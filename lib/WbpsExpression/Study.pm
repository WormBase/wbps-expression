use strict;
use warnings;
package WbpsExpression::Study;
use File::Path qw/make_path/;
use File::Basename;
use WbpsExpression::Study::Design;
use YAML qw/DumpFile LoadFile/;
use Carp qw/confess/;
use List::Util qw/all min pairmap pairs/;
use List::MoreUtils qw/duplicates uniq/;
use File::Slurp qw/write_file/;
use Regexp::Common qw/URI/;
use Data::Dumper;
use LWP;
use Log::Any '$log';

use open ':encoding(utf8)';
# use Smart::Comments '###';
sub new {
  my ($class, $study_id, $species, $brc4, $study_design, $study_config, $skipped_runs, $sources) = @_;
  return bless {
    study_id => $study_id,
    species => $species,
    brc4   => $brc4,
    design => $study_design,
    config => $study_config,
    skipped_runs => $skipped_runs,
    sources => $sources,
  }, $class;
}
#
# sub rnaseq_status {
#   my ($species, $wbps_assembly) = @_;
#   my $rnaseq_status_path = $ENV{PARASITE_CONF} . "/" . "WBPS".$ENV{PARASITE_VERSION}.".rnaseq_status.json";
#   open(FH,"<",$rnaseq_status_path) or die "$rnaseq_status_path file doesn't exist!\n";
#   my $data = do { local $/; <FH> };;
#   my $ret = JSON::decode_json( $data );
#   if($ret->{$species}) {
#     my $status = (exists($ret->{$species}) && $ret->{$species}->{$wbps_assembly}->{brc4}) ? "brc4" : "irap";
#     return $status;
#   } else {
#     return "irap"
#   }
# }

sub rnaseq_status {
    my ($brc4_dir, $species, $wbps_assembly) = @_;
    if (-d $brc4_dir . '/' . $species . '/' . $wbps_assembly) {
        return "brc4";
    }
    else {
        return "irap";
    }
}


sub from_folder {
  my ($class, $path, $species, $brc4_path) = @_;
  return unless -d $path;
  my $study_id = basename($path);
  my $design_path = sprintf("%s/%s.design.tsv", $path,$study_id);
  my $config_path = sprintf("%s/%s.config.yaml", $path,$study_id);
  my $skipped_runs_path = sprintf("%s/%s.skipped_runs.tsv", $path,$study_id);
  my $sources_path = sprintf("%s/%s.sources.tsv", $path,$study_id);
  my $brc4_study_path = ($brc4_path eq "") ? "" : "$brc4_path/$study_id";
  my $brc4 = ($brc4_path eq "") ? 0 : 1;
  return unless -f $config_path;
  return $class->new(
     $study_id,
     $species,
     $brc4,
     read_design($design_path),
     LoadFile($config_path),
     read_skipped_runs($skipped_runs_path),
     read_sources($species, $path, $brc4_study_path),
  );
}

sub to_folder {
  my ($self, $path) = @_;
  my $study_id = $self->{study_id};
  my $design_path = sprintf("%s/%s.design.tsv", $path,$study_id);
  my $config_path = sprintf("%s/%s.config.yaml", $path,$study_id);
  my $skipped_runs_path = sprintf("%s/%s.skipped_runs.tsv", $path,$study_id);
  my $sources_path = sprintf("%s/%s.sources.tsv", $path,$study_id);
  make_path $path;
  $self->{design}->to_tsv($design_path);
  write_design($design_path, $self->{design});
  DumpFile($config_path, $self->{config});
  write_skipped_runs($skipped_runs_path, $self->{skipped_runs});
  write_sources($sources_path, $self->{sources});
}

sub read_design {
  my ($path) = @_;
  return -s $path ? WbpsExpression::Study::Design::from_tsv($path) : WbpsExpression::Study::Design::empty;
}

sub write_design {
  my ($path, $design) = @_;
  unlink $path;
  return if $design->is_empty;
  $design->to_tsv($path);
}

sub read_skipped_runs {
  my ($path) = @_;
  return [] unless -s $path;
  open(my $fh, "<", $path) or die "$!: $path";
  my ($header, @runs) = <$fh>;
  die "No header?: $path" if $header ne "Run\n";
  chomp for @runs;
  return \@runs;
}

sub write_skipped_runs {
  my ($path, $runs) = @_;
  unlink $path && return unless @{$runs // []};
  write_file($path, join("\n", "Run", sort @{$runs})."\n");
}

sub read_sources {
  my ($species, $path, $brc4_path) = @_;
  return {} unless -s $path;
  if($brc4_path ne ""){
    my %sources = read_brc4_sources($brc4_path);
    return \%sources;
  }
  open(my $fh, "<", $path) or die "$!: $path";
  my ($header, @xs) = <$fh>;
  warn "Study::read_sources bad header at $path" if $header ne "Run\tLocation\tEnd\tQuality\n";
  my %sources = map {
    chomp;
    my ($run, $location, $end, $quality) = split "\t";
    $run => {
      location => $location,
      end => $end,
      quality => $quality,
    }
  } @xs;
  return \%sources;
}

sub read_brc4_sources {
  my ($study_path) = @_;
  return {} unless -s $study_path;
  opendir(my $dh, $study_path) || die "Can't open $study_path: $!";
  my @pre_directories = readdir $dh;
  @pre_directories = grep !m/^\./, @pre_directories;
  my @directories;
  for my $pre_dir (@pre_directories) {
    push @directories, $pre_dir if (-d "$study_path/$pre_dir");
  }
  closedir($dh);
  my %sources = map {
    my $run_dir = $_;
    my $run_id = $run_dir;

    # parse and pretty print mapping stats (mapping quality) from the brc4 pipeline output mappingStats.txt file.
    my $mapping_stats_file = "$study_path/$run_dir/mappingStats.txt";
    open(my $mapping_fh, "<", $mapping_stats_file) or die "$!: $mapping_stats_file";
    my ($header, @xs) = <$mapping_fh>;
    warn "Study::read_sources bad header at $mapping_stats_file" if $header ne "file\tcoverage\tmapped\tnumber_reads_mapped\taverage_read_length\tnumber_pairs_mapped\n";
    my @result_bam_lines = grep {/^results\.bam.+$/} @xs;
    my $result_bam_lines_count = @result_bam_lines;
    die "Study::read_sources There are multiple or 0 lines with the '/^results\.bam.+\$/' pattern in the $mapping_stats_file" unless $result_bam_lines_count == 1;
    my $result_bam_line = $result_bam_lines[0];
    my ($file, $coverage, $mapped, $number_reads_mapped, $average_read_length, $number_pairs_mapped) = split "\t", $result_bam_line;
    my $pretty_mapped = int($mapped * 100);
    my $location = "$study_path/$run_id";
    close($mapping_fh);

    # parse and pretty print paired end metadata info from the brc4 pipeline output metadata.json file.
    my $metadata_file = "$study_path/$run_dir/metadata.json";
    open(my $metadata_fh, "<", $metadata_file) or die "$!: $metadata_file";
    my $data = do {
      local $/;
      <$metadata_fh>
    };
    my $metadata = JSON::decode_json($data);
    my $end = $metadata->{hasPairedEnds};
    die "$metadata_file does not contain the hasPairedEnds field or the hasPairedEnds value is not true/false." unless ref($end) eq "JSON::PP::Boolean";
    my $pretty_end = ($end) ? "pe" : "se";
    close($metadata_fh);

    # add it all in a hash
    $run_id => { location => $location,
                 end      => $pretty_end,
                 quality  => $pretty_mapped }
  } @directories;

  return %sources;
}



sub write_sources {
  my ($path, $sources) = @_;
  write_file($path, join (
    "\n",
    "Run\tLocation\tEnd\tQuality", 
    map {
      join("\t", $_, $sources->{$_}{location}, $sources->{$_}{end}, $sources->{$_}{quality} )
    } sort keys %{$sources}
  )."\n"); 
}

sub list_of_contrasts_checks {
  my ($num_replicates_by_condition, $name, @values) = @_;
#### list_of_contrasts_checks: @_
  if ($#values > 55 ){
    return ("$name manageably many contrasts" => 0);
  }
  my @ds = duplicates map {$_->[2]} @values;
  if (@ds) {
    return ("contrast names should be all unique, duplicates: ".join (", ", @ds) => 0);
  }
  return map { 
    my ($reference, $test, $contrast_name) = @{$values[$_]};
    my $nrr = $num_replicates_by_condition->{$reference};
    my $nrt = $num_replicates_by_condition->{$test};
    $contrast_name && $nrr && $nrt
      ? ("$contrast_name should have enough replicates, is: $nrr ref $nrt test"
            => not ($nrr < 3 && $nrt < 3))
      : ("contrast $reference/$test/$contrast_name not matching the design" 
            => 0 )
  }  (0 .. $#values);
  
}
sub config_matches_design_checks {
  my ($config, $design) = @_;
#### contrasts: $config->{contrasts}
  my %num_replicates_by_condition = pairmap { $a => scalar @$b } %{$design->replicates_by_condition};
  my @conditions_config = keys %{$config->{condition_names}};
  my @conditions_match = map {("Condition $_ in config present in design" => $num_replicates_by_condition{$_})} @conditions_config; 
  my @contrasts_match = map {
    list_of_contrasts_checks(\%num_replicates_by_condition, $_->{name}, @{$_->{values}})
  } @{$config->{contrasts}};
#### @contrasts_match
  return @conditions_match, @contrasts_match;
}
sub sources_checks {
  my ($self) = @_;
  my @ftp_paths = map {
    my ($k, $v) = @{$_};
    my $location = $v->{location};
    "ISL FTP folder: $k => $location" => $location =~ m/$RE{URI}{FTP}/
  } pairs %{$self->{sources} // {}};
  my @bad_ends = grep {$_ ne 'pe' and $_ ne 'se'} map {$_->{end}} values %{$self->{sources} // {}};
  if (scalar @ftp_paths % 2 == 0) {
    return @ftp_paths, "Ends single/paired => " . (@bad_ends ? 0 : 1), "Ends single/paired => " . (@bad_ends ? 0 : 1);
  }
  else {
    return @ftp_paths, "Ends single/paired => " . (@bad_ends ? 0 : 1);
  }
}

sub all_runs {
  my ($self) = @_;
  my @runs_design = $self->{design}->all_runs;
  my @runs_skipped = @{$self->{skipped_runs} //[]};
  return sort @runs_design, @runs_skipped;
}
sub run_ids_checks {
  my ($self) = @_;
  my @runs_design = $self->{design}->all_runs;
  my @runs_skipped = @{$self->{skipped_runs} //[]};
  my @runs_sources = keys %{$self->{sources}};
  my @runs_design_no_sources = grep {not $self->{sources}{$_}} @runs_design;
  my @runs_skipped_no_sources = grep {not $self->{sources}{$_}} @runs_skipped;
  return (
     "All runs in design have a source dir " . join(", ", @runs_design_no_sources) => ! @runs_design || ! @runs_design_no_sources,
     "All skipped runs have a source dir " . join(", ", @runs_skipped_no_sources) => ! @runs_skipped || ! @runs_skipped_no_sources,
     "No more source dirs" => @runs_design + @runs_skipped == @runs_sources,
     "No duplicates between design and skipped" => not duplicates (@runs_design, @runs_skipped),
  );
  #todo: design + skipped_runs = results
}
sub consistency_checks {
  my ($self) = @_;
  return config_matches_design_checks($self->{config}, $self->{design});
}
sub design_checks {
  my ($self) = @_;
  return "Design empty" => 1 if $self->{design}->is_empty;
  return $self->{design}->data_quality_checks;
}

sub all_checks {
  my ($self) = @_;
  my %result;
  %result = (%result, $self->design_checks, $self->consistency_checks, $self->sources_checks, $self->run_ids_checks);
  return %result;
}
sub config_matches_design {
  my ($self) = @_;
  my %checks = $self->consistency_checks;
  return all {$_} values %checks; 
}
sub passes_checks {
  my ($self) = @_;
  my %checks = $self->all_checks;
  return all {$_} values %checks;
}

# use featurecounts results if they are available, if not use HTSeq2
sub quantification_method {
  my ($self) = @_;
  if($self->{brc4}){
    return "HTSeq2";
  }
  my $study_id          = $self->{study_id};
  my $species           = $self->{species}; 
  my $ftp               = 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena';

  my $featurecounts_raw = join("/", $ftp, $study_id, $species, 'genes.raw.featurecounts.tsv'); 
  my $response_fc_raw   = LWP::UserAgent->new->get($featurecounts_raw); 
  
  my $featurecounts_tpm = join("/", $ftp, $study_id, $species, 'genes.tpm.featurecounts.tsv');
  my $response_fc_tpm   = LWP::UserAgent->new->get($featurecounts_tpm);

  my $featurecounts_tpm_irap = join("/", $ftp, $study_id, $species, 'genes.tpm.featurecounts.irap.tsv');
  my $response_fc_tpm_irap   = LWP::UserAgent->new->get($featurecounts_tpm_irap);

  if ($response_fc_raw->is_success && $response_fc_tpm->is_success
      || $response_fc_raw->is_success && $response_fc_tpm_irap->is_success){
    return "FeatureCounts";
  }

  my $htseq_raw = join("/", $ftp, $study_id, $species, 'genes.raw.htseq2.tsv');
  my $response_htseq_raw   = LWP::UserAgent->new->get($htseq_raw);

  my $htseq_tpm = join("/", $ftp, $study_id, $species, 'genes.tpm.htseq2.tsv');
  my $response_htseq_tpm   = LWP::UserAgent->new->get($htseq_tpm); 

  my $htseq_tpm_irap = join("/", $ftp, $study_id, $species, 'genes.tpm.htseq2.irap.tsv');
  my $response_htseq_tpm_irap   = LWP::UserAgent->new->get($htseq_tpm_irap);

  if ($response_htseq_raw->is_success && $response_htseq_tpm->is_success
      || $response_htseq_raw->is_success && $response_htseq_tpm_irap->is_success){
    return "HTSeq2";
  } else{
    return "FeatureCounts"; }
}

sub source_counts {
  my ($self, $run_id) = @_;
  if ($self->{brc4}) {
    my $simple_stranded = join("/", $self->{sources}->{$run_id}{location}, "genes.htseq-union.stranded.sum.counts");
    my $simple_unstranded = join("/", $self->{sources}->{$run_id}{location}, "genes.htseq-union.unstranded.counts");
    if (-e $simple_stranded) {
      return $simple_stranded;
    } elsif (-e $simple_unstranded) {
      return $simple_unstranded;
    } else {
      die "Couldn't find a TPM file in $self->{sources}->{$run_id}{location} $simple_stranded or $simple_unstranded"
    };
  } else {
    my $m = lc $self->{quantification_method};
    return join("/", $self->{sources}->{$run_id}{location}, "$run_id.$self->{sources}->{$run_id}{end}.genes.raw.$m.tsv");
  }
}

# sometimes TPM files are named .irap.tsv, sometimes just .tsv. Try both.

sub source_tpm {
  my ($self, $run_id) = @_;
  if ($self->{brc4}) {
    my $simple_stranded = join("/", $self->{sources}->{$run_id}{location}, "genes.htseq-union.stranded.sum.counts.tpm2");
    my $simple_unstranded = join("/", $self->{sources}->{$run_id}{location}, "genes.htseq-union.unstranded.counts.tpm2");
    if (-e $simple_stranded) {
      return $simple_stranded;
    }
    elsif (-e $simple_unstranded) {
      return $simple_unstranded;
    }
    else {
      die "Couldn't find a count file in $self->{sources}->{$run_id}{location}"
    };
  } else {
    my $m = lc $self->{quantification_method};
    my $simple = join("/", $self->{sources}->{$run_id}{location}, "$run_id.$self->{sources}{$run_id}{end}.genes.tpm.$m.tsv");
    my $simple_response = LWP::UserAgent->new->get($simple);
    if ($simple_response->is_success) {
      return $simple;
    }
    return join("/", $self->{sources}->{$run_id}{location}, "$run_id.$self->{sources}->{$run_id}{end}.genes.tpm.$m.irap.tsv");
  }
}

sub source_bigwig {
  my ($self, $run_id) = @_;
  if ($self->{brc4}) {
    print "$run_id $self->{sources}->{$run_id}{location}\n";
    my $bigwig_path = join("/", $self->{sources}->{$run_id}{location}, "results.bw");
    if (! -e $bigwig_path) {
      die "Couldn't find a bigwig file in $self->{sources}->{$run_id}{location}"
    }
    return $bigwig_path;
  }
  return join("/", $self->{sources}->{$run_id}{location}, "$run_id.nospliced.bw");
}

sub mapping_quality {
  my ($self, $run_id) = @_;
  return $self->{sources}->{$run_id}{quality};
}

sub qc_issues_per_run {
  my ($self) = @_;
  my %result;
  for my $run_id ($self->{design}->all_runs){
    my $q = $self->{sources}->{$run_id}{quality};
    if ($q < 40){
      $result{$run_id} = ["Low mapping quality: $q"];
    }
  } 
  return \%result;
}

sub study_frontmatter {
  my ($self) = @_;
  my $quantification_method = $self->{quantification_method};
  my $aligner = $self->{brc4} ? "HISAT2" : "TopHat2";
  my $method = $self->{brc4} ? "by Ensembl Metazoa RNA-Seq pipeline" : "by RNASeq-er: https://www.ebi.ac.uk/fg/rnaseq/api";
  my @pubmed_lines = map {sprintf ("%s: https://www.ncbi.nlm.nih.gov/pubmed/%s", $_->[1][0], $_->[0]) } pairs %{$self->{config}{pubmed} // {}};
  return join ("\n",
    "",
    sprintf("Study %s: %s", $self->{study_id}, $self->{config}{title}),
    (@pubmed_lines
      ? "See ". join(", ", @pubmed_lines)
      : sprintf("Submitted to archives by %s", $self->{config}{submitting_centre})
    ),
    "Reads aligned with $aligner and quantified with $quantification_method, $method",
    "",
  );
}
sub analyses_required {
  my ($self) = @_;
#### $self
  my $study_id = $self->{study_id};
  my $min_reps = min map {scalar @{$_}} values %{$self->{design}->runs_by_condition};
  my $counts_file_name = "$study_id.counts_per_run.tsv";
  my %h = %{$self->{design}{replicates_by_run}};
  my $has_technical_replicates = keys %h > uniq values %h;
  return (
    {
      type => "study_design",
      file_name => "$study_id.metadata_per_run.tsv",
      title => "Characteristics and conditions per run",
    },
    {
      type => "counts_per_run",
      file_name => $counts_file_name,
      title => sprintf("Counts of aligned reads per run (%s)", $self->quantification_method),
    },
    {
      type => "tpms_per_run",
      file_name => "$study_id.tpm_per_run.tsv",
      title => "Gene expression (TPM) per run",
      description => join("\n",
        $self->study_frontmatter,
        "Values are transcripts per million units (TPMs) per gene for each run",
      ),
    }, 
    ($min_reps >= 2 ? {
      type => "average_by_condition",
      file_name => "$study_id.tpm.tsv",
      title => "Gene expression (TPM) per condition as median across ".($has_technical_replicates ? "replicates" : "runs"),
      description => join("\n",
        $self->study_frontmatter,
        "Values are transcripts per million units (TPMs) per gene, median across ". ($has_technical_replicates ? "technical, then biological replicates" : "runs"),
      ),
    } :()),
   (@{$self->{config}{contrasts}} ? {
      type => "differential_expression",
      title => "Differential expression",
      files => [
        (map {
         {
           file_name => sprintf("$study_id.de.%s.tsv" , $_->{name}),
           title => sprintf("Summary file: %s" , $_->{name} =~ tr/_/ /r),
           key => $_->{name},
           description => join("\n",
            $self->study_frontmatter,
            "Differential expression analysis, comparing pairs of conditions differing by " . $_->{name} =~ tr/_/ /r,
            ),
          }} @{$self->{config}{contrasts}}
         ),{
           file_name => "$study_id.de.contrasts.zip",
           title => sprintf("Full result files for %s contrasts (zipped)", scalar map {@{$_->{values}}} @{$self->{config}{contrasts}}),
           key => "all_contrasts",
         }
      ],
      source_file_name => $counts_file_name,
      contrasts => {map {$_->{name} => $_->{values}} @{$self->{config}{contrasts}}},
   } : ()), 
  );
}

sub to_hash {
  my ($self) = @_;
  my %result;
  $result{study_id} = $self->{study_id};
  $result{study_title} = $self->{config}{title};
  $result{study_category} = $self->{config}{category};
  $result{attributes} = {
    submitting_centre => $self->{config}{submitting_centre},
    "ENA first public" => $self->{config}{ena_first_public},
    "ENA last update" => $self->{config}{ena_last_update},
    "ENA study" => sprintf('<a href="https://www.ebi.ac.uk/ena/data/view/%s">Study page: %s</a>', $self->{study_id},$self->{study_id}),
    %{$self->{design}{values}{common}}
  };
  $result{attributes}{"Study description"} = $self->{config}{description}
    if $self->{config}{description};
  $result{attributes}{pubmed} = join (", ", sort {$b cmp $a} map {
      my ($pubmed_id, $xs) = @{$_};
      my (undef, $pubmed_description) = @{$xs};
      sprintf('<a href="https://www.ncbi.nlm.nih.gov/pubmed/%s">%s</a>: %s', $pubmed_id, $pubmed_id, $pubmed_description)
    } pairs %{$self->{config}{pubmed}})
    if %{$self->{config}{pubmed}};

  $result{attributes}{"ENA BioProject"} = sprintf('<a href="https://www.ebi.ac.uk/ena/data/view/%s">Study page: %s</a>', $self->{config}{bioproject},$self->{config}{bioproject})
    if $self->{config}{bioproject};
  for my $r (@{$self->{config}{resources}}){
    my ($property_name, $label, $url) = @{$r};
    $result{attributes}{"Linked resource: $property_name"} = sprintf('<a href="%s">%s</a>', $url, $label);
  }
  my @runs_curated = map {
    my ($condition, $replicate, $run_id) = @{$_};
    my %o;
    $o{condition} = $condition;
    $o{replicate} = $replicate if $replicate ne $run_id;
    $o{run_id} = $run_id;
    $o{bigwig} = $self->source_bigwig($run_id);
    my %attributes = map {
      $_ => $self->{design}->value_in_run($run_id, $_)
    } $self->{design}->characteristics_not_common;
    $o{attributes} = \%attributes;
    \%o
  } $self->{design}->condition_replicate_run_ordered_triples;

  my @runs_skipped = map {
    { run_id => $_, bigwig => $self->source_bigwig($_)}
  } @{$self->{skipped_runs}};
  $result{runs} = [@runs_curated, @runs_skipped];
  return \%result;
}
1;
