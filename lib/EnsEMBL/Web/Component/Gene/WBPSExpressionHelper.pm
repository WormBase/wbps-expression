package EnsEMBL::Web::Component::Gene::WBPSExpressionHelper;
use File::Basename;

sub new {
  my ($class, $args) = @_;
  return bless $args, $class;
}

sub from_folder {
  my ($self, $species, $studies_path) = @_;
  my ($spe, $cies, $bp) = split "_", $species;
  
  my $studies_file = "$studies_path/$spe"."_"."$cies.studies.tsv";
  my %metadata;
  open (my $fh, "<", $studies_file) or return;
  while(<$fh>){
    chomp;
    my ($accession, $category, $title) = split "\t";
    $metadata{$accession}{title} = $title;
    $metadata{$accession}{category} = $category;
  }
  my @studies;
  for my $study_path (grep {-d $_}  glob("$studies_path/*")){
    my $study_id = basename $study_path;
    my $study = {
      study_id => $study_id,
      study_title => $metadata{$study_id}{title},
      study_category => $metadata{$study_id}{category},
    };
    if ($metadata{$study_id}{category} eq "Response to treatment"){
      my @files = glob("$study_path/$study_id.de.*.tsv");
      next unless @files;
      for my $file_path (@files){
        my ($type) = $file_path =~ m{$study_path/$study_id\.de\.(.*)\.tsv};
        $study->{de}{$type} = $file_path;
      }
    } elsif ($metadata{$study_id}{category} eq "Other" ) {
      my $tpms_per_run = "$study_path/$study_id.tpm_per_run.tsv";
      next unless -s $tpms_per_run;
      $study->{tpms_per_run} = $tpms_per_run;
    } else {
      my $tpms_per_condition = "$study_path/$study_id.tpm.tsv";
      next unless -s $tpms_per_condition;
      $study->{tpms_per_condition} = $tpms_per_condition;
    }
    push @studies, $study;
  }
  return &new(__PACKAGE__, {
    studies => \@studies,
    species_display_name => sprintf("%s %s (%s)", ucfirst($spe), $cies, uc($bp)),
    species_url_name => join("_", map {lc($_)} ($spe, $cies, $bp)),
  });
}


sub render_page {
  my ($self, $gene_id, $category) = @_;
  my @studies = grep {$_->{study_category} eq $category} @{$self->{studies}};
  my @html_panes = response_as_html_panes($self->{species_url_name}, $gene_id, $category, \@studies);
  return (
   join("<br>", @html_panes) || html_no_results($self->{species_display_name}, $gene_id, $category)
  );
}
sub response_as_html_panes {
  my ($species, $gene_id, $category, $studies) = @_;
  return unless @$studies;
### response_as_html_panes: join("\t", $species , $gene_id , $category,  map {$_->{study_id}} @{$studies})
  if ($category eq  "Response to treatment"){
    my ($differential_expression_values, $studies_with_no_results) = list_of_differential_expression_values_in_studies_and_studies_with_no_results($species, $gene_id, $studies);
    return (@{$differential_expression_values} ? (html_differential_expression_values_table($differential_expression_values), @{$studies_with_no_results} ? html_studies_with_no_results($studies_with_no_results) : ()): ());
  } elsif ($category eq "Other") {
    return map {html_flat_horizontal_table($_->{column_headers}, $_->{values})} summary_stats_in_tables($species, $gene_id, $studies);
  } else {
    return map {html_flat_horizontal_table($_->{column_headers}, $_->{values})} tpms_in_tables($species, $gene_id, $studies);
  }

}
sub html_no_results {
  my ($species, $gene_id, $category) = @_;
  return "<p> $species: no results for gene $gene_id in category $category </p>";
}
sub html_header {
  my ($species, $gene_id, $category) = @_;
  return "<h2> Gene expression - $category </h2>"; 
}
sub html_differential_expression_values_table {
  my ($differential_expression_values) = @_;
  return (
       "<table>"
     . "<thead>"
     . "<tr>"
        . "<th>Study</th>"
        . "<th>Contrast</th>"
        . "<th>Log<sub>2</sub>-fold change</th>"
        . "<th>Adjusted p-value</th>"
     . "</tr>"
     . "</thead>"
     . "<tbody>"
     . join ("\n", map {
      "<tr>"
        . "<td>" . html_study_link($_) . "</td>"
        . "<td>" . $_->{contrast} . "</td>"
        . "<td>" . $_->{log2_fold_change} . "</td>"
        . "<td>" . $_->{adjusted_p_value} . "</td>"
      . "</tr>"
     } @{$differential_expression_values})
     . "</tbody>"
     . "</table>"
  );
}
sub html_flat_horizontal_table {
  my ($column_headers, $values) = @_;
  return (
       "<table>"
     . "<thead>"
     . "<tr>"
        . join("\n", map {
        "<th>$_</th>"
        } @{$column_headers})
     . "</tr>"
     . "</thead>"
     . "<tbody>"
     . "<tr>"
        . join("\n", map {
        "<td>$_</td>"
        } @{$values})
     . "</tr>"
     . "</tbody>"
     . "</table>"
  );
}
sub html_studies_with_no_results {
  my ($studies) = @_;
  return "<ul>" . join ("<br>",
    map { "<li>$_</li>" }
    map { html_study_link($_) }
    @{$studies}
   ) . "</ul>";
}
sub html_study_link {
  my ($study) = @_;
  return sprintf("<a href=\"%s\">%s</a>", $study->{study_url}, $study->{study_title});
}
sub study_url {
  my ($species, $study_id) = @_;
  return "/expression/$species/index.html#$study_id";
}
sub list_of_differential_expression_values_in_studies_and_studies_with_no_results {
  my ($species, $gene_id, $studies) = @_;
  my @differential_expression_values;
  my @studies_with_no_results;

  for my $study (@{$studies}){
     $study->{study_url} = study_url($species,$study->{study_id}),
     my @differential_expression_values_for_study;
     while (my ($type, $path) = each %{$study->{de}}) {
        my ($contrasts, $differential_expression_values) = search_in_file($path, $gene_id);
        C:
        for my $i (0..$#$contrasts){
           my $contrast = $contrasts->[$i];
           next C if $contrast =~ /^\!/; # Low replicates or failed QC
           my ($log2_fold_change, $adjusted_p_value) = split(" ", $differential_expression_values->[$i]);
           push @differential_expression_values_for_study, {
              study_url => $study->{study_url},
              study_title => $study->{study_title},
              contrast => $contrasts->[$i],
              log2_fold_change => $log2_fold_change,
              adjusted_p_value => $adjusted_p_value,
           }
        }
     }
     push @differential_expression_values, $_ for @differential_expression_values_for_study;
     push @studies_with_no_results, $study unless @differential_expression_values_for_study;
  }
  return \@differential_expression_values, \@studies_with_no_results;
}

# Adapted from: https://metacpan.org/source/SHLOMIF/Statistics-Descriptive-3.0612/lib/Statistics/Descriptive.pm#L620
sub quantile {
    my ( $data_sorted, $QuantileNumber ) = @_;
    return $data_sorted->[0] if ( $QuantileNumber == 0 );
    my $count = @{$data_sorted};
    return $data_sorted->[ $count - 1 ] if ( $QuantileNumber == 4 );
 
    my $K_quantile = ( ( $QuantileNumber / 4 ) * ( $count - 1 ) + 1 );
    my $F_quantile = $K_quantile - POSIX::floor($K_quantile);
    $K_quantile = POSIX::floor($K_quantile);
 
    my $aK_quantile     = $data_sorted->[ $K_quantile - 1 ];
    return $aK_quantile if ( $F_quantile == 0 );
    my $aKPlus_quantile = $data_sorted->[$K_quantile];
 
    my $quantile = $aK_quantile
      + ( $F_quantile * ( $aKPlus_quantile - $aK_quantile ) );
 
    return $quantile;
}
sub summary_stats_in_tables {
  my ($species, $gene_id, $studies) = @_;
### summary_stats_in_tables: $studies
  my @result;
  for my $study (@{$studies}){
     my ($runs, $expression_tpm) = search_in_file($study->{tpms_per_run}, $gene_id);
### $expression_tpm
     my @expression_tpm_sorted = sort @{$expression_tpm //[]};
     push @result, {
       study_url => study_url($species,$study->{study_id}),
       study_title => $study->{study_title},
       column_headers => ["N", "min", "Q1", "Q2", "Q3", "max"],
       values => [ 
          scalar @expression_tpm_sorted,
          quantile(\@expression_tpm_sorted, 0),
          sprintf("%.1f",quantile(\@expression_tpm_sorted, 1)),
          sprintf("%.1f",quantile(\@expression_tpm_sorted, 2)),
          sprintf("%.1f",quantile(\@expression_tpm_sorted, 3)),
          quantile(\@expression_tpm_sorted, 4),
       ],
     } if $runs and $expression_tpm;
  }
  return @result;
}
sub tpms_in_tables {
  my ($species, $gene_id, $studies) = @_;
### tpms_in_tables: join("\t", $species , $gene_id , map {$_->{study_id}} @{$studies})
  my @result;  
  for my $study (@{$studies}){
     my ($conditions, $expressions_tpm) = search_in_file($study->{tpms_per_condition}, $gene_id);
     my @conditions_warnings_as_text = map {s{^!\s*(.*)}{$1 <i>(has quality warnings)</i>}; $_} @{$conditions};
     push @result, {
       study_url => study_url($species,$study->{study_id}),
       study_title => $study->{study_title},
       column_headers => \@conditions_warnings_as_text,
       values => $expressions_tpm,
     } if @conditions_warnings_as_text and @{$expressions_tpm};
  }
  return @result;
}


sub search_in_file {
  my ($path, $gene_id) = @_;
### search_in_file: $path . " " . $gene_id
  my $l = `grep --max-count=1 "^$gene_id" $path`;
### $l
  chomp $l;
  return unless $l;
  my ($id, @xs) = split "\t", $l;
  return unless $id eq $gene_id and @xs;
  my $h = `grep --max-count=1 "^\t" $path`;
### $h
  chomp $h;
  return unless $h;
  my ($blank, @hs) = split "\t", $h;
  return unless not($blank) and @hs;
  return unless @hs == @xs;
  return \@hs, \@xs;
}

1;