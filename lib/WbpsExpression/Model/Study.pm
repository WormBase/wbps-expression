use strict;
use warnings;
package WbpsExpression::Model::Study;
use File::Path qw/make_path/;
use File::Basename;
use WbpsExpression::Model::Design;
use YAML qw/DumpFile LoadFile/;
use Carp qw/confess/;
use List::Util qw/all min pairmap/;
use List::MoreUtils qw/duplicates uniq/;
use open ':encoding(utf8)';
# use Smart::Comments '###';
sub new {
  my ($class, $study_id, $study_design, $study_config) = @_;
  return bless {
    study_id => $study_id,
    design => $study_design,
    config => $study_config,
  }, $class;
}

sub from_paths {
   my ($class, $study_id, $design_path, $config_or_path) = @_;
   return $class->new($study_id, 
      WbpsExpression::Model::Design::from_tsv($design_path), 
      ref $config_or_path eq "HASH" ? $config_or_path : LoadFile($config_or_path)
   );
}

sub from_folder {
   my ($class, $path) = @_;
   confess $path unless -d $path;
   my $study_id = basename($path);
   my $design_path = sprintf("%s/%s.tsv", $path,$study_id);
   my $config_path = sprintf("%s/%s.yaml", $path,$study_id);
   return unless -f $design_path and -f $config_path; 
   return $class->from_paths($study_id, $design_path, $config_path);
}

sub to_folder {
   my ($self, $path) = @_;
   make_path $path;
   $self->{design}->to_tsv(sprintf("%s/%s.tsv", $path, $self->{study_id}));
   DumpFile(sprintf("%s/%s.yaml", $path, $self->{study_id}),$self->{config});
}
sub config_matches_design {
  my ($self) = @_;
  my %checks = config_matches_design_checks($self->{config}, $self->{design});
  return all {$_} values %checks; 
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
sub consistency_checks {
  my ($self) = @_;
  return config_matches_design_checks($self->{config}, $self->{design});
}
sub passes_checks {
  my ($self) = @_;
  return $self->{design}->passes_checks && $self->config_matches_design;
}

# Currently not part of the config but wondrously inferred
# Logic: 
# - Everything should have counts per run
# - Everything should have TPM per run
# - Everything with three replicates somewhere should have a TPM per condition
#   + maybe with a "qc warning" 
#     #low replicates: egg(1), L3(2)
#     #low mean mapping quality: egg(56%)
#     !eggs, !L3 in the header
# - Every factor corresponds to a DE results file

# This needs to provide enough arguments to determine
# - which analysis to run (analysis code will have access to both the study and the data files)
# - how to display analysis results
# - how to link to the data files
sub analyses_required {
  my ($self) = @_;
  my $study_id = $self->{study_id};
  my $min_reps = min map {scalar @{$_}} values %{$self->{design}->runs_by_condition};
  my $counts_file_name = "$study_id.counts_per_run.tsv";
  my %h = %{$self->{design}{replicates_per_run}};
  my $has_technical_replicates = keys %h > uniq values %h;
  return (
    {
      type => "study_design",
      file_name => "$study_id.metadata_per_run.tsv",
      title => "Characteristics and conditions per run",
    },
    {
      type => "aggregate_by_run",
      file_name => $counts_file_name,
      title => "Raw data (counts of aligned reads) per run",
      description => "Raw data (counts of aligned reads) for study $study_id",
      source => "counts_htseq2",
      decorate_files => 0,
    },
    {
      type => "aggregate_by_run",
      file_name => "$study_id.tpm_per_run.tsv",
      title => "Gene expression (TPM) per run",
      description => "Gene expression in TPM for each run in study $study_id",
      source => "tpm_htseq2",
      decorate_files => 1,
    }, 
    ($min_reps >= 2 ? {
      type => "average_by_condition",
      file_name => "$study_id.tpm.tsv",
      title => "Gene expression (TPM) per condition as median across ".($has_technical_replicates ? "replicates" : "runs"),
      description => sprintf ("Gene expression in TPM - %s per condition for study $study_id", $has_technical_replicates ? "technical, then biological replicates" : "runs"),
      source => "tpm_htseq2",
      decorate_files => 1,
    } :()),
   (map {{
      type => "differential_expression",
      file_name => sprintf("$study_id.de%s%s.tsv" , $_->{name} ? ".": "" , $_->{name}),
      title => sprintf("Differential expression%s%s" , $_->{name} ? ": ": "", $_->{name} =~ tr/_/ /r),
      description => sprintf("Differential expression analysis for study $study_id%s. Values are base 2 log of maximum likelihood estimate of fold change and adjusted p-value by Wald test, given to 2 s.f., where adj_pval < 0.05 and abs(log2fc) > 0.5", ($_->{name} ? ", comparing across ".$_->{name} =~ tr/_/ /r : "")),
      source_file_name => $counts_file_name,
      contrasts => $_->{values}, 
   }} @{$self->{config}{contrasts}}),
  );
}
1;
