use strict;
use warnings;
package Model::Study;
use File::Path qw/make_path/;
use File::Basename;
use Model::Design;
use YAML qw/DumpFile LoadFile/;
use Carp qw/confess/;
use List::Util qw/all max/;
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
      Model::Design::from_tsv($design_path), 
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
sub config_matches_design_checks {
  my ($config, $design) = @_;
  my %conditions_design = map {$_=>1} $design->all_conditions;
  my @conditions_config = keys %{$config->{condition_names}};
  my @conditions_match = map {("Condition $_ in config present in design" => $conditions_design{$_})} @conditions_config; 
  my @keys_match = map {("Key in config".join("\t", %{$_}). " matches slice") => $design->lookup_slice($_)} @{$config->{slices}};
  return @conditions_match, @keys_match;
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
  my $max_reps = max map {scalar @{$_}} values %{$self->{design}->runs_by_condition};
  return (
    {
      type => "aggregate_by_run",
      file_name => "$study_id.counts_per_run.tsv",
      title => "Counts per run",
      description => "Raw data (counts of aligned reads) for study $study_id",
      source => "counts_htseq2",
    },
    {
      type => "aggregate_by_run",
      file_name => "$study_id.tpm_per_run.tsv",
      title => "Expression per run (TPM)",
      description => "Gene expression in TPM for each run in  study $study_id",
      source => "tpm_htseq2",
    }, 
    ($max_reps >= 3 ? {
      type => "average_by_condition",
      file_name => "$study_id.tpm.tsv",
      title => "Expression per condition (TPM)",
      description => "Gene expression in TPM - median across runs per condition for study $study_id",
      source => "tpm_htseq2",
    } :()),
  );
}
1;
