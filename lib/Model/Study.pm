use strict;
use warnings;
package Model::Study;
use File::Path qw/make_path/;
use File::Basename;
use Model::Design;
use YAML qw/DumpFile LoadFile/;
use Carp qw/confess/;
use List::Util qw/all/;
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
  my %checks = config_matches_design_checks(@_);
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
1;
