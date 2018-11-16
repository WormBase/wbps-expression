use strict;
use warnings;
package Model::Study;
use File::Path qw/make_path/;
use File::Basename;
use Model::Design;
use YAML qw/DumpFile LoadFile/;
use Carp qw/confess/;
sub new {
  my ($class, $study_id, $study_design, $study_config) = @_;
  return bless {
    study_id => $study_id,
    design => $study_design,
    config => $study_config,
  }, $class;
}

sub from_folder {
   my ($class, $path) = @_;
   confess $path unless -d $path;
   my $study_id = basename($path);
   my $design_path = sprintf("%s/%s.tsv", $path,$study_id);
   my $config_path = sprintf("%s/%s.yaml", $path,$study_id);
   
   return unless -f $design_path and -f $config_path; 
   return $class->new($study_id, 
      Model::Design::from_tsv($design_path), 
      LoadFile($config_path)
   );
}

sub to_folder {
   my ($self, $path) = @_;
   make_path $path;
   $self->{design}->to_tsv(sprintf("%s/%s.tsv", $path, $self->{study_id}));
   DumpFile(sprintf("%s/%s.yaml", $path, $self->{study_id}),$self->{config});
}

sub config_matches_design {
  my ($config, $design) = @_;
  my %all_conditions = map {$_=>1} $design->all_conditions;
  for my $condition (keys %{$config->{condition_names}}){
     return 0 unless $all_conditions{$condition};
  }
  for my $key (@{$config->{slices}}){
     return 0 unless $design->lookup_slice($key);
  }
  return 1;
}
1;
