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
   return $class->new($study_id, 
      Model::Design::from_tsv(sprintf("%s/%s.tsv", $path,$study_id)), 
      LoadFile(sprintf("%s/%s.yaml", $path,$study_id))
   );
}


sub to_folder {
   my ($self, $path) = @_;
   make_path $path;
   $self->{design}->to_tsv(sprintf("%s/%s.tsv", $path, $self->{study_id}));
   DumpFile(sprintf("%s/%s.yaml", $path, $self->{study_id}),$self->{config});
}
1;
