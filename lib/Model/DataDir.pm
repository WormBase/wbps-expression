use strict;
use warnings;
use File::Basename qw/basename/;
use File::Slurp qw/read_dir/;
package Model::DataDir;
# I am not sure what this can even encapsulate.
# Naming convention for files?
# Utility methods?
# - list done analyses
# - analysis code should know where to put stuff
# 
# - get output path for analysis / did analysis get done?


sub new {
  my ($class, $path) = @_;
  return bless {dir => $path, study_id => basename $path}, $class;
}

sub path {
  my ($self, $type) = @_;
  my $dir = $self->{dir};
  my $study_id = $self->{study_id};
  return "$dir/$study_id.$type.tsv";
}

sub list {
  my ($self) = @_;
  my $study_id = $self->{study_id};
  return map {
    (my $type = basename $_) =~ s/^$study_id\.(.*)\.tsv/$1/;
    ($type => "$_")
  } read_dir($self->{dir}, prefix => 1);
}
1;
