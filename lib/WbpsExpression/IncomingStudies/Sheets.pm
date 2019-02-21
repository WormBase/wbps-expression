use strict;
use warnings;
package WbpsExpression::IncomingStudies::Sheets;
use Text::CSV qw/csv/;
use File::Slurp qw/read_dir/;
use File::Path qw/remove_tree/;
use WbpsExpression::Model::Study;
use WbpsExpression::Model::SkippedRuns;

#use Smart::Comments;
sub new {
  my ($class, $src_dir) = @_;
  return bless {
     dir => "$src_dir/curation",
  }, $class;
}

sub read_directories {
  my ($self, $species) = @_;
  my @our_studies = map {WbpsExpression::Model::Study->from_folder($_)} $self->dir_content_paths("studies", $species);
  my @skipped_runs_in_our_studies;
  my @other_studies;
  for my $skipped_runs (map {WbpsExpression::Model::SkippedRuns->from_folder($_)} $self->dir_content_paths("skipped_runs", $species)){
    if(grep {$skipped_runs->{study_id} eq $_->{study_id}} @our_studies){
       push @skipped_runs_in_our_studies, $skipped_runs;
    } else {
       push @other_studies, $skipped_runs;
    }
  }
  return \@our_studies, \@skipped_runs_in_our_studies, \@other_studies; 
}

sub path {
  my ($self, @components) = @_;
  my $path = join("/", $self->{dir}, @components);
#### $path
  return $path;
}

sub dir_content_paths {
  my ($self, $name, $species) = @_;
  my $entry = $self->path($name, $species);
  return -d $entry ? map {"$entry/$_"} read_dir $entry : ();
}

sub tsvs_in_folders {
  my ($self, $name, $species) = @_;
  my $entry = $self->path($name, $species);
  my %result = -d $entry ? map {my $f = "$entry/$_/$_.tsv"; -f $f ? ($_, $f) : () } read_dir $entry : ();
  return \%result;
}

sub list {
  my ($self, $name, $species) = @_;
  my $entry = $self->path($name, $species);
  return -f "$entry.tsv" ? map {join "", @$_} @{csv(in=>"$entry.tsv")}: ();
}

sub write_list {
  my ($self, $list, $name, $species) = @_;
  my $entry = $self->path($name, $species);
  open(my $fh, ">:utf8", "$entry.tsv") or die "$entry.tsv: $!";
  map {print $fh "$_\n"} @{$list};
  close($fh);
}
sub remove_all {
  return remove_tree path(@_);
}
sub _aoa_to_double_hash {
  my ($aoa) = @_;
  my %result;
  for (@{$aoa}){
    my ($k1, $k2, $v1, @vs) = @$_;
    $result{$k1}{$k2} = @vs ? [$v1, @vs] : $v1;
  }
  return \%result;
}
sub double_hash {
  my ($self, $name, $species) = @_;
  my $entry = $self->path($name, $species);
  return &_aoa_to_double_hash(-f "$entry.tsv" ? csv(in=>"$entry.tsv", sep => "\t"): []);
}
1;
