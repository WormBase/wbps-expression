use strict;
use warnings;
package Curation::Sheets;
use File::Basename qw/dirname/;
use Text::CSV qw/csv/;
use File::Slurp qw/read_dir/;
#use Smart::Comments;
sub new {
  my ($class, $src_dir) = @_;
  return bless {
     dir => "$src_dir/curation",
  }, $class;
}

sub path {
   my ($self, @components) = @_;
   return join("/", $self->{dir}, @components);
}

sub list {
  my ($self, $name, $species) = @_;
  my $entry = $self->path($name, $species);
### $entry
  return -d $entry ? read_dir $entry : -f "$entry.tsv" ? map {join "", @$_} @{csv(in=>"$entry.tsv")}: ();
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
