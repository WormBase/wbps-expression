use strict;
use warnings;
package Production::Sheets;
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
  my %result = -d $entry ? map {my $f = "$entry/$_/$_.tsv"; return -f $f ? ($_, $f) : () } : ();
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
  open(my $fh, ">", "$entry.tsv") or die "$entry.tsv: $!";
  map {print $fh "$_\n"} @{$list};
  close($fh);
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
