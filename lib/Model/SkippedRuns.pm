use strict;
use warnings;
package Model::SkippedRuns;
use File::Path qw/make_path/;
use File::Basename;
use File::Slurp qw/read_file write_file/;
use YAML qw/DumpFile LoadFile/;
use Carp qw/confess/;
use open ':encoding(utf8)';

sub new {
  my ($class, $study_id, $study_config, $runs) = @_;
  return bless {
    study_id => $study_id,
    config => $study_config,
    runs => $runs,
  }, $class;
}

sub from_paths {
   my ($class, $study_id, $runs_path, $config_or_path) = @_;
   open(my $runs_fh, "<", $runs_path) or die "$!: $runs_path";
   my ($header, @runs) = <$runs_fh>;
   chomp for @runs;
   my $study_config = ref $config_or_path eq "HASH" ? $config_or_path : LoadFile($config_or_path);
   return $class->new($study_id, $study_config, \@runs);
}

sub from_folder {
   my ($class, $path) = @_;
   confess $path unless -d $path;
   my $study_id = basename($path);
   my $runs_path = sprintf("%s/%s.tsv", $path,$study_id);
   my $config_path = sprintf("%s/%s.yaml", $path,$study_id);
   return unless -f $runs_path;
   my $config_or_path = -f $config_path ? $config_path : {};
   return $class->from_paths($study_id, $runs_path, $config_or_path);
}

sub to_folder {
   my ($self, $path) = @_;
   make_path $path;
   my $tsv_path = sprintf("%s/%s.tsv", $path, $self->{study_id});
   my $yaml_path = sprintf("%s/%s.yaml", $path, $self->{study_id});
   write_file($tsv_path, join("\n", "Run", @{$self->{runs}})); 
   if(%{$self->{config}}){
     DumpFile($yaml_path,$self->{config});
   } else {
     unlink $yaml_path if -f $yaml_path;
   }
}
1;
