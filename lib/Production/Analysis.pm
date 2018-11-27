use strict;
use warnings;
use File::Slurp qw/read_dir/;
use Model::DataDir;
package Production::Analysis;
use Production::Analysis::DataFiles qw/aggregate average_and_aggregate/;
use Model::Design;
use List::Util qw/pairmap/;
#use Smart::Comments;
sub new {
  my ($class, $dir) = @_;
  return bless {
    dir => $dir,
  }, $class;
}

my %ANALYSES = (
  aggregate_by_run => sub { 
    my($study, $files, $output_path, %analysis_args) = @_; 
    my @name_to_path_pairs = map {[$_, $study->{files}{$_}{$analysis_args{source}}]} @{$study->{design}->all_runs};
    aggregate(\@name_to_path_pairs, $output_path, $analysis_args{description});
  },
  average_by_condition => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    my @name_to_pathlist_pairs = pairmap {[$a, [map {$study->{files}{$_}{$analysis_args{source}}} @{$b}]]} %{$study->{design}->runs_by_condition};  
    average_and_aggregate(\@name_to_pathlist_pairs, $output_path, $analysis_args{description});
  }
);
sub run {
  my ($self, $study, $files) = @_;
  my @analyses = $study->analyses_required;
  my $output_dir = Model::DataDir->new(join("/", $self->{dir}, $study->{study_id}));
  my %done = $output_dir->list;
  my @todo = grep {not $done{$_->{file_name}}} @analyses;
  for my $analysis (@todo){
    &{$ANALYSES{$analysis->{type}}}(
       $study,
       $files,
       $output_dir->path($analysis->{file_name}),
       %{$analysis},
    );
  }
  return @todo;
}

sub do_all {
  my ($self, $location_per_run_id, @studies) = @_;
  my %result;
  for my $study (@studies){
    my @analyses = $study->analyses_required;
    my $output_dir = Model::DataDir->new(join("/", $self->{dir}, $study->{study_id}));
    my %done = $output_dir->list;
# grep {not $output_dir->done($analysis)} using type and hm, name?
    for my $analysis (grep {not $done{$_->{type}}} @analyses){
# Hm, maybe not? How do I cope with the fact there will be more than one analysis per type?
# Reminds me of the system of URIs for Atlas resources I wrote
       my $out_path = $output_dir->path($analysis);
       &{$ANALYSES{$analysis->{type}}}($stu
       
    }
  }
}
1;
