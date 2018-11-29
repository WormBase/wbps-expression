use strict;
use warnings;
package Production::Analysis;
use File::Slurp qw/read_dir/;
use Production::Analysis::DataFiles qw/aggregate average_and_aggregate/;
use Model::Design;
use List::Util qw/pairmap/;
use File::Path qw/make_path/;
#use Smart::Comments;
sub new {
  my ($class, $dir) = @_;
  return bless {
    dir => $dir,
  }, $class;
}
sub low_qc_by_condition {
   my ($runs_by_condition, $qc_issues_per_run) = @_;
   my %result;
   my %d = %{$runs_by_condition};
    for my $c (keys %d){
      my @runs = @{$d{$c}};
      push @{$result{$c}}, sprintf("Low replicates (%s)", scalar @runs) if @runs < 3;
      my %qcs;
      for my $run (@runs){
         for my $qc (@{$qc_issues_per_run->{$run_id}}){
            push @{$qcs{$qc}}, $run;
         }
      }
      for my $qc (keys %qcs){
         push @{$result{$c}}, sprintf("%s: %s %s", $qc, (@{$qcs{$qc}} > 1 ? "runs" : "run"), join(", ", @{$qcs{$qc}}));
      }
   }
   return \%result;
}
my %ANALYSES = (
  aggregate_by_run => sub { 
    my($study, $files, $output_path, %analysis_args) = @_; 
    my @name_to_path_pairs = map {[$_, $files->{$_}{$analysis_args{source}}]} $study->{design}->all_runs;
    aggregate(\@name_to_path_pairs, $output_path, $analysis_args{description});
  },
  average_by_condition => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    my @frontmatter;
    my $runs_by_condition = $study->{design}->runs_by_condition;
    my %qc_issues_per_run = map {$_ => $files->{$_}{qc_issues}} map {@{$_}} values %{$runs_by_condition};
    my $low_qc_by_condition = low_qc_by_condition($runs_by_condition, \%qc_issues_per_run);
    my @frontmatter = sort pairmap {"!$a: ".join (". ", sort map ucfirst @{$b}) } %{$low_qc_by_condition};
    my @name_to_pathlist_pairs= pairmap {
       my $name = $low_qc_by_condition->{$a} ? "!$a" : $a;
       my @paths = map {$files->{$_}{$analysis_args{source}}} @{$b};
       [$name, \@paths]
    } %{$runs_by_condition};
    average_and_aggregate(\@name_to_pathlist_pairs, $output_path, $analysis_args{description}, @frontmatter);
  }
);
sub run {
  my ($self, $study, $files) = @_;
  my @analyses = $study->analyses_required;
  my $output_dir = join("/", $self->{dir}, $study->{study_id});
  make_path $output_dir;
  my %done = map { $_=>1 } read_dir $output_dir; 
  my @analyses_to_do = grep {not $done{$_->{file_name}}} @analyses;
  for my $analysis (@analyses_to_do){
    &{$ANALYSES{$analysis->{type}}}(
       $study,
       $files,
       join("/", $output_dir, $analysis->{file_name}),
       %{$analysis},
    );
  }
  return @analyses_to_do;
}
1;
