use strict;
use warnings;
package Production::Analysis;
use File::Slurp qw/read_dir/;
use Production::Analysis::DataFiles;
use Model::Design;
use List::Util qw/pairmap pairgrep/;
use File::Path qw/make_path/;
#use Smart::Comments '###';
sub new {
  my ($class, $dir) = @_;
  return bless {
    dir => $dir,
  }, $class;
}
sub low_qc_by_condition {
   my ($runs_by_condition, $qc_issues_per_run) = @_;
   my %result = ();
   my %d = %{$runs_by_condition};
    for my $c (keys %d){
      my @runs = @{$d{$c}};
      push @{$result{$c}}, sprintf("Low replicates (%s)", scalar @runs) if @runs < 3;
      my %qcs;
      for my $run_id (@runs){
         for my $qc (@{$qc_issues_per_run->{$run_id}}){
            push @{$qcs{$qc}}, $run_id;
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
    my @runs = $study->{design}->all_runs;
    my %qc_issues_per_run = pairgrep {$b} map {$_ => $files->{$_}{qc_issues}} @runs;
    my @qc_warnings = map {
       my $qcs = $qc_issues_per_run{$_};
       $qcs ? "!$_: ".join(". ", sort map ucfirst @{$qcs}) : ()
     } @runs;
    my @frontmatter = $analysis_args{decorate_files} ? ($analysis_args{description}, @qc_warnings) : ();
    my @name_to_path_pairs = map {
      my $name = $analysis_args{decorate_files} && $qc_issues_per_run{$_} ? "!$_" : $_;
      my $path = $files->{$_}{$analysis_args{source}};
      [$name, $path]
    } @runs;
    Production::Analysis::DataFiles::aggregate(\@name_to_path_pairs, $output_path, @frontmatter);
  },
  average_by_condition => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    my $runs_by_condition = $study->{design}->runs_by_condition;
    my @conditions_ordered = $study->{design}->all_conditions; 
    my %qc_issues_per_run = map {$_ => $files->{$_}{qc_issues}} map {@{$_}} values %{$runs_by_condition};
    my $low_qc_by_condition = low_qc_by_condition($runs_by_condition, \%qc_issues_per_run);
#### $low_qc_by_condition
    my @qc_warnings = map {
      $low_qc_by_condition->{$_} ? ("!$_: ".join(". ", sort map {ucfirst $_} @{$low_qc_by_condition->{$_}})): ()
    } @conditions_ordered;
    my @frontmatter = $analysis_args{decorate_files} ? ($analysis_args{description}, @qc_warnings) : ();
#### @frontmatter
    my @name_to_pathlist_pairs= map {
       my $name = $analysis_args{decorate_files} && $low_qc_by_condition->{$_} ? "!$_" : $_;
       my @paths = map {$files->{$_}{$analysis_args{source}}} @{$runs_by_condition->{$_}};
       [$name, \@paths]
    } @conditions_ordered;
    Production::Analysis::DataFiles::average_and_aggregate(\@name_to_pathlist_pairs, $output_path, @frontmatter);
  }
);
sub run {
  my ($self, $output_dir, $study, $files) = @_;
  my @analyses = $study->analyses_required;
  make_path join("/", $output_dir, $study->{study_id});
  my %done = map { $_=>1 } read_dir join("/", $output_dir, $study->{study_id}); 
  my @analyses_to_do = grep {not $done{$_->{file_name}}} @analyses;
  for my $analysis (@analyses_to_do){
    &{$ANALYSES{$analysis->{type}}}(
       $study,
       $files,
       join("/", $output_dir, $study->{study_id}, $analysis->{file_name}),
       %{$analysis},
    );
  }
  return @analyses_to_do;
}
sub run_all_and_produce_markdown_report {
  my ($self, %args) = @_;
  my $output_dir = join("/", $self->{dir}, $args{species} , $args{assembly});
  make_path $output_dir;
  open(my $fh, ">", "$output_dir/index.md") or die "$output_dir/index.md: $!";
  print $fh sprintf("# %s - public RNASeq studies\nAssembly: %s\n", do {
    my $species = $args{species};
    $species =~ s/_/ /g;
    ucfirst($species)
  }, $args{assembly});
  print $fh "## Analysed\n";
  for my $study (@{$args{studies}{todo}}){
     $self->run($output_dir, $study, $args{files}{$study->{study_id}});
     print $fh $study->to_markdown;
  }
  print $fh "## Failing curation checks\n";
  for my $study (@{$args{studies}{failed_checks}}){
     print $fh $study->to_markdown;
  }
  print $fh "## Skipped\n";
  print $fh sprintf("##### Total: %s\n", scalar @{$args{studies}{ids_skipped}});
  print $fh join(", ", @{$args{studies}{ids_skipped}})."\n";

  close $fh;
}
1;
