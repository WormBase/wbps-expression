use strict;
use warnings;
package WbpsExpression::Analysis;
use File::Slurp qw/read_dir write_file/;
use File::Basename;
use WbpsExpression::Analysis::DataFiles;
use WbpsExpression::Analysis::DESeq2;
use WbpsExpression::Model::Design;
use List::Util qw/pairmap pairs pairgrep uniq/;
use File::Path qw/make_path/;
# use Smart::Comments '###';

sub low_qc_by_condition {
   my ($runs_by_condition, $qc_issues_per_run) = @_;
   my %result = ();
   my %d = %{$runs_by_condition};
    for my $c (keys %d){
      my @runs = @{$d{$c}};
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
sub low_replicate_by_condition {
  my ($runs_by_condition) = @_;
  my %result = ();
  my %d = %{$runs_by_condition};
  for my $c (keys %d){
    my $number_runs = @{$d{$c}};
    $result{$c} = $number_runs if $number_runs < 3;
  }
  return \%result;
}

my %ANALYSES = (  
  study_design => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    $study->{design}->to_tsv($output_path);
  },
  aggregate_by_run => sub { 
    my($study, $files, $output_path, %analysis_args) = @_; 
    my @runs = $study->{design}->all_runs;
    my %qc_issues_per_run = pairgrep {$b} map {$_ => $files->{$_}{qc_issues}} @runs;
    my @qc_warnings = map {
       my $qcs = $qc_issues_per_run{$_};
       $qcs ? "!$_: ".join(". ", sort map {ucfirst $_} @{$qcs}) : ()
     } @runs;
    my @frontmatter = $analysis_args{decorate_files} ? ($analysis_args{description}, @qc_warnings) : ();
    my @name_to_path_pairs = map {
      my $name = $analysis_args{decorate_files} && $qc_issues_per_run{$_} ? "!$_" : $_;
      my $path = $files->{$_}{$analysis_args{source}};
      [$name, $path]
    } @runs;
    WbpsExpression::Analysis::DataFiles::aggregate(\@name_to_path_pairs, $output_path, @frontmatter);
  },
  average_by_condition => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    my $runs_by_condition = $study->{design}->runs_by_condition;
    my $runs_by_condition_then_replicate = $study->{design}->runs_by_condition_then_replicate;
    my @conditions_ordered = $study->{design}->all_conditions; 
    my %qc_issues_per_run = map {$_ => $files->{$_}{qc_issues}} map {@{$_}} map {values %{$_}} values %{$runs_by_condition_then_replicate};
    my $low_qc_by_condition = low_qc_by_condition($runs_by_condition, \%qc_issues_per_run);
    my $low_replicate_by_condition = low_replicate_by_condition($runs_by_condition);
### $low_qc_by_condition
### $low_replicate_by_condition
    my @qc_warnings = map {
      my @low_qc = sort map {ucfirst $_} @{$low_qc_by_condition->{$_} // []};
      my $low_rep = $low_replicate_by_condition->{$_} ? sprintf ("Low replicates (%s)", $low_replicate_by_condition->{$_}) : "";
      @low_qc || $low_rep ? "!$_: ".join(". ", @low_qc, $low_rep || ()) : ()
    } @conditions_ordered; 
    my @frontmatter = $analysis_args{decorate_files} ? ($analysis_args{description}, @qc_warnings) : ();
### @frontmatter
    my @name_to_pathlist_pairs= map {
       my $name = $_;
       
       if($analysis_args{decorate_files} && ($low_qc_by_condition->{$_} || $low_replicate_by_condition->{$_})){
          $name = "!$name";
       }
       
       my @paths = map {[map {$files->{$_}{$analysis_args{source}}} @{$_}]} values %{$runs_by_condition_then_replicate->{$_}};
       [$name, \@paths]
    } @conditions_ordered;
### @name_to_pathlist_pairs 
    WbpsExpression::Analysis::DataFiles::average_and_aggregate(\@name_to_pathlist_pairs, $output_path, @frontmatter);
  },
  differential_expression => sub {
    my($study, $files, $output_path, %analysis_args) = @_;
### %analysis_args
    my %runs_by_condition_all = %{$study->{design}->runs_by_condition};
    my @conditions_ordered = uniq map {($_->[0], $_->[1])} @{$analysis_args{contrasts}};
    my %runs_by_condition;
    @runs_by_condition{@conditions_ordered} = @runs_by_condition_all{@conditions_ordered};
    my %qc_issues_per_run = map {$_ => $files->{$_}{qc_issues}} map {@{$_}} values %runs_by_condition;
    my $low_qc_by_condition = low_qc_by_condition(\%runs_by_condition, \%qc_issues_per_run);
    my $low_replicate_by_condition = low_replicate_by_condition(\%runs_by_condition);
    my @qc_warnings;
    my @contrasts_amended_names;
    my $contrasts_low_replicates;
    for (@{$analysis_args{contrasts}}){
       my ($reference, $test, $name) = @{$_};
       my @qcs  = ((map {"$reference - $_"} @{$low_qc_by_condition->{$reference} //[]}),( map {"$test - $_"} @{$low_qc_by_condition->{$test} //[]}));
       my $low_replicates = $low_replicate_by_condition->{$reference} || $low_replicate_by_condition->{$test};
       push @qc_warnings, "!$name: ".join(". ", sort map {ucfirst $_} @qcs) if @qcs;
       $contrasts_low_replicates++ if $low_replicates;
       $name = "!$name" if @qcs || $low_replicates;
       push @contrasts_amended_names, [$reference, $test, $name ];
    }  
    # We don't allow contrasts with 1 replicate. So if low, then 2.
    if(keys %{$low_replicate_by_condition} == @conditions_ordered){
       push @qc_warnings, "! Contrasts based on conditions with low (2) replicates";
    } elsif (%{$low_replicate_by_condition}){
       push @qc_warnings, sprintf("! Conditions with low (2) replicates used in %s/%s contrasts:",  $contrasts_low_replicates, scalar @contrasts_amended_names)
         if $contrasts_low_replicates; 
       push @qc_warnings, "! - $_" for keys %{$low_replicate_by_condition} ;
    }

    my $source_file = join("/", dirname ($output_path), $analysis_args{source_file_name});
    die "Does not exist: $source_file" unless -f $source_file;
    my $design_dump_file = join("/", dirname($output_path), $study->{study_id}.".design.tsv.tmp");
    $study->{design}->to_tsv($design_dump_file);
    my @frontmatter = ($analysis_args{description}, @qc_warnings);
    eval {
      WbpsExpression::Analysis::DESeq2::do_analysis(
        $design_dump_file, $source_file, \@contrasts_amended_names, $output_path, @frontmatter
      );
    };
    unlink $design_dump_file;
    die $@ if $@;
  },
);
sub run {
  my ($output_dir, $study, $files) = @_;
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
sub run_all {
  my ($studies, $data_files, $output_dir) = @_;
  make_path $output_dir;
  for my $study (@{$studies}){
     print STDERR sprintf("Running: %s\n", $study->{study_id}) if $ENV{ANALYSIS_VERBOSE};
    run($output_dir, $study, $data_files->{$study->{study_id}}) unless $ENV{ANALYSIS_SKIP_ALL};
  }
}
1;
