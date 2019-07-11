use strict;
use warnings;
package WbpsExpression::Analysis;
use File::Slurp qw/read_dir write_file/;
use File::Basename;
use WbpsExpression::Analysis::QualityWarnings;
use WbpsExpression::Analysis::DataFiles;
use WbpsExpression::Analysis::DESeq2;
use WbpsExpression::Study::Design;
use List::Util qw/pairgrep pairs/;
use File::Path qw/make_path remove_tree/;
use Log::Any '$log';
use Archive::Zip;

# use Smart::Comments '###';

my %ANALYSES = (  
  study_design => sub {
    my($study, $output_path) = @_; 
    $study->{design}->to_tsv($output_path);
  },
  counts_per_run => sub {
    my($study, $output_path) = @_;
    my @runs = $study->{design}->all_runs;
    my @name_to_path_pairs = pairs map {
      $_ => $study->source_counts($_),
    } @runs;
    WbpsExpression::Analysis::DataFiles::aggregate(\@name_to_path_pairs, $output_path);
  },
  tpms_per_run => sub {
    my($study, $output_path, $description) = @_;
    my @runs = $study->{design}->all_runs;
    my %qc_issues_per_run = %{$study->qc_issues_per_run};

    my @warnings = map {
       my $qcs = $qc_issues_per_run{$_};
       $qcs ? "!$_: ".join(". ", sort map {ucfirst $_} @{$qcs}) : ()
     } @runs;

    my @frontmatter = (split("\n", $description), @warnings);
    
    my @name_to_path_pairs = map {
      my $name = $qc_issues_per_run{$_} ? "!$_" : $_;
      my $path = $study->source_tpm($_);
      [$name, $path]
    } @runs;
    WbpsExpression::Analysis::DataFiles::aggregate(\@name_to_path_pairs, $output_path, @frontmatter);
  },
  average_by_condition => sub {
    my($study, $output_path, $description) = @_; 
    my @conditions_ordered = $study->{design}->all_conditions;
    my %qc_issues_per_run = %{$study->qc_issues_per_run};
    my ($conditions_amended_names, $warnings) = WbpsExpression::Analysis::QualityWarnings::conditions_amended_names_and_warnings_for_per_condition_analysis(
      $study->{design},
      \%qc_issues_per_run,
      \@conditions_ordered,
    );

    my @frontmatter = (split("\n", $description), @{$warnings});
    my @name_to_pathlist_pairs = map {
       my $condition = $_;
       my $runs_by_replicate = $study->{design}->runs_by_condition_then_replicate->{$condition};
       my @paths = map {[map {$study->source_tpm($_)} @{$_}]} values %{$runs_by_replicate};
       [$conditions_amended_names->{$condition}, \@paths]
    } @conditions_ordered;

    WbpsExpression::Analysis::DataFiles::average_and_aggregate(\@name_to_pathlist_pairs, $output_path, @frontmatter);
  },
  differential_expression => sub {
    my($study, $output_paths, $descriptions, %analysis_args) = @_;
    my %contrasts = %{$analysis_args{contrasts}};

    my %frontmatters;
    my %warnings_by_contrast_name_by_key;
    my $qc_issues_per_run = $study->qc_issues_per_run;
    for my $key (keys %contrasts){
      my ($contrasts_with_amended_names, $warnings_per_contrast_name, $warnings_all_contrasts) = WbpsExpression::Analysis::QualityWarnings::amended_contrasts_and_warnings_for_per_contrast_analysis(
        $study->{design},
        $qc_issues_per_run,
        $contrasts{$key},
      );
      $contrasts{$key} = $contrasts_with_amended_names;
      $frontmatters{$key} = [split("\n", $descriptions->{$key}), @{$warnings_all_contrasts}];
      $warnings_by_contrast_name_by_key{$key} = $warnings_per_contrast_name;
    }

    (my $unfiltered_results_folder = $output_paths->{all_contrasts}) =~ s/\.zip$//;

    my $source_file = join("/", dirname ($unfiltered_results_folder), $analysis_args{source_file_name});
    die "Does not exist: $source_file" unless -f $source_file;
    my $design_dump_file = join("/", dirname($unfiltered_results_folder), $study->{study_id}.".design.tsv.tmp");
    $study->{design}->to_tsv($design_dump_file);


    mkdir $unfiltered_results_folder;
    eval {
      WbpsExpression::Analysis::DESeq2::do_analysis(
        $design_dump_file, $source_file, $unfiltered_results_folder, 
        \%contrasts, $output_paths, \%frontmatters, \%warnings_by_contrast_name_by_key,
      );
    };
    unlink $design_dump_file;
    die $@ if $@;
    my $zip = Archive::Zip->new();
    $zip->addTree($unfiltered_results_folder, basename($unfiltered_results_folder));
    $zip->writeToFileNamed($output_paths->{all_contrasts});
    remove_tree $unfiltered_results_folder;
  },
);
sub run {
  my ($output_dir, $study) = @_;
  my @analyses = $study->analyses_required;
  make_path join("/", $output_dir, $study->{study_id});
  my %done = map { $_=>1 } read_dir join("/", $output_dir, $study->{study_id}); 
  my @analyses_to_do = grep {
    my $file_name = $_->{file_name} // $_->{files}->[-1]->{file_name};
    not $done{$file_name}
  } @analyses;

  for my $analysis (@analyses_to_do){
    if($analysis->{files}){
       my %file_paths = map {$_->{key} => join("/", $output_dir, $study->{study_id}, $_->{file_name})} @{$analysis->{files}};
       my %descriptions = map {$_->{description} ? ( $_->{key} => $_->{description}) : ()} @{$analysis->{files}};

       &{$ANALYSES{$analysis->{type}}}($study, \%file_paths, \%descriptions, %{$analysis});
    } else {
        my $file_path = join("/", $output_dir, $study->{study_id}, $analysis->{file_name});
        my $description = $analysis->{description};

       &{$ANALYSES{$analysis->{type}}}($study, $file_path, $description, %{$analysis});
   }
  }
  return @analyses_to_do;
}
sub run_all {
  my ($studies, $output_dir) = @_;
  make_path $output_dir;
  $log->info(sprintf(__PACKAGE__ . " run_all %s studies, output_dir %s", scalar @{$studies}, $output_dir));
  for my $study (@{$studies}){
    run($output_dir, $study) unless $ENV{ANALYSIS_SKIP_ALL};
  }
}
1;
