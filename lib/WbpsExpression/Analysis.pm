use strict;
use warnings;
package WbpsExpression::Analysis;
use File::Slurp qw/read_dir write_file/;
use File::Basename;
use WbpsExpression::Analysis::QualityWarnings;
use WbpsExpression::Analysis::DataFiles;
use WbpsExpression::Analysis::DESeq2;
use WbpsExpression::Model::Design;
use List::Util qw/pairgrep/;
use File::Path qw/make_path/;
use Log::Any '$log';
# use Smart::Comments '###';

my %ANALYSES = (  
  study_design => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    $study->{design}->to_tsv($output_path);
  },
  aggregate_by_run => sub { 
    my($study, $files, $output_path, %analysis_args) = @_; 
    my $do_decorate_files = $analysis_args{description} ? 1 : 0;
    my @runs = $study->{design}->all_runs;
    my %qc_issues_per_run = pairgrep {$b} map {$_ => $files->{$_}{qc_issues}} @runs;
    my @qc_warnings = map {
       my $qcs = $qc_issues_per_run{$_};
       $qcs ? "!$_: ".join(". ", sort map {ucfirst $_} @{$qcs}) : ()
     } @runs;
    my @frontmatter = $do_decorate_files ? (split("\n", $analysis_args{description}), @qc_warnings) : ();
    my @name_to_path_pairs = map {
      my $name = $do_decorate_files && $qc_issues_per_run{$_} ? "!$_" : $_;
      my $path = $files->{$_}{$analysis_args{source}};
      [$name, $path]
    } @runs;
    WbpsExpression::Analysis::DataFiles::aggregate(\@name_to_path_pairs, $output_path, @frontmatter);
  },
  average_by_condition => sub {
    my($study, $files, $output_path, %analysis_args) = @_; 
    my @conditions_ordered = $study->{design}->all_conditions;
    my %qc_issues_by_run = map {$_ => $files->{$_}{qc_issues}} $study->{design}->all_runs;
    my ($conditions_amended_names, $warnings) = WbpsExpression::Analysis::QualityWarnings::conditions_amended_names_and_warnings_for_per_condition_analysis(
      $study->{design},
      \%qc_issues_by_run,
      \@conditions_ordered,
    );

    my @frontmatter = (split("\n", $analysis_args{description}), @{$warnings});
    my @name_to_pathlist_pairs = map {
       my $condition = $_;
       my $runs_by_replicate = $study->{design}->runs_by_condition_then_replicate->{$condition};
       my @paths = map {[map {$files->{$_}{$analysis_args{source}}} @{$_}]} values %{$runs_by_replicate};
       [$conditions_amended_names->{$condition}, \@paths]
    } @conditions_ordered;

    WbpsExpression::Analysis::DataFiles::average_and_aggregate(\@name_to_pathlist_pairs, $output_path, @frontmatter);
  },
  differential_expression => sub {
    my($study, $files, $output_path, %analysis_args) = @_;
    my %qc_issues_by_run = map {$_ => $files->{$_}{qc_issues}} $study->{design}->all_runs;
    my ($contrasts_amended_names, $warnings) = WbpsExpression::Analysis::QualityWarnings::amended_contrasts_and_warnings_for_per_contrast_analysis(
      $study->{design},
      \%qc_issues_by_run,
      $analysis_args{contrasts},
    );
    my @frontmatter = (split("\n", $analysis_args{description}), @{$warnings});
    my $source_file = join("/", dirname ($output_path), $analysis_args{source_file_name});
    die "Does not exist: $source_file" unless -f $source_file;
    my $design_dump_file = join("/", dirname($output_path), $study->{study_id}.".design.tsv.tmp");
    $study->{design}->to_tsv($design_dump_file);

    eval {
      WbpsExpression::Analysis::DESeq2::do_analysis(
        $design_dump_file, $source_file, $contrasts_amended_names, $output_path, @frontmatter,
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
  $log->info(sprintf(__PACKAGE__ . " run_all %s studies, output_dir %s", scalar @{$studies}, $output_dir));
  for my $study (@{$studies}){
    run($output_dir, $study, $data_files->{$study->{study_id}}) unless $ENV{ANALYSIS_SKIP_ALL};
  }
}
1;
