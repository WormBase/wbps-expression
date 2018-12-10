use strict;
use warnings;
package Production::Analysis;
use File::Slurp qw/read_dir write_file/;
use Text::MultiMarkdown qw/markdown/;
use File::Basename;
use Production::Analysis::DataFiles;
use Production::Analysis::DESeq2;
use Model::Design;
use List::Util qw/pairmap pairgrep/;
use File::Path qw/make_path/;
use View::StudiesPage;
# use Smart::Comments '###';
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
      push @{$result{$c}}, sprintf("low replicates (%s)", scalar @runs) if @runs < 3;
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
       $qcs ? "!$_: ".join(". ", sort map {ucfirst $_} @{$qcs}) : ()
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
  },
  differential_expression => sub {
    my($study, $files, $output_path, %analysis_args) = @_;
### %analysis_args
    my $runs_by_condition = $study->{design}->runs_by_condition;
    my @conditions_ordered = $study->{design}->all_conditions; 
    my %qc_issues_per_run = map {$_ => $files->{$_}{qc_issues}} map {@{$_}} values %{$runs_by_condition};
    my $low_qc_by_condition = low_qc_by_condition($runs_by_condition, \%qc_issues_per_run);
    my @qc_warnings;
    my @contrasts_amended_names;
    for (@{$analysis_args{contrasts}}){
       my ($reference, $test, $name) = @{$_};
       my @qcs  = ((map {"$reference - $_"} @{$low_qc_by_condition->{$reference} //[]}),( map {"$test - $_"} @{$low_qc_by_condition->{$test} //[]}));
       if(@qcs){
          push @qc_warnings, "!$name: ".join(". ", sort map {ucfirst $_} @qcs);
          $name = "!$name";
       }
       push @contrasts_amended_names, [$reference, $test, $name ];
    }

    my $source_file = join("/", dirname ($output_path), $analysis_args{source_file_name});
    die "Does not exist: $source_file" unless -f $source_file;
    my $design_dump_file = join("/", dirname($output_path), $study->{study_id}.".design.tsv.tmp");
    $study->{design}->to_tsv($design_dump_file);
    my @frontmatter = ($analysis_args{description}, @qc_warnings);
    eval {
      Production::Analysis::DESeq2::do_analysis(
        $design_dump_file, $source_file, \@contrasts_amended_names, $output_path, @frontmatter
      );
    };
    unlink $design_dump_file;
    die $@ if $@;
  },
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
sub run_all {
  my ($self, %args) = @_;
  my $output_dir = join("/", $self->{dir}, $args{species} , $args{assembly});
  make_path $output_dir;
  for my $study (@{$args{studies}}){
     print STDERR sprintf("Running: %s\n", $study->{study_id}) if $ENV{ANALYSIS_VERBOSE};
     $self->run($output_dir, $study, $args{files}{$study->{study_id}});
  }
  print STDERR "Writing page: $output_dir/index.html\n" if $ENV{ANALYSIS_VERBOSE};
  write_file("$output_dir/index.html", View::StudiesPage->new($args{species}, @{$args{studies}})->to_html);
}
1;