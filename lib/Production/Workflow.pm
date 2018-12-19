use strict;
use warnings;
package Production::Workflow;
use Production::Sheets;
use Production::CurationDefaults;
use Production::Analysis;
use PublicResources::Rnaseq;
use File::Basename qw/dirname/;
use File::Slurp qw/read_dir write_file/;
use File::Path qw/make_path/;
use List::Util qw/first/;
use List::MoreUtils qw/uniq/;
use Text::MultiMarkdown qw/markdown/;
use Model::Study;
use Model::SkippedRuns;
use View::StudiesPage;
#use Smart::Comments '###';
sub new {
  my ($class, $root_dir, $src_dir, $work_dir) = @_;
  my $sheets = Production::Sheets->new($src_dir);
  return bless {
     processing_path => "$root_dir/curation",
     sheets => $sheets,
     public_rnaseq_studies => PublicResources::Rnaseq->new($root_dir, $sheets),
     analysis => Production::Analysis->new($work_dir),
  }, $class;
}
sub get_studies_in_sheets {
  my ($self, $species) = @_;
  return map {Model::Study->from_folder($_)} $self->{sheets}->dir_content_paths("studies", $species);
}
sub get_skipped_runs_in_sheets {
  my ($self, $species) = @_;
  return map {Model::SkippedRuns->from_folder($_)} $self->{sheets}->dir_content_paths("skipped_runs", $species);
}

my %exceptions = (
  "ERP006987" => "C. sinensis: four runs, but the tpms could be meaningful",
  "SRP013211" => "O. viverrini: two runs, juvenile vs adult, from 2012",
);

sub should_reject_study {
  my ($self, $public_study_record) = @_;
  return not ($exceptions{$public_study_record->{study_id}}) && @{$public_study_record->{runs}//[]} < 6;
}

# I could be more complicated
# e.g. get the median number of characteristics in a study, reject if three times fewer?
sub should_reject_run {
  my ($study, $run) = @_;
  return not(grep {$_} values %{$run->{characteristics}});
}

sub fetch_incoming_studies {
  my ($self, $public_study_records, $current_studies, $current_skipped_runs_per_study_id) = @_;

  my %result = (NEW_SKIPPED_RUNS => [], SAVE => []); 
  for my $study (@{$public_study_records}){
    my $current_record = $current_studies->{$study->{study_id}};
    my $should_reject_all = $self->should_reject_study($study);
    my @should_reject_runs = map {
      $_->{run_id}
    } grep {
      $should_reject_all || should_reject_run($study, $_)
    } @{$study->{runs}};

    my @current_skipped_runs = @{ $current_skipped_runs_per_study_id->{$study->{study_id}}{runs}//[]};
    my %runs_to_skip = map {$_ => 1} (@should_reject_runs, @current_skipped_runs);
#### @current_skipped_runs
    if ($current_record and not $ENV{RECREATE_ALL_SKIPPED_RUNS}){
       delete $runs_to_skip{$_} for $current_record->{design}->all_runs;
    #### Remaining: sort keys %runs_to_skip
    }

    my @remaining_runs = grep {not $runs_to_skip{$_->{run_id}}} @{$study->{runs}};

    if (%runs_to_skip){
      my $config = @remaining_runs ? {} : &Production::CurationDefaults::config_base(%$study);
      push @{$result{NEW_SKIPPED_RUNS}}, Model::SkippedRuns->new($study->{study_id}, $config, [keys %runs_to_skip]);
    }
    if (@remaining_runs){
      $study->{runs} = \@remaining_runs; 
      $study = &Production::CurationDefaults::study(%$study);
      if ($current_record and Model::Study::config_matches_design_checks($current_record->{config}, $study->{design}) and not $ENV{RECREATE_ALL_CONFIGS}){
        $study->{config}{contrasts} = $current_record->{config}{contrasts};
        $study->{config}{condition_names} = $current_record->{config}{condition_names};
        # Additionally, characteristics in the current record were already reused
        # because they provided sources of attributes for the runs - see PublicResources::Rnaseq
      }
      push @{$result{SAVE}}, $study;
    }
  };
  return \%result;
}
sub run_checks {
  my ($self, @studies) = @_;
  my @failed_checks;
  my @passed_checks;
  for my $study (@studies){
    if($study->passes_checks){
       push @passed_checks, $study;
    } else {
       push @failed_checks, $study;
    }
  }
  return \@failed_checks, \@passed_checks;
}

sub do_everything {
  my ($self, $species, $assembly) = @_;
  unlink $self->{sheets}->path("ignore_studies", "$species.tsv");
  my @public_study_records = $self->{public_rnaseq_studies}->get($species, $assembly);  
  my %current_studies = map {$_->{study_id}=> $_} $self->get_studies_in_sheets($species);
  my %current_skipped_runs_per_study_id  = map {$_->{study_id}=> $_} $self->get_skipped_runs_in_sheets($species);

  my $incoming_studies = $self->fetch_incoming_studies(\@public_study_records, \%current_studies, \%current_skipped_runs_per_study_id);

  for my $study (@{$incoming_studies->{SAVE}}){
     $study->to_folder($self->{sheets}->path("studies", $species, $study->{study_id}));
  }
  for my $skipped_runs (@{$incoming_studies->{NEW_SKIPPED_RUNS}}){
     $skipped_runs->to_folder($self->{sheets}->path("skipped_runs", $species, $skipped_runs->{study_id}));
  }
  
  my ($studies_failing_checks, $studies_passing_checks) = $self->run_checks(@{$incoming_studies->{SAVE}});

  my %files;
  for my $study (@{$studies_passing_checks}){
     my $public_study_record = first {$_->{study_id} eq $study->{study_id}} @public_study_records; 
     for my $run (@{$public_study_record->{runs}}){
        $files{$study->{study_id}}{$run->{run_id}} = { 
          %{$run->{data_files}},
          qc_issues => $_->{qc_issues},
        };
     }
  }
  $self->{analysis}->run_all(
    species => $species,
    assembly => $assembly,
    studies => {
      passing_checks => $studies_passing_checks,
      failing_checks => $studies_failing_checks,
      skipped_runs => $incoming_studies->{NEW_SKIPPED_RUNS},
    },
    files => \%files,
  );
}
1;
