use strict;
use warnings;
package WbpsExpression::IncomingStudies;
use WbpsExpression::IncomingStudies::Sheets;
use WbpsExpression::IncomingStudies::CurationDefaults;
use WbpsExpression::Model::Study;
use WbpsExpression::Model::SkippedRuns;
use PublicResources::Rnaseq;
use File::Basename qw/dirname/;
use File::Slurp qw/read_dir write_file/;
use File::Path qw/make_path/;
use List::Util qw/first/;
use List::MoreUtils qw/uniq/;
# use Smart::Comments '###';
sub new {
  my ($class, $sheets, $public_rnaseq_studies) = @_;
  return bless {
     sheets => $sheets,
     public_rnaseq_studies => $public_rnaseq_studies,
  }, $class;
}

#Add to me for unit testing
our %exceptions = (
  "ERP006987" => "C. sinensis: four runs, but the tpms could be meaningful",
  "SRP013211" => "O. viverrini: two runs, juvenile vs adult, from 2012",
  "DRP003063" => "Emu study, four runs",
  "SRP131874" => "E. granulosus, four runs but a second study for the species",
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

sub fetch_studies {
  my ($self, $public_study_records, $current_studies, $current_skipped_runs_per_study_id) = @_;
  my @studies;
  my @skipped_runs;
  for my $study (@{$public_study_records}){
#### fetch_studies: $study
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
      my $config = @remaining_runs ? {} : &WbpsExpression::IncomingStudies::CurationDefaults::config_base(%$study);
      push @skipped_runs, WbpsExpression::Model::SkippedRuns->new($study->{study_id}, $config, [keys %runs_to_skip]);
    }
    if (@remaining_runs){
      $study->{runs} = \@remaining_runs; 
      $study = &WbpsExpression::IncomingStudies::CurationDefaults::study(%$study);
      if ($current_record and WbpsExpression::Model::Study::config_matches_design_checks($current_record->{config}, $study->{design}) and not $ENV{RECREATE_ALL_CONFIGS}){
        $study->{config}{contrasts} = $current_record->{config}{contrasts};
        $study->{config}{condition_names} = $current_record->{config}{condition_names};
        # Additionally, characteristics in the current record were already reused
        # because they provided sources of attributes for the runs - see PublicResources::Rnaseq
      }
      push @studies, $study;
    }
  };
  return \@studies, \@skipped_runs;
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

sub fetch_all {
  my ($self, $species, $assembly) = @_;
  unlink $self->{sheets}->path("ignore_studies", "$species.tsv");
  my @public_study_records = sort {$a->{study_id} cmp $b->{study_id} } $self->{public_rnaseq_studies}->get($species, $assembly);  
  my %data_files;
  for my $public_study_record (@public_study_records){
     for my $run (@{$public_study_record->{runs}}){
        $data_files{$public_study_record->{study_id}}{$run->{run_id}} = { 
          %{$run->{data_files}},
          qc_issues => $_->{qc_issues},
        };
     }
  }
  my ($our_studies, $skipped_runs_in_our_studies, $other_studies) = $self->{sheets}->read_directories($species);
  my %current_studies = map {$_->{study_id}=> $_} @{$our_studies};
  my %current_skipped_runs_per_study_id  = map {$_->{study_id}=> $_} (@{$skipped_runs_in_our_studies}, @{$other_studies});

  my ($selected_studies, $skipped_runs_in_studies)  = $self->fetch_studies(\@public_study_records, \%current_studies, \%current_skipped_runs_per_study_id);

  for my $study (@{$selected_studies}){
     $study->to_folder($self->{sheets}->path("studies", $species, $study->{study_id}));
  }
  for my $skipped_runs_in_study (@{$skipped_runs_in_studies}){
     $skipped_runs_in_study->to_folder($self->{sheets}->path("skipped_runs", $species, $skipped_runs_in_study->{study_id}));
  }

  my @other_studies = grep {
     my $skipped_study_id = $_->{study_id};
     my $found = grep {$skipped_study_id eq $_->{study_id}} @{$selected_studies};
     not $found
  } @{$skipped_runs_in_studies};

  return ($selected_studies, \@other_studies, \%data_files);
}
1;
