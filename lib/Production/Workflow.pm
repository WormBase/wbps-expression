use strict;
use warnings;
package Production::Workflow;
use Production::Sheets;
use Production::CurationDefaults;
use Production::Analysis;
use PublicResources::Rnaseq;
use File::Basename qw/dirname/;
use File::Slurp qw/read_dir/;
use File::Path qw/make_path/;
use List::Util qw/first/;
use List::MoreUtils qw/uniq/;
use Model::Study;
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

sub should_reject_study {
  my ($self, $study) = @_;
  return $study->{design}->all_runs < 6;
}

sub fetch_incoming_studies {
  my ($self, $public_study_records, $current_studies, $current_ignore_studies) = @_;

  my %result = (REJECT => [], SAVE => []); 
  for my $study ( map {&Production::CurationDefaults::study(%$_)} @{$public_study_records}){
    my $current_record = $current_studies->{$study->{study_id}};
    if ($current_ignore_studies->{$study->{study_id}} or $self->should_reject_study($study)){
      push @{$result{REJECT}}, $study;
    } else { 
      if ($current_record and Model::Study::config_matches_design_checks($current_record->{config}, $study->{design})){
        $study->{config}{slices} = $current_record->{config}{slices};
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
  my %result = (FAILED_CHECKS => [], PASSED_CHECKS => []);
  for my $study (@studies){
    push @{$result{$study->passes_checks ? "PASSED_CHECKS" : "FAILED_CHECKS"}},$study;
  }
  return \%result;
}
sub do_everything {
  my ($self, $species, $assembly) = @_;
  my @public_study_records = $self->{public_rnaseq_studies}->get($species, $assembly);  
  my %current_studies = map {$_->{study_id}=> $_} $self->get_studies_in_sheets($species);
  my %current_ignore_studies = map {$_=>1} $self->{sheets}->list("ignore_studies", $species);
  my $incoming_studies = $self->fetch_incoming_studies(\@public_study_records, \%current_studies, \%current_ignore_studies);
  for my $study (@{$incoming_studies->{SAVE}}){
     $study->to_folder($self->{sheets}->path("studies", $species, $study->{study_id}));
  }
  my @new_ignore_studies = uniq sort(keys %current_ignore_studies, map {$_->{study_id}} @{$incoming_studies->{REJECT}});
  if (@{$incoming_studies->{REJECT}}){
    $self->{sheets}->write_list( \@new_ignore_studies, "ignore_studies", $species);
  }
  
  my $todo_studies = $self->run_checks(@{$incoming_studies->{SAVE}});

  my %files;
  for my $study (@{$todo_studies->{PASSED_CHECKS}}){
     my $public_study_record = first {$_->{study_id} eq $study->{study_id}} @public_study_records; 
     for my $run (@{$public_study_record->{runs}}){
        $files{$study->{study_id}}{$run->{run_id}} = { 
          %{$run->{data_files}},
          qc_issues => $_->{qc_issues},
        };
     }
  }
  $self->{analysis}->run_all_and_produce_markdown_report(
    species => $species,
    assembly => $assembly,
    studies => {
      ids_skipped => \@new_ignore_studies,
      failed_checks => $todo_studies->{FAILED_CHECKS},
      todo => $todo_studies->{PASSED_CHECKS},
    },
    files => \%files,
  );
  # Future directions:
  # - Deployment directory - link between production directory results, a corner of FTP where they serve, and where the paths should go
  # - Instead of the markdown report, make something deployable
  # - Report on what just happened: I'm not sure actually, maybe it's better to always show state that was achieved after the run
}
1;
