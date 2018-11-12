
package Production::Workflow;
use Production::Sheets;
use Production::CurationDefaults;
use PublicResources::Rnaseq;
use File::Basename qw/dirname/;
use File::Slurp qw/read_dir/;
use File::Path qw/make_path/;
use List::MoreUtils qw/uniq/;
use Model::Study;
sub new {
  my ($class, $root_dir, $src_dir) = @_;
  my $sheets = Production::Sheets->new($src_dir);
  return bless {
     processing_path => "$root_dir/curation",
     sheets => $sheets,
     public_rnaseq_studies => PublicResources::Rnaseq->new($root_dir, $sheets),
  }, $class;
}
sub get_studies_in_sheets {
  my ($self, $species) = @_;
  return map {Model::Study->from_folder($_)} $self->{sheets}->dir_content_paths("studies", $species);
}

sub should_reject_study {
  my ($self, $study) = @_;
  return ($study->{design}->all_conditions < 2 || $study->{design}->all_runs < 6);
}

sub update_sheets_with_incoming_studies {
  my ($self, $species, @new_studies) = @_;
  my @new_studies_saved;
  my @new_studies_rejected;
  my @ignore_studies = $self->{sheets}->list("ignore_studies", $species);
  for my $study (@new_studies) {
    my $study_id = $study->{study_id};
    if (grep {$_ eq $study_id} @ignore_studies or $self->should_reject_study($study)){ 
      push @new_studies_rejected, $study->{study_id};
    } else {
      $study->to_folder($self->{sheets}->path("studies", $species, $study->{study_id}));
      push @new_studies_saved, $study->{study_id};
    }
  }
  if(@new_studies_rejected){
    $self->{sheets}->write_list( [sort uniq(@ignore_studies, @new_studies_rejected)], "ignore_studies", $species)
  }
  return \@new_studies_rejected , \@new_studies_saved;
}


sub list_queue {
  my ($self, $species, $assembly) = @_;
  return read_dir join("/", $self->{processing_path}, $species, $assembly);
}

sub do_everything {
  my ($self, $species, $assembly) = @_;
  my ($factors, $location_per_run_id, @public_study_records) = $self->{public_rnaseq_studies}->get($species, $assembly);  
  my @current_studies = grep {$_->{config}{public}} $self->get_studies_in_sheets($species);
  my %current_study_ids = map {($_->{study_id}, 1)} @current_studies;
  my @new_studies = map {&Production::CurationDefaults::study(%$_)} grep { not $current_study_ids{$_} } @public_study_records;
  my ($rejected_study_ids, $saved_study_ids) = $self->update_sheets_with_incoming_studies($species, @new_studies);

  # Remaining:
  # Work out analyses to do, and do them for all current studies, populating the production directory - needs input path locations and production directory 
  # Deployment directory - link between production directory results, a corner of FTP where they serve, and where the paths should go
  # Go through done / pending curation / rejected studies, and make a HTML page

  # I haven't yet worked out analyses ran for a study - this should be part of study config, with a file name and description
  # We need them to provide enough info for
  #   - what analyses should run?
  #   - production->deployment sync
  #   - HTML rendering with links pointing to correct places
}
1;
