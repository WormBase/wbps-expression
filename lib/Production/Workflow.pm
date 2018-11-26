
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

sub do_everything {
  my ($self, $species, $assembly) = @_;
  my @public_study_records = $self->{public_rnaseq_studies}->get($species, $assembly);  
  my %current_studies = map {$_->{study_id}=> $_} $self->get_studies_in_sheets($species);
  my %current_ignore_studies = map {$_=>1} $self->{sheets}->list("ignore_studies", $species);

  my @rejected;
  my @updated;
  my @saved;

  INCOMING_STUDIES:
  for my $study ( map {&Production::CurationDefaults::study(%$_)} grep { not $current_study_ids{$_} } @public_study_records){
    my $current_record = $current_studies{$study->{study_id}};
    if ($current_ignore_studies{$study->{study_id}} or $self->should_reject_study($study)){
      $self->{sheets}->remove_tree($species, $study->{study_id}) if $current_record; 
      push @rejected, $study->{study_id};
      next INCOMING_STUDIES;
    } 
    if ($current_record and Model::Study::config_matches_design($current_record->{config}, $study->{design})){
       $study->{config}{slices} = $current_record->{config}{slices};
       $study->{config}{condition_names} = $current_record->{config}{condition_names};
       # Additionally, characteristics in the current record were already reused
       # because they provided sources of attributes for the runs - see PublicResources::Rnaseq
       
       $study->to_folder($self->{sheets}->path("studies", $species, $study->{study_id}));
       push @updated, $study->{study_id};
       next INCOMING_STUDIES;
    }
    $study->to_folder($self->{sheets}->path("studies", $species, $study->{study_id}));
    push @saved, $study->{study_id};
    next INCOMING_STUDIES;
  };
  if (@rejected){
    $self->{sheets}->write_list( [uniq sort(keys %current_ignore_studies, @rejected)], "ignore_studies", $species)
  }


  # Remaining:
  # Run checks on everything
  # Work out analyses to do, and do them for all current studies, populating the production directory - needs input path locations and production directory 
  # Deployment directory - link between production directory results, a corner of FTP where they serve, and where the paths should go
  # Go through done / pending curation / rejected studies, and make a HTML page

  # I haven't yet worked out analyses ran for a study - this should be part of study config, with a file name and description
  # We need them to provide enough info for
  #   - what analyses should run?
  #   - production->deployment sync
  #   - HTML rendering with links pointing to correct places
  # Some internal report on what happened:
  #   - statuses: ignored / removed / unchanged / updated / failed_data_quality / failed_analysis 
  #   - more granular than html, which should only have successful and other
}
1;
