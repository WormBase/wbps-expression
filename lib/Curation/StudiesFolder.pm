
package Curation::StudiesFolder;
use PublicResources::Rnaseq;
use File::Basename qw/dirname/;
use File::Slurp qw/read_dir/;
use File::Path qw/make_path/;
my $IGNORE_STUDIES_CANNED = {
};
sub new {
  my ($class, $root_dir, $src_dir, $ignore_studies) = @_;
  $src_dir //= dirname(dirname(dirname(__FILE__)));
  $ignore_studies //= $IGNORE_STUDIES_CANNED; 
  return bless {
     processing_path => "$root_dir/curation",
     curated_studies_path => "$src_dir/curation",
     ignore_studies => $ignore_studies,
     public_rnaseq_studies => PublicResources::Rnaseq->new($root_dir),
  }, $class;
}

sub get_curated_studies {
  my ($self, $species) = @_;
  my $dir = join("/", $self->{curated_studies_path}, $species);
  return -d $dir ? read_dir $dir : ();
}

sub enqueue_all {
  my ($self, $species, $assembly) = @_;
  
  my @curated_studies = $self->get_curated_studies($species);

  my $dir = join("/", $self->{processing_path}, $species, $assembly);
  make_path $dir;

  my ($factors, $location_per_run_id, @studies) = 
    $self->{public_rnaseq_studies}->get($species, $assembly);

  
  for my $study (@studies){
     next if grep {$_ eq $study->{study_id}} @{$self->{ignore_studies}->{$species}};
     next if grep {$_ eq $study->{study_id}} @curated_studies;
     # make a nice TSV file
     # make a nice yaml file 
     make_path join("/", $dir, $study->{study_id});
  } 
}
sub list_queue {
  my ($self, $species, $assembly) = @_;
  return read_dir join("/", $self->{processing_path}, $species, $assembly);
}
1;
