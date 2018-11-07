
use Test::More;
use Test::MockModule;
use File::Path qw/make_path/;
use File::Temp qw/tempdir/;
use Curation::StudiesFolder;

my $rnaseq = Test::MockModule->new('PublicResources::Rnaseq');

my $studies_path_in_repo = Curation::StudiesFolder->new("temp dir")->{curated_studies_path};
ok(-d $studies_path_in_repo, "Curation folder in repository") or diag $studies_path_in_repo;

my $species = "schistosoma_mansoni";
my $assembly = "Schisto_v7";
my $study_id = "SRP001337";
sub enqueue_all_creates_directory_structure {
  my ($studies, %args) = @_;
  $args{expected} //= [map {$_->{study_id}} @$studies];
  $args{test_name} //= join ", ", @{$args{expected}};
  

  $rnaseq->mock('get', sub {return {}, {}, @$studies});
  my $subject = Curation::StudiesFolder->new(tempdir(CLEANUP => 1), $args{src_dir}, $args{ignore_studies});
  
  $subject->enqueue_all($species, $assembly); 
  
  is_deeply([$subject->list_queue($species, $assembly)], $args{expected}, $args{test_name});
}
enqueue_all_creates_directory_structure([], test_name => "Null case");
enqueue_all_creates_directory_structure([{study_id => $study_id}]);
enqueue_all_creates_directory_structure(
  [{study_id => $study_id}],
  expected => [],
  ignore_studies => {$species => [$study_id]},
  test_name => "Can ignore studies",
);
my $tmp_src_dir = tempdir(CLEANUP => 1);
make_path join ("/", $tmp_src_dir, "curation",$species, $study_id);
enqueue_all_creates_directory_structure(
  [{study_id => $study_id}],
  expected => [],
  src_dir => $tmp_src_dir, 
  test_name => "Can skip curated studies");

done_testing; 
