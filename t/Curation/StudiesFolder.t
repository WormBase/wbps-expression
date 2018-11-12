
use Test::More;
use Test::MockModule;
use File::Path qw/make_path/;
use File::Temp qw/tempdir/;
use File::Slurp qw/write_file/;

use Curation::StudiesFolder;

my $rnaseq = Test::MockModule->new('PublicResources::Rnaseq');

my $species = "schistosoma_mansoni";
my $assembly = "Schisto_v7";
my $study_id = "SRP001337";
sub enqueue_all_creates_directory_structure {
  my ($studies, %args) = @_;
  $args{expected} //= [map {$_->{study_id}} @$studies];
  $args{test_name} //= join ", ", @{$args{expected}};
  
  $rnaseq->mock('get', sub {return {}, {}, @$studies});
  my $tmp_dir = tempdir(CLEANUP => 1);
  make_path "$tmp_dir/curation/studies/$species/$_" for @{$args{curated_studies}};
  make_path "$tmp_dir/curation/ignore_studies";
  write_file("$tmp_dir/curation/ignore_studies/$species.tsv", join "\n", @{$args{ignore_studies}}) if $args{ignore_studies};  
  my $subject = Curation::StudiesFolder->new(tempdir(CLEANUP => 1), $tmp_dir);
  
  $subject->enqueue_all($species, $assembly); 
  
  is_deeply([$subject->list_queue($species, $assembly)], $args{expected}, $args{test_name});
}
enqueue_all_creates_directory_structure([], test_name => "Null case");
enqueue_all_creates_directory_structure([{study_id => $study_id}], test_name => "One study");
enqueue_all_creates_directory_structure(
  [{study_id => $study_id}],
  expected => [],
  ignore_studies => [$study_id],
  test_name => "Can ignore studies",
);
enqueue_all_creates_directory_structure(
  [{study_id => $study_id}],
  expected => [],
  curated_studies => [$study_id],
  test_name => "Can skip curated studies");

done_testing; 
