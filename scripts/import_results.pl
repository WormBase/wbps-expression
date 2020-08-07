use strict;
use warnings;
use File::Basename;
use File::Slurp qw/write_file/;
use File::Path qw/make_path/;
use PublicResources::Resources::RnaseqerMetadata;
use Smart::Comments;

for my $d (glob "Resources/*"){
  my $species = basename $d;
  my $api = PublicResources::Resources::RnaseqerMetadata->new("Resources", $species);
  for my $assembly (@{$api->access}){
     for my $study_id (@{$api->access($assembly)}){
        my $path = "/Users/wb4/dev/wbps-expression/curation/studies/$species/$study_id/$study_id.source_dirs.tsv";
        make_path (dirname $path);
        write_file($path,
        "Run\tLocation\n",
        map {
            join ("\t", $_, ($api->{location_per_run_id}{$_} =~ s{/[^/]*.nospliced.bw}{}r))."\n"
        } map {@{$api->access($assembly, $study_id, $_)}} @{$api->access($assembly, $study_id)});
     }
  }
}
