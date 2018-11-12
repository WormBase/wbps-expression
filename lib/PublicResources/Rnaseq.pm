
package PublicResources::Rnaseq;

use PublicResources::Resources::ArrayExpressMetadata;
use PublicResources::Resources::RnaseqerMetadata;
use PublicResources::Resources::EnaMetadata;
use PublicResources::Resources::GeoMetadata;
use PublicResources::Resources::Factors;
use PublicResources::Resources::RnaseqerStats;
use PublicResources::Resources::PubMed;
use PublicResources::Descriptions;
use PublicResources::Links;
sub new {
  my ($class, $root_dir, $src_dir) = @_;
  $src_dir //= dirname(dirname(dirname(__FILE__)));
  bless {
    root_dir => $root_dir,
    src_dir => $src_dir,
    links => 'PublicResources::Links',
  }, $class; 
}
sub get {
  my ($self, $species, $assembly) = @_;
  my $root_dir = $self->{root_dir};
  my $src_dir =  $self->{src_dir}; 
  $species = lc($species);
  $species =~ s/([a-z]*_[a-z]*).*/$1/;
  my $rnaseqer_metadata = PublicResources::Resources::RnaseqerMetadata->new($root_dir, $species);
  my $array_express_metadata = PublicResources::Resources::ArrayExpressMetadata->new($root_dir, $species);
  my $ena_metadata = PublicResources::Resources::EnaMetadata->new($root_dir, $species, $rnaseqer_metadata); 
  my $rnaseqer_stats = PublicResources::Resources::RnaseqerStats->new($root_dir, $species, $rnaseqer_metadata); 
  my $geo_metadata = PublicResources::Resources::GeoMetadata->new($root_dir, $species, $rnaseqer_metadata); 
  my $pubmed = PublicResources::Resources::PubMed->new($root_dir, $species, {
     rnaseqer=>$rnaseqer_metadata,
     array_express=>$array_express_metadata,
     ena=>$ena_metadata,
     geo=>$geo_metadata,
  });
  my $factors = PublicResources::Resources::Factors->new($root_dir, $species, $rnaseqer_metadata, $array_express_metadata);
  my $descriptions = PublicResources::Descriptions->create_for_species($src_dir, $species);
  my @studies;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    unless ($ena_metadata->{$assembly}{$study_id}){
       print STDERR "Study $study_id not in ENA, skipping\n";
       next;
    }
    my @runs;
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
       my $stats = $rnaseqer_stats->get_formatted_stats($run_id);
       my $links = $self->{links}->misc_links($study_id,$run_id, $rnaseqer_metadata->data_location($run_id),
         [keys %{$pubmed->{$assembly}{$study_id} || {}}]
       );
       my %attributes;
       for my $characteristic_type (@{$rnaseqer_metadata->access($assembly, $study_id, $run_id)}){
         $attributes{$characteristic_type} = $rnaseqer_metadata->access($assembly, $study_id, $run_id, $characteristic_type);
       }
       my ($run_description_short, $run_description_full) =
          $descriptions->run_description( $study_id, $run_id, $factors, \%attributes);
       push @runs, {
          run_id => $run_id,
          attributes => {%$stats, %$links, %attributes},
          run_description_short => $run_description_short,
          run_description_full => $run_description_full,
       };
    }
    my ($study_description_short, $study_description_full) =
         $descriptions->study_description($study_id, $ena_metadata->{$assembly}{$study_id});

    push @studies, {
      study_id => $study_id,
      runs => \@runs,
      study_description_short => $study_description_short,
      study_description_full => $study_description_full, 
      attributes => $ena_metadata->{$assembly}{$study_id}{attributes},
      pubmed => $pubmed->{$assembly}{$study_id},
    };
  }
  return $factors, $rnaseqer_metadata->{location_per_run_id}, @studies;
}
1;
