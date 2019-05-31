use strict;
use warnings;
package WbpsExpression::IncomingStudies::StudyMetadata;
use WbpsExpression::IncomingStudies::StudyMetadata::EnaClient;
use WbpsExpression::IncomingStudies::StudyMetadata::GeoClient;
use WbpsExpression::IncomingStudies::StudyMetadata::PubMedClient;
use List::MoreUtils qw/uniq/;
# delegates to EnaClient except it retrieves extra IDs from GEO, and looks up descriptions in PubMed

my $pubmed_curated_ids = {
  bursaphelenchus_xylophilus => {
    ERP010245 => [25981957],
    SRP115811 => [27224277],
  },
  hymenolepis_microstoma => {
    ERP004459 => [30455861],
  },
  heligmosomoides_polygyrus => {
    SRP157940 => [30349532],
  },
  heterorhabditis_bacteriophora => {
    SRP125059 => [28049427],
  },
  meloidogyne_incognita => {
    ERP009887 => [28594822],
    SRP109232 => [29230237],
    SRP152065 => [30486772],
  },
  schistosoma_mansoni => {
    ERP000427 => [22253936],
    ERP014584 => [28542189],
    ERP016356 => [27499125],
    SRP067194 => [27003592],
    SRP071285 => [27677173],
    SRP093920 => [30733716],
    SRP096638 => [30365505],
    SRP108901 => [28753630],
    SRP124650 => [29557781],
    SRP130864 => [29649665, 30029996],
  },
};

sub get {
  my ($species, $study_id) = @_;

  my $result = WbpsExpression::IncomingStudies::StudyMetadata::EnaClient::get_study_metadata($study_id);
  return unless $result;
  my $pubmed_ids_ena = delete $result->{pubmed_refs};

  my @pubmed_ids_geo = WbpsExpression::IncomingStudies::StudyMetadata::GeoClient::get_pubmed_ids($study_id);
  my $pubmed_ids_extra = $pubmed_curated_ids->{$species}{$study_id};
  my %publication_title_and_description_per_pubmed_id = map {
    my $xs = WbpsExpression::IncomingStudies::StudyMetadata::PubMedClient::title_and_description_for_pubmed_id($species, $_);
    $xs ? ($_ => $xs) : ()
  } uniq @{$pubmed_ids_ena}, @pubmed_ids_geo, @{ $pubmed_curated_ids->{$species}{$study_id} //[]};
  $result->{pubmed} = \%publication_title_and_description_per_pubmed_id;

  $result->{title} ||= do {
     (my $species_title = $species) =~ s/_/ /g;
     $species_title = ucfirst ($species_title);
     "$species_title study"
  };
  return $result;
}
1;
