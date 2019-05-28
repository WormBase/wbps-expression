use strict;
use warnings;

package WbpsExpression::IncomingStudies::StudyMetadata::GeoClient;

use YAML;
use LWP;
use XML::Simple;
use Log::Any '$log';

my $EUTILS_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

sub get_xml {
  my ($url) = @_;
  $log->info("GeoClient get_xml LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;
  return XMLin( $response->decoded_content );
}

sub get_pubmed_ids {
  my ($study_id) = @_;
  my ( $web, $key ) = session_bits_from_esearch_payload(
    get_xml(
      "$EUTILS_URL/esearch.fcgi?db=gds&term=${study_id}\[accn]&usehistory=y")
  );
  return unless $web and $key;
  return esummary_xml_to_pubmed_ids(
    get_xml("$EUTILS_URL/esummary.fcgi?db=gds&query_key=$key&WebEnv=$web") );
}

sub session_bits_from_esearch_payload {
  my $payload = shift;
  return if $payload->{WarningList}{OutputMessage} eq 'No items found.';

  my $web = $payload->{WebEnv};
  my $key = $payload->{QueryKey};

  die Dump($payload) unless $web and $key;
  return $web, $key;
}

sub esummary_xml_to_pubmed_ids {
  my $payload = shift;
  return map { $_->{content} || () }
    map { my $o = $_->{Item}; ref $o eq 'ARRAY' ? @$o : $o }
    grep ( { $_->{Name} eq 'PubMedIds' } @{ $payload->{DocSum}{Item} } );
}
