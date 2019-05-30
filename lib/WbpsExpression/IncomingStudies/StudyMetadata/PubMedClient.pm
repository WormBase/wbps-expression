use strict;
use warnings;

package WbpsExpression::IncomingStudies::StudyMetadata::PubMedClient;
use List::MoreUtils qw(uniq);
use LWP;
use XML::Simple;
use Log::Any '$log';

sub title_and_description_for_pubmed_id {
  my ( $species, $pubmed_id ) = @_;
  my $url =
    "https://www.ncbi.nlm.nih.gov/pubmed/$pubmed_id?report=xml&format=text";

  $log->info("PubMedClient get_xml LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;
  my $payload_string = XMLin( $response->decoded_content );
  return unless $payload_string; #rare but e.g. 30962434
  my $payload        = XMLin($payload_string);

# Use regex because XML::Simple is being too simple and creates too much structure
# E.g. 30049782: <ArticleTitle> Stuff in <i>Caenorhabditis elegans</i>.</ArticleTitle>
  my ($title) = $payload_string =~ m{<ArticleTitle>(.*)</ArticleTitle>};
  $title = italicise_species_in_title( $species, $title );

  my $authors = $payload->{MedlineCitation}{Article}{AuthorList}{Author};
  my @authors =
    $authors ? ref $authors eq 'ARRAY' ? @$authors : ($authors) : ();
  my $first_author = $authors[0]->{LastName};
  my $last_author  = $authors[-1]->{LastName};
  $authors =
      $first_author
    ? $last_author ne $first_author
      ? "$first_author .. $last_author"
      : $first_author
    : "";
  my $year =
    $payload->{MedlineCitation}{Article}{Journal}{JournalIssue}{PubDate}{Year};
  my $short_description = "$authors, $year";
  my $full_description =
    $title ? "$title ($authors, $year)" : $short_description;
  return [ $short_description, $full_description ];
}

sub italicise_species_in_title {
  my ( $species, $title ) = @_;
  return unless $title;
  my ( $spe, $cies ) = split / |_/, lc $species;
  my $s = substr( $spe, 0, 1 );
  $title =~
    s{(?:<i>)?($spe\.?[ _]$cies|$s\.?[ _]$cies)(?:(</i>)?)}{<i>$1</i>}gi;
  return $title;
}
1;
