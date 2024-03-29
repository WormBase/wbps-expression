use strict;
use warnings;

package WbpsExpression::IncomingStudies::StudyMetadata::PubMedClient;
use List::MoreUtils qw(uniq);
use LWP;
use XML::Simple;
use Log::Any '$log';
use Try::Tiny;
use Carp;
use Data::Dumper;
use File::Slurp qw/write_file/;

sub title_and_description_for_pubmed_id {
  my ( $species, $pubmed_id ) = @_;
  my $url =
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pubmed_id&retmode=xml";
  $log->info("PubMedClient get_xml LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;

  my $payload_string = $response->decoded_content;
  return unless $payload_string; #rare but e.g. 30962434
  my ($short_description, $full_description );

  try {
   my $payload = XMLin($payload_string);
   my $article_elt = $payload->{"PubmedArticle"}{"MedlineCitation"}{"Article"}
      || die "cannot find PubmedArticle|MedlineCitation|Article";
    
    my $title = $article_elt->{ArticleTitle}
      || die "cannot find PubmedArticle|MedlineCitation|Article|ArticleTitle";
    $title = italicise_species_in_title( $species, $title );
    print $title;
    my $authors = $article_elt->{AuthorList}{Author}
      || die "cannot find PubmedArticle|MedlineCitation|Article|AuthorList|Author";
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
    my $year = $article_elt->{Journal}{JournalIssue}{PubDate}{Year}
      || die "cannot find PubmedArticle|MedlineCitation|Article|Journal|JournalIssue|PubDate|Year";
    $short_description = "$authors, $year";
    $full_description =
      $title ? "$title ($authors, $year)" : $short_description;
  } catch {
    my $msg = $_;
    my $xmllog = "/homes/digri/tmp/pubmid-$pubmed_id-$$.xml";
    write_file($xmllog, $payload_string);
    confess "PubMed XML parsing error (XML saved as $xmllog): $msg $pubmed_id $url";
  };
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
