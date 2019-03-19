package PublicResources::Resources::PubMed;
use parent PublicResources::Resources::LocallyCachedResource;
use PublicResources::Resources::RnaseqerMetadata;
use List::MoreUtils qw(uniq);
use XML::Simple;

my $curation = {
  hymenolepis_microstoma => {
    ERP004459 => [30455861], 
  },
  heligmosomoides_polygyrus => {
    SRP157940 => [30349532],
  }, 
  meloidogyne_incognita => {
    ERP009887 => [28594822],
    SRP109232 => [29230237],
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
sub _fetch {
    my ( $class, $species, $metadata ) = @_;
    my %data;
    for my $assembly( @{$metadata->{rnaseqer}->access}){
       for my $study_id ( @{ $metadata->{rnaseqer}->access($assembly) } ) {
          my $ena_study_pubmed_ids = $metadata->{ena}{$assembly}{$study_id}{study_pubmed} // [];
          my $ena_bioproject_pubmed_ids = $metadata->{ena}{$assembly}{$study_id}{bioproject_pubmed} // [];
          my $geo_pubmed_ids = $metadata->{geo}{$assembly}{$study_id}{pubmed} // [];
          my $curated_ids = $curation->{$species}{$study_id} // [];
          for my $pubmed_id ( uniq(@$ena_study_pubmed_ids, @$ena_bioproject_pubmed_ids, @$geo_pubmed_ids, @$ae_pubmed_ids, @$curated_ids)){
              next if $pubmed_id eq '2971468'; # Summer 2019 or later? Check if PRJNA392315 still refers to this paper in error
              $data{$assembly}{$study_id}{$pubmed_id} = &_short_and_full_paper_description_from_payload($species, $class->get_xml(
                   "https://www.ncbi.nlm.nih.gov/pubmed/$pubmed_id?report=xml&format=text"
              ));
          } 
       }
    }
    return \%data;
}
sub italicise_species_in_title {
  my ($species, $title) = @_;
  return unless $title;
  my $species = lc $species;
  my ($spe, $cies) = split / |_/, $species;
  my $s = substr($spe, 0, 1);
  $title =~ s{(?:<i>)?($spe\.?[ _]$cies|$s\.?[ _]$cies)(?:(</i>)?)}{<i>$1</i>}gi;
  return $title;
}

sub _short_and_full_paper_description_from_payload {
    # PubMed formats this as string to encourage people to use their API
    # We are not encouraged enough, so we're going to parse twice.
    my ($species, $payload_string) = @_;
    my $payload = XMLin($payload_string);

    my @authors = @{$payload->{MedlineCitation}{Article}{AuthorList}{Author} || [] };
    # Use regex because XML::Simple is being too simple.
    # E.g. 30049782: <ArticleTitle> Stuff in <i>Caenorhabditis elegans</i>.</ArticleTitle>
    my ($title) = $payload_string =~ m{<ArticleTitle>(.*)</ArticleTitle>};
    my $title = italicise_species_in_title($species, $title);
    
    my $authors = $payload->{MedlineCitation}{Article}{AuthorList}{Author};
    my @authors = $authors ? ref $authors eq 'ARRAY' ? @$authors : ($authors) : ();
    my $first_author = @authors[0]->{LastName};
    my $last_author = @authors[-1]->{LastName};
    my $authors = $first_author ? $last_author ne $first_author ? "$first_author .. $last_author" : $first_author :  "";
    my $year = $payload->{MedlineCitation}{Article}{Journal}{JournalIssue}{PubDate}{Year};
    my $short_description = "$authors, $year";
    my $full_description = $title ? "$title ($authors, $year)" : $short_description;
    return [$short_description, $full_description];
}
1;
