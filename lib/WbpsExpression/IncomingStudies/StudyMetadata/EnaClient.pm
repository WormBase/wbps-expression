use strict;
use warnings;

package WbpsExpression::IncomingStudies::StudyMetadata::EnaClient;

use LWP;
use XML::Simple;
use Log::Any '$log';
use List::Util qw/uniq/;

# use Smart::Comments '###';

# Returns a hash:
# ena_first_public
# ena_last_update
# bioproject
# pubmed_refs // []
# resource_links // []
# title // ""
# description // ""
# submitting_centre

sub get_study_metadata {
  my ($study_id) = @_;
  my $data_for_study = &xml_to_data_for_study(
    get_xml("https://www.ebi.ac.uk/ena/data/view/$study_id&display=xml") );
#### $data_for_study
  if ( $data_for_study->{bioproject} ) {
    my ( $submitting_centre, $pubmed_refs, $resource_links) =
      &xml_to_data_for_bioproject(get_xml(
        sprintf( "https://www.ebi.ac.uk/ena/data/view/%s&display=xml",
          $data_for_study->{bioproject} )
      ));

    #TODO more useful stuff: maybe the descriptions?
    $data_for_study->{pubmed_refs} = [sort {$a cmp $b} uniq ( @{$pubmed_refs // []}, @{ $data_for_study->{pubmed_refs} //[] })];
    $data_for_study->{resource_links} = [sort {$a cmp $b } uniq (@{$resource_links//[]}, @{ $data_for_study->{resource_links} //[] })];
    $data_for_study->{submitting_centre} ||= $submitting_centre;

  }
  if ( $data_for_study->{attributes}{submitting_centre} and $data_for_study->{attributes}{submitting_centre} =~ /^null$/i ) {
    delete $data_for_study->{attributes}{submitting_centre};
  }
#### $data_for_study
  return $data_for_study;
}

sub get_xml {
  my ($url) = @_;
  $log->info("EnaClient get_xml LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;
  return XMLin( $response->decoded_content );
}

sub _url_link {
  my ( $label, $url ) = @_;

  # TODO can you rely on BioProjects to make ArrayExpress links?
  # Are there any more URL links?
  ( my $property_name = $label ) =~ s/.*ArrayExpress/ArrayExpress/;

  return ( $property_name, $label, $url );
}

sub extract_pubmed_refs_and_resource_links {
  my ($links) = @_;
#### $links
  my @pubmed_refs;
  my @resource_links;
  for my $link ( @{ $links // [] } ) {
    if ( uc( $link->{XREF_LINK}{DB} //"" ) eq 'PUBMED' ) {
      push @pubmed_refs, $link->{XREF_LINK}{ID};
    }
    elsif ( $link->{URL_LINK}{LABEL}
      and $link->{URL_LINK}{URL} )
    {
      my ( $property_name, $label, $url ) =
        _url_link( $link->{URL_LINK}{LABEL},
        $link->{URL_LINK}{URL} );
      push @resource_links, [ $property_name, $label, $url ];
    }
    else {
      #probably a link to an ENA something - skip
    }
  }
#### @pubmed_refs
#### @resource_links
  return \@pubmed_refs, \@resource_links;
}

sub xml_to_data_for_bioproject {
  my $payload = shift;
  return unless $payload;
  return $payload->{PROJECT}{center_name},
    extract_pubmed_refs_and_resource_links(
    $payload->{PROJECT}{PROJECT_LINKS}{PROJECT_LINK} );
}

# Found description as an empty hash in: SRP013211
sub to_string {
  my ($o) = @_;
  $o = join ", ", pairmap { "$a: $b" } %{$o} if ref $o eq 'HASH';
  $o = join ", ", @{$o} if ref $o eq 'ARRAY';
  return $o;
}

sub determine_bioproject {
  my $payload = shift;
  my @bioprojects;
  my $ids = $payload->{STUDY}{IDENTIFIERS}{EXTERNAL_ID};

  # XML::Simple is being a bit too simple
  # many ids -> array here: ERP016356
  # one id -> hash here: SRP093920
  my @ids = ref $ids eq 'ARRAY' ? @$ids : ref $ids eq 'HASH' ? ($ids) : ();
  for my $identifier (@ids) {
    push @bioprojects, $identifier->{content}
      if uc( $identifier->{namespace} ) eq "BIOPROJECT"
      and $identifier->{label} eq "primary";
  }
  my ( $bioproject, @other_bioprojects ) = @bioprojects;
  $bioproject //= $payload->{STUDY}{IDENTIFIERS}{SECONDARY_ID};
  die join( " ",
    $payload->{STUDY}{accession},
    ": could not determine BioProject",
    $bioproject, @other_bioprojects )
    if ( not $bioproject or @other_bioprojects );

# Working assumption: everything has a bioproject. Will soon find a counterexample. :)
  return $bioproject;
}

sub xml_to_data_for_study {
  my $payload = shift;
  return {} unless $payload;
  my $ena_first_public = join( " ",
    map { $_->{TAG} eq 'ENA-FIRST-PUBLIC' ? $_->{VALUE} : () }
      @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } );
  my $ena_last_update = join( " ",
    map { $_->{TAG} eq 'ENA-LAST-UPDATE' ? $_->{VALUE} : () }
      @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } );

  my $submitting_centre = (
    uc( $payload->{STUDY}{broker_name} // "" ) eq 'NCBI'
      and length( $payload->{STUDY}{center_name} ) < 10
      or $payload->{STUDY}{center_name} eq "BioProject"
  ) ? "" : $payload->{STUDY}{center_name};

  my ( $pubmed_refs, $resource_links ) =
    extract_pubmed_refs_and_resource_links(
    $payload->{STUDY}{STUDY_LINKS}{STUDY_LINK} );

  return {
    ena_first_public => $ena_first_public,
    ena_last_update  => $ena_last_update,
    bioproject       => determine_bioproject($payload),
    pubmed_refs      => $pubmed_refs,
    resource_links   => $resource_links,
    title            => clean_study_title(to_string( $payload->{STUDY}{DESCRIPTOR}{STUDY_TITLE} )),
    description =>
      clean_study_description(to_string( $payload->{STUDY}{DESCRIPTOR}{STUDY_DESCRIPTION} )),
    submitting_centre => $submitting_centre,
  };
}
sub clean_study_description {
  my ($study_description) = @_;
  return "" unless $study_description;
  if( $study_description =~ m{http://www.sanger.ac.uk/datasharing}){
     $study_description =~ s/\s+/ /g;
     $study_description =~ s{This data is part of a pre-publication release.*http://www.sanger.ac.uk/datasharing\W+}{};
  }
  return "" if length($study_description) > 500;
  return "" if length($study_description) < 5;
  die $study_description if $study_description =~ /This data is part of a pre-publication release/;
  return $study_description;
}
sub clean_study_title {
   my ($name) = @_;
   $name =~s/_+/ /g if scalar(split " ", $name) == 1;
   $name = ucfirst(lc($name)) if $name eq uc($name);
   $name =~ s/rejuvinate/rejuvenate/;
   $name =~ s/Some RNA-seq reads form/Some RNA-seq reads from/;
   return $name;
}
1;
