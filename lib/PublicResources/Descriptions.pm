use strict;
use warnings;
package PublicResources::Descriptions;

sub new {
   my ($class, $species, $curations) = @_;
   return bless {
      species => $species,
      curations => $curations, 
   }, $class;
}

# Duplicated with JBrowseDisplay.pm
my @types_blacklist = (
 "synonym",
 "bioproject_id",
  "species",
  "organism",
  "replicate",
  "sample_name",
  "batch",
  "barcode",
  "insdc_center_name",
  "insdc_first_public",
  "insdc_secondary_accession",
  "insdc_status",
  "insdc_last_update",
  "label",
  "model",
  "package",
  "ncbi_submission_model",
  "ncbi_submission_package",
  "sample_comment",
  "sample_title",
  "geo_accession",
  "biological_replicate",
  "block",
  "zone", #schmidtea mediterranea
  "repplicate",
  "in_house_sample_code",
  "collected_by",
  "biomaterial_provider",
  "description_title",
  "treatment_sources",
  "population",
  "sample_name",
  "agarosemigrationtemperature",
  "agarosemigrationttime",
  "baermanntemperature",
  "base_calling_software_version",
  "culturetemperature",
  "culturetime",
  "library_id",
  "library_preparation",
  "wash",
);

# if some of these types can be made nicer in labels,
# e.g. if they have _, add some code that prettifies them
my @types_that_help_explain_values = qw/strain isolate/;

sub run_description_from_attributes {
    my ( $attributes ) = @_;
    my @result_short;
    my @result_full;
    for my $t (sort keys %{$attributes}) {
        my $v = $attributes->{$t};
        next if grep {$_ eq $t} @types_blacklist;
        next unless $v;
        push @result_short, $v;
        if ( grep { $_ eq $t } @types_that_help_explain_values ) {
            push @result_full, "$t $v";
        }
        else {
            push @result_full, $v;
        }
    }
    return "" unless @result_short and @result_full;
    return [ join( ", ", @result_short ), join( ", ", @result_full ) ];
}

sub _get_run_description {
    my ( $self, $study_id, $run_id, $attributes ) = @_;
    return (
        $self->{curations}{$study_id}{$run_id}
          or $attributes->{sample_name}
          or run_description_from_attributes( $attributes )
          or ""
    );
}

sub run_description {
    my ( $self, $study_id, $run_id, $attributes ) = @_;
    my $species = $self->{species};
    $species =~ s/_/ /g;
    $species = ucfirst($species);
    my $r = _get_run_description(@_);
    $r = [$r, $r] unless ref $r eq 'ARRAY';
	my ($r1, $r2) = @{$r};
	return $r1 ? $r2 ? ($r1, $r2) : ($r1, $r1) : ( "", "sample from $species");
}
sub _clean_study_description {
  my ($study_description) = @_;
  if( $study_description =~ m{http://www.sanger.ac.uk/datasharing}){
     $study_description =~ s/\s+/ /g;
     $study_description =~ s{This data is part of a pre-publication release.*http://www.sanger.ac.uk/datasharing\W+}{};
  }
  return "" if length($study_description) > 500;
  return "" if length($study_description) < 5;
  die $study_description if $study_description =~ /This data is part of a pre-publication release/;
  return $study_description;
}
sub _clean_study_title {
   my $name = shift;
   $name =~s/_+/ /g if scalar(split " ", $name) == 1;
   $name = ucfirst(lc($name)) if $name eq uc($name);
   return $name;
}
sub study_description {
    my ( $self, $study_id, $study_metadata ) = @_;
    my $species = $self->{species};
    $species =~ s/_/ /g;
    $species = ucfirst($species);
    my $short_description = $study_metadata->{study_title} ? _clean_study_title($study_metadata->{study_title}) : "$species study";
    my $full_description = _clean_study_description($study_metadata->{study_description} // "");
    $full_description ||= $short_description;
    return $short_description, $full_description;
}
1;
