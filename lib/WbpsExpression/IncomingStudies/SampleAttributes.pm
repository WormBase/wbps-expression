use strict;
use warnings;

package WbpsExpression::IncomingStudies::SampleAttributes;
use List::Util qw/pairs/;
use List::MoreUtils qw(uniq);
use LWP;
use JSON;
use Log::Any '$log';

sub characteristics_by_run {
  my ($species, $study_id) = @_;

  my %result;
  for my $o (@{ get_json("https://www.ebi.ac.uk/fg/rnaseq/api/json/getSampleAttributesPerRunByStudy/$study_id") }){
    $result{$o->{RUN_ID}}{$o->{TYPE}} = $o->{VALUE};
  }
  for my $run_id (keys %result){
    $result{$run_id} = standardise_for_worm($species, $result{$run_id});
  }
  return \%result;
}

sub get_json {
  my ($url) = @_;
  $log->info("SampleAttributes get_json LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;
  return decode_json($response->decoded_content);
}

sub standardise_for_worm {
  my ($species, $sample_attributes) = @_;
  $species= lc($species);
  $species =~ s/[^a-zA-Z0-9]+/_/g;
  my %result;
  for my $p (pairs %{$sample_attributes}){
    my ($type, $value) = _normalise_type_and_value(@{$p});
    if($type eq "sample_name"){
 	  $value = fix_sample_name($species, $value);
    }
    $result{$type} = $value if $type and $value;
  }
  return put_sex_into_separate_characteristic_if_present_in_dev_stage_or_organism_part(\%result);
}

my @type_blacklist = (
 #Metadata
  "bioproject_id",
  "species",
  "organism",
 # Replicate - frequently 1,2,3
  "replicate",
  "repplicate",
  "replication",
  "biological_replicate",
 # GEO->ENA import artifacts
  "ncbi_submission_model",
  "ncbi_submission_package",
  "geo_accession", 
  "insdc_center_name",
  "insdc_first_public",
  "insdc_secondary_accession",
  "insdc_status",
  "insdc_last_update",
 # Software things
  "package",
  "library_id",
  "library_preparation",
  "base_calling_software_version",
);
sub _normalise_type_and_value {
  my ($type, $value) = @_;

  $type = lc($type);
  $type =~ s/^\s+|\s+$//g;
  $type =~ s/\W+/_/g;
  return "","" if grep {$_ eq $type} @type_blacklist;
# Reject accessions of sample names - regex originally from https://www.ebi.ac.uk/ena/submit/accession-number-formats
# SRS\d+ are NCBI sample names I think, although they don't commit to a regex
# Sometimes in "synonym" field
  return "","" if $value =~ /^(E|S|D)RS\d{6,}$/;
  return "","" if $value =~ /^SAM(E|D|N)[A-Z]?\d+$/;
  return "","" if $value =~ /^GSM\d+$/;
#Sometimes there's curation like: age+time unit
  return $type, $value if $type eq "age" and $value =~s/^\W+$//;
# But sometimes age means developmental stage
  $type =~ s/^age$/developmental_stage/;

# Try give developmental stages first-class support
  $type =~ s/^stage$/developmental_stage/;
  $type =~ s/life_cycle_stage/developmental_stage/;
  $type =~ s/dev_stage/developmental_stage/;
  $type =~ s/development_stage/developmental_stage/;
  if( $type eq "developmental_stage"){
     $value =~ s/^Adults?/adult/i;
     $value =~ s/^adult worms?$/adult/i;
  }
  $type =~ s/^tissue$/organism_part/;
  if($type eq "organism_part" ){
     $value = lc($value);
#http://purl.obolibrary.org/obo/UBERON_0000468
     $value =~ s/^whole$/whole organism/;
     $value =~ s/^whole body$/whole organism/;
     $value =~ s/^whole ?-?_?worms? ?-?_?(tissue)?$/whole organism/;
     $value =~ s/^whole ?-?_?organisi?m$/whole organism/;
     $value =~ s/^whole ?-?_?animals?$/whole organism/;
     $value =~ s/^intact ?-?_?animals?$/whole organism/;
  }
  if($type eq "host_disease"){
     $value = lc($value);
     $value =~ s/^schistosomiase$/schistosomiasis/;
  }
  if($type eq "genotype"){
     $value =~ s/^(wide|wild) ?-?_?type ?-?_?(genotype)?$/wild type/i;
     $value =~ s/^WT$/wild type/i;
  }
  $type =~ s/time_point/timepoint/;

  $value =~s/^not applicable.*$//i;
  $value =~s/^unknown$//i;
  $value =~s/^N\\?A$//i;
  $value =~s/^\W+$//;

  $value =~ s/^\s+|\s+$//g;
  return $type, $value;
}
# Sample name is important so give it special treatment
# We want it to be a run description and a meaningful guess at conditions
# Unfortunately someone out there has a default based on species name,
# also GEO announces stuff in the field, etc.
sub fix_sample_name {
  my ($species, $sample_name) = @_;
  return "" unless $sample_name;
  ( my $cies = $species ) =~ s/.*_//;
  return ""
    if (  $cies
      and $sample_name =~ /$cies/i
      and $sample_name =~ /sample from/i );
  return "" if $sample_name =~ /private and is scheduled to be released/;

  $sample_name =~ s/\[\w+ $cies RNAseq\]//;#https://www.ebi.ac.uk/ena/data/view/DRS026763&display=xml
  $sample_name =~ s/^\s+//;
  $sample_name =~ s/\s+$//;
  $sample_name =~ s/\s*(biological)? replicate \d+$//;
  return $sample_name;
}
sub put_sex_into_separate_characteristic_if_present_in_dev_stage_or_organism_part {
  my $o = shift;
  if(not $o->{sex} and $o->{developmental_stage} and $o->{developmental_stage} =~ /^adult ?-?_?males?$/){
     $o->{sex} = "male";
     $o->{developmental_stage} = "adult";
  }
  if(not $o->{sex} and $o->{developmental_stage} and $o->{developmental_stage} =~ /^adult ?-?_?females?$/){
     $o->{sex} = "female";
     $o->{developmental_stage} = "adult";
  }
  if(not $o->{sex} and $o->{developmental_stage} and $o->{developmental_stage} =~ /^adult ?-?_?males? ?-?_?and ?-?_?females?$/){
     $o->{sex} = "pooled male and female";
     $o->{developmental_stage} = "adult";
  }
  if(not $o->{sex} and $o->{organism_part} and $o->{organism_part} =~ /^whole ?-?_?males?$/){
     $o->{sex} = "male";
     $o->{organism_part} = "whole organism";
  }
  if(not $o->{sex} and $o->{organism_part} and $o->{organism_part} =~ /^whole ?-?_?females?$/){
     $o->{sex} = "female";
     $o->{organism_part} = "whole organism";
  }
  return $o;
}
1;
