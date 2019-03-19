use strict;
use warnings;
package PublicResources::Resources::RnaseqerMetadata;
use List::MoreUtils qw(uniq);
use File::Basename;
use parent 'PublicResources::Resources::LocallyCachedResource';

# Assembly -> study -> run -> type -> value

sub access {
  my $self = shift;
  my $h = $self->{metadata};
  while(@_ and ref $h){
    $h = $h->{shift @_};
  }
  return [sort keys %$h ] if ref $h;
  return [] unless $h;
  return $h;
}

sub access_characteristics {
  my ($self, $assembly, $study_id, $sample_id, $run_id) = @_;
  return $self->{metadata}{$assembly}{$study_id}{$sample_id}{$run_id};
}

sub data_location {
  my ($self, $run_id) = @_; 
  my $bigwig_location = $self->{location_per_run_id}{$run_id};
  return dirname $bigwig_location; 
}

sub _fetch {
  my ($class, $species) = @_;
  $species= lc($species);
  $species =~ s/[^a-zA-Z0-9]+/_/g;

  my $run_records = $class->_get_rnaseqer_runs_for_organism($species);

  my %run_attributes;
  for my $study_id (uniq (map {$_->{STUDY_ID}} @$run_records)){
     for my $attribute_record (@{ 
       $class->_get_rnaseqer_sample_attributes_per_run_for_study($study_id)
     }){
       my ($type, $value) = _normalise_type_and_value($attribute_record->{TYPE}, $attribute_record->{VALUE});
       if($type eq "sample_name"){
          $value = fix_sample_name($species, $value);
       }
       $run_attributes{$attribute_record->{RUN_ID}}{$type} = $value if $type and $value;
     }
     for my $run_id (keys %run_attributes){
        $run_attributes{$run_id} = _normalise_characteristics_hash($run_attributes{$run_id});
     }
  }
  my %data;
  for my $run_record (@$run_records){
      $data
        {$run_record->{ASSEMBLY_USED}}
        {$run_record->{STUDY_ID}}
        {$run_record->{SAMPLE_IDS}}
        {$run_record->{RUN_IDS}}
        = $run_attributes{$run_record->{RUN_IDS}}; 
  }

  my %location_per_run_id;
  for my $run_record (@$run_records){
  (my $bigwig_we_want = $run_record->{BIGWIG_LOCATION}) =~ s/.bw$/.nospliced.bw/;
      $location_per_run_id
         {$run_record->{RUN_IDS}}
         = $bigwig_we_want;
  }
  return {metadata => \%data, location_per_run_id => \%location_per_run_id};
}
# Stuff that won't be helpful to see, and also not helpful when curating the experiment
# Be conservative when putting stuff here, because how do we know?
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
  $type =~ s/\W+/_/g;
#each run has a sample and we can look it up in ENA but it's not a characteristic so filter it
  return "","" if grep {$_ eq $type} @type_blacklist;
# Reject accessions of sample names - regex originally from https://www.ebi.ac.uk/ena/submit/accession-number-formats
# SRS\d+ are NCBI sample names I think, although they don't commit to a regex
# Sometimes in "synonym" field
  return "","" if $value =~ /^(E|S)RS\d{6,}$/;
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
sub _normalise_characteristics_hash {
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
sub _get_rnaseqer_json {
  my $class = shift;
  return $class->get_json(
    join ("/",
      "https://www.ebi.ac.uk/fg/rnaseq/api/json",
      @_
    )
  );
}
# Ask for runs that had at least 20% of the reads mapped
# We want to exclude failures and queued entries
# RNASeq-er doesn't have complete records anyway!
sub _get_rnaseqer_runs_for_organism {
  my ($class, $species) = @_;
  my $payload = $class->_get_rnaseqer_json(
    "20", "getRunsByOrganism", $species
  );
  my @a = grep {$_->{BIGWIG_LOCATION} ne "NA"} @$payload;
  return \@a;
}
sub _get_rnaseqer_sample_attributes_per_run_for_study {
  my ($class, $study) = @_;
  return $class->_get_rnaseqer_json(
    "getSampleAttributesPerRunByStudy", $study
  );
}
1;
