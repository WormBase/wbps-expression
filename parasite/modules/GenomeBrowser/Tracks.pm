
package GenomeBrowser::Tracks;

use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::Factors;
use GenomeBrowser::RnaseqerMetadata;
# This is the convenient representation of tracks
# Fetches its dependent data and puts them on disk

# you can edit the dependent data by hand, like configs

# Convert me to JBrowse format
# Dump me to yaml to see what I'm like
# Use me to move files EBI -> Sanger

# TODO 
# Figure out what's up with hubs 
# Do we need more data fetched from somewhere?
# make yourself a FactorsCurationWizard
sub new {
  my ($class, $root_dir, $species, $assembly, $rnaseqer_metadata, $factors) = @_;
  my $rnaseqer_metadata = GenomeBrowser::RnaseqerMetadata->new($root_dir, $species);
  my $array_express_metadata = GenomeBrowser::ArrayExpressMetadata->new($root_dir, $species);
  my $factors = GenomeBrowser::Factors->new($root_dir, $species, $assembly, $rnaseqer_metadata, $array_express_metadata);
  my %data;
  for my $study_id ($rnaseqer_metadata->access($assembly)){
    for my $run_id ($rnaseqer_metadata->access($assembly, $study_id)){
       my %attributes = %{$rnaseqer_metadata->{$assembly}{$study_id}{$run}};
       $attributes{Path} = &_public_url($run_id);
       $attributes{Label} = "$study_id/$run_id\t".join (", ",
          map {$_ if $_} @attributes{@$factors}
       );
       $data{$run_id}=\%attributes;
    }
  }
  return bless \%data, $class;
}

sub list_names {
   return keys shift;
}

sub source_path {
  my ($self, $track_id) = @_;
  return "ebi:/arrayexpress_folder/$track_id.bw";
}

sub storage_path {
  my ($self, $track_id) = @_;
  return "sanger:/ngs_folder/$track_id.bw";
}


sub as_jbrowse_config {
  my $self = shift;
  return "Perl object that can be converted to tracks.conf, probably with IniFiles";
}

sub _public_url {
  my $track_id = shift;
  return "https://ngs.sanger.ac.uk/$track_id.bw";
}
1;
