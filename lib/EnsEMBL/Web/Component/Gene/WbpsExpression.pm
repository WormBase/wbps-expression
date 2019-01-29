use strict;
use warnings;

# Copied as a single file from https://github.com/wormbase/wbps-expression into WBPS web code
# Please no dependencies apart from those listed in ensembl-webcode cpanfile:
# https://github.com/Ensembl/ensembl-webcode/blob/master/cpanfile

package EnsEMBL::Web::Component::Gene::WbpsExpression;

sub new {
   # location of the data files to read from, passed in as config
}
sub from_folder {

}

sub get_data {
  my ($species, $gene_id) = @_;
  # look through the files, provide objects the UI can render
  

  # return format documented clearly
  return {};
}

1;
