use strict;
use warnings;

# Copied as a single file from https://github.com/wormbase/wbps-expression into WBPS web code
# Please no dependencies apart from those listed in ensembl-webcode cpanfile:
# https://github.com/Ensembl/ensembl-webcode/blob/master/cpanfile

package EnsEMBL::Web::Component::Gene::WbpsExpression;
use File::Basename;

sub new {
  my ($class, $args) = @_;
  return bless $args, $class;
}
sub from_folder {
  my ($species, $assembly, $studies_path) = @_;
  my $studies_file = "$studies_path/$species.$assembly.studies.tsv";
  my %metadata;
  open (my $fh, "<", $studies_file) or return;
  while(<$fh>){
    chomp;
    my ($accession, $title) = split "\t";
    $metadata{$accession}{title} = $title;
  }
  my @studies;
  for my $study_path (grep {-d $_}  glob("$studies_path/*")){
     my $study_id = basename $study_path;
     push @studies, {
        tpms_per_condition => "$study_path/$study_id.tpm.tsv",
        study_id => $study_id,
        study_title => $metadata{$study_id}{title},
     };
  }
  return &new(__PACKAGE__, {studies => \@studies});
}

sub get_data {
  my ($self, $gene_id) = @_;
  my @result;
  for my $study (@{$self->{studies}}){
     my ($conditions, $expression_tpm) = search_in_file($study->{tpms_per_condition}, $gene_id);
     push @result, {
       study_id => $study->{study_id},
       study_title => $study->{study_title},
       conditions => $conditions,
       expression_tpm => $expression_tpm,
     } if $conditions and $expression_tpm;
  }
  return \@result;;
}

sub search_in_file {
   my ($path, $gene_id) = @_;
   my $l = `grep --max-count=1 "^$gene_id" $path`;
   chomp $l;
   return unless $l;
   my ($id, @xs) = split "\t", $l;
   return unless $id eq $gene_id and @xs;
   my $h = `grep --max-count=1 "^\t" $path`;
   chomp $h;
   return unless $h;
   my ($blank, @hs) = split "\t", $h;
   return unless not($blank) and @hs;
   return \@hs, \@xs;
}
1;
