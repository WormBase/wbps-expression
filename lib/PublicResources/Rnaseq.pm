use strict;
use warnings;
package PublicResources::Rnaseq;

use PublicResources::Resources::RnaseqerMetadata;
use PublicResources::Resources::EnaMetadata;
use PublicResources::Resources::GeoMetadata;
use PublicResources::Resources::RnaseqerFtp;
use PublicResources::Resources::PubMed;
use PublicResources::Descriptions;
use PublicResources::Links;
use Production::Sheets;
use Model::Design;
use File::Basename qw/dirname/;
use List::Util qw/pairmap/;
use List::MoreUtils qw/duplicates/;
#use Smart::Comments;
sub new {
  my ($class, $root_dir, $sheets) = @_;
  $sheets //= Production::Sheets->new(dirname(dirname(dirname(__FILE__))));
  bless {
    root_dir => $root_dir,
    sheets => $sheets,
    links => 'PublicResources::Links',
  }, $class; 
}
sub get {
  my ($self, $species, $assembly) = @_;
  my $root_dir = $self->{root_dir};
  $species = lc($species);
  $species =~ s/([a-z]*_[a-z]*).*/$1/;
  my $rnaseqer_metadata = PublicResources::Resources::RnaseqerMetadata->new($root_dir, $species);
  my $ena_metadata = PublicResources::Resources::EnaMetadata->new($root_dir, $species, $rnaseqer_metadata); 
  my $rnaseqer_ftp = PublicResources::Resources::RnaseqerFtp->new($root_dir, $species, $rnaseqer_metadata); 
  my $geo_metadata = PublicResources::Resources::GeoMetadata->new($root_dir, $species, $rnaseqer_metadata); 
  my $pubmed = PublicResources::Resources::PubMed->new($root_dir, $species, {
     rnaseqer=>$rnaseqer_metadata,
     ena=>$ena_metadata,
     geo=>$geo_metadata,
  });
  my $descriptions = PublicResources::Descriptions->new($species, $self->{sheets}->double_hash('run_descriptions', $species));
  my %stored_characteristics = pairmap {$a => Model::Design::from_tsv($b)->characteristics_per_run } %{$self->{sheets}->tsvs_in_folders('studies', $species)};
  my @studies;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    unless ($ena_metadata->{$assembly}{$study_id}){
       print STDERR "Study $study_id not in ENA, skipping\n";
       next;
    }
    my %characteristics_per_run = map {my $run_id = $_; 
      $run_id => ($stored_characteristics{$study_id}{$run_id} // $rnaseqer_metadata->access_characteristics($assembly, $study_id, $run_id) // {})
    } @{$rnaseqer_metadata->access($assembly, $study_id)};
    
    my @characteristic_types_varying_in_study = do {
       my %flattened = map {pairmap {"$a\t$b" => $a} %{$_}} values %characteristics_per_run;
       duplicates values %flattened
    };

    my @runs;
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
### $run_id
       my $stats = $rnaseqer_ftp->get_formatted_stats($run_id);
       my $data_location = $rnaseqer_metadata->data_location($run_id);
       my $links = $self->{links}->misc_links($study_id,$run_id, $data_location);
       my %characteristics = %{$characteristics_per_run{$run_id}};
       my %characteristics_varying_in_study = map {$_ => $characteristics{$_}} @characteristic_types_varying_in_study;
       my ($run_description_short, $run_description_full) =
          $descriptions->run_description( $study_id, $run_id, \%characteristics_varying_in_study || \%characteristics);
       push @runs, {
          run_id => $run_id,
          qc_issues => $rnaseqer_ftp->get_qc_issues($run_id),
          data_files => $rnaseqer_ftp->{$run_id}{files},
          characteristics => \%characteristics,
          attributes => {%$stats, %$links, %characteristics},
          run_description_short => $run_description_short,
          run_description_full => $run_description_full,
       };
    }
    my ($study_description_short, $study_description_full) =
         $descriptions->study_description($study_id, $ena_metadata->{$assembly}{$study_id});

    push @studies, {
      study_id => $study_id,
      runs => \@runs,
      study_description_short => $study_description_short,
      study_description_full => $study_description_full, 
      attributes => $ena_metadata->{$assembly}{$study_id}{attributes},
      pubmed => $pubmed->{$assembly}{$study_id},
    };
  }
  return @studies;
}
1;
