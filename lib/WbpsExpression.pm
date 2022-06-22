use strict;
use warnings;
package WbpsExpression;
use WbpsExpression::IncomingStudies;
use WbpsExpression::Analysis;
use WbpsExpression::StudiesPage;
use File::Slurp qw/write_file read_dir/;
use File::Path qw/make_path/;
use FindBin;
my $studies_dir = "$ENV{PARASITE_REPOSITORIES}/studies";

our %exceptions = (
    "ERP006987"  => "C. sinensis: four runs, but the tpms could be meaningful",
    "SRP013211"  => "O. viverrini: two runs, juvenile vs adult, from 2012",
    "DRP003063"  => "Emu study, four runs",
    "SRP131874"  => "E. granulosus, four runs but a second study for the species",
    "SRP152065"  => "M. incognita IARI study, two controls two replicates",
    "SRP179824"  => "S. japonicum, three runs done for genome sequencing. Essentially a nice start",
    "SRP140458"  => "T. pseudospiralis, three runs for genome sequencing, different life stages",
    "SRP048819"  => "D. immitis, three runs, the community expressed interest in including this study",
    "SRP067884"  => "T. circ, drug resistant vs suseptible. Only annotated at our competitors', nematode.net.",
    "SRP140708" => "Couldn't find featurecounts file here http://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR702/000/SRR7026110/",
    "SRP064960" => "4 runs, but requested to be included",
    "ERP023175" => "They aligned the reads to other species. Will copy from WBPS16."
);

# use Smart::Comments '###';
use JSON;
sub to_lowercase_binomial {
  my ($species) = @_;
  $species =~ s/\s+/_/;
  $species = lc($species);
  my ($spe, $cies) = split "_", $species;
  return "${spe}_${cies}";
}
sub run {
  my ($species, $assembly, $output_dir) = @_;
  $species = to_lowercase_binomial($species);
  make_path $output_dir;
  my ($selected_studies, $other_studies) =
    WbpsExpression::IncomingStudies::update_studies($studies_dir, $species, $assembly);
  #### $selected_studies
  #### $other_studies
  #### $data_files
  return unless @$selected_studies or @$other_studies;

  WbpsExpression::Analysis::run_all($selected_studies, $output_dir);

  create_listing_and_webpage($species, $selected_studies, $other_studies, $output_dir);
}

sub run_no_updates {
  my ($species, $assembly, $output_dir) = @_;
  $species = to_lowercase_binomial($species);
  make_path $output_dir;
  return unless (-d "$studies_dir/$species");
  print $species."\n";
  my @paths = map {read_dir $_ } "$studies_dir/$species";
  ### @paths
  my @our_studies = map {WbpsExpression::Study->from_folder($_, $species)} map { "$studies_dir/$species/$_" } read_dir "$studies_dir/$species";
  ### @our_studies
  my @selected_studies;
  my @other_studies;
  for my $study (@our_studies){
     print $study->{study_id}."\n";
     if ($study->{design}->is_empty or $exceptions{$study->{study_id}}){
       push @other_studies, $study;
     } else {
       $study->{quantification_method} = $study->quantification_method;
       push @selected_studies, $study;
     }
  }
  my $ref_other_studies = \@other_studies;
  my $ref_selected_studies = \@selected_studies;
  #### $selected_studies
  #### $other_studies
  #### $data_files
  return unless @$ref_selected_studies or @$ref_other_studies;

  WbpsExpression::Analysis::run_all($ref_selected_studies, $output_dir);

  create_listing_and_webpage($species, $ref_selected_studies, $ref_other_studies, $output_dir);
}

sub run_web_only {
  my ($species, $output_dir) = @_;
  $species = to_lowercase_binomial($species);
  make_path $output_dir;
my @paths = map {read_dir $_ } "$studies_dir/$species";
### @paths
  my @our_studies = map {WbpsExpression::Study->from_folder($_, $species)} map { "$studies_dir/$species/$_" } read_dir "$studies_dir/$species";
### @our_studies
  my @selected_studies;
  my @other_studies;
  for my $study (@our_studies){
     if ($study->{design}->is_empty){
       push @other_studies, $study;
     } else {
       push @selected_studies, $study;
     }
  }
  return unless @selected_studies or @other_studies;
  create_listing_and_webpage($species, \@selected_studies, \@other_studies, $output_dir);
}
sub listing_tsv {
  my ($studies) = @_;
  return map {
    die "Study without design going into the listing? $_->{study_id}" if $_->{design}->is_empty;
    join("\t", $_->{study_id}, $_->{config}{category}, $_->{config}{title})."\n"
  } @{$studies};
}
sub create_listing_and_webpage {
  my($species, $selected_studies, $other_studies, $output_dir) = @_;
  return unless @{$selected_studies} or @{$other_studies};
  write_file("$output_dir/$species.studies.tsv", { binmode => ":utf8" }, listing_tsv($selected_studies)) or die "$!";
  write_file("$output_dir/$species.studies.json", { binmode => ":utf8" }, to_json([map {$_->to_hash} @$selected_studies, @$other_studies])) or die "$!";
  write_file("$output_dir/index.html", { binmode => ":utf8" }, WbpsExpression::StudiesPage::html($species, $selected_studies, $other_studies));
}
1;
