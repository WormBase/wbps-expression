use strict;
use warnings;
package WbpsExpression;
use WbpsExpression::IncomingStudies;
use WbpsExpression::Analysis;
use WbpsExpression::StudiesPage;
use File::Slurp qw/write_file read_dir/;
use File::Path qw/make_path/;
use FindBin;
my $studies_dir = "$FindBin::Bin/../studies";
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
sub run_web_only {
  my ($species, $output_dir) = @_;
  $species = to_lowercase_binomial($species);
  make_path $output_dir;
my @paths = map {read_dir $_ } "$studies_dir/$species";
### @paths
  my @our_studies = map {WbpsExpression::Study->from_folder($_)} map { "$studies_dir/$species/$_" } read_dir "$studies_dir/$species";
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
    join("\t", $_->{study_id}, $_->{config}{category}, $_->{config}{title})."\n"
  } @{$studies};
}
sub create_listing_and_webpage {
  my($species, $selected_studies, $other_studies, $output_dir) = @_;
  return unless @{$selected_studies} or @{$other_studies};
  write_file("$output_dir/$species.studies.tsv", { binmode => ":utf8" }, listing_tsv($selected_studies)) or die "$!";
  write_file("$output_dir/$species.studies.json", { binmode => ":utf8" }, JSON->new->encode([map {$_->to_hash} @$selected_studies, @$other_studies])) or die "$!";
  write_file("$output_dir/index.html", { binmode => ":utf8" }, WbpsExpression::StudiesPage::html($species, $selected_studies, $other_studies));
}
1;
