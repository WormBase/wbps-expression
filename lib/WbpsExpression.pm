use strict;
use warnings;
package WbpsExpression;
use WbpsExpression::IncomingStudies;
use WbpsExpression::Analysis;
use WbpsExpression::StudiesPage;
use File::Slurp qw/write_file/;
use File::Path qw/make_path/;
use Log::Any '$log';
# use Smart::Comments '###';
sub new {
  my ($class, $root_dir, $src_dir) = @_;

  my $sheets = WbpsExpression::IncomingStudies::Sheets->new($src_dir);
  my $public_rnaseq_studies = PublicResources::Rnaseq->new($root_dir, $sheets);

  return bless {
    sheets => $sheets,
    incoming_studies => WbpsExpression::IncomingStudies->new($sheets, $public_rnaseq_studies),
  }, $class;
}
sub to_lowercase_binomial {
  my ($species) = @_;
  $species =~ s/\s+/_/;
  $species = lc($species);
  my ($spe, $cies) = split "_", $species;
  return "${spe}_${cies}";
}
sub run {
  my ($self, $species, $assembly, $output_dir) = @_;
  $species = to_lowercase_binomial($species);
  make_path $output_dir;
  my ($selected_studies, $other_studies, $data_files) =$self->{incoming_studies}->fetch_all($species, $assembly);
  #### $selected_studies
  #### $other_studies
  #### $data_files
  return unless @$selected_studies or @$other_studies;

  my @selected_studies_passing_checks;
  for my $study (@$selected_studies){
    if ($study->passes_checks){
      push @selected_studies_passing_checks, $study;
    } else {
      $log->info(sprintf(__PACKAGE__ . " Study %s failing checks - see %s",$study->{study_id}, $self->{sheets}->path("studies", $study->{study_id})));
    }
  }
 
  WbpsExpression::Analysis::run_all(\@selected_studies_passing_checks, $data_files, $output_dir);

  create_listing_and_webpage($species, $selected_studies, $other_studies, $output_dir);
}

sub run_web_only {
  my ($self, $species, $output_dir) = @_;
  $species = to_lowercase_binomial($species);
  make_path $output_dir;
  my ($selected_studies, undef, $other_studies) = $self->{sheets}->read_directories($species);
  return unless @$selected_studies or @$other_studies;
  create_listing_and_webpage($species, $selected_studies, $other_studies, $output_dir);
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
  write_file("$output_dir/index.html", { binmode => ":utf8" }, WbpsExpression::StudiesPage::html($species, $selected_studies, $other_studies));
}
1;
