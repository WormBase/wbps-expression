use strict;
use warnings;
use Test::More;
use Test::MockModule;
use File::Temp qw/tempdir/;
use FindBin;

use PublicResources::Rnaseq;
use WbpsExpression;
use WbpsExpression::IncomingStudies;

my $tmp = tempdir(CLEANUP => 1);

my $species = "brugia_malayi";
my $assembly = "Bmal-4.0";

my $counts_htseq2 = "";
my $tpm_htseq2 = "";
my $bigwig = "";

subtest "run_web_only" => sub {
  WbpsExpression->new("$tmp/root_dir", "$FindBin::Bin/..")
    ->run_web_only($species, "$tmp/out");
  ok(not (-d "$tmp/root_dir"), "no data downloaded");
  ok(-d "$tmp/out", "make directory");
  ok(-s "$tmp/out/index.html", "make webpage");
  ok(-s "$tmp/out/$species.studies.tsv", "make listing");
  system("rm -rf $tmp/*");
  ok(not (-d "$tmp/out"), "cleanup directory");
};

subtest "run" => sub {
  my $public_resources = Test::MockModule->new("PublicResources::Rnaseq");
  $public_resources->mock("get" => ({
    study_id => "study_id",
    study_description_short => "study_description_short",
    study_description_full => "study_description_full",
    runs => [{
      run_id => "run_id",
      sample_id => "sample_id",
      characteristics => {"type" => "value"},
      data_files => {
        counts_htseq2 => \$counts_htseq2,
        tpm_htseq2 => \$tpm_htseq2,
        bigwig => \$bigwig,
      },
      run_description_short => "run_description_short",
      run_description_full => "run_description_full",
   }]}));

  my $subject = WbpsExpression->new("$tmp/root_dir", "$tmp/curation", "$tmp/work_dir");
  $WbpsExpression::IncomingStudies::exceptions{study_id} = "Less than usual amount of runs, still okay for testing";

  $subject->run($species, $assembly, "$tmp/out");
  
  ok(not (-d "$tmp/root_dir"), "no data downloaded - mocked out");

  ok(-d "$tmp/out", "make directory");
  ok(-s "$tmp/out/study_id/study_id.metadata_per_run.tsv", "metadata_per_run.tsv");
  ok(-s "$tmp/out/study_id/study_id.tpm_per_run.tsv", "tpm_per_run.tsv");
  ok(-s "$tmp/out/study_id/study_id.counts_per_run.tsv", "counts_per_run.tsv");

  ok(-s "$tmp/out/index.html", "make webpage");
  ok(-s "$tmp/out/$species.studies.tsv", "make listing");
  system("rm -rf $tmp/*");
  ok(not (-d "$tmp/out"), "cleanup directory");
};
done_testing;
