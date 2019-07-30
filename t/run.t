use strict;
use warnings;
use Test::More;
use File::Temp qw/tempdir/;
use JSON;
use File::Slurp qw/read_file/;

use WbpsExpression;

my $tmp = tempdir(CLEANUP => 1);

my $species = "brugia_malayi";
my $assembly = "Bmal-4.0";

my $counts_htseq2 = "";
my $tpm_htseq2 = "";
my $bigwig = "";

subtest "run_web_only" => sub {
  WbpsExpression::run_web_only($species, "$tmp/out");
  ok(-d "$tmp/out", "make directory");
  ok(-s "$tmp/out/index.html", "make webpage");
  ok(-s "$tmp/out/$species.studies.tsv", "make listing");
  ok(-s "$tmp/out/$species.studies.json", "make JSON for tracks");
  my @studies = @{from_json(read_file("$tmp/out/$species.studies.json", {binmode => ":utf8"}))};
  ok (scalar @studies, "Some studies");
  system("rm -rf $tmp/*");
  ok(not (-d "$tmp/out"), "cleanup directory");
};

done_testing;
