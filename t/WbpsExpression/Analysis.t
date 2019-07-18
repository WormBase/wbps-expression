use strict;
use warnings;
use Test::More;
use File::Temp qw/tempdir/;
use File::Path qw/make_path/;
use File::Slurp qw/read_file write_file/;

use WbpsExpression::Analysis;
use WbpsExpression::Study;
use WbpsExpression::Study::Design;

my $tmp = tempdir(CLEANUP => 1);

make_path "$tmp/in/SRR3209257";

write_file "$tmp/in/SRR3209257/SRR3209257.pe.genes.raw.htseq2.tsv", <<"EOF";
mk4.000000.02	0
mk4.000000.03	207
mk4.000000.08	0
EOF

write_file "$tmp/in/SRR3209257/SRR3209257.pe.genes.tpm.htseq2.irap.tsv", <<"EOF";
Gene	V2
mk4.000000.02	26.95
mk4.000000.03	21.14
mk4.000000.08	0
EOF

my $study = bless {
  study_id => "study_id",
  design => WbpsExpression::Study::Design::from_tsv(\"Run\tCondition\torganism part\nSRR3209257\thead\thead\n"),
  config => {
    title => "config.title",
    submitting_centre => "config.submitting_centre",
    pubmed => {},
    rnaseqer_last_update => "2017-01-01",
    contrasts => [],
  },
  sources => {
	SRR3209257 => {
	  location => "$tmp/in/SRR3209257",
      end => "pe",
      quality => 30,
	}
  },
}, 'WbpsExpression::Study';

WbpsExpression::Analysis::run("$tmp/out", $study);

my $expected = <<"EOF";
gene_id	SRR3209257
mk4.000000.02	0
mk4.000000.03	207
mk4.000000.08	0
---
Run	Condition	organism part
SRR3209257	head	head
---
# 
# Study study_id: config.title
# Submitted to archives by config.submitting_centre
# Reads aligned with TopHat2 and quantified with HTSeq, by RNASeq-er: https://www.ebi.ac.uk/fg/rnaseq/api
# 
# Values are transcripts per million units (TPMs) per gene for each run
# !SRR3209257: Low mapping quality: 30
gene_id	!SRR3209257
mk4.000000.02	26.95
mk4.000000.03	21.14
mk4.000000.08	0
EOF
my ($counts_expected, $metadata_expected, $tpms_expected) = map {s/^\s*//; $_} split "---", $expected;
my $counts = read_file("$tmp/out/study_id/study_id.counts_per_run.tsv");
my $metadata = read_file("$tmp/out/study_id/study_id.metadata_per_run.tsv");
my $tpms = read_file("$tmp/out/study_id/study_id.tpm_per_run.tsv");
is_deeply($counts, $counts_expected, "counts");
is_deeply($metadata, $metadata_expected , "metadata");
is_deeply($tpms, $tpms_expected, "tpms");

done_testing;
