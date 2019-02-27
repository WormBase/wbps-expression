use strict;
use warnings;
use EnsEMBL::Web::Component::Gene::WbpsExpression;
use Test::More;
use File::Temp qw/tempdir/;
use File::Slurp qw/write_file/;
use File::Path qw/make_path/;

my $dir = tempdir(CLEANUP => 1);
my $file = <<EOF;
#temp file
	heads	tails
g1	1.1	1.2
g2	2.1	2.2
EOF
my $study_id = "SRP071241";
my $category = "Organism parts";
my $study_title = "Comparison of gene expression between female Schistosoma mansoni heads and tails";
my $file_name = "$study_id.tpm.tsv";
my $file_path = join("/", $dir, $study_id, $file_name);

my $species = "schistosoma_mansoni_prjea36577";
make_path(join("/", $dir, $study_id));
write_file($file_path, $file);

my $second_file = <<EOF;
# DESeq2 version: ‘1.22.1’
	5-AzaC vs untreated
g1	1.1 0.04
g3	-2.3 1.2e-3
EOF

my $second_study_id ="SRP130864";
my $second_category = "Response to treatment";
my $second_study_title = "5-AzaC effect on Schistosoma mansoni Transcriptome";
my $second_file_name = "$second_study_id.de.treatment.tsv";
my $second_file_path = join("/", $dir, $second_study_id, $second_file_name);
make_path(join("/", $dir, $second_study_id));
write_file($second_file_path, $second_file);

my $third_file = <<EOF;
	SRR5664530	SRR5664533	SRR5664534	SRR5664531	SRR5664532	SRR5664529	SRR5664535
g1	0.8	1.3	3.54	35.85	41.03	41.95	48.82
g2	11.61	11.32	13.32	120.31	128.8	140.44	148.89
EOF

my $third_study_id ="SRP108901";
my $third_category = "Other";
my $third_study_title = "Schistosoma mansoni strain LE - Transcriptome or Gene expression";
my $third_file_name = "$third_study_id.tpm_per_run.tsv";
my $third_file_path = join("/", $dir, $third_study_id, $third_file_name);
make_path(join("/", $dir, $third_study_id));
write_file($third_file_path, $third_file);

my $studies_file = <<"EOF";
$study_id\t$category\t$study_title
$second_study_id\t$second_category\t$second_study_title
$third_study_id\t$third_category\t$third_study_title
EOF

write_file(join("/", $dir, "$species.studies.tsv"), $studies_file);

my $subject = EnsEMBL::Web::Component::Gene::WbpsExpression::from_folder(
   $species, $dir
);
is( scalar @{$subject->{studies}}, 3, "Read in three studies");

sub is_empty_response {
  my ($payload, $test_name) = @_;
  subtest $test_name => sub {
     like($payload, qr{<html>.*</html>}, "is html");
     like($payload, qr{no results}i,  "no results");
  };
}
sub is_tables {
  my ($payload, $num_tables, $num_rows_per_table, $num_columns_per_row, $test_name) = @_;
  subtest $test_name => sub {
    my @table_tags = $payload =~ /<table>/g;
    my @th_tags = $payload =~ /<th>/g;
    my @tr_tags = $payload =~ /<tr>/g;
    my @td_tags = $payload =~ /<td>/g;
    is(scalar @table_tags, $num_tables, "num tables");
    is(scalar @th_tags, $num_tables, "each table has header");
    is(scalar @tr_tags, $num_tables * ($num_rows_per_table) , "num rows per table");
    is(scalar @td_tags, $num_tables * ($num_rows_per_table +1) * $num_columns_per_row, "num columns per row");
  } or diag explain $payload;
}

is_empty_response($subject->render_page("g0", $_), "Invalid gene, category $_") for ($category, $second_category, $third_category, "invalid category");
is_empty_response($subject->render_page("g1", "invalid category"), "g1 invalid category");
is_tables($subject->render_page("g1", $category), 1, 1, 2, "$category - one row in a flat table");
is_tables($subject->render_page("g1", $second_category), 1, 1, 5, "$second_category ");
is_tables($subject->render_page("g1", $third_category), 1, 1, 6, "$third_category - six stats");
done_testing;
