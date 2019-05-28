use strict;
use warnings;

use Test::More;
use Test::Mock::LWP;

use WbpsExpression::IncomingStudies::RnaseqerResults;

$Mock_ua->mock(get => sub { return $Mock_response });
$Mock_response->mock( decoded_content => sub { return ""; } );
is_deeply(WbpsExpression::IncomingStudies::RnaseqerResults::get_results_by_study("other_species"), {}, "Null case");

my $species = "ancylostoma_ceylanicum";
my $study_id = "SRP035476";
my $assembly_used = "Acey_2013.11.30.genDNA";

my $out = <<"EOF";
[{
  "STUDY_ID": "$study_id",
  "SAMPLE_IDS": "SAMN02585476",
  "BIOREP_ID": "SRR1124909",
  "RUN_IDS": "SRR1124909",
  "ORGANISM": "$species",
  "REFERENCE_ORGANISM": "$species",
  "STATUS": "Complete",
  "ASSEMBLY_USED": "$assembly_used",
  "ENA_LAST_UPDATED": "Fri Jun 19 2015 18:20:11",
  "LAST_PROCESSED_DATE": "Sun Jan 21 2018 17:51:55",
  "CRAM_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909/SRR1124909.cram",
  "BEDGRAPH_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909/SRR1124909.bedgraph",
  "BIGWIG_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909/SRR1124909.bw",
  "MAPPING_QUALITY": 96
}, {
  "STUDY_ID": "$study_id",
  "SAMPLE_IDS": "SAMN02585474",
  "BIOREP_ID": "SRR1124913",
  "RUN_IDS": "SRR1124913",
  "ORGANISM": "$species",
  "REFERENCE_ORGANISM": "$species",
  "STATUS": "Complete",
  "ASSEMBLY_USED": "$assembly_used",
  "ENA_LAST_UPDATED": "Fri Jun 19 2015 18:20:11",
  "LAST_PROCESSED_DATE": "Mon Jan 22 2018 05:39:26",
  "CRAM_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/003/SRR1124913/SRR1124913.cram",
  "BEDGRAPH_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/003/SRR1124913/SRR1124913.bedgraph",
  "BIGWIG_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/003/SRR1124913/SRR1124913.bw",
  "MAPPING_QUALITY": 89
}, {
  "STUDY_ID": "$study_id",
  "SAMPLE_IDS": "SAMN02585481",
  "BIOREP_ID": "SRR1124900",
  "RUN_IDS": "SRR1124900",
  "ORGANISM": "$species",
  "REFERENCE_ORGANISM": "$species",
  "STATUS": "Complete",
  "ASSEMBLY_USED": "$assembly_used",
  "ENA_LAST_UPDATED": "Fri Jun 19 2015 18:20:10",
  "LAST_PROCESSED_DATE": "Sun Jan 21 2018 19:52:12",
  "CRAM_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/000/SRR1124900/SRR1124900.cram",
  "BEDGRAPH_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/000/SRR1124900/SRR1124900.bedgraph",
  "BIGWIG_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/000/SRR1124900/SRR1124900.bw",
  "MAPPING_QUALITY": 96
}]
EOF

$Mock_response->mock( decoded_content => sub { return $out; } );
is_deeply(
  WbpsExpression::IncomingStudies::RnaseqerResults::get_results_by_study($species),
  { $study_id => {
  assembly_used => $assembly_used,
  source_dirs_by_run => {
   SRR1124909 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909",
   SRR1124913 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/003/SRR1124913",
   SRR1124900 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/000/SRR1124900",
  },
  replicates_by_run => {SRR1124909 => "SRR1124909", "SRR1124913" => "SRR1124913", SRR1124900 => "SRR1124900"},
}},
  "Test case",
);

$Mock_response->mock( decoded_content => sub { my $x = $out; $x =~ s/"ASSEMBLY_USED": "$assembly_used"/"ASSEMBLY_USED": "pancake"/; return $x; });
is_deeply(
  WbpsExpression::IncomingStudies::RnaseqerResults::get_results_by_study($species),
  { $study_id => {
  assembly_used => "pancake",
  source_dirs_by_run => {
   SRR1124909 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909",
   SRR1124913 => "",
   SRR1124900 => "",
  },
  replicates_by_run => {SRR1124909 => "SRR1124909", "SRR1124913" => "SRR1124913", SRR1124900 => "SRR1124900"},
}}, 
  "Newest assembly",
);
$Mock_response->mock( decoded_content => sub { my $x = $out; $x =~ s/SAMN02585481/SAMN02585474/g; return $x; });
is_deeply(
  WbpsExpression::IncomingStudies::RnaseqerResults::get_results_by_study($species),
  { $study_id => {
  assembly_used => $assembly_used,
  source_dirs_by_run => {
   SRR1124909 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909",
   SRR1124913 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/003/SRR1124913",
   SRR1124900 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/000/SRR1124900",
  },
  replicates_by_run => {SRR1124909 => "SRR1124909", "SRR1124913" => "SAMN02585474", SRR1124900 => "SAMN02585474"},
}},
  "Aggregate replicates when there is more than one",
);
$Mock_response->mock( decoded_content => sub { my $x = $out; $x =~ s/SAMN\d+/SAMN0000001/g; return $x; });
is_deeply(
  WbpsExpression::IncomingStudies::RnaseqerResults::get_results_by_study($species),
  { $study_id => {
  assembly_used => $assembly_used,
  source_dirs_by_run => {
   SRR1124909 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/009/SRR1124909",
   SRR1124913 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/003/SRR1124913",
   SRR1124900 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR112/000/SRR1124900",
  },
  replicates_by_run => {SRR1124909 => "SRR1124909", "SRR1124913" => "SRR1124913", SRR1124900 => "SRR1124900"},
}},
  "Samples are unreliable replicates when they're uniform across the study",
);
done_testing;
