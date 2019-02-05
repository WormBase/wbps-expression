use strict;
use warnings;
use EnsEMBL::Web::Component::Gene::WbpsExpression;
use Test::More;
use File::Temp qw/tempdir/;
use File::Slurp qw/write_file/;
use File::Path qw/make_path/;

my $dir = tempdir(CLEANUP => 1);
my $f1 = <<EOF;
#temp file
	heads	tails
g1	1.1	1.2
g2	2.1	2.2
EOF
my $study_id = "SRP071241";
my $study_title = "Comparison of gene expression between female Schistosoma mansoni heads and tails";
my $f1_name = "$study_id.tpm.tsv";
my $f1_path = join("/", $dir, $study_id, $f1_name);

my $species = "schistosoma_mansoni";
my $assembly = "Smansoni_v7";
make_path(join("/", $dir, $study_id));
write_file($f1_path, $f1);
write_file(join("/", $dir, "$species.$assembly.studies.tsv"), "$study_id\t$study_title\n");

my $subject = EnsEMBL::Web::Component::Gene::WbpsExpression::from_folder(
   $species, $assembly, $dir
);

is_deeply($subject, bless({
  studies => [{
     study_id => "$study_id",
     study_title => $study_title,
     tpms_per_condition => $f1_path, 
  }]
}, 'EnsEMBL::Web::Component::Gene::WbpsExpression'), "Create reads in the config");

is_deeply($subject->get_data("invalid ID"), [], "Null case");

is_deeply($subject->get_data("g2"), [{
  study_id => $study_id,
  study_title => $study_title,
  conditions => ["heads", "tails"],
  expression_tpm => [2.1, 2.2],
}] , "One line");


done_testing;
