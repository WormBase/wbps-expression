
use EnsEMBL::Web::Component::Gene::WbpsExpression;
use Test::More skip_all => "Not implemented yet";
use File::Temp qw/tempdir/;
use File::Slurp qw/write_file/;

my $dir = tempdir(CLEANUP => 1);
my $f1 = <<EOF
#temp file
	heads	tails
g1	1.1	1.2
g2	2.1	2.1
EOF
my $study_id = "SRP071241";
my $f1_name = "$study_id.tpm.tsv";
write_file(join("/", $dir, $study_id, $f1_name), $f1);

my $subject = EnsEMBL::Web::Component::Gene::WbpsExpression::from_folder(
   $dir
);

is_deeply($subject, bless({
  studies => [{
     accession => "$study_id",
     title => "Comparison of gene expression between female Schistosoma mansoni heads and tails",
     tpms_per_condition => join("/", $study_id, $f1_name),
  }]
}, 'EnsEMBL::Web::Component::Gene::WbpsExpression'), "Create reads in the config");


done_testing;
