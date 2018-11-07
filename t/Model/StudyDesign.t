
use Test::More;
use Test::MockModule;
use File::Temp qw/tempdir/;
use Model::StudyDesign;

sub test_preserve_format {
   my ($tsv, $test_name) = @_;
   unless($test_name){
      my ($header, @lines) = split "\n", $tsv;
      $test_name = sprintf("%s + %s lines", $header , 0+@lines);
   }
   my $subject = Model::StudyDesign::from_tsv(\$tsv);
   my $tmp = "";
   $subject->to_tsv(\$tmp);
   is($tmp, $tsv, $test_name) or diag explain $subject;
}

test_preserve_format("Condition\n");
test_preserve_format("Condition\torganism part\n");
my $tsv = <<EOF;
Condition	organism	developmental stage	sex	organism part
head	Schistosoma mansoni	adult	female	head
tail	Schistosoma mansoni	adult	female	tail
EOF
test_preserve_format($tsv);
done_testing;
