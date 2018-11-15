use Test::More;
use Test::MockModule;
use File::Temp qw/tempdir/;
use Model::Design;
use List::Util qw/sum pairs/;

sub describe_tsv {
  my ($header, @lines) = split "\n", shift;
  return sprintf("%s + %s lines", $header , 0+@lines);
}
sub test_preserve_format {
   my ($tsv, $test_name) = @_;
   $test_name //= describe_tsv($tsv);
   my $subject = Model::Design::from_tsv(\$tsv);
   my $tmp = "";
   $subject->to_tsv(\$tmp);
   is($tmp, $tsv, "from_tsv then to_tsv preserves tsv: $test_name") or diag explain $subject;
   my @a = ($subject->{conditions_per_run}, $subject->characteristics_per_run, $subject->{characteristics_in_order});
   
   is_deeply(Model::Design::from_data_by_run(@a), $subject, "characteristics_per_run stable: $test_name") or diag explain @a;
}
sub fails_data_quality_checks {
   my ($tsv, $test_name) = @_;
   $test_name //= "fails_data_quality_checks: ".describe_tsv($tsv);
   my %checks = Model::Design::from_tsv(\$tsv)->data_quality_checks;
   my $failed = sum map {not $_} values %checks;
   ok($failed, sprintf("Failed %s/%s checks: %s", $failed,scalar %checks, $test_name));
}
sub passes_data_quality_checks {
   my ($tsv, $test_name) = @_;
   $test_name //= "passes_data_quality_checks: ".describe_tsv($tsv);
   subtest $test_name => sub {
     my @pairs = pairs Model::Design::from_tsv(\$tsv)->data_quality_checks;
     plan tests => scalar @pairs;
     for my $pair (@pairs) {
       ok($pair->value, $pair->key);
     }
   }
}
sub well_formatted_but_fails_checks {
   test_preserve_format(@_);
   fails_data_quality_checks(@_);
}
sub design_ok {
   test_preserve_format(@_);
   passes_data_quality_checks(@_);
}

my $in = -t STDIN ? "" : do {local $/; <>};
if($in){
   design_ok($in, "test stdin");
   done_testing;
   exit;
}
well_formatted_but_fails_checks("Run\tCondition\n");
well_formatted_but_fails_checks("Run\tCondition\torganism part\n");
design_ok("Run\tCondition\torganism part\nSRR3209257\thead\thead\n");
my $tsv = <<EOF;
Run	Condition	organism	developmental stage	sex	organism part
SRR3209257	head	Schistosoma mansoni	adult	female	head
SRR3209258	head	Schistosoma mansoni	adult	female	head
SRR3209259	head	Schistosoma mansoni	adult	female	head
SRR3209260	tail	Schistosoma mansoni	adult	female	tail
SRR3209261	tail	Schistosoma mansoni	adult	female	tail
SRR3209262	tail	Schistosoma mansoni	adult	female	tail
EOF
design_ok($tsv);
(my $tsv_no_factors = $tsv) =~ s/head/tail/g;
design_ok($tsv_no_factors, "same as before, head -> tail ie no factors at all");
(my $tsv_extra_condition = $tsv) =~ s/head/2_head/;
design_ok($tsv_extra_condition, "same, pointless extra condition\n$tsv_extra_condition\n");
(my $tsv_incoherent_per_run = $tsv) =~ s/Schistosoma mansoni/other wormie/;
well_formatted_but_fails_checks($tsv_incoherent_per_run, "data not assembling by condition\n$tsv_incoherent_per_run\n");
done_testing;
