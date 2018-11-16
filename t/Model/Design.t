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
   my $tabs = $tsv =~ tr/\t//;
   my $lines = $tsv =~ tr/\n//;
   my $subject = Model::Design::from_tsv(\$tsv);
   subtest $test_name => sub {
	 is(%{$subject->{conditions_per_run}}+1, $lines,
       "Dimension check: header + one line per run expected to be number of newlines"
     );
	 is($lines * (@{$subject->{characteristics_in_order}}+1), $tabs,
       "Dimension check: number of newlines * ( num characteristics + two header fields - one for no tab at EOL) expected to be number of tabs"
     );
	 my $tmp = "";
	 $subject->to_tsv(\$tmp);
	 is($tmp, $tsv, "from_tsv then to_tsv preserves tsv") or diag explain $subject;
	 my @a = ($subject->{conditions_per_run}, $subject->characteristics_per_run, $subject->{characteristics_in_order});
	 
	 is_deeply(Model::Design::from_data_by_run(@a), $subject, "characteristics_per_run stable") or diag explain \@a;
   }
}
sub fails_data_quality_checks {
   my ($tsv, $test_name) = @_;
   $test_name //= "fails_data_quality_checks: ".describe_tsv($tsv);
   my %checks = Model::Design::from_tsv(\$tsv)->data_quality_checks;
   my $failed = sum map {not $_} values %checks;
   ok($failed, sprintf("%s inadequate data quality: caught by %s/%s checks", $test_name, $failed,scalar %checks));
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
design_ok("Run\tCondition\torganism part\nSRR3209257\thead\thead\nSRR3209258\ttail\ttail\n");
well_formatted_but_fails_checks("Run\tCondition\tstrain\nSRR3209257\thead\tLE\nSRR3209258\ttail\tLE\n", "condition head/tail but organism part missing");

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
well_formatted_but_fails_checks($tsv_extra_condition, "same, pointless extra condition - should fail because it makes improper slices");
(my $tsv_incoherent_per_run = $tsv) =~ s/Schistosoma mansoni/other wormie/;
well_formatted_but_fails_checks($tsv_incoherent_per_run, "data not assembling by condition");


my $two_factor_tsv = <<EOF;
Run	Condition	organism	developmental stage	sex	organism part
SRR3209257	hf	Schistosoma mansoni	adult	female	head
SRR3209258	hf	Schistosoma mansoni	adult	female	head
SRR3209259	hf	Schistosoma mansoni	adult	female	head
SRR0000007	hm	Schistosoma mansoni	adult	male	head
SRR0000008	hm	Schistosoma mansoni	adult	male	head
SRR0000009	hm	Schistosoma mansoni	adult	male	head
SRR3209260	tf	Schistosoma mansoni	adult	female	tail
SRR3209261	tf	Schistosoma mansoni	adult	female	tail
SRR3209262	tf	Schistosoma mansoni	adult	female	tail
SRR0000000	tm	Schistosoma mansoni	adult	male	tail
SRR0000001	tm	Schistosoma mansoni	adult	male	tail
SRR0000002	tm	Schistosoma mansoni	adult	male	tail
EOF
design_ok($two_factor_tsv);

sub test_choose {
  my ($n, $expected) = @_;
  my $actual = [Model::Design::n_choose_one(0..$n)];
  is_deeply($actual, $expected, "n_choose_one 0..$n") or diag explain $actual;
}
test_choose(-1, []);
test_choose(0, [[0,[]]]);
test_choose(1, [[0,[1]], [1, [0]]]);
test_choose(2, [[0,[1,2]], [1,[0,2]], [2,[0,1]]]);

sub test_subsets {
  my @actual = Model::Design::subsets(@_);
  ok(@actual == 2**@_, "Subsets: @_") or diag explain \@actual;
}
test_subsets();
test_subsets(1);
test_subsets("apple", "banana");
test_subsets("apple", "banana", "pear");
done_testing;
