use Test::More;
use WbpsExpression::IncomingStudies::CurationDefaults;
use WbpsExpression::Study::Design;
sub test_subsets {
  my @actual = WbpsExpression::IncomingStudies::CurationDefaults::subsets(@_);
  ok(@actual == 2**@_, "Subsets: @_") or diag explain \@actual;
}
test_subsets();
test_subsets(1);
test_subsets("apple", "banana");
test_subsets("apple", "banana", "pear");

sub contrasts_as_expected {
  my ($design_tsv, $expected, $test_name) = @_;
  my $actual = WbpsExpression::IncomingStudies::CurationDefaults::contrasts(WbpsExpression::Study::Design::from_tsv(\$design_tsv));
  is_deeply($actual, $expected, $test_name) or warn explain $actual, $expected;
}
contrasts_as_expected("Run\tCondition\n", [], "null case");
my $tsv = <<'EOF';
Run	Condition	developmental_stage
r11	head	head
r12	head	head
r13	head	head
r21	tail	tail
r22	tail	tail
r23	tail	tail
EOF
contrasts_as_expected($tsv, [{name => "developmental_stage", values => [["head", "tail", "head vs tail"]]}], "minimal example");
(my $no_characteristics_tsv = $tsv) =~ s/head\thead/head\ttail/;
contrasts_as_expected($no_characteristics_tsv, [{name => "", values => [["head", "tail", "head vs tail"]]}], "no characteristics needed");
(my $no_factors_tsv = $tsv) =~ s/head/tail/g;
contrasts_as_expected($no_factors_tsv, [], "no factors example");
my $two_factor_tsv = <<'EOF';
Run	Condition	sax	developmental_stage
r11	hf	female	head
r12	hf	female	head
r13	hf	female	head
r21	tf	female	tail
r22	tf	female	tail
r23	tf	female	tail
r31	hm	male	head
r32	hm	male	head
r33	hm	male	head
r41	tm	male	tail
r42	tm	male	tail
r43	tm	male	tail
EOF
contrasts_as_expected($two_factor_tsv, [
  {
    'name' => 'developmental_stage',
    'values' => [
      [
        'hf',
        'tf',
        'hf vs tf'
      ],
      [
        'hm',
        'tm',
        'hm vs tm'
      ]
    ]
  },
  {
    'name' => 'sax',
    'values' => [
      [
        'hf',
        'hm',
        'hf vs hm'
      ],
      [
        'tf',
        'tm',
        'tf vs tm'
      ]
    ]
  }
], "two characteristics example");

is(WbpsExpression::IncomingStudies::CurationDefaults::trim_values("3rd stage dispersal juvenile, pooled male and female", "3rd stage dispersal juvenile", "pooled male and female"), "", "trim values");

done_testing;
