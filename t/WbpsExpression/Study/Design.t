use strict;
use warnings;
use Test::More;
use WbpsExpression::Study::Design;
use List::Util qw/sum pairs pairmap/;
use List::MoreUtils qw/uniq/;

sub describe_tsv {
  my ($tsv) = @_;
  my $num_colums  = 1 + scalar $tsv =~ tr/\t//;
  my $num_rows = -1 + scalar $tsv =~ tr/\n//;
  my ( $header, @lines ) = split "\n", $tsv;
  return "$num_rows rows $num_colums columns header:$header";
}

sub dimension_check_based_on_input_tsv {
  my ( $subject, $tsv, $test_name ) = @_;
  my $tabs  = scalar $tsv =~ tr/\t//;
  my $lines = scalar $tsv =~ tr/\n//;
  subtest $test_name => sub {

    my $has_replicates = grep {$_} pairmap {$a ne $b } %{ $subject->{replicates_by_run} } ;
    is(
      1 + keys %{ $subject->{replicates_by_run} },
      $lines,
"Dimension check: header + one line per run expected to be number of newlines"
    );
    is(
      $lines * ( @{ $subject->{characteristics_in_order} } + ($has_replicates ? 2 : 1)),
      $tabs,
"Dimension check: number of newlines * ( num characteristics + two or three header fields - one for no tab at EOL) expected to be number of tabs"
    );
    my $tmp = "";
    $subject->to_tsv( \$tmp );
    is( $tmp, $tsv, "from_tsv then to_tsv preserves tsv" )
      or diag explain $subject;
  };
}

sub design_remakes_itself_from_data_by_run {
  my ( $subject, $test_name ) = @_;
  my @a = (
    $subject->{replicates_by_run},
    $subject->conditions_by_run,
    $subject->characteristics_by_run,
    $subject->{characteristics_in_order}
  );
  my $subject_remade = WbpsExpression::Study::Design::from_data_by_run(@a);
  is_deeply( $subject_remade, $subject, $test_name ) or diag explain $subject_remade, $subject;
}

sub some_checks_fail {
  my ( $test_name, %checks ) = @_;
  return unless %checks;
  my $failed = sum map { not $_ } values %checks;
  ok( $failed,
    sprintf( "%s fails %s/%s checks", $test_name, $failed, scalar %checks ) );
}

sub all_checks_pass {
  my ( $test_name, @checks ) = @_;
  return unless @checks;
  subtest $test_name => sub {
    my @pairs = pairs @checks;
    plan tests => scalar @pairs;
    for my $pair (@pairs) {
      ok( $pair->[1], $pair->[0] );
    }
  };
}

sub checks_as_expected {
  my ( $expected, $test_name, @checks ) = @_;
}

sub test_design_formatting_and_checks {
  my ( $expect_checks_pass, $tsv, $test_name ) = @_;
  my $design = WbpsExpression::Study::Design::from_tsv(\$tsv);
  subtest $test_name => sub {
	dimension_check_based_on_input_tsv( $design, $tsv, "dimension_check_based_on_input_tsv" );
	design_remakes_itself_from_data_by_run( $design,
	  "design_remakes_itself_from_data_by_run" );
	if ($expect_checks_pass) {
	  all_checks_pass( "Checks pass as expected", $design->data_quality_checks );
	}
	else {
	  some_checks_fail( "Some checks fail as expected", $design->data_quality_checks );
	}
  }
}
sub design_ok {
  my ($tsv, $test_name) = @_;
  $test_name //= describe_tsv($tsv);
  test_design_formatting_and_checks(1, $tsv, $test_name);
}
sub tsv_well_formatted_but_fails_checks {
  my ($tsv, $test_name) = @_;
  $test_name //= "Formally incorrect: " . describe_tsv($tsv);
  test_design_formatting_and_checks(0, $tsv, $test_name);
}

tsv_well_formatted_but_fails_checks("Run\tCondition\n");
tsv_well_formatted_but_fails_checks("Run\tCondition\torganism part\n");
design_ok("Run\tCondition\torganism part\nSRR3209257\thead\thead\n");
design_ok("Run\tSample\tCondition\torganism part\nSRR3209257\tSRS1327198\thead\thead\n");
design_ok(
"Run\tCondition\torganism part\nSRR3209257\thead\thead\nSRR3209258\ttail\ttail\n"
);
design_ok(
"Run\tSample\tCondition\torganism part\nSRR3209257\tSRS1327198\thead\thead\nSRR3209258\tSRS1327197\ttail\ttail\n"
);
tsv_well_formatted_but_fails_checks(
  "Run\tCondition\tstrain\nSRR3209257\thead\tLE\nSRR3209258\ttail\tLE\n",
  "condition head/tail but organism part missing" );
tsv_well_formatted_but_fails_checks(
  "Run\tSample\tCondition\tstrain\nSRR3209257\tSRS1327198\thead\tLE\nSRR3209258\tSRS1327197\ttail\tLE\n",
  "condition head/tail but organism part missing" );
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
( my $tsv_no_factors = $tsv ) =~ s/head/tail/g;
tsv_well_formatted_but_fails_checks( $tsv_no_factors,
  "same as before, head -> tail ie no factors at all" );
( my $tsv_extra_condition = $tsv ) =~ s/head/2_head/;
tsv_well_formatted_but_fails_checks( $tsv_extra_condition,
"same, pointless extra condition - should fail because it makes improper slices"
);
( my $tsv_incoherent_by_run = $tsv ) =~ s/Schistosoma mansoni/other wormie/;
tsv_well_formatted_but_fails_checks( $tsv_incoherent_by_run,
  "data not assembling by condition" );

my $replicates_tsv = <<EOF;
Run	Sample	Condition	organism	developmental stage	sex	organism part
SRR3209257	head_1	head	Schistosoma mansoni	adult	female	head
SRR3209258	head_1	head	Schistosoma mansoni	adult	female	head
SRR3209259	head_2	head	Schistosoma mansoni	adult	female	head
SRR3209260	tail_1	tail	Schistosoma mansoni	adult	female	tail
SRR3209261	tail_2	tail	Schistosoma mansoni	adult	female	tail
SRR3209262	tail_3	tail	Schistosoma mansoni	adult	female	tail
EOF
design_ok($replicates_tsv);
( my $replicates_tsv_incoherent_by_run = $replicates_tsv ) =~ s/Schistosoma mansoni/other wormie/;
tsv_well_formatted_but_fails_checks( $replicates_tsv_incoherent_by_run,
  "data not assembling by replicate" );
( my $replicates_tsv_conditions_not_assembling_by_replicate = $replicates_tsv ) =~ s/tail_1/head_2/;
tsv_well_formatted_but_fails_checks($replicates_tsv_conditions_not_assembling_by_replicate,
  "conditions not grouping by replicate" );

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
done_testing;
