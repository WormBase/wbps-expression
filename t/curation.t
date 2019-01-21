use strict;
use warnings;
use Test::More;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use File::Find;
use File::Slurp qw/read_file/;
use FindBin;
use Model::Design;
use Model::Study;
use Model::SkippedRuns;
use List::Util qw/sum pairs pairmap/;
use List::MoreUtils qw/uniq/;

sub describe_tsv {
  my ( $header, @lines ) = split "\n", shift;
  return sprintf( "%s + %s lines", $header, 0 + @lines );
}

sub design_dumps_to_tsv {
  my ( $subject, $tsv, $test_name ) = @_;
  $test_name //= describe_tsv($tsv);
  my $tabs  = $tsv =~ tr/\t//;
  my $lines = $tsv =~ tr/\n//;
  subtest $test_name => sub {

    my $has_replicates = grep {$_} pairmap {$a ne $b } %{ $subject->{replicates_per_run} } ;
    is(
      1 + keys %{ $subject->{replicates_per_run} },
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
    $subject->{replicates_per_run},
    $subject->conditions_per_run,
    $subject->characteristics_per_run,
    $subject->{characteristics_in_order}
  );
  my $subject_remade = Model::Design::from_data_by_run(@a);
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
  if ($expected) {
    all_checks_pass( "expect pass checks: $test_name", @checks );
  }
  else {
    some_checks_fail( "expect fail some checks: $test_name", @checks );
  }
}

sub test_design_formatting_and_checks {
  my ( $expect_checks_pass, $design, $tsv, $test_name ) = @_;
  design_dumps_to_tsv( $design, $tsv, "design_dumps_to_tsv: $test_name" );
  design_remakes_itself_from_data_by_run( $design,
    "design_remakes_itself_from_data_by_run: $test_name" );
  checks_as_expected(
    $expect_checks_pass,
    "design in $test_name",
    $design->data_quality_checks
  );
}

sub test_study {
  my ( $subject, $tsv, $test_name, $design_checks_should_pass,
    $config_checks_should_pass )
    = @_;
  $test_name //= describe_tsv($tsv);
  test_design_formatting_and_checks( $design_checks_should_pass,
    $subject->{design}, $tsv, $test_name );
  checks_as_expected(
    $config_checks_should_pass,
    "config in $test_name",
    $subject->consistency_checks
  );
}

sub test_tsv_and_config {
  my ( $tsv, $config, $test_name, $design_checks_should_pass,
    $config_checks_should_pass )
    = @_;
  test_study( Model::Study->from_paths( "study_id", \$tsv, $config ),
    $tsv, $test_name, $design_checks_should_pass, $config_checks_should_pass, );
}

sub test_folder {
  my ($path) = @_;
  my $study_id = basename $path;
  ok( -f "$path/$study_id.tsv",  "$path design exists" );
  ok( -f "$path/$study_id.yaml", "$path config exists" );
  my $tsv = read_file "$path/$study_id.tsv";
  test_study(
    Model::Study->from_folder($path),
    $tsv, "test_folder $path",
    1, 1
  );
}

sub test_tsv_with_empty_config {
  test_tsv_and_config( shift, {}, shift, shift, 1 );
}

sub tsv_well_formatted_but_fails_checks {
  test_tsv_with_empty_config( shift, shift, 0 );
}

sub design_ok {
  test_tsv_with_empty_config( shift, shift, 1 );
}

sub ok_design_fails_with_wrong_condition_names {
  test_tsv_and_config( shift, { condition_names => { wrong_condition => "" } },
    shift, 1, 0 );
}

sub ok_design_fails_with_wrong_contrasts {
  test_tsv_and_config( shift, { contrasts => [{ name => "contrasts.name", values => [["wrong_reference", "wrong_test", "wrong_reference vs wrong_test : contrast name"]]}] } ,
    shift, 1, 0 );
}

sub skipped_runs_ok {
  my ($skipped_runs_path, $study_path) = @_;
  my $subject = Model::SkippedRuns->from_folder($skipped_runs_path);
  subtest "$skipped_runs_path" => sub {
    cmp_ok (scalar @{$subject->{runs}}, '>', 0 , "Some runs");
    if($study_path){
      ok (not (keys %{$subject->{config}}), "Study config should go with curated runs");
      my @curated_runs = Model::Study->from_folder($study_path)->{design}->all_runs;
      my @overlap = grep {my $x = $_; grep {$x eq $_ } @curated_runs } @{$subject->{runs}};
      is_deeply(\@overlap, [], "Skipped and curated runs shouldn't overlap");
    } else {
      ok (scalar keys %{$subject->{config}}, "Study skipped whole - should have config with study description"); 
    }
  }
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
ok_design_fails_with_wrong_condition_names($tsv);
ok_design_fails_with_wrong_contrasts($tsv);
( my $tsv_no_factors = $tsv ) =~ s/head/tail/g;
tsv_well_formatted_but_fails_checks( $tsv_no_factors,
  "same as before, head -> tail ie no factors at all" );
( my $tsv_extra_condition = $tsv ) =~ s/head/2_head/;
tsv_well_formatted_but_fails_checks( $tsv_extra_condition,
"same, pointless extra condition - should fail because it makes improper slices"
);
( my $tsv_incoherent_per_run = $tsv ) =~ s/Schistosoma mansoni/other wormie/;
tsv_well_formatted_but_fails_checks( $tsv_incoherent_per_run,
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
ok_design_fails_with_wrong_condition_names($replicates_tsv);
ok_design_fails_with_wrong_contrasts($replicates_tsv);
( my $replicates_tsv_incoherent_per_run = $replicates_tsv ) =~ s/Schistosoma mansoni/other wormie/;
tsv_well_formatted_but_fails_checks( $replicates_tsv_incoherent_per_run,
  "data not assembling by replicate" );
( my $replicates_tsv_conditions_not_assembling_by_replicate = $replicates_tsv ) =~ s/tail_1/head_2/;
eval {
  Model::Design::from_tsv($replicates_tsv_conditions_not_assembling_by_replicate)
};
like($@, "/head_2/", "conditions not grouping by replicate - death");

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
ok_design_fails_with_wrong_condition_names($two_factor_tsv);
ok_design_fails_with_wrong_contrasts($two_factor_tsv);
my %studies;
my $arg = $ARGV[0] // "";
my $study_folder_pattern = qr/$arg.*\w+\d+$/;
find(
  sub {
    $studies{basename $File::Find::name} = $File::Find::name if -d $File::Find::name and $File::Find::name =~ $study_folder_pattern;
  },
  "$FindBin::Bin/../curation/studies"
);
my %skipped_runs;
find(
  sub {
    $skipped_runs{basename $File::Find::name} = $File::Find::name if -d $File::Find::name and $File::Find::name =~ $study_folder_pattern;
  },
  "$FindBin::Bin/../curation/skipped_runs"
);

test_folder($_) for values %studies;
for my $study_id (keys %skipped_runs){
  skipped_runs_ok($skipped_runs{$study_id}, $studies{$study_id}); 
}


done_testing;
