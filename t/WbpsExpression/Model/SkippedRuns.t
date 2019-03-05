use strict;
use warnings;
use Test::More;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;

use WbpsExpression::Model::SkippedRuns;
sub test_create {
  my ($runs, $config, $test_name) = @_;
  my $runs_path = join("\n", "Run", @$runs);
  my $tmpdir = tempdir(CLEANUP => 1);
  my $study_id = basename $tmpdir;
  my $subject = WbpsExpression::Model::SkippedRuns->from_paths($study_id, \$runs_path, $config);
  is_deeply($subject->{config}, $config, "$test_name - preserves config");
  is_deeply($subject->{runs}, $runs, "$test_name - preserves runs");
  $subject->to_folder($tmpdir);
  is_deeply(WbpsExpression::Model::SkippedRuns->from_folder($tmpdir), $subject, "$test_name - to and from folder") or diag `ls $tmpdir`;
}

test_create([], {}, "Null case");
test_create(["run_id"], {}, "One run");
test_create([], {description => "study description"},  "One property");

done_testing;