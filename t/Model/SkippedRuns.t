use strict;
use warnings;
use Test::More;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;

use Model::SkippedRuns;
sub skipped_runs_formally_correct {
  my ($runs, $config, $test_name) = @_;
  my $runs_path = join("\n", "Run", @$runs);
  my $tmpdir = tempdir(CLEANUP => 1);
  my $study_id = basename $tmpdir;
  my $subject = Model::SkippedRuns->from_paths($study_id, \$runs_path, $config);
  is_deeply($subject->{config}, $config, "$test_name - preserves config");
  is_deeply($subject->{runs}, $runs, "$test_name - preserves runs");
  $subject->to_folder($tmpdir);
  is_deeply(Model::SkippedRuns->from_folder($tmpdir), $subject, "$test_name - to and from folder") or diag `ls $tmpdir`;
}

test_create([], {}, "Null case");
test_create(["run id"], {}, "One run");
test_create([], {description => "study description"},  "One property");

done_testing;
