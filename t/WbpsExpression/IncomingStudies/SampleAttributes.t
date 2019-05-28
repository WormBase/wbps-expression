use strict;
use warnings;

use Test::More;
use WbpsExpression::IncomingStudies::SampleAttributes;

my $species = "schistosoma mansoni";

sub test_case {
  my ($input, $expected, $test_name) = @_;
  my $actual = WbpsExpression::IncomingStudies::SampleAttributes::standardise_for_worm(
     $species, $expected
  );
  is_deeply($actual, $expected, $test_name);
}

sub unchanged {
  my ($o, $test_name) = @_;
  test_case($o, $o, $test_name);
}
unchanged({}, "null case");
unchanged({unrelated_type => "unrelated_value"}, "domain agnostic case");

sub reject_sample_name {
  my ($sample_name) = @_;
  test_case({sample_name => $sample_name}, {}, "Reject $sample_name");
}

reject_sample_name("Sample from $species");
reject_sample_name("DRS026763");

test_case({dev_stage => "adult females"}, {developmental_stage => "adult", sex => "female"}, "Standardise dev stage");
done_testing;
