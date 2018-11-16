use Test::More;
use Model::Design;
# Test utility code here
# More tests for designs are bundled in t/Model/curation.t
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
