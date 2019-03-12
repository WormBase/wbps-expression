use strict;
use warnings;
use Test::More;
use WbpsExpression::Analysis::QualityWarnings;
use WbpsExpression::Model::Design;
use List::MoreUtils qw/zip/;
use List::Util qw/pairs/;

sub test_warnings_per_condition {
  my ($design_tsv, $qc_issues_by_run, $expected_num_issues, $test_name) = @_;
  
  my $design = WbpsExpression::Model::Design::from_tsv(\$design_tsv);
  my $conditions_ordered = [$design->all_conditions];

  my ($conditions_amended_names, $warnings) = 
     WbpsExpression::Analysis::QualityWarnings::conditions_amended_names_and_warnings_for_per_condition_analysis(
     $design, $qc_issues_by_run, $conditions_ordered,
  );
  subtest $test_name => sub {
    is_deeply([sort keys %$conditions_amended_names], $conditions_ordered, "All conditions as keys") or diag explain $conditions_amended_names;
    my $conditions_changed_names = grep {$_->[0] ne $_->[1]} pairs %$conditions_amended_names;
    is ($conditions_changed_names, $expected_num_issues, "Each expected issue is a changed condition name") or diag explain $conditions_amended_names;
    is (scalar @{$warnings}, $expected_num_issues, "Each expected issue is a line of warning") or diag explain $warnings;
  };
}
sub test_warnings_per_contrast {
  my ($design_tsv, $qc_issues_by_run, $expected_num_issues,$expected_warning_lines, $test_name) = @_;

  my $design = WbpsExpression::Model::Design::from_tsv(\$design_tsv);
  my $conditions_ordered = [$design->all_conditions];
  my @contrasts = map {
    my ($r, $t) = @{$_};
    $r && $t ? ([$r, $t, "$r vs $t"]) : () 
  } pairs @{$conditions_ordered};
  my ($amended_contrasts, $warnings) = 
    WbpsExpression::Analysis::QualityWarnings::amended_contrasts_and_warnings_for_per_contrast_analysis(
      $design, $qc_issues_by_run, \@contrasts 
  );
  subtest $test_name => sub {
    is_deeply([map {$_->[0]} @{$amended_contrasts}], [ map {$_->[0]} @contrasts], "references in contrasts not changed");
    is_deeply([map {$_->[1]} @{$amended_contrasts}], [ map {$_->[1]} @contrasts], "tests in contrasts not changed");
    my @z = pairs zip @contrasts, @{$amended_contrasts};
    my $contrasts_changed_names = grep {$_->[0][2] ne $_->[1][2]} @z;
    is($contrasts_changed_names, $expected_num_issues, "Each expected issue is a changed contrast name");
    is(scalar @{$warnings}, $expected_warning_lines, "Warning lines as expected") or diag explain $warnings;
  };
}
sub test_tsv {
  my $result = "Run\tCondition\n";
  my $c = 0;
  for my $x (@_){
    $c++;
    for my $i (1..$x){
      $result .= sprintf("run_%s%s\tcondition_%s\n", $c, $i, $c);
    }
  }
  return $result;
}
test_warnings_per_contrast(test_tsv(), {}, 0, 0, "Null case contrasts");
test_warnings_per_contrast(test_tsv(3,3), {}, 0, 0, "All is good - one contrast");
test_warnings_per_contrast(test_tsv(3,3,3,3), {}, 0, 0, "All is good - two contrasts");
test_warnings_per_contrast(test_tsv(2,3), {}, 1, 1, "One low rep - one contrast");
test_warnings_per_contrast(test_tsv(2,3,3,3), {}, 1, 1, "One low rep - two contrasts");
test_warnings_per_contrast(test_tsv(2,3,2,3), {}, 2, 1, "All low reps - two contrasts");
test_warnings_per_contrast(test_tsv(3,3), {"run_11" => ["Low qc run_11"]}, 1,1, "One low rep - one contrast");
test_warnings_per_contrast(test_tsv(3,3), {"run_11" => ["Low qc run_11"], "run_12" => ["Low qc run_12"]}, 1,1, "One low rep - one contrast");
test_warnings_per_contrast(test_tsv(3,3,3,3) , {"run_11" => ["Low qc run_11"]}, 1,1, "One low rep - two contrasts");

test_warnings_per_condition(test_tsv(), {}, 0, "Null case conditions");
test_warnings_per_condition(test_tsv(3), {}, 0, "All good - one condition");
test_warnings_per_condition(test_tsv(3, 3), {}, 0, "All good - two conditions");
test_warnings_per_condition(test_tsv(2), {}, 1, "One low rep - one condition");
test_warnings_per_condition(test_tsv(2, 3), {}, 1, "One low rep - two conditions");
test_warnings_per_condition(test_tsv(3), {"run_11" => ["Low qc run_11"]}, 1, "One low qc - one condition");
test_warnings_per_condition(test_tsv(3, 3), {"run_11" => ["Low qc run_11"]}, 1, "One low qc - two conditions");

done_testing;
