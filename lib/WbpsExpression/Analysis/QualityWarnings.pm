use strict;
use warnings;

package WbpsExpression::Analysis::QualityWarnings;
use List::Util qw/pairmap pairs pairgrep uniq/;
use WbpsExpression::Model::Design;

sub conditions_amended_names_and_warnings_for_per_condition_analysis {
    my ( $design, $qc_issues_by_run, $conditions_ordered ) = @_;
    my ( $low_qc_by_condition, $low_replicate_by_condition ) =
      low_qc_and_replicate_by_condition( $design, $qc_issues_by_run,
        $conditions_ordered );
    my %conditions_amended_names;
    my @warnings;
    for my $condition ( @{$conditions_ordered} ) {
        my @low_qcs =
          sort map { ucfirst $_ } @{ $low_qc_by_condition->{$condition} // [] };
        my $low_reps = $low_replicate_by_condition->{$condition};
        my $warning =
          join( ". ", @low_qcs, $low_reps ? "Low replicates ($low_reps)" : () );
        push @warnings, "!$condition: $warning" if $warning;
        $conditions_amended_names{$condition} =
          ( @low_qcs || $low_reps ) ? "!$condition: " : $condition;
    }
    return \%conditions_amended_names, \@warnings;
}

sub amended_contrasts_and_warnings_for_per_contrast_analysis {
    my ( $design, $qc_issues_by_run, $contrasts ) = @_;
    my @conditions_ordered = uniq map { ( $_->[0], $_->[1] ) } @{$contrasts};
    my ( $low_qc_by_condition, $low_replicate_by_condition ) =
      low_qc_and_replicate_by_condition( $design, $qc_issues_by_run,
        \@conditions_ordered );

    my @amended_contrasts;
    my @warnings;
    my $contrasts_low_replicates;
    for my $contrast ( @{$contrasts} ) {
        my ( $reference, $test, $contrast_name ) = @{$contrast};
        my @low_qcs_reference =
          map { "$reference - $_" }
          @{ $low_qc_by_condition->{$reference} // [] };
        my @low_qcs_test =
          map { "$test - $_" } @{ $low_qc_by_condition->{$test} // [] };
        my $low_reps = $low_replicate_by_condition->{$reference}
          || $low_replicate_by_condition->{$test};
        my $contrast_amended_name =
          ( @low_qcs_reference || @low_qcs_test || $low_reps )
          ? "!$contrast_name"
          : $contrast_name;
        push @amended_contrasts, [ $reference, $test, $contrast_amended_name ];
        my $warning = join( ". ",
            sort map { ucfirst $_ } @low_qcs_reference,
            @low_qcs_test );
        push @warnings, $warning if $warning;
        $contrasts_low_replicates++ if $low_reps;
    }

    # We don't allow contrasts with 1 replicate. So if low, then 2.
    if ($contrasts_low_replicates) {
        push @warnings,
          $contrasts_low_replicates == @{$contrasts}
          ? "!Data quality warning: all contrasts based on conditions with low (2) replicates"
          : sprintf(
            "!Data quality warning: %s out of %s contrasts (marked with !) based on conditions with low (2) replicates",
            $contrasts_low_replicates, scalar @{$contrasts} );
    }
    return \@amended_contrasts, \@warnings;
}

sub low_qc_and_replicate_by_condition {
    my ( $design, $qc_issues_per_run, $conditions ) = @_;
    my @conditions = @{$conditions};

    my %runs_by_condition_all = %{ $design->runs_by_condition };
    my %runs_by_condition;
    @runs_by_condition{@conditions} = @runs_by_condition_all{@conditions};

    my %replicates_by_condition_all = %{ $design->replicates_by_condition };
    my %replicates_by_condition;
    @replicates_by_condition{@conditions} =
      @replicates_by_condition_all{@conditions};

    return low_qc_by_condition( \%runs_by_condition, $qc_issues_per_run ),
      low_replicate_by_condition( \%replicates_by_condition );
}

sub low_qc_by_condition {
    my ( $runs_by_condition, $qc_issues_per_run ) = @_;
    my %result = ();
    my %d      = %{$runs_by_condition};
    for my $c ( keys %d ) {
        my @runs = @{ $d{$c} };
        my %qcs;
        for my $run_id (@runs) {
            for my $qc ( @{ $qc_issues_per_run->{$run_id} } ) {
                push @{ $qcs{$qc} }, $run_id;
            }
        }
        for my $qc ( keys %qcs ) {
            push @{ $result{$c} },
              sprintf( "%s: %s %s",
                $qc,
                ( @{ $qcs{$qc} } > 1 ? "runs" : "run" ),
                join( ", ", @{ $qcs{$qc} } ) );
        }
    }
    return \%result;
}

sub low_replicate_by_condition {
    my ($replicates_by_condition) = @_;
    my %result                    = ();
    my %d                         = %{$replicates_by_condition};
    for my $c ( keys %d ) {
        my $number_replicates = @{ $d{$c} };
        $result{$c} = $number_replicates if $number_replicates < 3;
    }
    return \%result;
}
1;
