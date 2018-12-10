use strict;
use warnings;

package Production::CurationDefaults;
use PublicResources::Rnaseq;
use List::MoreUtils qw/uniq/;
use Model::Study;
use Model::Design;

#use Smart::Comments '###';
sub design_from_runs {
  my (@runs) = @_;
  my %conditions_per_run =
    map { ( $_->{run_id}, $_->{run_description_short} ) } @runs;
  my %characteristics_per_run =
    map { ( $_->{run_id}, $_->{characteristics} ) } @runs;
  my @characteristics_in_order =
    uniq sort { $a cmp $b } map { keys %{ $_->{characteristics} } } @runs;
  return Model::Design::from_data_by_run( \%conditions_per_run,
    \%characteristics_per_run, \@characteristics_in_order );
}

sub condition_names {
  my (@runs) = @_;
  my %result;
  for my $run (@runs) {
    my ( $rs, $rf ) =
      ( $run->{run_description_short}, $run->{run_description_full} );
    $result{$rs} = $rf unless $rs eq $rf;
  }
  return \%result;
}

# adapted from Math::Subsets::List 1.008
sub subsets {
  my @result;
  my $n = scalar(@_);    # Size of list to be subsetted
  my $l = 0;             # Current item
  my @p = ();            # Current subset
  my @P = @_;            # List to be subsetted

  my $p;
  $p = sub {
    if ( $l < $n ) {
      ++$l;
      &$p();
      push @p, $P[ $l - 1 ];
      &$p();
      pop @p;
      --$l;
    }
    else {
      push @result, [@p];
    }
  };
  &$p;
  $p = undef;
  return @result;
}

sub partition_by_characteristics_values {
  my ( $design, $chs, $conditions ) = @_;
  my %result;
  for my $condition ( @{$conditions} ) {
    my $k = join "\t",
      map { $design->value_in_condition( $condition, $_ ) } @{$chs};
    push @{ $result{$k} }, $condition;
  }
  return values %result;
}

sub is_different_values_on_all_characteristics {
  my ( $design, $chs, $conditions ) = @_;
  for my $ch ( @{$chs} ) {
    return 0
      if @{$conditions} > uniq map { $design->value_in_condition( $_, $ch ) }
      @{$conditions};
  }
  return 1;
}

sub n_choose_two {
  my @result;
  for my $i ( 0 .. $#_ ) {
    for my $j ( $i + 1 .. $#_ ) {
      push @result, [ $_[$i], $_[$j] ];
    }
  }
  return @result;
}

sub contrast_name {
  my ( $design, $factors, $reference, $test ) = @_;
### I could be more complicated
  return "$reference vs $test";
}

sub contrasts {
  my ($design)   = @_;
  my %replicates = %{$design->runs_by_condition};
  my @conditions = grep {@{$replicates{$_}} >= 3 } $design->all_conditions;
  return [] unless @conditions;
  my @chs        = $design->characteristics_varying_by_condition;
#### All characteristics varying by conditions somewhere in the design: @chs
  my @result;
  die "Too many characteristics to iterate through subsets: @chs" if @chs > 10;
  for my $s ( grep { @{$_} < 4 } subsets(@chs) ) {
    my @subset_chs = @{$s};
#### Characteristics - chosen: @subset_chs
    my @other_chs = grep {
      my $ch = $_;
      not( grep { $_ eq $ch } @subset_chs )
    } @chs;
    my @partition = grep { @{$_} > 1 }
      partition_by_characteristics_values( $design, \@other_chs, \@conditions );
#### Conditions that share values on non-chosen characteristics: @partition
    @partition = grep { @{$_} > 1 } map {
      my $p = $_;
      my @partitioned_p =
        partition_by_characteristics_values( $design, \@subset_chs, $p );
      [ map { @{$_} == 1 ? @{$_} : () } @partitioned_p ]
    } @partition if @subset_chs;
#### Filter to those that also differ in values on chosen characteristics: @partition
    @partition = grep {
      is_different_values_on_all_characteristics( $design, \@subset_chs, $_ )
    } @partition;
#### Filter further to those that differ in all values chosen characteristics: @partition
    next unless @partition;
    my @contrasts =
      map {
      [ $_->[0], $_->[1], contrast_name( $design, \@subset_chs, @{$_} ) ]
      }
      map { n_choose_two( @{$_} ) }
      sort { join( "", @{$a} ) cmp join( "", @{$b} ) }
      map { [ sort @{$_} ] } @partition;
    push @result,
      {
      name   => join( "+", @subset_chs ),
      values => \@contrasts,
      };
  }
#### @result
  return \@result;
}

sub study {
  my (%args) = @_;
  my $design = design_from_runs( @{ $args{runs} } );
  return Model::Study->new(
    $args{study_id},
    $design,
    {
      condition_names => condition_names( @{ $args{runs} } ),
      title           => $args{study_description_short},
      description     => $args{study_description_full},
      pubmed          => $args{pubmed},
      contrasts       => contrasts($design),
    }
  );
}
1;
