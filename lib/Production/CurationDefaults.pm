use strict;
use warnings;

package Production::CurationDefaults;
use PublicResources::Rnaseq;
use List::MoreUtils qw/uniq/;
use Model::Study;
use Model::Design;
use Data::Compare;

#use Smart::Comments '###';
sub design_from_runs {
  my (@runs) = @_;
  my %conditions_per_run =
    map { ( $_->{run_id}, $_->{run_description_short} ) } @runs;
  my %replicates_per_run = 
    map { ( $_->{run_id}, $_->{sample_id})} @runs;
  # If no extra information in the samples, skip them
  if (Compare([keys %replicates_per_run] , [ values %replicates_per_run ]){
    %replicates_per_run = map { ( $_->{run_id}, $_->{run_id}) } @runs;
  }
  my %characteristics_per_run =
    map { ( $_->{run_id}, $_->{characteristics} ) } @runs;
  my @characteristics_in_order =
    uniq sort { $a cmp $b } map { keys %{ $_->{characteristics} } } @runs;
  return Model::Design::from_data_by_run( \%replicates_per_run, \%conditions_per_run,
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

sub all_characteristics_vary {
  my ( $design, $chs, $conditions ) = @_;
  for my $ch ( @{$chs} ) {
    return 0
      if 1 == uniq map { $design->value_in_condition( $_, $ch ) }
      @{$conditions};
  }
  return 1;
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
sub trim_values {
  my ($result, @values) = @_;
  $result =~ s/([,\s]*)$_(s|es)?([,\s]*)/$1/ for map {quotemeta $_ } @values;
  $result =~ s/^\W+|[,\s]+$//;
  return $result;
}
sub venn {
  my ($xs, $ys) = @_;
  my @xxs;
  my @is;
  my @yys;
  for my $x (@{$xs}){
    if( grep {$_ eq $x} @{$ys} ){
      push @is, $x;
    } else {
      push @xxs, $x;
    }
  }
  for my $y (@{$ys}){
    if (not grep {$_ eq $y} @$xs){
       push @yys, $y;
    }
  }
  return \@xxs, \@is, \@yys;
}

sub contrast_name {
  my ( $design, $factors, $reference, $test ) = @_;
  my @reference_values = map {$design->value_in_condition($reference, $_) } @{$factors};
  my $reference_short = trim_values($reference, @reference_values);
  my @test_values = map {$design->value_in_condition($test, $_) } @{$factors}; 
  my $test_short = trim_values($test,@test_values); 
  if($reference_short eq $test_short){
    my ($ref_only, $intersection, $test_only) = venn(\@reference_values, \@test_values);
    my @ref_only = grep {$_} @{$ref_only};
    my @test_only = grep {$_} @{$test_only};
    my $common = $reference_short || join(", ", grep {$_} @{$intersection});
    $common = "$common: " if $common;
    return $common
      . (@ref_only ? join (", ", @ref_only ) : "''")
      . " vs "
      . (@test_only ? join(", ", @test_only) : "''");
  } else {
    return "$reference vs $test";
  }
}
sub is_contrast_across_characteristic {
  my ( $characteristic, $design, $reference, $test ) = @_;
  my $r = $design->value_in_condition($reference, $characteristic) // "";
  my $t = $design->value_in_condition($test, $characteristic) // "";
#### is_contrast_across_characteristic : "$characteristic $reference - $r vs $test - $t"
  return $r ne $t;
}
sub if_value_defined_then_empty {
  my ( $characteristic, $design, $reference, $test ) = @_;
  my $r = $design->value_in_condition($reference, $characteristic);
  my $t = $design->value_in_condition($test, $characteristic);
  return (! defined $r || not $r) && (! defined $t || not $t);
}
sub is_life_cycle {
  my ($x, $y, @others) = sort @_;
  return $x && $y && $x eq "developmental_stage" && $y eq "sex" && not @others; 
}
sub is_drug_assay {
  my ($x, $y, @others) = sort @_;
  return $x && $y && $x eq "timepoint" && $y eq "treatment" && not @others; 
}
sub contrasts {
  my ($design)   = @_;
  my %replicates = %{$design->runs_by_condition};
  my @conditions = grep {@{$replicates{$_}} >= 2 } $design->all_conditions;
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
    next unless all_characteristics_vary($design, \@subset_chs, [ map {@$_} @partition]); 
    @partition = grep {
      is_different_values_on_all_characteristics( $design, \@subset_chs, $_ )
    } @partition unless is_life_cycle(@subset_chs) || is_drug_assay(@subset_chs);
#### Filter further to those that differ in all values chosen characteristics: @partition
    next unless @partition;
    my @contrasts =
      map {
      [ $_->[0], $_->[1], contrast_name( $design, \@subset_chs, @{$_} ) ]
      }
      grep {
        ! is_life_cycle(@subset_chs) || if_value_defined_then_empty("treatment", $design, @{$_})
      }
      grep {
        ! is_drug_assay(@subset_chs) || is_contrast_across_characteristic("treatment", $design, @{$_}) 
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
  if(grep {$_->{name} eq "developmental_stage+sex"} @result){
     @result = grep {$_->{name} ne "developmental_stage" && $_->{name} ne "sex"} @result;
  }
  if(grep {$_->{name} eq "timepoint+treatment"} @result){
     @result = grep {$_->{name} ne "timepoint" && $_->{name} ne "treatment"} @result;
  }
#### @result
  return \@result;
}
sub config_base {
  my (%args) = @_;
  return {
      title             => $args{study_description_short},
      description       => $args{study_description_full},
      pubmed            => $args{pubmed},
      submitting_centre => $args{attributes}{submitting_centre},
  };
}
sub study {
  my (%args) = @_;
  my $design = design_from_runs( @{ $args{runs} } );
  return Model::Study->new(
    $args{study_id},
    $design,
    {
      %{ config_base(%args) },
      condition_names => condition_names( @{ $args{runs} } ),
      contrasts       => contrasts($design),
    }
  );
}
1;
