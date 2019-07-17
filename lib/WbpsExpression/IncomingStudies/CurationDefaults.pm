use strict;
use warnings;

package WbpsExpression::IncomingStudies::CurationDefaults;
use List::Util qw/pairmap/;
use List::MoreUtils qw/uniq all any/;
use WbpsExpression::Study;
use WbpsExpression::Study::Design;
use Data::Compare;
# use Smart::Comments '###';

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
  my ( $design, $subset_chs, $reference, $test ) = @_;
  my @factors = grep {$design->value_in_condition($reference, $_) ne $design->value_in_condition($test, $_)} @{$subset_chs};
  my @reference_values = map {$design->value_in_condition($reference, $_) } @factors;
  my $reference_short = trim_values($reference, @reference_values);
  my @test_values = map {$design->value_in_condition($test, $_) } @factors; 
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
  return $x && $y && ($x eq "concentration" || $x eq "timepoint")  && $y eq "treatment" && not @others; 
}
sub try_make_time_series {
   my ($subset_chs, @conditions) = @_;
   return () unless @conditions >=5;
   my ($ch, @other_chs) = @{$subset_chs};
   return () unless $ch eq 'timepoint' and not @other_chs;
   return () unless all {$_ =~ /.*, ?\d+ \w+$/} @conditions;
   return () unless 1 == scalar uniq map {$_ =~ s/^(.*, )?\d+ (\w+)s?$/$1\t$2/ and $_} map {$_} @conditions;
   my ($first_c, @other_cs) = sort {
     (my $aa = $a) =~ s/^.*, ?(\d+) \w+$/$1/;
     (my $bb = $b) =~ s/^.*, ?(\d+) \w+$/$1/;
     $aa <=> $bb
   } @conditions;
   my @result = map {[$first_c, $_] } @other_cs;
#### try_make_time_series in: $subset_chs, @conditions
#### try_make_time_series result: @result
   return @result;
}
sub clear_control_characteristic_in_drug_assay {
  my ($design, $c1, $c2) = @_;
  my ($control_value, @other_control_values) = uniq grep {not $_ or $_ =~ /control/i } map {$design->value_in_condition($_, 'treatment') } $design->all_conditions;
  return @other_control_values ? "" : $control_value;
}

sub if_clear_control_value_then_use_as_reference {
  my ($subset_chs, $design, $c1, $c2) = @_;
  if (is_drug_assay(@{$subset_chs})){
    my $control_value = clear_control_characteristic_in_drug_assay($design, $c1, $c2);
    if ($control_value){
      if ($design->value_in_condition($c1, 'treatment') eq $control_value) {
        return [$c1, $c2];
      } elsif ($design->value_in_condition($c2, 'treatment') eq $control_value) {
        return [$c2, $c1];
      } else {
        return ();
      } 
    }
  }
  if ($c2 =~ /control|no treatment/i and $c1 !~ /control|no treatment/i){
     return [$c2, $c1];
  }
  return [$c1, $c2];
}
sub subset_is_rnai_in_study_with_varying_timepoints {
  my ($subset_chs, $chs) = @_;
  my ($x,@xs) = @$subset_chs;
  my $result = (grep {$_ eq "timepoint"} @$chs) && ($x eq "rnai") && not @xs;
  return $result;
}
sub contrasts {
  my ($design)   = @_;
  return [] unless $design->passes_checks;
  my %replicates = pairmap {$a => scalar @$b } %{$design->replicates_by_condition};
  # Soft minimum of replicates is 3, but we also allow 2 vs 3
  my @conditions = grep {$replicates{$_} >= 2 } $design->all_conditions;
  return [] unless grep { $replicates{$_} >= 3 } @conditions;
  my @chs        = $design->characteristics_not_common;
#### All characteristics varying by conditions somewhere in the design: @chs
  my @result;
  die "Too many characteristics to iterate through subsets: @chs" if @chs > 10;
  for my $s ( grep { @{$_} < 4 and @{$_} > 0 } subsets(@chs) ) {
    next if subset_is_rnai_in_study_with_varying_timepoints($s, \@chs);
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
    my @pairs =
       map { 
        my @as = try_make_time_series(\@subset_chs, @{$_});
        my @bs = n_choose_two( @{$_} );
        @as ? @as : @bs
      }
      sort { join( "", @{$a} ) cmp join( "", @{$b} ) }
      map { [ sort @{$_} ] } @partition;
    my @contrasts =
      map {
      [ $_->[0], $_->[1], contrast_name( $design, \@subset_chs, @{$_} ) ]
      }
      grep {
        ! is_life_cycle(@subset_chs) || if_value_defined_then_empty("treatment", $design, @{$_})
      }
      map {
        if_clear_control_value_then_use_as_reference(\@subset_chs, $design, @{$_})
      }
      grep {
        ! is_drug_assay(@subset_chs) || is_contrast_across_characteristic("treatment", $design, @{$_}) 
      }
      grep {
         not ($replicates{$_->[0]} < 3 and $replicates{$_->[1]} < 3)
      } @pairs;
    
    push @result,
      {
      name   => join( "+", @subset_chs ),
      values => \@contrasts,
      } if @contrasts;
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
my %life_stage_categories = (
  developmental_stage => 1,
  age => 1,
  sex => 1,
  host_infection => 1,
);
my %treatment_categories = (
  treatment => 1,
  rnai => 1,
  irradiation => 1,
  plane_of_amputation => 1,
  rnai_feedings => 1,
  culture_conditions => 1,
);
sub category {
  my ($design, $contrasts) = @_;
  my @chs        = $design->characteristics_not_common;
  return "Other" if (
    any {@{$_} < 2 } values %{ $design->replicates_by_condition}
  );
  return "Life cycle" if (
     all { $life_stage_categories{$_}} @chs and scalar $design->all_conditions > 2
  );
  return "Organism parts" if (
     all { $_ eq "organism_part" || $life_stage_categories{$_}} @chs and any { $_ eq "organism_part"} @chs
  );
  return "Variation within species" if (
     all { $_ eq "isolate" || $_ eq "strain" || $life_stage_categories{$_}} @chs and any { $_ eq "isolate" || $_ eq "strain" } @chs
  );
  return "Other" if not @{$contrasts};
  my $mentions_treatment = grep {$treatment_categories{$_}} @chs;
  my $mentions_cell_type = grep { $_ eq "cell_type"} @chs;
  return "Response to treatment" if (
    $mentions_treatment && ! $mentions_cell_type
  );
  return "Cell types" if (
   $mentions_cell_type && ! $mentions_treatment
  );
  return "Other";
}
1;
