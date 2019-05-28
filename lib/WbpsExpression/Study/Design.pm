use strict;
use warnings;
package WbpsExpression::Study::Design;
use Text::CSV qw/csv/;
use List::Util qw/pairmap all/;
use List::MoreUtils qw/uniq /;
use open ':encoding(utf8)';
# use Smart::Comments '###';
sub new {
  my ($class, %args) = @_;
  return bless \%args, $class; 
}

sub from_tsv {
  my ($path) = @_;
  my %characteristics_by_run = ();
  my %replicates_by_run = ();
  my %conditions_by_run = ();

  my ($header, @lines) = @{csv(in=>$path, sep_char => "\t")};
  my $runs_are_by_sample; 
  my ($_1,$_2,$_3, @characteristics) = @{$header};
  if ($_1 eq 'Run' and $_2 eq 'Sample' and $_3 eq 'Condition'){
    $runs_are_by_sample = 1;
  } elsif ($_1 eq 'Run' and $_2 eq 'Condition') {
    $runs_are_by_sample = '';
    @characteristics = $_3 ? ($_3, @characteristics) : ();
  } else {
    die "Bad header at $path: @$header";
  }
  for my $line (@lines){
     my @line = @{$line};
     my ($run_id, $replicate_id, $condition, @values) = $runs_are_by_sample ? @line : ($line[0], @line); 
### require: scalar @values == scalar @characteristics
     $_ =~ s/^\s+|\s+$//g for @values;
     $replicates_by_run{$run_id} = $replicate_id;
     $conditions_by_run{$run_id} = $condition;
     my %h;
     @h{@characteristics} = @values;
     $characteristics_by_run{$run_id} = \%h;
  }
  return from_data_by_run(\%replicates_by_run, \%conditions_by_run, \%characteristics_by_run, \@characteristics);
}

sub reverse_hoa {
  my $original = shift;
  my %inverse;
  push @{ $inverse{ $original->{$_} } }, $_ for keys %$original;
### require: keys %inverse <= keys %$original
  return \%inverse;
}

sub empty {
  return from_data_by_run({},{},{},[]);
}

sub is_empty {
  my ($self) = @_;
  return not $self->all_runs;
}

sub from_data_by_run {
  my ($replicates_by_run, $conditions_by_run, $characteristics_by_run, $characteristics_in_order) = @_;

  my %runs_by_replicate = %{reverse_hoa($replicates_by_run)};
  my $conditions_by_replicate = {};
  for my $replicate (keys %runs_by_replicate){
     my ($condition, @others) = uniq map {$conditions_by_run->{$_}} @{$runs_by_replicate{$replicate}};
     die "$replicate ambiguous condition: $condition @others" if @others; 
     $conditions_by_replicate->{$replicate} = $condition;
  }
  my %replicates_by_condition =  %{reverse_hoa($conditions_by_replicate)};

  my %characteristics_varying_across_conditions;
  my %characteristics_varying_within_any_condition;
  my %characteristics_varying_within_any_replicate;
  for my $characteristic (@{$characteristics_in_order}){
     $characteristics_varying_across_conditions{$characteristic} 
          += -1 + uniq map {$characteristics_by_run->{$_}{$characteristic}} keys %{$replicates_by_run};
     for my $replicate (keys %runs_by_replicate){
        $characteristics_varying_within_any_replicate{$characteristic} 
          += -1 + uniq map {$characteristics_by_run->{$_}{$characteristic}} @{$runs_by_replicate{$replicate}};
     }
     for my $condition (keys %replicates_by_condition){
        $characteristics_varying_within_any_condition{$characteristic} 
          += -1 + uniq map {$characteristics_by_run->{$_}{$characteristic}} map {@{$runs_by_replicate{$_}}} @{$replicates_by_condition{$condition}};
     }
  }
  my %characteristics = (by_run => {}, by_replicate => {}, by_condition => {}, common => {});
#### $characteristics_in_order
#### %characteristics_varying_within_any_condition
#### %characteristics_varying_within_any_replicate
#### %characteristics_varying_across_conditions
#### $characteristics_by_run

  for my $characteristic (@{$characteristics_in_order}){
    if ($characteristics_varying_within_any_replicate{$characteristic}){
      for my $run_id (keys %{$characteristics_by_run}){
         $characteristics{by_run}{$run_id}{$characteristic} = $characteristics_by_run->{$run_id}{$characteristic} // ""; 
      }
    } elsif ($characteristics_varying_within_any_condition{$characteristic}){
      for my $replicate (keys %runs_by_replicate){
         my ($value, @others) = uniq map {$characteristics_by_run->{$_}{$characteristic}} @{$runs_by_replicate{$replicate}};
### require: not @others
         $characteristics{by_replicate}{$replicate}{$characteristic} = $value // "";
      }
    } elsif($characteristics_varying_across_conditions{$characteristic}) {
      for my $condition (keys %replicates_by_condition){
         my ($value, @others) = uniq map {$characteristics_by_run->{$_}{$characteristic}} map {@{$runs_by_replicate{$_}}} @{$replicates_by_condition{$condition}};
### require: not @others
         $characteristics{by_condition}{$condition}{$characteristic} = $value // ""; 
      }
    } else {
      my ($value, @others) = uniq map {$characteristics_by_run->{$_}{$characteristic}} keys %{$conditions_by_run};
### require: defined $value
### require: not @others
      $characteristics{common}{$characteristic} = $value;
    }
  }
  return new(__PACKAGE__, conditions_by_replicate => $conditions_by_replicate, replicates_by_run => $replicates_by_run, values => \%characteristics, characteristics_in_order => $characteristics_in_order);
}
sub common_value {
   my ($self, $characteristic) = @_;
#### common value: $characteristic
   my $common_value = $self->{values}{common}{$characteristic};
#### $common_value
   return $common_value;
}
sub value_in_condition {
   my ($self, $condition, $characteristic) = @_;
#### value in condition: $condition, $characteristic
   my $value_in_condition = $self->common_value($characteristic);
   $value_in_condition //= $self->{values}{by_condition}{$condition}{$characteristic} if $self->{values}{by_condition}{$condition};
#### $value_in_condition
   return $value_in_condition;
}
sub value_in_replicate {
   my ($self, $replicate, $characteristic) = @_;
#### value in replicate: $replicate, $characteristic
   my $value_in_replicate = $self->value_in_condition($self->{conditions_by_replicate}{$replicate}, $characteristic);
   $value_in_replicate //= $self->{values}{by_replicate}{$replicate}{$characteristic} if $self->{values}{by_replicate}{$replicate};
#### $value_in_replicate
   return $value_in_replicate;
}
sub value_in_run {
   my ($self, $run_id, $characteristic) = @_;
#### value in run: $run_id, $characteristic
   my $value_in_run = $self->value_in_replicate($self->{replicates_by_run}{$run_id}, $characteristic);
   $value_in_run //= $self->{values}{by_run}{$run_id}{$characteristic} if $self->{values}{by_run}{$run_id};
#### $value_in_run
   return $value_in_run;
}
sub runs_by_condition_then_replicate {
  my ($self) = @_;
  my $h = reverse_hoa($self->{replicates_by_run});
  my $H = reverse_hoa($self->{conditions_by_replicate});
  my %result;
  for my $condition (keys %{$H}){
     for my $replicate (@{$H->{$condition}}){
        $result{$condition}{$replicate} = $h->{$replicate}; 
     }
  }
  return \%result;
}
sub replicates_by_condition {
  my %result = pairmap { $a => [ keys %$b ] } %{shift->runs_by_condition_then_replicate};
  return \%result;
}
sub runs_by_replicate {
  my %result = map {%$_} values %{shift->runs_by_condition_then_replicate};
  return \%result;
}
sub runs_by_condition {
  my %result = pairmap {$a => [ map { @$_ } values %{$b} ]} %{shift->runs_by_condition_then_replicate};
  return \%result;
}
sub conditions_by_run {
  my ($self) = @_;
  my %result;
  while (my ($run, $replicate) = each %{$self->{replicates_by_run}}){
    $result{$run} = $self->{conditions_by_replicate}{$replicate};
  } 
  return \%result;
}
sub condition_replicate_run_ordered_triples {
  my %d = %{shift->runs_by_condition_then_replicate};
  my @result;
  for my $condition (sort keys %d){
     for my $replicate (sort keys %{$d{$condition}}){
	   for my $run_id (sort @{$d{$condition}{$replicate}}){
		  push @result, [$condition, $replicate, $run_id];
	   }
     }
  }
  return @result;
}
sub all_runs {
  my @a = map {$_->[2]} shift->condition_replicate_run_ordered_triples;
  return @a;
}
sub all_replicates {
  my @a = uniq map {$_->[1]} shift->condition_replicate_run_ordered_triples;
  return @a;
}
sub all_conditions {
  my @a = uniq map {$_->[0]} shift->condition_replicate_run_ordered_triples;
  return @a;
}
sub characteristics_by_run {
  my ($self) = @_;
  my @all_runs = $self->all_runs;
  my %result;
  for my $run_id (@all_runs){
    for my $ch (@{$self->{characteristics_in_order}}){
      $result{$run_id}{$ch} = $self->value_in_run($run_id, $ch);
    }
  }
  return \%result;
}

sub to_tsv {
  my($self, $path) = @_;
  open(my $fh, ">", $path) or die "$path: $!";
#### $path
  my @a = @{$self->{characteristics_in_order}};
#### characteristics in order: @a
  my @ts = $self->condition_replicate_run_ordered_triples;
  my $needs_three_columns = @ts && not all {$_->[1] eq $_->[2]} @ts;
  print $fh join("\t", ($needs_three_columns ? qw/Run Sample Condition/ : qw/Run Condition/), @a)."\n";
#### $self
  for my $t (@ts){
     my ($condition, $replicate, $run) = @{$t};
     print $fh join ("\t", $run, ($needs_three_columns ? $replicate : ()), $condition, map {$self->value_in_run($run, $_)} @a)."\n";
  }
  close $fh; 
}
sub characteristics_varying_by_condition {
  my ($self) = @_;
  my %h = map {$_ => 1 } map {keys %{$_} } values %{$self->{values}{by_condition}};
  return grep {$h{$_}} @{$self->{characteristics_in_order}};
}
sub data_quality_checks {
  my ($self) = @_;
  my %characteristics_by_replicate = %{$self->{values}{by_replicate}};
  my %characteristics_by_run = %{$self->{values}{by_run}};
  my @replicates_well_defined = pairmap {
    "Replicate $a should have uniform characteristics in runs @{$b}"
      => 2 > uniq map {$_ ? join("", sort %$_): ()} @characteristics_by_run{@{$b}} 
  } %{$self->runs_by_replicate};

  my @conditions_well_defined = pairmap {
    "Condition $a should have uniform characteristics in replicates @{$b}"
      => 2 > uniq map {$_ ? join("", sort %$_): ()} @characteristics_by_replicate{@{$b}}
  } %{$self->replicates_by_condition};
#### replicates by condition: $self->replicates_by_condition
#### %characteristics_by_run
#### @conditions_well_defined
  my @characteristics_varying_by_condition = $self->characteristics_varying_by_condition;
  my @conditions_unique = pairmap {
    "Characteristics $a should define precisely one condition"
     => (@{$b} == 1 )
  } %{reverse_hoa({ map {
      my $c = $_; 
      my $s = join("\t", map {$self->value_in_condition($c, $_)} @characteristics_varying_by_condition );
      $c => $s
   } $self->all_conditions})};
  my $num_runs = scalar $self->all_runs;
  my $num_replicates = scalar $self->all_replicates;
  my $num_conditions = scalar $self->all_conditions;
  return (
    "Study should have some runs"
       => $num_runs,
    "Study should have some characteristics" 
       => 0+@{$self->{characteristics_in_order}},
    "Characteristics should have non-blank names"
       => ( 0 ==  grep {not $_} @{$self->{characteristics_in_order}}),
    "Conditions should have non-blank names",
       => ( 0 ==  grep {not $_} $self->all_conditions),
    "Conditions should have reasonably short names - below 60 chars"
       => ( 0 ==  grep {length $_ > 60 } $self->all_conditions),
    "Some characteristics should vary by condition"
       => (2 > $self->all_runs or 0 < $self->characteristics_varying_by_condition ),
     "If the study has fewer samples than replicates, it should have fewer conditions than samples"
       => ($num_runs < 2*$num_replicates || $num_conditions < $num_replicates ),
     @replicates_well_defined,
     @conditions_well_defined,
     @conditions_unique,
  ); 
}
sub passes_checks {
  my %checks = shift->data_quality_checks;
  return all {$_} values %checks;
}
1;
