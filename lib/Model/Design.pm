use strict;
use warnings;
package Model::Design;
use Text::CSV qw/csv/;
use List::Util qw/pairmap all/;
use List::MoreUtils qw/uniq /;
use open ':encoding(utf8)';
#use Smart::Comments '###';
sub new {
  my ($class, %args) = @_;
  return bless \%args, $class; 
}

sub from_tsv {
  my ($path) = @_;
  my %characteristics_per_run = ();
  my %conditions_per_run = ();

  my ($header, @lines) = @{csv(in=>$path, sep_char => "\t")};
  my ($_1,$_2, @characteristics) = @{$header};
  for my $line (@lines){
### require: @$line == 2 + @characteristics 
     my ($run_id,$condition, @values) = @$line;
     $_ =~ s/^\s+|\s+$//g for @values;
     $conditions_per_run{$run_id} = $condition;
     my %h;
     @h{@characteristics} = @values;
     $characteristics_per_run{$run_id} = \%h;
  }
  return from_data_by_run(\%conditions_per_run, \%characteristics_per_run, \@characteristics);
}

sub reverse_hoa {
  my $original = shift;
  my %inverse;
  push @{ $inverse{ $original->{$_} } }, $_ for keys %$original;
### require: keys %inverse <= keys %$original
  return \%inverse;
}

sub from_data_by_run {
  my ($conditions_per_run,$characteristics_per_run, $characteristics_in_order) = @_;  
  my %runs_by_condition = %{reverse_hoa($conditions_per_run)};

  my %characteristics_varying_across_conditions;
  my %characteristics_varying_within_any_condition;
  for my $characteristic (@{$characteristics_in_order}){
     $characteristics_varying_across_conditions{$characteristic} += -1 + uniq map {$characteristics_per_run->{$_}{$characteristic}} keys %{$conditions_per_run};
     for my $condition (keys %runs_by_condition){
        $characteristics_varying_within_any_condition{$characteristic} += -1 + uniq map {$characteristics_per_run->{$_}{$characteristic}} @{$runs_by_condition{$condition}};
     }
  }
  my %characteristics = (by_run => {}, by_condition => {}, common => {});
#### $characteristics_in_order
#### %characteristics_varying_within_any_condition
#### %characteristics_varying_across_conditions
#### $characteristics_per_run

  for my $characteristic (@{$characteristics_in_order}){
    if ($characteristics_varying_within_any_condition{$characteristic}){
      for my $run_id (keys %{$characteristics_per_run}){
         $characteristics{by_run}{$run_id}{$characteristic} = $characteristics_per_run->{$run_id}{$characteristic} // ""; 
      }
    } elsif($characteristics_varying_across_conditions{$characteristic}) {
      for my $condition (keys %runs_by_condition){
         my ($value, @others) = uniq map {$characteristics_per_run->{$_}{$characteristic}} @{$runs_by_condition{$condition}};
### require: not @others
         $characteristics{by_condition}{$condition}{$characteristic} = $value // ""; 
      }
    } else {
      my ($value, @others) = uniq map {$characteristics_per_run->{$_}{$characteristic}} keys %{$conditions_per_run};
### require: defined $value
### require: not @others
      $characteristics{common}{$characteristic} = $value;
    }
  }
  return new(__PACKAGE__, conditions_per_run => $conditions_per_run, values => \%characteristics, characteristics_in_order => $characteristics_in_order);
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
   my $value_in_condition = $self->common_value($characteristic) // $self->{values}{by_condition}{$condition}{$characteristic};
#### $value_in_condition
   return $value_in_condition;
}
sub value_in_run {
   my ($self, $run_id, $characteristic) = @_;
#### value in run: $run_id, $characteristic
   my $value_in_run = $self->value_in_condition($self->{conditions_per_run}{$run_id}, $characteristic) // $self->{values}{by_run}{$run_id}{$characteristic};
#### $value_in_run
   return $value_in_run;
}
sub runs_by_condition {
  return reverse_hoa(shift->{conditions_per_run});
}
sub condition_run_ordered_pairs {
  my %d = %{shift->runs_by_condition};
  my @result;
  for my $condition (sort keys %d){
     for my $run_id (sort @{$d{$condition}}){
        push @result, [$condition, $run_id];
     }
  }
  return @result;
}
sub all_runs {
  my @a = map {$_->[1]} shift->condition_run_ordered_pairs;
  return @a;
}
sub all_conditions {
  my @a = uniq map {$_->[0]} shift->condition_run_ordered_pairs;
  return @a;
}
sub characteristics_per_run {
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
  print $fh join("\t", "Run", "Condition", @a)."\n";
#### $self
  for my $p ($self->condition_run_ordered_pairs){
     print $fh join ("\t", $p->[1], $p->[0], map {$self->value_in_run($p->[1], $_)} @a)."\n";
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
  my %runs_by_condition = %{$self->runs_by_condition};
  my %characteristics_by_run = %{$self->{values}{by_run}};
  my @conditions_well_defined = pairmap {
    "Condition $a should have uniform characteristics in runs @{$b}"
      => 2 > uniq map {$_ ? join("", sort %$_): ()} @characteristics_by_run{@{$b}}
  } %runs_by_condition;
  my @characteristics_varying_by_condition = $self->characteristics_varying_by_condition;
  my @conditions_unique = pairmap {
    "Characteristics $a should define precisely one condition"
     => (@{$b} == 1 )
  } %{reverse_hoa({ map {
      my $c = $_; 
      my $s = join("\t", map {$self->value_in_condition($c, $_)} @characteristics_varying_by_condition );
      $c => $s
   } $self->all_conditions})};
  return (
    "Study should have some runs",
       => scalar %runs_by_condition,
    "Study should have some characteristics" 
       => 0+@{$self->{characteristics_in_order}},
    "Characteristics should have non-blank names"
       => ( 0 ==  grep {not $_} @{$self->{characteristics_in_order}}),
    "Conditions should have non-blank names",
       => ( 0 ==  grep {not $_} $self->all_conditions),
    "Conditions should have reasonably short names - below 60 chars",
       => ( 0 ==  grep {length $_ > 60 } $self->all_conditions),
    "Some characteristics should vary by condition"
       => (2 > $self->all_runs or 0 < $self->characteristics_varying_by_condition ),
     @conditions_well_defined,
     @conditions_unique,
  ); 
}
sub passes_checks {
  my %checks = shift->data_quality_checks;
  return all {$_} values %checks;
}
1;
