use strict;
use warnings;
package Model::Design;
use Text::CSV qw/csv/;
use List::Util qw/pairmap/;
use List::MoreUtils qw/uniq/;
#use Smart::Comments '###';

sub new {
  my ($class, %args) = @_;
  return bless \%args, $class; 
}

sub from_tsv {
  my ($path) = @_;
  my %characteristics_per_run;
  my %conditions_per_run;

  my ($header, @lines) = @{csv(in=>$path, sep_char => "\t")};
  my ($_1,$_2, @characteristics) = @{$header};
  for my $line (@lines){
### require: @$line == 2 + @characteristics 
     my ($run_id,$condition, @values) = @$line;
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
### require: %inverse <= %$original
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
  my %characteristics;
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
### require: $value
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
sub all_runs {
  return keys %{shift->{conditions_per_run}};
}
sub all_conditions {
  return keys %{reverse_hoa(shift->{conditions_per_run})};
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
  my %runs_by_condition = %{reverse_hoa($self->{conditions_per_run})};
  for my $condition (sort keys %runs_by_condition){
     for my $run_id (sort @{$runs_by_condition{$condition}}){
        print $fh join ("\t", $run_id, $condition, map {$self->value_in_run($run_id, $_)} @a)."\n";
     }
  }  
  close $fh; 
}

sub data_quality_checks {
  my ($self) = @_;
  my %runs_by_condition = %{reverse_hoa($self->{conditions_per_run})};
  my %characteristics_by_run = %{$self->{values}{by_run} || {}};
  return (
    "Some runs",
       => scalar %runs_by_condition,
    "Some characteristics" 
       => 0+@{$self->{characteristics_in_order}},
    "If there are multiple conditions, then some characteristics vary by condition"
       => (%runs_by_condition < 2 or %{$self->{values}{by_condition}} > 0 ),
     pairmap { 
       "Condition $a is well-defined: uniform characteristics in runs @{$b}" 
         => not scalar grep {$_} @characteristics_by_run{@{$b}}
     } %runs_by_condition
  ); 
}
1;
