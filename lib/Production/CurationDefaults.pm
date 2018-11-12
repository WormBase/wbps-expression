use strict;
use warnings;
package Production::CurationDefaults;
use PublicResources::Rnaseq;
use List::MoreUtils qw/uniq/;
use Model::Study;
use Model::Design;
#use Smart::Comments;
sub design_from_runs {
  my (@runs) = @_;
  my %conditions_per_run = map {($_->{run_id}, $_->{run_description_short})} @runs;
  my %characteristics_per_run = map {($_->{run_id}, $_->{characteristics})} @runs;
  my @characteristics_in_order = sort {$a cmp $b} uniq map {keys %{$_->{characteristics}}} @runs;
  return Model::Design::from_data_by_run(\%conditions_per_run, \%characteristics_per_run, \@characteristics_in_order);
}
sub condition_names {
  my (@runs) = @_;
  my %result;
  for my $run (@runs){
    my ($rs, $rf) = ($run->{run_description_short}, $run->{run_description_full});
    $result{$rs} = $rf unless $rs eq $rf;
  }
  return \%result;
}

# If two conditions vary by exactly one characteristic, they're probably a contrast
# This won't work if  multiple characteristics imply each other
sub auto_contrasts {
   my($design) = @_;
### require: ref $design eq 'Model::Design'
   my %vs = %{$design->{values}{by_condition} || {}};
   my @characteristics_varying_by_condition = uniq map {keys %{$_}} values %vs;
   my @contrasts;
   for my $c1 (keys %vs){
      for my $c2 (keys %vs){
         next if $c1 ge $c2;
         my ($ch, @others) = grep {$vs{$c1}{$_} ne $vs{$c2}{$_}} @characteristics_varying_by_condition;
         push @contrasts, {
            conditions => [$c1,$c2],
            characteristics => [$ch],
            name => sprintf("%s: %s vs %s", $ch, $vs{$c1}{$ch}, $vs{$c2}{$ch}),
         } if $ch and not @others;
      }
   }
   return \@contrasts;
}

sub study {
  my (%args) = @_;
  my $design =  design_from_runs(@{$args{runs}});
  return Model::Study->new($args{study_id}, $design , {
     condition_names => condition_names(@{$args{runs}}),
     contrasts => auto_contrasts($design),
     title => $args{study_description_short},
     description => $args{study_description_full},
     pubmed => $args{pubmed},
     public => 0,
  });
}
1;
