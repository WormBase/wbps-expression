use strict;
use warnings;
package Production::CurationDefaults;
use PublicResources::Rnaseq;
use List::MoreUtils qw/uniq/;
use Model::Study;
use Model::Design;

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

sub study {
  my (%args) = @_;
  my $design = design_from_runs(@{$args{runs}});
  return Model::Study->new($args{study_id}, $design , {
     condition_names => condition_names(@{$args{runs}}),
     slices => $design->slices_keys,
     title => $args{study_description_short},
     description => $args{study_description_full},
     pubmed => $args{pubmed},
     public => 0,
  });
}
1;
