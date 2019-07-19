use strict;
use warnings;

package WbpsExpression::IncomingStudies;
use WbpsExpression::IncomingStudies::CurationDefaults;
use WbpsExpression::IncomingStudies::RnaseqerResults;
use WbpsExpression::IncomingStudies::RnaseqerFtp;
use WbpsExpression::IncomingStudies::StudyMetadata;
use WbpsExpression::IncomingStudies::SampleAttributes;
use WbpsExpression::Study;
use File::Basename qw/dirname/;
use File::Slurp qw/read_dir write_file/;
use File::Path qw/make_path/;
use List::Util qw/first pairmap/;
use List::MoreUtils qw/uniq singleton duplicates all/;
use Log::Any qw/$log/;

use Data::Compare;

# use Smart::Comments '###';

#Add to me for unit testing
our %exceptions = (
  "ERP006987" => "C. sinensis: four runs, but the tpms could be meaningful",
  "SRP013211" => "O. viverrini: two runs, juvenile vs adult, from 2012",
  "DRP003063" => "Emu study, four runs",
  "SRP131874" => "E. granulosus, four runs but a second study for the species",
  "SRP152065" => "M. incognita IARI study, two controls two replicates",
  "SRP179824" => "S. japonicum, three runs done for genome sequencing. Essentially a nice start",
  "SRP140458" => "T. pseudospiralis, three runs for genome sequencing, different life stages",
  "SRP048819" => "D. immitis, three runs, the community expressed interest in including this study",
  "SRP067884" => "T. circ, drug resistant vs suseptible. Only annotated at our competitors', nematode.net.",
);

sub same_runs {
  my ( $xs, $ys ) = @_;
  return Compare( uniq sort $xs, uniq sort $ys);
}

sub design_and_skipped_runs_from_data_by_run_and_previous_values {
  my ($characteristics_by_run, $replicates_by_run, $design_now, $skipped_runs_now) = @_;

  # Skipped runs stay skipped
  delete $characteristics_by_run->{$_} for @{$skipped_runs_now};
  
  my @skipped_runs;
  for my $run_id (keys %{$replicates_by_run}){
     unless ($characteristics_by_run->{$run_id}){
        delete $replicates_by_run->{$run_id};
        push @skipped_runs, $run_id;
     }
  }
  
  my %flattened = map {pairmap {"$a\t$b" => $a} %{$_}} values %{$characteristics_by_run};
  my @characteristic_types_varying_in_study = sort {$a cmp $b} duplicates values %flattened;
  my @characteristic_types_constant_across_study = sort {$a cmp $b} singleton values %flattened;
  my $characteristics_in_order = [@characteristic_types_varying_in_study, @characteristic_types_constant_across_study];
  for my $run_id (keys %{$replicates_by_run}){
    for my $ch (@{$characteristics_in_order}){
       $characteristics_by_run->{$run_id}{$ch} //= "";
    }
  }

  my $conditions_by_run_now = $design_now->conditions_by_run;
  my %conditions_by_run;
  for my $run_id (keys %{$replicates_by_run}){
    $conditions_by_run{$run_id} = 
      $conditions_by_run_now->{$run_id}
      || $characteristics_by_run->{$run_id}{sample_name}
      || join(", ", grep {$_} map {
            my $k = $_;
            my $v = $characteristics_by_run->{$run_id}{$k};
            $v && $v =~ /^[A-Z0-9]{1,6}$/ && length $k < 10 ? "$k $v" : $v
         } (@characteristic_types_varying_in_study 
            || qw/developmental_stage age sex host_infection organism_part strain isolate treatment rnai irradiation plane_of_amputation rnai_feedings/
         ))
      || $replicates_by_run->{$run_id};
  }

  my $design = keys %{$replicates_by_run} ? WbpsExpression::Study::Design::from_data_by_run($replicates_by_run, \%conditions_by_run, $characteristics_by_run, $characteristics_in_order): WbpsExpression::Study::Design::empty;
  return $design, \@skipped_runs;
}
sub update_study_with_results {
  my ( $path, $species, $study_id, $rnaseqer_last_update, $location_by_run, $quality_by_run, $replicates_by_run ) = @_;
  my @all_runs = sort keys %$location_by_run;
  my ( $design, $skipped_runs, $sources);
  my $study_now = WbpsExpression::Study->from_folder($path);
  if ( $study_now and same_runs(
      [ $study_now->all_runs],
      \@all_runs,
    )) {
    $design       = $study_now->{design};
    $skipped_runs = $study_now->{skipped_runs};
# Original curation used 30 as minimum quality. Take them.
  } elsif ($study_now and same_runs([$study_now->{design}->all_runs], grep {$quality_by_run->{$_} >= 30 } @all_runs)) {
    $design = $study_now->{design};
    $skipped_runs = [grep {my $run= $_; not (grep {$run eq $_} $study_now->{design}->all_runs)} @all_runs]; 
# Skipped studies stay skipped
  } elsif  ( $study_now and $study_now->{design}->is_empty){
    $design = WbpsExpression::Study::Design::empty;
    $skipped_runs = \@all_runs;
  } elsif  ( keys %{$location_by_run} < 6 and not $exceptions{$study_id} ){
    $design = WbpsExpression::Study::Design::empty;
    $skipped_runs = \@all_runs;
  } else {
    # New study, or more rarely there were runs added to the study
    my $design_now = $study_now ? $study_now->{design} : WbpsExpression::Study::Design::empty;
    my $skipped_runs_now = $study_now ? $study_now->{skipped_runs} : [];
    my $characteristics_by_run = WbpsExpression::IncomingStudies::SampleAttributes::characteristics_by_run($species, $study_id);

    ( $design, $skipped_runs ) = design_and_skipped_runs_from_data_by_run_and_previous_values(
      $characteristics_by_run,
      $replicates_by_run,
      $design_now, $skipped_runs_now
    );
  }

  my %sources = map {
    $_ => {
      quality => $quality_by_run->{$_},
      location => $location_by_run->{$_},
      end =>
        ( $study_now and $study_now->{config}{rnaseqer_last_update} // "" eq $rnaseqer_last_update and $study_now->{sources}{$_}{end} )
        || WbpsExpression::IncomingStudies::RnaseqerFtp::get_end_for_run($_, $location_by_run->{$_}, $rnaseqer_last_update),
    }
  } @all_runs;


  my $study_metadata =
    $study_now && ( all {$study_now->{config}{$_}} qw/title submitting_centre ena_first_public ena_last_update/ )
    ? $study_now->{config}
    : WbpsExpression::IncomingStudies::StudyMetadata::get( $species, $study_id );

  return unless $study_metadata;

  my $contrasts =
    WbpsExpression::IncomingStudies::CurationDefaults::contrasts($design);
  my $category =
    WbpsExpression::IncomingStudies::CurationDefaults::category( $design, $contrasts );

  my $config = {
    %{$study_metadata},
    contrasts => $contrasts,
    category => $category,
    rnaseqer_last_update => $rnaseqer_last_update,
  };

  my $study =
    WbpsExpression::Study->new($study_id, $design, $config, $skipped_runs, \%sources);
#### $study
  $study->to_folder($path);
  return $study;
}

sub update_studies {
  my ( $studies_dir, $species, $assembly ) = @_;
  my @studies;
  my @other_studies;
  my $rnaseqer_results_by_study_id = WbpsExpression::IncomingStudies::RnaseqerResults::get_results_by_study($species);
#### $rnaseqer_results_by_study_id
  for my $study_id (sort keys %{$rnaseqer_results_by_study_id}){
    next unless $rnaseqer_results_by_study_id->{$study_id}{assembly_used} eq $assembly; 
    next if $study_id eq "DRP002615";
    $log->info( __PACKAGE__ . " processing $study_id");
    my $study_path = join("/", $studies_dir, $species, $study_id);
    my $study = update_study_with_results($study_path, $species, $study_id,
        $rnaseqer_results_by_study_id->{$study_id}{rnaseqer_last_update},
       $rnaseqer_results_by_study_id->{$study_id}{location_by_run},
       $rnaseqer_results_by_study_id->{$study_id}{quality_by_run},
       $rnaseqer_results_by_study_id->{$study_id}{replicates_by_run},
    );
     $log->info( __PACKAGE__ . "Could not construct $study_id - skipping") unless $study;
    next unless $study;
    my $passes_checks = $study->passes_checks;
    $log->info( __PACKAGE__ . ": Study failing checks - see $study_path") unless $passes_checks;
    if ( not $study->{design}->is_empty and $passes_checks) {
       push @studies, $study;
    } else {
      push @other_studies, $study;
    }
  }
  return ( \@studies, \@other_studies );
}
1;
