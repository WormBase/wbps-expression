use strict;
use warnings;

package WbpsExpression::IncomingStudies::RnaseqerResults;
use Log::Any '$log';
use LWP;
use JSON;
use File::Basename;
use DateTime;
use DateTime::Format::Strptime;
use DateTime::Format::ISO8601::Format;
use List::Util qw/pairs unpairs/;
use File::Slurp qw(read_file read_dir);
# use Smart::Comments '###';

# 
# { study_id -> {
#     assembly_used,
#     rnaseqer_last_update,
#     location_by_run :: {Run -> URL},
#     quality_by_run :: {Run -> Int},
#     replicates_by_run :: {Run -> Replicate}
#   }
# }

sub get_json {
  my ($url) = @_;
  $log->info("RnqseqerResults::get_json LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;
  return $response->decoded_content && decode_json($response->decoded_content);
}

sub get_results_by_study {
  my ($species) = @_;
  my %runs_per_study;
  for my $o (@{get_json ("https://www.ebi.ac.uk/fg/rnaseq/api/json/20/getRunsByOrganism/$species") || []}){
    push @{$runs_per_study{$o->{STUDY_ID}}}, $o if $o->{STATUS} eq "Complete";
  };
  my %result;
  for my $study_id (keys %runs_per_study){
    my @runs = @{$runs_per_study{$study_id}};
    my ($rnaseqer_last_update, $assembly_used) = date_and_value_for_newest_key(map {
       $_->{LAST_PROCESSED_DATE} => $_->{ASSEMBLY_USED}
    } @runs);
    if ($rnaseqer_last_update lt "2019-07-15" && $rnaseqer_last_update gt "2019-04-15" && grep {$study_id eq $_} qw/SRP174213 SRP175031 SRP140708 SRP179824 SRP140458 SRP187776/){
       # On 2019-07-08 we've copied all data for which the last update date is after 2019-04-15, but still has htseq rather than feature counts.
       # If it's year 2020, check the two studies and remove this.
       $rnaseqer_last_update = "2019-04-01";
    }
    my %location_by_run;
    my %quality_by_run;
    for my $run (@runs){
#### $run
      my $run_id = $run->{RUN_IDS};
    
# the JSON retrieved just above may contain URLs that don't exist
# e.g. ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR100/0045/SRR10003845
# this appears to be because '0045' is used in the path where it should be '045'
# => try a HEAD request, and if it fails, attempt a workaround by removing a leading '0' from a directory in the path
{
   my $bwurl = $run->{BIGWIG_LOCATION};
   my $bwresponse = LWP::UserAgent->new->head($bwurl);
   if('404' eq $bwresponse->code()) {
      my $olddir = dirname $bwurl;
      $log->warn("JSON retrieved from https://www.ebi.ac.uk/fg/rnaseq/api/json/20/getRunsByOrganism/$species referes to $olddir which doesn't exist!");
      $bwurl =~ s~(/rnaseq/[A-Z]+\d+/)00(\d\d+)(/[A-Z]+\d+/)~${1}0${2}${3}~;
      my $newdir = dirname $bwurl;
      $bwresponse = LWP::UserAgent->new->head($bwurl);
      if('200' eq $bwresponse->code()) {
         $log->warn("$newdir appears to exist, so will try using this instead");
         $run->{BIGWIG_LOCATION} = $bwurl;
      } else {
         # we know the URL is bad, but it's not clear if data from this location will be required; so it's not yet a fatal error but requires a warning.
         $log->warn("No alternative to $olddir found;  if data from here are required, an error will occur.");
      }
   }
}
      
      $location_by_run{$run_id} = ($run->{ASSEMBLY_USED} //"") eq $assembly_used ? dirname ($run->{BIGWIG_LOCATION}) : "";
      $quality_by_run{$run_id} = $run->{MAPPING_QUALITY};
    }
    my %runs_by_sample;
    for (@runs){
       push @{$runs_by_sample{$_->{SAMPLE_IDS}}} , $_->{RUN_IDS};
    }
    my %replicates_by_run;
    for my $sample_id (keys %runs_by_sample){
       my @run_ids_for_sample_id = @{$runs_by_sample{$sample_id}};
       my $sample_has_much_more_than_half_ids = (scalar @run_ids_for_sample_id - 1) * 2 > scalar @runs;
       for my $run_id (@run_ids_for_sample_id){
          $replicates_by_run{$run_id} = (( $sample_has_much_more_than_half_ids ) || 1 == @run_ids_for_sample_id ? $run_id : $sample_id);
       }
    }
    $result{$study_id} = {
       assembly_used => $assembly_used,
       rnaseqer_last_update => $rnaseqer_last_update,
       location_by_run=> \%location_by_run,
       quality_by_run=> \%quality_by_run,
       replicates_by_run => \%replicates_by_run,
    };
  }
  return \%result;
}


sub date_and_value_for_newest_key {
  my %h = @_;
  my $format_rnaseqer = DateTime::Format::Strptime->new(pattern=> "%a %b %e %Y %T", strict=>1); # Fri Jun 19 2015 18:20:10
  my $format_iso = DateTime::Format::ISO8601::Format->new;
  my ($k, $v, @xs) = unpairs sort {DateTime->compare($b->[0], $a->[0])} map {[$format_rnaseqer->parse_datetime($_->[0]) , $_->[1]]} pairs %h;
  return $format_iso->format_date($k), ($v//"");
}
1;
