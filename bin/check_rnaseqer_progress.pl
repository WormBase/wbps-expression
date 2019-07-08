#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
# Our code
use lib "$FindBin::Bin/../lib";
# Local modules - installed with `carton`
use lib "$FindBin::Bin/../local/lib/perl5";

use List::MoreUtils qw/uniq/;
use JSON;
use LWP;
use feature 'say';
use ProductionMysql;
use DateTime;
use DateTime::Format::Strptime;
use DateTime::Format::ISO8601::Format;
my $format_rnaseqer = DateTime::Format::Strptime->new(pattern=> "%a %b %e %Y %T", strict=>1); # Fri Jun 19 2015 18:20:10
my $format_iso = DateTime::Format::ISO8601::Format->new;

my $most_recent_complete_date;
my $most_recent_complete_message;

my $runs_complete;
my %studies_complete;
my %species_complete;

my $num_failed;

my %in_progress;

my %not_yet_on_ftp;
for my $species (uniq map {/([a-z]*)_([a-z]*)/ ? "$1_$2" : ()} ProductionMysql->staging->species(@ARGV ? @ARGV : ("core_$ENV{PARASITE_VERSION}"))){
   my $url = "https://www.ebi.ac.uk/fg/rnaseq/api/json/0/getRunsByOrganism/$species"; 
   my $response = LWP::UserAgent->new->get($url);
   die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;

  for (@{decode_json($response->decoded_content) // []}){
    my $date = $format_rnaseqer->parse_datetime($_->{LAST_PROCESSED_DATE});
    my $date_str = $format_iso->format_date($date);
    if ($_->{STATUS} eq "Complete"){
      $runs_complete++;
      $studies_complete{$_->{STUDY_ID}}++;
      $species_complete{$species}++;
      if (not $most_recent_complete_date or $most_recent_complete_date < $date){
         $most_recent_complete_date = $date;
         $most_recent_complete_message = sprintf("%s study %s, run %s - %s", $species, $_->{STUDY_ID}, $_->{RUN_IDS}, $_->{LAST_PROCESSED_DATE});
      }
    } elsif($_->{STATUS} eq "Not_yet_on_ftp"){
      push @{$not_yet_on_ftp{$date_str}{$species}{$_->{STUDY_ID}}}, $_->{RUN_IDS};
    } elsif ($_->{STATUS} eq "In_progress") {
      push @{$in_progress{$species}{$_->{STUDY_ID}}}, $_->{RUN_IDS};
    } elsif ($_->{STATUS} =~ /failed/) {
      $num_failed++;
    } else {
      die "Have not seen this status for WormBase ParaSite data before: " . JSON->new->pretty(1)->encode($_);
    }
  }
}

my $studies_complete = keys %studies_complete;
my $species_complete = keys %species_complete;
say "Complete: $runs_complete runs among $studies_complete studies in $species_complete species" if $runs_complete;
say "Most recently completed: $most_recent_complete_message" if $most_recent_complete_message;
say "Failed: $num_failed runs" if $num_failed;
if(%in_progress){
  say "In progress: ";
  for my $species (sort keys %in_progress){
    say "  $species\t". join("; ", map {
      $_ . ": ".join (", ", @{$in_progress{$species}{$_}})
    }  keys %{$in_progress{$species}});
  }
}
if (%not_yet_on_ftp){
  say "Not yet on ftp:";
  for my $date_str (sort keys %not_yet_on_ftp){
    say "  $date_str\t" . join ("; ", map {
       my $species = $_;
       my $studies = keys %{$not_yet_on_ftp{$date_str}{$species}};
      "$species: $studies"
    } sort keys %{$not_yet_on_ftp{$date_str}});
  }
}

