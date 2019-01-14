#!/usr/bin/env perl
use File::Find;
use File::Slurp qw/read_file/;
use FindBin;
#unless (-d "$FindBin::Bin/../local"){
#  print STDERR `carton`;
#}
#use lib "$FindBin::Bin/../local";
use lib "$FindBin::Bin/../lib";
use View::StudiesPage;
use Model::Study;
use Model::SkippedRuns;

my ($species) = @ARGV;
die "Usage: $0 species" unless $species;

my @studies;
find(
  sub {
    push @studies, Model::Study->from_folder($File::Find::name)
      if -d $File::Find::name and $File::Find::name =~ /\w+\d+$/;
  },
  "$FindBin::Bin/../curation/studies/$species"
);
my @skipped;
find(
  sub {
    push @skipped, Model::SkippedRuns->from_folder($File::Find::name)
      if -d $File::Find::name and $File::Find::name =~ /\w+\d+$/;
  },
  "$FindBin::Bin/../curation/skipped_runs/$species"
);


print View::StudiesPage->new($species, {
  passing_checks => \@studies,
  failing_checks => [],
  skipped_runs => \@skipped,
})->to_html;
