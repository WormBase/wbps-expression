#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw/make_path/;
use FindBin;
use lib "$FindBin::Bin/../lib";
use WbpsExpression::Model::Design;
my %chs;
my %samples;
my %conditions;
while(<>){
  chomp;
  my ($run, $sample, $study, $life_stage, $sex, $species, $strain, $pmid, $organism_part) = split "\t";
  $species =~ s/ /_/;
  $species = lc $species;
  $life_stage = "adult" if $life_stage =~ /adult/;
  $chs{$species}{$study}{$run}{developmental_stage} = $life_stage;
  $chs{$species}{$study}{$run}{sex} = $sex;
  $chs{$species}{$study}{$run}{strain} = $strain;
  $chs{$species}{$study}{$run}{organism_part} = $organism_part;
  $samples{$species}{$study}{$run} = $sample;
  $conditions{$species}{$study}{$run} = "";
}
my @characteristics_in_order = ("developmental_stage", "sex", "strain", "organism_part");
for my $species (keys %chs) {
  for my $study (keys %{$chs{$species}}){
    make_path "curation/studies/$species/$study";
    my $design = WbpsExpression::Model::Design::from_data_by_run(
       $samples{$species}{$study},
       $conditions{$species}{$study},
       $chs{$species}{$study},
       \@characteristics_in_order,
    );
    $design->to_tsv("curation/studies/$species/$study/$study.tsv");
  }
}
