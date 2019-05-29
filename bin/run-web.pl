#!/usr/bin/env perl
# Create a studies page from curation data in the repository
use FindBin;
#unless (-d "$FindBin::Bin/../local"){
#  print STDERR `carton`;
#}
use lib "$FindBin::Bin/../local/lib/perl5";
use lib "$FindBin::Bin/../lib";
use WbpsExpression;

my ($species, $output_dir) = @ARGV;
die "Usage: $0 species output_dir" unless $species and $output_dir;

WbpsExpression::run_web_only($species, $output_dir);
