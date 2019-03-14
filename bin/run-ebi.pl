#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
# Our code
use lib "$FindBin::Bin/../lib";
# Local modules - installed with `carton`
use lib "$FindBin::Bin/../local/lib/perl5";
# Local R installation - downloaded the long way from https://www.r-project.org/ 
# Also, DESeq2 installed the long way from Bioconductor
# Check this: `R --slave --no-restore --file=- <<< 'installed.packages()' | grep ^DESeq2`
# Statistics::R depends on $PATH to find the R binary
$ENV{PATH} = "$FindBin::Bin/../local/R-3.5.1/bin:$ENV{PATH}";
use WbpsExpression;

use Log::Any::Adapter;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
Log::Any::Adapter->set('Log4perl');

### This part is WormBase ParaSite specific
use ProductionMysql;
my @core_dbs = ProductionMysql->staging->core_databases(@ARGV);
my $data_dir = "/nfs/nobackup/ensemblgenomes/wormbase/parasite/production/jbrowse/WBPS$ENV{PARASITE_VERSION}";
my $src_dir = "$FindBin::Bin/..";
my $work_dir = "$data_dir/Production-".`whoami`;
chomp $work_dir;

my $processing = WbpsExpression->new("$data_dir/Resources",$src_dir);
for my $core_db (@core_dbs) {
  my ($spe, $cies, $bp ) = split "_", $core_db;
  next if $bp eq 'core';
  next unless "$data_dir/Resources/${spe}_${cies}";
  my $assembly = ProductionMysql->staging->meta_value($core_db, "assembly.name");
  my $species = join ("_", $spe, $cies, $bp);
  $processing->run("${spe}_${cies}", $assembly, "$work_dir/$species");
  $ENV{DO_DEPLOY_WEB} and print `sudo -u wormbase rsync --delete -av $work_dir/$species/  /ebi/ftp/pub/databases/wormbase/parasite/web_data/rnaseq_studies/releases/next/$species/`;
} 
