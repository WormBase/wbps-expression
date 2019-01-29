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

# Turn on logs
$ENV{ANALYSIS_VERBOSE} //= 1;
$ENV{LOCALLY_CACHED_RESOURCE_VERBOSE} //= 1;

use ProductionMysql;
my @core_dbs = ProductionMysql->staging->core_databases(@ARGV);

use Production::Workflow;
my $data_dir = "/nfs/nobackup/ensemblgenomes/wormbase/parasite/production/jbrowse/WBPS$ENV{PARASITE_VERSION}";
my $src_dir = "$FindBin::Bin/..";
my $work_dir = "$data_dir/Production-".`whoami`;
chomp $work_dir;

my $processing = Production::Workflow->new("$data_dir/Resources",$src_dir, $work_dir );
for my $core_db (@core_dbs) {
  my ($spe, $cies, $bp ) = split "_", $core_db;
  next if $bp eq 'core';
  next unless "$data_dir/Resources/${spe}_${cies}";
  my $assembly = ProductionMysql->staging->meta_value($core_db, "assembly.name");
  $processing->do_everything ("${spe}_${cies}", $assembly);
  $ENV{DO_DEPLOY_WEB} and print `sudo -u wormbase rsync --delete -av $work_dir/${spe}_${cies}/$assembly/  /ebi/ftp/pub/databases/wormbase/parasite/web_data/rnaseq_studies/releases/next/${spe}_${cies}_${bp}/`;
} 
