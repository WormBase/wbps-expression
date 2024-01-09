#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use JSON;
# Our code
use lib "$FindBin::Bin/../lib";
# Local modules - installed with `carton`
use lib "$FindBin::Bin/../local/lib/perl5";
# Local R installation - downloaded the long way from https://www.r-project.org/ 
# Also, DESeq2 installed the long way from Bioconductor
# Check this: `R --slave --no-restore --file=- <<< 'installed.packages()' | grep ^DESeq2`
# Statistics::R depends on $PATH to find the R binary
$ENV{PATH} = "$FindBin::Bin/../local/R-3.6.3/bin:$ENV{PATH}";
use WbpsExpression;
$SIG{__WARN__} = sub { die @_ } if $ENV{DIE_ON_WARNINGS};

use Log::Any::Adapter;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
Log::Any::Adapter->set('Log4perl');

### This part is WormBase ParaSite specific
use ProductionMysql;
my @core_dbs = ProductionMysql->staging->core_databases(@ARGV);
my $work_dir = "$ENV{PARASITE_SCRATCH}/jbrowse/WBPS$ENV{PARASITE_VERSION}/WbpsExpression";
my $assemblies_rename_json = "$ENV{PARASITE_CONF}/brc4_rnaseq.assemblies-rename.json";
my $studies_dir = "$ENV{EXPRESSION_CODE}/studies";
my $prev_release_dir = "$ENV{PARASITE_SCRATCH}/old_jbrowse/WbpsExpression/";
my $brc4_dir = "$ENV{PARASITE_SCRATCH}/brc4rnaseq/WBPS$ENV{PARASITE_VERSION}/gene_expression";

my $json_data;
{
    local $/; # enable slurp mode
    open my $fh, "<", $assemblies_rename_json
        or die "could not open $assemblies_rename_json: $!";
    $json_data = <$fh>;
    close $fh;
}
my $assembly_rename_dict = decode_json($json_data);

if (@core_dbs > 5) {
  $ENV{DO_THROTTLE_GEO} //=1;
}



my $greenlight = 1;
for my $core_db (@core_dbs) {
  my ($spe, $cies, $bp ) = split "_", $core_db;
  # if ($core_db eq "ancylostoma_duodenale_prjna72581_core_17_105_1") {$greenlight = 0};
  # next if $bp eq 'core';
  # next if $greenlight == 0;
  # print $core_db."\n";
  my $wbps_assembly = ProductionMysql->staging->meta_value($core_db, "assembly.default");
  my $assembly = $wbps_assembly;
  $wbps_assembly =~ s/$_/$assembly_rename_dict->{$_}/g for keys %$assembly_rename_dict;
  $assembly =~ s/N__americanus_v1/N_americanus_v1/;
  $assembly =~ s/A_simplex_0011_upd/A_simplex_v1_5_4/;
  $assembly =~ s/S_solidus_NST_G2_0011_upd/S_solidus_NST_G2_v1_5_4/;
  $assembly =~ s/T_regenti_v1_0_4_001_upd/T_regenti_v1_0_4/;
  $assembly =~ s/T_canis_Equador_0011_upd/T_canis_Equador_v1_5_4/;
  $assembly =~ s/N_brasiliensis_RM07_v1_5_4_0011_upd/N_brasiliensis_RM07_v1_5_4/;
  $assembly = "" if $assembly eq "Heterorhabditis_bacteriophora-7.0";
  my $species = join ("_", $spe, $cies, $bp);
  WbpsExpression::run_brc4("$species", $assembly, $wbps_assembly, "$work_dir/$species", $studies_dir, "$prev_release_dir/$species", $brc4_dir);
  # $ENV{DO_DEPLOY_WEB} and print `sudo -u wormbase rsync --delete -av $work_dir/$species/  $ENV{PARASITE_FTP}/web_data/rnaseq_studies/releases/next/$species/`;
}
