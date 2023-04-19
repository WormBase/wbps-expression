use strict;
use warnings;
use feature 'state';
package WbpsExpression::Analysis::DESeq2;
use Statistics::R; 
use WbpsExpression::Analysis::Common;
use List::MoreUtils qw/zip zip_unflatten/;
use Log::Any '$log';
use File::Basename;
#use Smart::Comments '###';

sub R_dds_for_design_path_and_counts_path {
  my ($R, $design_path, $counts_path) = @_;
  $R->set('colPath', $design_path);
  $R->run(q`colData = read.csv(colPath, header=TRUE, sep="\t", row.names=1)`);
  $R->set('countsPath', $counts_path);
  $R->run(q`countData = read.csv(countsPath, header=TRUE, sep="\t", row.names=1)`);
  $log->info("DESeqDataSetFromMatrix");
  $log->info($counts_path);
  $R->run(q`dataSet = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Condition)`);
# If there are no samples, do not collapse the replicates
  $R->run(q`if (exists("dataSet.Sample")) { collapseReplicates(dataSet, dataSet.Sample) }`);
  $R->run(q`dds = DESeq(dataSet)`);
}

sub R_results_for_reference_and_test {
  my ($R, $reference, $test) = @_;
  $R->set('reference', $reference);
  $R->set('test', $test);
  $R->run(q`res = results(dds, contrast = c("Condition", reference, test), alpha = 0.05)`);
}
sub R_append_unfiltered_results_to_file {
  my ($R, $path) = @_;
  $R->set('target_path', $path);
  $R->run(q`write.table(
     data.frame(
       gene_id=rownames(res),
       baseMean=round(res$baseMean, digits=2),
       log2FoldChange=round(res$log2FoldChange, digits=2),
       lfcSE=round(res$lfcSE, digits=2),
       stat=signif(res$stat, digits=2),
       pvalue=signif(res$pvalue, digits=2),
       padj=signif(res$padj, digits=2)
     ),
     file = target_path, sep = "\t", row.names = FALSE, quote=FALSE, na="NA", append=TRUE
  )`);
}
sub R_get_values_as_hash {
  my ($R) = @_;
  $R->run(
    q`res = subset(res, padj < 0.05)`,
    q`res = subset(res, abs(log2FoldChange) > 0.5 )`,
    q`rn = rownames(res)`,
    q`fc = res$log2FoldChange`,
    q`fc = signif(fc, digits=2)`,
    q`pv = res$padj`,
    q`pv = signif(pv, digits=2)`,
  );
  my $fc = $R->get('fc');
  my $rn = $R->get('rn');
  my $pv = $R->get('pv');
  my @fold_changes = ref $fc eq 'ARRAY' ? @{$fc} : ();
  my @row_names = ref $rn eq 'ARRAY' ? @{$rn} : ();
  my @adjusted_p_values = ref $pv eq 'ARRAY' ? @{$pv} : ();
  die unless @fold_changes == @row_names and @fold_changes == @adjusted_p_values;
  my @values = map {join("\t", @$_)} zip_unflatten @fold_changes, @adjusted_p_values;
  my %h = zip @row_names, @values;
  return \%h;
}
sub R_has_objects {
  my ($R, @object_names) = @_;
  return unless $R;
  return unless $R->is_started();
  for my $object_name (@object_names){
    my $o_exists = $R->run("exists(\"$object_name\")");
    return unless $o_exists =~ /TRUE/;
  }
  return 1;
}

our $r_version;
our $deseq_version;
our $R_GLOBAL;

sub establish_r_session_and_set_r_and_deseq_versions {
  state $attempts_to_get_r;
  die "R doesn't seem to work today" if $attempts_to_get_r++ > 3;
  $log->info(__PACKAGE__." starting a new R session with a DESeq2 package");
  eval {
    local $SIG{ALRM} = sub { die "GET_R_SESSION alarm\n" };
    alarm 180;
    $R_GLOBAL = Statistics::R->new();
    $R_GLOBAL->run(q`suppressPackageStartupMessages(library(DESeq2))`);
    $r_version //= ($R_GLOBAL->get('getRversion()') =~ s/[^\d\.]//rg);
    $deseq_version //= ($R_GLOBAL->get('packageVersion("DESeq2")') =~ s/[^\d\.]//rg);
    alarm 0;
  };
  if($@){
    die $@ unless $@ eq "GET_R_SESSION alarm\n"; 
    $R_GLOBAL->stop if $R_GLOBAL;
    establish_r_session_and_set_r_and_deseq_versions();
  }
  $attempts_to_get_r=0;
}

sub values_for_contrast {
  my ($design_path, $counts_path, $reference, $test, $unfiltered_results_path) = @_;
  state $current_dds;

  TRY_REUSE_DDS:
  if($current_dds && $current_dds eq "$design_path.$counts_path"){
    goto SET_RESULTS;
  }
  $current_dds = "$design_path.$counts_path";

  TRY_REUSE_R:
  if($R_GLOBAL && $R_GLOBAL->is_started){
    goto LOAD_DDS;
  }

  GET_R_SESSION:
  establish_r_session_and_set_r_and_deseq_versions();

  LOAD_DDS:
  unless (R_has_objects($R_GLOBAL, "DESeqDataSetFromMatrix", "collapseReplicates", "DESeq")){
     goto GET_R_SESSION;
  }
  $log->info(__PACKAGE__."::R_dds_for_design_path_and_counts_path");
  R_dds_for_design_path_and_counts_path($R_GLOBAL, $design_path, $counts_path);

  SET_RESULTS:
  unless (R_has_objects($R_GLOBAL, "results", "dds")){
     goto GET_R_SESSION;
  }
  $log->info(__PACKAGE__."::R_results_for_reference_and_test $reference $test");
  R_results_for_reference_and_test($R_GLOBAL, $reference, $test);

  GET_VALUES:
  unless (R_has_objects($R_GLOBAL, "res")){
     goto GET_R_SESSION;
  }
  R_append_unfiltered_results_to_file($R_GLOBAL, $unfiltered_results_path);
  return R_get_values_as_hash($R_GLOBAL);
}
sub do_analysis {
  my ($design_path, $counts_path, $unfiltered_results_folder, $contrasts, $output_paths, $frontmatters, $warnings_by_contrast_name_by_key) = @_;

  for my $key (keys %{$contrasts}){
  $log->info(__PACKAGE__."::do_analysis " . (basename (dirname $design_path)));
    do_analysis_for_key($design_path, $counts_path, $unfiltered_results_folder,
       $contrasts->{$key}, $output_paths->{$key}, $frontmatters->{$key}, $warnings_by_contrast_name_by_key->{$key}
    );
  }
}

sub contrast_name_to_file_name {
  my ($contrast_name) = @_;
# Inspired by:
# https://metacpan.org/release/Text-Slugify/source/lib/Text/Slugify.pm
  $contrast_name =~ s/[^a-z0-9]+/-/gi;
  $contrast_name =~ s/^-?(.+?)-?$/$1/;
  $contrast_name =~ s/^(.+)$/\L$1/;
  return $contrast_name;
}
sub do_analysis_for_key {
  my ($design_path, $counts_path, $unfiltered_results_folder,
     $contrasts, $out_path, $frontmatter, $warnings_by_contrast_name) = @_;
  my @frontmatter = @{$frontmatter};
  my @name_to_data_pairs;
  my @analysis_warnings;
  establish_r_session_and_set_r_and_deseq_versions();
  for my $contrast (@{$contrasts}){
     my ($reference, $test, $contrast_name) = @{$contrast};
     push @analysis_warnings, @{$warnings_by_contrast_name->{$contrast_name}} if $warnings_by_contrast_name->{$contrast_name};
     my $unfiltered_results_path = join ("/", $unfiltered_results_folder, contrast_name_to_file_name($contrast_name).".tsv");
     WbpsExpression::Analysis::Common::write_frontmatter_only($unfiltered_results_path,
        @frontmatter,
        "Using R version $r_version, DESeq2 version $deseq_version",
        "Results for contrast $contrast_name",
        "Values as in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#more-information-on-results-columns",
        @{$warnings_by_contrast_name->{$contrast_name}}
     );
     my $h = values_for_contrast($design_path, $counts_path,$reference, $test, $unfiltered_results_path);
     die "Missing: $unfiltered_results_path" unless -f $unfiltered_results_path;
     if(%$h){
       push @name_to_data_pairs, ["$contrast_name\t", $h];
     } else {
       push @analysis_warnings, "No differentially expressed genes found in contrast $contrast_name";
     }
  }
  WbpsExpression::Analysis::Common::write_named_hashes(\@name_to_data_pairs, $out_path,
    @frontmatter, 
      "Using R version $r_version, DESeq2 version $deseq_version",
      "Values are base 2 logarithm of maximum likelihood estimate of fold change and adjusted p-value by Wald test per gene for each contrast",
      "Values rounded to 2.s.f., filtered past a significance threshold of adj_pval < 0.05 and abs(log2fc) > 0.5",
    @analysis_warnings,
  );
}
1;
