use strict;
use warnings;
use feature 'state';
package WbpsExpression::Analysis::DESeq2;
use Statistics::R; 
use WbpsExpression::Analysis::Common;
use List::MoreUtils qw/zip zip_unflatten/;
use Log::Any '$log';
#use Smart::Comments '###';

sub R_dds_for_design_path_and_counts_path {
  my ($R, $design_path, $counts_path) = @_;
  $R->set('colPath', $design_path);
  $log->info("Reading design $design_path");
  $R->run(q`colData = read.csv(colPath, header=TRUE, sep="\t", row.names=1)`);
  $R->set('countsPath', $counts_path);
  $log->info("Reading counts $counts_path");
  $R->run(q`countData = read.csv(countsPath, header=TRUE, sep="\t", row.names=1)`);
  $log->info("DESeqDataSetFromMatrix");
  $R->run(q`dataSet = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Condition)`);
# If there are no samples, do not collapse the replicates
  $R->run(q`if (exists("dataSet.Sample")) { collapseReplicates(dataSet, dataSet.Sample) }`);
  $log->info("DESeq");
  $R->run(q`dds = DESeq(dataSet)`);
}

sub R_results_for_reference_and_test {
  my ($R, $reference, $test) = @_;
  $R->set('reference', $reference);
  $R->set('test', $test);
  $R->run(q`res = results(dds, contrast = c("Condition", reference, test), alpha = 0.05)`);
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
  my @values = map {join(" ", @$_)} zip_unflatten @fold_changes, @adjusted_p_values;
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

sub values_for_contrast {
  my ($design_path, $counts_path, $reference, $test) = @_;
  my $attempts_to_get_r_for_this_contrast;
  state $current_dds;

  TRY_REUSE_DDS:
  if($current_dds && $current_dds eq "$design_path.$counts_path"){
    goto SET_RESULTS;
  }
  $current_dds = "$design_path.$counts_path";

  GET_R_SESSION:
  die "values_for_contrast $current_dds $reference $test: R doesn't seem to work today" if $attempts_to_get_r_for_this_contrast++ > 3;
  $log->info("starting new R session with DESeq2");
  $R_GLOBAL = Statistics::R->new();
  $R_GLOBAL->run(q`suppressPackageStartupMessages(library(DESeq2))`);
  $r_version //= $R_GLOBAL->get('getRversion()');
  $deseq_version //= $R_GLOBAL->get('packageVersion("DESeq2")');

  LOAD_DDS:
  unless (R_has_objects($R_GLOBAL, "DESeqDataSetFromMatrix", "collapseReplicates", "DESeq")){
     goto GET_R_SESSION;
  }
  $log->info("R_dds_for_design_path_and_counts_path");
  R_dds_for_design_path_and_counts_path($R_GLOBAL, $design_path, $counts_path);

  SET_RESULTS:
  unless (R_has_objects($R_GLOBAL, "results", "dds")){
     goto GET_R_SESSION;
  }
  $log->info("R_results_for_reference_and_test $reference $test");
  R_results_for_reference_and_test($R_GLOBAL, $reference, $test);

  GET_VALUES:
  unless (R_has_objects($R_GLOBAL, "res")){
     goto GET_R_SESSION;
  }
  return R_get_values_as_hash($R_GLOBAL);
}
sub do_analysis {
  my ($design_path, $counts_path, $contrasts, $out_path, @frontmatter) = @_;
  $log->info("DESeq2::do_analysis $out_path");
  my @name_to_data_pairs;
  my @analysis_warnings;
  for my $contrast (@{$contrasts}){
     my ($reference, $test, $contrast_name) = @{$contrast};
     my $h = values_for_contrast($design_path, $counts_path,$reference, $test);
     if(%$h){
       push @name_to_data_pairs, [$contrast_name, $h];
     } else {
       push @analysis_warnings, "No differentially expressed genes found in contrast $contrast_name";
     }
  }
  WbpsExpression::Analysis::Common::write_named_hashes(\@name_to_data_pairs, $out_path,
    shift @frontmatter,
    "R version: $r_version", "DESeq2 version: $deseq_version",
    @frontmatter, @analysis_warnings,
  );
}
1;
