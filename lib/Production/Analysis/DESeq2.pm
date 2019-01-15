use strict;
use warnings;
package Production::Analysis::DESeq2;
use Statistics::R; 
use Production::Analysis::Common;
use List::MoreUtils qw/zip/;
#use Smart::Comments '###';

sub do_analysis {
  my ($design_path, $counts_path, $contrasts, $out_path, @frontmatter) = @_;
  my @name_to_data_pairs;
  my @analysis_warnings;
  my $R = Statistics::R->new();
  my $r_version = $R->get('getRversion()');

  $R->set('colPath', $design_path);
  $R->run(q`colData = read.csv(colPath, header=TRUE, sep="\t", row.names=1)`);
  $R->set('countsPath', $counts_path);
  $R->run(q`countData = read.csv(countsPath, header=TRUE, sep="\t", row.names=1)`);

  $R->run(q`suppressPackageStartupMessages(library(DESeq2))`);
  my $deseq_version = $R->get('packageVersion("DESeq2")');
  $R->run(q`dataSet = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Condition)`);
# TODO this should handle the two-or-three runs case
# If there are no samples, do not collapse the replicates
  $R->run(q`collapseReplicates(dataSet, ~ Sample)`);
  $R->run(q`dds = DESeq(dataSet)`);
  for my $contrast (@{$contrasts}){
     my ($reference, $test, $contrast_name) = @{$contrast};
     $R->set('reference', $reference);
     $R->set('test', $test);
     $R->run(
      q`res = results(dds, contrast = c("Condition", reference, test), alpha = 0.05)`,
      q`res = subset(res, padj < 0.05)`,
      q`rn = rownames(res)`,
      q`fc = res$log2FoldChange`,
      q`fc = round(fc, digits=2)`,
     );
     my $fc = $R->get('fc');
     my $rn = $R->get('rn');
     my @fold_changes = ref $fc eq 'ARRAY' ? @{$fc} : ();
     my @row_names = ref $rn eq 'ARRAY' ? @{$rn} : ();
     die unless @fold_changes == @row_names;
     my %h = zip @row_names, @fold_changes;
     if(%h){
       push @name_to_data_pairs, [$contrast_name, \%h];
     } else {
       push @analysis_warnings, "No differentially expressed genes found in contrast: $contrast_name";
     }
  }
  $R->stop;
  Production::Analysis::Common::write_named_hashes(\@name_to_data_pairs, $out_path,
    shift @frontmatter,
    "R version: $r_version", "DESeq2 version: $deseq_version",
    @frontmatter, @analysis_warnings,
  );
}
1;
