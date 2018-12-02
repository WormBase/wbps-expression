use strict;
use warnings;

package Production::Analysis::DESeq2;


sub new {
  my ($class, $rscript_path) = @_;
  $rscript_path //= `which Rscript`;
  die "Please put Rscript in PATH" unless $rscript_path;
  bless {
    rscript => $rscript_path,
  }, $class;
}
my $r_program = <<'EOF';
suppressMessages( library( DESeq2 ) )
args = commandArgs(trailingOnly=TRUE)
colData = read.csv(args[0], header=TRUE, sep="\t", row.names=1)
countData = read.csv(args[1], header=TRUE, sep="\t", row.names=1)
formula = args[2]
outPath = args[3]
contrasts = args[4:] # TODO

dds = DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = formula)
dds = DESeq(dds)
allResults = *new data frame*
for c in contrasts:
  res = results(dds, contrast = c, alpha = 0.05)
  res2 = subset(res, padj < 0.05)
  res3 = as.data.frame(res2$log2FoldChange)
  rownames(res3) = rownames(res2)
  # Pad with zeros?
  addResultTODO(allResults, res3, c)

write.csv(allResults, file=outPath)
EOF

sub do_deseq2 {
  my ($self, $design_path, $counts_path, $batch_effects, $contrasts, $out_path) = @_;
  my $design_formula = "~ " . join ("+", @{$batch_effects}, "Condition");
  system("<<<", $r_program, $self{rscript}, "--vanilla", "-", 
      $design_path, $counts_path, $design_formula, $out_path, @{$contrasts}, 
  ) and do {
     unlink $out_path;
     die "Failed to create $out_path: R exitted with $?";
  };
  die "Failed to create $out_path: file not found" unless -f $out_path;
  
}
