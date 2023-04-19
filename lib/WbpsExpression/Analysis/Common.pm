use strict;
use warnings;
package WbpsExpression::Analysis::Common;
use Log::Any '$log';

sub write_frontmatter_only {
  my($out_path, @frontmatter) = @_;
  open(my $fh, ">", $out_path) or die "$out_path: $!";
  print $fh "# $_\n" for @frontmatter;
  close $fh;
}
sub write_named_hashes {
  my ($name_to_data_pairs, $out_path, @frontmatter) = @_;
  $log->info(sprintf("write_named_hashes %s -> %s\n", scalar @{$name_to_data_pairs}, $out_path));
  my %row_labels;
  for my $p (@{$name_to_data_pairs}){
    for my $label(keys %{$p->[1]}){
      $row_labels{$label}++;
    }
  }
  open(my $fh, ">:utf8", $out_path) or die "$out_path: $!";
  print $fh "# $_\n" for @frontmatter;
  if(@{$name_to_data_pairs}){
    print $fh join ("\t", "gene_id", map {$_->[0]} @{$name_to_data_pairs})."\n";
    for my $row (sort keys %row_labels){
       print $fh join ("\t",$row, map {$_->[1]->{$row} // 0} @{$name_to_data_pairs})."\n";
    }
  }
  close $fh;
}
1;
