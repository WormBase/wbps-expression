use strict;
use warnings;
package Production::Analysis::Common;

sub write_named_hashes {
  my ($name_to_data_pairs, $out_path, @frontmatter) = @_;
  print STDERR sprintf("write_named_hashes %s -> %s\n", scalar @{$name_to_data_pairs}, $out_path) if $ENV{ANALYSIS_VERBOSE};
  my %row_labels;
  for my $p (@{$name_to_data_pairs}){
    for my $label(keys %{$p->[1]}){
      $row_labels{$label}++;
    }
  }
  open(my $fh, ">", $out_path) or die "$out_path: $!";
  print $fh "# $_\n" for @frontmatter;
  print $fh join ("\t", "", map {$_->[0]} @{$name_to_data_pairs})."\n";
  for my $row (sort keys %row_labels){
     print $fh join ("\t",$row, map {$_->[1]->{$row} // ""} @{$name_to_data_pairs})."\n";
  }
  close $fh;
}
1;
