use strict;
use warnings;
package View::Design;
use Text::MultiMarkdown qw/markdown/;
sub new {
  my ($class, $design) = @_;
  return bless {design => $design}, $class; 
}

# TODO remove, create nice html instead
sub to_markdown {
  my ($design) = @_;
  my $table = ""; 
  open(my $fh, ">:utf8", \$table);
  my @a = @{$design->{characteristics_in_order}};
  print $fh join(" | ","" , "Run", "Condition", @a, "")."\n";
  print $fh join(" | ","", "---", "---",(map {"---"} @a), "")."\n";
  for my $p ($design->condition_run_ordered_pairs){
     print $fh join (" | ","", $p->[1], $p->[0], (map {$design->value_in_run($p->[1], $_)} @a), "")."\n";
  }
  close $fh;
  my $summary = sprintf("#### Design: %s conditions across %s runs\n", scalar $design->all_conditions, scalar $design->all_runs);
  return "$summary\n$table\n";
}

sub to_html {
  my ($self) = @_;
  return markdown(to_markdown($self->{design}));
}
1;
