use strict;
use warnings;
package View::StudiesPage;
use View::Study;
use View::SkippedRuns;
sub new {
  my ($class, $species, $studies) = @_;
  return bless {species => $species, studies => $studies}, $class;
}

sub to_html {
  my ($self) = @_;
  my $result = "";
  open(my $fh, ">:utf8", \$result);
  print $fh sprintf("<h2> %s - public RNASeq studies</h2>\n", do {
    my $species = $self->{species};
    $species =~ s/_/ /g;
    ucfirst($species)
  });

  
  print $fh "<h3>Studies</h3>\n";
  print $fh "<h4>Analysed</h4>\n";
  for my $study (@{$self->{studies}{passing_checks}}){
     print $fh View::Study->new($study)->to_html . "\n";
  }
  for my $study (@{$self->{studies}{failing_checks}}){
     print $fh View::Study->new($study)->to_html . "\n";
  }
  print $fh "<h4>Other</h4>\n";
  for my $skipped_runs (@{$self->{studies}{skipped_runs}}){
     print $fh View::SkippedRuns->new($skipped_runs)->to_html . "\n";
  }
  close $fh; 
  return sprintf('
   <html>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css" />
    <body>
      %s
    </body>
   </html>
',
   $result);
}
1;
