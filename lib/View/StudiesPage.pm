use strict;
use warnings;
package View::StudiesPage;
use View::Study;
sub new {
  my ($class, $species, @studies) = @_;
  return bless {species => $species, studies => \@studies}, $class;
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

  my $toc = join "\n", map {"<li>$_</li>\n"} map {$_->{study_id}} @{$self->{studies}};
  print $fh "<h3>TOC</h3>\n<ul>\n$toc</ul>\n";
  
  print $fh "<h3>Studies</h3>\n";
  for my $study (@{$self->{studies}}){
     print $fh View::Study->new($study)->to_html . "\n";
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
