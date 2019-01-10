use strict;
use warnings;

sub new {
  my ($class, $design) = @_;
  return bless {design => $design}, $class; 
}


1;
