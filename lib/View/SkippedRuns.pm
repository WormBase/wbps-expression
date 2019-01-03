use strict;
use warnings;
package View::SkippedRuns;
use Text::MultiMarkdown qw/markdown/;
use View::Design;
sub new {
  my ($class, $study) = @_;
  return bless {study => $study}, $class;
}
# I'm copypasted from View::Study
# Extract the part drawing configs into a common UI?
# Could also have: links to ENA and possibly other places
# Could also have: run descriptions
sub to_markdown {
   my ($study) = @_;
   my $result = "";
   open (my $fh, ">:utf8", \$result);
   print $fh sprintf("### %s: %s\n", $study->{study_id}, $study->{config}{title} // "<title>" );
   print $fh sprintf("**%s**\n\n", $study->{config}{submitting_centre}) if $study->{config}{submitting_centre};
   print $fh sprintf("*%s*\n\n", $study->{config}{description}) if $study->{config}{description} and $study->{config}{description} ne $study->{config}{title};
   while (my ($k, $v) = each %{$study->{config}{pubmed} //{}}){
     print $fh sprintf("[%s](https://www.ncbi.nlm.nih.gov/pubmed/%s)\n\n", $v->[1], $k);
   }
   close $fh;
   return $result;
}
sub to_html {
  my ($self) = @_;
  return (%{$self->{study}{config}//{}} ? markdown(to_markdown($self->{study})) : sprintf("<h3>%s</h3>\n", $self->{study}{study_id}))  
     . "\n<h5>Runs</h5>"
     . "\n<ul>"
     . join ("\n", map {"<li>$_</li>"} @{$self->{study}{runs}})
     . "\n</ul>";
}
1;

