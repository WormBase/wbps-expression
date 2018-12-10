use strict;
use warnings;
package View::Study;
use Text::MultiMarkdown qw/markdown/;
use View::Design;
sub new {
  my ($class, $study) = @_;
  return bless {study => $study}, $class;
}
sub to_markdown {
   my ($study) = @_;
   my $result = "";
   open (my $fh, ">", \$result);
   print $fh sprintf("### %s: %s\n", $study->{study_id}, $study->{config}{title} // "<title>" );
   print $fh sprintf("*%s*\n", $study->{config}{description}) if $study->{config}{description} and $study->{config}{description} ne $study->{config}{title};
   while (my ($k, $v) = each %{$study->{config}{pubmed} //{}}){
     print $fh sprintf("[%s](https://www.ncbi.nlm.nih.gov/pubmed/%s)\n", $v->[1], $k);
   }
   my $fails_checks = not $study->passes_checks;
   print $fh "*Not analysed - needs curation*\n" if $fails_checks;
   print $fh "#### Data\n";
   for my $analysis($study->analyses_required){
     print $fh sprintf("##### %s\n", $analysis->{title});
     if($fails_checks){
        print $fh sprintf( "<strike>%s</strike>\n", $analysis->{description});
     } else {
        print $fh  sprintf("[%s](%s/%s)\n", $analysis->{description}, $study->{study_id}, $analysis->{file_name});
     }
   }
  # print $fh $study->{design}->to_markdown;
   close $result;
   return $result;
}
sub to_html {
  my ($self) = @_;
  return markdown(to_markdown($self->{study})) . "\n" . View::Design->new($self->{study}->{design})->to_html;
}
1;

