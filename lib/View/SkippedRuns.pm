use strict;
use warnings;
package View::SkippedRuns;
use HTML::Template;
use File::Basename;
sub new {
  my ($class, $study) = @_;
  return bless {study => $study}, $class;
}

sub to_html {  
  my ($self) = @_;
  my $study = $self->{study};

  my $skip_tmpl = HTML::Template->new(path => [dirname (__FILE__) ], filename => 'templates/skip.tmpl');

  $skip_tmpl->param(STUDYID =>  $study->{study_id});
  $skip_tmpl->param(STUDYTITLE =>  $study->{config}{title} // "" );
  $skip_tmpl->param(STUDYCENTRE => $study->{config}{submitting_centre}) if $study->{config}{submitting_centre};
  

  my @description;
  if ($study->{config}{description} and $study->{config}{description} ne $study->{config}{title}){
      push @description, $study->{config}{description};
  }
  while (my ($k, $v) = each %{$study->{config}{pubmed} //{}}){
     push @description, sprintf ("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/%s\">%s</a> ", $k, $v->[1]);
  }

  $skip_tmpl->param(STUDYDESCRIPTION => join("\n<br>", @description));
  $skip_tmpl->param(NUMRUNS => scalar @{$study->{runs}});
  
  return $skip_tmpl->output;
}

1;

