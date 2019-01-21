use strict;
use warnings;
package View::Study;
use HTML::Template;
use Data::Dumper;
use File::Basename;
sub new {
  my ($class, $study) = @_;
  return bless {study => $study}, $class;
}

sub to_html {
  my ($self) = @_;
  my $study = $self->{study};
  my $study_tmpl = HTML::Template->new(path => [dirname (__FILE__) ], filename => 'templates/study.tmpl');
   
   
  $study_tmpl->param(STUDYID =>  $study->{study_id});
  $study_tmpl->param(STUDYTITLE =>  $study->{config}{title} // "" );
  $study_tmpl->param(STUDYCENTRE => $study->{config}{submitting_centre}) if $study->{config}{submitting_centre};

  my @description;
  if ($study->{config}{description} and $study->{config}{description} ne $study->{config}{title}){
      push @description, $study->{config}{description};
  }
  while (my ($k, $v) = each %{$study->{config}{pubmed} //{}}){
     push @description, sprintf ("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/%s\">%s</a> ", $k, $v->[1]);
  }

  $study_tmpl->param(STUDYDESCRIPTION => join("\n<br>", @description));

  my @analyses;
  my $fails_checks = not $study->passes_checks;
  for my $analysis($study->analyses_required){
    my %item;
    if($fails_checks){
        $item{FAILEDANALYSIS} = $analysis->{title};
        push(@analyses, \%item);
    } else {
        $item{ANALYSIS}     = $analysis->{title};
        $item{ANALYSISLINK} = join("/", $study->{study_id}, $analysis->{file_name});
        push(@analyses, \%item);
    }
  }

  if($fails_checks) {    
    $study_tmpl->param(FAILEDANALYSES => \@analyses);
  } else {
    $study_tmpl->param(ANALYSES => \@analyses);
  }

  my $design = $study->{design};
  my $design_summary = sprintf("Design: %s conditions\n", scalar $design->all_conditions);
  $study_tmpl->param(DESIGNSUMMARY => $design_summary);

  my @a = ('Condition', 'num. replicates');
  push (@a, map { ucfirst $_ } @{$design->{characteristics_in_order}});
  my @columns_name = map {
     my $label = $_;
     $label =~ s/_/ /;
     # Also: ucfirst?
     {'COLUMN' => $label }
  } @a;

  $study_tmpl->param(DESIGNCOLUMNS => \@columns_name);
  
  my @rows;
  my %h = %{$design->runs_by_condition_then_replicate};
  for my $c (sort keys %h){
     my @b = ($c);
     push @b, scalar keys $h{$c};
     push @b, map {$design->value_in_condition($c, $_) // "" } @{$design->{characteristics_in_order}};
     my @row;
     push (@row, map {{'ROW' => $_}} @b);
     push (@rows, {'DESIGNROW' => \@row});
  }   
  $study_tmpl->param(DESIGNROWS => \@rows);
  
  return $study_tmpl->output;   
}

1;

