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
  $study_tmpl->param(STUDYTITLE =>  $study->{config}{title} // "NO TITLE" );
  $study_tmpl->param(STUDYCENTRE => $study->{config}{submitting_centre}) if $study->{config}{submitting_centre};

  my $description = "";
  $description .= $study->{config}{description} if $study->{config}{description} and $study->{config}{description} ne $study->{config}{title};
  while (my ($k, $v) = each %{$study->{config}{pubmed} //{}}){
     $description .= sprintf ("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/%s\">%s</a> ", $k, $v->[1]);
  }
  $study_tmpl->param(STUDYDESCRIPTION => $description);

  my @analyses;
  my $fails_checks = not $study->passes_checks;
  for my $analysis($study->analyses_required){
    my %item;
    if($fails_checks){
        $item{FAILEDANALYSIS} = $analysis->{description};
        push(@analyses, \%item);
    } else {
        $item{ANALYSIS}     = $analysis->{description};
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
  my $design_summary = sprintf("Design: %s conditions across %s runs\n", scalar $design->all_conditions, scalar $design->all_runs);
  $study_tmpl->param(DESIGNSUMMARY => $design_summary);

  my @a = ('Run', 'Condition');
  push (@a, map { ucfirst $_ } @{$design->{characteristics_in_order}});
  my @columns_name = map {{'COLUMN' => ucfirst $_}} @a;
  $study_tmpl->param(DESIGNCOLUMNS => \@columns_name);
  

  my @rows;
  for my $p ($design->condition_run_ordered_pairs){
     my @b = ($p->[1], $p->[0]);
     push (@b, map {$design->value_in_run($p->[1], $_)} @{$design->{characteristics_in_order}});
     my @row;
     push (@row, map {{'ROW' => $_}} @b);
     push (@rows, {'DESIGNROW' => \@row});
  }   
  $study_tmpl->param(DESIGNROWS => \@rows);
  
  if ($study->{SKIPPED_RUNS}) {
    $study_tmpl->param(HASSKIPS => 1);
    $study_tmpl->param(SKIPPEDCOUNT => scalar @{$study->{SKIPPED_RUNS}});
    my $skr = join ', ' , @{$study->{SKIPPED_RUNS}};
    $study_tmpl->param(SKIPPEDRUNS => $skr);
  }
  
  return $study_tmpl->output;   
}

1;

