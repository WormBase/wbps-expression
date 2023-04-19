use strict;
use warnings;
package WbpsExpression::StudiesPage;
use open ':std', ':encoding(UTF-8)';
use File::Basename;
use HTML::Template;
use Data::Dumper;

sub template {
  my ($name) = @_;
  return HTML::Template->new(path => [dirname (__FILE__) ], filename => "StudiesPage/$name.tmpl");
}

sub html {
  my ($species, $our_studies, $other_studies) = @_;
  
  my $template = template("studies");

  (my $species_param = $species) =~ s/_/ /g;
  $species_param = ucfirst($species_param);
  $template->param(SPECIES => $species_param);
  $template->param(MORE_WELCOMING_MESSAGE_FOR_SKIPPED_STUDIES => ($species_param =~ /Caenorhabditis elegans|Pristionchus pacificus|Globodera pallida/));

  my @toc = map {{
    TOC_ITEM_ID => $_->{study_id},
    TOC_ITEM_NAME => $_->{config}{title},
  }} @{$our_studies};
  push @toc, {
    TOC_ITEM_ID => "wbps_expression_other",
    TOC_ITEM_NAME => @{$other_studies} == 1 ? "1 other study" : sprintf("%s other studies", scalar @{$other_studies}),
  } if @{$other_studies};
  $template->param(TOC => \@toc);

  $template->param(RUNSTUDIES => join("\n", map {study_html($_)} @{$our_studies} ));
  $template->param(OTHERSTUDIES => join ("\n", map {other_study_html($_)} @{$other_studies}));

  return $template->output;
}

sub description_html {
  my ($study) = @_;

  my @description;
  if ($study->{config}{description} and $study->{config}{description} ne $study->{config}{title}){
      push @description, sprintf("<em>%s</em>", $study->{config}{description});
  }
  for my $pubmed_id (sort keys %{$study->{config}{pubmed}}){
     push @description, sprintf ("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/%s\">%s</a> ", $pubmed_id, $study->{config}{pubmed}{$pubmed_id}[1]);
  }
  for my $resource_link (@{$study->{config}{resource_links}}){
     my ($resource_type, $resource_title, $resource_url) = @{$resource_link};
     push @description, sprintf ("<a href=\"%s\">%s</a>",$resource_url, $resource_title);
  }
  return join("<br>\n", @description);
}

sub downloads_html {
  my ($study) = @_;
  my $template = template("downloads");
  my @analyses;
  my $fails_checks = $study->passes_checks;
  for my $analysis($study->analyses_required){
    my %item;
    if($fails_checks){
        $item{FAILEDANALYSIS} = $analysis->{title};
        push(@analyses, \%item);
    } else {
        $item{ANALYSIS}     = $analysis->{title};
        if ($analysis->{file_name}){
            $item{ANALYSISLINK} = join("/", $study->{study_id}, $analysis->{file_name});
        } else {
            $item{ANALYSISLINKS} = [map {
                {SUBANALYSIS=> $_->{title}, SUBANALYSISLINK => join("/", $study->{study_id}, $_->{file_name})}
            } @{$analysis->{files}}];
        }
        push(@analyses, \%item);
    }
  }
  if($fails_checks) {    
    $template->param(FAILEDANALYSES => \@analyses);
  } else {
    $template->param(ANALYSES => \@analyses);
  }
  return $template->output;
}
sub design_html {
  my ($design) = @_;
  my $template = template("design");
  my $design_summary = sprintf("Design: %s conditions\n", scalar $design->all_conditions);
  $template->param(DESIGNSUMMARY => $design_summary);

  my @a = ('Condition', 'num. replicates');
  push (@a, map { ucfirst $_ } @{$design->{characteristics_in_order}});
  my @columns_name = map {
     my $label = $_;
     $label =~ s/_/ /g;
     $label =~ s/Rnai/RNAi/;
     {'COLUMN' => $label }
  } @a;

  $template->param(DESIGNCOLUMNS => \@columns_name);
  
  my @rows;
  my %h = %{$design->runs_by_condition_then_replicate};
  for my $c (sort keys %h){
     my @b = ($c);
     push @b, scalar keys %{$h{$c}};
     push @b, map {$design->value_in_condition($c, $_) // "" } @{$design->{characteristics_in_order}};
     my @row;
     push (@row, map {{'ROW' => $_}} @b);
     push (@rows, {'DESIGNROW' => \@row});
  }   
  $template->param(DESIGNROWS => \@rows);
  return $template->output;
} 
sub study_html {
  my ($study) = @_;
  my $template = template("study");
  $template->param(STUDYID =>  $study->{study_id});
  $template->param(STUDYTITLE =>  $study->{config}{title});
  $template->param(STUDYCENTRE => $study->{config}{submitting_centre}) if $study->{config}{submitting_centre};
  $template->param(STUDYDESCRIPTION => description_html($study));
  $template->param(STUDYDOWNLOADS => downloads_html($study));
  $template->param(STUDYDESIGN => design_html($study->{design}));
  return $template->output;   
}

sub other_study_html {
  my ($study) = @_;
  my $template = template("skip");
  $template->param(STUDYID =>  $study->{study_id});
  $template->param(STUDYTITLE =>  $study->{config}{title} // "" );
  $template->param(STUDYCENTRE => $study->{config}{submitting_centre}) if $study->{config}{submitting_centre};
  $template->param(STUDYDESCRIPTION => description_html($study));
  my $num_runs = scalar keys %{$study->{sources}};
  $template->param(RUNS => $num_runs eq 1 ? "1 run in ENA": "$num_runs runs in ENA");
  
  return $template->output;
}

1;

