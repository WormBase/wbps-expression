use strict;
use warnings;
package WbpsExpression::StudiesPage;
use WbpsExpression::StudiesPage::Study;
use WbpsExpression::StudiesPage::SkippedRuns;
use open ':std', ':encoding(UTF-8)';
use File::Basename;

sub html {
  my ($species, $our_studies, $other_studies) = @_;
  
  my $studies_tmpl = HTML::Template->new(path => [dirname (__FILE__) ], filename => "StudiesPage/templates/studies.tmpl");

  (my $species_param = $species) =~ s/_/ /g;
  $species_param = ucfirst($species_param);
  $studies_tmpl->param(SPECIES => $species_param);
  $studies_tmpl->param(MORE_WELCOMING_MESSAGE_FOR_SKIPPED_STUDIES => ($species_param =~ /Caenorhabditis elegans|Pristionchus pacificus/));

  my @toc = map {{
    TOC_ITEM_ID => $_->{study_id},
    TOC_ITEM_NAME => $_->{config}{title},
  }} @{$our_studies};
  push @toc, {
    TOC_ITEM_ID => "wbps_expression_other",
    TOC_ITEM_NAME => @{$our_studies} == 1 ? "1 other study" : sprintf("%s other studies", scalar @{$our_studies}),
  } if @{$our_studies};
  $studies_tmpl->param(TOC => \@toc);

  $studies_tmpl->param(RUNSTUDIES => join("\n", map {WbpsExpression::StudiesPage::Study->new($_)->to_html } @{$our_studies} ));
  $studies_tmpl->param(OTHERSTUDIES => join ("\n", map {WbpsExpression::StudiesPage::SkippedRuns->new($_)->to_html} @{$other_studies}));

  return $studies_tmpl->output;
}

1;
