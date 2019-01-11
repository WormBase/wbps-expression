use strict;
use warnings;
package View::StudiesPage;
use View::Study;
use View::SkippedRuns;
use open ':std', ':encoding(UTF-8)';

sub new {
  my ($class, $species, $studies) = @_;
  return bless {species => $species, studies => $studies}, $class;
}

sub to_html {
  my ($self) = @_;
  my $studies_tmpl = HTML::Template->new(filename => 'studies.tmpl');

  $studies_tmpl->param(SPECIES => do {
    my $species = $self->{species};
    $species =~ s/_/ /g;
    ucfirst($species)
  });

  my $studies;

  for my $study (@{$self->{studies}{passing_checks}}){
    my ($skip) = grep {$_->{study_id} eq $study->{study_id}} @{$self->{studies}{skipped_runs}};
    if ($skip) {
      $study->{SKIPPED_RUNS} = $skip->{runs};
    }
     $studies .= View::Study->new($study)->to_html . "\n";
  }

  for my $study (@{$self->{studies}{failing_checks}}){
     $studies .= View::Study->new($study)->to_html . "\n";
  }

  $studies_tmpl->param(RUNSTUDIES => $studies);

  my $skip_studies;
  for my $skipped_r (@{$self->{studies}{skipped_runs}}){
    my $found = grep {$_->{study_id} eq $skipped_r->{study_id}} @{$self->{studies}{passing_checks}};
    if (!$found) {
      $skip_studies .= View::SkippedRuns->new($skipped_r)->to_html . "\n";
    }
  }

  $studies_tmpl->param(OTHERSTUDIES => $skip_studies);

  return $studies_tmpl->output;
}

1;
