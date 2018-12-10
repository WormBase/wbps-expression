use PublicResources::Links;
use Test::More;
my $links = 'PublicResources::Links';
my $result = $links->misc_links("run_id", "study_id", "<data_location_url>");
is_deeply($result,
   {
      "ENA study" => '<a href="https://www.ebi.ac.uk/ena/data/view/run_id">Study page: run_id</a>',
      "Mapping results" => '<a href="<data_location_url>">RNASeq-er processing directory: study_id</a>',
  },
   "unit tests are good"
) or diag explain $result;
done_testing();
