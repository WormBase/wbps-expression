use Test::More;

use Curation::Sheets;
is_deeply(Curation::Sheets::_aoa_to_double_hash([]), {}, "Null case");
is_deeply(Curation::Sheets::_aoa_to_double_hash([
  [  "study", "run", "desc"],
]), {"study" => {"run"=>"desc"}});
is_deeply(Curation::Sheets::_aoa_to_double_hash([
  [  "study", "run", "desc", "long_description"],
]), {"study" => {"run"=>["desc", "long_description"]}});

is_deeply(Curation::Sheets::_aoa_to_double_hash([
[  "study", "run", "desc", "long_description"],
[  "study", "run_2", "desc_2"],
]), {"study" => {"run"=>["desc", "long_description"], "run_2"=>"desc_2" }});
done_testing;
