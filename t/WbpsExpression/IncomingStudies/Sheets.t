use Test::More;

use WbpsExpression::IncomingStudies::Sheets;
is_deeply(WbpsExpression::IncomingStudies::Sheets::_aoa_to_double_hash([]), {}, "Null case");
is_deeply(WbpsExpression::IncomingStudies::Sheets::_aoa_to_double_hash([
  [  "study", "run", "desc"],
]), {"study" => {"run"=>"desc"}});
is_deeply(WbpsExpression::IncomingStudies::Sheets::_aoa_to_double_hash([
  [  "study", "run", "desc", "long_description"],
]), {"study" => {"run"=>["desc", "long_description"]}});

is_deeply(WbpsExpression::IncomingStudies::Sheets::_aoa_to_double_hash([
[  "study", "run", "desc", "long_description"],
[  "study", "run_2", "desc_2"],
]), {"study" => {"run"=>["desc", "long_description"], "run_2"=>"desc_2" }});
done_testing;
