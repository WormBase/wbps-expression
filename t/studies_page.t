use WbpsExpression::StudiesPage;
use Test::More;
 
like(WbpsExpression::StudiesPage::html("brugia_malayi", [],[]), qr/Brugia malayi.*studies/);
done_testing;
