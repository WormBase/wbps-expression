use View::StudiesPage;
use Test::More;

my @studies = ();
my $html = View::StudiesPage->new("brugia_malayi", @studies)->to_html;
like($html, qr/Brugia malayi.*studies/);
done_testing;
