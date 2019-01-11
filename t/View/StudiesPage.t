use View::StudiesPage;
use Test::More;
use FindBin;
 
$ENV{HTML_TEMPLATE_ROOT} = "$FindBin::Bin/../../lib/View/template";
my @studies = ();
my $html = View::StudiesPage->new("brugia_malayi", @studies)->to_html;
like($html, qr/Brugia malayi.*studies/);
done_testing;
