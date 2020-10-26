use FindBin;
use lib "$FindBin::Bin/../local/lib/perl5";
#$ENV{PATH} = "$FindBin::Bin/../local/R-3.6.1/bin:$ENV{PATH}";
$ENV{R_LIBS_USER} = "$FindBin::Bin/../local/R-3.6.1/library"
use Statistics::R;

my $R = Statistics::R->new();
my $input_value = 1;
$R->set('x', $input_value);
$R->run(q`y <- x^2`);
my $output_value = $R->get('y');
print "y = $output_value\n";


