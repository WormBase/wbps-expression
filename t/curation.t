use strict;
use warnings;
use Test::More;
use File::Basename qw/basename dirname/;
use File::Find;
use File::Slurp qw/read_file/;
use FindBin;
use WbpsExpression::Study;
use List::Util qw/pairs/;

sub test_folder {
  my ($path) = @_;
  subtest $path => sub {
    my $study_id = basename $path;
	ok( -f "$path/$study_id.design.tsv" or -f "$path/$study_id.skipped_runs.tsv",  "$path design/runs exist" );
	ok( -f "$path/$study_id.config.yaml", "$path config exists" );
	ok( -f "$path/$study_id.sources.tsv", "$path sources TSV exists" );
	my $study;
	eval {
	  $study = WbpsExpression::Study->from_folder($path);
	};
	ok ($study, "Construct study from folder") or return;
    my @checks = $study->all_checks;
	  subtest "Study checks" => sub {
	  my @pairs = pairs @checks;
	  plan tests => scalar @pairs;
	  for my $pair (@pairs) {
		ok( $pair->[1], $pair->[0] );
	  }
	}
  }
}
my @studies;
my $arg = $ARGV[0] // "";
my $study_folder_pattern = $arg =~ /^[A-Z]+\d+$/ ? qr/$arg$/ : qr/$arg.*[A-Z]+\d+$/;
find(
  sub {
    my $bn = basename $File::Find::name;
    my $bbn = basename (dirname $File::Find::name);
    push @studies, $File::Find::name if -d $File::Find::name and $File::Find::name =~ $study_folder_pattern;
  },
  "$FindBin::Bin/../studies"
);

test_folder($_) for @studies;

done_testing;
