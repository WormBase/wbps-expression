use strict;
use warnings;

package WbpsExpression::IncomingStudies::RnaseqerFtp;
use Log::Any '$log';
use LWP;
use File::Basename;
use List::Util qw/pairs unpairs/;
use File::Slurp qw(read_file read_dir);
# use Smart::Comments '###';

my $CAN_SEE_EBI_FILESYSTEM = -d "/nfs/ftp";

sub get_ftp_dir {
  my ($url) = @_;
 if ($CAN_SEE_EBI_FILESYSTEM and $url =~ m{ftp://ftp.ebi.ac.uk}){
    # Use a shortcut
    # Not exactly the same: EBI's ftp server replies with ls -l output
	(my $local = $url) =~ s{ftp://ftp.ebi.ac.uk}{/nfs/ftp};
	$log->info("RnaseqerFtp::get_ftp_dir read_dir $local");
	return join "\n", read_dir $local;
  } else {
	$log->info("RnaseqerFtp::get_ftp_dir LWP::get $url");
    my $response = LWP::UserAgent->new->get($url);
	die "$url error:" . $response->status_line . "\n"
	  unless $response->is_success;
	return $response->decoded_content;
  }
}

sub get_end_for_run {
  my ($run_id, $ftp_path) = @_;
  my $listing = get_ftp_dir($ftp_path);
  my @files = map {basename $_} $listing =~ /($run_id\..*)/g;
  my $pe = grep { /^$run_id.pe/ } @files;
  my $se = grep { /^$run_id.se/ } @files;
  my $end = ($pe xor $se ) 
    ? $pe ? "pe" : "se"
    : die $listing;
#### @files
#### $end
  die "No counts: $ftp_path" unless grep {$_ eq "$run_id.$end.genes.raw.htseq2.tsv" } @files;  
  die "No TPMs: $ftp_path" unless grep {$_ eq "$run_id.$end.genes.tpm.htseq2.irap.tsv"} @files;  
  die "No bigwig: $ftp_path" unless grep {$_ eq "$run_id.nospliced.bw"} @files;
  return $end;
}
1;
