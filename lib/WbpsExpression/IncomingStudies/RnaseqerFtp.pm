use strict;
use warnings;

package WbpsExpression::IncomingStudies::RnaseqerFtp;
use Log::Any '$log';
use LWP;
use File::Basename;
use List::Util qw/pairs unpairs/;
use File::Slurp qw(read_file read_dir);
# use Smart::Comments '###';
use Carp;
use Try::Tiny;

my $CAN_SEE_EBI_FILESYSTEM = -d "/nfs/ftp";

sub get_ftp_dir {
  my ($url) = @_;
 if ($CAN_SEE_EBI_FILESYSTEM and $url =~ m{ftp://ftp.ebi.ac.uk}){
    # Use a shortcut
    # Not exactly the same: EBI's ftp server replies with ls -l output
	(my $local = $url) =~ s{ftp://ftp.ebi.ac.uk/pub}{/nfs/ftp/public};
	$log->info("RnaseqerFtp::get_ftp_dir read_dir $local");
	my @local_dir;
	try {
     @local_dir = read_dir $local;
   } catch {
     my $msg = $_;
     confess "Can't read local directory $local: $msg"; 
   };
	return join "\n", @local_dir;
  } else {
	$log->info("RnaseqerFtp::get_ftp_dir LWP::get $url");
    my $response = LWP::UserAgent->new->get($url);
	die "$url error:" . $response->status_line . "\n"
	  unless $response->is_success;
	return $response->decoded_content;
  }
}

sub get_end_for_run {
my $this = (caller(0))[3];
  my ($run_id, $ftp_path, $rnaseqer_last_update) = @_;
  $run_id               || confess "$this called without positional arg 0 \$run_id";
  $ftp_path             || confess "$this called without positional arg 1 \$ftp_path";
  $rnaseqer_last_update || confess "$this called without positional arg 2 \$rnaseqer_last_update";
  my $listing = get_ftp_dir($ftp_path);
  my @files = map {basename $_} $listing =~ /($run_id\..*)/g;
  my $pe = grep { /^$run_id.pe/ } @files;
  my $se = grep { /^$run_id.se/ } @files;
  my $end = ($pe xor $se ) 
    ? $pe ? "pe" : "se"
    : die $listing;
#### @files
#### $end
  if ($rnaseqer_last_update ge "2019-04-15"){
    $log->error("ERROR: No counts: $ftp_path")  unless grep {$_ eq "$run_id.$end.genes.raw.featurecounts.tsv" } @files;
    $log->error("ERROR: No TPMs: $ftp_path")    unless grep {$_ eq "$run_id.$end.genes.tpm.featurecounts.irap.tsv"} @files;  
  } else {
    $log->error("ERROR: No counts: $ftp_path")  unless grep {$_ eq "$run_id.$end.genes.raw.htseq2.tsv" } @files;  
    $log->error("ERROR: No TPMs: $ftp_path")    unless grep {$_ eq "$run_id.$end.genes.tpm.htseq2.irap.tsv"} @files;  
  }
  $log->error("ERROR: No bigwig: $ftp_path")    unless grep {$_ eq "$run_id.nospliced.bw"} @files;
  return $end;
}
1;
