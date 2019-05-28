use strict;
use warnings;

package WbpsExpression::IncomingStudies::RnaseqerResults;
use Log::Any '$log';
use LWP;
use JSON;
use File::Basename;
use DateTime;
use DateTime::Format::Strptime;
use List::Util qw/pairs unpairs/;
use File::Slurp qw(read_file read_dir);
# use Smart::Comments '###';

# 
# { study_id -> {
#     assembly_used,
#     source_dirs_by_run :: {location -> URL, end -> "pe|se", quality -> Int},
#     replicates_by_run :: {Run -> Replicate}
#   }
# }

sub get_json {
  my ($url) = @_;
  $log->info("RnqseqerResults::get_json LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:" . $response->status_line . "\n"
    unless $response->is_success;
  return $response->decoded_content && decode_json($response->decoded_content);
}

my $CAN_SEE_EBI_FILESYSTEM = -d "/nfs/ftp";

sub get_ftp_dir {
  my ($url) = @_;
 if ($CAN_SEE_EBI_FILESYSTEM and $url =~ m{ftp://ftp.ebi.ac.uk}){
    # Use a shortcut
    # Not exactly the same: EBI's ftp server replies with ls -l output
	(my $local = $url) =~ s{ftp://ftp.ebi.ac.uk}{/nfs/ftp};
	$log->info("RnaseqerResults::get_ftp_dir read_dir $local");
	return join "\n", read_dir $local;
  } else {
	$log->info("RnaseqerResults::get_ftp_dir LWP::get $url");
    my $response = LWP::UserAgent->new->get($url);
	die "$url error:" . $response->status_line . "\n"
	  unless $response->is_success;
	return $response->decoded_content;
  }
}

sub get_end_for_run {
  my ($ftp_path, $run_id) = @_;
  my $listing = get_ftp_dir($ftp_path);
  my @files = map {basename $_} $listing =~ /($run_id\..*)/g;
  my $pe = grep { /^$run_id.pe/ } @files;
  my $se = grep { /^$run_id.se/ } @files;
  my $end = ($pe xor $se ) 
    ? $pe ? "pe" : "se"
    : die $listing;
#### @files
#### $end
  die "No counts: $ftp_path" unless grep {$_ eq "$run_id.$end.genes.raw.htseq2.tsv"} @files;  
  die "No TPMs: $ftp_path" unless grep {$_ eq "$run_id.$end.genes.tpm.htseq2.irap.tsv"} @files;  
  die "No bigwig: $ftp_path" unless grep {$_ eq "$run_id.nospliced.bw"} @files;
  return $end;
}

sub get_results_by_study {
  my ($species) = @_;
  my %runs_per_study;
  for my $o (@{get_json ("https://www.ebi.ac.uk/fg/rnaseq/api/json/20/getRunsByOrganism/$species") || []}){
    push @{$runs_per_study{$o->{STUDY_ID}}}, $o if $o->{STATUS} eq "Complete";
  };
  my %result;
  for my $study_id (keys %runs_per_study){
    my @runs = @{$runs_per_study{$study_id}};
    my $assembly_used = value_for_newest_key(map {
       $_->{LAST_PROCESSED_DATE} => $_->{ASSEMBLY_USED}
    } @runs);
    my %source_dirs_by_run;
    for my $run (@runs){
#### $run
      my $run_id = $run->{RUN_IDS};
      my $ftp_path = dirname ($run->{BIGWIG_LOCATION});
      $source_dirs_by_run{$run_id} = $run->{ASSEMBLY_USED} eq $assembly_used ? {
        location => $ftp_path,
        end => get_end_for_run($ftp_path, $run_id),
        quality => $run->{MAPPING_QUALITY},
      } : {};
    } 
    my %runs_by_sample;
    for (@runs){
       push @{$runs_by_sample{$_->{SAMPLE_IDS}}} , $_->{RUN_IDS};
    }
    my %replicates_by_run;
    for my $sample_id (keys %runs_by_sample){
       my @run_ids_for_sample_id = @{$runs_by_sample{$sample_id}};
       my $sample_has_more_than_half_ids = scalar @run_ids_for_sample_id * 2 > scalar @runs;
       for my $run_id (@run_ids_for_sample_id){
          $replicates_by_run{$run_id} = ($sample_has_more_than_half_ids || 1 == @run_ids_for_sample_id ? $run_id : $sample_id);
       }
    }
    $result{$study_id} = {
       assembly_used => $assembly_used,
       source_dirs_by_run => \%source_dirs_by_run,
       replicates_by_run => \%replicates_by_run,
    };
  }
  return \%result;
}


sub value_for_newest_key {
  my %h = @_;
  my $f = DateTime::Format::Strptime->new(pattern=> "%a %b %e %Y %T", strict=>1); # Fri Jun 19 2015 18:20:10
  my ($k, $v, @xs) = unpairs sort {DateTime->compare($f->parse_datetime($a->[0]), $f->parse_datetime($b->[0]))} pairs %h;
  return $v;
}
1;
