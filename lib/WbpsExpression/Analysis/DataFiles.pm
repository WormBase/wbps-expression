use strict;
use warnings;
package WbpsExpression::Analysis::DataFiles;
use File::Slurp qw/read_dir/;
use WbpsExpression::Analysis::Common;
use LWP;
use Log::Any '$log';
use Scalar::Util qw/looks_like_number/;
# use Smart::Comments '###';
my $CAN_SEE_EBI_FILESYSTEM = -d "/nfs/ftp";

sub open_read_fh {
  my ($path) = @_;
  $log->info("open_read_fh $path");
  my $fh;
  if (not ref $path and $path =~ m{ftp://ftp.ebi.ac.uk}) {
    if( $CAN_SEE_EBI_FILESYSTEM ) {
       $path =~ s{ftp://ftp.ebi.ac.uk}{/nfs/ftp};
       open ($fh, "<", $path) or die "$path: $!";
    } else {
       my $response = LWP::UserAgent->new->get($path);
       die "$path error:".$response->status_line."\n" unless $response->is_success;
       my $body = $response->decoded_content;
       open ($fh, "<", \$body) or die "$path: $!";
    }
  } else {
    open ($fh, "<", $path) or die "$path: $!";
  }
  return $fh;
}

sub read_file_into_hash {
  my ($path) = @_;
#### read_file_into_hash: $path
  my %result;
  my $fh = open_read_fh($path);
  my $header = <$fh>;
  chomp $header;
  my ($k, $v) = split "\t", $header;
  $result{$k} = $v if $k and $k ne 'Gene' and looks_like_number $v;
  while(<$fh>){
    chomp;
    my ($k, $v) = split "\t";
    $result{$k} = $v;
  }
  close $fh;
  return \%result;
}

sub read_files_into_averaged_hash {
  my (@paths) = @_;
  my %result;
  for my $i (0 .. $#paths){
    for my $path (@{$paths[$i]}) {
	my $fh = open_read_fh($path);
	my $header = <$fh>; 
#### $header
	while(<$fh>){
      chomp;
#### $_
	  my ($k, $v) = split "\t";
	  push @{$result{$k}[$i]},$v;
	}
      close $fh;
    }
  }
#### %result
  for my $k (keys %result){
    my @vs = map { calculate_median( $_)} @{$result{$k}};
    $result{$k} = sprintf("%.1f", calculate_median(\@vs)); 
  }
  return \%result;
}

#Adapted from: https://metacpan.org/source/SHLOMIF/Statistics-Descriptive-3.0612/lib/Statistics/Descriptive.pm#L237
sub calculate_median {
    my ( $expressions ) = @_;
    my @expressionsSorted = sort {$a <=> $b} @$expressions;
    my $count = @expressionsSorted;
    ##Even or odd
    if ($count % 2){
        return @expressionsSorted[($count-1)/2];
    } else {
        return (
            (@expressionsSorted[($count)/2] + @expressionsSorted[($count-2)/2] ) / 2
        );
    }
}
sub aggregate {
  my ($name_to_path_pairs, $out_path, @frontmatter) = @_;
  my @name_to_data_pairs = map {
   my $name = $_->[0];
   my $data = read_file_into_hash($_->[1]);
   [$name, $data]
  }  @{$name_to_path_pairs};
  WbpsExpression::Analysis::Common::write_named_hashes(\@name_to_data_pairs, $out_path, @frontmatter);
}
sub average_and_aggregate {
  my ($name_to_pathlist_pairs, $out_path, @frontmatter) = @_;
  my @name_to_data_pairs =  map {
    my $name = $_->[0];
    my $data = read_files_into_averaged_hash(@{$_->[1]});
    [$name, $data]
  } @{$name_to_pathlist_pairs};
  WbpsExpression::Analysis::Common::write_named_hashes(\@name_to_data_pairs, $out_path, @frontmatter);
}
1;
