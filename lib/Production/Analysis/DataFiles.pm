use strict;
use warnings;
use File::Slurp qw/read_dir/;
use Model::DataDir;
package Production::Analysis::DataFiles;

sub read_file_into_hash {
  my ($path) = @_;
  my %result;
  open(my $fh, "<", $path) or die "$path: $!";
  my $header = <$fh>;
  while(<$fh>){
    my ($k, $v) = split "\t";
    $result{$k} = $v;
  }
  close $fh;
  return \%result;
}

sub read_files_into_averaged_hash {
  my (@paths) = @_;
  my %result;
  for my $path (@paths) {
	open(my $fh, "<", $path) or die "$path: $!";
	my $header = <$fh>;
	while(<$fh>){
	  my ($k, $v) = split "\t";
	  push @{$result{$k}},$v;
	}
    close $fh;
  }
  for my $k (keys %result){
    $result{$k} = sprintf("%.1f", calculate_median($result{$k})); 
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
sub write_named_hashes {
  my ($name_to_data_pairs, $out_path) = @_;
  my %row_labels;
  for my $p (@{$name_to_data_pairs}){
    for my $label(keys %{$p->[1]}){
      $row_labels{$label}++;
    }
  }
  open(my $fh, ">", $out_path) or die "$out_path: $!";
  print $fh join ("\t", "", map {$_->[0]} @{$name_to_data_pairs})."\n";
  for my $row (sort keys %row_labels){
     print $fh join ("\t",$row, map {$_->[1]->{$row} // ""} @{$name_to_data_pairs})."\n";
  }
  close $fh;
}
sub aggregate {
  my ($name_to_path_pairs, $out_path) = @_;
  my @name_to_data_pairs = map {[$_->[0], read_file_into_hash($_->[1])]}  @{$name_to_path_pairs};
  write_named_hashes(\@name_to_data_pairs, $out_path);
}
sub average_and_aggregate {
  my ($name_to_pathlist_pairs, $out_path) = @_;
  my @pairs =  map {[$_->[0], read_files_into_averaged_hash(@{$_->[1]})]} @{$name_to_pathlist_pairs};
  write_named_hashes(\@name_to_data_pairs, $out_path);
}
1;
