use strict;
use warnings;
package PublicResources::Resources::LocallyCachedResource;
use LWP;
use YAML;
use File::Path qw(make_path);
use JSON;
use XML::Simple;
use Text::CSV qw(csv);
use File::Slurp qw(read_file read_dir);
use Log::Any '$log', default_adapter => 'Stderr';

my $CAN_SEE_EBI_FILESYSTEM = -d "/nfs/ftp";

#use Smart::Comments;
sub new {
    my ($class,$root_dir,$species, @other_args) = @_;
### require: $species
    my $dir = "$root_dir/$species";
    my $path_to_local_copy = "$dir/$class";
    make_path $dir;

    YAML::DumpFile($path_to_local_copy, 
       $class->_fetch($species, @other_args)
    ) unless -s $path_to_local_copy;

    return bless YAML::LoadFile($path_to_local_copy), $class;
}

sub get_text {
  my ($class, $url) = @_;
  my $errors;
  if ($CAN_SEE_EBI_FILESYSTEM and $url =~ m{ftp://ftp.ebi.ac.uk}){
     (my $local = $url) =~ s{ftp://ftp.ebi.ac.uk}{/nfs/ftp};
     if( -d $local){
        #Not exactly the same: EBI's ftp server replies with ls -l output
        $log->info(__PACKAGE__." get_text: read_dir $local");
        return join "\n", read_dir $local;
     }elsif(-f $local){
        $log->info(__PACKAGE__." get_text: read_file $local");
        return read_file $local;
     }
  }
  $log->info(__PACKAGE__." get_text LWP::get $url");
  my $response = LWP::UserAgent->new->get($url);
  if($response->is_success){
    return $response->decoded_content;
  } else {
    die "$url error:".$response->status_line."\n";
  }
}
sub get_csv {
  my $class = shift;
  my $text = $class->get_text(@_);
  return csv ({ allow_whitespace => 1, in=>\$text});
}
sub get_json { 
  my $class = shift;
  return decode_json($class->get_text(@_));
}

sub get_xml { 
  my $class= shift;
  return XMLin($class->get_text(@_));
}
1;
