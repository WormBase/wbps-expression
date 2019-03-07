# Produces a folder of TSV files
# Format:
# study run location  num_reads_approx  mapping_quality_approx  @factors
use strict;
use warnings;
package PublicResources::Report;
use PublicResources::Rnaseq;
use File::Path qw(make_path);
use File::Slurp qw(write_file);
use File::Basename qw(dirname);
use Log::Any '$log', default_adapter => 'Stderr';

sub new {
    my ( $class, %args ) = @_;
    $args{root_dir} //=
      "$ENV{PARASITE_SCRATCH}/jbrowse/WBPS$ENV{PARASITE_VERSION}";
    return bless {
        dir       => "$args{root_dir}/tsv",
        resources => PublicResources::Rnaseq->new("$args{root_dir}/Resources"),
    }, $class;
}

sub path {
    my ( $self, @args ) = @_;
    my $result = join "/", $self->{dir}, @args;
    $log->info(__PACKAGE__." path $result");
    return $result;
}

sub make_report {
    my ($self, $species, $assembly, %opts) = @_;

    my @studies = sort {$a->{study_id} cmp $b->{study_id}} $self->{resources}->get( $species, $assembly );
    my $runs_path = $self->path("$species.runs.tsv");
    unlink $runs_path if -f $runs_path;
    
    return unless @studies;
    make_path dirname $runs_path;
    open(my $runs_fh, '>:utf8', $runs_path ) or die $runs_path;

    print $runs_fh join ("\t",
       "Study", "Track",
       "Library size (reads)", "Mapping quality (reads uniquely mapped)", 
    ) . "\n";
    for my $study (@studies) {
        my $study_id = $study->{study_id};
        for my $run ( sort {$a->{run_id} cmp $b->{run_id}} @{ $study->{runs} } ) {
            my $run_id = $run->{run_id};
            print $runs_fh join ("\t",
                 $study_id, join(": ", grep {$_} $run_id, $run->{run_description_short}),
                 $run->{attributes}{library_size_reads_approximate}, $run->{attributes}{fraction_of_reads_mapping_uniquely_approximate},
            ). "\n";     
        }
    }
    close $runs_fh;
    
    my $studies_path = $self->path("$species.studies.tsv");
    unlink $studies_path if -f $studies_path;

    open(my $studies_fh, '>:utf8', $studies_path) or die $studies_path;
    print $studies_fh join ("\t",
      "Study", "Submitting centre", "ENA first public", "ENA last update",
      "Description",
      "PubMed id", "PubMed",
    )."\n";
    for my $study (@studies){
      print $studies_fh join ("\t",
        $study->{study_id},
        $study->{attributes}{submitting_centre}, $study->{attributes}{"ENA first public"}, $study->{attributes}{"ENA last update"},
        $study->{study_description_short},
        join(", " , keys %{$study->{pubmed}|| {}}), join(", " , map {$_->[1]} values %{$study->{pubmed}|| {}}),
      )."\n";
    }
    close $studies_fh; 
}
1;
