use strict;
use warnings;

# Copied as a single file from https://github.com/wormbase/wbps-expression into WBPS web code
# Please no dependencies apart from those listed in ensembl-webcode cpanfile:
# https://github.com/Ensembl/ensembl-webcode/blob/master/cpanfile

package EnsEMBL::Web::Component::Gene::WbpsExpression;
use File::Basename;

sub new {
  my ($class, $args) = @_;
  return bless $args, $class;
}
sub from_folder {
  my ($species, $studies_path) = @_;
  my $studies_file = "$studies_path/$species.studies.tsv";
  my %metadata;
  open (my $fh, "<", $studies_file) or return;
  while(<$fh>){
    chomp;
    my ($accession, $category, $title) = split "\t";
    $metadata{$accession}{title} = $title;
    $metadata{$accession}{category} = $category;
  }
  my @studies;
  for my $study_path (grep {-d $_}  glob("$studies_path/*")){
    my $study_id = basename $study_path;
    my $study = {
      study_id => $study_id,
      study_title => $metadata{$study_id}{title},
      study_category => $metadata{$study_id}{category},
    };
    if ($metadata{$study_id}{category} eq "Response to treatment"){
      my @files = glob("$study_path/$study_id.de.*.tsv");
      next unless @files;
      for my $file_path (@files){
        my ($type) = m{$study_path/$study_id\.de\.(.*)\.tsv};
        $study->{de}{$type} = $file_path;
      }
    } elsif ($metadata{$study_id}{category} eq "Other" ) {
      my $tpms_per_run = "$study_path/$study_id.tpm_per_run.tsv";
      next unless -s $tpms_per_run;
      $study->{tpms_per_run} = $tpms_per_run;
    } else {
      my $tpms_per_condition = "$study_path/$study_id.tpm.tsv";
      next unless -s $tpms_per_condition;
      $study->{tpms_per_condition} = $tpms_per_condition;
    }
    push @studies, $study;
  }
  return &new(__PACKAGE__, {studies => \@studies});
}


sub render_page {
  my ($self, $gene_id, $category) = @_;
  my $species = "TODO - from request, or from \$self?";
  my @studies = grep {$_->{study_category} eq $category} @{$self->{studies}};
  return html_no_studies_for_category($gene_id, $category) unless @studies;
  my $header = html_header($gene_id, $category);
  return (
    "<html><body>"
    . html_header($gene_id, $category)
    . "<br>"
    . join("<br>", response_as_html_panes($gene_id, $category, \@studies)) 
    . "</body></html>"
  );
}
sub response_as_html_panes {
  my ($gene_id, $category, $studies) = @_;
  if ($category eq  "Response to treatment"){
    my ($fold_changes, $studies_with_no_results) = list_of_fold_changes_in_studies_and_studies_with_no_results($gene_id, $studies);
    return (
       (@{$fold_changes} ? html_fold_changes_table($fold_changes) : ()),
       (@{$studies_with_no_results} ? html_studies_with_no_results($studies_with_no_results) : ()),
    );
  } elsif ($category eq "Other") {
    return map {html_flat_horizontal_table($_)} summary_stats_in_tables($gene_id, $studies);
  } else {
    return map {html_flat_horizontal_table($_)} tpms_in_tables($gene_id, $studies);
  }
}
sub html_no_studies_for_category {
  my ($gene_id, $category) = @_;
  return "<p> No results for gene $gene_id in category $category </p>";
}
sub html_header {
  my ($gene_id, $category) = @_;
  return "<h2> Gene expression - $category </h2>"; 
}
sub html_fold_changes_table {
  my ($fold_changes) = @_;
  return (
     "<table>"
     . "<th>"
       . "<td>Study</td>"
       . "<td>Contrast type</td>"
       . "<td>Contrast name</td>"
       . "<td>Log<sub>2</sub>-fold change</td>"
     . "</th>"
     . join ("\n", map {
       "<tr>"
         . "<td>" . html_study_link($_) . "</td>"
         . "<td>" . ucfirst($_->{contrast_type}) . "</td>"
         . "<td>" . $_->{contrast_name} . "</td>"
         . "<td>" . $_->{fold_change} . "</td>"
       . "</tr>"
     } @{$fold_changes})
     . "</table>"
  );
}
sub list_of_fold_changes_in_studies_and_studies_with_no_results {
  my ($species, $gene_id, $studies) = @_;
  my @fold_changes;
  my @studies_with_no_results;

  for my $study (@{$studies}){
     my @fold_changes_for_study;
     while (my ($type, $path) = each %{$study->{de}}) {
        my ($contrasts, $fold_changes) = search_in_file($path, $gene_id);
        for my $i (0..$#$contrasts){
           push @fold_changes_for_study, {
              species => $species,
              study_id => $study->{study_id},
              study_title => $study->{study_title},
              contrast_type => $type,
              contrast_name => $contrasts->[$i],
              fold_change => $fold_changes->[$i],
           }
        }
     }
     push @fold_changes, $_ for @fold_changes_for_study;
     push @studies_with_no_results, $study unless @fold_changes_for_study;
  }
  return \@fold_changes, \@studies_with_no_results;
}

# Adapted from: https://metacpan.org/source/SHLOMIF/Statistics-Descriptive-3.0612/lib/Statistics/Descriptive.pm#L620
sub quantile {
    my ( $data_sorted, $QuantileNumber ) = @_;
 
    return $data_sorted->[0] if ( $QuantileNumber == 0 );
    my $count = @{$data_sorted};
    return $data_sorted->[ $count - 1 ] if ( $QuantileNumber == 4 );
 
    my $K_quantile = ( ( $QuantileNumber / 4 ) * ( $count - 1 ) + 1 );
    my $F_quantile = $K_quantile - POSIX::floor($K_quantile);
    $K_quantile = POSIX::floor($K_quantile);
 
    my $aK_quantile     = $data_sorted->[ $K_quantile - 1 ];
    return $aK_quantile if ( $F_quantile == 0 );
    my $aKPlus_quantile = $data_sorted->[$K_quantile];
 
    my $quantile = $aK_quantile
      + ( $F_quantile * ( $aKPlus_quantile - $aK_quantile ) );
 
    return $quantile;
}
sub summary_stats_in_tables {
  my ($species, $gene_id, $studies) = @_;
  my @result;
  for my $study (@{$studies}){
     my ($runs, $expression_tpm) = search_in_file($study->{tpms_per_run}, $gene_id);
     my @expression_tpm_sorted = sort @{$expression_tpm //[]};
     push @result, {
       species => $species,
       study_id => $study->{study_id},
       study_title => $study->{study_title},
       column_headers => ["N", "min", "Q1", "Q2", "Q3", "max"],
       values => [ scalar @expression_tpm_sorted, map {quantile(\@expression_tpm_sorted, $_)} (0..4)],
     } if $runs and $expression_tpm;
  }
  return \@result;
}
sub tpms_in_tables {
  my ($species, $gene_id, $studies) = @_;
  my @result;  
  for my $study (@{$studies}){
     my ($conditions, $expression_tpm) = search_in_file($study->{tpms_per_condition}, $gene_id);
     push @result, {
       species => $species,
       study_id => $study->{study_id},
       study_title => $study->{study_title},
       column_headers => $conditions,
       values => $expression_tpm,
     } if $conditions and $expression_tpm;
  }
  return \@result;
}


sub search_in_file {
  my ($path, $gene_id) = @_;
  my $l = `grep --max-count=1 "^$gene_id" $path`;
  chomp $l;
  return unless $l;
  my ($id, @xs) = split "\t", $l;
  return unless $id eq $gene_id and @xs;
  my $h = `grep --max-count=1 "^\t" $path`;
  chomp $h;
  return unless $h;
  my ($blank, @hs) = split "\t", $h;
  return unless not($blank) and @hs;
  return unless @hs == @xs;
  return \@hs, \@xs;
}
1;
