package PublicResources::Links;

sub _link {
   my ($url, $name) = @_;
   return "<a href=\"$url\">$name</a>";
}

sub misc_links {
   my ($self, $study_id, $run_id, $data_location) = @_;
   
   my $result = {
      "Mapping results" =>
          _link($data_location, "RNASeq-er processing directory: $run_id"),
      "ENA study" =>
          _link("https://www.ebi.ac.uk/ena/data/view/$study_id", "Study page: $study_id"),
   };
   return $result;
}

1;
