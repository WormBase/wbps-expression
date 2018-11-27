use PublicResources::Descriptions;
use Test::More;
use Test::MockModule;
my $species = "schistosoma_mansoni";
my $study_id = "SRP1234000";
my $run_id = "SRR1230321";
my $default_short_desc = "";
my $default_long_desc = "sample from Schistosoma mansoni";

sub assert_run_description {
  my ($attributes, $curated, $expected_description_short, $expected_description_full, $desc) = @_;
   my ($description_short, $description_full) =
     PublicResources::Descriptions->new($species, {$study_id=>$curated})->run_description($study_id, $run_id, $attributes);
   is_deeply($description_short, $expected_description_short, "description_short $desc");
   is_deeply($description_full, $expected_description_full, "description_full $desc");
}

assert_run_description({},{}, "", $default_long_desc, "Default long description ok in null case");
assert_run_description({},{$run_id => "curated value"}, "curated value", "curated value", "curated value");
assert_run_description({},{$run_id => ["curated value", "curated value long"]}, "curated value", "curated value long", "curated value long");
assert_run_description({organism_part=>"head"},{}, "head", "head","Use attributes"); 
assert_run_description({organism_part=>"head", replicate =>"1"},{}, "head", "head","Use attributes - skip blacklist"); 
assert_run_description({organism_part=>"head", strain => "NMRI"},{}, "head, NMRI", "head, strain NMRI","Use multiple attributes"); 
sub sample_name_ok {
  my ($sample_name, $desc) = @_;
  assert_run_description({sample_name => $sample_name}, {}, "$sample_name", "$sample_name", 
    $desc // "sample_name_ok $sample_name");
}
sub sample_name_rejected {
  my ($sample_name, $desc) = @_;
  assert_run_description({sample_name => $sample_name}, {}, $default_short_desc, $default_long_desc,
    $desc // "sample_name_rejected $sample_name");
}
sample_name_ok("words words words");
sample_name_rejected("");
sample_name_rejected("123");
sample_name_rejected("short SN");
sample_name_rejected("onlyOneWordIsBadSampleName");
sample_name_rejected("something about being private in GEO");
sample_name_rejected("sample from S. mansoni");

sub assert_study_description {
  my ($study_attributes, $expected, $desc) = @_;
   $desc //= "assert_study_description $study_id -> $expected";
   is_deeply(
     [PublicResources::Descriptions->new($species, {})->study_description($study_id, $study_attributes)],
     $expected, $desc);
}

assert_study_description(
    {},
    ["Schistosoma mansoni study","Schistosoma mansoni study"],
    "Default value",
);
assert_study_description(
    {study_description => "Desc"},
    ["Schistosoma mansoni study","Desc"],
);
assert_study_description(
    {"study_title" => "Study title"},
    [("Study title") x 2 ],
);
assert_study_description(
    {"study_title"=> "Study title", "study_description" => "Desc",},
    ["Study title","Desc"],
);
sub assert_names_get_better {
  my ($from, $to) = @_;
  assert_study_description(
    {"study_title" => $from },
    [($to) x 2],
    "Names get better: $from -> $to",
  );
}

assert_names_get_better("Globodera_pallida_transcriptomics", "Globodera pallida transcriptomics");
assert_names_get_better("MIYAZAKI", "Miyazaki");



done_testing();
