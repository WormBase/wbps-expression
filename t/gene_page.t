use strict;
use warnings;
use EnsEMBL::Web::Component::Gene::WbpsExpression;
use Test::More;
use File::Temp qw/tempdir/;
use File::Slurp qw/write_file/;
use File::Path qw/make_path/;

my $dir = tempdir(CLEANUP => 1);
my $f1 = <<EOF;
#temp file
	heads	tails
g1	1.1	1.2
g2	2.1	2.2
EOF
my $study_id = "SRP071241";
my $category = "Organism parts";
my $study_title = "Comparison of gene expression between female Schistosoma mansoni heads and tails";
my $f1_name = "$study_id.tpm.tsv";
my $f1_path = join("/", $dir, $study_id, $f1_name);

my $species = "schistosoma_mansoni";
my $assembly = "Smansoni_v7";
make_path(join("/", $dir, $study_id));
write_file($f1_path, $f1);
write_file(join("/", $dir, "$species.$assembly.studies.tsv"), "$study_id\t$category\t$study_title\n");

my $subject = EnsEMBL::Web::Component::Gene::WbpsExpression::from_folder(
   $species, $assembly, $dir
);

is_deeply($subject, bless({
  studies => [{
     study_id => "$study_id",
     study_title => $study_title,
     study_category => $category,
     tpms_per_condition => $f1_path,
  }]
}, 'EnsEMBL::Web::Component::Gene::WbpsExpression'), "Create reads in the config");

is_deeply($subject->expression_for_gene_in_category("invalid ID", $category), [], "Null case - gene");
is_deeply($subject->expression_for_gene_in_category("g1", "Different category"), [], "Null case - category");

is_deeply($subject->expression_for_gene_in_category("g1", $category), [{
  study_id => $study_id,
  study_title => $study_title,
  conditions => ["heads", "tails"],
  expression_tpm => [1.1, 1.2],
}] , "One line");
# 
my $spec = << '/SPEC';
studies have a category - known from the tsv

We will show a section at a time, among:
- Life stages
- Organism parts
- Variation within species
- Cell types
- Response to treatment
- Other
The first four sections will have the same, 'baseline', display format.

### Left hand side menu content
The content should be determined during 'packing' the species, as an expandable menu, with the parent node saying 'Expression' in gray, and as its children the subset of above six categories for which there is data, in blue, as links to individual section.
As a common special case, 'Expression -> Other' should appear as one-level non-expandable 'Expression'.

### Interface description
Each page will have a section title (e.g. "Gene expression - life stages") and one or more panes.
In the 'baseline' format, each pane will have a title:
Expression (TPM) - [link to study page](title)
Under a title, there should be a table with the data. Usually, this should be a "flat", one row table listing the values in TPMs, and column headers listing the conditions.
Where possible to construct, in the case of multi-dimensional designs, this should be a two-dimensional table with both row and column headers, with column headers listing the largest dimension of the design, and row headers listing the remaining dimensions.

The 'Other' category should be like above, except the data in a one row table should be statistical aggregates, and column headers listing the statistics : n, low, Q1, Q2, Q3, high.

The 'Response to treatment' category should present results from each study where significant expression was found, with a title:
Differential expression (log2foldchange, where adjusted p-value < 0.05) - [link to study page](title)
The values should be presented in a one column vertical table, with row headers being the comparison type + contrast name.

If there was no data for some studies, then under the results a separate pane should appear, titled "No significant results" and a bullet point style list of studies, of the format [link to study page](title).

### Presentation layer
The display code will create styled HTML from a plain Perl payload of the format: 
```
{
  section_title, # saying "Gene expression" and category
  panes : [{
     pane_title,
     study_title,
     study_id,
     type = "table_flat_horizontal",
     column_headers : [String],
     values : [Number],
   } | {
     pane_title,
     study_title,
     study_id,
     type = "table_rows_and_columns",
     column_headers : [String],
     row_headers : [String],
     values: [[Number]],
  } | {
     pane_title,
     study_title,
     study_id,
     type = "table_flat_vertical",
     row_headers : [String],
     values : [Number],
  } | {
     pane_title,
     type = "list_of_studies",
     values : [{
       study_title,
       study_id,
     }]
  }]
}
Studies should be linked to sections in the study explorer page, as URLs following pattern `/expression/$species/#$study_id`.
```
### Architecture
We are confident that this functionality can be achieved without storing the data in a relational database, or introducing library depencencies to the web server code. The web server code will gain one Perl module tasked with retrieving the data, with a class capturing the desired behaviour.

An instance of that class should be created for each species. The construction parameters should be species name, assembly, and a data folder, location of which will be provided via config. Upon construction, the module will parse a listing of study id | category | title, and note the locations of all data files it will need access to. The instance should be cached in web server memory.

Then on request, the module should get the data it needs from each file using a `grep` on the file twice: once for the data values, once for the conditions from the data header. If the approach of table scanning the whole file will not be fast enough, we will try storing the data in binary files as B-trees, leveraging a BerkeleyDB type component already installed on the web servers. We anticipate that this may not be necessary - the files are of low quantity and modest size.

/SPEC

done_testing;
