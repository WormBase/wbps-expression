[![Build Status](https://travis-ci.org/WormBase/wbps-expression.svg?branch=master)](https://travis-ci.org/WormBase/wbps-expression)
# wbps-expression
Transcriptomic data for WormBase ParaSite

## Introduction
This is the pipeline for providing WormBase ParaSite with RNASeq data. It encompasses a curation platform, data retrieval and analysis, and a UI oriented around static pages and files.
### Curation
On a fresh run for a species of interest, the pipeline retrieves:
 - run metadata from [RNASeq-er](https://www.ebi.ac.uk/fg/rnaseq/api/), who retrieve it from [ENA](http://www.ebi.ac.uk/ena)
 - study and publication metadata from [ENA](http://www.ebi.ac.uk/ena), and [GEO](https://www.ncbi.nlm.nih.gov/geo/)
 - publication details from [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/)
 
and preserves the information as YAML, in a format corresponding to [pipeline's internal Perl modules](https://github.com/WormBase/wbps-expression/tree/master/lib/PublicResources/Resources) such that it can then be used on subsequent runs.

This is then used to update the [curation](https://github.com/WormBase/wbps-expression/tree/master/curation) folder, preserved together with source code in this repository. The pipeline makes some guesses on what to accept and reject - with a few exceptions, we only allow studies of at least six runs - and runs consistency checks on the annotation. A curator then amends the files and re-runs iteratively, until the checks pass and they are satisfied with the results.

### Data retrieval and analysis
#### Source files
We leave alignment, quantification, etc. to RNASeq-er. Their data, assembled by study but without any extra interpretation, is available for every study as "counts of aligned reads per run" and "TPMs per run".

#### TPMs by condition
Where there is enough replicates, we provide median per condition. Not many studies have both technical and biological replicates, but where they do, we take median of technical replicates for each biological replicates, and then take median of biological replicates.

#### Differential expression analysis
For studies where it makes sense to form contrasts from appropriate pairs of conditions, we run differential expression analysis. The pipeline picks contrasts automatically, through a number of  heuristics:
- compare everything to a reference condition, if there is a clear reference
- in the general case, try make contrasts from all pairs of conditions
- except if conditions differ by multiple types of characteristics, then pick only the pairs that differ by one characteristic
- except a parasite-specific curation of life stages as two characteristics (developmental stage and sex) should still work
- except drug treatment assays curated as treatment+concentration or treatment+timepoint should still work
It works slightly better than it sounds.

The analysis uses DESeq2 in a very standard way, with fold changes and p-values extracted and filtered past a significance threshold.

### Outputs
The analysis results are returned as files within a per-study directory structure, together with a listing of metadata for programmatic use else where, and a static HTML page presenting the content.

This HTML page is intended as a primary point of reference for people interested in the data. It lists the studies, with metadata for each, and links to the analysis results.

WormBase ParaSite deploys the content by syncing the data to a particular place in a file system within web servers' environment, with web servers configured to read content from there.

### Gene page
The repository also contains a module capable of searching in text files with `grep` and formatting results as HTML pages, which forms the basis of our gene page. Due to a small volume of the files we had no need for a database to store the data.

### Tracks
The pipeline also supports WormBase ParaSite track hubs and JBrowse displays - [code somewhere else](https://github.com/wormbase/wormbase-pipeline/), it uses the `PublicResources` modules and YAML files only.

## Curation tools
### Tracking progress
The curation files go together with the source code, and `git` is really good at tracking what happened, when, and why. `git status` will show you what new files appeared after a run. Very convenient!

### What to edit
Primarily edit TSV files in study folders, to fix the per-run metadata.

The "Condition" field is auto-generated, add custom conditions to files in [curation/run_descriptions](https://github.com/WormBase/wbps-expression/tree/master/curation/run_descriptions) folder.

Skip runs by adding them to the [curation/skipped_runs](https://github.com/WormBase/wbps-expression/tree/master/curation/skipped_runs) folder.

Don't edit YAMLs apart from removing contrasts if there are too many: other changes will be lost.

There are also a few places with essentially corner-case curation, scattered around the source code:
- characteristics get standardised after retrieving them from RNASeq-er through a bunch of regex-based heuristics centered around parasite specific stuff like life stages
- PubMed ids not in ENA or GEO
- Sample accessions that do not correspond to biological replicates (because e.g. runs with different conditions have been incorrectly marked with the same sample ID by the submitter)
- Studies with fewer than six runs that are nevertheless worth including, e.g. when that's the best data for the species

### Editing tsv files
#### Command line environment
See the [scripts](https://github.com/WormBase/wbps-expression/tree/master/scripts) folder for utilities that pick out specific columns (good with `paste`) or transpose the rows and columns (good with `grep` or `sed`). Rerun the pipeline changing the file to get it in a more "standard" format.

#### Online - GitHub's editing tools
You can use GitHub's editing tools. Then make a commit and open a PR - Travis will run the build with the checks.
This is great for very small changes, but the editor doesn't help you much, it's just a text box. You can paste it somewhere more convenient, Google Spreadsheets didn't like tabs when I tried but Excel worked well. 

## Development
### Install
Not automatised, but should present no difficulties. This software can run both in a personalized or cluster computing environment.

#### Instructions
Install R, and DESeq2. Clone the repository and install all the Perl modules. Write a wrapper script similar to [bin/run-ebi.pl](https://github.com/WormBase/wbps-expression/blob/master/bin/run-ebi.pl), hooking up libraries, inputs, and outputs.

#### Perl modules
You could do this:
```
cpanm -v --installdeps --notest .
```

#### R
Unfortunately DESeq2 pulls down a lot of dependencies: interfacing C++ code, plotting, etc. Install BioConductor, and then install DESeq2 using BioConductor.
```
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2"))
```
