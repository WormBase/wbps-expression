[![Build Status](https://travis-ci.org/WormBase/wbps-expression.svg?branch=master)](https://travis-ci.org/WormBase/wbps-expression)
# wbps-expression
Transcriptomic data for WormBase ParaSite

## Install
Clone the repository and install all the Perl modules. Nothing else is required, apart from one bit that also uses R.

### Perl
```
cpanm -v --installdeps --notest .
```
### R - DESeq2
Unfortunately DESeq2 pulls down a lot of dependencies (interfacing C++ code, plotting, etc. )
- Install R
- Install DESeq2 using BioConductor
```
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2"))
```
### Running the code
Not standardised yet. 
There's a wrapper, which defines example locations to folders the code needs:
- folder to download public resources to - from ENA, RNASeq-er, and so on
- where the curation is
- where to write files to
```
bin/run-ebi.sh <species> <assembly>
bin/checks-ebi.sh # shows what's wrong with studies that got skipped
```
### Curation
#### what gets created automatically, and what gets preserved
Characteristics for runs: `curation/studies/$study/$study.tsv`
What to name the contrasts: `curation/studies/$study/$study.yaml` unless you set an environment variable `RECREATE_ALL_CONFIGS=1`
Conditions: these are automatically created, but you can override them in `curation/run_descriptions/$species.tsv`

#### Example workflow
Run the code for a species, to get an automatic attempt at curation and some analyses. Commit it so that it's easy to refer to what was done.
Fix the characteristics. Then commit (maybe) and run the code again - it will fix condition names, and if all is good some more analyses will run.

If it's not clear how to curate a study, you can ignore it by not committing it. It will come up next time around when the code is ran. You can also exclude a study forever, by adding the accession to `curation/ignore_studies/$species.tsv`. 

When you're done with a species, make sure the state of Git is what you want it to be, push it to a branch and open a pull request. It will then be easy to compare what was done etc. You can also run the checks yourself locally and push to the master branch.
#### What the characteristics should be
Broadly, they should agree on conditions. The samples that represent the same condition should have the same characteristics, or else the analyses won't work.

#### How to conveniently edit the files
Copy the design files somewhere where you can edit them - Excel works well, so does online editing via Github.  You can also use programs - `scripts/transpose` swaps rows and columns, and then you can use `grep -v` to remove unwanted columns.
