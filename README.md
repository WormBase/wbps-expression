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
### Example setup
Have a look at
```
bin/run-ebi.sh
```
This pulls in the depenencies using Carton, and refers to a local R installation. I installed R manually - it takes so long, and it's used for so little (one analysis) that there's no point in writing a provisioning script.
