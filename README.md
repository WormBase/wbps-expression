[![Build Status](https://travis-ci.org/WormBase/wbps-expression.svg?branch=master)](https://travis-ci.org/WormBase/wbps-expression)
# wbps-expression
Transcriptomic data for WormBase ParaSite

## Install

### Perl
```
cpanm -v --installdeps --notest .
```
Retrieving resources also requires `LWP::UserAgent` making requests over https - should be included in core.
### R - DESeq2
Unfortunately DESeq2 pulls down a lot of dependencies (interfacing C++ code, plotting, etc. )
- Install R
- Install DESeq2 using BioConductor
```
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2"))
```
