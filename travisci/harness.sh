#!/bin/bash
export PERL5LIB=$PWD/lib

echo "Running test suite"
echo "Using $PERL5LIB"

perl $PWD/ensembl-test/scripts/runtests.pl $PWD/t/
