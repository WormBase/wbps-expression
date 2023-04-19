#!/usr/bin/bash
set -euo pipefail

ROOT_DIR=${EXPRESSION_CODE}

( cd $ROOT_DIR && [ -d "local" ] || carton )

if [ -d "${PARASITE_PRODUCTION}/data/conda_envs/" ] ; then
  PATH=${PARASITE_PRODUCTION}/data/conda_envs/wbps-expression/bin:$PATH
  export PATH
fi

if [ 0 -eq $(R --slave --no-restore --file=- <<< 'installed.packages()' | grep -c ^DESeq2) ]; then
  echo "Your R doesn't have DESeq2 installed: " $(which R)
fi

if [ "$#" -gt 0 ]; then
  PERL5LIB="$ROOT_DIR/lib:$ROOT_DIR/local/lib/perl5:${PARASITE_SOFTWARE}/perl_modules/lib/perl5:${PERL5LIB}" prove $( find "$ROOT_DIR/t/" -type f $(perl -e 'print join(" -o ", map {"-name *$_*"} @ARGV )' "$@")  )
else
  PERL5LIB="$ROOT_DIR/lib:$ROOT_DIR/local/lib/perl5:${PARASITE_SOFTWARE}/perl_modules/lib/perl5:${PERL5LIB}" prove -r $ROOT_DIR/t/
fi