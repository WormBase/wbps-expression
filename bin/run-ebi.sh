#!/usr/bin/bash
set -euo pipefail

if [ $# -lt 2 -o "$PARASITE_VERSION" == "" ] ; then
  echo "Usage: PARASITE_VERSION=<> $0 species assembly"
  exit 1
fi
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

( cd $ROOT_DIR && [ -d "local" ] || carton )

if [ -d "$ROOT_DIR/local/R-3.5.1/bin" ] ; then
  PATH=$ROOT_DIR/local/R-3.5.1/bin:$PATH
  export PATH
fi

hasDESeq2=$(R --slave --no-restore --file=- <<< 'installed.packages()' | grep ^DESeq2)
if [ 0 -eq $(R --slave --no-restore --file=- <<< 'installed.packages()' | grep -c ^DESeq2) ]; then
  echo "Your R doesn't have DESeq2 installed: " $(which R)
fi

ROOT_DIR=$ROOT_DIR PERL5LIB="$ROOT_DIR/lib:$ROOT_DIR/local/lib/perl5" ANALYSIS_VERBOSE=1 LOCALLY_CACHED_RESOURCE_VERBOSE=1 perl -e '
   use Production::Workflow;
   my $data_dir = "/nfs/nobackup/ensemblgenomes/wormbase/parasite/production/jbrowse/WBPS$ENV{PARASITE_VERSION}";
   my $src_dir = $ENV{ROOT_DIR};
   Production::Workflow->new("$data_dir/Resources", $src_dir, "$data_dir/Production")->do_everything(@ARGV);
' "$@"
