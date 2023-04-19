ROOT_DIR=${EXPRESSION_CODE}

if [ -d "${PARASITE_PRODUCTION}/data/conda_envs/" ] ; then
  PATH=${PARASITE_PRODUCTION}/data/conda_envs/wbps-expression/bin:$PATH
  export PATH
else echo "Couldn't find ${PARASITE_PRODUCTION}/data/conda_envs/";
fi;

if [ 0 -eq $(R --slave --no-restore --file=- <<< 'installed.packages()' | grep -c ^DESeq2) ]; then
  echo "Your R doesn't have DESeq2 installed: " $(which R)
fi

export PERL5LIB="$ROOT_DIR/lib:$ROOT_DIR/local/lib/perl5:${PARASITE_SOFTWARE}/perl_modules/lib/perl5:${PERL5LIB}"