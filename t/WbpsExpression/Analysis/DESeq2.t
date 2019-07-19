use strict;
use warnings;
use Test::More;
use WbpsExpression::Analysis::DESeq2;

WbpsExpression::Analysis::DESeq2::establish_r_session_and_set_r_and_deseq_versions();
ok($WbpsExpression::Analysis::DESeq2::R_GLOBAL->is_started, "Can start R");
is($WbpsExpression::Analysis::DESeq2::R_GLOBAL->get(1), 1, "Can talk to R");
done_testing;
