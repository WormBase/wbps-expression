requires 'Data::Compare';
requires 'File::Basename';
requires 'File::Path';
requires 'File::Slurp';
requires 'JSON';
requires 'LWP';
requires 'List::MoreUtils';
requires 'List::Util';
requires 'Text::CSV';
requires 'XML::Simple';
requires 'YAML';
requires 'HTML::Template';
requires 'Statistics::R';

test_requires 'File::Temp';
test_requires 'Test::MockModule';
test_requires 'Test::More';


on 'develop' => sub {
    suggests 'Smart::Comments';
};
