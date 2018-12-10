requires 'Data::Dumper';
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
requires 'Text::MultiMarkdown';

test_requires 'File::Temp';
test_requires 'Test::MockModule';
test_requires 'Test::More';

on 'develop' => sub {
    requires 'Statistics::R';
    suggests 'Smart::Comments';
};
