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
requires 'Log::Any';
requires 'Log::Any::Adapter::Log4perl';
requires 'Log::Log4perl';
requires 'DateTime::Format::Strptime';
requires 'DateTime::Format::ISO8601::Format';
requires 'Sub::Throttle';
requires 'Archive::Zip';

test_requires 'File::Temp';
test_requires 'Test::Mock::LWP';
test_requires 'Test::More';


on 'develop' => sub {
    suggests 'Smart::Comments';
};
