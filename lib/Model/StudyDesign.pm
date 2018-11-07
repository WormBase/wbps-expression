package Model::StudyDesign;
use Text::CSV qw/csv/;

sub new {
  my ($class, $characteristics_in_order, $characteristics_per_condition) = @_;
  return bless {
     characteristics_in_order => $characteristics_in_order // [],
     characteristics_per_condition => $characteristics_per_condition // {},
  }, $class;
}

sub from_tsv {
  my ($path) = @_;
  my %characteristics_per_condition;

  my ($header, @lines) = @{csv(in=>$path, sep_char => "\t")};
  my ($__, @characteristics) = @{$header};
  for my $line (@lines){
     my ($condition, @values) = @$line;
     my %h;
     @h{@characteristics} = @values;
     $characteristics_per_condition{$condition} = \%h;
  }
  return new(__PACKAGE__, \@characteristics, \%characteristics_per_condition);
}

sub to_tsv {
  my($self, $path) = @_;
  open(my $fh, ">", $path) or die $path;
  my @a = @{$self->{characteristics_in_order}};
  my %h = %{$self->{characteristics_per_condition}};
  print $fh join("\t", "Condition", @a)."\n";
  for my $condition (sort keys %h){
     my %cs = %{$h{$condition}};
     print $fh join("\t", $condition, @cs{@a})."\n";
  }
  close $fh; 
}
1;
