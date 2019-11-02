#!/usr/bin/env perl
#    generate-negative.pl -- Generate a set of negative sequences based on
#      the proportion of bases in a positive sequence.
#    Copyright (C) 2016-2019  Raymond Wan
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FindBin qw ($Bin);
use lib $FindBin::Bin;  ##  Search the directory where the script is located

use diagnostics;
use strict;
use warnings;

use File::Temp;

use AppConfig;
use AppConfig::Getopt;
use Pod::Usage;


########################################
##  Important constants
########################################

my $WINDOW_SIZE = 50;


########################################
##  Important variables
########################################

##  Arguments provided by the user
my $window_size_arg = $WINDOW_SIZE;
my $gccount_arg = 0;

##  Important data structure
my %raw_counts;  ##  Hash of nucleotide counts
my @nucleotides_seen;  ##  Nucleotides seen
my @cumm_prob;  ##  Cummulative probabilities


########################################
##  Process arguments
########################################

##  Create AppConfig and AppConfig::Getopt objects
my $config = AppConfig -> new ({
           GLOBAL => {
             DEFAULT => undef,     ##  Default value for new variables
           }
  });

my $getopt = AppConfig::Getopt -> new ($config);

$config -> define ("winsize", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=i",
  DEFAULT => $WINDOW_SIZE
});                        ##  Window size
$config -> define ("gccount!", {
  ARGCOUNT => AppConfig::ARGCOUNT_NONE
});                        ##  Perform a GC count and exit
$config -> define ("verbose!", {
  ARGCOUNT => AppConfig::ARGCOUNT_NONE
});                        ##  Verbose output
$config -> define ("help!", {
  ARGCOUNT => AppConfig::ARGCOUNT_NONE
});                        ##  Help screen

##  Process the command-line options
$config -> getopt ();


########################################
##  Validate the settings
########################################
if ($config -> get ("help")) {
  pod2usage (-verbose => 0);
  exit (1);
}

if (!defined ($config -> get ("winsize"))) {
  printf STDERR "EE\tThe option --winsize is required with a window size.\n";
  exit (1);
}
$window_size_arg = $config -> get ("winsize");

$gccount_arg = 0;
if ($config -> get ("gccount")) {
  $gccount_arg = 1;
}

printf STDERR "II\tWindow size:  %u\n", $window_size_arg;
printf STDERR "II\tGC count only:  ";
if ($gccount_arg == 0) {
  printf STDERR "No\n";
}
else {
  printf STDERR "Yes\n";
}


########################################
##  Initialize hash for keeping track of base frequencies
########################################

##  Add bases here to include them in the raw_counts
$raw_counts{"A"} = 0;
$raw_counts{"C"} = 0;
$raw_counts{"G"} = 0;
$raw_counts{"T"} = 0;


########################################
##  Read in FASTA sequences
########################################

my $records_total_lines = 0;
my $dropped_bases = 0;
my $total_bases = 0;
while (<STDIN>) {
  my $header = $_;
  chomp $header;
  my $seq = <STDIN>;  
  chomp $seq;

  ##  Validate the header is in FASTA format
  if ($header !~ /^>/) {
    printf STDERR "EE\tProblem with FASTA header on line %u.\n", $.;
    exit (1);
  }
  
  ##  Count the number of records seen
  $records_total_lines++;

  ##  Unmask the sequence
#   $seq = uc ($seq);

  my @bases = split //, $seq;
  
  ##  Ensure the value to --winsize is the same as the sequence length
  if (scalar (@bases) != $window_size_arg) {
    printf STDERR "EE\tArgument provided with --winsize differs in value from the actual sequence length of %u on line %u!\n", scalar (@bases), $.;
    exit (1);
  }
  
  ##  Count each base.  Bases which are not in the hash are dropped
  for (my $i = 0; $i < scalar (@bases); $i++) {
    if (!defined $raw_counts{$bases[$i]}) {
      $dropped_bases++;
    }
    else {
      $raw_counts{$bases[$i]}++;
    }
    $total_bases++;
  }
}

my $kept_bases = $total_bases - $dropped_bases;
printf STDERR "II\tRecords read:  %u\n", $records_total_lines;
printf STDERR "II\tTotal bases:  %u\n", $total_bases;
printf STDERR "II\t  Bases dropped:  %u\n", $dropped_bases;
printf STDERR "II\t  Bases kept:  %u\n", $kept_bases;

my $k = 0;  ##  Nucleotides seen so far (array position)
foreach my $key (sort (keys %raw_counts)) {
  $nucleotides_seen[$k] = $key;  ##  Nucleotides seen so far
  $cumm_prob[$k] = $raw_counts{$key} / $kept_bases;  ##  Probability of that base
  if ($config -> get ("verbose")) {
    printf STDERR "II\t  %s (%u):  %u (%.7f)\n", $key, $k, $raw_counts{$key}, $cumm_prob[$k];
  }
  $k++;
}

##  Simplify the probabilities so that P("C") = P("G") and P("A") = P("T")
if (($nucleotides_seen[0] eq "A") && ($nucleotides_seen[3] eq "T")) {
  $cumm_prob[0] = ($cumm_prob[0] + $cumm_prob[3]) / 2;
  $cumm_prob[3] = $cumm_prob[0];
}
else {
  printf STDERR "EE\tBases A and T are not in positions 0 and 3 of the array!\n";
  exit (1);
}

if (($nucleotides_seen[1] eq "C") && ($nucleotides_seen[2] eq "G")) {
  $cumm_prob[1] = ($cumm_prob[1] + $cumm_prob[2]) / 2;
  $cumm_prob[2] = $cumm_prob[1];
}
else {
  printf STDERR "EE\tBases C and G are not in positions 1 and 2 of the array!\n";
  exit (1);
}

for ($k = 0; $k < scalar (@cumm_prob); $k++) {
  printf STDERR "II\t  Simplified %s (%u):  %.7f\n", $nucleotides_seen[$k], $k, $cumm_prob[$k];
}

##  Convert the probabilities to cummulative probabilities
for ($k = 1; $k < scalar (@cumm_prob); $k++) {
  $cumm_prob[$k] = $cumm_prob[$k] + $cumm_prob[$k - 1];
}

##  Due to rounding, it is possible that $cumm_prob[scalar (@cumm_prob) - 1] is not exactly 1
$cumm_prob[scalar (@cumm_prob) - 1] = 1;

if ($config -> get ("verbose")) {
  for ($k = 0; $k < scalar (@cumm_prob); $k++) {
    printf STDERR "II\t%u\t%.7f\t%s\n", $k, $cumm_prob[$k], $nucleotides_seen[$k];
  }
  printf STDERR "\n";
}

##  If just counting GC content, then exit now and don't print out FASTA sequences
if ($gccount_arg == 1) {
  exit;
}


########################################
##  Create random FASTA sequences
########################################

##  Create the same number of records as the input
for (my $i = 0; $i < $records_total_lines; $i++) {
  printf STDOUT ">Negative-%u\n", $i;
  for (my $j = 0; $j < $window_size_arg; $j++) {
    my $base = "";
    my $value = rand ();  ##  Generate a random value in the range [0, 1]
    for ($k = 0; $k < scalar (@cumm_prob); $k++) {
      if ($value < $cumm_prob[$k]) {
        $base = $nucleotides_seen[$k];
        last;
      }
    }
#     printf STDERR "%u\t%.7f\t%s\n", $i, $value, $base;

    ##  Try to catch the very small chance that $base has no value because we exited the loop
    ##    with a ($value < 1) but ($value > $cumm_prob[scalar (@cumm_prob) - 1]).  This would be very 
    ##    close to 1.  Since we set the maximum probability above to exactly 1, this shouldn't happen.
    ##    So, we should never exit here.
    if (length ($base) == 0) {
      printf STDERR "EE\tUnexpected error -- \$base has no value!\n";
      exit (1);
    }
    printf STDOUT "%s", $base;
  }
  printf STDOUT "\n";
}


=pod

=head1 NAME

generate-negative.pl -- Generate random negative sequences in FASTA format.

=head1 SYNOPSIS

B<generate-negative.pl> <input.fasta >output.fasta

=head1 DESCRIPTION

Generate random negative sequences with the same GC content as the positive sequences.  The length of the sequences generated are all equal to those of the positive sequences.  If masked sequences are to be handled, then:

=over 5

=item  1. The call to uc () should be removed.

=item  2. Lower cased versions of {A, C, G, T} should be added to the hash of counts.

=back

=head1 ALGORITHM

This script executes the following steps:

=over 5

=item  1. Initialize the hash of counts (if other bases are to be recorded, then add them here).

=item  2. Read in the input positive sequences in FASTA format.  Those that have not been previously initialized in the hash are dropped.

=item  3. Convert the base probabilities to cummulative probabilities in the range [0, 1].  (The largest probability is forcibly set to 1 exactly to prevent rounding errors.  i.e., just in case the value 0.9999999... is less than 1 but greater than the largest cummulative probability.)

=item  4. If we're just counting GC content, exit here.

=item  5. Create random FASTA sequences.

=back

=head1 OPTIONS

=over 5

=item --winsize I<number>

Size of the window of sequences to extract.

=item --gccount

Print out the GC content statistics and exit.

=item --verbose

Display verbose information about the execution of the program.

=item --help

Display this help message.

=back

=head1 EXAMPLE

=over 5

cat input.fasta | ./generate-negative.pl >output.fasta

=back

=head1 LIMITATIONS

=over 5

=item * Bases that have been hard-masked to "N" are ignored from both the counting and the total.

=item * Bases that have been soft-masked to {"a", "c", "g", "t"} are also ignored.  A call to uc () is needed if they should be counted.

=back

=head1 AUTHOR

Raymond Wan <rwan.work@gmail.com>

=head1 COPYRIGHT

Copyright (C) 2016-2019, Raymond Wan, All rights reserved.


