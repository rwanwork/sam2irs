#!/usr/bin/env perl
#    score-irs.pl -- Score the introns in an IRS file, while taking the top N.
#    Copyright (C) 2019  Raymond Wan
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

use AppConfig;
use AppConfig::Getopt;
use Pod::Usage;


########################################
##  Important constants
########################################


########################################
##  Important functions
########################################



########################################
##  Important variables
########################################

##  Arguments provided by the user
my $verbose_arg = 0;
my $first_arg = 0;

my %duplicate_check;


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

$config -> define ("first", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  How many to keep
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

$verbose_arg = 0;
if ($config -> get ("verbose")) {
  $verbose_arg = 1;
}

if (!defined ($config -> get ("first"))) {
  printf STDERR "EE\tThe option --first requires a number.\n";
  exit (1);
}
$first_arg = $config -> get ("first");


########################################
##  Summarize the settings
########################################

if ($verbose_arg) {
  printf STDERR "II\tFirst file:  %s\n", $first_arg;
}


########################################
##  Read in the first file
########################################

open (my $fp, "<", $first_arg) or die "EE\tCould not open $first_arg for reading.\n";
while (<$fp>) {
  my $line = $_;
  chomp $line;
  
  my ($chr_tmp, $source_gtf, $type_gtf, $start_gtf, $end_gtf, $score_gtf, $strand_gtf, $phase_gtf, $attributes_gtf) = split /\t/, $line;

  if ($type_gtf ne "intron") {
    next;
  }

  my $key = $chr_tmp."-".$start_gtf."-".$end_gtf.$strand_gtf;
  
  if (defined ($duplicate_check{$key})) {
    printf STDERR "EE\tUnexpected error -- possible duplicates in input file (%s).\n", $key;
    exit (1);
  }
  $duplicate_check{$key} = 1;
}
close ($fp);


########################################
##  Read in the second file from STDIN
########################################

while (<STDIN>) {
  my $line = $_;
  chomp $line;
  
  my ($chr_tmp, $source_gtf, $type_gtf, $start_gtf, $end_gtf, $score_gtf, $strand_gtf, $phase_gtf, $attributes_gtf) = split /\t/, $line;

  if ($type_gtf ne "intron") {
    next;
  }

  my $key = $chr_tmp."-".$start_gtf."-".$end_gtf.$strand_gtf;
  
  if (defined ($duplicate_check{$key})) {
    $score_gtf = 0;
    $attributes_gtf = "";
  
    my $new_line = join ("\t", $chr_tmp, $source_gtf, $type_gtf, $start_gtf, $end_gtf, $score_gtf, $strand_gtf, $phase_gtf, $attributes_gtf);
    printf STDOUT "%s\n", $new_line;
  }
}



=pod

=head1 NAME

intersect-irs.pl -- Scores an IRS file and takes the top N records (the records themselves are unchanged).


=head1 AUTHOR

Raymond Wan <rwan.work@gmail.com>

=head1 COPYRIGHT

Copyright (C) 2016-2019, Raymond Wan, All rights reserved.


