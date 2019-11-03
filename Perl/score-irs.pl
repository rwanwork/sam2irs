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

##  Constant for normalization
my $NORMALIZER_CONSTANT = 100000000000;


########################################
##  Important functions
########################################

sub GetCount {
  my ($gtf_attributes) = @_;
  my $count = 0;
  
  if ($gtf_attributes =~ /count=(\d+);/) {
    $count = $1;
  }
  else {
    printf STDERR "EE\tCould not obtain the count in %s.\n", $gtf_attributes;
    exit (1);
  }

  return ($count);
}

sub GetWidth {
  my ($gtf_attributes) = @_;
  my $width = 0;
  
  if ($gtf_attributes =~ /width=(\d+);/) {
    $width = $1;
  }
  else {
    printf STDERR "EE\tCould not obtain the width in %s.\n", $gtf_attributes;
    exit (1);
  }

  return ($width);
}

sub GetAlignedBases {
  my ($gtf_attributes) = @_;
  my $aligned_bases = 0;
  
  if ($gtf_attributes =~ /aligned_bases=(\d+)$/) {
    $aligned_bases = $1;
  }
  elsif ($gtf_attributes =~ /aligned_bases=(\d+);/) {
    $aligned_bases = $1;
  }
  else {
    printf STDERR "EE\tCould not obtain the aligned bases in %s.\n", $gtf_attributes;
    exit (1);
  }

  return ($aligned_bases);
}


sub ScoreCmpFunction {
    my %record_a = %$a;
    my %record_b = %$b;
    my $first = $record_a{"score"};
    my $second = $record_b{"score"};
    if ($first < $second) { 
      return 1; 
    }
    elsif ($first == $second) { 
      return 0; 
    }
    elsif ($first > $second) { 
      return -1; 
    }
}


########################################
##  Important variables
########################################

##  Arguments provided by the user
my $verbose_arg = 0;
my $topn_arg = 0;

my @gtf_array;  ##  GTF file in an array
my %chr_list;  ##  List of chromosomes to process


########################################
##  Important global counter variables
########################################

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

$config -> define ("topn", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=i"
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

if (!defined ($config -> get ("topn"))) {
  printf STDERR "EE\tThe option --topn requires a number.\n";
  exit (1);
}
$topn_arg = $config -> get ("topn");


########################################
##  Summarize the settings
########################################

if ($verbose_arg) {
  printf STDERR "II\tTop N:  %s\n", $topn_arg;
  printf STDERR "II\tNormalization constant:  %s\n", $NORMALIZER_CONSTANT;
}


########################################
##  Read in GTF file
########################################

my $gtf_pos = 0;
my $gtf_intron_total = 0;
my $gtf_total = 0;
while (<STDIN>) {
  my $line = $_;
  chomp $line;
  
  $gtf_total++;

  my ($chr_tmp, $source_gtf, $type_gtf, $start_gtf, $end_gtf, $score_gtf, $strand_gtf, $phase_gtf, $attributes_gtf) = split /\t/, $line;

  if ($type_gtf ne "intron") {
    next;
  }
  $gtf_intron_total++;

  my $count = GetCount ($attributes_gtf);
  my $width = GetWidth ($attributes_gtf);
  my $aligned_bases = GetAlignedBases ($attributes_gtf);
  my $score = $count / $width / $aligned_bases * $NORMALIZER_CONSTANT;
  
  ##  Maintaining a list of chromosomes to process (and how many times it was encountered)
  if (!defined $chr_list{$chr_tmp}) {
    $chr_list{$chr_tmp} = 0;
  }
  $chr_list{$chr_tmp}++;
  
  ##  Store the GTF record in an array
  my $r = {score => $score, record => $line};
  push (@gtf_array, $r);
  
  ##  Move to the next array position
  $gtf_pos++;
}

my $num_chromosomes = keys %chr_list;
if ($verbose_arg) {
  printf STDERR "II\tGTF records read:  %u\n", $gtf_total;
  printf STDERR "II\t  Intron records:  %u\n", $gtf_intron_total;
  printf STDERR "II\tNumber of chromosomes in GTF file:  %u\n", $num_chromosomes;
}


my @sorted_gtf_array = sort ScoreCmpFunction @gtf_array;

for (my $i = 0; $i < $topn_arg; $i++) {
#   printf STDERR "%f\t[%s]\n", $sorted_gtf_array[$i]{"score"}, $sorted_gtf_array[$i]{"record"};
  printf STDOUT "%s\n", $sorted_gtf_array[$i]{"record"};
}


=pod

=head1 NAME

score-irs.pl -- Scores an IRS file and takes the top N records (the records themselves are unchanged).


=head1 AUTHOR

Raymond Wan <rwan.work@gmail.com>

=head1 COPYRIGHT

Copyright (C) 2016-2019, Raymond Wan, All rights reserved.


