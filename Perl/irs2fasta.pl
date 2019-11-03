#!/usr/bin/env perl
#    irs2fasta -- Takes the set of introns and extracts a set of FASTA
#      sequences surrounding them.
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

my $WINDOW_SIZE_CONST = 200;


########################################
##  Important functions
########################################

##  Reverse complement a sequence that may or may not be masked.
sub DNAReverseComplement {
  my $dna = shift;
  my $revcomp = reverse ($dna);

  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}


########################################
##  Important variables
########################################

##  Arguments provided by the user
my $verbose_arg = 0;
my $genome_dir_arg = "";
my $genome_fn_arg = "";
my $region_arg = "";
my $hardmask_arg = 0;
my $unmask_arg = 0;
my $window_size_arg = $WINDOW_SIZE_CONST;

my @gtf_array;  ##  GTF file in an array
my %chr_list;  ##  List of chromosomes to process
my %all_chrs;  ##  All chromosomes


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

$config -> define ("genomedir", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  Genome directory
$config -> define ("genomefn", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  Genome file
$config -> define ("region", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  Indicate what to extract
$config -> define ("winsize", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=i"
});                        ##  The window size to use
$config -> define ("unmask!", {
  ARGCOUNT => AppConfig::ARGCOUNT_NONE
});                        ##  Indicate sequence should not be masked
$config -> define ("hardmask!", {
  ARGCOUNT => AppConfig::ARGCOUNT_NONE
});                        ##  Indicate sequence should be hard-masked
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

if ((defined ($config -> get ("genomedir"))) && (defined ($config -> get ("genomefn")))) {
  printf STDERR "EE\tYou cannot use both the --genomedir and --genomefn options.\n";
  exit (1);
}

if ((!defined ($config -> get ("genomedir"))) && (!defined ($config -> get ("genomefn")))) {
  printf STDERR "EE\tYou need to use either the --genomedir and --genomefn options.\n";
  exit (1);
}

if (defined ($config -> get ("genomedir"))) {
  $genome_dir_arg = $config -> get ("genomedir");
}

if (defined ($config -> get ("genomefn"))) {
  $genome_fn_arg = $config -> get ("genomefn");
}

if (!defined ($config -> get ("region"))) {
  printf STDERR "EE\tPlease indicate the region to be extracted using the --region option.  Choose one of 'intron', 'upstream', or 'downstream' (case-sensitive).\n";
  exit (1);
}
$region_arg = $config -> get ("region");

if (($region_arg ne "intron") &&
    ($region_arg ne "upstream") &&
    ($region_arg ne "downstream")) {
  printf STDERR "EE\tFor the region to be extracted, choose one of 'intron', 'upstream', or 'downstream' (case-sensitive).\n";
  exit (1);
}

if (defined ($config -> get ("winsize"))) {
  $window_size_arg = $config -> get ("winsize");
  if ($window_size_arg <= 0) {
    printf STDERR "EE\tThe window size cannot be 0 or negative.\n";
    exit (1);
  }
}

$unmask_arg = 0;
if ($config -> get ("unmask")) {
  $unmask_arg = 1;
}

$hardmask_arg = 0;
if ($config -> get ("hardmask")) {
  $hardmask_arg = 1;
}


########################################
##  Summarize the settings
########################################

if ($verbose_arg) {
  if (length ($genome_dir_arg) != 0) {
    printf STDERR "II\tGenome directory:  %s\n", $genome_dir_arg;
  }
  if (length ($genome_fn_arg) != 0) {
    printf STDERR "II\tGenome file:  %s\n", $genome_fn_arg;
  }
  printf STDERR "II\tRegions to extract:  %s\n", $region_arg;
  if ($region_arg ne "intron") {
    printf STDERR "II\tWindow size to extract:  %u\n", $window_size_arg;
  }
  if ($unmask_arg == 0) {
    printf STDERR "II\tUnmask sequence:  No\n";
  }
  else {
    printf STDERR "II\tUnmask sequence:  Yes\n";
  }
  if ($hardmask_arg == 0) {
    printf STDERR "II\tHard-mask sequence:  No\n";
  }
  else {
    printf STDERR "II\tHard-mask sequence:  Yes\n";
  }
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
  
  ##  Maintaining a list of chromosomes to process (and how many times it was encountered)
  if (!defined $chr_list{$chr_tmp}) {
    $chr_list{$chr_tmp} = 0;
  }
  $chr_list{$chr_tmp}++;
  
  ##  Store the GTF record in an array
  $gtf_array[$gtf_pos] = $line;
  
  ##  Move to the next array position
  $gtf_pos++;
}

my $num_chromosomes = keys %chr_list;
if ($verbose_arg) {
  printf STDERR "II\tGTF records read:  %u\n", $gtf_total;
  printf STDERR "II\t  Intron records:  %u\n", $gtf_intron_total;
  printf STDERR "II\tNumber of chromosomes in GTF file:  %u\n", $num_chromosomes;
}


########################################
##  If genome given as a single compressed file, then read all of it in
########################################

##  Create a temporary file, which is deleted automatically when it goes out of scope
if (defined ($config -> get ("genomefn"))) {
  ##  Define the path to the compressed genome
  my $curr_genome_fn = $genome_fn_arg;
  if (! -e $curr_genome_fn) {
    printf STDERR "EE\tThe input file %s could not be found!\n", $curr_genome_fn;
    exit (1);
  }
   
  ##  Uncompress the chromosome
  my $tmp_fn = "";
  if (($curr_genome_fn =~ /\.gz$/) && ($verbose_arg)) {
    printf STDERR "II\tUncompressing %s to %s...\n", $curr_genome_fn, $tmp_fn;
    my $tmp_fh = File::Temp -> new ();
    $tmp_fn = $tmp_fh -> filename;
    `gunzip -c $curr_genome_fn >$tmp_fn`;
  }
  else {
    $tmp_fn = $curr_genome_fn;
  }
   
  ##  Read in the chromosome
  my $chr = "";
  my $chr_name = "";
  my $chr_count = 0;
  open (my $fp, "<", $tmp_fn) or die "EE\tCould not open $tmp_fn for input!\n";
  while (<$fp>) {
    my $line = $_;
    chomp $line;
    
    if ($line =~ /^>(.+)$/) {
      if ($. != 1) {
        $all_chrs{$chr_name} = $chr;
        $chr = "";
        $chr_count++;
      }
      $chr_name = $1;
    }
    else {
      $chr = $chr.$line;
    }
  }
  close ($fp);
  $all_chrs{$chr_name} = $chr;
  $chr_count++;
  
  if ($verbose_arg) {
    printf STDERR "II\tNumber of chromosomes in genome file:  %u\n", $chr_count;
  
#     foreach my $key (sort (keys %all_chrs)) {
#       printf STDERR "II\t%s\t%u\n", $key, length ($all_chrs{$key});
#     }
  }
}


########################################
##  Main program body -- applied to each chromosome independently
########################################

my $curr_genome_fn = "";
foreach my $key (sort (keys %chr_list)) {
#   printf STDERR "II\tProcessing %s...\n", $key;

  ##  Current chromosome
  my $chr = "";

  ##  Since Perl arrays are 0-based but coordinates are 1-based, we add a dummy value in the beginning
#   $chr = "X";  
  
  if (length ($genome_dir_arg) != 0) {
    ##  Create a temporary file, which is deleted automatically when it goes out of scope
    my $tmp_fh = File::Temp -> new ();
    my $tmp_fn = $tmp_fh -> filename;

    ##  Define the path to the compressed genome
    $curr_genome_fn = $genome_dir_arg."/".$key.".fa.gz";
    if (! -e $curr_genome_fn) {
      printf STDERR "EE\tThe input file %s could not be found!\n", $curr_genome_fn;
      exit (1);
    }
   
    ##  Uncompress the chromosome
    if ($verbose_arg) {
      printf STDERR "II\tUncompressing %s to %s...\n", $curr_genome_fn, $tmp_fn;
    }
    `gunzip -c $curr_genome_fn >$tmp_fn`;
   
    ##  Read in the chromosome
    open (my $fp, "<", $tmp_fn) or die "EE\tCould not open $tmp_fn for input!\n";
    my $header = <$fp>;
    while (<$fp>) {
      my $line = $_;
      chomp $line;
    
      if ($line =~ /^>/) {
        printf STDERR "EE\tUnexpected FASTA header within chromosome on line %u!  There should be only one chromosome per FASTA file!\n", $.;
        exit (1);
      }
    
      $chr = $chr.$line;
    }
    close ($fp);
    printf STDERR "II\tLength of chromosome %s:  %u\n", $key, length ($chr);
  }
  elsif (length ($genome_fn_arg) != 0) {
    if (!defined ($all_chrs{$key})) {
      printf STDERR "EE\tThe sequence for chromosome %s was not read in from file!\n", $key;
      exit (1);
    }
    
    $chr = $chr.$all_chrs{$key};
  }
  else {
    printf STDERR "EE\tYou need to use either the --genomedir and --genomefn options.\n";
    exit (1);
  }
  
  ##  Process each GTF record
  for (my $i = 0; $i < $gtf_pos; $i++) {
    my $line = $gtf_array[$i];

    my ($seqid_gtf, $source_gtf, $type_gtf, $start_gtf, $end_gtf, $score_gtf, $strand_gtf, $phase_gtf, $attributes_gtf) = split /\t/, $line;
    
    if ($seqid_gtf ne $key) {
      next;
    }

    ##  Create a name for the FASTA record
    my $name = sprintf (">%s_%u_%u", $seqid_gtf, $start_gtf, $end_gtf);

    ##  Shift the positions over by 1 since Perl arrays are 0-based
    $start_gtf = $start_gtf - 1;
    $end_gtf = $end_gtf - 1;
    
    my $window_start = 0;
    my $window_end = 0;
    my $window_seq = "";

    if ($region_arg eq "intron") {
      $window_start = $start_gtf;
      $window_end = $end_gtf;
    }
    elsif ($region_arg eq "upstream") {
      $window_start = $start_gtf - $window_size_arg;
      $window_end = $start_gtf - 1;
      if ($strand_gtf eq "-") {
        $window_start = $end_gtf + 1;
        $window_end = $end_gtf + $window_size_arg;
      }
    }
    elsif ($region_arg eq "downstream") {
      $window_start = $end_gtf + 1;
      $window_end = $end_gtf + $window_size_arg;
      if ($strand_gtf eq "-") {
        $window_start = $start_gtf - $window_size_arg;
        $window_end = $start_gtf - 1;
      }
    }
    
    if ($window_start < 0) {
      printf STDERR "WW\tChromosome %s exceeded position 0.\n", $key;
      printf STDERR "WW\t  [%s]\n", $line;
      $window_start = 0;
    }
    if ($window_end > length ($chr)) {
      printf STDERR "WW\tChromosome %s exceeded last position.\n", $key;
      printf STDERR "WW\t  [%s]\n", $line;
      $window_end = length ($chr);
    }
    
    $window_seq = substr ($chr, $window_start, $window_end - $window_start + 1);
    if ($verbose_arg) {
#       printf STDERR "II\t[%s]  Extracted from %u until %u [%s]\n", $strand_gtf, $window_start, $window_end, $window_seq;
    }
    
    ##  Reverse complement the sequence
    if ($strand_gtf eq "-") {
      $window_seq = DNAReverseComplement ($window_seq);
    }
    
    ##  Hard-mask the sequence, if requested
    if ($hardmask_arg) {
      $window_seq =~ s/a/N/gs;
      $window_seq =~ s/c/N/gs;
      $window_seq =~ s/g/N/gs;
      $window_seq =~ s/t/N/gs;
    }

    ##  Unmask the sequence, if requested
    if ($unmask_arg) {
      $window_seq = uc ($window_seq);
    }

    printf STDOUT "%s\n", $name;
    printf STDOUT "%s\n", $window_seq;
  }
}



=pod

=head1 NAME

irs2fasta.pl -- Extract the DNA sequence corresponding to each cluster.  Input format is a tab-separated GTF formatted file; output format is FASTA.

=head1 SYNOPSIS

B<irs2fasta.pl> <input.txt >output.fasta

=head1 DESCRIPTION

Extract the sequences.

=head1 OPTIONS

=over 5

=item --genomedir I<string>

Path to the directory of genomes.  Assumes they are B<compressed> and that their filenames matches the names used in the cluster file.

=item --unmask

Unmask the sequence by converting it all to upper case.  Of course, this implies that the original sequence was already masked.  There is no harm unmasking already masked sequences.

=item --verbose

Display verbose information about the execution of the program.

=item --help

Display this help message.

=back

=head1 EXAMPLE

=over 5

cat input.txt | ./irs2fasta.pl --genomedir ~/some-dir/

=back

=head1 AUTHOR

Raymond Wan <rwan.work@gmail.com>

=head1 COPYRIGHT

Copyright (C) 2016-2019, Raymond Wan, All rights reserved.


