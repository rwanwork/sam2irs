#!/usr/bin/env perl
#    sam2irs.pl -- Converts a set of mapped reads in SAM format to a 
#      table of intron retention scores
#
#    Copyright (C) 2016-2020  Raymond Wan
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

my $UNINIT_FREQUENCY = 0;
my $UNINIT_GENE = "";

##  An initialized base is an exon -- that is, a non-negative integer [0, ...]
##    A base associated with an exon is indicated by a non-negative integer 
##    pointing to the GTF record
my $UNINIT_BASE = -2;
my $INTRON_BASE = -1;

my $CLASSIFY_SAM_UNKNOWN = 0;
my $CLASSIFY_SAM_UNALIGNED = 1;
my $CLASSIFY_SAM_EXON = 2;
my $CLASSIFY_SAM_INTRON = 3;
my $CLASSIFY_SAM_BOTH = 4;

my $BOOLEAN_TRUE = 1;
my $BOOLEAN_FALSE = 0;

my $VERBOSE_NONE = 0;
my $VERBOSE_SUMMARY = 1;
my $VERBOSE_CHR = 2;
my $VERBOSE_ALL = 3;


########################################
##  Important functions
########################################

sub GetGeneName {
  my ($gtf_attributes) = @_;
  my $gene_name = "";
  
  if ($gtf_attributes =~ /gene_id \"([^"]+)\"/) {
    $gene_name = $1;
  }
  else {
    printf STDERR "EE\tCould not match gene name in %s.\n", $gtf_attributes;
    exit (1);
  }

  return ($gene_name);
}


########################################
##  Important variables
########################################

##  Arguments provided by the user
my $verbose_level_arg = $VERBOSE_NONE;
my $chrlist_arg = "";
my $gtf_arg = "";
my $sam_record_arg = "";
my $exons_arg = 0;
my $log_fn_arg = "";

##  Arrays of input files
my @sam_input;  ##  All SAM records (= # of lines in the SAM file)
my @classify_sam_input;  ##  Classification of the SAM records (equal in length to @sam_input)
my @gtf_input;  ##  All GTF records (= # of lines in the GTF file)

##  Arrays for outputing
my @introns_output;
my $num_introns_output = 0;
my @exons_output;
my $num_exons_output = 0;

##  Hashes of input files
my %chrlist_hash_input;  ##  Hash of chromosomes to process (= # of chromosomes)
my %chrlist_sam_acc;
my %chrlist_gtf_acc;

my $total_bases = 0;


########################################
##  Important global counter variables
########################################

my $total_sam = 0;  ##  Number of records in the SAM file

my $overlap_gtf = 0;
my $not_overlap_gtf = 0;
my $exon_gtf = 0;
my $not_exon_gtf = 0;

my $unavailable_cigar_sam = 0;  ##  Number of unavailable SAM records (ignored)
my $n_cigar_sam = 0;  ##  SAM records with 'N' (ignored)
my $not_n_cigar_sam = 0;  ##  SAM records without 'N' (used)


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

$config -> define ("chrlist", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});
$config -> define ("gtf", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  Indicate the GTF file to read
$config -> define ("samrecord", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  The record file for SAM
$config -> define ("exons!", {
  ARGCOUNT => AppConfig::ARGCOUNT_NONE
});                        ##  Flag to include exons
$config -> define ("verbose", {
  DEFAULT  => 0,
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=i"
});                        ##  Verbose output
$config -> define ("log", {
  ARGCOUNT => AppConfig::ARGCOUNT_ONE,
  ARGS => "=s"
});                        ##  The log file for verbose levels >= $VERBOSE_CHR
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

$verbose_level_arg = $VERBOSE_NONE;
if (defined ($config -> get ("verbose"))) {
  if ($config -> get ("verbose") == $VERBOSE_NONE) {
    $verbose_level_arg = $VERBOSE_NONE;
  }
  elsif ($config -> get ("verbose") == $VERBOSE_SUMMARY) {
    $verbose_level_arg = $VERBOSE_SUMMARY;
  }
  elsif ($config -> get ("verbose") == $VERBOSE_CHR) {
    $verbose_level_arg = $VERBOSE_CHR;
  }
  elsif ($config -> get ("verbose") == $VERBOSE_ALL) {
    $verbose_level_arg = $VERBOSE_ALL;
  }
  else {
    printf STDERR "EE\tInvalid verbose level provided.  Only one of {0, 1, 2, 3} is allowed.\n";
    exit (1);
  }
}

if (!defined ($config -> get ("chrlist"))) {
  printf STDERR "EE\tThe option --chrlist requires a list of chromosomes (i.e., format is genome sizes file).\n";
  exit (1);
}
$chrlist_arg = $config -> get ("chrlist");

if (!defined ($config -> get ("gtf"))) {
  printf STDERR "EE\tThe option --gtf requires a GTF filename.\n";
  exit (1);
}
$gtf_arg = $config -> get ("gtf");

if (! (-e $gtf_arg)) {
  printf STDERR "EE\tThe GTF file %s does not exist!\n", $gtf_arg;
  exit (1);
}

if (defined ($config -> get ("samrecord"))) {
  $sam_record_arg = $config -> get ("samrecord");
}

if (defined ($config -> get ("log"))) {
  $log_fn_arg = $config -> get ("log");
  if (($config -> get ("verbose") == $VERBOSE_NONE) || ($config -> get ("verbose") == $VERBOSE_SUMMARY)) {
    $verbose_level_arg = $VERBOSE_CHR;
  }
}
else {
  if ($verbose_level_arg >= $VERBOSE_CHR) {
    printf STDERR "EE\tFilename for log file required with a verbose level >= %u.\n", $VERBOSE_CHR;
    exit (1);
  }
}

$exons_arg = 0;
if ($config -> get ("exons")) {
  $exons_arg = 1;
}


########################################
##  Summarize the settings
########################################

if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "II\tFilenames\n";
  printf STDERR "II\t  Input chromosome list:  %s\n", $chrlist_arg;
  printf STDERR "II\t  Input GTF file:  %s\n", $gtf_arg;
  printf STDERR "II\t  Record of SAM file:  ";
  if (defined ($config -> get ("samrecord"))) {
      printf STDERR "%s\n", $sam_record_arg;
  }
  else {
      printf STDERR "Not provided\n";
  }
  printf STDERR "II\t  Log file:  ";
  if (defined ($config -> get ("log"))) {
      printf STDERR "%s\n", $log_fn_arg;
  }
  else {
      printf STDERR "Not provided\n";
  }
  if ($exons_arg) {
    printf STDERR "II\tOutput exons:  Yes\n";
  }
  else {
    printf STDERR "II\tOutput exons:  No\n";
  }
  printf STDERR "\n";
}


########################################
##  Create log file
########################################

my $log_fp;
if (length ($log_fn_arg) != 0) {
  open ($log_fp, ">", $log_fn_arg) or die "EE\tCould not open $log_fn_arg for output!\n";
}


########################################
##  Read in the entire GTF file, discarding records that are not relevant 
##    to our data set, as well as predicted genes
########################################

my $total_gtf = 0;
my $gtf_input_pos = 0;
my $gm_gtf = 0;

my $exon_count = 0;
my $cds_count = 0;
my $start_codon_count = 0;
my $stop_codon_count = 0;
my $other_count = 0;

open (my $gtf_fp, "<", $gtf_arg) or die "EE\tCould not open $gtf_arg for input!\n";
while (<$gtf_fp>) {
  my $line = $_;
  chomp $line;
  
  ##  Handle the rare chance that there are blank lines at the end of the GTF file
  if (length ($line) == 0) {
    next;
  }
  
  $total_gtf++;
  
  my ($seqname_tmp, $source_tmp, $feature_tmp, $start_tmp, $end_tmp, $score_tmp, $strand_tmp, $frame_tmp, $attribute_tmp) = split /\t/, $line;
 
  ##  Drop predicted genes
  if ($attribute_tmp =~ /gene_id "Gm\d+"/) {
    $gm_gtf++;
    next;
  }
  
  ##  Count the number of features taken; only mark exons on the genome
  if ($feature_tmp eq "exon") {
    $exon_count++;
    
    ##  Store GTF in an array
    $gtf_input[$gtf_input_pos] = join ("\t", $seqname_tmp, $source_tmp, $feature_tmp, $start_tmp, $end_tmp, $score_tmp, $strand_tmp, $frame_tmp, $attribute_tmp);
    $gtf_input_pos++;  
  }
  elsif ($feature_tmp eq "CDS") {
    $cds_count++;
  }
  elsif ($feature_tmp eq "start_codon") {
    $start_codon_count++;
  }
  elsif ($feature_tmp eq "stop_codon") {
    $stop_codon_count++;
  }
  else {
    $other_count++;
  }
}
close ($gtf_fp);

if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "II\tProcessed GTF file (%s, %u records):\n", $gtf_arg, $total_gtf;
  printf STDERR "II\t  Retained ('exon'):\t\t%u\n", $gtf_input_pos;
  printf STDERR "II\t  Discarded:\n";
  printf STDERR "II\t    Gene names with 'Gm':\t%u\n", $gm_gtf;
  printf STDERR "II\t    # of 'CDS':\t\t\t%u\n", $cds_count;
  printf STDERR "II\t    # of 'start_codon':\t\t%u\n", $start_codon_count;
  printf STDERR "II\t    # of 'stop_codon':\t\t%u\n", $stop_codon_count;
  printf STDERR "II\t    # of Others:\t\t%u\n", $other_count;
  printf STDERR "\n";
}


########################################
##  Read in the entire SAM file, keeping every record since 
##    we may be using the --samrecord option
########################################

my $sam_header_lines = 0;
my $sam_total_lines = 0;
my $sam_invalid_lines = 0;

$total_sam = 0;
while (<STDIN>) {
  my $line = $_;
  chomp $line;
  
  $sam_total_lines++;
  
  ##  Ignore the SAM header
  if ($line =~ /^@/) {
    $sam_header_lines++;
    next;
  }
  
  ##  Handle the rare chance that there are blank lines at the end of the SAM file
  if (length ($line) == 0) {
    next;
  }
  
  my @tmp_array = split /\t/, $line;
  if (scalar (@tmp_array) < 11) {
    printf STDERR "WW\tLine %u of the SAM file is incomplete with only %u columns [%s].\n", $sam_total_lines, scalar (@tmp_array), $line;
    $sam_invalid_lines++;
    next;
  }

  my ($qname_tmp, $flag_tmp, $rname_tmp, $pos_tmp, $mapq_tmp, $cigar_tmp, $rnext_tmp, $pnext_tmp, $tlen_tmp, $seq_tmp, $qual_tmp, $others_fields_tmp) = split /\t/, $line;

  ##  Store only the alignment in an array
  $sam_input[$total_sam] = join ("\t", $flag_tmp, $rname_tmp, $pos_tmp, $cigar_tmp);
  
  $total_sam++;
}

if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "II\tSAM file (STDIN, %u records):\n", $sam_total_lines;
  printf STDERR "II\t  Number of header lines:\t\t%u\n", $sam_header_lines;
  printf STDERR "II\t  Number of incomplete records:\t\t%u\n", $sam_invalid_lines;
  printf STDERR "II\t  Number of records stored:\t\t%u\n", $total_sam;
}

##  After knowing the size of the SAM file, make another pass over it, 
##    classifying the records.
for (my $i = 0; $i < $total_sam; $i++) {
  my $line = $sam_input[$i];
    
  my ($flag_tmp, $rname_tmp, $pos_tmp, $cigar_tmp) = split /\t/, $line;
  
  ##  Skip features that are not aligned to any chromosome
  if ($rname_tmp eq "*") {
    $classify_sam_input[$i] = $CLASSIFY_SAM_UNALIGNED;
  }
  else {
    $classify_sam_input[$i] = $CLASSIFY_SAM_UNKNOWN;
  }
}


########################################
##  Make a pass over the SAM and then the GTF files:
##    1)  Use the SAM file to determine which chromosomes to process, as well as the tenative lengths
##    2)  Update the lengths using the GTF file
########################################

for (my $i = 0; $i < $total_sam; $i++) {  
  my $line = $sam_input[$i];

  if ($classify_sam_input[$i] == $CLASSIFY_SAM_UNALIGNED) {
    next;
  }
  
  my ($flag_tmp, $rname_tmp, $pos_tmp, $cigar_tmp) = split /\t/, $line;

  if (!defined ($chrlist_hash_input{$rname_tmp})) {
    $chrlist_hash_input{$rname_tmp} = $pos_tmp;
    $chrlist_sam_acc{$rname_tmp} = 0;
    $chrlist_gtf_acc{$rname_tmp} = 0;
  }
  else {
    if ($pos_tmp > $chrlist_hash_input{$rname_tmp}) {
      $chrlist_hash_input{$rname_tmp} = $pos_tmp;
    }
  }
}

for (my $i = 0; $i < $gtf_input_pos; $i++) {  
  my $line = $gtf_input[$i];

  my ($seqname_tmp, $source_tmp, $feature_tmp, $start_tmp, $end_tmp, $score_tmp, $strand_tmp, $frame_tmp, $attribute_tmp) = split /\t/, $line;
  
  if (!defined ($chrlist_hash_input{$seqname_tmp})) {
    next;
  }
  else {
    if ($end_tmp > $chrlist_hash_input{$seqname_tmp}) {
      $chrlist_hash_input{$seqname_tmp} = $end_tmp;
    }
  }
}


########################################
##  Read in the chromosome list to "confirm" the lengths are correct
########################################

if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "II\tProcessing chromosome lengths (%s)...\n", $chrlist_arg;
}

open (my $chr_fp, "<", $chrlist_arg) or die "EE\tCould not open $chrlist_arg for input!\n";
while (<$chr_fp>) {
  my $line = $_;
  chomp $line;

  my ($chr_tmp, $size_tmp) = split /\t/, $line;
  
  if (!defined ($chrlist_hash_input{$chr_tmp})) {
    next;
  }
  else {
    if ($size_tmp >= $chrlist_hash_input{$chr_tmp}) {
      $chrlist_hash_input{$chr_tmp} = $size_tmp;
    }
    elsif ($size_tmp < $chrlist_hash_input{$chr_tmp}) {
      printf STDERR "WW\tThe chromosome length for %s (%u) is too short!  Needs to be more than %u!\n", $chr_tmp, $size_tmp, $chrlist_hash_input{$chr_tmp};
      exit (1);
    }
  }
}
close ($chr_fp);

if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "\n";
}


########################################
##  Process a chromosome at a time
########################################

my $num_chrs_processed = 0;
my $total_chrs_processed = keys %chrlist_hash_input;
if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "II\tTotal chromosomes to process:  %u\n", $total_chrs_processed;
}

foreach my $this_chr (sort (keys %chrlist_hash_input)) {
  ##  Increase the length of the chromosome (1-based, so add 1)
  my $chr_length = $chrlist_hash_input{$this_chr} + 1;

  ########################################
  ##  Initialize the accumulator and gtf arrays
  ########################################

  ##  Arrays that are equal to the size of the chromosome
  my @freq_chr_array;  ##  Accumulator of number of reads that align;
                       ##  size equal to (length of the chromosome) + 1
  my @exon_chr_array;  ##  Array of indexes to the GTF array so that each
                       ##    position points to an exon
                       ##  Size equal to (length of the chromosome) + 1
  my @gene_chr_array;  ##  Array of indexes to the GTF array so that each
                       ##    position points to a gene
                       ##  Size equal to (length of the chromosome) + 1

  ##  Initialize each position in the chromosome.
  ##    Position 0 is unused since genomes are numbered from position 1.
  for (my $k = 0; $k <= $chr_length; $k++) {
    $freq_chr_array[$k] = $UNINIT_FREQUENCY;
    $exon_chr_array[$k] = $UNINIT_BASE;
    $gene_chr_array[$k] = $UNINIT_GENE;
  }


  ########################################
  ##  Find the end points of each gene
  ########################################

  ##  Beginning and ending position of the gene, as well as gene name
  my %gene_begin;
  my %gene_end;
  my %gene_name;

  ##  Find the features that are at either end of a gene (i.e., not necessarily start/end codon)
  for (my $i = 0; $i < $gtf_input_pos; $i++) {
    my $line = $gtf_input[$i];

    my ($seqname_tmp, $source_tmp, $feature_tmp, $start_tmp, $end_tmp, $score_tmp, $strand_tmp, $frame_tmp, $attribute_tmp) = split /\t/, $line;

    ##  Skip features that are not part of this chromosome
    if ($seqname_tmp ne $this_chr) {
      next;
    }
    
    ##  Get the gene name
    my $gene_id_tmp = GetGeneName ($attribute_tmp);

    ##  Record the gene
    if (!defined ($gene_name{$gene_id_tmp})) {
      $gene_name{$gene_id_tmp} = $gene_id_tmp;
    }

    ##  Handle the beginning of the gene
    if (!defined ($gene_begin{$gene_id_tmp})) {
      $gene_begin{$gene_id_tmp} = $start_tmp;
    }
    if ($start_tmp < $gene_begin{$gene_id_tmp}) {
      $gene_begin{$gene_id_tmp} = $start_tmp;
    }
    
    ##  Handle the end of the gene
    if (!defined ($gene_end{$gene_id_tmp})) {
      $gene_end{$gene_id_tmp} = $end_tmp;
    }
    if ($end_tmp > $gene_end{$gene_id_tmp}) {
      $gene_end{$gene_id_tmp} = $end_tmp;
    }
  }
  

  ########################################
  ##  Remove overlapping genes
  ########################################

  my %overlaps;  ##  Keep track of overlapping gene names

  ##  $key is the name of the gene
  foreach my $key (sort (keys %gene_name)) {
    ##  Check that we have both a beginning and an end to the gene
    if (!defined ($gene_begin{$key})) {
      printf STDERR "EE\tBeginning of gene missing for:  %s\n", $key;
      exit (1);
    }
    if (!defined ($gene_end{$key})) {
      printf STDERR "EE\tEnd of gene missing for:  %s\n", $key;
      exit (1);
    }
    
#     printf STDERR "%s\t%s\t%s\n", $key, $gene_begin{$key}, $gene_end{$key};

    my $gene_begin = $gene_begin{$key};
    my $gene_end = $gene_end{$key};
  
    ##  Should never happen
    if ($gene_begin > $gene_end) {
      printf STDERR "EE\tStart coordinate greater than end coordinate!\n";
      exit (1);
    }
  
    ##  Mark the region of this gene
    for (my $k = $gene_begin; $k <= $gene_end; $k++) {
      if (length ($gene_chr_array[$k]) == 0) {
        ##  No gene occupying this position yet
        $gene_chr_array[$k] = $key;
      }
      elsif ($gene_chr_array[$k] eq $key) {
        ##  Do nothing, already marked
      }
      elsif ($gene_chr_array[$k] ne $key) {
        ##  Overlapping gene position -- keep track of both genes that contributed 
        ##  to the overlap
        $overlaps{$key} = $key;
        $overlaps{$gene_chr_array[$k]} = $gene_chr_array[$k];
      }
    }
  }

  my $num_genes_tmp = keys (%gene_name);
  $overlap_gtf = keys (%overlaps);
  $not_overlap_gtf = $num_genes_tmp - $overlap_gtf;

  if (($verbose_level_arg >= $VERBOSE_CHR) || ($verbose_level_arg >= $VERBOSE_ALL)) {
    printf $log_fp "II\tTotal number of genes:  %u\n", $num_genes_tmp;
    printf $log_fp "II\t  Overlapping:  %u\n", $overlap_gtf;
    printf $log_fp "II\t  Not overlapping:  %u\n", $not_overlap_gtf;
    
    if ($verbose_level_arg >= $VERBOSE_ALL) {
      printf $log_fp "II\tList of genes:\n";
  
      foreach my $key (sort (keys %overlaps)) {
        printf $log_fp "GG\t  %s\n", $key;
      }
    }
  }


  ########################################
  ##  Mark the genes on @exon_chr_array as all introns from end-to-end
  ########################################

  ##  $key is the name of the gene
  foreach my $key (sort (keys %gene_name)) {
    my $gene_begin = $gene_begin{$key};
    my $gene_end = $gene_end{$key};
  
    ##  Skip this gene if it is part of an overlapping gene
    if (defined ($overlaps{$key})) {
      next;
    }

    ##  Should never happen
    if ($gene_begin > $gene_end) {
      printf STDERR "EE\tStart coordinate greater than end coordinate!\n";
      exit (1);
    }
      
    ##  Mark the entire region of this gene as an intron, as a form of initialization
    for (my $k = $gene_begin; $k <= $gene_end; $k++) {
      $exon_chr_array[$k] = $INTRON_BASE;
    }
  }
  
  ########################################
  ##  Mark the exons on @exon_chr_array
  ########################################

  ##  Mark @exon_chr_array with the position of the GTF record in @gtf_pos
  for (my $i = 0; $i < $gtf_input_pos; $i++) {
    my $line = $gtf_input[$i];
  
    my ($seqname_tmp, $source_tmp, $feature_tmp, $start_tmp, $end_tmp, $score_tmp, $strand_tmp, $frame_tmp, $attribute_tmp) = split /\t/, $line;

    ##  Skip features that are not part of this chromosome
    if ($seqname_tmp ne $this_chr) {
      next;
    }
    
    ##  Get the gene name; skip this exon if it is part of an overlapping gene
    my $gene_id_tmp = GetGeneName ($attribute_tmp);
    if (defined ($overlaps{$gene_id_tmp})) {
      next;
    }

    ##  Check that start position is less than end position before proceeding.
    if ($start_tmp > $end_tmp) {
      printf STDERR "EE\tUnexpected error in GTF:  %s\n", $line;
      printf STDERR "EE\t  Start position greater than end position!\n";
      exit (1);
    }
    
    ##  Mark the exon area bounded by start_tmp and end_tmp
    for (my $k = $start_tmp; $k <= $end_tmp; $k++) {
      $exon_chr_array[$k] = $i;
    }
  }

  ##  Number of introns, exons, etc. positions in this chromosome
  my $intron_pos_count = 0;
  my $not_exon_pos_count = 0;
  my $exon_pos_count = 0;

  if ($verbose_level_arg >= $VERBOSE_CHR) {
    for (my $k = 0; $k <= $chr_length; $k++) {
      if ($exon_chr_array[$k] == $INTRON_BASE) {
        $intron_pos_count++;
      }
      elsif ($exon_chr_array[$k] == $UNINIT_BASE) {
        $not_exon_pos_count++;
      }
      else {
        $exon_pos_count++;
      }
    }

    printf $log_fp "II\tSize of chromosome:  %u\n", $chr_length;
    printf $log_fp "II\t  Intron positions:  %u (%.2f%%)\n", $intron_pos_count, $intron_pos_count / ($chr_length) * 100;
    printf $log_fp "II\t  Exon positions:  %u (%.2f%%)\n", $exon_pos_count, $exon_pos_count / ($chr_length) * 100;
    printf $log_fp "II\t  Non-intron/exon positions:  %u (%.2f%%)\n", $not_exon_pos_count, $not_exon_pos_count / ($chr_length) * 100;
  }


  ########################################
  ##  Check each stored alignment
  ########################################

  my $total_bases_chr = 0;
  for (my $i = 0; $i < $total_sam; $i++) {  
    my $line = $sam_input[$i];
    
    my ($flag_tmp, $rname_tmp, $pos_tmp, $cigar_tmp) = split /\t/, $line;

    ##  Skip alignments that are not part of this chromosome
    if ($rname_tmp ne $this_chr) {
      next;
    }

    my $contains_n = 0;
    if ($cigar_tmp eq "*") {
      $unavailable_cigar_sam++;
      $classify_sam_input[$i] = $CLASSIFY_SAM_UNKNOWN;
      next;
    }
    elsif ($cigar_tmp =~ /(\d+)N/) {
      $n_cigar_sam++;
      $classify_sam_input[$i] = $CLASSIFY_SAM_EXON;  ##  If there is an 'N', assume it is part of an exon
      $contains_n = 1;
    }
    else {
      $not_n_cigar_sam++;
    }
  
    ##  Flags for indicating if the SAM record overlaps with intron or exon or both
    my $is_intron = $BOOLEAN_FALSE;
    my $is_exon = $BOOLEAN_FALSE;
  
    my $k = $pos_tmp;
    my $end_k = $k;
    while (length ($cigar_tmp) != 0) {
      my $len = 0;
      my $instruction = "";
      if ($cigar_tmp =~ /^(\d+)(\D)(.*)$/) {
        $len = $1;
        $instruction = $2;
        $cigar_tmp = $3;
      }
      else {
        printf STDERR "EE\tCould not match the cigar string %s of %s!\n", $cigar_tmp, $line;
        exit (1);
      }

      if ($instruction eq "M") {
        $end_k = $k + $len;
        $total_bases_chr += $len;
      }
      elsif ($instruction eq "I") {
        ##  No movement on reference -- insertion into the reference
      }
      elsif ($instruction eq "D") {
        ##  Deletion from the reference; need to advance
        $end_k = $k + $len;  
        $total_bases_chr += $len;
      }
      elsif ($instruction eq "N") {
        ##  Skipped region on the reference, with no bases being aligned to this area
        $k = $k + $len;
        $end_k = $k;
      }
      elsif ($instruction eq "S") {
        ##  No movement on reference
      }
      elsif ($instruction eq "H") {
        ##  No movement on reference
      }
      elsif ($instruction eq "P") {
        ##  No movement on reference
      }
      elsif ($instruction eq "=") {
        $end_k = $k + $len;
        $total_bases_chr += $len;
      }
      elsif ($instruction eq "X") {
        $end_k = $k + $len;
        $total_bases_chr += $len;
      }
      else {
        printf STDERR "EE\tUnknown instruction %s in %s!\n", $instruction, $line;
        exit (1);
      }
    
      while ($k < $end_k) {
        if ($exon_chr_array[$k] == $INTRON_BASE) {
          $is_intron = $BOOLEAN_TRUE;
        }
        elsif ($exon_chr_array[$k] >= 0) {
          $is_exon = $BOOLEAN_TRUE;
        }
        
        if ($contains_n == 0) {
          $freq_chr_array[$k]++;
        }
        $k++;
      }
    }
    
    if (($is_intron == $BOOLEAN_TRUE) && ($is_exon == $BOOLEAN_TRUE)) {
      $classify_sam_input[$i] = $CLASSIFY_SAM_BOTH;
    }
    elsif (($is_intron == $BOOLEAN_TRUE) && ($is_exon == $BOOLEAN_FALSE)) {
      $classify_sam_input[$i] = $CLASSIFY_SAM_INTRON;
    }
    elsif (($is_intron == $BOOLEAN_FALSE) && ($is_exon == $BOOLEAN_TRUE)) {
      $classify_sam_input[$i] = $CLASSIFY_SAM_EXON;
    }
    else {
      $classify_sam_input[$i] = $CLASSIFY_SAM_UNKNOWN;
    }
  }

  $total_bases += $total_bases_chr;
  if ($verbose_level_arg >= $VERBOSE_CHR) {
    printf $log_fp "II\tTotal bases:  %u\n", $total_bases_chr;
    printf $log_fp "II\tSAM entries:\n";
    printf $log_fp "II\t  Unavailable CIGAR string (skipped):  %u\n", $unavailable_cigar_sam;
    printf $log_fp "II\t  With N (skipped):  %u\n", $n_cigar_sam;
    printf $log_fp "II\t  Without N (used):  %u\n", $not_n_cigar_sam;
  }


  ########################################
  ##  Output locations of intron retention, scanning from beginning of chromosome
  ########################################

  ##  At the end of the while () loop:
  ##
  ##  EEEEEEEEIIIIIIIIEEEEEEEE
  ##         i        j
  ##
  ##    E = exon
  ##    I = intron
  ##

  my $intron_retention = 0;  ##  Number of intron retentions
  my $not_intron_retention = 0;  ##  Number of non-intron retentions
  my %unique_genes;

  for (my $i = 0; $i < $chr_length; $i++) {
    my $j = $i + 1;
    if ($j == $chr_length) {
      last;
    }
  
    ##  $i continually moves until $i and $j mark the border between the end of an 
    ##  exon and the start of a non-exon region
    if (($exon_chr_array[$i] >= 0) && ($exon_chr_array[$j] == $INTRON_BASE)) {
      ##  Need to advance $j to the next exon
      while ($exon_chr_array[$j] == $INTRON_BASE) {
        $j++;
        if ($j == $chr_length) {
          ##  Set an exit condition
          $j = $i;
          last;
        }
      }

      ##  We've reached the end of the chromosome
      if ($i == $j) {
        last;
      }
    
      ##  Take the two exon records at the opposite ends of the intron
      my $exon_i = $gtf_input[$exon_chr_array[$i]];
      my $exon_j = $gtf_input[$exon_chr_array[$j]];
    
      ##  Make sure both are from the same gene, since that is the definition of an intron
      my ($seqname_i, $source_i, $feature_i, $start_i, $end_i, $score_i, $strand_i, $frame_i, $attribute_i) = split /\t/, $exon_i;
      my ($seqname_j, $source_j, $feature_j, $start_j, $end_j, $score_j, $strand_j, $frame_j, $attribute_j) = split /\t/, $exon_j;
    
      my $gene_id_i = GetGeneName ($attribute_i);
      my $gene_id_j = GetGeneName ($attribute_j);

      if ($gene_id_i ne $gene_id_j) {
        printf STDERR "EE\tThis should never happen since we're guaranteed to be searching between two exons of the same gene!\n";
        printf STDERR "EE\t  Genes:  %s and %s\n", $gene_id_i, $gene_id_j;
        exit (1);
      }
    
      if ($strand_i ne $strand_j) {
        printf STDERR "WW\tThis should never happen since we're guaranteed to be searching between two exons of the same gene!\n";
        printf STDERR "WW\t  Genes:  %s and %s\n", $gene_id_i, $gene_id_j;
        printf STDERR "WW\t  Strand:  %s and %s\n", $strand_i, $strand_j;
        exit (1);
      }
      
      my $intron_width = $j - $i - 1;
      my $intron_count = 0;
      for (my $k = $i + 1; $k < $j; $k++) {
        $intron_count += $freq_chr_array[$k];
      }

      if ($intron_count == 0) {
        ##  Not a case of intron retention
        $not_intron_retention++;
      }
    
      if (!defined ($unique_genes{$gene_id_i})) {
        $unique_genes{$gene_id_i} = $gene_id_i;
      }

      my $out_start = $i + 1;
      my $out_end = $j - 1;
      my $out_attr = sprintf ("name=%s;count=%u;width=%u", $gene_id_i, $intron_count, $intron_width);
      my $outline = join ("\t", $this_chr, "sam2irs", "intron", $out_start, $out_end, $intron_count, $strand_i, "0", $out_attr);
      
      $introns_output[$num_introns_output] = $outline;
      $num_introns_output++;
      
      $intron_retention++;
    }
  }

  
  ########################################
  ##  Perform a sanity check across the chromosome
  ##    * There shouldn't be two different genes next to each other
  ########################################

  if ($exons_arg) {
    for (my $i = 0; $i < $chr_length; $i++) {
      my $j = $i + 1;
      if ($j == $chr_length) {
        last;
      }
  
      if (($exon_chr_array[$i] >= 0) && ($exon_chr_array[$j] >= 0)) {
        ##  Take both exon records
        my $exon_i = $gtf_input[$exon_chr_array[$i]];
        my $exon_j = $gtf_input[$exon_chr_array[$j]];
    
        my ($seqname_i, $source_i, $feature_i, $start_i, $end_i, $score_i, $strand_i, $frame_i, $attribute_i) = split /\t/, $exon_i;
        my ($seqname_j, $source_j, $feature_j, $start_j, $end_j, $score_j, $strand_j, $frame_j, $attribute_j) = split /\t/, $exon_j;
    
        my $gene_id_i = GetGeneName ($attribute_i);
        my $gene_id_j = GetGeneName ($attribute_j);
      
        ##  In between two different genes.  
        ##    Not an intron, so go to next loop iteration
        if ($gene_id_i ne $gene_id_j) {
          printf STDERR "WW\t* Two exons of two different genes are connected to each other!\n";
          printf STDERR "WW\t  Positions:  %u vs %u\n", $i, $j;
          printf STDERR "WW\t  [%s] vs [%s]\n", $gene_id_i, $gene_id_j;
          printf STDERR "WW\t  [%s]\n", $exon_i;
          printf STDERR "WW\t  [%s]\n\n", $exon_j;
          exit (1);
        }
      }
    }
  }

  
  ########################################
  ##  Output locations of exons, if requested, scanning from beginning of chromosome
  ########################################

  ##  At the end of the while () loop:
  ##
  ##  XXXXXXXXEEEEEEEEXXXXXXXX
  ##         i        j
  ##
  ##    E = exon
  ##    X = non-exon (i.e., intron, etc.); value of UNINIT_BASE or INTRON_BASE
  ##

  ##  Change <= to < so that we stay within the array; perform a single pass over $exon_chr_array
  if ($exons_arg) {
    for (my $i = 0; $i < $chr_length; $i++) {
      my $j = $i + 1;
      if ($j == $chr_length) {
        last;
      }
  
      ##  $i continually moves until $i and $j mark the border at the beginning of an exon
      if (($exon_chr_array[$i] < 0) && ($exon_chr_array[$j] >= 0)) {
        my $exon_i = $gtf_input[$exon_chr_array[$i + 1]];
        my ($seqname_i, $source_i, $feature_i, $start_i, $end_i, $score_i, $strand_i, $frame_i, $attribute_i) = split /\t/, $exon_i;
        my $gene_id_i = GetGeneName ($attribute_i);

        my $exon_j = "";
        my $gene_id_j = "";
      
        ##  Need to advance $j to the next non-exon
        while ($exon_chr_array[$j] >= 0) {
          $j++;
        
          ##  Unfortunately, we also have to check that we're on the same gene since two genes
          ##  may be adjacent to each other
          $exon_j = $gtf_input[$exon_chr_array[$j - 1]];
    
          my ($seqname_j, $source_j, $feature_j, $start_j, $end_j, $score_j, $strand_j, $frame_j, $attribute_j) = split /\t/, $exon_j;
    
          $gene_id_j = GetGeneName ($attribute_j);

          ##  Two different genes next to each other, so leave this loop
          if ($gene_id_i ne $gene_id_j) {
            last;
          }
        
          if ($j == $chr_length) {
            ##  Set an exit condition
            $j = $i;
            last;
          }
        }

        ##  We've reached the end of the chromosome
        if ($i == $j) {
          last;
        }

        ##  Take the ends of the exons, and check the gene
        my ($seqname_j, $source_j, $feature_j, $start_j, $end_j, $score_j, $strand_j, $frame_j, $attribute_j) = split /\t/, $exon_j;
      
        ##  In between two different genes.  
        ##    Not an intron, so go to next loop iteration
        if ($gene_id_i ne $gene_id_j) {
          printf STDERR "WW\tTwo exons of two different genes are connected to each other!\n";
          printf STDERR "WW\t  Positions:  %u vs %u\n", $i, $j;
          printf STDERR "WW\t  [%s] vs [%s]\n", $gene_id_i, $gene_id_j;
          printf STDERR "WW\t  [%s]\n", $exon_i;
          printf STDERR "WW\t  [%s]\n\n", $exon_j;
#           for (my $k = $i; $k < $i; $k++) {
#             my ($seqname_k, $source_k, $feature_k, $start_k, $end_k, $score_k, $strand_k, $frame_k, $attribute_k) = split /\t/, $gtf_input[$exon_chr_array[$k]];
#             my $gene_id_k = GetGeneName ($attribute_k);
#             printf STDERR "EE\t%u\t%s\n", $k, $gene_id_k;
#           }
          exit (1);
        }
    
        if ($strand_i ne $strand_j) {
          printf STDERR "WW\tTwo exons of the same gene are on opposite strands!\n";
          printf STDERR "WW\t  [%s] vs [%s]\n", $gene_id_i, $gene_id_j;
          printf STDERR "WW\t  Strand:  %s and %s\n", $strand_i, $strand_j;
          exit (1);
        }
      
        my $exon_width = $j - $i - 1;
        my $out_start = $i + 1;
        my $out_end = $j - 1;
        my $out_attr = sprintf ("name=%s;count=0;width=%u", $gene_id_i, $exon_width);
        my $outline = join ("\t", $this_chr, "sam2irs", "exon", $out_start, $out_end, "0", $strand_i, "0", $out_attr);
      
        $exons_output[$num_exons_output] = $outline;
        $num_exons_output++;
      }
    }
  }
    
  my $intron_retention_genes = keys (%unique_genes);

  if ($verbose_level_arg >= $VERBOSE_CHR) {
    printf $log_fp "II\tIntron retention genes:  %u\n", $intron_retention_genes;
    printf $log_fp "II\tIntron retentions:  %u\n", $intron_retention;
    printf $log_fp "II\tNon-intron retentions:  %u\n", $not_intron_retention;
    printf $log_fp "\n";
  }

  $num_chrs_processed++;

  if ($num_chrs_processed % 10 == 0) {
    if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
      printf STDERR "II\t%u of %u chromosomes processed...\n", $num_chrs_processed, $total_chrs_processed;
    }
  }
}  ##  End outer loop


########################################
##  Classify each SAM record
########################################

my $classify_sam_unaligned_count = 0;
my $classify_sam_unknown_count = 0;
my $classify_sam_intron_count = 0;
my $classify_sam_exon_count = 0;
my $classify_sam_both_count = 0;
for (my $i = 0; $i < $total_sam; $i++) {
  if ($classify_sam_input[$i] == $CLASSIFY_SAM_INTRON) {
    $classify_sam_intron_count++;
  }
  elsif ($classify_sam_input[$i] == $CLASSIFY_SAM_EXON) {
    $classify_sam_exon_count++;
  }
  elsif ($classify_sam_input[$i] == $CLASSIFY_SAM_BOTH) {
    $classify_sam_both_count++;
  }
  elsif ($classify_sam_input[$i] == $CLASSIFY_SAM_UNKNOWN) {
    $classify_sam_unknown_count++;
  }
  elsif ($classify_sam_input[$i] == $CLASSIFY_SAM_UNALIGNED) {
    $classify_sam_unaligned_count++;
  }
}

if (defined ($config -> get ("samrecord"))) {
  open (my $sam_record_fp, ">", $sam_record_arg) or die "EE\tCould not create $sam_record_arg!\n";
  printf $sam_record_fp "%u\n", $total_sam;
  printf $sam_record_fp "%u\n", $classify_sam_intron_count;
  printf $sam_record_fp "%u\n", $classify_sam_exon_count;
  printf $sam_record_fp "%u\n", $classify_sam_both_count;
  printf $sam_record_fp "%u\n", $classify_sam_unaligned_count;
  printf $sam_record_fp "%u\n", $classify_sam_unknown_count;
  for (my $i = 0; $i < $total_sam; $i++) {
    printf $sam_record_fp "%u", $classify_sam_input[$i];
  }
  printf $sam_record_fp "\n";
  close ($sam_record_fp);
}

my $aligned = $total_sam - $classify_sam_unaligned_count;

for (my $i = 0; $i < $num_introns_output; $i++) {
#   printf "%s\n", $introns_output[$i];
  printf "%s;aligned_reads=%u;aligned_bases=%u;\n", $introns_output[$i], $aligned, $total_bases;
}

if ($exons_arg) {
  for (my $i = 0; $i < $num_exons_output; $i++) {
    printf "%s;aligned_reads=%u;aligned_bases=%u;\n", $exons_output[$i], $aligned, $total_bases;
  }
}

if ($verbose_level_arg >= $VERBOSE_SUMMARY) {
  printf STDERR "\n";
  printf STDERR "II\tTotal SAM records processed:  %u\n", $total_sam;
  printf STDERR "II\t  Intron records:  %u\n", $classify_sam_intron_count;
  printf STDERR "II\t  Exon records:  %u\n", $classify_sam_exon_count;
  printf STDERR "II\t  Intron+Exon records:  %u\n", $classify_sam_both_count;
  printf STDERR "II\t  Unaligned records:  %u\n", $classify_sam_unaligned_count;
  printf STDERR "II\t  Unknown records:  %u\n", $classify_sam_unknown_count;
  printf STDERR "\n";
}


########################################
##  Close log file
########################################

if (length ($log_fn_arg) != 0) {
  close ($log_fp);
}



=pod

=head1 NAME

sam2irs.pl -- Converts a set of mapped reads in SAM format to a table of intron retention scores (regions between exons with aligned reads).

Input format is SAM; output format is GFF3.  Use perldoc to properly format this documentation.

=head1 SYNOPSIS

B<sam2irs.pl> <input.sam >output.txt

More detailed help:

  perldoc sam2irs.pl

=head1 DESCRIPTION

This script scores the amount of aligned bases within introns by reading in a set of aligned reads.

=head1 ALGORITHM

This script executes the following steps:

=over 5

=item 1. Initialization of variables.

=item 2. Process the command-line arguments.

=item 3. Read in the gene feature annotations (i.e., the GTF file provided with the I<--gtf> argument).

=over 2

=item a. Ignore the annotations for chromosomes (or scaffolds) whose chromosome lengths were not provided.

=item b. Ignore predicted genes (i.e., those of the form "Gm\d+").

=item c. Store gene features that are "exons".  For all other gene features, just count how many there are.

=back

=item 4. Read in the SAM file (i.e., provided via standard in) and store all of the alignments in an array.

=item 5. Make a pass over the stored SAM file and indicate whether the read is:

=over 2

=item a. Unaligned or

=item b. Unknown (i.e., will determine later)

=back

=item 6. Make a pass over the SAM file to determine which chromosomes are to be processed, and how long each one is.

=item 7. Use the GTF file to "correct" these chromosome lengths since some genes may exist after where reads align to.

=item 8. For each chromosome, do the following:

=over 2

=item a. Initialize the chromosome arrays.

=item b. For each gene, find its beginning and end.  This would be the exons with the lowest and highest chromosome positions.

=item c. For each gene, mark its location on the chromosome array B<gene_chr_array>.  Using this array, if two genes overlap, then mark both of them as overlapping.

=item d. For each gene, initialize the entire gene's region as an intron.

=item e. For each exon, do the following:

=over 2

=item i. Skip features that are not associated with the chromosome being processed.

=item ii. Get the gene associated with this exon and skip it if it is part of an overlapping gene.

=item iii. Mark area from the start to the end position of the exon.

=item iv. Because of the initialization step previously, any unmarked regions are introns.

=back

=back

=item f. For each alignment, do the following:

=over 2

=item i. Skip alignments not associated with the chromosome being processed.

=item ii. Process each instruction in the CIGAR string, incrementing the number of aligned bases a position at a time.

=item iii. Record whether the alignment is on both an intronic region, an exonic region, or both.

=back

=item g. Scan the entire chromosome for locations where introns occur and output each intron's score.

=item h. If the I<--samrecord> option was provided, then also output a record of where each alignment mapped to (i.e., an intron, an exon, both, or somewhere else).

=back

=head1 OPTIONS

=over 5

=item --chrlist I<file>

A file indicating the list of chromosomes and their lengths.

=item --gtf I<file>

Indicate the GTF file of genomic features to use.  Only exons will be used; all other features will be ignored.
  
=item --samrecord I<file>

Create a record of what type of region each read in the SAM file aligns to.
  
=item --exons

By default, intron regions are output.  This flag indicates that exons should be output as well.
  
=item --verbose I<level>

Display verbose debugging information about the execution of this script based on the "level" provided to standard error.  Note that a lot of output is provided if levels 2 or 3 are selected.  Because of this, the --log option is also required.

=over 2

=item 0.  Do not output any debugging information.

=item 1.  Provide summary information as the script starts and as it finishes to STDERR.

=item 2.  Provide information as each chromosome is processed.

=item 3.  Provide all debugging information.

=back

=item --log I<file>

Filename of the log file to store verbose information at levels 2 and above.
  
=item --help

Display this help message.

=back

=head1 EXAMPLE

Go into the Examples/ subdirectory and type the following:

cat test1.sam | ../Perl/sam2irs.pl --verbose --chrlist test1.genome --gtf test1.gtf 2>/dev/null`

=head1 NOTES

1. Co-ordinates output from this script are B<1-based>.

2. The co-ordinates of each intron are inclusive of the end-points.  Note that the original co-ordinates of the gene are not retained in the output, even if introns are "merged".  This is because only the beginning of the first intron and the end of the last intron contained within the output file.

=head1 AUTHOR

Raymond Wan <rwan.work@gmail.com>

=head1 COPYRIGHT

Copyright (C) 2016-2020, Raymond Wan, All rights reserved.


