sam2irs
=======


Introduction
------------

This repository, **sam2irs**, consists of a set of programs for calculating intron retention scores (IRS') using aligned sequences in the sequence alignment/map (SAM) format as input.  While the main script is called sam2irs, additional scripts have been included to supplement it.  An accompanying script called irs2fasta takes the intronic regions and extracts the up and downstream sequences.

Besides this document, this repository contains the following:

  * source code in Perl,
  * `environment.yml` file for creating a [conda](https://docs.conda.io/en/latest/) environment,
  * licensing information, and
  * a small data file for testing.

The manuscript that makes use of this software has been published as:

    L. Yue, R. Wan, S. Luan, W. Zeng, T. H. Cheung.  Dek Modulates Global Intron Retention during Muscle Stem Cells Quiescence Exit, Developmental Cell, 53(6), pg. 661--676, 2020.

It is available at this [link](https://doi.org/10.1016/j.devcel.2020.05.006).  An earlier version (which should not be cited) appears in Cell Press' Sneak Peak [service](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3385139).


Requirements
------------

| Software                | Version | Required? | Web site                                    |
|:-----------------------:|:-------:|:---------:|:-------------------------------------------:|
| Perl                    | 5.22.0  | Yes       | https://www.perl.org                        |
| AppConfig (Perl module) | 1.71    | Yes       | https://metacpan.org/pod/AppConfig          |
| conda                   | 4.7.11  | Yes       | https://docs.conda.io/en/latest/            |
| Dreme (Meme suite)      | 5.0.5   | Yes       | http://meme-suite.org/doc/dreme.html        |

Experiments with this software have been successfully run on Linux systems running:
  * CentOS 6.9 64-bit and
  * Ubuntu 18.04 and 19.04

The versions above represent the tools used during software development or when running the experiments in the paper. They do not represent the minimum requirements; it is possible that lower versions can be used.

We highly recommend using `conda`.  If you are, then you can create an environment using `conda env create -f environment.yml`, as explained in the online conda [instructions](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).  This will install the required programs above (with the exception of `conda` itself, of course).


Organization
------------

After cloning this repository from GitHub using the `git clone` command, the following directory structure is obtained:

    .
    ├── environment.yml
    ├── Examples                     Examples to test sam2irs against
    ├── LICENSE                      Software license (GNU GPL v3)
    ├── Perl                         Perl scripts
    │   ├── generate-negative.pl     Generate a set of negative sequences
    │   ├── intersect-irs.pl         Intersect two sets of IRS files
    │   ├── irs2fasta.pl             Extract FASTA sequences upstream / downstream of introns
    │   ├── sam2irs.pl               Main sam2irs script
    │   └── score-irs.pl             Script for calculating IRS and taking the top N IRS'
    └── README.md                    This README file


Input/Output Format
-------------------

### Input

As input, sam2irs requires:

1.  A list of chromosomes;
2.  A list of the exons of each gene (in the gene transfer format, [GTF](https://asia.ensembl.org/info/website/upload/gff.html)); and
3.  The aligned reads (in the Sequence Alignment/Map format, [SAM](http://www.htslib.org/doc/sam.html)).

The list of chromosomes is a tab-separated file, with each row corresponding to a chromosome.  Each row has two fields.  The first is the name of the chromosome and the second is the length of the chromosome in base pairs.  Only chromosomes in this list are processed; chromosomes mentioned in the SAM file that are **not** in this list are therefore skipped.

The input format for `irs2fasta` is the output format of `sam2irs`.


### Output

The output is in [GTF](https://asia.ensembl.org/info/website/upload/gff.html).  The second field is the "source" -- that is, the program that performed the calculation.  In our case, it is always "sam2irs".  The third field is the feature type, which is always "intron" (a non-standard feature type which we have used).

The ninth and last field is a set of attributes separated by semi-colons.  In these attributes, we list the following information:

  * name -- name of the gene associated with this intron
  * count -- number of bases aligned to this intron
  * width -- width of the intron
  * aligned_reads -- number of reads that were aligned (used for normalisation)
  * aligned_bases -- number of bases that were aligned (used for normalisation)

The number of aligned reads and bases is the same for every feature.  While inefficient in terms of disk space, this allows each row to be normalised independently.

The output format for irs2fasta is the FASTA format.  That is, it is a text file of records.  Each record consists of two lines.  The first line is the name of the sequence preceded by a ">" character.  The second line is a DNA sequence.


Method
------

The sam2irs script proceeds as follows:

  1.  Read in the list of exons:
    * Exclude genes from chromosomes that are are not to be considered (i.e., usually unassigned scaffolds, etc.).
    * Exclude predicted genes that start with "Gm".
    * Exclude all features (as indicated by the third field of the GTF format) that are not exons.  More specifically, the features "CDS", "start_codon", and "stop_codon" are dropped.
  2.  Read in the aligned reads.
  3.  Using the aligned reads and the list of exons, determine the "length" of each chromosome.
  4.  Process each chromosome of the genome one-by-one (see below).
  5.  Output the intron retention scores.
  6.  If the --samrecord argument is provided, then print to file the "classification" of each read.

In essence, the script creates multiple arrays, with each array equal in size to the length of the chromosome under consideration.  These arrays might indicate the gene associated with a particular location (thus, if more than one gene is associated with a single location, then it is an overlapping gene); whether the position is an intron or an exon; and the number of bases aligned to that position.  In the Perl version of `sam2irs`, the names of these 3 arrays are `gene_chr_array`, `exon_chr_array`, and `freq_chr_array`, respectively.

The main part of the script is step #4, which processes each chromosome one-by-one.  The main steps this part of the script is as follows:

  1.  Record the starting and ending position of each gene (since the GTF file does not have to be sorted by location) in a set of dictionaries.
  2.  For each gene, mark the area between its starting and ending position.  Overlapping areas imply overlapping genes.  If this happen, the genes in question (i.e., potentially 2 or more) are discarded.
  3.  For each remaining gene, mark the entire area between the starting and ending position as an "intron" as a form of initialisation.  Afterwards, go back and indicate the location of the exons.  Thus, this "collapses" all of the exons of a gene's isoforms so that intronic regions are those regions which are introns for all isoforms of the gene.
  4.  For each aligned read, use the CIGAR string to indicate how many bases aligned at each position of the chromosome.  Also record whether the read aligns to an exon, an intron, both intron and exon (i.e., it straddles the boundary between the two), or neither.


Additional documentation
------------------------

In addition to this `README.md`, documentation has been placed at the end of the script in Perl's [PerlPod](https://perldoc.perl.org/perlpod) format.  To view it, you can type `perldoc sam2irs.pl`.

Note that the PerlDoc software has to be installed in order to view it.  Alternatively, you can just go to the end of the script with your favourite editor.


Running example
---------------

###  Calculating the IRS

In the `Examples/` directory there is a simple example that depicts how IRS is calculated.  This example consists of one chromosome, one gene, and 6 aligned reads of 6 base pairs in length each.  The gene has two exons.  Bases that lie within an intron are in blue; all other bases are in red.  This is illustrated in the following image (the first position is 1):

![Examples/test1.svg](Examples/test1.svg)

Type the following to process this example:

  * `cat test.sam | ../Perl/sam2irs.pl --verbose 0 --chrlist test.sizes --gtf test1.gtf 1>output.gtf 2>error.txt`

Any errors will be sent to the file `error.txt`.  The output will be sent to the file `output.gtf`, which in our example is the following single line in GTF format:

chr1    sam2irs        intron  11      17      15      -       0       name=Test;count=15;width=7;aligned_reads=6;aligned_bases=36;

This output indicates that the chromosome "chr1" has an intron from position 11 to 17 (where the first base is position 1 and endpoints are included in the interval), which has 15 intronic bases.  The width of this region is 7 base pairs.  The total number of aligned bases is 36 for the entire data set.

The intron retention score is:  (15 / 7) / 36 = 0.0595


###  Obtain the top-N IRS regions

The above toy example only has only one intronic region.  In a real data set, you might have many regions and you want to have the top **N** sites.  To do that, take the `output.gtf` above and use the `Perl/score-irs.pl` as follows:

  * `cat output.gtf | ../Perl/score-irs.pl --topn 100 >output-topN.gtf`

By default, only the top N regions are printed.  The normalised score (i.e., 0.0595 in the above example) it not output.  If you want them output, then use the `--showscore` option.  The score will be shown as the first column.

In practice, the number of aligned bases will be very large and multiplying by a constant across all data sets is needed to provide "sensible" values.  In our pulished work, we arbitrarily chose a constant of `100000000000`.  You will need to find a constant that gives you the range of values that fits for your data set, and then divide it through for your entire file.  The value of this constant does not matter; what matters is that you use the same constant for all of your samples.  You can modify this constant in the code, or you can provide a new value using the `--normaliser` option.

For example, combining these two options, you can do the following:

  * `cat output.gtf | ../Perl/score-irs.pl --topn 100 --showscore --normaliser 100000 >output-topN.gtf`
  
Note that without the `--showscore` option, the normaliser constant is not too important since the regions are printed out based on relative values.


Problems
--------

Hopefully, you do not encounter problems while using this program.  However, if you do, then you can obtain additional information by increasing the level of verbosity (i.e., the `--verbose` option).  Instead of the `0` in the above example, consider other values:

0.  Do not output any debugging information.
1.  Provide summary information as the script starts and as it finishes to standard error.
2.  Provide information as each chromosome is processed to a log file (--log).
3.  Provide all debugging information to a log file (--log).

Note that fatal errors will cause the program to terminate immediately with an error message.  This will happen regardless of the verbose level provided.


Hints
-----

1.  `sam2irs` does not directly support compressed SAM files in BAM format.  However, since it accepts the SAM file via standard in, you could simply do this instead:

    `samtools view test.bam | ../Perl/sam2irs.pl --verbose 0 --gtf test1.gtf 1>output.gtf 2>error.txt`

2.  The amount of memory used by `sam2irs` is dominated by the length of the longest chromosome being processed, the size of the SAM file, and the size of the GTF file.  If the length of the longest chromosome is `n`, then the amount of memory used is `3 * 4n`.  So, if the longest chromosome is 1000 base pairs, then the amount of memory used by this program is *approximately* 12000 bytes plus the sizes of the SAM and GTF files.

3.  While each chromosome is processed independently, the entire SAM file is read in and kept in memory.  If you are having problems processing your data, then you can separate the SAM file so that each chromosome is a separate file (i.e., by using `grep`).  Then, process each chromsome one-by-one.


About sam2irs
-------------

This software was implemented while I was at the Hong Kong University of Science and Technology.  My contact details:

     E-mail:  rwan.work@gmail.com 

My homepage is [here](http://www.rwanwork.info/).

The latest version of sam2irs can be downloaded from [GitHub](https://github.com/rwanwork/sam2irs).

If you have any information about bugs, suggestions for the documentation or just have some general comments, feel free to contact me via e-mail or as a GitHub issue.  (Of the two, I prefer GitHub.)


Copyright and License
---------------------

     sam2irs (SAM to intron retention score calculator)
     Copyright (C) 2016-2020 by Raymond Wan

sam2irs is distributed under the terms of the GNU General Public License (GPL, version 3 or later) -- see the file LICENSE for details.

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts and no Back-Cover Texts. A copy of the license is included with the archive as LICENSE.


