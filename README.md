sam2irs
=======


Introduction
------------

sam2irs is a program for calculating intron retention scores (IRS) using input in
the sequence alignment/map (SAM) format.  It was developed using Perl but a C++
version is in development.

This document accompanies:
  * source code in Perl,
  * licensing information, and
  * a small data file for testing.
  
The manuscript that makes use of this software is currently in preparation.


Requirements
------------

|Software                |Version  |Required?  |Ubuntu package   |Anaconda package |
|:----------------------:|:-------:|:---------:|:---------------:|:---------------:|
|Perl                    | 5.22.0  | Yes       |perl             |perl             |
|AppConfig (Perl module) | 1.71    | Yes       |libappconfig-perl|perl-appconfig   |

Experiments with this software have been successfully run on Linux systems running:
  * CentOS 6.9 64-bit and
  * Ubuntu 18.04

The versions represent the tools used during software development or when running the experiments in the paper. They do not represent the minimum requirements; it is possible that lower versions can be used.


Organization
------------

After cloning this repository from GitHub, the following directory structure is obtained:

    ├── Cpp                     C++ version of sam2irs (in development).
    ├── Examples                Examples to test sam2irs against.
    │   ├── test1.*             Files related to a very simple test case.
    ├── LICENSE                 Software license (GNU GPL v3).
    ├── Perl                    Perl version of sam2irs.
    └── README.md               This README file.


Running examples
----------------

### test1

In the Examples/ directory type:

  * `cat test1.sam | ../Perl/sam2irs.pl --verbose --debug --chrlist test1.genome --gtf test1.gtf 2>/dev/null`

With standard error sent to /dev/null, the output will be a single line in GTF format:

chr1    sam2irs        intron  11      17      15      -       0       name=Test;count=15;width=7;aligned_reads=6;aligned_bases=36;

This output indicates that the chromosome "chr1" has an intron at from position 11 to 17 (endpoints are 
included in the interval), which has 15 intronic bases.  The width of this region is 7 base pairs.
The total number of aligned bases is 36 for the entire data set.  

The intron retention score is:  (15 / 7) / 36.


Future Work
-----------

A hopefully faster implementation in C++ is currently being developed.


About sam2irs
-------------

This software was implemented by Raymond Wan as part of my employment at the 
Hong Kong University of Science and Technology.  My contact details:

     E-mails:  rwan.work@gmail.com 
               OR 
               raymondwan@ust.hk

My homepage is [here](http://www.rwanwork.info/).

The latest version of sam2irs can be downloaded from [GitHub](https://github.com/rwanwork/sam2irs).

If you have any information about bugs, suggestions for the documentation or just have some general 
comments, feel free to write to one of the above addresses.


Copyright and License
---------------------

     sam2irs (SAM to intron retention score calculator)
     Copyright (C) 2016-2018 by Raymond Wan

sam2irs is distributed under the terms of the GNU General
Public License (GPL, version 3 or later) -- see the file COPYING for details.

Permission is granted to copy, distribute and/or modify this document under the
terms of the GNU Free Documentation License, Version 1.3 or any later version
published by the Free Software Foundation; with no Invariant Sections, no
Front-Cover Texts and no Back-Cover Texts. A copy of the license is included
with the archive as COPYING.DOC.


Wednesday, June 27, 2018


