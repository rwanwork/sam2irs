#####################################################################
##  other-rules.smk
##
##  Raymond Wan (raymondwan@cuhk.edu.hk)
##
##  Department of Medicine and Therapeutics, 
##  Division of Endocrinology & Diabetes,
##  The Chinese University of Hong Kong, Hong Kong
##
##  Copyright (C) 2019, Raymond Wan, All rights reserved.
#####################################################################

rule Score_IRS:
  input:
    input_fn1 = INPUT_DIR + "/{sample}/{sample}-1.txt",
    input_fn2 = INPUT_DIR + "/{sample}/{sample}-2.txt"
  output:
    output_fn1 = OUTPUT_DIR + "/scores/{sample}/{topn}/score-1.txt",
    output_fn2 = OUTPUT_DIR + "/scores/{sample}/{topn}/score-2.txt",
    output_fn1 = OUTPUT_DIR + "/scores/{sample}/{topn}/score-1-IRS.txt",
    output_fn2 = OUTPUT_DIR + "/scores/{sample}/{topn}/score-2-IRS.txt"
  params:
    topn = "{topn}"
  shell:
    """
    cat {input.input_fn1} | Perl/score-irs.pl --topn {params.topn} >{output.output_fn1}
    cat {input.input_fn2} | Perl/score-irs.pl --topn {params.topn} >{output.output_fn2}

    cat {input.input_fn1} | Perl/score-irs.pl --showscore --topn {params.topn} >{output.output_fn1}
    cat {input.input_fn2} | Perl/score-irs.pl --showscore --topn {params.topn} >{output.output_fn2}
    """

rule Calculate_Intersect:
  input:
    input_fn1 = OUTPUT_DIR + "/scores/{sample}/{topn}/score-1.txt",
    input_fn2 = OUTPUT_DIR + "/scores/{sample}/{topn}/score-2.txt"
  output:
    output_fn1 = OUTPUT_DIR + "/intersect/{sample}/{topn}/intersect.txt"
  shell:
    """
    cat {input.input_fn1} | Perl/intersect-irs.pl --first {input.input_fn2} >{output.output_fn1}
    """

rule Get_Positive_Seqs:
  input:
    input_fn1 = OUTPUT_DIR + "/intersect/{sample}/{topn}/intersect.txt",
    input_fn2 = "genome-masked.fna"
  output:
    output_fn1 = OUTPUT_DIR + "/sequences/{sample}/{topn}/{winsize}/{method}/positive.fna"
  params:
    winsize = "{winsize}",
    method = "{method}"
  shell:
    """
    cat {input.input_fn1} | Perl/irs2fasta.pl --verbose --genomefn {input.input_fn2} --region {params.method} >{output.output_fn1}
    """


rule Get_Negative_Seqs:
  input:
    input_fn1 = OUTPUT_DIR + "/sequences/{sample}/{topn}/{winsize}/{method}/positive.fna"
  output:
    output_fn1 = OUTPUT_DIR + "/sequences/{sample}/{topn}/{winsize}/{method}/negative.fna"
  params:
    winsize = "{winsize}"
  shell:
    """
    cat {input.input_fn1} | Perl/generate-negative.pl >{output.output_fn1}
    """

rule Run_Dreme:
  input:
    input_fn1 = OUTPUT_DIR + "/sequences/{sample}/{topn}/{winsize}/{method}/positive.fna",
    input_fn2 = OUTPUT_DIR + "/sequences/{sample}/{topn}/{winsize}/{method}/negative.fna"
  output:
    output_fn1 = OUTPUT_DIR + "/motif/{sample}/{topn}/{winsize}/{method}/dreme.html"
  log:
    log_fn1 = OUTPUT_DIR + "/motif/{sample}/{topn}/{winsize}/{method}/log.txt"
  params:
    sample = "{sample}",
    winsize = "{winsize}",
    topn = "{topn}",
    method = "{method}"
  shell:
    """
    time dreme-py3 -p {input.input_fn1} -n {input.input_fn2} -oc {OUTPUT_DIR}/motif/{params.sample}/{params.topn}/{params.winsize}/{params.method}/ -mink 6 -maxk 8 2>&1 >{log.log_fn1}
    """


