#####################################################################
##  Snakefile.smk
##
##  Raymond Wan (raymondwan@cuhk.edu.hk)
##
##  Department of Medicine and Therapeutics, 
##  Division of Endocrinology & Diabetes,
##  The Chinese University of Hong Kong, Hong Kong
##
##  Copyright (C) 2019, Raymond Wan, All rights reserved.
#####################################################################

from snakemake.utils import validate, min_version

##  Set the minimum snakemake version
min_version("5.7")

##  Define global constraints on wildcards
wildcard_constraints:
  topn = "\d+",
  winsize = "\d+",
  sample = "\w+",
  method = "upstream|downstream|intron"

##  Include additional functions and rules
include:  "rules/global-vars.smk"
include:  "rules/other-rules.smk"


##  Set the shell and the prefix to run before "shell" commands -- activate the conda environment
shell.executable ("/bin/bash")
shell.prefix ("source /opt/miniconda3/etc/profile.d/conda.sh; conda activate sam2irs; ")


##  Top-level rule
rule all:
  input: 
    OUTPUT_DIR + "/Complete/ASC_50_1000.done",
    OUTPUT_DIR + "/Complete/DEK_50_1000.done",
    OUTPUT_DIR + "/Complete/fSC_50_1000.done",
    OUTPUT_DIR + "/Complete/pSC_50_1000.done",
    OUTPUT_DIR + "/Complete/uSC_50_1000.done",
    OUTPUT_DIR + "/Complete/YFP_50_1000.done",
    OUTPUT_DIR + "/Complete/ASC_100_1000.done",
    OUTPUT_DIR + "/Complete/DEK_100_1000.done",
    OUTPUT_DIR + "/Complete/fSC_100_1000.done",
    OUTPUT_DIR + "/Complete/pSC_100_1000.done",
    OUTPUT_DIR + "/Complete/uSC_100_1000.done",
    OUTPUT_DIR + "/Complete/YFP_100_1000.done",
    OUTPUT_DIR + "/Complete/ASC_200_1000.done",
    OUTPUT_DIR + "/Complete/DEK_200_1000.done",
    OUTPUT_DIR + "/Complete/fSC_200_1000.done",
    OUTPUT_DIR + "/Complete/pSC_200_1000.done",
    OUTPUT_DIR + "/Complete/uSC_200_1000.done",
    OUTPUT_DIR + "/Complete/YFP_200_1000.done",
    OUTPUT_DIR + "/Complete/ASC_400_1000.done",
    OUTPUT_DIR + "/Complete/DEK_400_1000.done",
    OUTPUT_DIR + "/Complete/fSC_400_1000.done",
    OUTPUT_DIR + "/Complete/pSC_400_1000.done",
    OUTPUT_DIR + "/Complete/uSC_400_1000.done",
    OUTPUT_DIR + "/Complete/YFP_400_1000.done",
    OUTPUT_DIR + "/Complete/ASC_50_500.done",
    OUTPUT_DIR + "/Complete/DEK_50_500.done",
    OUTPUT_DIR + "/Complete/fSC_50_500.done",
    OUTPUT_DIR + "/Complete/pSC_50_500.done",
    OUTPUT_DIR + "/Complete/uSC_50_500.done",
    OUTPUT_DIR + "/Complete/YFP_50_500.done",
    OUTPUT_DIR + "/Complete/ASC_100_500.done",
    OUTPUT_DIR + "/Complete/DEK_100_500.done",
    OUTPUT_DIR + "/Complete/fSC_100_500.done",
    OUTPUT_DIR + "/Complete/pSC_100_500.done",
    OUTPUT_DIR + "/Complete/uSC_100_500.done",
    OUTPUT_DIR + "/Complete/YFP_100_500.done",
    OUTPUT_DIR + "/Complete/ASC_200_500.done",
    OUTPUT_DIR + "/Complete/DEK_200_500.done",
    OUTPUT_DIR + "/Complete/fSC_200_500.done",
    OUTPUT_DIR + "/Complete/pSC_200_500.done",
    OUTPUT_DIR + "/Complete/uSC_200_500.done",
    OUTPUT_DIR + "/Complete/YFP_200_500.done",
    OUTPUT_DIR + "/Complete/ASC_400_500.done",
    OUTPUT_DIR + "/Complete/DEK_400_500.done",
    OUTPUT_DIR + "/Complete/fSC_400_500.done",
    OUTPUT_DIR + "/Complete/pSC_400_500.done",
    OUTPUT_DIR + "/Complete/uSC_400_500.done",
    OUTPUT_DIR + "/Complete/YFP_400_500.done"


rule Combine_Results:
  input:
    input_fn1 = OUTPUT_DIR + "/motif/{sample}/{topn}/{winsize}/upstream/dreme.html",
    input_fn2 = OUTPUT_DIR + "/motif/{sample}/{topn}/{winsize}/downstream/dreme.html",
    input_fn3 = OUTPUT_DIR + "/motif/{sample}/{topn}/{winsize}/intron/dreme.html"
  output:
    output_fn1 = OUTPUT_DIR + "/Complete/{sample}_{winsize}_{topn}.done"
  shell:
    """
    touch {output.output_fn1}
    """


