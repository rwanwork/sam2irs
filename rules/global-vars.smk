#####################################################################
##  global-vars.smk
##
##  Raymond Wan (raymondwan@cuhk.edu.hk)
##
##  Department of Medicine and Therapeutics, 
##  Division of Endocrinology & Diabetes,
##  The Chinese University of Hong Kong, Hong Kong
##
##  Copyright (C) 2019, Raymond Wan, All rights reserved.
#####################################################################


#####################################################################
##  Libraries
#####################################################################

from pathlib import Path
from getpass import getuser


#####################################################################
##  Specifications about the server
#####################################################################

##  Determine the home directory
MY_HOMEDIR = str (Path.home())

##  Get the current user name
MY_USERNAME = getuser ()

##  Set the number of threads
NUM_THREADS = 12


#####################################################################
##  Local directories
#####################################################################

INPUT_DIR = "Input"
OUTPUT_DIR = "Output"


