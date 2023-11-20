#!/bin/bash

#BSUB -o $PWD/%J.log
#BSUB -e $PWD/%J.err
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 8000

/software/R-4.1.3/bin/Rscript check_env.R
