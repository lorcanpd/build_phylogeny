#!/bin/bash

#BSUB -o $PWD/%J.log
#BSUB -e $PWD/%J.err
#BSUB -q normal
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 2000
#BSUB -n 1
#BSUB -J test_singularity
#BSUB -G team274-grp

echo "Running job"
/software/singularity/v3.10.0/bin/singularity exec --bind $PWD:/mnt /nfs/users/nfs_l/lp23/sifs/manas-pipe.sif Rscript /mnt/check_env.R
echo "Job finished"


