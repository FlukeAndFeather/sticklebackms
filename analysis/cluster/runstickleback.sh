#!/bin/bash
#
#SBATCH --job-name=runstickleback
#SBATCH -p hns
#
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL

# All the stickleback (Py) and sticklbackms (R) dependencies, plus devtools (R)
# have to be installed first!

cd $SCRATCH/repos/sticklebackms
ml R/4.1
ml python/3.9

Rscript analysis/cluster/runstickleback.R breath 60 3 8 2 8 4 8 200 20 20
