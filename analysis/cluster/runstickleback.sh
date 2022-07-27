#!/bin/bash
#
#SBATCH --job-name=runstickleback
#
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL

# All the stickleback (Py) and sticklbackms (R) dependencies, plus devtools (R)
# have to be installed first!

cd $HOME/repos/sticklebackms
ml R/4.1
ml python/3.9

Rscript analysis/cluster/runstickleback.R lunge 50 5 8 2 8 4 8 200 20 20
