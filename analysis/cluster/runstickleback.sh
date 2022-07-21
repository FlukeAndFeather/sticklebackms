#!/bin/bash
#
#SBATCH --job-name=runstickleback
#
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL

cd $HOME/repos/sticklebackms
ml R
ml python/3.9

SIF=$GROUP_HOME/$USER/sing/sb-test_latest.sif
SCRIPT=stickleback_test.py
DATA=$SCRATCH/stickleback/bw_breaths_sensors_events.pkl
WINDOW=50
FOLDS=2
TREES=64

ml python/3.9
cd $HOME/repos/run-stickleback
echo WINDOW=$WINDOW FOLDS=$FOLDS TREES=$TREES
python3 $SCRIPT $DATA $WINDOW $FOLDS $TREES $SLURM_CPUS_PER_TASK

# Verify writable directories exist
# mkdir -p $SCRATCH/.numba_cache
# MPLCONFIGDIR=$SCRATCH/.mplconfig

# singularity exec --env-file stickleback_env $SIF python3 $SCRIPT $DATA $WINDOW $FOLDS $TREES $SLURM_CPUS_PER_TASK
