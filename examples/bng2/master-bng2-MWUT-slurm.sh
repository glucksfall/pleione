#!/bin/sh

#SBATCH --no-requeue
#SBATCH --partition=cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --job-name=pleione-bng2
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt

export PYTHONPATH="$PYTHONPATH:$HOME/opt/github-repositories/glucksfall.pleione"

MODEL=pysbmodel-example6-bng2.bngl

PARTITION=$SLURM_JOB_PARTITION
DATA=../exp-data/bng2/data-*.txt

NUM_ITER=100
NUM_SIMS=10
POP_SIZE=100
POP_BEST=0

SWAP=0.5
RATE=0.5
ERROR="MWUT"
UTABLE=./ucrit.txt

~/bin/python3 -m pleione.bng2 --model=$MODEL \
--iter=$NUM_ITER --inds=$POP_SIZE --sims=$NUM_SIMS --best=$POP_BEST \
--data=$DATA --rate=$RATE --swap=$SWAP --error=$ERROR --crit=$UTABLE \
--slurm=$PARTITION --seed=0
