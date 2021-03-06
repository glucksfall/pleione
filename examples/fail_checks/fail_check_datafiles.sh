#!/bin/sh

#SBATCH --no-requeue
#SBATCH --partition=spica

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --job-name=pleione-kasim
#SBATCH --output=stdout_python_fail.txt
#SBATCH --error=stderr_python_fail.txt

export PYTHONPATH="$PYTHONPATH:$HOME/opt/github-repositories/glucksfall.pleione"

MODEL=pysbmodel-example6-kasim.kappa
FINAL=660
STEPS=10

PARTITION=$SLURM_JOB_PARTITION
DATA=./exp-data/kasim/data-*.txt
# DATA="data1 data2 data3"

NUM_ITER=100
NUM_SIMS=10
POP_SIZE=100
POP_BEST=0

SWAP=0.5
RATE=0.5
ERROR="MWUT"
UTABLE=./ucrit.txt

python3 -m pleione.kasim --model=$MODEL --final=$FINAL --steps=$STEPS \
--iter=$NUM_ITER --inds=$POP_SIZE --sims=$NUM_SIMS --best=$POP_BEST \
--data $DATA --rate=$RATE --swap=$SWAP --error=$ERROR --crit=$UTABLE \
--syntax=3 --seed=0 --python='python3'
