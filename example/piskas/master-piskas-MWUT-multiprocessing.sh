#!/bin/sh
export PYTHONPATH="$PYTHONPATH:$HOME/opt/github-repositories/glucksfall.pleione"

MODEL=pysbmodel-example6-piskas.kappa
FINAL=660
STEPS=66

DATA=../exp-data/piskas/data-*.txt

NUM_ITER=100
NUM_SIMS=10
POP_SIZE=100
POP_BEST=0

SWAP=0.5
RATE=0.5
ERROR="MWUT"
UTABLE=./ucrit.txt

~/bin/python3 -m pleione.piskas --model=$MODEL --final=$FINAL --steps=$STEPS \
--iter=$NUM_ITER --inds=$POP_SIZE --sims=$NUM_SIMS --best=$POP_BEST \
--data=$DATA --rate=$RATE --swap=$SWAP --error=$ERROR --crit=$UTABLE --seed=0
