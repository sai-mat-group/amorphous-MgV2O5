#!/bin/bash
cp pot.mtp curr.mtp
mv pot.mtp init.mtp
cp train.cfg train_init.cfg
export OMP_NUM_THREADS=1

/home/ns/Lmlip/mlip-2/bin/mlp calc-grade curr.mtp train.cfg train.cfg out.cfg --als-filename=state.als

