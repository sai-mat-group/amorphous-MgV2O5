#!/bin/bash

export OMP_NUM_THREADS=1

#/home/vijaychoyal/Lmlip/mlip-2/bin/mlp convert-cfg relax.cfg input.pos --output-format=lammps-datafile
/home/ns/Lmlip/interface-lammps-mlip-2/lmp_mpi -in md.lmp
