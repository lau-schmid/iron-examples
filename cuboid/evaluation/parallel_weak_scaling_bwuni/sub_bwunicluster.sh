#!/bin/bash
#MSUB -l nodes=16:ppn=16
#MSUB -l walltime=12:00:00
#MSUB -l mem=32gb
#MSUB -N weak_scaling_on_16x16
#MSUB -o parallel_weak_scaling_16x16-bwd.txt
#MSUB -q multinode

export OMP_NUM_THREADS=${MOAB_PROCCOUNT}
#module load compiler/gnu/4.9

while true; do
  #/home/st/st_st/st_smt86072/opencmiss/examples/cuboid/evaluation/parallel_weak_scaling_bwuni/run.sh
  /home/st/st_st/st_smt86072/opencmiss/examples/cuboid/evaluation/parallel_weak_scaling_bwuni/run_backwards.sh
done
