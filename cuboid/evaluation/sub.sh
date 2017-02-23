SCENARIO=parallel_weak_scaling

echo "OPENCMISS_EVALUATION_DIR=$OPENCMISS_EVALUATION_DIR"
echo "SCENARIO=$SCENARIO"


#qsub -q fast.q $OPENCMISS_EVALUATION_DIR/$SCENARIO/run.sh -j eo -o $OPENCMISS_EVALUATION_DIR/$SCENARIO
export nproc=96
qsub -q full.q -pe full_mpi_rr $nproc $OPENCMISS_EVALUATION_DIR/$SCENARIO/run.sh -o $OPENCMISS_EVALUATION_DIR/$SCENARIO

#qsub /home/maierbn/opencmiss/examples/cuboid/evaluation/testmpi.sh


# full.q : alle Rechenknoten
# big.q : leadx01 mit 64 Cores und 2.3 GHz
# fast.q : leadx02 bis leadx07 mit 32 Cores und 3.2GHz
# s_x2.q bis s_x7.q: einzelne Rechenknoten


# view queues in GUI:
# module load tools/1; tempo
# view loads and jobs in shell:
# qstat -f

# SGE manual: http://gridscheduler.sourceforge.net/htmlman/manuals.html

# programming environments: 
# qconf -spl
# abaqus_big
# abaqus_ff
# abaqus_rr
# abaqus_small
# big_open_mpi   (only 1 node with 64 cores)
# full_mpi_ff    (whole fillup)
# full_mpi_rr    (whole round robin)
# make
# matlab
# matlab_rr
# mpi_ff
# mpi_rr

# queues
# $ qconf -sql
# all.q
# big.q
# fast.q
# full.q
# s_x2.q
# s_x3.q
# s_x4.q
# s_x5.q
# s_x6.q
# s_x7.q
# single.q

# arguments of qsub: -pe <pe> <n>-<m>,  -q <queue>    
