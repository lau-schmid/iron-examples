#!/bin/bash 
# multi process execution with parameters

. ~/.bashrc
#module add openmpi


# read in command line parameters
nproc=1
x=3
y=4
z=1
f=1
a=1
SCENARIO=parallel_weak_scaling

if [[ $# -eq 0 ]]; then
  echo ""
  echo "usage: program <nproc> <x> <y> <z> <f> <a> <scenario>"
  echo ""
  echo ""
fi

if [[ $# -gt 0 ]]; then
  nproc=$1
  shift
fi
if [[ $# -gt 0 ]]; then
  x=$1
  shift
fi
if [[ $# -gt 0 ]]; then
  y=$1
  shift
fi
if [[ $# -gt 0 ]]; then
  z=$1
  shift
fi
if [[ $# -gt 0 ]]; then
  f=$1
  shift
fi
if [[ $# -gt 0 ]]; then
  a=$1
  shift
fi
if [[ $# -gt 0 ]]; then
  SCENARIO=$1
  shift
fi


#echo "[run_with_parameters] OPENCMISS_REL_DIR=$OPENCMISS_REL_DIR, OPENCMISS_INPUT_DIR=$OPENCMISS_INPUT_DIR"
echo "[run_with_parameters] nproc=$nproc, x=$x, y=$y, z=$z, f=$f, a=$a, SCENARIO=$SCENARIO"

export command="mpiexec -n $nproc $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR $x $y $z $f $a | tee $OPENCMISS_EVALUATION_DIR/$SCENARIO/out_$nproc.txt"
echo "$command"
ulimit -a
echo "modules:"
module list
echo "$(date) run_with_parameters.sh nproc=$nproc, x=$x, y=$y, z=$z, f=$f, a=$a, SCENARIO=$SCENARIO memlock limit=$(ulimit -l)" >> $OPENCMISS_REL_DIR/log.txt
$command
