#!/bin/bash
# test cases for example-0404 (1D problem with Hodgkin-Huxley)

echo "compiling and running example $(pwd)"

folder=$1

mkdir -p $folder

echo "  compiling $folder"
cd $folder
cmake -DCMAKE_BUILD_TYPE=$folder -DOPENCMISS_BUILD_TYPE=$folder ..
make
cd ..
echo "  running $folder"
# <number elements X> <interpolation type> <solver type> <PDE step size> <stop time> <output frequency> <CellML Model URL> <slow-twitch> <ODE time-step>

#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt005 && ./$folder/src/example 64 1 0 0.005 10 10 hodgkin_huxley_1952.cellml T 0.001 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt005
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt001 && ./$folder/src/example 64 1 0 0.001 10 50 hodgkin_huxley_1952.cellml T 0.0002 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt001
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt0005 && ./$folder/src/example 64 1 0 0.0005 10 100 hodgkin_huxley_1952.cellml T 0.0001 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt0005
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt00025 && ./$folder/src/example 64 1 0 0.00025 10 200 hodgkin_huxley_1952.cellml T 0.00005 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt00025


mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt005 && ./$folder/src/example 64 1 0 0.005 10 10 hodgkin_huxley_1952.cellml T 0.005 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt005
mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt001 && ./$folder/src/example 64 1 0 0.001 10 50 hodgkin_huxley_1952.cellml T 0.001 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt001
mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt0005 && ./$folder/src/example 64 1 0 0.0005 10 100 hodgkin_huxley_1952.cellml T 0.0005 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt0005
mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt00025 && ./$folder/src/example 64 1 0 0.00025 10 200 hodgkin_huxley_1952.cellml T 0.00025 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt00025


#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt01 && ./$folder/src/example 64 1 0 0.01 10 5 hodgkin_huxley_1952.cellml F 0.01 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt01
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt005 && ./$folder/src/example 64 1 0 0.005 10 10 hodgkin_huxley_1952.cellml F 0.005 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt005
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt001 && ./$folder/src/example 64 1 0 0.001 10 50 hodgkin_huxley_1952.cellml F 0.001 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt001
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt0005 && ./$folder/src/example 64 1 0 0.0005 10 100 hodgkin_huxley_1952.cellml F 0.0005 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt0005
#mkdir -p results/current_run/l1x1_n64_i1_s0_O1_dt00025 && ./$folder/src/example 64 1 0 0.00025 10 200 hodgkin_huxley_1952.cellml F 0.00025 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O1_dt00025
