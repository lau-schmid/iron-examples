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

#mkdir -p results/current_run/l1x1_n4_i1_s0_01 && ./$folder/src/example 4 1 0 0.01 10 1 hodgkin_huxley_1952.cellml F 0.005 && mv *.ex* results/current_run/l1x1_n4_i1_s0_01
#mkdir -p results/current_run/l1x1_n4_i1_s0_05 && ./$folder/src/example 4 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n4_i1_s0_05
#mkdir -p results/current_run/l1x1_n6_i1_s0_05 && ./$folder/src/example 6 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n6_i1_s0_05
#mkdir -p results/current_run/l1x1_n7_i1_s0_05 && ./$folder/src/example 7 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n7_i1_s0_05
mkdir -p results/current_run/l1x1_n8_i1_s0_05 && ./$folder/src/example 8 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n8_i1_s0_05
#mkdir -p results/current_run/l1x1_n16_i1_s0_05 && ./$folder/src/example 16 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n16_i1_s0_05
#mkdir -p results/current_run/l1x1_n18_i1_s0_05 && ./$folder/src/example 18 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n18_i1_s0_05
#mkdir -p results/current_run/l1x1_n32_i1_s0_05 && ./$folder/src/example 32 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.025 && mv *.ex* results/current_run/l1x1_n32_i1_s0_05
#mkdir -p results/current_run/l1x1_n2048_i1_s0_005 && ./$folder/src/example 2048 1 0 0.005 5 1 hodgkin_huxley_1952.cellml F 0.0025 && mv *.ex* results/current_run/l1x1_n2048_i1_s0_005
