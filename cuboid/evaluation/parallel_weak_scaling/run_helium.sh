# weak scaling

while true; do
  echo "run quietly, output to out*.txt"
  mpirun -n 1  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 3 4  1 1 2>&1 | tee -a out1.txt
  mpirun -n 2  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 3 4  2 1 2>&1 | tee -a out2.txt
  mpirun -n 4  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 3 8  2 1 2>&1 | tee -a out4.txt
  mpirun -n 8  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 6 8  2 1 2>&1 | tee -a out8.txt
  mpirun -n 12 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 6 8  3 1 2>&1 | tee -a out12.txt
  mpirun -n 16 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 6 8  4 1 2>&1 | tee -a out16.txt
  mpirun -n 24 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 6 12 4 1 2>&1 | tee -a out24.txt
  mpirun -n 32 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 6 16 4 1 2>&1 | tee -a out32.txt
  mpirun -n 48 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 9 16 4 1 2>&1 | tee -a out48.txt
done
