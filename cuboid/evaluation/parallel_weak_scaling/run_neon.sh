# weak scaling

  echo "run quietly, output to out*.txt"
  #mpirun -n 128 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 8 16 8 1 2>&1 | tee -a out128.txt
  #mpirun -n 96 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 8 16 6 1 > out96.txt
  mpirun -n 64 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 8 16 4 1 2>&1 | tee -a out064.txt
  mpirun -n 48 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 6 16 4 1 2>&1 | tee -a out48.txt
  mpirun -n 32 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 4 16 4 1 2>&1 | tee -a out32.txt
  mpirun -n 24 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 4 12 4 1 2>&1 | tee -a out24.txt
  mpirun -n 1  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 2 4  1 1 2>&1 | tee -a out1.txt   &
  mpirun -n 2  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 2 4  2 1 2>&1 | tee -a out2.txt   &
  mpirun -n 4  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 2 8  2 1 2>&1 | tee -a out4.txt   &
  mpirun -n 16 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 4 8  4 1 2>&1 | tee -a out16.txt
  mpirun -n 8  $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 4 8  2 1 2>&1 | tee -a out8.txt   &
  mpirun -n 12 $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 4 8  3 1 2>&1 | tee -a out12.txt  &
