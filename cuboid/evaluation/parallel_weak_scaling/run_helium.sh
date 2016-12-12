# weak scaling

while true; do
  mpirun -n 1 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 3 4 1 1 > out1.txt
  mpirun -n 2 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 3 4 2 1 > out2.txt
  mpirun -n 4 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 3 8 2 1 > out4.txt
  mpirun -n 8 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 8 2 1 > out8.txt
  mpirun -n 12 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 8 3 1 > out12.txt
  mpirun -n 16 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 8 4 1 > out16.txt
  mpirun -n 24 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 12 4 1 > out24.txt
  mpirun -n 32 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 16 4 1 > out32.txt
  mpirun -n 48 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 9 16 4 1 > out48.txt
done
