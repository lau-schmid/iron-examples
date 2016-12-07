# strong scaling

mpirun -n 64 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 32 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 16 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 12 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 8 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 4 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 2 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1
mpirun -n 1 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 8 8 8 1

