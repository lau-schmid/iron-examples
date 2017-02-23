# weak scaling to 256


mpirun -n 256 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 12 32 8 1
mpirun -n 192 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 12 24 8 1
mpirun -n 128 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 12 16 8 1
mpirun -n 96 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 12 16 6 1
mpirun -n 64 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 12 16 4 1
mpirun -n 48 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 9 16 4 1
mpirun -n 32 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 16 4 1
mpirun -n 24 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 12 4 1
mpirun -n 16 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 8 4 1
mpirun -n 12 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 8 3 1
mpirun -n 8 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 6 8 2 1
mpirun -n 4 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 3 8 2 1
mpirun -n 2 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 3 4 2 1
mpirun -n 1 $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR 3 4 1 1
