# CMake generated Testfile for 
# Source directory: /usr/local/home/kraemer/software/opencmiss_examples/classicalfield_laplace_simple/Fortran
# Build directory: /usr/local/home/kraemer/software/opencmiss_examples/classicalfield_laplace_simple/build_release/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Laplace_Fortran "/usr/local/home/kraemer/software/opencmiss_examples/classicalfield_laplace_simple/build_release/Fortran/laplace_fortran")
set_tests_properties(Laplace_Fortran PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/home/kraemer/software/opencmiss/iron/install/x86_64_linux/gnu-4.9-F4.9/openmpi_release/release/bin:")
