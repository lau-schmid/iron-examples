# CMake generated Testfile for 
# Source directory: /usr/local/home/kraemer/software/iron-examples/PoissonExample/Fortran
# Build directory: /usr/local/home/kraemer/software/iron-examples/PoissonExample/release_build/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Poisson "/usr/local/home/kraemer/software/iron-examples/PoissonExample/release_build/Fortran/poisson")
set_tests_properties(Poisson PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/home/kraemer/software/OpenCMISS/iron/install/x86_64_linux/gnu-4.9-F4.9/openmpi_release/release/bin:")
