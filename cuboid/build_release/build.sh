alias ccc='rm -rf CMakeFiles/ CMakeCache.txt cmake_install.cmake Makefile CTestTestfile.cmake export mpi_verification packaging Tests support'
ccc && cmake -DOPENCMISS_BUILD_TYPE=RELEASE -DCMAKE_BUILD_TYPE=RELEASE  .. && make clean && make all
