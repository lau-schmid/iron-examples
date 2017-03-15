echo "start"
gcc --version
g++ --version 

gcc -shared -o libfunc.so -fPIC func.c
g++ main.cpp -ldl -o main && ./main ./libfunc.so
