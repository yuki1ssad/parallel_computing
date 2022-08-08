
# Martix_Multiplication with openmp
The examples are written for Linux, not Windows.

# How-to
## Compilation
The simplest way to configure is to run command:
```
g++ -shared Matrix.cpp -o libMatrix.so -I$PWD -fPIC -fopenmp
g++ -fopenmp -o test_omp test_omp.cpp -I$PWD -L$PWD -l Matrix
```

## Run
The following command runs the compiled codes.

```
LD_LIBRARY_PATH=$PWD ./test_omp
```
