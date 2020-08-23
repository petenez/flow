# compile pfc-flow.c
gcc pfc-flow.c -fopenmp -lfftw3_omp -lfftw3 -lm -O3 -Wall -o pfc-flow

# run pfc-flow
./pfc-flow
