#!/bin/bash

# compile fields.c
mpicc fields.c -lfftw3_mpi -lfftw3 -lm -O3 -Wall -o fields

# iterate
for (( i = 0 ; i <= 1000000 ; i += 100 )) ; do
	
	# visualize crystals and flow
	if [[ -e test-$i.n && -e test-$i.v && !( -e flow-$i.png ) ]] ; then
		java -jar flowplotter.jar test-$i.n test-$i.v flow-$i.png
	fi
	
	# get the orientation field
	if [[ -e test-$i.n && !( -e test-$i.phi ) ]] ; then
		mpirun -np 4 fields test-$i.n test-$i.phi 6 0.2 0.15
	fi
	
	# visualize crystals and lattice orientation
	if [[ -e test-$i.n && -e test-$i.phi && !( -e phi-$i.png ) ]] ; then
		java -jar orienter.jar test-$i.n test-$i.phi phi-$i.png
	fi
	
done
