# MFEM-HiOp example

This example, `adv.cpp`, illustrates the integration of HiOp with MFEM and it is
based on the parallel HiOp-extended version of example 9 from MFEM.

For a description of the problem being solved, see the documentation in the
beginning of the source file.

Sample runs:
```
mpirun -np 4 ./adv
mpirun -np 4 ./adv -m ../data/periodic-segment.mesh -rs 3 -p 0 -o 2 -dt 0.002
mpirun -np 4 ./adv -m ../data/disc-nurbs.mesh -p 1 -rs 2 -dt 0.005 -tf 9
mpirun -np 4 ./adv -m ../data/periodic-square.mesh -p 3 -rs 3 -dt 0.0025 -tf 9
mpirun -np 4 ./adv -m ../data/periodic-cube.mesh -p 0 -rs 2 -o 2 -dt 0.02 -tf 8
```

For a full list of options, see
```
./adv -h
```
