# MFEM-PUMI example

This example, `adapt.cpp`, illustrates the integration of PUMI with MFEM and it
is based on the parallel PUMI-based version of example 6 from MFEM.

For a description of the problem being solved, see the documentation in the
beginning of the source file.

Sample runs:
```
mpirun -np 8 ./adapt
```

For a full list of options, see
```
./adapt -h
```
