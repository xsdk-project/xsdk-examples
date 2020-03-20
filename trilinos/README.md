# Trilinos examples

Example codes demonstrating the use of
[Trilinos](https://trilinos.github.io) using other XSDK
packages.

## Trilinos + SuperLU

'SimpleSolve_WithParameters.cpp` is an example of sparse linear system direct solvers
allocated to single MPI rank. The code uses Amesos2 (interface for sparse 
direct solver libraries) APIs to call
[SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist).

The source code ptovides the examples of solver parameters.
    Teuchos::ParameterList superludist_params = amesos2_params.sublist("SuperLU_DIST");
    superludist_params.set("Trans","No","Whether to solve with A^T");
    superludist_params.set("npcol",1,"Number of Processor Columns");
    superludist_params.set("nprow",1,"Number of Processor Rows");
    superludist_params.set("ColPerm","NATURAL","Use 'natural' ordering of columns");


**Usage**

```
mpirun -np 1 ./SimpleSolve_WithParameters
```
