# MFEM-PETSc example
This example solves the classical obstacle problem which models an edge clamped 
elastic membrane pulled over a rigid obstacle.  MFEM is used to discretize the
underlying Poisson equation and PETSC-TAO is used the solve to optimization
problem.  This proble also demonstrates the how MFEM and PETSc can share 
the data from vectors of each type.

This example is built to run in serial, so launch it wth your desired options:
```
./obstacle --order 2
```

Useful non-default options:
|   Flag                | Meaning                                               |
|:----------------------| :-----------------------------------------------------|
| --order n             | Set the polynomial order of the discretization.       |
| --visit               | Output VisIt files for visualation of the solution.   |