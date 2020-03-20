# xSDK Examples

Example codes demonstrating the use of of various xSDK libraries in tandem to solve problems of 
interest.  Each of the library folders has one or more examples codes that are built of that library 
and utilize the code integrations with the other xSDK libraries.  Running these example codes and
examining the output is a ood way to better understand how these libraries can work together, and the
code samples are a good place to start for new projects.  More details about the examples can be found 
in the README.md files in these folders.

## Code Example Summary

|   Example                                  | Libraries                   | Description                                       |
|:-------------------------------------------|:----------------------------|:--------------------------------------------------|
|  hypre/ij_laplacian.c                      | HYPRE+SuperLU_Dist          | 2D Laplacian problem                              |
|  libensemble/test_persistent_aposmm_tao.py | libEnsemble+PETSc           | 2D constrained optimization problem               |
|  mfem/ginkgo/mfem_ex1_gko.cpp              | MFEM+Ginkgo                 | 2D Poisson problem with Ginko solver              |
|  mfem/hypre-superlu/convdiff.cpp           | MFEM+HYPRE+SuperLU_Dist     | 2D steady state convective diffusion              |
|  mfem/petsc/obstacle.cpp                   | MFEM+PETSc                  | Membrane obstacle problem (min energy functional) |
|  mfem/sundials/transient-heat.cpp          | MFEM+Sundials               | 2D Transient nonlinear heat conduction            |
|  petsc/ex19.c                              | PETSc+HYPRE+SuperLU_Dist    | 2D nonlinear driven cavity problem                |
|  sundials/ark_brusselator1D_FEM_sludist.cpp| SUNDIALS+SuperLU_Dist       | Chemical kinetics brusselator problem             |
|  sundials/cv_petsc_ex7.c                   | SUNDIALS+PETSc              | 2D nonlinear PDE solution                         |
|  trilinos/SimpleSolve_WithParameters.cpp   | Trilinos+SuperLU_Dist       | Small linear system direct solution               |

## Install the code samples

The examples can be installed along with the xSDK utilizing the spack package.
```
spack install xsdk
spack install xsdk-examples
```

Further details on how to run each example can be found in each exampl folder's README.md file.