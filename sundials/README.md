# SUNDIALS examples

Example codes demonstrating the use of
[SUNDIALS](https://computing.llnl.gov/projects/sundials) using other XSDK
packages.

## SUNDIALS + SuperLU

The `ark_brusselator1D_FEM_sludist.cpp` example simulates a brusselator problem
from chemical kinetics. This program solves the problem with the diagonally
implicit Runge--Kutta method from the SUNDIALS ARKode time integrator package.
It uses the SUNDIALS serial vector for the solution data, the SUNDIALS Newton
nonlinear solver to solve a nonlinear system at every time step, and the
[SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist) parallel sparse-direct
linear solver to solve the resulting linear system. Jacobian data is stored in
a SuperLU_DIST SuperMatrix. 100 outputs are printed to an output file at equal time
intervals, and run statistics are printed to stdout at the end.

**Usage**

```
mpirun -np 1 ./ark_brusselator1D_FEM_sludist
```

## SUNDIALS + PETSc

The ``cv_petsc_ex7.c`` example solves a nonlinear, time-dependent PDE in 2d. The
example is based on the [PETSc](https://www.mcs.anl.gov/petsc/) TS ex7.c, but
uses the BDF method from the SUNDIALS CVODE time integrator package. It
interfaces with the PETSc Vec for the solution data, the PETSc SNES nonlinear
solvers to solve a nonlinear system at every time step, an the PETSc KSP linear
solvers to solve the resulting linear systems. Output, nonlinear/linear solver
and vector options can be controlled using the typical PETSc runtime arguments.

**Usage**

To run with the default options:

```
./cv_petsc_ex7
```

To monitor the KSP linear solver progress:

```
./cv_petsc_ex7 -ksp_monitor
```

To use the Anderson method from SNES instead of the default Newton line search: 

```
./cv_petsc_ex7 -snes_type anderson
```

To view all the options available:

```
./cv_petsc_ex7 -h
```

## SUNDIALS + MAGMA

The ``cvRoberts_blockdiag_magma.cpp`` example solves a group of independent ODEs
simulating a problem from chemical kinetics. The example groups the independent ODEs together
to form a larger system which is evolved with the CVODE implicit BDF time integration
method. The linear system that arises is block-diagonal and there is no coupling
between the blocks. Thus, a batched LU method from the [MAGMA](https://icl.utk.edu/magma/) linear algebra library
is utilized. This example is CUDA-enabled.

To run:

```
./cvRoberts_blockdiag_magma [number of groups]
```


## License

SUNDIALS is released under the BSD 3-clause license. See the
[LICENSE](./LICENSE) and [NOTICE](./NOTICE) files for details. All new
contributions must be made under the BSD 3-clause license.

**Please Note** If you are using SUNDIALS with any third party libraries linked
in (e.g., LAPACK, KLU, SuperLU_MT, PETSc, or *hypre*), be sure to review the
respective license of the package as that license may have more restrictive
terms than the SUNDIALS license.

```text
SPDX-License-Identifier: BSD-3-Clause

LLNL-CODE-667205  (ARKODE)
UCRL-CODE-155951  (CVODE)
UCRL-CODE-155950  (CVODES)
UCRL-CODE-155952  (IDA)
UCRL-CODE-237203  (IDAS)
LLNL-CODE-665877  (KINSOL)
```
