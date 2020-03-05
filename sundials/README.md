# SUNDIALS examples

Example codes demonstrating the use of
[SUNDIALS](https://computing.llnl.gov/projects/sundials) using other XSDK
packages.

## SUNDIALS + SuperLU

The `ark_brusselator1D_FEM_sludist.cpp` example simulates a brusselator problem
from chemical kinetics. This program solves the problem with the diagonally
implicit Runge--Kutta method from the SUNDIALS ARKode time integrator package.
It uses the SUNDIALS OpenMP vector for the solution data, the SUNDIALS Newton
nonlinear solver to solve a nonlinear system at every time step, and the
[SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist)
parallel sparse-direct linear solver to solve the resulting linear system.

## SUNDIALS + PETSc

The ``cv_petsc_ex7.c`` example solves a nonlinear, time-dependent PDE in 2d. The
example is based on the [PETSc](https://www.mcs.anl.gov/petsc/) TS ex7.c, but
uses the BDF method from the SUNDIALS CVODE time integrator package. It
interfaces with the PETSc Vec for the solution data, the PETSc SNES nonlinear
solvers to solve a nonlinear system at every time step, an the PETSc KSP linear
solvers to solve the resulting linear systems.

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
