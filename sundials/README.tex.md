# SUNDIALS examples

Example codes demonstrating the use of SUNDIALS using other XSDK packages.


## SUNDIALS + SuperLU

The `ark_brusselator1D_FEM_sludist.cpp` example simulates a brusselator
problem from chemical kinetics.
This is a PDE system with 3 components, $Y = [u,v,w]$, satisfying the equations,
$$
\begin{aligned}
   u_t &= du \cdot u_{xx} + a - (w+1)u + vu^2, \\
   v_t &= dv \cdot v_{xx} + wu - vu^2, \\
   w_t &= dw \cdot w_{xx} + (b-w)/ep - wu,
\end{aligned}
$$
for $t \in [0, 80]$, $x \in [0, 1]$, with initial conditions
$$
\begin{aligned}
   u(0,x) &=  a  + 0.1\sin(\pi x), \\
   v(0,x) &= b/a + 0.1\sin(\pi x), \\
   w(0,x) &=  b  + 0.1\sin(\pi x),
\end{aligned}
$$
and with stationary boundary conditions, i.e.
$$
\begin{aligned}
   u_t(t,0) &= u_t(t,1) = 0, \\
   v_t(t,0) &= v_t(t,1) = 0, \\
   w_t(t,0) &= w_t(t,1) = 0.
\end{aligned}
$$
Here, we use a piecewise linear Galerkin finite element
discretization in space, where all element-wise integrals are
computed using 3-node Gaussian quadrature (since we will have
quartic polynomials in the reaction terms for the $u_t$ and $v_t$
equations, including the test function).  The time derivative
terms for this system will include a mass matrix, giving rise
to an ODE system of the form
$$
     M y_t = L y + R(y),
$$
where $M$ is the block mass matrix for each component, $L$ is
the block Laplace operator for each component, and $R(y)$ is
a 3x3 block comprised of the nonlinear reaction terms for
each component.  Since it it highly inefficient to rewrite
this system as
$$
     y_t = M^{-1}(L y + R(y)),
$$
we solve this system using ARKStep, with a user-supplied mass
matrix.  We therefore provide functions to evaluate the ODE RHS
$$
   f(t,y) = L y + R(y),
$$
its Jacobian
$$
   J(t,y) = L + \frac{dR}{dy},
$$
and the mass matrix, $M$.
This program solves the problem with the DIRK method, using a
Newton iteration with the SuperLU_DIST SUNLinearSolver.
100 outputs are printed at equal time intervals, and run
statistics are printed at the end.

## SUNDIALS + PETSc

``cv_petsc_ex7.c``

## License

SUNDIALS is released under the BSD 3-clause license. See the [LICENSE](./LICENSE)
and [NOTICE](./NOTICE) files for details. All new contributions must be made
under the BSD 3-clause license.

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
