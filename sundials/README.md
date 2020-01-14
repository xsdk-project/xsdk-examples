# SUNDIALS examples

Example codes demonstrating the use of SUNDIALS using other XSDK packages.


## SUNDIALS + SuperLU

The `ark_brusselator1D_FEM_sludist.cpp` example simulates a brusselator
problem from chemical kinetics.
This is a PDE system with 3 components, <img alt="$Y = [u,v,w]$" src="svgs/959bb0ca5827460670264d2b3fa169d4.svg" align="middle" width="88.7766pt" height="24.56552999999997pt"/>, satisfying the equations,
<p align="center"><img alt="$$&#10;\begin{aligned}&#10;   u_t &amp;= du \cdot u_{xx} + a - (w+1)u + vu^2, \\&#10;   v_t &amp;= dv \cdot v_{xx} + wu - vu^2, \\&#10;   w_t &amp;= dw \cdot w_{xx} + (b-w)/ep - wu,&#10;\end{aligned}&#10;$$" src="svgs/4e3a4548d84ce8f1402d4960eb2b2da8.svg" align="middle" width="255.41174999999996pt" height="69.622905pt"/></p>
for <img alt="$t \in [0, 80]$" src="svgs/738a9094e7816e91b64dfe50cef848ee.svg" align="middle" width="66.916905pt" height="24.56552999999997pt"/>, <img alt="$x \in [0, 1]$" src="svgs/2510e5860f95e80cadf9cf45baa50227.svg" align="middle" width="62.174310000000006pt" height="24.56552999999997pt"/>, with initial conditions
<p align="center"><img alt="$$&#10;\begin{aligned}&#10;   u(0,x) &amp;=  a  + 0.1\sin(\pi x), \\&#10;   v(0,x) &amp;= b/a + 0.1\sin(\pi x), \\&#10;   w(0,x) &amp;=  b  + 0.1\sin(\pi x),&#10;\end{aligned}&#10;$$" src="svgs/b456576c21632aaaa540d7b365d201bb.svg" align="middle" width="196.12394999999998pt" height="65.69194499999999pt"/></p>
and with stationary boundary conditions, i.e.
<p align="center"><img alt="$$&#10;\begin{aligned}&#10;   u_t(t,0) &amp;= u_t(t,1) = 0, \\&#10;   v_t(t,0) &amp;= v_t(t,1) = 0, \\&#10;   w_t(t,0) &amp;= w_t(t,1) = 0.&#10;\end{aligned}&#10;$$" src="svgs/7145dfb04640484f0e3ddfd91571774b.svg" align="middle" width="159.92129999999997pt" height="65.69194499999999pt"/></p>
Here, we use a piecewise linear Galerkin finite element
discretization in space, where all element-wise integrals are
computed using 3-node Gaussian quadrature (since we will have
quartic polynomials in the reaction terms for the <img alt="$u_t$" src="svgs/e6897b8647f3bd38144535d3f40078e2.svg" align="middle" width="14.322330000000001pt" height="14.102549999999994pt"/> and <img alt="$v_t$" src="svgs/3e3c6ee78813607a4d976d92c19dd36e.svg" align="middle" width="12.885510000000002pt" height="14.102549999999994pt"/>
equations, including the test function).  The time derivative
terms for this system will include a mass matrix, giving rise
to an ODE system of the form
<p align="center"><img alt="$$&#10;     M y_t = L y + R(y),&#10;$$" src="svgs/feb6dcc9921528ee025c6b144635b061.svg" align="middle" width="131.642775pt" height="16.376943pt"/></p>
where <img alt="$M$" src="svgs/fb97d38bcc19230b0acd442e17db879c.svg" align="middle" width="17.67348pt" height="22.381919999999983pt"/> is the block mass matrix for each component, <img alt="$L$" src="svgs/ddcb483302ed36a59286424aa5e0be17.svg" align="middle" width="11.145420000000001pt" height="22.381919999999983pt"/> is
the block Laplace operator for each component, and <img alt="$R(y)$" src="svgs/4051c5cf4a2c287d8c463da35cb695a5.svg" align="middle" width="33.91872pt" height="24.56552999999997pt"/> is
a 3x3 block comprised of the nonlinear reaction terms for
each component.  Since it it highly inefficient to rewrite
this system as
<p align="center"><img alt="$$&#10;     y_t = M^{-1}(L y + R(y)),&#10;$$" src="svgs/7e24cd0234b30362b1a26b71378daf1f.svg" align="middle" width="162.028845pt" height="18.269295pt"/></p>
we solve this system using ARKStep, with a user-supplied mass
matrix.  We therefore provide functions to evaluate the ODE RHS
<p align="center"><img alt="$$&#10;   f(t,y) = L y + R(y),&#10;$$" src="svgs/48173febf3555b003e1f02241d2be7e0.svg" align="middle" width="144.49248pt" height="16.376943pt"/></p>
its Jacobian
<p align="center"><img alt="$$&#10;   J(t,y) = L + \frac{dR}{dy},&#10;$$" src="svgs/2ab039b6508ad5aa0e3a4083cb6f3b58.svg" align="middle" width="127.93934999999999pt" height="36.953894999999996pt"/></p>
and the mass matrix, <img alt="$M$" src="svgs/fb97d38bcc19230b0acd442e17db879c.svg" align="middle" width="17.67348pt" height="22.381919999999983pt"/>.
This program solves the problem with the DIRK method, using a
Newton iteration with the SuperLU_DIST SUNLinearSolver.
100 outputs are printed at equal time intervals, and run
statistics are printed at the end.

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