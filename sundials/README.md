# SUNDIALS examples

Example codes demonstrating the use of SUNDIALS using other XSDK packages.


## SUNDIALS + SuperLU

The `ark_brusselator1D_FEM_sludist.cpp` example simulates a brusselator
problem from chemical kinetics.
This is a PDE system with 3 components, <img src="/sundials/svgs/959bb0ca5827460670264d2b3fa169d4.svg?invert_in_darkmode" align=middle width=88.7766pt height=24.56552999999997pt/>, satisfying the equations,
<p align="center"><img src="/sundials/svgs/4e3a4548d84ce8f1402d4960eb2b2da8.svg?invert_in_darkmode" align=middle width=255.41174999999996pt height=69.622905pt/></p>
for <img src="/sundials/svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.9139630000000025pt height=20.14650000000001pt/> in <img src="/sundials/svgs/309241d5d19086f1455c7cbb72bc3e21.svg?invert_in_darkmode" align=middle width=40.95267pt height=24.56552999999997pt/>, <img src="/sundials/svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.359955000000003pt height=14.102549999999994pt/> in <img src="/sundials/svgs/e88c070a4a52572ef1d5792a341c0900.svg?invert_in_darkmode" align=middle width=32.764215pt height=24.56552999999997pt/>, with initial conditions
<p align="center"><img src="/sundials/svgs/b456576c21632aaaa540d7b365d201bb.svg?invert_in_darkmode" align=middle width=196.12394999999998pt height=65.69194499999999pt/></p>
and with stationary boundary conditions, i.e.
<p align="center"><img src="/sundials/svgs/7145dfb04640484f0e3ddfd91571774b.svg?invert_in_darkmode" align=middle width=159.92129999999997pt height=65.69194499999999pt/></p>
Here, we use a piecewise linear Galerkin finite element
discretization in space, where all element-wise integrals are
computed using 3-node Gaussian quadrature (since we will have
quartic polynomials in the reaction terms for the <img src="/sundials/svgs/e6897b8647f3bd38144535d3f40078e2.svg?invert_in_darkmode" align=middle width=14.322330000000001pt height=14.102549999999994pt/> and <img src="/sundials/svgs/3e3c6ee78813607a4d976d92c19dd36e.svg?invert_in_darkmode" align=middle width=12.885510000000002pt height=14.102549999999994pt/>
equations, including the test function).  The time derivative
terms for this system will include a mass matrix, giving rise
to an ODE system of the form
<p align="center"><img src="/sundials/svgs/feb6dcc9921528ee025c6b144635b061.svg?invert_in_darkmode" align=middle width=131.642775pt height=16.376943pt/></p>
where <img src="/sundials/svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.67348pt height=22.381919999999983pt/> is the block mass matrix for each component, <img src="/sundials/svgs/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode" align=middle width=11.145420000000001pt height=22.381919999999983pt/> is
the block Laplace operator for each component, and <img src="/sundials/svgs/4051c5cf4a2c287d8c463da35cb695a5.svg?invert_in_darkmode" align=middle width=33.91872pt height=24.56552999999997pt/> is
a 3x3 block comprised of the nonlinear reaction terms for
each component.  Since it it highly inefficient to rewrite
this system as
<p align="center"><img src="/sundials/svgs/7e24cd0234b30362b1a26b71378daf1f.svg?invert_in_darkmode" align=middle width=162.028845pt height=18.269295pt/></p>
we solve this system using ARKStep, with a user-supplied mass
matrix.  We therefore provide functions to evaluate the ODE RHS
<p align="center"><img src="/sundials/svgs/48173febf3555b003e1f02241d2be7e0.svg?invert_in_darkmode" align=middle width=144.49248pt height=16.376943pt/></p>
its Jacobian
<p align="center"><img src="/sundials/svgs/2ab039b6508ad5aa0e3a4083cb6f3b58.svg?invert_in_darkmode" align=middle width=127.93934999999999pt height=36.953894999999996pt/></p>
and the mass matrix, <img src="/sundials/svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.67348pt height=22.381919999999983pt/>.
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

```
SPDX-License-Identifier: BSD-3-Clause

LLNL-CODE-667205  (ARKODE)
UCRL-CODE-155951  (CVODE)
UCRL-CODE-155950  (CVODES)
UCRL-CODE-155952  (IDA)
UCRL-CODE-237203  (IDAS)
LLNL-CODE-665877  (KINSOL)
```
