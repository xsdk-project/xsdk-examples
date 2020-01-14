# SUNDIALS examples

Example codes demonstrating the use of SUNDIALS using other XSDK packages.


## SUNDIALS + SuperLU

The `ark_brusselator1D_FEM_sludist.cpp` example simulates a brusselator
problem from chemical kinetics.
This is a PDE system with 3 components, Y = [u,v,w], satisfying the equations,
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/6cedb29b4a78a44f92b4c86acd199b2c.svg?invert_in_darkmode" align=middle width=698.24205pt height=18.269295pt/></p>
for t in [0, 80], x in [0, 1], with initial conditions
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/cd1636e383e0d49f72e8dfed5df42de1.svg?invert_in_darkmode" align=middle width=645.7703999999999pt height=16.376943pt/></p>
and with stationary boundary conditions, i.e.
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/72de083d1e05cb26f5f6c7e5eb65a4e7.svg?invert_in_darkmode" align=middle width=457.90965pt height=16.376943pt/></p>
Here, we use a piecewise linear Galerkin finite element
discretization in space, where all element-wise integrals are
computed using 3-node Gaussian quadrature (since we will have
quartic polynomials in the reaction terms for the u_t and v_t
equations, including the test function).  The time derivative
terms for this system will include a mass matrix, giving rise
to an ODE system of the form
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/feb6dcc9921528ee025c6b144635b061.svg?invert_in_darkmode" align=middle width=131.642775pt height=16.376943pt/></p>
where M is the block mass matrix for each component, L is
the block Laplace operator for each component, and R(y) is
a 3x3 block comprised of the nonlinear reaction terms for
each component.  Since it it highly inefficient to rewrite
this system as
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/7e24cd0234b30362b1a26b71378daf1f.svg?invert_in_darkmode" align=middle width=162.028845pt height=18.269295pt/></p>
we solve this system using ARKStep, with a user-supplied mass
matrix.  We therefore provide functions to evaluate the ODE RHS
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/48173febf3555b003e1f02241d2be7e0.svg?invert_in_darkmode" align=middle width=144.49248pt height=16.376943pt/></p>
its Jacobian
<p align="center"><img src="https://rawgit.com/in	git@github.com:xsdk-project/xsdk-examples/master/svgs/af613be3db3346dc98a03b655d73560a.svg?invert_in_darkmode" align=middle width=149.24712pt height=16.376943pt/></p>
and the mass matrix, M.
This program solves the problem with the DIRK method, using a
Newton iteration with the SuperLU_DIST SUNLinearSolver.
100 outputs are printed at equal time intervals, and run
statistics are printed at the end.