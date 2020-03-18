# MFEM-hypre-superlu example

This example code utilizes MFEM to discretize a Convection-Diffusion equation.  This problem is 
models the steady-state concentration of a material as it diffuses and is carried through a flowing 
medium.  From the command line the user can vary the velocity of the flow field.  The linear system 
that arises from this discretization can optionally be solved utilizing HYPRE's
BoomerAMG solver or the SuperLU solver.

This example can be used to show the advantages and disadvantages of HYPRE's iterative solvers,
and SuperLU's direct solvers.  At lower velocities the HYPRE solver/preconditioner is 
sufficient to solve the system and does so very quickly.  At higher velocities the BoomerAMG solver
is no longer able to handle the increasingly ill-conditioned system and a direct solver is
required.  For these high velocities selecting the SuperLU solver will allow the system to 
be solved.

This example is built to run in parallel, so launch it with mpirun with your desired options:
mpirun -np 4 convdiff --velocity 100.0 --no-superlu

Useful non-default options:
|   Flag                | Meaning                                               |
|:----------------------| :-----------------------------------------------------|
| --refine n            | Number of times to uniformly refine the initial mesh. |
| --order n             | Set the polynomial order of the discretization.       |
| --velocity n          | Velocity of the flow field.                           |
| --superlu             | Use the SuperLU direct solver.                        |
| --no-superlu          | Use the interative HYPRE BoomerAMG solver.            |
| --visit               | Output VisIt files for visualation of the solution.   |
