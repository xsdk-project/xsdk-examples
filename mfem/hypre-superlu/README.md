# MFEM-hypre-superlu example

This example code utilizes MFEM to discretize a Convection-Diffusion equation.  This problem is 
models the steady-state concentration of a material as it diffuses and is carried through a flowing 
medium.  From the command line the user can vary the velocity of the flow field.  The linear system 
that arises from this discretization can optionally be solved utilizing HYPRE's
BoomerAMG solver or the SuperLU_Dist solver.

This example can be used to show the advantages and disadvantages of HYPRE's iterative solvers,
and SuperLU's direct solvers.  At lower velocities the HYPRE solver/preconditioner is 
sufficient to solve the system and does so very quickly.  At higher velocities the BoomerAMG solver
is no longer able to handle the increasingly ill-conditioned system and a direct solver is
required.  For these high velocities selecting the SuperLU solver will allow the system to 
be solved.

This example is built to run in parallel, so launch it with mpirun and your desired options:
```
mpirun -np 4 ./convdiff --velocity 100.0 --no-superlu
mpirun -np 4 ./convdiff --velocity 100.0 --superlu
```

Available options:
|   Flag                | Meaning                                               |
|:----------------------| :-----------------------------------------------------|
| --help                | Print a help message and exit.                        |
| --refine n            | Number of times to uniformly refine the initial mesh. |
| --order n             | Set the polynomial order of the discretization.       |
| --velocity n          | Velocity of the flow field.                           |
| --superlu             | Use the SuperLU_Disy direct solver.                   |
| --no-superlu          | Use the GMRES + HYPRE BoomerAMG solver. (default)     |
| --slu-colperm n       | Set the SuperLU Column permutation algorithm:         |
|                       |    0 - natural                                        |
|                       |    1 - mmd-ata                                        |
|                       |    2 - mmd_at_plus_a                                  |
|                       |    3 - colamd                                         |
|                       |    4 - metis_at_plus_a (default)                      |
|                       |    5 - parmetis                                       |
|                       |    6 - zoltan                                         |
| --slu-rowperm n       | Set the SuperLU Row permutation algorithm:            |
|                       |    0 - NoRowPerm                                      |
|                       |    1 - LargeDiag (default)                            |
|                       |    2 - MyPermR                                        |
| --one-matrix          | Solve with one matrix. (default)                      |
| --two-matrix          | Solve with two different matrices.                    |
| --one-matrix          | Solve with one rhs. (default)                         |
| --two-matrix          | Solve with two different rhs.                         |
| --visit               | Output VisIt files for visualation of the solution.   |
| --no-visit            | Do not output VisIt files. (default)                  |
