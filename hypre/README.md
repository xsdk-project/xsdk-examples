<!--
Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
HYPRE Project Developers. See the top-level COPYRIGHT file for details.

SPDX-License-Identifier: (Apache-2.0 OR MIT)
-->
# HYPRE examples

Example codes demonstrating the use of [HYPRE](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) using other XSDK packages (SuperLU_DIST).

## HYPRE + SuperLU_DIST

The `ij_laplacian.c` example demonstrates HYPRE's Linear-Algebraic (IJ) interface
for solving the 2-D Laplacian problem with zero boundary conditions on an n x n grid.
The number of unknowns is N=n^2.
The standard 5-point stencil is used, and we solve for the interior nodes only.

Available solvers are AMG (default), PCG, PCG with AMG or Parasails
preconditioners, flexible GMRES with AMG. 
Within HYPRE's AMG V-cycle there is an option to use a 
distributed LU factorization from
SuperLU_DIST as the coarse level solver.

### Usage

**Sample run**:   `mpirun -np 4 ij_laplacian -dslu_th 50`

 - This run solves a system corresponding to a discretization 
    of the Laplace equation -Delta u = 1 with zero boundary
    conditions on a 33 x 33 grid, using HYPRE's BoomerAMG with SuperLU_DIST's
    distributed LU factorization for the coarse level solver where `50` 
    is a maximum number of coarse level degrees of freedom using 4 MPI tasks. 
 - The output of the example is various information regarding the  
    solver and solver performance. 

By specifying the command line parameter `-dslu_th` to be
the maximum coarse level number of degrees of freedom, the
coarse level solver within BoomerAMG will be changed from the
default Gaussian Elimination, to a sparse LU decomposition and
triangular solution from SuperLU_DIST.

Optionally running with `-vis` command-line option, the solution is saved for GLVis visualization (see `vis/glvis-ij-laplacian.sh` for example usage of GLVis visualization).

To see the full list of command-line options, run: `ij_laplacian -help`

## License
HYPRE is distributed under the terms of both the MIT license and the Apache
License (Version 2.0). Users may choose either license, at their option.

All new contributions must be made under both the MIT and Apache-2.0 licenses.

See [LICENSE-MIT](./LICENSE-MIT), [LICENSE-APACHE](./LICENSE-APACHE),
[COPYRIGHT](./COPYRIGHT), and [NOTICE](./NOTICE) for details.

SPDX-License-Identifier: (Apache-2.0 OR MIT)

LLNL-CODE-778117
