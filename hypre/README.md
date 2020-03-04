<!--
Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
HYPRE Project Developers. See the top-level COPYRIGHT file for details.

SPDX-License-Identifier: (Apache-2.0 OR MIT)
-->
# HYPRE examples

Example codes demonstrating the use of [HYPRE](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) using other XSDK packages (SuperLU_DIST).

## HYPRE + SuperLU_DIST

The `ij_laplacian.c` example demonstrates HYPRE's Linear-Algebraic (IJ) interface
for solving the 2-D Laplacian problem with zero boundary conditions on an n x n grid
using various solvers (AMG, PCG, PCG with AMG preconditioner, Flexible GMRES with AMG
preconditioner) from HYPRE with the option to use a distributed LU factorization from
SuperLU_DIST as the coarse level solver in BoomerAMG.  The number of unknowns is N=n^2.
The standard 5-point stencil is used, and we solve for the interior nodes only.

Optionally, the solution may be saved for GLVis visualization (see 
`vis/glvis-ij-laplacian.sh`).

## License
HYPRE is distributed under the terms of both the MIT license and the Apache
License (Version 2.0). Users may choose either license, at their option.

All new contributions must be made under both the MIT and Apache-2.0 licenses.

See [LICENSE-MIT](./LICENSE-MIT), [LICENSE-APACHE](./LICENSE-APACHE),
[COPYRIGHT](./COPYRIGHT), and [NOTICE](./NOTICE) for details.

SPDX-License-Identifier: (Apache-2.0 OR MIT)

LLNL-CODE-778117
