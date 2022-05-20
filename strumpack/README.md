<!--
STRUMPACK -- STRUctured Matrix PACKage, Copyright (c) 2014-2022, The
Regents of the University of California, through Lawrence Berkeley
National Laboratory (subject to receipt of any required approvals from
the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this
software, please contact Berkeley Lab's Technology Transfer Department
at TTD@lbl.gov.

NOTICE. This software is owned by the U.S. Department of Energy. As
such, the U.S. Government has been granted for itself and others
acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
license in the Software to reproduce, prepare derivative works, and
perform publicly and display publicly.  Beginning five (5) years after
the date permission to assert copyright is obtained from the
U.S. Department of Energy, and subject to any subsequent five (5) year
renewals, the U.S. Government is granted for itself and others acting
on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
in the Software to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so.
-->
# STRUMPACK examples

Example codes demonstrating the use of [STRUMPACK](https://portal.nersc.gov/project/sparse/strumpack/) using other XSDK packages ([ButterflyPACK](https://github.com/liuyangzhuan/ButterflyPACK)).

## STRUMPACK + ButterflyPACK

The `sparse.cpp` example demonstrates STRUMPACK's algebraic sparse
direct solver and preconditioners for solving the 3-D Laplacian
problem with zero boundary conditions on an n x n x n grid.  The
number of unknowns is N=n^3.  The standard 7-point stencil is used,
and we solve for the interior nodes only.

STRUMPACK implements multifrontal sparse LU factorization, with the
option of compression of the larger frontal matrices. Compression
algorithms include HODLR (Hierarchically Off-Diagonal Low Rank) and
HODBF (Hierarchically Off-Diagonal Butterfly), BLR (Block Low Rank),
HSS (Hierarchically Semi Separable), lossy or lossless.

After factorization, linear system can be solved using forward and
backward substitution with the lower and upper triangular factors
respectively. Without compression, the solver behaves as a sparse
direct method. The sparse direct solver still uses iterative
refinement, but typically only needs a single iteration.  When
compression is enabled, the LU factorization is only approximate, and
the solver is used as a preconditioner for GMRES (or BiCGStab).


### Usage

**Sample run**:   `OMP_NUM_THREADS=4 mpirun -np 4 ./sparse 100 --sp_compression hodlr --hodlr_butterfly_levels 10`

 - This run solves a system corresponding to a discretization
    of the Laplace equation -Delta u = f with zero boundary
    conditions on a 100 x 100 x 100 grid, using STRUMPACK's
    distributed LU factorization, and with both HODBF compression
    for the largest fronts (dense sub-blocks in the sparse factors).
 - The output of the example is various information regarding the
    solver and solver performance.

Options for the compression algorithm include `none` (direct solver),
`blr`, `lossy`/`lossless` (needs the ZFP library), `hodlr` (needs
ButterflyPACK), `blr_hodlr` (needs ButterflyPACK). `blr_hodlr`
combines both BLR (on medium sized blocks) and HODLR (on the largest
blocks).

The thresholds for using the compression schemes can be set using
`--sp_compression_min_sep_size 1000` or, for specific formats
``--sp_blr_min_sep_size 1000``.  The compression tolerance can be
tuned using `--blr_rel_tol 1e-6` or `--hodlr_rel_tol 1e-6`


To see the full list of command-line options, run: `./sparse --help`

For more information on how to tune the preconditioners, see
[here](https://portal.nersc.gov/project/sparse/strumpack/master/prec.html).

## License
STRUMPACK -- STRUctured Matrix PACKage, Copyright (c) 2014-2022, The
Regents of the University of California, through Lawrence Berkeley
National Laboratory (subject to receipt of any required approvals from
the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this
software, please contact Berkeley Lab's Technology Transfer Department
at TTD@lbl.gov.

NOTICE. This software is owned by the U.S. Department of Energy. As
such, the U.S. Government has been granted for itself and others
acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
license in the Software to reproduce, prepare derivative works, and
perform publicly and display publicly.  Beginning five (5) years after
the date permission to assert copyright is obtained from the
U.S. Department of Energy, and subject to any subsequent five (5) year
renewals, the U.S. Government is granted for itself and others acting
on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
in the Software to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so
