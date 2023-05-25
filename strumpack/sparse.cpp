/*
 * STRUMPACK -- STRUctured Matrix PACKage, Copyright (c) 2014-2022,
 * The Regents of the University of California, through Lawrence
 * Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov.
 *
 * NOTICE. This software is owned by the U.S. Department of Energy. As
 * such, the U.S. Government has been granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable,
 * worldwide license in the Software to reproduce, prepare derivative
 * works, and perform publicly and display publicly.  Beginning five
 * (5) years after the date permission to assert copyright is obtained
 * from the U.S. Department of Energy, and subject to any subsequent
 * five (5) year renewals, the U.S. Government is granted for itself
 * and others acting on its behalf a paid-up, nonexclusive,
 * irrevocable, worldwide license in the Software to reproduce,
 * prepare derivative works, distribute copies to the public, perform
 * publicly and display publicly, and to permit others to do so.
 *
 * Developers: Pieter Ghysels, and others.
 *             (Lawrence Berkeley National Lab, Computational Research
 *             Division).
 *
 */
#include <iostream>
#include "StrumpackSparseSolverMPIDist.hpp"
#include "sparse/CSRMatrix.hpp"
#include "misc/TaskTimer.hpp"

typedef double scalar;
// typedef int64_t integer;  // to use 64 bit integers
typedef int integer;

using namespace strumpack;

int main(int argc, char* argv[]) {
  int thread_level, myrank;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_level);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (thread_level != MPI_THREAD_MULTIPLE && myrank == 0)
    std::cout << "MPI implementation does not support MPI_THREAD_MULTIPLE"
              << std::endl;
  {
    int n = 30, nrhs = 1;
    if (argc > 1) n = atoi(argv[1]); // get grid size
    else std::cout << "# please provide grid size" << std::endl;
    // get number of right-hand sides
    // if (argc > 2) nrhs = std::max(1, atoi(argv[2]));
    if (!myrank)
      std::cout << "solving 3D " << n << "^3 Poisson problem"
                << " with " << nrhs << " right hand sides" << std::endl;

    // Create the main solver object, using an MPI communicator.
    StrumpackSparseSolverMPIDist<scalar,integer> spss(MPI_COMM_WORLD);

    // The matching phase finds a column permutation that maximizes
    // the diagonal elements. Since the 3D Poisson problem is already
    // diagonally dominant, we can disable this matching.
    spss.options().set_matching(MatchingJob::NONE);

    // A fill reducing ordering ordering is a symmtric permutation of
    // the matrix that minimizes the fill in the sparse triangular
    // factors. Since the problem here is defined on a regular mesh,
    // we can use a simple geometric nested dissection algorithm (see
    // also below spss.reorder(n, n, n) where the mesh dimensions are
    // specified). For general sparse matrices, use any other option,
    // such as the default ReorderingStrategy::METIS, or the parallel
    // ReorderingStrategy::PARMETIS (if STRUMPACK was configured with
    // PARMETIS support).
    spss.options().set_reordering_method(ReorderingStrategy::GEOMETRIC);
    spss.options().set_from_command_line(argc, argv);

    // construct a sparse matrix from a simple 7 point stencil, with
    // zero boundary conditions
    CSRMatrix<scalar,integer> A;
    if (!myrank) {
      int n2 = n * n;
      int N = n * n2;
      int nnz = 7 * N - 6 * n2;
      A = CSRMatrix<scalar,integer>(N, nnz);
      integer* col_ptr = A.ptr();
      integer* row_ind = A.ind();
      scalar* val = A.val();

      nnz = 0;
      col_ptr[0] = 0;
      for (integer xdim=0; xdim<n; xdim++)
        for (integer ydim=0; ydim<n; ydim++)
          for (integer zdim=0; zdim<n; zdim++) {
            integer ind = zdim+ydim*n+xdim*n2;
            val[nnz] = 6.0;
            row_ind[nnz++] = ind;
            if (zdim > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-1; } // left
            if (zdim < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+1; } // right
            if (ydim > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-n; } // front
            if (ydim < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+n; } // back
            if (xdim > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-n2; } // up
            if (xdim < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+n2; } // down
            col_ptr[ind+1] = nnz;
          }
      A.set_symm_sparse();
    }
    // This scatters the sparse matrix A from the root over all the
    // ranks, using a 1D block row distribution, see
    // https://portal.nersc.gov/project/sparse/strumpack/master/sparse_example_usage.html#autotoc_md9
    CSRMatrixMPI<scalar,integer> Adist(&A, MPI_COMM_WORLD, true);
    // delete sequential sparse matrix (on the root)
    A = CSRMatrix<scalar,integer>();

    auto n_local = Adist.local_rows();
    DenseMatrix<scalar> b(n_local, nrhs), x(n_local, nrhs),
      x_exact(n_local, nrhs);

    // construct a random exact solution
    x_exact.random();

    // compute a right hand-side corresponding to exact solution
    // x_exact
    Adist.spmv(x_exact, b);


    // One can also directly pass the CSR rowptr, colind, and value
    //   arrays, see
    //   https://portal.nersc.gov/project/sparse/strumpack/master/classstrumpack_1_1SparseSolver.html
    spss.set_matrix(Adist);

    // For geometric nested dissection, the the mesh dimensions n x n
    // x n (and separator width if not 1) need to be provided. For
    // other fill-reducing orderings, such as the default
    // ReorderingStrategy::GEOMETRIC, just call spss.reorder();
    spss.reorder(n, n, n);

    // the actual numerical factorization phase. If reorder() was not
    // already called, it will be called by factor internally.
    spss.factor();

    // solve a linear system Ax=b for x. If factor was not already
    // called, then it will be called by solve internally.
    spss.solve(b, x);

    auto scaled_res = Adist.max_scaled_residual(x, b);
    x.scaled_add(-1., x_exact);
    auto relerr = x.normF() / x_exact.normF();
    if (!myrank) {
      std::cout << "# COMPONENTWISE SCALED RESIDUAL = "
                << scaled_res << std::endl;
      std::cout << "# relative error = ||x-x_exact||_F/||x_exact||_F = "
                << relerr << std::endl;
    }
  }
  scalapack::Cblacs_exit(1);
  MPI_Finalize();
  return 0;
}
