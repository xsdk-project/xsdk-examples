## libEnsemble with PETSc/TAO

This example script demonstrates interoperability between libEnsemble and
PETSc/TAO; another xSDK package.

#### Note: The supplied scripts work with libEnsemble v0.5.2 (As used in xSDK release v0.5.0).

This script uses libEnsemble to optimize a two-dimensional algebraic function, the [six-hump camel](https://www.sfu.ca/~ssurjano/camel6.html) function, over a bound-constrained domain.

The optimization is governed by the [APOSMM](https://www.mcs.anl.gov/~jlarson/APOSMM/) generating function, which performs multistart local optimization.

The local optimization instances are governed by a quasi-Newton method in
PETSc/TAO, [BLMVM](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Tao/TAOBLMVM.html).

### Usage

Execute via one of the following commands (e.g. 3 workers).

Using MPI (mpi4py):

    mpiexec -np 4 python3 test_persistent_aposmm_tao.py

Using multiprocessing:

    python3 test_persistent_aposmm_tao.py --nworkers 3 --comms local
