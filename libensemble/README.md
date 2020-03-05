A libEnsemble example that intferfaces with the PETSc/TAO blmvm solver.

This is example works with libEnsemble v0.5.2 (As used in xSDK release v0.5.0).

The example finds minima for a 6-hump camel function using the limited memory
variable metric method for bound constraind minimization.

Execute via one of the following commands (e.g. 3 workers).

Using MPI (mpi4py):

    mpiexec -np 4 python3 test_persistent_aposmm_tao.py

Using multiprocessing:

    python3 test_persistent_aposmm_tao.py --nworkers 3 --comms local
