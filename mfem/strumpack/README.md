# MFEM-STRUMPACK-SuperLU-HYPRE example

This example code utilizes MFEM to discretize a diffusion eigenvalue problem
-Delta u = lambda u, with homogeneous Dirichlet boundary conditions. The problem
is discretized using continuous finite elements of arbitrary specified order on
any given mesh. The discretized eigenvalue problem is solved using HYPRE's
LOBPCG eigenvalue solver with linear system preconditioner/solver using
STRUMPACK, SuperLU, or HYPRE's BoomerAMG.

This example is built to run in parallel, so launch it with mpirun and your
desired options, e.g. using STRUMPACK:
```
mpirun -np 4 ./diffusion-eigen -m ../data/star.mesh --strumpack
```
SuperLU:
```
mpirun -np 4 ./diffusion-eigen -m ../data/star.mesh --superlu
```
or HYPRE:
```
mpirun -np 4 ./diffusion-eigen -m ../data/star.mesh --no-strumpack
```

For a full list of options, see
```
./diffusion-eigen -h
```
Note that when using STRUMPACK (the default) command-line parameters are also
passed to STRUMPACK to support STRUMPACK-specifiic parameters.
