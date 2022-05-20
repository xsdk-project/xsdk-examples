# MFEM-HYPRE example

This example demonstrates the integration of MFEM with HYPRE on both CPUs and
GPUs.

This example code utilizes MFEM to discretize a steady state magnetic diffusion
problem, curl curl E + E = f, with suitable boundary condition. The problem is
discretized using H(curl)-conforming finite elements of arbitrary specified
order on any given mesh. The discretized problem is then solved using PCG with
HYPRE's AMS preconditioner.

Sample CPU runs:
```
mpirun -np 4 ./magnetic-diffusion -m ../data/star.mesh
mpirun -np 4 ./magnetic-diffusion -m ../data/beam-hex.mesh
mpirun -np 4 ./magnetic-diffusion -m ../data/beam-hex.mesh -o 2 -pa
```

Sample GPU runs, replace `<dev>` with `cuda` or `hip`:
```
mpirun -np 4 ./magnetic-diffusion -m ../data/star.mesh -pa -d <dev>
mpirun -np 4 ./magnetic-diffusion -m ../data/star.mesh -no-pa -d <dev>
mpirun -np 4 ./magnetic-diffusion -m ../data/beam-hex.mesh -pa -d <dev>
```

For a full list of options, see
```
./magnetic-diffusion -h
```
