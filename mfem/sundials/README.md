# MFEM-SUNDIALS examples

## transient-heat.cpp

This example solves a time dependent heat equation with a temperature dependent
conductivity. MFEM is used for the spatial discretization and SUNDIALS is 
used for the ODE time integration.  This problem demonstrates the MFEM 
integrations with the SUNDIALS CVODE and ARKODE solvers which give MFEM
access to an array of advance ODE solvers.

This example is built to run in parallel, so launch it with mpirun and your desired options:
```
mpirun -np 4 ./transient-heat --kappa 0.5 --alpha 0.01 --ode-solver 8
```

Useful non-default options:
|   Flag                | Meaning                                               |
|:----------------------| :-----------------------------------------------------|
| --order n             | Set the polynomial order of the discretization.       |
| --kappa n             | Set up the conductivity model C(T) = kappa + alpha T. |
| --alpha n             | Set up the conductivity model C(T) = kappa + alpha T. |
| --ode-solver n        | Pick the ODE solver used for the time integration.    |
| ..............        | 1  - MFEM (Forward Euler).                            |
| ..............        | 2  - MFEM (RK2).                                      |
| ..............        | 3  - MFEM (RK3 SSP).                                  |
| ..............        | 4  - MFEM (RK4).                                      |
| ..............        | 5  - MFEM (Backward Euler).                           |
| ..............        | 6  - MFEM (SDIRK 2).                                  |
| ..............        | 7  - MFEM (SDIRK 3).                                  |
| ..............        | 8  - CVODE (implicit Adams).                          |
| ..............        | 9  - CVODE (implicit BDF).                            |
| ..............        | 10 - ARKODE (default explicit).                       |
| ..............        | 11 - ARKODE (explicit Fehlberg-6-4-5).                |
| ..............        | 12 - ARKODE (default impicit).                        |
| --mesh s              | Mesh file to use.                                     |


## advection.cpp

This example code solves the time-dependent advection equation. MFEM is used for the
spatial discretization and SUNDIALS is used for the ODE time integration. This problem
demonstrates MFEM integration with the SUNDIALS CVODE and ARKODE solvers for CUDA.

This example is built to run in parallel, so launch it with mpirun and your desired options:
```
mpirun -np 4 ./advection --device cuda --partial-assembly --ode-solver 9
mpirun -np 4 ./advection --device cpu --partial-assembly --ode-solver 9
```

Useful non-default options:
|   Flag                | Meaning                                               |
|:----------------------| :-----------------------------------------------------|
| --device s            | Device configuration string. Choose cpu or cuda.      |
| --order n             | Set the polynomial order of the discretization.       |
| --partial-assembly    | Enable Partial Assembly.                              |
| --element-assembly    | Enable Element Assembly.                              |
| --full-assembly       | Enable Full Assembly.                                 |
| --ode-solver n        | Pick the ODE solver used for the time integration.    |
| ..............        | 1  - MFEM (Forward Euler).                            |
| ..............        | 2  - MFEM (RK2).                                      |
| ..............        | 3  - MFEM (RK3 SSP).                                  |
| ..............        | 4  - MFEM (RK4).                                      |
| ..............        | 6  - MFEM (RK6).                                      |
| ..............        | 7  - CVODE (implicit Adams).                          |
| ..............        | 8  - ARKODE (default 4th order explicit).             |
| ..............        | 9  - ARKODE (RK8).                                    |
| --mesh s              | Mesh file to use.                                     |
