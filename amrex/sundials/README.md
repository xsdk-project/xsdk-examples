# AMReX + SUNDIALS examples

Example codes demonstrating the use of [AMReX](https://amrex-codes.github.io/)
and [SUNDIALS](https://computing.llnl.gov/projects/sundials).

## Advection-Diffusion Example

This is an example of a scalar-valued advection-diffusion problem for chemical
transport. The governing equation is:

$$ \frac{\partial u}{\partial t} + \vec{a} \cdot \nabla u -  \nabla \cdot ( D \nabla u ) = 0 $$

where $u = u(t,x,y)$ is the chemical concentration, $\vec{a}$ is the advection
vector, and $D$ is a diagonal matrix of diffusion coefficients. The problem is
solved on the square domain $[-1, 1]^2$ with periodic boundary conditions and
the initial condition is

$$ u(0,x,y) = u_0(x,y) = \frac{10}{\sqrt{2\pi}} e^{-50 (x^2 + y^2)} $$

## Problem Options

The problem inputs are listed below and may be specified on the command line
e.g., `./amrex_sundials_advection_diffusion help=1` or by supplying an input
file of values e.g., `./amrex_sundials_advection_diffusion inputs` where
`inputs` is a text file with `option = value` lines.

| Option               | Type   | Description                                        | Default  |
|----------------------|--------|----------------------------------------------------|----------|
| `help`               | `int`  | print input options and exit (1)                   | 0        |
| `n_cell`             | `int`  | number of cells on each side of the square domain  | 128      |
| `max_grid_size`      | `int`  | max size of boxes in box array                     | 64       |
| `plot_int`           | `int`  | enable (1) or disable (0) plots                    | 0        |
| `arkode_order`       | `int`  | ARKStep method order                               | 4        |
| `nls_max_iter`       | `int`  | maximum number of nonlinear iterations             | 3        |
| `ls_max_iter`        | `int`  | maximum number of linear iterations                | 5        |
| `rhs_adv`            | `int`  | advection: implicit (0) or explicit (1)            | 1        |
| `rtol`               | `Real` | relative tolerance                                 | 1e-4     |
| `atol`               | `Real` | absolute tolerance                                 | 1e-9     |
| `fixed_dt`           | `Real` | use a fixed time step size (if `fixed_dt` > 0.0)   | -1.0     |
| `tfinal`             | `Real` | final integration time                             | 1e4      |
| `dtout`              | `Real` | output frequency                                   | `tfinal` |
| `max_steps`          | `int`  | maximum number of steps between outputs            | 10000    |
| `advCoeffx`          | `Real` | advection speed in the x-direction                 | 5.0e-4   |
| `advCoeffy`          | `Real` | advection speed in the y-direction                 | 2.5e-4   |
| `diffCoeffx`         | `Real` | diffusion coefficient in the x-direction           | 1.0e-6   |
| `diffCoeffy`         | `Real` | diffusion coefficient in the y-direction           | 1.0e-6   |
| `use_preconditioner` | `int`  | use preconditioning (1) or not (0)                 | 0        |

If preconditioning is enabled, then additional options may be set (see AMReX
documentation of the `MLMG` solver for descriptions):

| Option                      | Type   | Default |
|-----------------------------|--------|---------|
| `mlmg.agglomeration`        | `int`  | 1       |
| `mlmg.consolidation`        | `int`  | 1       |
| `mlmg.max_coarsening_level` | `int`  | 1000    |
| `mlmg.linop_maxorder`       | `int`  | 2       |
| `mlmg.max_iter`             | `int`  | 1000    |
| `mlmg.max_fmg_iter`         | `int`  | 1000    |
| `mlmg.mg_fixed_iter`        | `int`  | 0       |
| `mlmg.verbose`              | `int`  | 0       |
| `mlmg.bottom_verbose`       | `int`  | 0       |
| `mlmg.use_hypre`            | `int`  | 1       |
| `mlmg.hypre_interface`      | `int`  | 3       |
| `mlmg.use_petsc`            | `int`  | 0       |
| `mlmg.tol_rel`              | `Real` | 1.0e-6  |
| `mlmg.tol_abs`              | `Real` | 1.0e-6  |
