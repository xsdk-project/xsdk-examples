# AMReX + SUNDIALS examples

Example codes demonstrating the use of [AMReX](https://amrex-codes.github.io/)
and [SUNDIALS](https://computing.llnl.gov/projects/sundials).

## Advection-Diffusion Example

This is an example of a scalar-valued advection-diffusion problem for chemical
transport. The governing equation is:

u_t + a \cdot \nabla u - \nabla \cdot ( D \nabla u ) = 0

where u(t,x,y) is the chemical concentration, a is the advection vector, and D
is diagonal matrix containing anisotropic diffusion coefficients. The problem is
solved on the unit square domain centered at the origin and evolved from time 0
to 10^3.

### Problem Options

The problem inputs are listed below and may be specified on the command line
e.g., `./amrex_sundials_advection_diffusion help=1` or by supplying an input
file of values e.g., `./amrex_sundials_advection_diffusion inputs` where
`inputs` is a text file with `option = value` lines.

| Option               | Type   | Description                                        | Default  |
|:---------------------|:-------|:---------------------------------------------------|:---------|
| `n_cell`             | `int`  | number of cells on each side of the square domain  | 128      |
| `max_grid_size`      | `int`  | max size of boxes in box array                     | 64       |
| `plot_int`           | `int`  | enable (1) or disable (0) plots                    | -1       |
| `arkode_order`       | `int`  | ARKStep method order                               | 4        |
| `nls_max_iter`       | `ìnt`  | maximum number of nonlinear iterations             | 3        |
| `ls_max_iter`        | `int`  | maximum number of linear iterations                | 5        |
| `rhs_adv`            | `int`  | advection: disable (0), implicit (1), explicit (2) | 2        |
| `rhs_diff`           | `int`  | diffusion: disable (0), implicit (1), explicit (2) | 1        |
| `rtol`               | `Real` | relative tolerance                                 | 1e-4     |
| `atol`               | `Real` | absolute tolerance                                 | 1e-9     |
| `fixed_dt`           | `Real` | use a fixed time step size (if `fixed_dt` > 0.0)   | -1.0     |
| `tfinal`             | `Real` | final integration time                             | 1e3      |
| `dtout`              | `Real` | output frequency                                   | `tfinal` |
| `max_steps`          | `int`  | maximum number of steps between outputs            | 10000    |
| `advCoeffx`          | `Real` | advection speed in the x-direction                 | 5e-4     |
| `advCoeffy`          | `Real` | advection speed in the y-direction                 | 5e-4     |
| `diffCoeffx`         | `Real` | diffusion coefficient in the x-direction           | 2e-5     |
| `diffCoeffy`         | `Real` | diffusion coefficient in the y-direction           | 2e-5     |
| `use_preconditioner` | `int`  | use preconditioning (1) or not (0)                 | 0        |

If preconditioning is enabled, then additional options may be set (see AMReX
documentation of the `MLMG` solver for descriptions):

| Option                    | Type   | Default |
|:--------------------------|:-------|:--------|
| mlmg.agglomeration        | `int`  | 1       |
| mlmg.consolidation        | `int`  | 1       |
| mlmg.max_coarsening_level | `int`  | 1000    |
| mlmg.linop_maxorder       | `int`  | 2       |
| mlmg.max_iter             | `int`  | 1000    |
| mlmg.max_fmg_iter         | `int`  | 1000    |
| mlmg.verbose              | `int`  | 0       |
| mlmg.bottom_verbose       | `int`  | 0       |
| mlmg.use_hypre            | `int`  | 1       |
| mlmg.hypre_interface      | `int`  | 3       |
| mlmg.use_petsc            | `int`  | 0       |
| mlmg.tol_rel              | `Real` | 1.0e-6  |
