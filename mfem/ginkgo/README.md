# MFEM-Ginkgo Example

This example modifies MFEM's `example 22` to demonstrate the use of Ginkgo from
within MFEM.  The code can solve three variations of a damped harmonic oscillator.
MFEM handles the discretization, linear operator (assembled as a `SparseMatrix` or 
applied directly by MFEM within Ginkgo as a matrix-free operator), right hand side 
vector, and preconditioner; Ginkgo's GMRES solver is used in place of MFEM's.



Useful non-default options:
|   Flag                | Meaning                                           |
|:----------------------| :-------------------------------------------------|
|  -m [file]            | Mesh file to use.                                 |
| -d cuda               | Use the MFEM cuda backend and Ginkgo CudaExecutor |
| -d hip                | Use the MFEM hip backend and Ginkgo HipExecutor   |
| -no-pa,  --no-partial-assembly | Don't use partial assembly, but build the full SparseMatrix for the system  |
| -no-gko, --no-ginkgo  | Don't use Ginkgo (for easy comparison with "all-MFEM" case) |

## More about the MFEM-Ginkgo integration

Modifying MFEM code to use Ginkgo via MFEM's wrappers is simple.  First, we create a 
`GinkgoExecutor` wrapper (for Ginkgo's `Executor` class) to match the device configuration of MFEM:

```
 Ginkgo::GinkgoExecutor exec(device); //device is the MFEM Device configuration
```

The Ginkgo wrappers provided in MFEM provide mix-and-match interoperability: 
Ginkgo solvers can use MFEM preconditioners, and MFEM can use Ginkgo preconditioners.
To allow Ginkgo to use the MFEM block diagonal preconditioner, `BDP`, we create
a wrapped object:

```
Ginkgo::MFEMPreconditioner gko_precond(exec, BDP);
```

Then, we use this preconditioner when creating the Ginkgo GMRES solver:

```
Ginkgo::GMRESSolver gmres(exec, gko_precond, 50); // 50 is GMRES restart value
```

From this point, the configuration and application of the GMRES solver proceed with
code identical to that for the MFEM GMRES solver:

```
// Finish configuring the solver
gmres.SetOperator(*A.Ptr());
gmres.SetRelTol(1e-12);
gmres.SetMaxIter(1000);
gmres.SetPrintLevel(1);
// Apply the solver with rhs B and solution U
gmres.Mult(B, U);
```

The `Operator`, `A`, can be a `SparseMatrix` or a matrix-free operator. For `SparseMatrix`
and `Vector` objects, MFEM and Ginkgo operate on the same data, whether on host or device,
without unnecessary copying between the two libraries. 

An example of MFEM using a Ginkgo preconditioner can be found in MFEM's Ginkgo-enabled example 1,
(`mfem/examples/ginkgo/ex1.cpp`).
For more about Ginkgo, see the Ginkgo [documentation](https://ginkgo-project.github.io/ginkgo/doc/develop/). 
