# MFEM-Ginkgo Example

This example discretizes a Poisson problem with the MFEM library and uses the 
the CG solver from the Ginkgo library to solve the linear system.  To do 
this, it uses two special wrapper classes, for MFEM Operators and Vectors.  The 
MFEMOperatorWrapper class allows the use of MFEM Operators that are not built into
sparse matrices (e.g., MFEM partial assembly mode) with Ginkgo, provided the vectors
given to the solver are of MFEMVectorWrapper type.

The example also uses an MFEMOperatorWrapper to wrap a basic MFEM preconditioner for 
Ginkgo's use.

Currently, MFEM's CUDA backend can also be used, except for the full assembly + preconditioning
combination.

This example is built to run in serial, so launch it with your desired options:
```
./mfem-gko --no-partial-assembly
```

Useful non-default options:
|   Flag                | Meaning                                           |
|:----------------------| :-------------------------------------------------|
|  -m [file]            | Mesh file to use.                                 |
| -no-pa,  --no-partial-assembly | Don't use partial assembly, but build the full SparseMatrix for the system  |
| --no-partial-assembly | build a full SparseMatrix for the system          |
| -pc, --preconditioner | Use preconditioning                               |
| -d "cuda"             | Use the MFEM cuda backend and Ginkgo CudaExecutor.|


## More about the MFEMVectorWrapper class

The MFEMVectorWrapper class inherits from Ginkgo's `matrix::Dense<double>` class,
which is what Ginkgo uses for double-precision vectors.  It stores a 
`std::unique_ptr` to an MFEM Vector, and the private data of the Ginkgo Dense
object points to the MFEM Vector data.  Thus, the same data is being used 
through the Ginkgo Dense and MFEM Vector interfaces, and either library
can modify the data through its usual routines.  This also avoids data copying.

As an example, consider the MFEMVectorWrapper's constructor:

```
class MFEMVectorWrapper : public gko::matrix::Dense<double> {
public:
    MFEMVectorWrapper(std::shared_ptr<const gko::Executor> exec,
                      gko::size_type size, mfem::Vector *mfem_vec,
                      bool on_device = true, bool ownership = false)
        : gko::matrix::Dense<double>(
              exec, gko::dim<2>{size, 1},
              gko::Array<double>::view(exec, size,
                                       mfem_vec->ReadWrite(on_device)),  // Here we use Ginkgo's 
              1)                                                         // `view` option, which 
                                                                         // creates a non-owning Array
                                                                         // through data pointers only.
                                                                         //  MFEM's `ReadWrite` access 
                                                                         //  will give either host or 
                                                                         //  device pointer to the 
                                                                         //  Vector's data.
    {
        // This controls whether or not we want Ginkgo to own its MFEM Vector.
        // Normally, when we are wrapping an MFEM Vector created outside 
        // Ginkgo, we do not want ownership to be true. However, Ginkgo
        // creates its own temporary vectors as part of its solvers, and 
        // these will be owned (and deleted) by Ginkgo. 
        if (ownership) {
            using deleter = mfem_destroy<mfem::Vector>;
            mfem_vec_ = std::unique_ptr<mfem::Vector,
                                        std::function<void(mfem::Vector *)>>(
                mfem_vec, deleter{});
        } else {
            using deleter = gko::null_deleter<mfem::Vector>;
            mfem_vec_ = std::unique_ptr<mfem::Vector,
                                        std::function<void(mfem::Vector *)>>(
                mfem_vec, deleter{});
        }
    }
```

When we need to directly access the MFEM Vector object, we can call one of two 
functions:

```
// Return reference to MFEM Vector object
mfem::Vector &get_mfem_vec_ref() const { return *(this->mfem_vec_.get()); }

// Return const reference to MFEM Vector object
const mfem::Vector &get_mfem_vec_const_ref() const
{
    return const_cast<const mfem::Vector &>(*(this->mfem_vec_.get()));
}
```

## More about the MFEMOperatorWrapper class

The MFEMOperatorWrapper is a custom Ginkgo LinOp.  It stores a pointer (as private
member `mfem_oper_`) to an MFEM OperatorHandle, which can point to any MFEM Operator type.  

When creating custom LinOps in Ginkgo, the key requirement is the implementation
of the two `apply` functions:

```
void apply_impl(const gko::LinOp *b, gko::LinOp *x) const override;

void apply_impl(const gko::LinOp *alpha, const gko::LinOp *b,
                const gko::LinOp *beta, gko::LinOp *x) const override;
```

To see how this wrapper works, let's look at the simpler `apply_impl` in more detail,
from [mfem_wrapper.cpp](./mfem_wrapper.cpp):

```
void MFEMOperatorWrapper::apply_impl(const gko::LinOp *b, gko::LinOp *x) const
{

    // Cast to MFEMVectorWrapper; only accept this type for this impl
    const MFEMVectorWrapper *mfem_b = gko::as<const MFEMVectorWrapper>(b);
    MFEMVectorWrapper *mfem_x = gko::as<MFEMVectorWrapper>(x);

    // Access the MFEM Operator's Mult function and apply to the 
    //  MFEM Vectors stored in the MFEMVectorWrapper objects 
    this->mfem_oper_->Mult(mfem_b->get_mfem_vec_const_ref(),
                           mfem_x->get_mfem_vec_ref());
}
```

In other words, the combination of MFEMOperatorWrapper and MFEMVectorWrapper allows 
Ginkgo to let MFEM handle the application of its operators without needing a matrix
(or to really know anything about what the MFEM Operator is doing in its `Mult()` function).

## Using the wrappers to solve MFEM example 1

There are three main sections where Ginkgo-related code has been added to MFEM's [ex1.cpp](https://github.com/mfem/mfem/blob/master/examples/ex1.cpp).

The first is to create a custom Ginkgo logger to store and output the residuals during
the CG iterations.  This code is modified from Ginkgo's `custom-logger` [example](https://github.com/ginkgo-project/ginkgo/tree/develop/examples/custom-logger).

The second is to create Ginkgo Executors as necessary for the `cpu` or `cuda` backend options.

The third section replaces MFEM's solver call with a Ginkgo solver. First, we wrap the 
right hand side `B` and solution vector `X` for Ginkgo:

```
// Create MFEM Vector wrappers.
auto gko_rhs = MFEMVectorWrapper::create(executor, B.Size(), &B, on_device);
auto gko_x = MFEMVectorWrapper::create(executor, X.Size(), &X, on_device);
```

We also need to wrap the MFEM Operator `A`:

```
auto oper_wrapper = MFEMOperatorWrapper::create(executor, B.Size(), A);
```

After some more set-up for the Ginkgo solver, we generate a solver for our 
operator and apply it:

```
// Generate CG solver for our specific wrapper MFEM Operator
auto solver = solver_gen->generate(gko::give(oper_wrapper));
    
// Solve system
solver->apply(gko::lend(gko_rhs), gko::lend(gko_x));
```

The `gko::give` and `gko::lend` functions are part of the smart pointer
management used by Ginkgo.  For more about this and the building of 
Ginkgo solver factories, see the Ginkgo [documentation](https://ginkgo-project.github.io/ginkgo/doc/develop/). 
