//                          MFEM Example 1, modified
//
// This code has been modified from `ex1.cpp` provided in the examples of 
// MFEM.  The sections not marked as related to Ginkgo are largely unchanged
// from the version provided by MFEM.  
//
// The linear system is now solved using the CG solver provided by 
// the Ginkgo library.  It requires the use of the wrapper functions in
// `mfem_wrapper.hpp`.  Together, the wrappers and this modified example
// show how a user can combine an MFEM preconditioner with a Ginkgo solver,
// even without building a full sparse matrix for the system (partial 
// assembly mode).
//
// The default mesh option is "beam-hex.mesh", provided by MFEM.  
// Important non-default options:
//  -m [file] : Mesh file.  
//  -no-pa, --no-partial-assembly : Don't use partial assembly, build a full
//                                  SparseMatrix for the system
//  -pc, --preconditioner : Use preconditioning.  
//  -d "cuda" : Use the MFEM cuda backend and Ginkgo CudaExecutor.
// 
//  NOTE: Currently the combination of "-no-pa -pc -d cuda" is not allowed.
//
// MFEM's provided information about `ex1.cpp`:
// Description:  This example code demonstrates the use of MFEM to define a
//               simple finite element discretization of the Laplace problem
//               -Delta u = 1 with homogeneous Dirichlet boundary conditions.
//               Specifically, we discretize using a FE space of the specified
//               order, or if order < 1 using an isoparametric/isogeometric
//               space (i.e. quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)
//
//               The example highlights the use of mesh refinement, finite
//               element grid functions, as well as linear and bilinear forms
//               corresponding to the left-hand side and right-hand side of the
//               discrete linear system. We also cover the explicit elimination
//               of essential boundary conditions, static condensation, and the
//               optional connection to the GLVis tool for visualization.

#include "mfem.hpp"

#include "mfem_wrapper.hpp"

#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// -------------------------------------------------------------------
// -- Start Ginkgo custom logger definition and auxiliary functions --
//
// Adapted from the Ginkgo `custom-logger` example, see 
//   `ginkgo/examples/custom-logger/custom-logger.cpp`

// Utility function which gets the scalar value of a Ginkgo gko::matrix::Dense
// matrix representing the norm of a vector.
template <typename ValueType>
double get_norm(const gko::matrix::Dense<ValueType> *norm)
{
    // Put the value on CPU thanks to the master executor
    auto cpu_norm = clone(norm->get_executor()->get_master(), norm);
    // Return the scalar value contained at position (0, 0)
    return cpu_norm->at(0, 0);
}

template <typename ValueType>
double compute_norm(const gko::matrix::Dense<ValueType> *b)
{
    // Get the executor of the vector
    auto exec = b->get_executor();
    // Initialize a result scalar containing the value 0.0.
    auto b_norm = gko::initialize<gko::matrix::Dense<ValueType>>({0.0}, exec);
    // Use the dense `compute_norm2` function to compute the norm.
    b->compute_norm2(lend(b_norm));
    // Use the other utility function to return the norm contained in `b_norm``
    return get_norm(lend(b_norm));
}


template <typename ValueType>
struct ResidualLogger : gko::log::Logger {

    using gko_dense = gko::matrix::Dense<ValueType>;
    // Customize the logging hook which is called everytime an iteration is
    // completed
    void on_iteration_complete(const gko::LinOp *,
                               const gko::size_type &iteration,
                               const gko::LinOp *residual,
                               const gko::LinOp *solution,
                               const gko::LinOp *residual_norm) const override
    {

        // If the solver shares a residual norm, log its value
        if (residual_norm) {
            auto dense_norm = gko::as<gko_dense>(residual_norm);
            // Add the norm to the `recurrent_norms` vector
            recurrent_norms.push_back(get_norm(dense_norm));
            // Otherwise, use the recurrent residual vector
        } else {
            auto dense_residual = gko::as<gko_dense>(residual);
            // Compute the residual vector's norm
            auto norm = compute_norm(gko::lend(dense_residual));
            // Add the computed norm to the `recurrent_norms` vector
            recurrent_norms.push_back(norm);
        }

        // Add the current iteration number to the `iterations` vector
        iterations.push_back(iteration);
        std::cout << "Iteration  " << iteration << " : residual norm " 
             << recurrent_norms[iteration] << std::endl;

    }

    // Construct the logger and store the system matrix and b vectors
    ResidualLogger(std::shared_ptr<const gko::Executor> exec,
                   const gko::LinOp *matrix, const gko_dense *b)
        : gko::log::Logger(exec, gko::log::Logger::iteration_complete_mask),
          matrix{matrix},
          b{b}
    {}

private:
    // Pointer to the system matrix
    const gko::LinOp *matrix;
    // Pointer to the right hand sides
    const gko_dense *b;
    // Vector which stores all the recurrent residual norms
    mutable std::vector<ValueType> recurrent_norms{};
    // Vector which stores all the iteration numbers
    mutable std::vector<std::size_t> iterations{};
};
// ------------- End Ginkgo custom logger definition -----------------
// -------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // 1. Parse command-line options.
    const char *mesh_file = "../data/beam-hex.mesh";
    int order = 2;
    bool static_cond = false;
    bool pa = true;
    const char *device_config = "cpu";
    bool visualization = true;
    bool pc = false;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for"
                   " isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                   "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                   "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&device_config, "-d", "--device",
                   "Device configuration string, see Device::Configure().");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&pc, "-pc", "--preconditioner",
                   "-no-pc", "--no-preconditioner", 
                   "Enable OperatorJacobi preconditioner");
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    // ---------------------------------------------------------------
    // -------------------- Start Ginkgo set-up ----------------------

    // Create Ginkgo executor.

    // This will point to the selected executor default executor
    std::shared_ptr<gko::Executor> executor;

    // We will always need an OpenMP executor.
    auto ref_executor = gko::ReferenceExecutor::create();

    // If the user has requested to use CUDA, then build a 
    // CudaExecutor and set `executor` to it; otherwise,
    // use the OmpExecutor
    bool on_device = false;
    if (!strcmp(device_config, "cuda")) {
        auto cuda_executor =
            gko::CudaExecutor::create(0, ref_executor);
        executor = cuda_executor;
        on_device = true;
        // Check we don't want full assembly + pc with this backend:
        if (!pa && pc) {
            cout << "Error: currently cannot use -no-pa, -pc, and -d cuda together" << endl;
            return 1;
        }
    } else {
        executor = ref_executor;
    }
  
    // --------------------- End Ginkgo set-up -----------------------
    // ---------------------------------------------------------------

    // 2. Enable hardware devices such as GPUs, and programming models such as
    //    CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(device_config);
    device.Print();

    // 3. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
    //    the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 4. Refine the mesh to increase the resolution. In this example we do
    //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    //    largest number that gives a final mesh with no more than 50,000
    //    elements.
    {
        int ref_levels =
            (int)floor(log(50000. / mesh->GetNE()) / log(2.) / dim);
        for (int l = 0; l < ref_levels; l++) {
            mesh->UniformRefinement();
        }
    }

    // 5. Define a finite element space on the mesh. Here we use continuous
    //    Lagrange finite elements of the specified order. If order < 1, we
    //    instead use an isoparametric/isogeometric space.
    FiniteElementCollection *fec;
    if (order > 0) {
        fec = new H1_FECollection(order, dim);
    } else if (mesh->GetNodes()) {
        fec = mesh->GetNodes()->OwnFEC();
        cout << "Using isoparametric FEs: " << fec->Name() << endl;
    } else {
        fec = new H1_FECollection(order = 1, dim);
    }
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
         << endl;

    // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking all
    //    the boundary attributes from the mesh as essential (Dirichlet) and
    //    converting them to a list of true dofs.
    Array<int> ess_tdof_list;
    if (mesh->bdr_attributes.Size()) {
        Array<int> ess_bdr(mesh->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    // 7. Set up the linear form b(.) which corresponds to the right-hand side
    // of
    //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
    //    the basis functions in the finite element fespace.
    LinearForm *b = new LinearForm(fespace);

    ConstantCoefficient one(1.0);
    b->AddDomainIntegrator(new DomainLFIntegrator(one));
    b->Assemble();

    // 8. Define the solution vector x as a finite element grid function
    //    corresponding to fespace. Initialize x with initial guess of zero,
    //    which satisfies the boundary conditions.
    GridFunction x(fespace);
    x = 0.0;

    // 9. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the Laplacian operator -Delta, by adding the
    //    Diffusion domain integrator.
    BilinearForm *a = new BilinearForm(fespace);
    if (pa) {
        a->SetAssemblyLevel(AssemblyLevel::PARTIAL);
    }
    a->AddDomainIntegrator(new DiffusionIntegrator(one));

    // 10. Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.
    if (static_cond) {
        a->EnableStaticCondensation();
    }
    a->Assemble();

    OperatorPtr A;
    Vector B, X;

    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

    // 11. Solve the linear system A X = B.

    // ----------------------------------------------------------------
    // ----------------- Start Ginkgo solve section -------------------

    // Create MFEM Vector wrappers.
    auto gko_rhs = MFEMVectorWrapper::create(executor, B.Size(), &B, on_device);
    auto gko_x = MFEMVectorWrapper::create(executor, X.Size(), &X, on_device);

    // Create MFEM Operator Wrapper for Ginkgo.
    // Note we must set A's ownership to false or else the operator will be 
    // deleted when doing a shallow copy to set mfem_oper_ in the
    // MFEMOperatorWrapper object.  
    A.SetOperatorOwner(false);
    auto oper_wrapper = MFEMOperatorWrapper::create(executor, B.Size(), A);

    // Create Ginkgo solver and solve the system.
    using cg = gko::solver::Cg<double>;

    if (!pc) {  // No preconditioner
      // Build Ginkgo CG solver factory
      auto solver_gen =
          cg::build()
              .with_criteria(
                  gko::stop::Iteration::build().with_max_iters(X.Size()).on(
                      executor),
                  gko::stop::ResidualNormReduction<>::build().with_reduction_factor(
                      1.e-6).on(executor))
              .on(executor);

      // Create and add a custom ResidualLogger to the solver factory 
      auto logger = std::make_shared<ResidualLogger<double>>(executor, 
                                                           gko::lend(oper_wrapper),
                                                           gko::lend(gko_rhs));
      solver_gen->add_logger(logger);

      // Generate CG solver for our specific wrapper MFEM Operator
      auto solver = solver_gen->generate(gko::give(oper_wrapper));
    
      // Solve system
      solver->apply(gko::lend(gko_rhs), gko::lend(gko_x));

   } else {

      // Create MFEM preconditioner
      OperatorPtr M_ptr;
      if (pa) {
        // Partial assembly version
        OperatorJacobiSmoother* M = new OperatorJacobiSmoother(*a, ess_tdof_list);
        M_ptr.Reset(M, false);
      } else {
        // Full assembly version (A is a SparseMatrix)
        DSmoother *M = new DSmoother();
        M->SetOperator(*A);
        M_ptr.Reset(M, false);
      }

      // Wrap preconditioner for Ginkgo
      auto prec_wrap = MFEMOperatorWrapper::create(executor, B.Size(), M_ptr);

      // Build Ginkgo CG solver factory with MFEM preconditioner
      auto solver_gen =
        cg::build()
            .with_criteria(
                gko::stop::Iteration::build().with_max_iters(X.Size()).on(
                    executor),
                gko::stop::ResidualNormReduction<>::build().with_reduction_factor(
                    1.e-6).on(executor))
            .with_generated_preconditioner(gko::share(prec_wrap))
            .on(executor);

      // Create and add a custom ResidualLogger to the solver factory 
      auto logger = std::make_shared<ResidualLogger<double>>(executor, 
                                                           gko::lend(oper_wrapper),
                                                           gko::lend(gko_rhs));
      solver_gen->add_logger(logger);

      // Generate CG solver for our specific wrapper MFEM Operator
      auto solver = solver_gen->generate(gko::give(oper_wrapper));
  
      // Solve system
      solver->apply(gko::lend(gko_rhs), gko::lend(gko_x));
    } 

    // ------------------ End Ginkgo solve section -------------------
    // ---------------------------------------------------------------

    // 12. Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, x);

    // 13. Save the refined mesh and the solution. This output can be viewed
    // later
    //     using GLVis: "glvis -m refined.mesh -g sol.gf".
    ofstream mesh_ofs("refined.mesh");
    mesh_ofs.precision(8);
    mesh->Print(mesh_ofs);
    ofstream sol_ofs("sol.gf");
    sol_ofs.precision(8);
    x.Save(sol_ofs);

    // 14. Send the solution by socket to a GLVis server.
    if (visualization) {
        char vishost[] = "localhost";
        int visport = 19916;
        socketstream sol_sock(vishost, visport);
        sol_sock.precision(8);
        sol_sock << "solution\n" << *mesh << x << flush;
    }

    // 15. Free the used memory.
    delete a;
    delete b;
    delete fespace;
    if (order > 0) {
        delete fec;
    }
    delete mesh;
}
