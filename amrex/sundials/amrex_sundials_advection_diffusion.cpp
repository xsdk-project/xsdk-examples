/* -----------------------------------------------------------------------------
 * AMReX + SUNDIALS xSDK 2D Advection-Diffusion example code
 *
 * Based on Hands-on Lessons with SUNDIALS + AMReX from the Argonne Training
 * Program in Extreme-Scale Computing (ATPESC) written by (alphabetical):
 *   David Gardner (gardner48@llnl.gov)
 *   John Loffeld (loffeld1@llnl.gov)
 *   Daniel Reynolds (reynolds@smu.edu)
 *   Donald Willcox (dewillcox@lbl.gov)
 * ---------------------------------------------------------------------------*/

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Sundials.H>

#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "amrex_sundials_advection_diffusion.h"

using namespace amrex;

void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data)
{
  // Extract problem data and options
  Geometry* geom         = prob_data->geom;
  int       plot_int     = prob_opt->plot_int;
  int       arkode_order = prob_opt->arkode_order;
  int       nls_max_iter = prob_opt->nls_max_iter;
  int       ls_max_iter  = prob_opt->ls_max_iter;
  int       rhs_adv      = prob_opt->rhs_adv;
  int       rhs_diff     = prob_opt->rhs_diff;
  Real      rtol         = prob_opt->rtol;
  Real      atol         = prob_opt->atol;
  Real      fixed_dt     = prob_opt->fixed_dt;
  Real      tfinal       = prob_opt->tfinal;
  Real      dtout        = prob_opt->dtout;
  int       max_steps    = prob_opt->max_steps;
  int use_preconditioner = prob_opt->use_preconditioner;

  // initial time, number of outputs, and error flag
  Real time = 0.0;
  int  nout = ceil(tfinal/dtout);
  int  ier  = 0;

  // Write a plotfile of the initial data
  if (plot_int > 0)
  {
    const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
    MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
    WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
  }

  // Create the ARK stepper
  void* arkode_mem = nullptr;

  if (rhs_adv > 0 && rhs_diff > 0)
  {
    if (rhs_adv > 1 && rhs_diff > 1)
    {
      // explicit advection and diffusion
      arkode_mem = ARKStepCreate(ComputeRhsAdvDiff, nullptr, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
    else if (rhs_adv > 1)
    {
      // explicit advection and implicit diffusion
      arkode_mem = ARKStepCreate(ComputeRhsAdv, ComputeRhsDiff, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
    else if (rhs_diff > 1)
    {
      // implicit advection and explicit diffusion
      arkode_mem = ARKStepCreate(ComputeRhsDiff, ComputeRhsAdv, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
    else
    {
      // implicit advection and diffusion
      arkode_mem = ARKStepCreate(nullptr, ComputeRhsAdvDiff, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
  }
  else if (rhs_adv > 0)
  {
    if (rhs_adv > 1)
    {
      // explicit advection
      arkode_mem = ARKStepCreate(ComputeRhsAdv, nullptr, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
    else
    {
      // implicit advection
      arkode_mem = ARKStepCreate(nullptr, ComputeRhsAdv, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
  }
  else if (rhs_diff > 0)
  {
    if (rhs_diff > 1)
    {
      // explicit diffusion
      arkode_mem = ARKStepCreate(ComputeRhsDiff, nullptr, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
    else
    {
      // implicit diffusion
      arkode_mem = ARKStepCreate(nullptr, ComputeRhsDiff, time, nv_sol,
                                 *amrex::sundials::The_Sundials_Context());
    }
  }
  else
  {
    amrex::Print() << "Invalid RHS options for ARKode" << std::endl;
    return;
  }

  // Attach the user data structure to ARKStep
  ARKStepSetUserData(arkode_mem, prob_data);

  // Set the method order
  ARKStepSetOrder(arkode_mem, arkode_order);

  // Set the time step size or integration tolerances
  if (fixed_dt > 0.0)
    ARKStepSetFixedStep(arkode_mem, fixed_dt);
  else
    ARKStepSStolerances(arkode_mem, atol, rtol);

  // Set the max number of steps between outputs
  ARKStepSetMaxNumSteps(arkode_mem, max_steps);

  // Attach linear solver (if needed)
  if (rhs_adv == 1 || rhs_diff == 1)
  {
    // Create and attach GMRES linear solver for Newton
    SUNLinearSolver LS;
    if (use_preconditioner)
      LS = SUNLinSol_SPGMR(nv_sol, PREC_LEFT, ls_max_iter,
                           *amrex::sundials::The_Sundials_Context());
    else
      LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter,
                           *amrex::sundials::The_Sundials_Context());

    ier = ARKStepSetLinearSolver(arkode_mem, LS, nullptr);
    if (ier != ARKLS_SUCCESS)
    {
      amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
      return;
    }

    if (use_preconditioner)
    {
      // Attach preconditioner setup/solve functions
      ier = ARKStepSetPreconditioner(arkode_mem, precondition_setup, precondition_solve);
      if (ier != ARKLS_SUCCESS)
      {
        amrex::Print() << "Attachment of preconditioner unsuccessful" << std::endl;
        return;
      }
    }

    // Set max number of nonlinear iterations
    ier = ARKStepSetMaxNonlinIters(arkode_mem, nls_max_iter);
    if (ier != ARK_SUCCESS)
    {
      amrex::Print() << "Error setting max number of nonlinear iterations" << std::endl;
      return;
    }
  }

  // Advance the solution in time
  Real tout = time + dtout; // first output time
  Real tret;                // return time
  for (int iout=0; iout < nout; iout++)
  {
    ier = ARKStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
    if (ier < 0)
    {
      amrex::Print() << "Error in ARKStepEvolve" << std::endl;
      return;
    }

    // Get integration stats
    long nfe_evals, nfi_evals;
    ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
    amrex::Print() << "t = " << std::setw(5) << tret
                   << "  explicit evals = " << std::setw(7) << nfe_evals
                   << "  implicit evals = " << std::setw(7) << nfi_evals
                   << std::endl;

    // Write output
    if (plot_int > 0)
    {
      const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
      MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
    }

    // Update output time
    tout += dtout;
    if (tout > tfinal) tout = tfinal;
  }

  // Output final solution statistics
  amrex::Print() << "\nFinal Solver Statistics:\n" << std::endl;
  ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
}


void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data)
{
  // ParmParse is way of reading inputs from the inputs file
  ParmParse pp;

  // --------------------------------------------------------------------------
  // Problem options
  // --------------------------------------------------------------------------

  // Enable (>0) or disable (<0) writing output files
  prob_opt.plot_int = -1; // plots off
  pp.query("plot_int", prob_opt.plot_int);

  // Specify the ARKode method order
  prob_opt.arkode_order = 4; // 4th order
  pp.query("arkode_order", prob_opt.arkode_order);

  // Specify the max number of nonlinear iterations
  prob_opt.nls_max_iter = 3;
  pp.query("nls_max_iter", prob_opt.nls_max_iter);

  // Specify the max number of linear iterations
  prob_opt.ls_max_iter = 5;
  pp.query("ls_max_iter", prob_opt.ls_max_iter);

  // Specify RHS functions/splitting
  prob_opt.rhs_adv  = 2; // explicit advection
  prob_opt.rhs_diff = 1; // implicit diffusion
  pp.query("rhs_adv", prob_opt.rhs_adv);
  pp.query("rhs_diff", prob_opt.rhs_diff);

  // Specify relative and absolute tolerances
  prob_opt.rtol = 1.0e-4;
  prob_opt.atol = 1.0e-9;
  pp.query("rtol", prob_opt.rtol);
  pp.query("atol", prob_opt.atol);

  // Specify a fixed time step size
  prob_opt.fixed_dt = -1.0; // diabled by default (use adaptive steps)
  pp.query("fixed_dt", prob_opt.fixed_dt);

  // Specify final time for integration
  prob_opt.tfinal = 1.0e3;
  pp.query("tfinal", prob_opt.tfinal);

  // Specify output frequency
  prob_opt.dtout = prob_opt.tfinal;
  pp.query("dtout", prob_opt.dtout);

  // Specify maximum number of steps between outputs
  prob_opt.max_steps = 10000;
  pp.query("max_steps", prob_opt.max_steps);

  // Decide whether to use a preconditioner or not
  prob_opt.use_preconditioner = 0;
  pp.query("use_preconditioner", prob_opt.use_preconditioner);

  // --------------------------------------------------------------------------
  // Problem data
  // --------------------------------------------------------------------------

  // The number of cells on each side of a square domain.
  prob_data.n_cell = 128;
  pp.query("n_cell", prob_data.n_cell);

  // The domain is broken into boxes of size max_grid_size
  prob_data.max_grid_size = 64;
  pp.query("max_grid_size", prob_data.max_grid_size);

  // Advection coefficients
  prob_data.advCoeffx = 5.0e-4;
  prob_data.advCoeffy = 2.5e-4;
  pp.query("advCoeffx", prob_data.advCoeffx);
  pp.query("advCoeffy", prob_data.advCoeffy);

  // Diffusion coefficients
  prob_data.diffCoeffx = 1.0e-6;
  prob_data.diffCoeffy = 1.0e-6;
  pp.query("diffCoeffx", prob_data.diffCoeffx);
  pp.query("diffCoeffy", prob_data.diffCoeffy);

  // MLMG options
  ParmParse ppmg("mlmg");
  prob_data.mg_agglomeration = 1;
  ppmg.query("agglomeration", prob_data.mg_agglomeration);
  prob_data.mg_consolidation = 1;
  ppmg.query("consolidation", prob_data.mg_consolidation);
  prob_data.mg_max_coarsening_level = 1000;
  ppmg.query("max_coarsening_level", prob_data.mg_max_coarsening_level);
  prob_data.mg_linop_maxorder = 2;
  ppmg.query("linop_maxorder", prob_data.mg_linop_maxorder);
  prob_data.mg_max_iter = 1000;
  ppmg.query("max_iter", prob_data.mg_max_iter);
  prob_data.mg_max_fmg_iter = 1000;
  ppmg.query("max_fmg_iter", prob_data.mg_max_fmg_iter);
  prob_data.mg_verbose = 0;
  ppmg.query("verbose", prob_data.mg_verbose);
  prob_data.mg_bottom_verbose = 0;
  ppmg.query("bottom_verbose", prob_data.mg_bottom_verbose);
  prob_data.mg_use_hypre = 1;
  ppmg.query("use_hypre", prob_data.mg_use_hypre);
  prob_data.mg_hypre_interface = 3;
  ppmg.query("hypre_interface", prob_data.mg_hypre_interface);
  prob_data.mg_use_petsc = 0;
  ppmg.query("use_petsc", prob_data.mg_use_petsc);
  prob_data.mg_tol_rel = 1.0e-6;
  ppmg.query("tol_rel", prob_data.mg_tol_rel);

  // Ouput problem options and parameters
  amrex::Print()
    << "n_cell        = " << prob_data.n_cell        << std::endl
    << "max_grid_size = " << prob_data.max_grid_size << std::endl
    << "plot_int      = " << prob_opt.plot_int       << std::endl
    << "arkode_order  = " << prob_opt.arkode_order   << std::endl
    << "rhs_adv       = " << prob_opt.rhs_adv        << std::endl
    << "rhs_diff      = " << prob_opt.rhs_diff       << std::endl;
  if (prob_opt.fixed_dt > 0.0)
    amrex::Print()
      << "fixed_dt      = " << prob_opt.fixed_dt << std::endl;
  else
    amrex::Print()
      << "rtol          = " << prob_opt.rtol << std::endl
      << "atol          = " << prob_opt.atol << std::endl;
  amrex::Print()
    << "tfinal        = " << prob_opt.tfinal << std::endl
    << "dtout         = " << prob_opt.dtout  << std::endl;
  if (prob_opt.rhs_adv > 0)
    amrex::Print()
      << "advCoeffx     = " << prob_data.advCoeffx << std::endl
      << "advCoeffy     = " << prob_data.advCoeffy << std::endl;
  if (prob_opt.rhs_diff > 0)
    amrex::Print()
      << "diffCoeffx    = " << prob_data.diffCoeffx << std::endl
      << "diffCoeffy    = " << prob_data.diffCoeffy << std::endl;
  if ((prob_opt.rhs_adv > 0) && (prob_opt.rhs_diff > 0) &&
      (prob_opt.rhs_adv != prob_opt.rhs_diff))
    if (prob_opt.rhs_diff > 1)
      amrex::Print() << "ImEx treatment: implicit advection and explicit diffusion" << std::endl;
    else
      amrex::Print() << "ImEx treatment: implicit diffusion and explicit advection" << std::endl;
  if (prob_opt.use_preconditioner)
    amrex::Print()
      << "preconditioning enabled" << std::endl
      << "  mlmg.agglomeration        = " << prob_data.mg_agglomeration << std::endl
      << "  mlmg.consolidation        = " << prob_data.mg_consolidation << std::endl
      << "  mlmg.max_coarsening_level = " << prob_data.mg_max_coarsening_level << std::endl
      << "  mlmg.linop_maxorder       = " << prob_data.mg_linop_maxorder << std::endl
      << "  mlmg.max_iter             = " << prob_data.mg_max_iter << std::endl
      << "  mlmg.max_fmg_iter         = " << prob_data.mg_max_fmg_iter << std::endl
      << "  mlmg.verbose              = " << prob_data.mg_verbose << std::endl
      << "  mlmg.bottom_verbose       = " << prob_data.mg_bottom_verbose << std::endl
      << "  mlmg.use_hypre            = " << prob_data.mg_use_hypre << std::endl
      << "  mlmg.hypre_interface      = " << prob_data.mg_hypre_interface << std::endl
      << "  mlmg.use_petsc            = " << prob_data.mg_use_petsc << std::endl
      << "  mlmg.tol_rel              = " << prob_data.mg_tol_rel << std::endl;
}


int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  DoProblem();

  amrex::Finalize();
  return 0;
}


void DoProblem()
{
  // What time is it now?  We'll use this to compute total run time.
  Real strt_time = amrex::second();

  // Set problem data and options
  ProblemData prob_data;
  ProblemOpt  prob_opt;
  ParseInputs(prob_opt, prob_data);

  // Make BoxArray and Geometry
  BoxArray ba;
  Geometry geom;
  SetUpGeometry(ba, geom, prob_data);

  // How Boxes are distrubuted among MPI processes
  DistributionMapping dm(ba);
  prob_data.dmap = &dm;

  // Allocate the solution MultiFab
  int nGhost = 1;  // number of ghost cells for each array
  int nComp  = 1;  // number of components for each array
  MultiFab sol(ba, dm, nComp, nGhost);

  // Allocate the linear solver coefficient MultiFabs
  MultiFab acoef(ba, dm, nComp, nGhost);
  MultiFab bcoef(ba, dm, nComp, nGhost);
  acoef = 1.0;
  bcoef = 1.0;
  prob_data.acoef = &acoef;
  prob_data.bcoef = &bcoef;

  // Build the flux MultiFabs
  Array<MultiFab, AMREX_SPACEDIM> flux;
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
  {
    // flux(dir) has one component, zero ghost cells, and is nodal in
    // direction dir
    BoxArray edge_ba = ba;
    edge_ba.surroundingNodes(dir);
    flux[dir].define(edge_ba, dm, 1, 0);
  }
  prob_data.flux = &flux;

  // Create an N_Vector wrapper for the solution MultiFab
  sunindextype length = nComp * prob_data.n_cell * prob_data.n_cell;
  N_Vector nv_sol     = amrex::sundials::N_VMake_MultiFab(length, &sol);

  // Set the initial condition
  FillInitConds2D(sol, geom);

  // Integrate in time
  ComputeSolutionARK(nv_sol, &prob_opt, &prob_data);

  // Call the timer again and compute the maximum difference between the start
  // time and stop time over all processors
  Real stop_time = amrex::second() - strt_time;
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

  // Tell the I/O Processor to write out the "run time"
  amrex::Print() << "Run time = " << stop_time << std::endl;
}


void FillInitConds2D(MultiFab& sol, const Geometry& geom)
{
  const auto dx = geom.CellSize();
  const auto prob_lo = geom.ProbLo();
  const auto prob_hi = geom.ProbHi();

  Real sigma = 0.1;
  Real a = 1.0/(sigma*sqrt(2*M_PI));
  Real b = -0.5/(sigma*sigma);

  for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& fab = sol[mfi].array();

    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         Real y = prob_lo[1] + (((Real) j) + 0.5) * dx[1];
         Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];
         Real r = x * x + y * y;
         fab(i,j,k,n) = a * exp(b * r);
       });
  }
}

void SetUpGeometry(BoxArray& ba, Geometry& geom, ProblemData& prob_data)
{
  // Extract problem options
  int n_cell = prob_data.n_cell;
  int max_grid_size = prob_data.max_grid_size;

  IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
  IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
  Box domain(dom_lo, dom_hi); // cell-centered

  // Initialize the boxarray "ba" from the single box "domain"
  ba.define(domain);

  // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a
  // direction
  ba.maxSize(max_grid_size);

  // This defines the physical box, [-1,1] in each direction.
  RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                   {AMREX_D_DECL(1.0, 1.0, 1.0)});

  // This defines a Geometry object
  Vector<int> is_periodic(AMREX_SPACEDIM, 1);  // periodic in all direction
  geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

  prob_data.geom = &geom;
  prob_data.grid = &ba;
}


/* ---------------------------------------------------------------------------
 * SUNDIALS RHS functions
 * ---------------------------------------------------------------------------*/

int ComputeRhsAdv(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  // extract MultiFabs
  MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_rhs);

  // extract problem data
  ProblemData* prob_data = (ProblemData*) data;
  Geometry* geom = prob_data->geom;
  Real advCoeffx = prob_data->advCoeffx;
  Real advCoeffy = prob_data->advCoeffy;

  // clear the RHS
  *rhs = 0.0;

  // fill ghost cells
  sol->FillBoundary(geom->periodicity());

  // compute advection
  ComputeAdvectionUpwind(*sol, *rhs, *geom, advCoeffx, advCoeffy);

  return 0;
}

int ComputeRhsDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  // extract MultiFabs
  MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_rhs);

  // extract problem data
  ProblemData *prob_data = (ProblemData*) data;
  Geometry* geom = prob_data->geom;
  Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
  Real diffCoeffx = prob_data->diffCoeffx;
  Real diffCoeffy = prob_data->diffCoeffy;

  // fill ghost cells
  sol->FillBoundary(geom->periodicity());

  // clear the RHS
  *rhs = 0.0;

  // compute diffusion
  ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom,
                   diffCoeffx, diffCoeffy);

  return 0;
}

int ComputeRhsAdvDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  // extract MultiFabs
  MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_rhs);

  // extract problem data
  ProblemData* prob_data = (ProblemData*) data;
  Geometry* geom = prob_data->geom;
  Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
  Real advCoeffx = prob_data->advCoeffx;
  Real advCoeffy = prob_data->advCoeffy;
  Real diffCoeffx = prob_data->diffCoeffx;
  Real diffCoeffy = prob_data->diffCoeffy;

  // clear the RHS
  *rhs = 0.0;

  // fill ghost cells
  sol->FillBoundary(geom->periodicity());

  // compute advection
  ComputeAdvectionUpwind(*sol, *rhs, *geom, advCoeffx, advCoeffy);

  // compute diffusion
  ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom,
                   diffCoeffx, diffCoeffy);

  return 0;
}

/* ---------------------------------------------------------------------------
 * Advection RHS functions
 * ---------------------------------------------------------------------------*/

// Assumes ghost cells already filled
// Adds result to adv_mf MultiFab
void ComputeAdvectionUpwind(MultiFab& sol_mf, MultiFab& adv_mf, Geometry& geom,
                            Real advCoeffx, Real advCoeffy)
{
  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
  Real sideCoeffx = advCoeffx * dxInv;
  Real sideCoeffy = advCoeffy * dyInv;

  for (MFIter mfi(sol_mf,TilingIfNotGPU); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& sol_fab = sol_mf[mfi].array();
    Array4<Real> const& adv_fab = adv_mf[mfi].array();

    // x-direction
    if (advCoeffx > 0)
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffx *
             (sol_fab(i,j,k,n) - sol_fab(i-1,j,k,n));
         });
    }
    else
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffx *
             (sol_fab(i+1,j,k,n) - sol_fab(i,j,k,n));
         });
    }

    // y-direction
    if (advCoeffy > 0)
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffy *
             (sol_fab(i,j,k,n) - sol_fab(i,j-1,k,n));
         });
    }
    else
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffy *
             (sol_fab(i,j+1,k,n) - sol_fab(i,j,k,n));
         });
    }
  }
}

/* ---------------------------------------------------------------------------
 * Diffusion RHS functions
 * ---------------------------------------------------------------------------*/

// Assumes ghots cells are already filled
// Adds result to diff_mf
void ComputeDiffusion(MultiFab& sol, MultiFab& diff_mf, MultiFab& fx_mf,
                      MultiFab& fy_mf, Geometry& geom,
                      Real diffCoeffx, Real diffCoeffy)
{
  ComputeDiffFlux(sol, fx_mf, fy_mf, geom, diffCoeffx, diffCoeffy);
  ComputeDivergence(diff_mf, fx_mf, fy_mf, geom);
}

// Assumes ghost cells already filled
// Overwrites fx_mf and fy_mf MultiFabs
void ComputeDiffFlux(MultiFab& sol_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                     Geometry& geom, Real diffCoeffx, Real diffCoeffy)
{
  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
  Real coeffX = diffCoeffx * dxInv;
  Real coeffY = diffCoeffy * dyInv;

  for (MFIter mfi(sol_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& sol = sol_mf[mfi].array();
    Array4<Real> const& fx = fx_mf[mfi].array();
    Array4<Real> const& fy = fy_mf[mfi].array();

    // x-flux
    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         // always use zero component for flux
         fx(i,j,k,0) = coeffX * (sol(i,j,k,n) - sol(i-1,j,k,n));
       });

    // y-flux
    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         // always use zero component for flux
         fy(i,j,k,0) = coeffY * (sol(i,j,k,n) - sol(i,j-1,k,n));
       });
  }
}

// Assumes ghost cells already filled
// Adds result to div_mf MultiFab
void ComputeDivergence(MultiFab& div_mf, MultiFab& fx_mf,
                       MultiFab& fy_mf, Geometry& geom)
{
  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh

  for (MFIter mfi(div_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& div = div_mf[mfi].array();
    Array4<Real> const& fx = fx_mf[mfi].array();
    Array4<Real> const& fy = fy_mf[mfi].array();

    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         // always use zero component for flux
         div(i,j,k,n) += (dxInv * (fx(i+1,j,k,0) - fx(i,j,k,0)) +
                          dyInv * (fy(i,j+1,k,0) - fy(i,j,k,0)));
       });
  }
}

/* ---------------------------------------------------------------------------
 * Preconditioning routines
 * ---------------------------------------------------------------------------*/

int precondition_setup(realtype tn, N_Vector u, N_Vector fu,
                       booleantype jok, booleantype *jcurPtr,
                       realtype gamma, void *user_data)
{
  return 0;
}

int precondition_solve(realtype tn, N_Vector u, N_Vector fu,
                       N_Vector r, N_Vector z,
                       realtype gamma, realtype delta,
                       int lr, void *user_data)
{
  ProblemData *prob_data = (ProblemData*) user_data;

  auto geom = *(prob_data->geom);
  auto grid = *(prob_data->grid);
  auto dmap = *(prob_data->dmap);
  auto& acoef = *(prob_data->acoef);
  auto& bcoef = *(prob_data->acoef);

  MultiFab* solution = amrex::sundials::N_VGetVectorPointer_MultiFab(z);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(r);

  LPInfo info;
  info.setAgglomeration(prob_data->mg_agglomeration);
  info.setConsolidation(prob_data->mg_consolidation);
  info.setMaxCoarseningLevel(prob_data->mg_max_coarsening_level);

  const Real tol_abs = 0.0;
  const Real ascalar = 1.0;
  const Real bscalar = gamma;

  MLABecLaplacian mlabec({geom}, {grid}, {dmap}, info);

  mlabec.setMaxOrder(prob_data->mg_linop_maxorder);

  // Set periodic BC
  mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                   LinOpBCType::Periodic,
                                   LinOpBCType::Periodic)},
    {AMREX_D_DECL(LinOpBCType::Periodic,
                  LinOpBCType::Periodic,
                  LinOpBCType::Periodic)});

  mlabec.setLevelBC(0, nullptr);

  mlabec.setScalars(ascalar, bscalar);

  mlabec.setACoeffs(0, acoef);

  Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
  {
    const BoxArray& ba = amrex::convert(bcoef.boxArray(),
                                        IntVect::TheDimensionVector(idim));
    face_bcoef[idim].define(ba, bcoef.DistributionMap(), 1, 0);

    switch (idim)
    {
    case 0:
      face_bcoef[idim] = prob_data->diffCoeffx;
    case 1:
      face_bcoef[idim] = prob_data->diffCoeffy;
    }
  }

  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

  MLMG mlmg(mlabec);
  mlmg.setMaxIter(prob_data->mg_max_iter);
  mlmg.setMaxFmgIter(prob_data->mg_max_fmg_iter);
  mlmg.setVerbose(prob_data->mg_verbose);
  mlmg.setBottomVerbose(prob_data->mg_bottom_verbose);
#ifdef AMREX_USE_HYPRE
  if (prob_data->mg_use_hypre)
  {
    mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    if (prob_data->mg_hypre_interface == 1)
      mlmg.setHypreInterface(amrex::Hypre::Interface::structed);
    else if (prob_data->mg_hypre_interface == 2)
      mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
    else
      mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
  }
#endif
#ifdef AMREX_USE_PETSC
  if (prob_data->mg_use_petsc)
  {
    mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
  }
#endif

  mlmg.solve({solution}, {rhs}, prob_data->mg_tol_rel, tol_abs);

  return 0;
}
