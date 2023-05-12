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

#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "amrex_sundials_advection_diffusion.h"

using namespace amrex;

void ComputeSolution(N_Vector nv_sol, ProblemOpt* prob_opt,
                     ProblemData* prob_data)
{
  BL_PROFILE("ComputeSolution()");

  // Extract problem data and options
  Geometry* geom         = prob_data->geom;
  int       plot_int     = prob_opt->plot_int;
  int       arkode_order = prob_opt->arkode_order;
  int       nls_max_iter = prob_opt->nls_max_iter;
  int       ls_max_iter  = prob_opt->ls_max_iter;
  int       rhs_adv      = prob_opt->rhs_adv;
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
  if (rhs_adv)
  {
    // explicit advection and implicit diffusion
    arkode_mem = ARKStepCreate(ComputeRhsAdv, ComputeRhsDiff, time, nv_sol,
                               *amrex::sundials::The_Sundials_Context());
  }
  else
  {
    // implicit advection and diffusion
    arkode_mem = ARKStepCreate(nullptr, ComputeRhsAdvDiff, time, nv_sol,
                               *amrex::sundials::The_Sundials_Context());
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

  // Create and attach GMRES linear solver for Newton
  SUNLinearSolver LS = nullptr;
  if (use_preconditioner)
    LS = SUNLinSol_SPGMR(nv_sol, PREC_LEFT, ls_max_iter,
                         *amrex::sundials::The_Sundials_Context());
  else
    LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter,
                         *amrex::sundials::The_Sundials_Context());

  ier = ARKStepSetLinearSolver(arkode_mem, LS, nullptr);
  if (ier != ARKLS_SUCCESS)
  {
    amrex::Print() << "Creating linear solver failed" << std::endl;
    return;
  }

  if (use_preconditioner)
  {
    // Attach preconditioner setup/solve functions
    ier = ARKStepSetPreconditioner(arkode_mem, nullptr, precondition_solve);
    if (ier != ARKLS_SUCCESS)
    {
      amrex::Print() << "Attaching preconditioner failed" << std::endl;
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

  // Advance the solution in time
  Real tout = time + dtout; // first output time
  Real tret;                // return time
  for (int iout=0; iout < nout; iout++)
  {
    BL_PROFILE_VAR("ARKStepEvolve()", pevolve);
    ier = ARKStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(pevolve);
    if (ier < 0)
    {
      amrex::Print() << "Error in ARKStepEvolve" << std::endl;
      return;
    }

    // Get integration stats
    long nfe_evals, nfi_evals;
    ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
    if (nfe_evals > 0)
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  explicit evals = " << std::setw(7) << nfe_evals
                     << "  implicit evals = " << std::setw(7) << nfi_evals
                     << std::endl;
    else
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  RHS evals = " << std::setw(7) << nfi_evals
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
  if (amrex::ParallelDescriptor::IOProcessor())
    ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  return;
}

// -----------------------------------------------------------------------------
// Advection-Diffusion main
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  BL_PROFILE_VAR("main()", pmain);

  // Set problem data and options
  ProblemData prob_data;
  ProblemOpt  prob_opt;
  if (ParseInputs(prob_opt, prob_data))
  {
    amrex::Finalize();
    return 1;
  }
  PrintSetup(prob_opt, prob_data);

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
  ComputeSolution(nv_sol, &prob_opt, &prob_data);

  BL_PROFILE_VAR_STOP(pmain);
  amrex::Finalize();

  return 0;
}


// -----------------------------------------------------------------------------
// Parse inputs
// -----------------------------------------------------------------------------


int ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data)
{
  // ParmParse is way of reading inputs from the inputs file
  ParmParse pp;

  pp.query("help", prob_opt.help);
  if (prob_opt.help)
  {
    PrintHelp();
    return 1;
  }

  // Problem options
  pp.query("plot_int", prob_opt.plot_int);
  pp.query("arkode_order", prob_opt.arkode_order);
  pp.query("rtol", prob_opt.rtol);
  pp.query("atol", prob_opt.atol);
  pp.query("fixed_dt", prob_opt.fixed_dt);
  pp.query("tfinal", prob_opt.tfinal);
  prob_opt.dtout = prob_opt.tfinal;
  pp.query("dtout", prob_opt.dtout);
  pp.query("max_steps", prob_opt.max_steps);
  pp.query("nls_max_iter", prob_opt.nls_max_iter);
  pp.query("ls_max_iter", prob_opt.ls_max_iter);
  pp.query("rhs_adv", prob_opt.rhs_adv);
  pp.query("use_preconditioner", prob_opt.use_preconditioner);

  // Grid options
  pp.query("n_cell", prob_data.n_cell);
  pp.query("max_grid_size", prob_data.max_grid_size);

  // Advection and diffusion coefficient values
  pp.query("advCoeffx", prob_data.advCoeffx);
  pp.query("advCoeffy", prob_data.advCoeffy);
  pp.query("diffCoeffx", prob_data.diffCoeffx);
  pp.query("diffCoeffy", prob_data.diffCoeffy);

  // ParmParse options prefixed with mlmg.
  ParmParse ppmg("mlmg");

  // MLMG Preconditioner options
  ppmg.query("agglomeration", prob_data.mg_agglomeration);
  ppmg.query("consolidation", prob_data.mg_consolidation);
  ppmg.query("max_coarsening_level", prob_data.mg_max_coarsening_level);
  ppmg.query("linop_maxorder", prob_data.mg_linop_maxorder);
  ppmg.query("max_iter", prob_data.mg_max_iter);
  ppmg.query("max_fmg_iter", prob_data.mg_max_fmg_iter);
  ppmg.query("fixed_iter", prob_data.mg_fixed_iter);
  ppmg.query("verbose", prob_data.mg_verbose);
  ppmg.query("bottom_verbose", prob_data.mg_bottom_verbose);
  ppmg.query("use_hypre", prob_data.mg_use_hypre);
  ppmg.query("hypre_interface", prob_data.mg_hypre_interface);
  ppmg.query("use_petsc", prob_data.mg_use_petsc);
  ppmg.query("tol_rel", prob_data.mg_tol_rel);
  ppmg.query("tol_abs", prob_data.mg_tol_abs);

  return 0;
}


// -----------------------------------------------------------------------------
// Print help message
// -----------------------------------------------------------------------------


void PrintHelp()
{
  amrex::Print()
    << std:: endl
    << "Usage: amrex_sundials_advection_diffusion [fname] [options]" << std::endl
    << "Options:" << std::endl
    << "  help=1" << std::endl
    << "    Print this help message and exit." << std::endl
    << "  plot_int=<int>" << std::endl
    << "    enable (1) or disable (0) plots [default=0]." << std::endl
    << "  arkode_order=<int>" << std::endl
    << "    ARKStep method order [default=4]." << std::endl
    << "  nls_max_iter=<int>" << std::endl
    << "    maximum number of nonlinear iterations [default=3]." << std::endl
    << "  ls_max_iter=<int>" << std::endl
    << "    maximum number of linear iterations [default=5]." << std::endl
    << "  rhs_adv=<int>" << std::endl
    << "    treat advection implicitly (0) or explicitly (1) [default=1]." << std::endl
    << "  fixed_dt=<float>" << std::endl
    << "    use a fixed time step size (if value > 0.0) [default=-1.0]." << std::endl
    << "  rtol=<float>" << std::endl
    << "    relative tolerance for time step adaptivity [default=1e-4]." << std::endl
    << "  atol=<float>" << std::endl
    << "    absolute tolerance for time step adaptivity [default=1e-9]." << std::endl
    << "  tfinal=<float>" << std::endl
    << "    final integration time [default=1e4]." << std::endl
    << "  dtout=<float>" << std::endl
    << "    time between outputs [default=tfinal]." << std::endl
    << "  max_steps=<int>" << std::endl
    << "    maximum number of internal steps between outputs [default=10000]." << std::endl
    << "  use_preconditioner=<int>"  << std::endl
    << "    use preconditioning (1) or not (0) [default=0]." << std::endl
    << "  n_cell=<int>" << std::endl
    << "    number of cells on each side of the square domain [default=128]." << std::endl
    << "  max_grid_size=<int>" << std::endl
    << "    max size of boxes in box array [default=64]." << std::endl
    << "  advCoeffx=<float>" << std::endl
    << "    advection speed in the x-direction [default=5e-4]." << std::endl
    << "  advCoeffy=<float>" << std::endl
    << "    advection speed in the y-direction [default=2.5e-4]." << std::endl
    << "  diffCoeffx=<float>" << std::endl
    << "    diffusion coefficient in the x-direction [default=1e-6]." << std::endl
    << "  diffCoeffy=<float>" << std::endl
    << "    diffusion coefficient in the y-direction [default=1e-6]." << std::endl << std::endl
    << "If preconditioning is enabled, then additional options may be set" << std::endl
    << "(see AMReX documentation of the MLMG solver for descriptions)" << std::endl
    << "  mlmg.agglomeration=<int> [default=1]" << std::endl
    << "  mlmg.consolidation=<int> [default=1]" << std::endl
    << "  mlmg.max_coarsening_level=<int> [default=1000]" << std::endl
    << "  mlmg.linop_maxorder=<int> [default=2]" << std::endl
    << "  mlmg.max_iter=<int> [default=1000]" << std::endl
    << "  mlmg.max_fmg_iter=<int> [default=1000]" << std::endl
    << "  mlmg.fixed_iter=<int> [default=0]" << std::endl
    << "  mlmg.verbose=<int> [default=0]" << std::endl
    << "  mlmg.bottom_verbose=<int> [default=0]" << std::endl
    << "  mlmg.use_hypre=<int> [default=1]" << std::endl
    << "  mlmg.hypre_interface=<int> [default=3]" << std::endl
    << "  mlmg.use_petsc=<int> [default=0]" << std::endl
    << "  mlmg.tol_rel=<float> [default=1e-6]" << std::endl << std::endl
    << "  mlmg.tol_abs=<float> [default=1e-6]" << std::endl << std::endl
    << "If a file name 'fname' is provided, it will be parsed for each of the above" << std::endl
    << "options.  If an option is specified in both the input file and on the" << std::endl
    << "command line, then the command line option takes precedence." << std::endl << std::endl;
  return;
}


// -----------------------------------------------------------------------------
// Print problem setup
// -----------------------------------------------------------------------------


void PrintSetup(ProblemOpt& prob_opt, ProblemData& prob_data)
{
  // Ouput problem options and parameters
  amrex::Print()
    << "n_cell        = " << prob_data.n_cell        << std::endl
    << "max_grid_size = " << prob_data.max_grid_size << std::endl
    << "plot_int      = " << prob_opt.plot_int       << std::endl
    << "arkode_order  = " << prob_opt.arkode_order   << std::endl;
  if (prob_opt.rhs_adv)
    amrex::Print()
      << "ImEx treatment (implicit diffusion, explicit advection)" << std::endl;
  else
    amrex::Print()
      << "fully implicit treatment" << std::endl;
  if (prob_opt.fixed_dt > 0.0)
    amrex::Print()
      << "fixed_dt      = " << prob_opt.fixed_dt << std::endl;
  else
    amrex::Print()
      << "rtol          = " << prob_opt.rtol << std::endl
      << "atol          = " << prob_opt.atol << std::endl;
  amrex::Print()
    << "tfinal        = " << prob_opt.tfinal      << std::endl
    << "dtout         = " << prob_opt.dtout       << std::endl
    << "advCoeffx     = " << prob_data.advCoeffx  << std::endl
    << "advCoeffy     = " << prob_data.advCoeffy  << std::endl
    << "diffCoeffx    = " << prob_data.diffCoeffx << std::endl
    << "diffCoeffy    = " << prob_data.diffCoeffy << std::endl;
  amrex::Print()
    << "Newton nonlinear solver:" << std::endl
    << "  max_iter    = " << prob_opt.nls_max_iter << std::endl
    << "  ls_max_iter = " << prob_opt.ls_max_iter  << std::endl;
  if (prob_opt.use_preconditioner)
    amrex::Print()
      << "Preconditioning enabled:" << std::endl
      << "  mlmg.agglomeration        = " << prob_data.mg_agglomeration        << std::endl
      << "  mlmg.consolidation        = " << prob_data.mg_consolidation        << std::endl
      << "  mlmg.max_coarsening_level = " << prob_data.mg_max_coarsening_level << std::endl
      << "  mlmg.linop_maxorder       = " << prob_data.mg_linop_maxorder       << std::endl
      << "  mlmg.max_iter             = " << prob_data.mg_max_iter             << std::endl
      << "  mlmg.max_fmg_iter         = " << prob_data.mg_max_fmg_iter         << std::endl
      << "  mlmg.fixed_iter           = " << prob_data.mg_fixed_iter           << std::endl
      << "  mlmg.verbose              = " << prob_data.mg_verbose              << std::endl
      << "  mlmg.bottom_verbose       = " << prob_data.mg_bottom_verbose       << std::endl
      << "  mlmg.use_hypre            = " << prob_data.mg_use_hypre            << std::endl
      << "  mlmg.hypre_interface      = " << prob_data.mg_hypre_interface      << std::endl
      << "  mlmg.use_petsc            = " << prob_data.mg_use_petsc            << std::endl
      << "  mlmg.tol_rel              = " << prob_data.mg_tol_rel              << std::endl
      << "  mlmg.tol_abs              = " << prob_data.mg_tol_abs              << std::endl;
  return;
}


// -----------------------------------------------------------------------------
// Set initial state
// -----------------------------------------------------------------------------


void FillInitConds2D(MultiFab& sol, const Geometry& geom)
{
  BL_PROFILE("FillInitConds2D()");

  const auto dx = geom.CellSizeArray();
  const auto prob_lo = geom.ProbLoArray();
  const auto prob_hi = geom.ProbHiArray();

  Real sigma = 0.1;
  Real a = 1.0/(sigma*sqrt(2*M_PI));
  Real b = -0.5/(sigma*sigma);
  Real dx0 = dx[0];
  Real dx1 = dx[1];
  Real pl0 = prob_lo[0];
  Real pl1 = prob_lo[1];

  for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& fab = sol[mfi].array();

    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         Real y = pl1 + (((Real) j) + 0.5) * dx1;
         Real x = pl0 + (((Real) i) + 0.5) * dx0;
         Real r = x * x + y * y;
         fab(i,j,k,n) = a * exp(b * r);
       });
  }
}


// -----------------------------------------------------------------------------
// Setup domain
// -----------------------------------------------------------------------------


void SetUpGeometry(BoxArray& ba, Geometry& geom, ProblemData& prob_data)
{
  BL_PROFILE("SetUpGeometry()");

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


// -----------------------------------------------------------------------------
// User-supplied ODE RHS functions for SUNDIALS
// -----------------------------------------------------------------------------


int ComputeRhsAdv(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  BL_PROFILE("ComputeRhsAdv()");

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
  BL_PROFILE("ComputeRhsDiff()");

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
  BL_PROFILE("ComputeRhsAdvDiff()");

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


// -----------------------------------------------------------------------------
// Utility functions to compute ODE RHS functions
// -----------------------------------------------------------------------------


// Assumes ghost cells already filled, adds result to adv_mf MultiFab
void ComputeAdvectionUpwind(MultiFab& sol_mf, MultiFab& adv_mf, Geometry& geom,
                            Real advCoeffx, Real advCoeffy)
{
  BL_PROFILE("ComputeAdvectionUpwind()");

  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
  Real sideCoeffx = advCoeffx * dxInv;
  Real sideCoeffy = advCoeffy * dyInv;

  for (MFIter mfi(sol_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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


// Assumes ghots cells are already filled, adds result to diff_mf
void ComputeDiffusion(MultiFab& sol, MultiFab& diff_mf, MultiFab& fx_mf,
                      MultiFab& fy_mf, Geometry& geom,
                      Real diffCoeffx, Real diffCoeffy)
{
  BL_PROFILE("ComputeDiffusion()");
  ComputeDiffFlux(sol, fx_mf, fy_mf, geom, diffCoeffx, diffCoeffy);
  ComputeDivergence(diff_mf, fx_mf, fy_mf, geom);
}


// Assumes ghost cells already filled, overwrites fx_mf and fy_mf MultiFabs
void ComputeDiffFlux(MultiFab& sol_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                     Geometry& geom, Real diffCoeffx, Real diffCoeffy)
{
  BL_PROFILE("ComputeDiffFlux()");

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


// Assumes ghost cells already filled, adds result to div_mf MultiFab
void ComputeDivergence(MultiFab& div_mf, MultiFab& fx_mf,
                       MultiFab& fy_mf, Geometry& geom)
{
  BL_PROFILE("ComputeDivergence()");

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


// -----------------------------------------------------------------------------
// Preconditioner
// -----------------------------------------------------------------------------


int precondition_solve(amrex::Real tn, N_Vector u, N_Vector fu, N_Vector r,
                       N_Vector z, amrex::Real gamma, amrex::Real delta, int lr,
                       void *user_data)
{
  BL_PROFILE("precondition_solve()");

  ProblemData *prob_data = (ProblemData*) user_data;

  auto geom = *(prob_data->geom);
  auto grid = *(prob_data->grid);
  auto dmap = *(prob_data->dmap);
  auto& acoef = *(prob_data->acoef);
  auto& bcoef = *(prob_data->bcoef);

  MultiFab* solution = amrex::sundials::N_VGetVectorPointer_MultiFab(z);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(r);

  LPInfo info;
  info.setAgglomeration(prob_data->mg_agglomeration);
  info.setConsolidation(prob_data->mg_consolidation);
  info.setMaxCoarseningLevel(prob_data->mg_max_coarsening_level);

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
  mlmg.setFixedIter(prob_data->mg_fixed_iter);
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

  mlmg.solve({solution}, {rhs}, prob_data->mg_tol_rel, prob_data->mg_tol_abs);

  return 0;
}
