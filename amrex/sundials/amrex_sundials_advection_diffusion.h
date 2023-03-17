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

#ifndef ADVECTION_DIFFUSION_H
#define ADVECTION_DIFFUSION_H

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Sundials.H>

// ------------------------
// User defined structures
// ------------------------

// Structure for problem options
struct ProblemOpt
{
  int help               = 0;      // Print help message and exit (1)
  int plot_int           = -1;     // Enable (>0) or disable (<0) output files
  int arkode_order       = 4;      // ARKODE method order
  amrex::Real rtol       = 1.0e-4; // Relative tolerance
  amrex::Real atol       = 1.0e-9; // Absolute tolerance
  amrex::Real fixed_dt   = -1.0;   // Fixed time step size (-1 = adaptive)
  amrex::Real tfinal     = 1.0e4;  // Final integration time
  amrex::Real dtout      = 1.0e4;  // Output frequency
  int max_steps          = 10000;  // Max steps between outputs
  int nls_max_iter       = 3;      // Max number of nonlinear iterations
  int ls_max_iter        = 5;      // Max number of linear iterations
  int rhs_adv            = 1;      // Implicit (0) or Explicit (1) advection
  int use_preconditioner = 0;      // Enable (>0) preconditioner
};


// User-defined structure passed through SUNDIALS to RHS functions
struct ProblemData
{
  // Grid options
  int n_cell        = 128; // Number of cells on each side of a square domain.
  int max_grid_size = 64;  // The domain split into boxes of size max_grid_size

  // AMReX grid data structures
  amrex::Geometry* geom = nullptr;
  amrex::BoxArray* grid = nullptr;
  amrex::DistributionMapping* dmap = nullptr;

  // AMReX MLMG preconditioner data and parameters
  amrex::MultiFab* acoef = nullptr;
  amrex::MultiFab* bcoef = nullptr;

  // MLMG Preconditioner options
  int mg_agglomeration        = 1;
  int mg_consolidation        = 1;
  int mg_max_coarsening_level = 1000;
  int mg_linop_maxorder       = 2;
  int mg_max_iter             = 1000;
  int mg_max_fmg_iter         = 1000;
  int mg_fixed_iter           = 0;
  int mg_verbose              = 0;
  int mg_bottom_verbose       = 0;
  int mg_use_hypre            = 0;
  int mg_hypre_interface      = 3;
  int mg_use_petsc            = 0;
  amrex::Real mg_tol_rel      = 1.0e-6;
  amrex::Real mg_tol_abs      = 1.0e-6;

  // Advection coefficients
  amrex::Real advCoeffx = 5.0e-4;
  amrex::Real advCoeffy = 2.5e-4;

  // Diffusion coefficients
  amrex::Real diffCoeffx = 1.0e-6;
  amrex::Real diffCoeffy = 1.0e-6;

  // Scratch space for flux computation
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* flux = nullptr;
};

// -----------------------
// User defined functions
// -----------------------

// Read the problem inputs, print help info, and problem setup
int ParseInputs(ProblemOpt& prob_opt,
                ProblemData& prob_data);

// Print help information message and exit
void PrintHelp();

// Print the problem setup
void PrintSetup(ProblemOpt& prob_opt,
                ProblemData& prob_data);

// Decompose the problem in space
void SetUpGeometry(amrex::BoxArray& ba,
                   amrex::Geometry& geom,
                   ProblemData& prob_data);

// Set the problem initial state
void FillInitConds2D(amrex::MultiFab& sol,
                     const amrex::Geometry& geom);

// Setup the time integrator and advance the solution with SUNDIALS
void ComputeSolution(N_Vector nv_sol,
                     ProblemOpt* prob_opt,
                     ProblemData* prob_data);

// ---------------------------------------------
// User-supplied ODE RHS functions for SUNDIALS
// ---------------------------------------------

// Advection RHS function
int ComputeRhsAdv(amrex::Real t,
                  N_Vector nv_sol,
                  N_Vector nv_rhs,
                  void* data);

// Diffusion RHS function
int ComputeRhsDiff(amrex::Real t,
                   N_Vector nv_sol,
                   N_Vector nv_rhs,
                   void* data);

// Advection-Diffusion RHS function
int ComputeRhsAdvDiff(amrex::Real t,
                      N_Vector nv_sol,
                      N_Vector nv_rhs,
                      void* data);

// -----------------------------------------------
// Utility functions to compute ODE RHS functions
// -----------------------------------------------

// Advective portion of ODE RHS
void ComputeAdvectionUpwind(amrex::MultiFab& sol,
                            amrex::MultiFab& advection,
                            amrex::Geometry& geom,
                            amrex::Real advCoeffx,
                            amrex::Real advCoeffy);

// Diffusive portion of ODE RHS
void ComputeDiffusion(amrex::MultiFab& sol,
                      amrex::MultiFab& diff_mf,
                      amrex::MultiFab& fx_mf,
                      amrex::MultiFab& fy_mf,
                      amrex::Geometry& geom,
                      amrex::Real diffCoeffx,
                      amrex::Real diffCoeffy);

// Utility functions for computing diffusion
void ComputeDiffFlux(amrex::MultiFab& sol,
                     amrex::MultiFab& fx,
                     amrex::MultiFab& fy,
                     amrex::Geometry& geom,
                     amrex::Real diffCoeffx,
                     amrex::Real diffCoeffy);

void ComputeDivergence(amrex::MultiFab& div,
                       amrex::MultiFab& fx,
                       amrex::MultiFab& fy,
                       amrex::Geometry& geom);

// ---------------------------------------------------
// User-supplied preconditioner function for SUNDIALS
// ---------------------------------------------------

int precondition_solve(amrex::Real tn,
                       N_Vector u,
                       N_Vector fu,
                       N_Vector r,
                       N_Vector z,
                       amrex::Real gamma,
                       amrex::Real delta,
                       int lr,
                       void *user_data);

#endif
