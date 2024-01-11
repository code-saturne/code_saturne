/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscdraw.h>
#include <petscviewer.h>
#include <petscksp.h>
#endif

#if defined(HAVE_HYPRE)
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-linear_solvers.c
 *
 * \brief Linear solvers examples.
 *
 * See \ref parameters for examples.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

#if defined(HAVE_PETSC)

/*----------------------------------------------------------------------------
 * User function example for setup options of a PETSc KSP solver.
 *
 * This function is called at the end of the setup stage for a KSP solver.
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp_p   <-> pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

/* Conjugate gradient with Jacobi preconditioning */
/*------------------------------------------------*/

/*! [sles_petsc_hook_1] */
static void
_petsc_p_setup_hook(void  *context,
                    void  *ksp_p)
{
  CS_UNUSED(context);
  KSP ksp = (KSP)ksp_p;
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCJACOBI);  /* Jacobi (diagonal) preconditioning */
}
/*! [sles_petsc_hook_1] */

/* Conjugate gradient with GAMG preconditioning */
/*----------------------------------------------*/

/*! [sles_petsc_hook_gamg] */
static void
_petsc_p_setup_hook_gamg(void  *context,
                         void  *ksp_p)
{
  CS_UNUSED(context);
  KSP ksp = (KSP)ksp_p;
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCGAMG);  /* GAMG (geometric-algebraic multigrid)
                             preconditioning */
}
/*! [sles_petsc_hook_gamg] */

/* Conjugate gradient with HYPRE BoomerAMG preconditioning */
/*---------------------------------------------------------*/

/*! [sles_petsc_hook_bamg] */
static void
_petsc_p_setup_hook_bamg(void  *context,
                         void  *ksp_p)
{
  CS_UNUSED(context);
  KSP ksp = (KSP)ksp_p;
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCHYPRE);  /* HYPRE BoomerAMG preconditioning */
}
/*! [sles_petsc_hook_bamg] */

/*----------------------------------------------------------------------------
 * User function example for setup options of a PETSc KSP solver.
 *
 * This example outputs the matrix structure and values, based on several
 * options.
 *
 * This function is called the end of the setup stage for a KSP solver.
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp_p   <-> pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

/*! [sles_petsc_hook_view] */
static void
_petsc_p_setup_hook_view(void  *context,
                         void  *ksp_p)
{
  CS_UNUSED(context);

  KSP ksp = (KSP)ksp_p;

  const char *p = getenv("CS_USER_PETSC_MAT_VIEW");

  if (p != NULL) {

    /* Get system and preconditioner matrices */

    Mat a, pa;
    KSPGetOperators(ksp, &a, &pa);

    /* Output matrix in several ways depending on
       CS_USER_PETSC_MAT_VIEW environment variable */

    if (strcmp(p, "DEFAULT") == 0)
      MatView(a, PETSC_VIEWER_DEFAULT);

    else if (strcmp(p, "DRAW_WORLD") == 0)
      MatView(a, PETSC_VIEWER_DRAW_WORLD);

    else if (strcmp(p, "DRAW") == 0) {

      PetscViewer viewer;
      PetscDraw draw;
      PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, "PETSc View",
                          0, 0, 600, 600, &viewer);
      PetscViewerDrawGetDraw(viewer, 0, &draw);
      PetscViewerDrawSetPause(viewer, -1);
      MatView(a, viewer);
      PetscDrawPause(draw);

      PetscViewerDestroy(&viewer);

    }

  }
}
/*! [sles_petsc_hook_view] */

/*----------------------------------------------------------------------------
 * Function pointer for user settings of a PETSc KSP solver setup.
 *
 * This function is called the end of the setup stage for a KSP solver.
 *
 * Note that using the advanced KSPSetPostSolve and KSPSetPreSolve functions,
 * this also allows setting furthur function pointers for pre and post-solve
 * operations (see the PETSc documentation).
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp_p   <-> pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

void
cs_user_sles_petsc_hook(void  *context,
                        void  *ksp_p)
{
  CS_UNUSED(ksp_p);

  /*! [sles_petsc_cdo_hook] */
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;

  if (slesp == NULL)
    return;

  /* Usually the name of the equation or the field id of the associated
     variable */
  if (strcmp(slesp->name, "Name_Of_The_System") == 0) {

    /* Assume a PETSc version greater or equal to 3.7.0 */
    if (slesp->precond == CS_PARAM_PRECOND_AMG) {
      if (slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_V) {

        PetscOptionsSetValue(NULL,
                             "-pc_hypre_boomeramg_strong_threshold", "0.7");

      }
    }

  }
  /*! [sles_petsc_cdo_hook] */
}

#endif /* defined(HAVE_PETSC) */

#if defined(HAVE_HYPRE)

/*----------------------------------------------------------------------------
 * User function example for setup options of a Hypre KSP solver.
 *
 * This function is called at the end of the setup stage for a KSP solver.
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * Check HYPRE documentation for available options:
 * https://hypre.readthedocs.io/en/latest/index.html
 *
 * parameters:
 *   verbosity <-- verbosity level
 *   context   <-> pointer to optional (untyped) value or structure
 *   solver    <->  handle to HYPRE solver
 *----------------------------------------------------------------------------*/

/* Conjugate gradient with BoomerAMG preconditioning */
/*---------------------------------------------------*/

/*! [sles_hypre_hook_1] */
static void
_hypre_p_setup_hook(int    verbosity,
                    void  *context,
                    void  *solver_p)
{
  CS_NO_WARN_IF_UNUSED(verbosity);
  CS_NO_WARN_IF_UNUSED(context);

  HYPRE_Solver  solver = solver_p;
  HYPRE_Solver  precond = NULL;

  /* Get pointer to preconditioner, based on solver type (here for PCG) */
  HYPRE_PCGGetPrecond(solver, &precond);

  /* Assuming the preconditioner is BoomerAMG, set options */
  HYPRE_BoomerAMGSetCoarsenType(precond, 8) ;        /* HMIS */
  HYPRE_BoomerAMGSetAggNumLevels(precond, 2);
  HYPRE_BoomerAMGSetPMaxElmts(precond, 4);
  HYPRE_BoomerAMGSetInterpType(precond, 7);          /* extended+i */
  HYPRE_BoomerAMGSetStrongThreshold(precond, 0.5);   /* 2d=>0.25 3d=>0.5 */
  HYPRE_BoomerAMGSetRelaxType(precond, 6);   /* Sym G.S./Jacobi hybrid */
  HYPRE_BoomerAMGSetRelaxOrder(precond, 0);
}
/*! [sles_hypre_hook_1] */

#endif /* defined(HAVE_HYPRE) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *
 * Available native iterative linear solvers include conjugate gradient,
 * Jacobi, BiCGStab, BiCGStab2, and GMRES. For symmetric linear systems,
 * an algebraic multigrid solver is available (and recommended).
 *
 * External solvers may also be setup using this function, the cs_sles_t
 * mechanism allowing such through user-define functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void)
{
  /* Available native iterative linear solvers are:
   *
   *  CS_SLES_PCG                 (preconditioned conjugate gradient)
   *  CS_SLES_JACOBI              (Jacobi)
   *  CS_SLES_BICGSTAB            (Bi-conjugate gradient stabilized)
   *  CS_SLES_BICGSTAB2           (BiCGStab2)
   *  CS_SLES_GMRES               (generalized minimal residual)
   *  CS_SLES_P_GAUSS_SEIDEL      (process-local Gauss-Seidel)
   *  CS_SLES_P_SYM_GAUSS_SEIDEL  (process-local symmetric Gauss-Seidel)
   *  CS_SLES_PCR3                (3-layer conjugate residual)
   *
   *  The multigrid solver uses the conjugate gradient as a smoother
   *  and coarse solver by default, but this behavior may be modified. */

  /* Example: use multigrid for wall distance computation */
  /*------------------------------------------------------*/

  /*! [sles_wall_dist] */
  cs_multigrid_define(-1, "wall_distance", CS_MULTIGRID_V_CYCLE);
  /*! [sles_wall_dist] */

  /* Example: use BiCGStab2 for user variable (named user_1) */
  /*---------------------------------------------------------*/

  /*! [sles_user_1] */
  cs_field_t *cvar_user_1 = cs_field_by_name_try("user_1");
  if (cvar_user_1 != NULL) {
    cs_sles_it_define(cvar_user_1->id,
                      NULL,   /* name passed is NULL if field_id > -1 */
                      CS_SLES_BICGSTAB2,
                      1,      /* polynomial precond. degree (default 0) */
                      10000); /* n_max_iter */
  }
  /*! [sles_user_1] */

  /* Example: increase verbosity parameters for pressure */
  /*-----------------------------------------------------*/

  /*! [sles_verbosity_1] */
  {
    cs_sles_t *sles_p = cs_sles_find_or_add(CS_F_(p)->id, NULL);
    cs_sles_set_verbosity(sles_p, 4);
  }
  /*! [sles_verbosity_1] */

  /* Example: visualize local error for velocity and pressure */
  /*----------------------------------------------------------*/

  /*! [sles_viz_1] */
  {
    cs_sles_t *sles_p = cs_sles_find_or_add(CS_F_(p)->id, NULL);
    cs_sles_set_post_output(sles_p, CS_POST_WRITER_DEFAULT);

    cs_sles_t *sles_u = cs_sles_find_or_add(CS_F_(vel)->id, NULL);
    cs_sles_set_post_output(sles_u, CS_POST_WRITER_DEFAULT);
  }
  /*! [sles_viz_1] */

  /* Example: change multigrid parameters for pressure */
  /*---------------------------------------------------*/

  /*! [sles_mgp_1] */
  {
    cs_multigrid_t *mg = cs_multigrid_define(CS_F_(p)->id,
                                             NULL,
                                             CS_MULTIGRID_V_CYCLE);

    cs_multigrid_set_coarsening_options(mg,
                                        3,    /* aggregation_limit (default 3) */
                                        0,    /* coarsening_type (default 0) */
                                        10,   /* n_max_levels (default 25) */
                                        30,   /* min_g_cells (default 30) */
                                        0.95, /* P0P1 relaxation (default 0.95) */
                                        20);  /* postprocessing (default 0) */

    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_JACOBI, /* descent smoother type (default: CS_SLES_PCG) */
       CS_SLES_JACOBI, /* ascent smoother type (default: CS_SLES_PCG) */
       CS_SLES_PCG,    /* coarse solver type (default: CS_SLES_PCG) */
       50,             /* n max cycles (default 100) */
       5,              /* n max iter for descent (default 2) */
       5,              /* n max iter for asscent (default 10) */
       1000,           /* n max iter coarse solver (default 10000) */
       0,              /* polynomial precond. degree descent (default 0) */
       0,              /* polynomial precond. degree ascent (default 0) */
       1,              /* polynomial precond. degree coarse (default 0) */
       -1.0,           /* precision multiplier descent (< 0 forces max iters) */
       -1.0,           /* precision multiplier ascent (< 0 forces max iters) */
       0.1);           /* requested precision multiplier coarse (default 1) */

  }
  /*! [sles_mgp_1] */

  /* Set parallel grid merging options for all multigrid solvers */
  /*-------------------------------------------------------------*/

  /*! [sles_mg_parall] */
  {
    cs_multigrid_t *mg = cs_multigrid_define(CS_F_(p)->id,
                                             NULL,
                                             CS_MULTIGRID_V_CYCLE);

    cs_multigrid_set_merge_options(mg,
                                   4,    /* # of ranks merged at a time */
                                   300,  /* mean # of cells under which we merge */
                                   500); /* global # of cells under which we merge */
  }
  /*! [sles_mg_parall] */

  /* Example: conjugate gradient preconditioned by multigrid for pressure */
  /*----------------------------------------------------------------------*/

  /*! [sles_mgp_2] */
  {
    cs_sles_it_t *c = cs_sles_it_define(CS_F_(p)->id,
                                        NULL,
                                        CS_SLES_FCG,
                                        -1,
                                        10000);
    cs_sles_pc_t *pc = cs_multigrid_pc_create(CS_MULTIGRID_V_CYCLE);
    cs_multigrid_t *mg = cs_sles_pc_get_context(pc);
    cs_sles_it_transfer_pc(c, &pc);

    assert(strcmp(cs_sles_pc_get_type(cs_sles_it_get_pc(c)), "multigrid") == 0);

    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_P_GAUSS_SEIDEL, /* descent smoother (CS_SLES_P_SYM_GAUSS_SEIDEL) */
       CS_SLES_P_GAUSS_SEIDEL, /* ascent smoother (CS_SLES_P_SYM_GAUSS_SEIDEL) */
       CS_SLES_PCG,            /* coarse solver (CS_SLES_P_GAUSS_SEIDEL) */
       1,              /* n max cycles (default 1) */
       1,              /* n max iter for descent (default 1) */
       1,              /* n max iter for asscent (default 1) */
       500,            /* n max iter coarse solver (default 1) */
       0,              /* polynomial precond. degree descent (default) */
       0,              /* polynomial precond. degree ascent (default) */
       0,              /* polynomial precond. degree coarse (default 0) */
       -1.0,           /* precision multiplier descent (< 0 forces max iters) */
       -1.0,           /* precision multiplier ascent (< 0 forces max iters) */
       1.0);           /* requested precision multiplier coarse (default 1) */

  }
  /*! [sles_mgp_2] */

  /* Example: conjugate gradient preconditioned by K-cycle multigrid in the *
   *          the saddle-point system for coupled velocity-pressure relying *
   *          on CDO face-based schemes. One considers this solver for the  *
   *          velocity block (i.e. the momentum equation). Case of a Stokes *
   *          equations                                                     */
  /*------------------------------------------------------------------------*/

  /*! [sles_kamg_momentum] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
    cs_param_sles_t  *slesp = eqp->sles_param;
    assert(slesp->field_id > -1);

    /* In case of an in-house K-cylcle multigrid as a preconditioner of a
       linear iterative solver */

    if (slesp->precond == CS_PARAM_PRECOND_AMG) { /* If multigrid is chosen as
                                                     preconditioner */

      if (slesp->amg_type == CS_PARAM_AMG_HOUSE_K) { /* If this is a K-cycle
                                                        algebraic multigrid */

        /* Retrieve the different context structures to apply additional
           settings */

        cs_sles_t  *sles = cs_sles_find_or_add(slesp->field_id, NULL);
        cs_sles_it_t  *itsol = cs_sles_get_context(sles);
        cs_sles_pc_t  *pc = cs_sles_it_get_pc(itsol);
        cs_multigrid_t  *mg = NULL;

        if (itsol == NULL) { /* Not already defined */

          if (slesp->solver == CS_PARAM_ITSOL_CG ||
              slesp->solver == CS_PARAM_ITSOL_FCG)
            itsol =  cs_sles_it_define(slesp->field_id, NULL,
                                       CS_SLES_IPCG, -1,
                                       slesp->cvg_param.n_max_iter);
          else
            bft_error(__FILE__, __LINE__, 0,
                      " %s: Case not treated.\n", __func__);

        }
        assert(itsol != NULL);

        if (pc == NULL) { /* Not already defined */

          pc = cs_multigrid_pc_create(CS_MULTIGRID_K_CYCLE);
          mg = cs_sles_pc_get_context(pc);
          cs_sles_it_transfer_pc(itsol, &pc);

        }
        else
          mg = cs_sles_pc_get_context(pc);

        assert(mg != NULL && pc != NULL);

        /* Available settings:
         * - max. number of elements in an aggregation
         * - type of algorithm to perform the aggregation
         * - max. number of levels (i.e. grids)
         * - max global number of rows at the coarsest level
         * - type of relaxation (weighting between a P_0 and P_1). For K-cycle,
         * this should be equal to 0.
         * - Activation of the postprocessing for the aggregation if > 0.
         * Aggregation set is numbered by its coarse row number modulo this
         * value
         */

        cs_multigrid_set_solver_options
          (mg,
           CS_SLES_TS_F_GAUSS_SEIDEL,
           CS_SLES_TS_B_GAUSS_SEIDEL,
           CS_SLES_PCG,       /* coarse solver */
           1,                 /* n_max_cycles */
           1,                 /* n_max_iter_descent, */
           1,                 /* n_max_iter_ascent */
           200,               /* n_max_iter_coarse */
           0,                 /* poly_degree_descent */
           0,                 /* poly_degree_ascent */
           0,                 /* poly_degree_coarse */
           -1.0,              /* precision_mult_descent */
           -1.0,              /* precision_mult_ascent */
           1.0);              /* precision_mult_coarse */

        cs_multigrid_set_coarsening_options(mg,
                                            8,    /* aggregation_limit*/
                                            CS_GRID_COARSENING_SPD_PW,
                                            10,   /* n_max_levels */
                                            500,  /* min_g_cells */
                                            0.,   /* P0P1 relaxation */
                                            0);   /* postprocess */

      } /* K-cycle */
    }   /* Multigrid as preconditioner */
  }
  /*! [sles_kamg_momentum] */

  /* Set a non-default linear solver for DOM radiation. */
  /*----------------------------------------------------*/

  /* The solver must be set for each direction; here, we assume
     a quadrature with 32 directions is used */

  /*! [sles_rad_dom_1] */
  {
    for (int i = 0; i < 32; i++) {
      char name[16];
      sprintf(name, "radiation_%03d", i+1);
      cs_sles_it_define(-1,
                        name,
                        CS_SLES_JACOBI,
                        0,      /* poly_degree */
                        1000);  /* n_max_iter */

    }
  }
  /*! [sles_rad_dom_1] */


  /* Example: activate convergence plot for pressure */
  /*-------------------------------------------------*/

  /*! [sles_plot_1] */
  {
    const cs_field_t *f = CS_F_(p);
    cs_sles_t *sles_p = cs_sles_find_or_add(f->id, NULL);

    bool use_iteration = true; /* use iteration or wall clock time for axis */

    if (strcmp(cs_sles_get_type(sles_p), "cs_sles_it_t") == 0) {
      cs_sles_it_t *c = cs_sles_get_context(sles_p);
      cs_sles_it_set_plot_options(c, f->name, use_iteration);
    }
    else if (strcmp(cs_sles_get_type(sles_p), "cs_multigrid_t") == 0) {
      cs_multigrid_t *c = cs_sles_get_context(sles_p);
      cs_multigrid_set_plot_options(c, f->name, use_iteration);
    }

  }
  /*! [sles_plot_1] */

#if defined(HAVE_PETSC)

  /* Setting global options for PETSc */
  /*----------------------------------*/

  /*! [sles_petsc_1] */
  {
    /* Initialization must be called before setting options;
       it does not need to be called before calling
       cs_sles_petsc_define(), as this is handled automatically. */

    PETSC_COMM_WORLD = cs_glob_mpi_comm;
    PetscInitializeNoArguments();

    /* See the PETSc documentation for the options database */
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsSetValue(NULL, "-ksp_type", "cg");
    PetscOptionsSetValue(NULL, "-pc_type", "jacobi");
#else
    PetscOptionsSetValue("-ksp_type", "cg");
    PetscOptionsSetValue("-pc_type", "jacobi");
#endif
  }
  /*! [sles_petsc_1] */

  /* Setting pressure solver with PETSc */
  /*------------------------------------*/

  /*! [sles_petsc_2] */
  {
    cs_sles_petsc_define(CS_F_(p)->id,
                         NULL,
                         MATSHELL,
                         _petsc_p_setup_hook,
                         NULL);

  }
  /*! [sles_petsc_2] */

  /* Setting global options for PETSc with GAMG preconditioner */
  /*-----------------------------------------------------------*/

  /*! [sles_petsc_gamg_1] */
  {
    /* Initialization must be called before setting options;
       it does not need to be called before calling
       cs_sles_petsc_define(), as this is handled automatically. */

    PETSC_COMM_WORLD = cs_glob_mpi_comm;
    PetscInitializeNoArguments();

    /* See the PETSc documentation for the options database */
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsSetValue(NULL, "-ksp_type", "cg");
    PetscOptionsSetValue(NULL, "-pc_type", "gamg");
    PetscOptionsSetValue(NULL, "-pc_gamg_agg_nsmooths", "1");
    PetscOptionsSetValue(NULL, "-mg_levels_ksp_type", "richardson");
    PetscOptionsSetValue(NULL, "-mg_levels_pc_type", "sor");
    PetscOptionsSetValue(NULL, "-mg_levels_ksp_max_it", "1");
    PetscOptionsSetValue(NULL, "-pc_gamg_threshold", "0.02");
    PetscOptionsSetValue(NULL, "-pc_gamg_reuse_interpolation", "TRUE");
    PetscOptionsSetValue(NULL, "-pc_gamg_square_graph", "4");
#else
    PetscOptionsSetValue("-ksp_type", "cg");
    PetscOptionsSetValue("-pc_type", "gamg");
    PetscOptionsSetValue("-pc_gamg_agg_nsmooths", "1");
    PetscOptionsSetValue("-mg_levels_ksp_type", "richardson");
    PetscOptionsSetValue("-mg_levels_pc_type", "sor");
    PetscOptionsSetValue("-mg_levels_ksp_max_it", "1");
    PetscOptionsSetValue("-pc_gamg_threshold", "0.02");
    PetscOptionsSetValue("-pc_gamg_reuse_interpolation", "TRUE");
    PetscOptionsSetValue("-pc_gamg_square_graph", "4");
#endif
  }
  /*! [sles_petsc_gamg_1] */

  /* Setting pressure solver with PETSc and GAMG preconditioner */
  /*------------------------------------------------------------*/

  /*! [sles_petsc_gamg_2] */
  {
    cs_sles_petsc_define(CS_F_(p)->id,
                         NULL,
                         MATMPIAIJ,
                         _petsc_p_setup_hook_gamg,
                         NULL);

  }
  /*! [sles_petsc_gamg_2] */

  /* Setting global options for PETSc with HYPRE BoomerAMG preconditioner */
  /*----------------------------------------------------------------------*/

  /*! [sles_petsc_bamg_1] */
  {

    /* Initialization must be called before setting options;
       it does not need to be called before calling
       cs_sles_petsc_define(), as this is handled automatically. */

    PETSC_COMM_WORLD = cs_glob_mpi_comm;
    PetscInitializeNoArguments();

    /* See the PETSc documentation for the options database */
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsSetValue(NULL, "-ksp_type", "cg");
    PetscOptionsSetValue(NULL, "-pc_type", "hypre");
    PetscOptionsSetValue(NULL, "-pc_hypre_type","boomeramg");
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", "HMIS");
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_interp_type", "ext+i-cc");
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_agg_nl","2");
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_P_max","4");
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_strong_threshold", "0.5");
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_no_CF","");
#else
    PetscOptionsSetValue("-ksp_type", "cg");
    PetscOptionsSetValue("-pc_type", "hypre");
    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");
    PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type", "HMIS");
    PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type", "ext+i-cc");
    PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl", "2");
    PetscOptionsSetValue("-pc_hypre_boomeramg_P_max", "4");
    PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.5");
    PetscOptionsSetValue("-pc_hypre_boomeramg_no_CF", "");
#endif
  }
  /*! [sles_petsc_bamg_1] */

  /* Setting pressure solver with PETSc and BoomerAMG preconditioner */
  /*-----------------------------------------------------------------*/

  /*! [sles_petsc_bamg_2] */
  {
    cs_sles_petsc_define(CS_F_(p)->id,
                         NULL,
                         MATMPIAIJ,
                         _petsc_p_setup_hook_bamg,
                         NULL);

  }
  /*! [sles_petsc_bamg_2] */

#endif /* defined(HAVE_PETSC) */

#if defined(HAVE_HYPRE)

  /* Setting global options for HYPRE */
  /*----------------------------------*/

  /*! [sles_hypre_1] */
  {
    /* Initialization must be called before setting options;
       it does not need to be called before calling
       cs_sles_hypre_define(), as this is handled automatically. */

    /* No global options set yet... */
  }
  /*! [sles_hypre_1] */

  /* Setting pressure solver with hypre with Default PCG+BoomerAMG options */
  /*-----------------------------------------------------------------------*/

  /*! [sles_hypre_2] */
  {
    cs_sles_hypre_define(CS_F_(p)->id,
                         NULL,
                         CS_SLES_HYPRE_PCG,            /* solver type */
                         CS_SLES_HYPRE_BOOMERAMG,      /* preconditioner type */
                         NULL,
                         NULL);

  }
  /*! [sles_hypre_2] */

  /* Setting pressure solver with hypre on GPU and  user-defined options */
  /*---------------------------------------------------------------------*/

  /*! [sles_hypre_3] */
  {
    cs_sles_hypre_t *sc
      = cs_sles_hypre_define(CS_F_(p)->id,
                             NULL,
                             CS_SLES_HYPRE_PCG,            /* solver type */
                             CS_SLES_HYPRE_BOOMERAMG,      /* preconditioner type */
                             _hypre_p_setup_hook,
                             NULL);

    cs_sles_hypre_set_host_device(sc, 1);  /* run on GPU */
  }
  /*! [sles_hypre_3] */

#endif /* defined(HAVE_HYPRE) */

  /* Setting pressure solver with AMGX */
  /*-----------------------------------*/

#if defined(HAVE_AMGX)
  /*! [sles_amgx] */
  {
    cs_sles_amgx_t *amgx_p = cs_sles_amgx_define(CS_F_(p)->id, NULL);

    cs_sles_amgx_set_config_file(amgx_p, "PCG_CLASSICAL_V_JACOBI.json");
  }
  /*! [sles_amgx] */
#endif /* defined(HAVE_AMGX) */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
