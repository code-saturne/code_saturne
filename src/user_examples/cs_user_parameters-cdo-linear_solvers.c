/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#if defined(HAVE_MUMPS)
#include <dmumps_c.h>
#include <smumps_c.h>
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

#if defined(HAVE_MUMPS)
#include "cs_sles_mumps.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t    *domain)
{
  /*
     Example: Use MUMPS to solve the saddle-point problem arising from CDO
     -------  schemes for (Navier-)Stokes
  */

  /*! [cdo_sles_navsto_full_mumps] */
  {
    /* Parameters related to the Navier-Stokes settings. General strategy. */

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    cs_navsto_param_set(nsp, CS_NSKEY_SLES_STRATEGY, "mumps");

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_navsto_full_mumps] */

  /*! [cdo_sles_navsto_alu_mumps] */
  {
    /* Parameters related to the Navier-Stokes settings. General strategy. */

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    cs_navsto_param_set(nsp, CS_NSKEY_SLES_STRATEGY, "alu");
    cs_navsto_param_set(nsp, CS_NSKEY_GD_SCALE_COEF, "5e3");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_RTOL, "1e-8");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_ATOL, "1e-14");

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_navsto_alu_mumps] */

  /*! [cdo_sles_navsto_gkb_mumps] */
  {
    /* Parameters related to the Navier-Stokes settings. General strategy. */

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    cs_navsto_param_set(nsp, CS_NSKEY_SLES_STRATEGY, "gkb");
    cs_navsto_param_set(nsp, CS_NSKEY_GD_SCALE_COEF, "1e3");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_RTOL, "1e-8");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_ATOL, "1e-14");

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_navsto_gkb_mumps] */

  /*! [cdo_sles_navsto_gkb_kcycle] */
  {
    /* Parameters related to the Navier-Stokes settings. General strategy. */

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    cs_navsto_param_set(nsp, CS_NSKEY_SLES_STRATEGY, "gkb");
    cs_navsto_param_set(nsp, CS_NSKEY_GD_SCALE_COEF, "0");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_RTOL, "1e-8");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_ATOL, "1e-14");

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");

    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND_BLOCK_TYPE, "none");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");

    /* Tolerance for the inner solver */

    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL_EPS, "1e-5");
  }
  /*! [cdo_sles_navsto_gkb_kcycle] */

  /*! [cdo_sles_navsto_uzacg] */
  {
    /* Parameters related to the Navier-Stokes settings. General strategy. */

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    cs_navsto_param_set(nsp, CS_NSKEY_SLES_STRATEGY, "uzawa_cg");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_RTOL, "1e-6");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_ATOL, "1e-14");

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");

    /* Set the inner solver for the velocity block */

#if defined(HAVE_HYPRE)
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "fgmres");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL_EPS, "1e-1");
#else
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL_EPS, "1e-4");
#endif

  /*===============================
    Set the Schur complement solver
    =============================== */

    /* Available approximations are:
     *
     *  CS_PARAM_SCHUR_DIAG_INVERSE,
     *  CS_PARAM_SCHUR_LUMPED_INVERSE,
     *  CS_PARAM_SCHUR_MASS_SCALED,
     */

    nsp->sles_param->schur_approximation = CS_PARAM_SCHUR_MASS_SCALED;

  }
  /*! [cdo_sles_navsto_uzacg] */

  /*! [cdo_sles_navsto_minres] */
  {
    /* Parameters related to the Stokes settings.
       MINRES is not possible with a non-symmetric saddle-point system
       General strategy. */

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    cs_navsto_param_set(nsp, CS_NSKEY_SLES_STRATEGY, "minres");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_RTOL, "1e-9");
    cs_navsto_param_set(nsp, CS_NSKEY_IL_ALGO_ATOL, "1e-14");

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");

    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");

#if defined(HAVE_PETSC)
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND_BLOCK_TYPE, "diag");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL_EPS, "1e-1");

    /* Must be set after the previous line to switch to PETSC in order to be
       able to use a block preconditioning for the velocity block */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");

#else  /* PETSc not installed */
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL_EPS, "1e-4");
#endif

  /*===============================
    Set the Schur complement solver
    =============================== */

    /* Available approximations are:
     *
     *  CS_PARAM_SCHUR_DIAG_INVERSE,
     *  CS_PARAM_SCHUR_LUMPED_INVERSE,
     *  CS_PARAM_SCHUR_MASS_SCALED  --> Good choice for the Stokes eq.
     */

    nsp->sles_param->schur_approximation = CS_PARAM_SCHUR_MASS_SCALED;
  }
  /*! [cdo_sles_navsto_minres] */
}

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
  /*! [linear_solver_immediate_exit] */
  /* Redefine the threshold under which an immediate exit occurs */

  cs_sles_set_epzero(1e-15);
  /*! [linear_solver_immediate_exit] */

  /*! [linear_solver_kcycle] */
  {
    /* Retrieve the set of SLES parameters for the "momentum equation */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("momentum");
    cs_param_sles_t  *slesp = eqp->sles_param;
    assert(slesp->field_id > -1);

    /* In case of an in-house K-cylcle multigrid as a preconditioner of a
       linear iterative solver */

    if (slesp->precond == CS_PARAM_PRECOND_AMG) {

      /* If multigrid is chosen as preconditioner */

      if (slesp->amg_type == CS_PARAM_AMG_HOUSE_K) {

        /* If this is a K-cycle multigrid */

        /* Retrieve the different context structures to modify/apply additional
           settings */

        cs_sles_t  *sles = cs_sles_find_or_add(slesp->field_id, NULL);
        cs_sles_it_t  *itsol = cs_sles_get_context(sles);
        cs_sles_pc_t  *pc = cs_sles_it_get_pc(itsol);
        cs_multigrid_t  *mg = NULL;

        if (itsol == NULL) { /* Not defined yet. */

          if (slesp->solver == CS_PARAM_ITSOL_CG ||
              slesp->solver == CS_PARAM_ITSOL_FCG)
            itsol =  cs_sles_it_define(slesp->field_id, NULL,
                                       CS_SLES_IPCG, -1, slesp->n_max_iter);
          else
            bft_error(__FILE__, __LINE__, 0,
                      " %s: Case not treated.\n", __func__);

        }
        assert(itsol != NULL);

        if (pc == NULL) { /* Not defined yet */

          pc = cs_multigrid_pc_create(CS_MULTIGRID_K_CYCLE);
          mg = cs_sles_pc_get_context(pc);
          cs_sles_it_transfer_pc(itsol, &pc);

        }
        else
          mg = cs_sles_pc_get_context(pc);

        assert(mg != NULL && pc != NULL);

        cs_multigrid_set_solver_options
          (mg,
           CS_SLES_TS_F_GAUSS_SEIDEL,
           CS_SLES_TS_B_GAUSS_SEIDEL,
           CS_SLES_PCG,           /* coarse solver */
           1,                     /* n_max_cycles */
           1,                     /* n_max_iter_descent, */
           1,                     /* n_max_iter_ascent */
           200,                   /* n_max_iter_coarse */
           -1,                    /* poly_degree_descent */
           -1,                    /* poly_degree_ascent */
           1,                     /* poly_degree_coarse */
           -1.0,                  /* precision_mult_descent */
           -1.0,                  /* precision_mult_ascent */
           1.0);                  /* precision_mult_coarse */

        /* Available settings:
         * - max. number of elements in an aggregation
         * - type of algorithm to perform the aggregation
         * - max. number of levels (i.e. grids)
         * - max global number of rows at the coarsest level
         * - type of relaxation (weighting between a P_0 and P_1). For K-cycle,
         *   this should be equal to 0.
         * - Activation of the postprocessing for the aggregation if > 0.
         * Aggregation set is numbered by its coarse row number modulo this
         * value
         */

        cs_multigrid_set_coarsening_options(mg,
                                            8,    /* aggregation_limit*/
                                            CS_GRID_COARSENING_SPD_PW,
                                            10,   /* n_max_levels */
                                            500,  /* min_g_cells */
                                            0.,   /* P0P1 relaxation */
                                            0);   /* postprocess */

      } /* K-cycle */

    } /* Multigrid as preconditioner */
  }
  /*! [linear_solver_kcycle] */
}

#if defined(HAVE_MUMPS)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for advanced user settings of a MUMPS solver.
 *        This function is called two times during the setup stage.
 *        1. Before the analysis step
 *        2. Before the factorization step
 *
 * One can recover the MUMPS step through the "job" member.
 * MUMPS_JOB_ANALYSIS or MUMPS_JOB_FACTORIZATION
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that structure should
 * not be temporary (i.e. local);
 *
 * \param[in]      slesp    pointer to the related cs_param_sles_t structure
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] pmumps   pointer to DMUMPS_STRUC_C or SMUMPS_STRUC_C struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_sles_mumps_hook(const cs_param_sles_t   *slesp,
                        void                    *context,
                        void                    *pmumps)
{
  CS_UNUSED(slesp);
  CS_UNUSED(context);

  /*! [sles_mumps_advanced_hook] */

  /* If MUMPS in single-precision is used, one has to modify the declaration as
     follows:

     SMUMPS_STRUC_C  *mumps = pmumps;

     All the remaining settings apply to the mumps variable in the single or
     double-precision case.
   */

  DMUMPS_STRUC_C  *mumps = pmumps;
  assert(mumps != NULL);

  int  ir_steps = 0;
  bool  use_parmetis = false;
  bool  use_openmp = true;
  bool  static_pivoting = true;
  double  blr_rate = -1; /* If <= 0 -> no BLR */

  if (use_parmetis) {

    mumps->ICNTL(28) = 2; // Analysis parallel
    mumps->ICNTL(29) = 2; // Choix de parmetis comme renum parall

  }
  else {

    /* Choose the way numbering is performed inside MUMPS.
     * This option may have a strong effect on the elapsed time
     * 0: AMD
     * 3: Scotch (need a MUMPS library compiled with Scotch)
     * 4: PORD
     * 5: METIS (need a MUMPS library compiled with METIS)
     * 6: QAMD
     * 7: automatic choice done by MUMPS
     */

     mumps->ICNTL(7) = 6;

  }

  /* Clustering for the analysis (Only for vector-valued system) */

  mumps->ICNTL(15) = -3; // Number of terms to be clustered

  /* Advanced settings for openMP */

  if (use_openmp && cs_glob_n_threads > 1) {

    mumps->KEEP(401) = 1; // Activate agressive openMP (comment when omp=1)
    mumps->KEEP(370) = 1; // Useful for hybrid MPI/openMP
    mumps->KEEP(371) = 1; // Useful for hybrid MPI/openMP

  }

  /* Iterative refinement */

  if (ir_steps < 0)
    mumps->ICNTL(10) = ir_steps;

  /* BLR compression */

  if (blr_rate > 0) {

    mumps->ICNTL(35) = 2; // Activate BLR algo (facto + solve)
    mumps->ICNTL(36) = 1; // Variante de BLR (0 is also a good choice)
    mumps->ICNTL(37) = 0; // Memory compression (= 1) but time consumming
    mumps->ICNTL(40) = 0; // Memory compression (= 1) mixed precision
    mumps->CNTL(7) = blr_rate; // Compression rate

  }

  if (static_pivoting)
    mumps->CNTL(4) = 0.0;
  /*! [sles_mumps_advanced_hook] */
}
#endif  /* HAVE_MUMPS */

/*----------------------------------------------------------------------------*/

END_C_DECLS
