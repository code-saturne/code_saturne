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

#if defined(HAVE_MUMPS)
#include <dmumps_c.h>
#include <smumps_c.h>
#endif

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
  CS_NO_WARN_IF_UNUSED(domain);

  /*! [cdo_sles_solver_family] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEqName");

    /* Specify that this is an in-house solvers (default choice) */

    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "saturne");

    /* Use solver/preconditioner available in the MUMPS library */

#if defined(HAVE_MUMPS)
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif

    /* Use solver/preconditioner available in the PETSc library */

#if defined(HAVE_PETSC)
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: PETSc is not available\n", __func__);
#endif

    /* Use solver/preconditioner available in the HYPRE library */

#if defined(HAVE_HYPRE)
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "hypre");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: HYPRE is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_solver_family] */

  /*! [cdo_sles_user_simple] */
  {
    /*
      Example: Use a BiCGstab iterative solver with a 1st order Neumann
               polynomial preconditioner is used to solve a user-defined
               equation
    */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEqName");

    cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "bicgs");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "poly1");
  }
  /*! [cdo_sles_user_simple] */

  /*! [cdo_sles_allow_no_op] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEqName");

    /* In case of a in-house solver, allow one to skip the solve in some
       specific situation rhs and solution nearly equal to zero */

    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_NO_OP, "true");
  }
  /*! [cdo_sles_allow_no_op] */

  /*! [cdo_sles_user_mumps] */
  {
    /*
      Example: Use MUMPS to solve a user-defined equation
    */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");

#if defined(HAVE_MUMPS)
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_user_mumps] */

  /*! [cdo_sles_user_mumps_advanced] */
  {
    /* Parameters related to a user-defined equation */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");
    cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(eqp);

    cs_param_sles_mumps(slesp,
                        false,                     /* single-precision ? */
                        CS_PARAM_MUMPS_FACTO_LU);  /* type of factorization */

    cs_param_sles_mumps_advanced(slesp,
                                 CS_PARAM_MUMPS_ANALYSIS_QAMD,
                                 1,  /* size of the block for analysis */
                                 -1, /* pct memory increase < 0 --> not used */
                                 0,  /* BLR compression: 0 --> not used */
                                 1,  /* iterative refinement steps */
                                 CS_PARAM_MUMPS_MEMORY_AUTO,
                                 true); /* advanced optimizations */
  }
  /*! [cdo_sles_user_mumps_advanced] */

  /*
     Example: Use MUMPS to solve the saddle-point problem arising from CDO
     -------  schemes for (Navier-)Stokes
  */

  /*! [cdo_sles_navsto_full_mumps] */
  {
    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");

#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "mumps");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_navsto_full_mumps] */

  /*! [cdo_sles_navsto_alu_mumps] */
  {
    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "alu");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_AUGMENT_SCALING, "1e3");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-8");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "mumps");

    /* More advanced usage */

    cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(mom_eqp);

    cs_param_sles_mumps(slesp,
                        false,  /* single-precision ? */
                        CS_PARAM_MUMPS_FACTO_LU);

    cs_param_sles_mumps_advanced(slesp,
                                 CS_PARAM_MUMPS_ANALYSIS_AUTO,
                                 3,     /* size of the block for analysis */
                                 -1,    /* pct memory increase < 0 = not used */
                                 0,     /* BLR compression:  0 = not used */
                                 0,     /* iterative refinement steps */
                                 CS_PARAM_MUMPS_MEMORY_AUTO, /* memory usage */
                                 true); /* advanced optimizations */
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_navsto_alu_mumps] */

  /*! [cdo_sles_navsto_gkb_mumps] */
  {
    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings for the saddle-point system */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "gkb");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_AUGMENT_SCALING, "10");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-8");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

    /* Linear algebra settings for the (1,1) block */

#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [cdo_sles_navsto_gkb_mumps] */

  /*! [cdo_sles_navsto_gkb_kcycle] */
  {
    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings for the saddle-point system */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "gkb");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_AUGMENT_SCALING, "0");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-8");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

    /* Linear algebra settings for the (1,1)-block */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-5");
  }
  /*! [cdo_sles_navsto_gkb_kcycle] */

  /*! [cdo_sles_navsto_uzacg] */
  {
    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings for the saddle-point system */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "uzawa_cg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-6");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

    /* Linear algebra settings for the (1,1)-block */

#if defined(HAVE_HYPRE)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "fgmres");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-1");
#else
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-4");
#endif

    /* Linear algebra settings for the Schur complement approximation */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SCHUR_APPROX, "mass_scaled");
  }
  /*! [cdo_sles_navsto_uzacg] */

  /*! [cdo_sles_navsto_minres] */
  {
    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");
    cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(mom_eqp);

    /* Parameters related to the Stokes settings.
       MINRES is not possible with a non-symmetric saddle-point system */

    /* Linear algebra settings for the saddle-point system */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SLES_VERBOSITY, "2");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "minres");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-9");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

    /* Linear algebra settings for the (1,1)-block */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");

#if defined(HAVE_PETSC) /* One assumes that PETSc is installed with Hypre */
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-1");

    /* Must be set after the previous line to switch to PETSC in order to be
       able to use a block preconditioning for the velocity block */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND_BLOCK_TYPE, "diag");

    /* Set the main parameters for BoomerAMG (Please refer to the Hypre
       documentation for more details) */

    cs_param_sles_boomeramg(slesp,
                            /* n_down_iter, down smoother */
                            1, CS_PARAM_AMG_BOOMER_FORWARD_L1_GS,
                            /* n_up_iter, up smoother */
                            1, CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS,
                            /* coarse solver */
                            CS_PARAM_AMG_BOOMER_GAUSS_ELIM,
                            /* coarsening algorithm */
                            CS_PARAM_AMG_BOOMER_COARSEN_HMIS);

    /* Set advanced parameters for BoomerAMG */

    cs_param_sles_boomeramg_advanced(slesp,
                                     0.5, /* strong threshold */
                                     CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_CC,
                                     8,   /* Pmax */
                                     2,   /* n_agg_levels */
                                     2);  /* n_agg_paths */
#else  /* PETSc not installed */
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-4");

    cs_param_sles_amg_inhouse(slesp,
                              /* n_down_iter, down smoother, down poly. deg. */
                              1, CS_PARAM_AMG_INHOUSE_FORWARD_GS, -1,
                              /* n_up_iter, up smoother, up poly. deg. */
                              1, CS_PARAM_AMG_INHOUSE_BACKWARD_GS, -1,
                              /* coarse solver, coarse poly. deg. */
                              CS_PARAM_AMG_INHOUSE_CG, 1,
                              /* coarsen algo, aggregation limit */
                              CS_PARAM_AMG_INHOUSE_COARSEN_SPD_PW, 8);

    /* Set advanced parameters for an in-house AMG */

    cs_param_sles_amg_inhouse_advanced(slesp,
                                       CS_CDO_KEEP_DEFAULT, // max. levels
                                       100,                 // min. n_g_rows
                                       CS_CDO_KEEP_DEFAULT, // p0p1 relax.
                                       CS_CDO_KEEP_DEFAULT, // coarse max. iter
                                       1e-2);               // coarse rtol mult.
#endif

    /* Linear algebra settings for the Schur complement approximation */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SCHUR_APPROX, "mass_scaled");
  }
  /*! [cdo_sles_navsto_minres] */

  /*! [cdo_sles_boomer] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");

    cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "boomer");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RTOL, "1e-1");

    /* Set the main parameters of the BoomerAMG algorithm */

    cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(eqp);

    cs_param_sles_boomeramg(slesp,
                            /* n_down_iter, down smoother */
                            2, CS_PARAM_AMG_BOOMER_HYBRID_SSOR,
                            /* n_up_iter, up smoother */
                            2, CS_PARAM_AMG_BOOMER_HYBRID_SSOR,
                            /* coarse solver */
                            CS_PARAM_AMG_BOOMER_GAUSS_ELIM,
                            /* coarsening algorithm */
                            CS_PARAM_AMG_BOOMER_COARSEN_PMIS);

    /* Set advanced parameters for BoomerAMG */

    cs_param_sles_boomeramg_advanced(slesp,
                                     0.5, /* strong threshold */
                                     CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_CC,
                                     8,   /* Pmax */
                                     2,   /* n_agg_levels */
                                     2);  /* n_agg_paths */
  }
  /*! [cdo_sles_boomer] */
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
  {
    /* Redefine the threshold under which an immediate exit occurs */

    cs_sles_set_epzero(1e-15);
  }
  /*! [linear_solver_immediate_exit] */
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

/*! [mumps_user_hook] */
void
cs_user_sles_mumps_hook(const cs_param_sles_t   *slesp,
                        void                    *context,
                        void                    *pmumps)
{
  CS_UNUSED(slesp);
  CS_UNUSED(context);

  DMUMPS_STRUC_C  *mumps = (DMUMPS_STRUC_C *)pmumps;
  assert(mumps != nullptr);

  /* If MUMPS is used in single-precision, one has to modify the declaration as
     follows:

     SMUMPS_STRUC_C  *mumps = pmumps;

     All the remaining settings are identical in the single or double-precision
     case.
   */

  mumps->CNTL(4) = 0.0;    /* Static pivoting */

  mumps->ICNTL(58) = 2;    /* Symbolic factorization {0, 1, or 2}*/
}
/*! [mumps_user_hook] */
#endif  /* HAVE_MUMPS */

/*----------------------------------------------------------------------------*/

END_C_DECLS
