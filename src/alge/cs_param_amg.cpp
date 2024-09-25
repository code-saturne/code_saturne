/*============================================================================
 * Routines to handle the set of parameters for algebraic multigrids (AMG)
 *============================================================================*/

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
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_amg.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the name of the type of algebraic multigrid (AMG)
 *
 * \param[in] type     type of AMG
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_amg_get_type_name(cs_param_amg_type_t  type)
{
  switch (type) {

  case CS_PARAM_AMG_NONE:
    return  "None";

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
    return  "Boomer V-cycle (Hypre)";
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    return  "Boomer W-cycle (Hypre)";
  case CS_PARAM_AMG_INHOUSE_K:
    return  "In-house (K-cycle)";
  case CS_PARAM_AMG_INHOUSE_V:
    return  "In-house (V-cycle)";
  case CS_PARAM_AMG_PETSC_GAMG_V:
    return  "GAMG V-cycle (PETSc)";
  case CS_PARAM_AMG_PETSC_GAMG_W:
    return  "GAMG W-cycle (PETSc)";
  case CS_PARAM_AMG_PETSC_HMG_V:
    return  "HMG V-cycle (PETSc)";
  case CS_PARAM_AMG_PETSC_HMG_W:
    return  "HMG W-cycle (PETSc)";
  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of AMG. Stop execution.", __func__);
  }

  return "";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the related solver class from the amg type
 *
 * \param[in] amg_type    type of AMG to consider
 *
 * \return the related solver class or CS_PARAM_SOLVER_CLASS_CS
 */
/*----------------------------------------------------------------------------*/

cs_param_solver_class_t
cs_param_amg_get_class(cs_param_amg_type_t  amg_type)
{
  switch (amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    return CS_PARAM_SOLVER_CLASS_HYPRE;

  case CS_PARAM_AMG_PETSC_GAMG_V:
  case CS_PARAM_AMG_PETSC_GAMG_W:
  case CS_PARAM_AMG_PETSC_HMG_V:
  case CS_PARAM_AMG_PETSC_HMG_W:
    return CS_PARAM_SOLVER_CLASS_PETSC;

  default:
    return CS_PARAM_SOLVER_CLASS_CS;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new structure storing a set of parameters used when calling
 *        boomerAMG. Set default values for all parameters.
 *
 * \return a pointer to a new set of boomerAMG parameters
 */
/*----------------------------------------------------------------------------*/

cs_param_amg_boomer_t *
cs_param_amg_boomer_create(void)
{
  cs_param_amg_boomer_t  *bamgp = nullptr;

  BFT_MALLOC(bamgp, 1, cs_param_amg_boomer_t);

  /* Main options */

  bamgp->coarsen_algo = CS_PARAM_AMG_BOOMER_COARSEN_HMIS;
  bamgp->coarse_solver = CS_PARAM_AMG_BOOMER_GAUSS_ELIM;

  /* From the HYPRE documentation: "There are further parameter choices for the
   * individual smoothers, which are described in the reference manual. The
   * default relaxation type is l1-Gauss-Seidel, using a forward solve on the
   * down cycle and a backward solve on the up-cycle, to keep symmetry. Note
   * that if BoomerAMG is used as a preconditioner for conjugate gradient, it
   * is necessary to use a symmetric smoother. Other symmetric options are
   * weighted Jacobi or hybrid symmetric Gauss-Seidel."
   */

  bamgp->n_down_iter = 1;
  bamgp->down_smoother = CS_PARAM_AMG_BOOMER_FORWARD_L1_GS;
  /* CS_PARAM_AMG_BOOMER_HYBRID_SSOR is also a good choice */

  bamgp->n_up_iter = 1;
  bamgp->up_smoother = CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS;
  /* CS_PARAM_AMG_BOOMER_HYBRID_SSOR is also a good choice */

  /* Advanced options */

  /* For best performance, it might be necessary to set certain parameters,
   * which will affect both coarsening and interpolation. One important
   * parameter is the strong threshold.  The default value is 0.25, which
   * appears to be a good choice for 2-dimensional problems and the low
   * complexity coarsening algorithms. For 3-dimensional problems a better
   * choice appears to be 0.5, when using the default coarsening
   * algorithm. However, the choice of the strength threshold is problem
   * dependent.
   */

  bamgp->strong_threshold = 0.5;
  bamgp->interp_algo = CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_CC;
  bamgp->p_max = 8;
  bamgp->n_agg_levels = 2;
  bamgp->n_agg_paths = 2;  /* HYPRE default is 1 */

  return bamgp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy the given set of parameters used when calling boomerAMG into a
 *        new structure
 *
 * \param[in] bamgp   reference set of boomerAMG parameters
 *
 * \return a pointer to a new set of boomerAMG parameters
 */
/*----------------------------------------------------------------------------*/

cs_param_amg_boomer_t *
cs_param_amg_boomer_copy(const cs_param_amg_boomer_t  *bamgp)
{
  cs_param_amg_boomer_t  *cpy = cs_param_amg_boomer_create();

  cpy->coarsen_algo = bamgp->coarsen_algo;
  cpy->coarse_solver = bamgp->coarse_solver;

  cpy->n_down_iter = bamgp->n_down_iter;
  cpy->down_smoother = bamgp->down_smoother;

  cpy->n_up_iter = bamgp->n_up_iter;
  cpy->up_smoother = bamgp->up_smoother;

  cpy->strong_threshold = bamgp->strong_threshold;
  cpy->interp_algo = bamgp->interp_algo;
  cpy->p_max = bamgp->p_max;
  cpy->n_agg_levels = bamgp->n_agg_levels;
  cpy->n_agg_paths = bamgp->n_agg_paths;

  return cpy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the name of the smoother used with BoomerAMG (HYPRE library)
 *
 * \param[in] smoother  smoother type
 *
 * \return name of the given smoother type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_amg_get_boomer_smoother_name(cs_param_amg_boomer_smoother_t  smoother)
{
  switch (smoother) {

  case CS_PARAM_AMG_BOOMER_JACOBI:
    return "Jacobi (0)";
  case CS_PARAM_AMG_BOOMER_FORWARD_GS:
    return "Forward Gauss-Seidel (3)";
  case CS_PARAM_AMG_BOOMER_BACKWARD_GS:
    return "Backward Gauss-Seidel (4)";
  case CS_PARAM_AMG_BOOMER_HYBRID_SSOR:
    return "Hybrid symmetric SOR (6)";
  case CS_PARAM_AMG_BOOMER_L1_SGS:
    return "L1 symmetric Gauss-Seidel (8)";
  case CS_PARAM_AMG_BOOMER_GAUSS_ELIM:
    return "Gauss elimination (9)";
  case CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS:
    return "Backward l1 Gauss-Seidel (13)";
  case CS_PARAM_AMG_BOOMER_FORWARD_L1_GS:
    return "Forward l1 Gauss-Seidel (14)";
  case CS_PARAM_AMG_BOOMER_CG:
    return "Conjugate Gradient (15)";
  case CS_PARAM_AMG_BOOMER_CHEBYSHEV:
    return "Chebyshev (16)";
  case CS_PARAM_AMG_BOOMER_FCF_JACOBI:
    return "FCF Jacobi (17)";
  case CS_PARAM_AMG_BOOMER_L1_JACOBI:
    return "L1 Jacobi (18)";

  default:
    return "Undefined";
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the set of parameters used for setting BoomerAMG
 *
 * \param[in] name      name related to the current SLES
 * \param[in] bamgp     set of boomerAMG parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_param_amg_boomer_log(const char                  *name,
                       const cs_param_amg_boomer_t  *bamgp)
{
  if (bamgp == nullptr)
    return;

  char  *prefix = nullptr;
  int  len = strlen(name) + strlen("  *  |") + 1;
  BFT_MALLOC(prefix, len, char);
  sprintf(prefix, "  * %s |", name);

  cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_down_smoothing: %1d it.| %s\n",
                prefix, bamgp->n_down_iter,
                cs_param_amg_get_boomer_smoother_name(bamgp->down_smoother));
  cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_up_smoothing:   %1d it.| %s\n",
                prefix, bamgp->n_up_iter,
                cs_param_amg_get_boomer_smoother_name(bamgp->up_smoother));
  cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarse_solver:  %s\n",
                prefix,
                cs_param_amg_get_boomer_smoother_name(bamgp->coarse_solver));

  switch (bamgp->coarsen_algo) {

  case CS_PARAM_AMG_BOOMER_COARSEN_FALGOUT:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarsening:    %s\n",
                  prefix, "Falgout (6)");
    break;

  case CS_PARAM_AMG_BOOMER_COARSEN_PMIS:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarsening:     %s\n",
                  prefix, "PMIS (8)");
    break;

  case CS_PARAM_AMG_BOOMER_COARSEN_HMIS:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarsening:     %s\n",
                  prefix, "HMIS (10)");
    break;

  case CS_PARAM_AMG_BOOMER_COARSEN_CGC:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarsening:     %s\n",
                  prefix, "CGC (21)");
    break;

  case CS_PARAM_AMG_BOOMER_COARSEN_CGC_E:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarsening:     %s\n",
                  prefix, "CGC-E (22)");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_coarsening:     %s\n",
                  prefix, "Unknown");
    break;

  } /* Coarsening algorithm */

  /* Advanced parameters */

  cs_log_printf(CS_LOG_SETUP, "%s   strong_threshold:       %f\n",
                prefix, bamgp->strong_threshold);

  cs_log_printf(CS_LOG_SETUP,
                "%s   aggressive_coarsening:  %d lv. | %d paths\n",
                prefix, bamgp->n_agg_levels, bamgp->n_agg_paths);

  cs_log_printf(CS_LOG_SETUP, "%s   Pmax set:               %d\n",
                prefix, bamgp->p_max);

  switch (bamgp->interp_algo) {

  case CS_PARAM_AMG_BOOMER_INTERP_HYPERBOLIC:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "For hyperbolic PDEs (2)");
    break;

  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_CC:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "extended+i (ext+i-cc) (6)");
    break;

  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "extended+i (ext+i) (7)");
    break;

  case CS_PARAM_AMG_BOOMER_INTERP_FF1:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "FF1 (13)");
    break;

  case CS_PARAM_AMG_BOOMER_INTERP_EXTENDED:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "extended (14)");
    break;

  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_MATRIX:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "extended+i in matrix form (ext+i-mm) (17)");
    break;

  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_E_MATRIX:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "extended+e in matrix form (ext+e-mm) (18)");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s BoomerAMG_interpolation:  %s\n",
                  prefix, "Unknown");
    break;

  } /* Interpolation algorithm */

  BFT_FREE(prefix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new structure storing a set of parameters used when calling
 *        the in-house AMG algo. Set default values for all parameters.
 *
 * \param[in] used_as_solver   true or false
 * \param[in] used_as_k_cycle  true or false
 *
 * \return a pointer to a new set of parameters
 */
/*----------------------------------------------------------------------------*/

cs_param_amg_inhouse_t *
cs_param_amg_inhouse_create(bool  used_as_solver,
                            bool  used_as_k_cycle)
{
  cs_param_amg_inhouse_t  *amgp = nullptr;

  BFT_MALLOC(amgp, 1, cs_param_amg_inhouse_t);

  /* Options shared among all configurations */

  amgp->max_levels = 15;
  amgp->min_n_g_rows = 256;

  amgp->coarse_rtol_mult = 1.0;
  amgp->coarse_max_iter = 2500;
  amgp->coarse_poly_degree = 0;
  amgp->coarse_solver = CS_PARAM_AMG_INHOUSE_CG;

  /* Down smoother options */

  amgp->down_poly_degree = 0;
  amgp->down_smoother = CS_PARAM_AMG_INHOUSE_CG;
  amgp->n_down_iter = 2;

  /* Up smoother options */

  amgp->up_poly_degree = 0;
  amgp->up_smoother = CS_PARAM_AMG_INHOUSE_CG;
  amgp->n_up_iter = 2;

  if (used_as_k_cycle) {

    /* Coarsening options */

    amgp->coarsen_algo = CS_PARAM_AMG_INHOUSE_COARSEN_SPD_PW;
    amgp->p0p1_relax = 0.;

    if (used_as_solver) {

      amgp->aggreg_limit = 4;

      /* Down smoother options */

      amgp->down_poly_degree = -1;
      amgp->down_smoother = CS_PARAM_AMG_INHOUSE_PROCESS_SGS;
      amgp->n_down_iter = 3;

      /* Up smoother options */

      amgp->up_poly_degree = -1;
      amgp->up_smoother = CS_PARAM_AMG_INHOUSE_PROCESS_SGS;
      amgp->n_up_iter = 3;

    }
    else { /* Used as a preconditioner */

      amgp->aggreg_limit = 8; /* More aggresive */

    }

  }
  else { /* This is a V-cycle AMG algorithm */

    /* Coarsening options */

    amgp->aggreg_limit = 3;
    amgp->coarsen_algo = CS_PARAM_AMG_INHOUSE_COARSEN_SPD_MX;
    amgp->p0p1_relax = 0.95;

    if (used_as_solver) {

      /* Down smoother options */

      amgp->n_down_iter = 5;

      /* Up smoother options */

      amgp->n_up_iter = 5;

    }

  } /* V-cycle */

  return amgp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy the given set of parameters used when calling in-house AMG algo.
 *        into a new structure
 *
 * \param[in] amgp  reference set of in-house AMG parameters
 *
 * \return a pointer to a new set of in-house AMG parameters
 */
/*----------------------------------------------------------------------------*/

cs_param_amg_inhouse_t *
cs_param_amg_inhouse_copy(const cs_param_amg_inhouse_t  *amgp)
{
  cs_param_amg_inhouse_t  *cpy = cs_param_amg_inhouse_create(true, true);

  cpy->max_levels = amgp->max_levels;
  cpy->min_n_g_rows = amgp->min_n_g_rows;
  cpy->p0p1_relax = amgp->p0p1_relax;

  cpy->aggreg_limit = amgp->aggreg_limit;
  cpy->coarsen_algo = amgp->coarsen_algo;

  cpy->coarse_solver = amgp->coarse_solver;
  cpy->coarse_max_iter = amgp->coarse_max_iter;
  cpy->coarse_poly_degree = amgp->coarse_poly_degree;
  cpy->coarse_rtol_mult = amgp->coarse_rtol_mult;

  cpy->n_down_iter = amgp->n_down_iter;
  cpy->down_smoother = amgp->down_smoother;
  cpy->down_poly_degree = amgp->down_poly_degree;

  cpy->n_up_iter = amgp->n_up_iter;
  cpy->up_smoother = amgp->up_smoother;
  cpy->up_poly_degree = amgp->up_poly_degree;

  return cpy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the name of the solver used with in-house AMG algo.
 *
 * \param[in] solver  solver type
 *
 * \return name of the given solver type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_amg_get_inhouse_solver_name(cs_param_amg_inhouse_solver_t  solver)
{
  switch (solver) {

  case CS_PARAM_AMG_INHOUSE_FORWARD_GS:
    return "Forward Gauss-Seidel smoother";
  case CS_PARAM_AMG_INHOUSE_BACKWARD_GS:
    return "Backward Gauss-Seidel smoother";
  case CS_PARAM_AMG_INHOUSE_JACOBI:
    return "Jacobi";
  case CS_PARAM_AMG_INHOUSE_PROCESS_GS:
    return "Process local Gauss-Seidel/Jacobi";
  case CS_PARAM_AMG_INHOUSE_PROCESS_SGS:
    return "Process local Symm. Gauss-Seidel/Jacobi";
  case CS_PARAM_AMG_INHOUSE_CG:
    return "Conjugate Gradient";
  case CS_PARAM_AMG_INHOUSE_CR3:
    return "3-layer Conjugate Residual (CR3)";
  case CS_PARAM_AMG_INHOUSE_GCR:
    return "Generalized Conjugate Residual (GCR)";
  case CS_PARAM_AMG_INHOUSE_GMRES:
    return "GMRES";

  default:
    return "Undefined";

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the set of parameters used for setting in-house AMG algorithms
 *
 * \param[in] name  name related to the current SLES
 * \param[in] amgp  set of in-house AMG parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_param_amg_inhouse_log(const char                    *name,
                         const cs_param_amg_inhouse_t  *amgp)
{
  if (amgp == nullptr)
    return;

  char  *prefix = nullptr;
  int  len = strlen(name) + strlen("  *  |") + 1;
  BFT_MALLOC(prefix, len, char);
  sprintf(prefix, "  * %s |", name);

  cs_log_printf(CS_LOG_SETUP,
                "%s in-house.AMG_down_smoothing: %1d it.| poly.degree %1d |"
                " %s\n",
                prefix, amgp->n_down_iter, amgp->down_poly_degree,
                cs_param_amg_get_inhouse_solver_name(amgp->down_smoother));
  cs_log_printf(CS_LOG_SETUP,
                "%s in-house.AMG_up_smoothing:   %1d it.| poly.degree %1d |"
                " %s\n",
                prefix, amgp->n_up_iter, amgp->up_poly_degree,
                cs_param_amg_get_inhouse_solver_name(amgp->up_smoother));
  cs_log_printf(CS_LOG_SETUP,
                "%s in-house.AMG_coarse_solver:  poly.degree %1d | %s\n",
                prefix, amgp->coarse_poly_degree,
                cs_param_amg_get_inhouse_solver_name(amgp->coarse_solver));
  cs_log_printf(CS_LOG_SETUP,
                "%s                              rtol.mult   %.2e\n",
                prefix, amgp->coarse_rtol_mult);

  switch (amgp->coarsen_algo) {

  case CS_PARAM_AMG_INHOUSE_COARSEN_SPD_DX:
    cs_log_printf(CS_LOG_SETUP, "%s in-house.AMG_coarsening:    %s\n",
                  prefix, "SPD, diag/extradiag ratio-based");
    break;
  case CS_PARAM_AMG_INHOUSE_COARSEN_SPD_MX:
    cs_log_printf(CS_LOG_SETUP, "%s in-house.AMG_coarsening:    %s\n",
                  prefix, "SPD, max diag/extradiag ratio-based");
    break;
  case CS_PARAM_AMG_INHOUSE_COARSEN_SPD_PW:
    cs_log_printf(CS_LOG_SETUP, "%s in-house.AMG_coarsening:    %s\n",
                  prefix, "SPD, pairwise aggregation");
    break;
  case CS_PARAM_AMG_INHOUSE_COARSEN_CONV_DIFF_DX:
    cs_log_printf(CS_LOG_SETUP, "%s in-house.AMG_coarsening:    %s\n",
                  prefix, "Conv-Diff, diag/extradiag ratio-based");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s in-house.AMG_coarsening:     %s\n",
                  prefix, "Unknown");
    break;

  } /* Coarsening algorithm */

  /* Advanced parameters related to the coarsening */

  cs_log_printf(CS_LOG_SETUP, "%s   max_levels:             %d\n",
                prefix, amgp->max_levels      );
  cs_log_printf(CS_LOG_SETUP, "%s   min_n_g_rows:           %lu\n",
                prefix, amgp->min_n_g_rows);
  cs_log_printf(CS_LOG_SETUP, "%s   p0p1_relax:             %f\n",
                prefix, amgp->p0p1_relax);
  cs_log_printf(CS_LOG_SETUP, "%s   aggregation_limit:      %d\n",
                prefix, amgp->aggreg_limit);

  BFT_FREE(prefix);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
