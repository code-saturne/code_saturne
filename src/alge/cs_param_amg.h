#ifndef __CS_PARAM_AMG_H__
#define __CS_PARAM_AMG_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param_types.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_param_amg.h

  \brief Routines to handle the set of parameters for algebraic multigrids
         (AMG) like boomerAMG of the HYPRE library, GAMG of the PETSc library
         or the Notay's K-cycle (in-house implementation) for instance

*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \enum cs_param_amg_type_t
 * Type of AMG (Algebraic MultiGrid) algorithm to use (either as a
 * preconditioner with or a solver). There are different choices of
 * implementation and of type of cycle
 */

typedef enum {

  CS_PARAM_AMG_NONE,            /*!< No specified algorithm */

  CS_PARAM_AMG_HYPRE_BOOMER_V,  /*!< V-cycle Boomer algorithm (Hypre lib.) */
  CS_PARAM_AMG_HYPRE_BOOMER_W,  /*!< W-cycle Boomer algorithm (Hypre lib.) */
  CS_PARAM_AMG_INHOUSE_K,       /*!< In-house algorithm with K-cycle */
  CS_PARAM_AMG_INHOUSE_V,       /*!< In-house algorithm with V-cycle */
  CS_PARAM_AMG_PETSC_GAMG_V,    /*!< V-cycle GAMG algorithm (PETSc lib.) */
  CS_PARAM_AMG_PETSC_GAMG_W,    /*!< W-cycle GAMG algorithm (PETSc lib.) */
  CS_PARAM_AMG_PETSC_HMG_V,     /*!< V-cycle HMG (hybrid AMG) from PETSc */
  CS_PARAM_AMG_PETSC_HMG_W,     /*!< W-cycle HMG (hybrid AMG) from PETSc */

  CS_PARAM_N_AMG_TYPES

} cs_param_amg_type_t;


/*! \enum cs_param_amg_boomer_coarsen_algo_t
 *  \brief Type of algorithm used in boomerAMG to coarsen a level. Only a
 *  selection of algorithms is available here. Values are those given in HYPRE
 */

typedef enum {

  CS_PARAM_AMG_BOOMER_COARSEN_FALGOUT = 6,
  CS_PARAM_AMG_BOOMER_COARSEN_PMIS = 8,
  CS_PARAM_AMG_BOOMER_COARSEN_HMIS = 10, /* (default) */
  CS_PARAM_AMG_BOOMER_COARSEN_CGC = 21,
  CS_PARAM_AMG_BOOMER_COARSEN_CGC_E = 22,

  CS_PARAM_AMG_BOOMER_N_COARSEN_ALGOS

} cs_param_amg_boomer_coarsen_algo_t;

/*! \enum cs_param_amg_boomer_interp_algo_t
 *  \brief Type of algorithm used in boomerAMG to coarsen a level. Only a
 *  selection of algorithms is available here. Values are those given in HYPRE
 */

typedef enum {

  CS_PARAM_AMG_BOOMER_INTERP_HYPERBOLIC = 2,
  CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_CC = 6, /* (default) Also for GPU */
  CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I = 7,
  CS_PARAM_AMG_BOOMER_INTERP_FF1 = 13,
  CS_PARAM_AMG_BOOMER_INTERP_EXTENDED = 14,     /* Also for GPU */
  CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_MATRIX = 17,
  CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_E_MATRIX = 18,

  CS_PARAM_AMG_BOOMER_N_INTERP_ALGOS

} cs_param_amg_boomer_interp_algo_t;

/*! \enum cs_param_amg_boomer_smoother_t
 *  \brief Type of algorithm used in boomerAMG to smooth a level. Only a
 *  selection of algorithms is available here. Values are those used in HYPRE.
 */

typedef enum {

  CS_PARAM_AMG_BOOMER_JACOBI = 0,
  CS_PARAM_AMG_BOOMER_FORWARD_GS = 3,
  CS_PARAM_AMG_BOOMER_BACKWARD_GS = 4,
  CS_PARAM_AMG_BOOMER_HYBRID_SSOR = 6,
  CS_PARAM_AMG_BOOMER_L1_SGS = 8,
  CS_PARAM_AMG_BOOMER_GAUSS_ELIM = 9,      /* for the coarsest level only */
  CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS = 13, /* (default) */
  CS_PARAM_AMG_BOOMER_FORWARD_L1_GS = 14,  /* (default) */
  CS_PARAM_AMG_BOOMER_CG = 15,
  CS_PARAM_AMG_BOOMER_CHEBYSHEV = 16,
  CS_PARAM_AMG_BOOMER_FCF_JACOBI = 17,
  CS_PARAM_AMG_BOOMER_L1_JACOBI = 18,

  CS_PARAM_AMG_BOOMER_N_SMOOTHERS

} cs_param_amg_boomer_smoother_t;

/*! \struct cs_param_amg_boomer_t
 *  \brief Set of the main parameters to setup the algebraic multigrid
 *         BoomerAMG belonging to the HYPRE library. These parameters are used
 *         to define this AMG directly in HYPRE or through the PETSc library
 *         according to the settings and the installed dependencies. Please
 *         refer to the HYPRE documentation for more details.
 */

typedef struct {

  /* Read the function \ref _petsc_pchypre_hook or \ref _hypre_boomeramg_hook
     for more details and read the HYPRE user guide */

  double                               strong_threshold;
  cs_param_amg_boomer_coarsen_algo_t   coarsen_algo;
  cs_param_amg_boomer_interp_algo_t    interp_algo;
  int                                  p_max;
  int                                  n_agg_levels;
  int                                  n_agg_paths;

  cs_param_amg_boomer_smoother_t       down_smoother;
  cs_param_amg_boomer_smoother_t       up_smoother;
  cs_param_amg_boomer_smoother_t       coarse_solver;

  int                                  n_down_iter;
  int                                  n_up_iter;

} cs_param_amg_boomer_t;

/* In-house AMG algorithms */
/* ----------------------- */

/*! \enum cs_param_amg_inhouse_solver_t
 *  \brief Type of algorithm used in the in-house algorithm for smoothing each
 *  level or solving the coarse level. Only the most relevant algorithms are
 *  available here.
 *
 * This enum avoids using \ref cs_sles_it_type_t which includes higher level
 * headers leading to loop cycle issues.
 */

 typedef enum {

   CS_PARAM_AMG_INHOUSE_FORWARD_GS = 1,   /* smoother only */
   CS_PARAM_AMG_INHOUSE_BACKWARD_GS = 2,  /* smoother only */

   CS_PARAM_AMG_INHOUSE_JACOBI = 3,
   CS_PARAM_AMG_INHOUSE_PROCESS_GS = 4,
   CS_PARAM_AMG_INHOUSE_PROCESS_SGS = 5,

   CS_PARAM_AMG_INHOUSE_CG = 6,
   CS_PARAM_AMG_INHOUSE_CR3 = 7,
   CS_PARAM_AMG_INHOUSE_GCR = 8,
   CS_PARAM_AMG_INHOUSE_GMRES = 9,

   CS_PARAM_AMG_INHOUSE_N_SOLVERS

 } cs_param_amg_inhouse_solver_t;

/*! \enum cs_param_amg_inhouse_coarsen_t
 *  \brief Type of algorithm used in the in-house algorithm to coarsen each
 *  level. This enum avoids using the associated \ref cs_grid_coarsening_t type
 *  which includes higher level headers.
 */

 typedef enum {

   /* For symmetric positive definite matrices (SPD) */

   CS_PARAM_AMG_INHOUSE_COARSEN_SPD_DX = 1 ,  /*!< SPD, diag/extradiag ratio
                                                based */
   CS_PARAM_AMG_INHOUSE_COARSEN_SPD_MX = 2,   /*!< SPD, diag/extradiag ratio
                                                based */
   CS_PARAM_AMG_INHOUSE_COARSEN_SPD_PW = 3 ,  /*!< SPD, pairwise aggregation */

   CS_PARAM_AMG_INHOUSE_COARSEN_CONV_DIFF_DX = 4,  /*!< for general matrices */

   CS_PARAM_AMG_INHOUSE_N_COARSENINGS

 } cs_param_amg_inhouse_coarsen_t;

/*! \struct cs_param_amg_inhouse_t
 *  \brief Set of the main parameters used to setup the algebraic multigrid
 *         available natively in code_saturne (in-house implementations). These
 *         parameters are the most impacting ones. For a more advanced
 *         usage, this is still possible to consider the function \ref
 *         cs_user_linear_solvers
 */

typedef struct {

  /* Coarsening algorithm */

  int                             max_levels;     /* advanced settings */
  cs_gnum_t                       min_n_g_rows;   /* advanced settings */
  double                          p0p1_relax;     /* advanced settings */

  int                             aggreg_limit;
  cs_param_amg_inhouse_coarsen_t  coarsen_algo;

  /* Down smoother */

  int                             n_down_iter;
  cs_param_amg_inhouse_solver_t   down_smoother;
  int                             down_poly_degree;

  /* Up smoother */

  int                             n_up_iter;
  cs_param_amg_inhouse_solver_t   up_smoother;
  int                             up_poly_degree;

  /* Coarse solver */

  double                          coarse_rtol_mult; /* advanced settings */
  int                             coarse_max_iter;  /* advanced settings */
  cs_param_amg_inhouse_solver_t   coarse_solver;
  int                             coarse_poly_degree;

} cs_param_amg_inhouse_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return true if the settings rely on boomeramg, otherwise false
 *
 * \param[in] solver   type of SLES solver
 * \param[in] precond  type of preconditioner
 * \param[in] amg      type of AMG
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_param_amg_boomer_is_needed(cs_param_solver_type_t   solver,
                              cs_param_precond_type_t  precond,
                              cs_param_amg_type_t      amg)
{
  if (precond == CS_PARAM_PRECOND_AMG || solver == CS_PARAM_SOLVER_AMG)
    if (amg == CS_PARAM_AMG_HYPRE_BOOMER_V ||
        amg == CS_PARAM_AMG_HYPRE_BOOMER_W)
      return true;

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return true if the settings rely on the in-house implementation,
 *        otherwise false
 *
 * \param[in] solver   type of SLES solver
 * \param[in] precond  type of preconditioner
 * \param[in] amg      type of AMG
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_param_amg_inhouse_is_needed(cs_param_solver_type_t   solver,
                               cs_param_precond_type_t  precond,
                               cs_param_amg_type_t      amg)
{
  if (precond == CS_PARAM_PRECOND_AMG || solver == CS_PARAM_SOLVER_AMG)
    if (amg == CS_PARAM_AMG_INHOUSE_K || amg == CS_PARAM_AMG_INHOUSE_V)
      return true;

  return false;
}

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
cs_param_amg_get_type_name(cs_param_amg_type_t  type);

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
cs_param_amg_get_class(cs_param_amg_type_t  amg_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new structure storing a set of parameters used when calling
 *        boomerAMG. Set default values for all parameters.
 *
 * \return a pointer to a new set of boomerAMG parameters
 */
/*----------------------------------------------------------------------------*/

cs_param_amg_boomer_t *
cs_param_amg_boomer_create(void);

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
cs_param_amg_boomer_copy(const cs_param_amg_boomer_t  *bamgp);

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
cs_param_amg_get_boomer_smoother_name(cs_param_amg_boomer_smoother_t  smoother);

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
                       const cs_param_amg_boomer_t  *bamgp);

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
                            bool  used_as_k_cycle);

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
cs_param_amg_inhouse_copy(const cs_param_amg_inhouse_t  *amgp);

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
cs_param_amg_get_inhouse_solver_name(cs_param_amg_inhouse_solver_t  solver);

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
                         const cs_param_amg_inhouse_t  *amgp);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_AMG_H__ */
