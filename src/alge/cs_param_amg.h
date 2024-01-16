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
  CS_PARAM_AMG_PETSC_GAMG_V,    /*!< V-cycle GAMG algorithm (PETSc lib.) */
  CS_PARAM_AMG_PETSC_GAMG_W,    /*!< W-cycle GAMG algorithm (PETSc lib.) */
  CS_PARAM_AMG_PETSC_PCMG,      /*!< preconditioned MG algorithm from PETSc */
  CS_PARAM_AMG_HOUSE_V,         /*!< In-house algorithm with V-cycle */
  CS_PARAM_AMG_HOUSE_K,         /*!< In-house algorithm with K-cycle */

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

/*! \enum cs_param_amg_boomer_interp_type_t
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

/*! \struct cs_param_boomer_amg_t

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

/*============================================================================
 * Global variables
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return true if the settings rely on boomeramg, otherwise false
 *
 * \param[in] solver     type of SLES solver
 * \param[in] precond    type of preconditioner
 * \param[in] amg        type of AMG
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_param_amg_boomer_is_needed(cs_param_itsol_type_t    solver,
                              cs_param_precond_type_t  precond,
                              cs_param_amg_type_t      amg)
{
  if (precond == CS_PARAM_PRECOND_AMG) {

    if (amg == CS_PARAM_AMG_HYPRE_BOOMER_V)
      return true;
    else if (amg == CS_PARAM_AMG_HYPRE_BOOMER_W)
      return true;

  }

  if (solver == CS_PARAM_ITSOL_AMG) {

    if (amg == CS_PARAM_AMG_HYPRE_BOOMER_V)
      return true;
    else if (amg == CS_PARAM_AMG_HYPRE_BOOMER_W)
      return true;

  }

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
 * \return the related solver class or CS_PARAM_SLES_CLASS_CS
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_class_t
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

END_C_DECLS

#endif /* __CS_PARAM_AMG_H__ */
