/*============================================================================
 * Routines/structure to handle the settings to solve a saddle-point problem
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

#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_saddle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_param_saddle.c

 * \brief Handle the settings of saddle-point systems.
 *        These systems arise from the monolithic coupling of the Navier-Stokes
 *        equations or in mixed formulation of scalar-valued equations.
 */

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
 * \brief Create a cs_param_saddle_t structure and assign default settings
 *
 * \param[in] block11_slesp   set of parameters for the (1,1) block
 *
 * \return a pointer to the new cs_param_saddle_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_saddle_t *
cs_param_saddle_create(const cs_param_sles_t   *block11_slesp)
{
  cs_param_saddle_t  *saddlep = NULL;

  BFT_MALLOC(saddlep, 1, cs_param_saddle_t);

  saddlep->verbosity = 0;
  saddlep->name = NULL;

  saddlep->solver = CS_PARAM_SADDLE_SOLVER_NONE;  /* Not used by default */
  saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;

  saddlep->cvg_param =  (cs_param_sles_cvg_t) {
    .n_max_iter = 100,
    .atol = 1e-12,       /* absolute tolerance */
    .rtol = 1e-6,        /* relative tolerance */
    .dtol = 1e3 };       /* divergence tolerance */

  saddlep->block11_sles_param = block11_slesp; /* shared */

  saddlep->schur_approximation = CS_PARAM_SCHUR_NONE;
  saddlep->schur_sles_param = NULL;

  return saddlep;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the structure storing the parameter settings for a saddle-point
 *        system
 *
 * \param[in, out] p_saddlep    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_free(cs_param_saddle_t  **p_saddlep)
{
  if (p_saddlep == NULL)
    return;

  cs_param_saddle_t  *saddlep = *p_saddlep;

  cs_param_sles_free(&(saddlep->schur_sles_param));

  /* block11_sles_param is shared */

  BFT_FREE(saddlep);
  *p_saddlep = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the name of the saddle-point system.
 *
 * \param[in]      basename   prefix for the naming of the Schur system
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_name(const char         *name,
                         cs_param_saddle_t  *saddlep)
{
  if (saddlep == NULL)
    return;

  size_t  len = strlen(name);
  BFT_MALLOC(saddlep->name, len + 1, char);
  strncpy(saddlep->name, name, len + 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a \ref cs_param_sles_t structure for the Schur
 *        approximation nested inside a \ref cs_param_saddle_t structure. By
 *        default, this member is not allocated. Do nothing if the related
 *        structure is already allocated.
 *
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_init_schur_sles(cs_param_saddle_t  *saddlep)
{
  if (saddlep == NULL)
    return;

  if (saddlep->schur_sles_param != NULL)
    return; /* Initialization has already been performed */

  char  *schur_name = NULL;
  const char  *basename = (saddlep->name == NULL) ?
    saddlep->block11_sles_param->name : saddlep->name;

  size_t  len = sizeof(basename) + sizeof("_schur_system");

  BFT_MALLOC(schur_name, len + 1, char);
  sprintf(schur_name, "%s_schur_system", basename);

  cs_param_sles_t  *schur_slesp = cs_param_sles_create(-1, schur_name);

  schur_slesp->precond = CS_PARAM_PRECOND_AMG;   /* preconditioner */
  schur_slesp->solver = CS_PARAM_ITSOL_FCG;      /* iterative solver */
  schur_slesp->amg_type = CS_PARAM_AMG_HOUSE_K;  /* no predefined AMG type */
  schur_slesp->cvg_param.rtol = 1e-3;            /* relative tolerance to stop
                                                    the iterative solver */

  BFT_FREE(schur_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a cs_param_saddle_t structure from ref to dest
 *
 * \param[in]      ref     reference structure to be copied
 * \param[in, out] dest    destination structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_copy(const cs_param_saddle_t  *ref,
                     cs_param_saddle_t        *dest)
{
  if (ref == NULL)
    return;

  dest->solver = ref->solver;
  dest->precond = ref->precond;
  dest->schur_approximation = ref->schur_approximation;

  dest->cvg_param.rtol = ref->cvg_param.rtol;
  dest->cvg_param.atol = ref->cvg_param.atol;
  dest->cvg_param.dtol = ref->cvg_param.dtol;
  dest->cvg_param.n_max_iter = ref->cvg_param.n_max_iter;

  dest->block11_sles_param = ref->block11_sles_param;

  if (ref->schur_sles_param != NULL) {

    if (dest->name == NULL)
      cs_param_saddle_set_name("automatic", dest); /* Avoid using the same
                                                      name */

    if (dest->schur_sles_param == NULL)
      cs_param_saddle_init_schur_sles(dest);

    cs_param_sles_copy_from(ref->schur_sles_param, dest->schur_sles_param);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup information for the given cs_param_saddle_t structure
 *
 * \param[in] saddlep     pointer to the structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_log(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == NULL)
    return;

  if (saddlep->solver == CS_PARAM_SADDLE_SOLVER_NONE)
    return;

  const char  *basename = (saddlep->name == NULL) ?
    saddlep->block11_sles_param->name : saddlep->name;

  char  *prefix = NULL;
  int  len = strlen(basename) + strlen("  *  |") + 1;
  BFT_MALLOC(prefix, len, char);
  sprintf(prefix, "  * %s |", basename);

  /* Log */

  cs_log_printf(CS_LOG_SETUP, "%s Verbosity: %d\n", prefix, saddlep->verbosity);

  /* Solver */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_GCR:
    cs_log_printf(CS_LOG_SETUP,
                  "%s Solver: Generalized Conjugate Residual (GCR)\n", prefix);
    break;

  case CS_PARAM_SADDLE_SOLVER_MINRES:
    cs_log_printf(CS_LOG_SETUP, "%s Solver: MINRES\n", prefix);
    break;

  case CS_PARAM_SADDLE_SOLVER_MUMPS:
    cs_log_printf(CS_LOG_SETUP, "%s Solver: MUMPS\n", prefix);
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s Solver: Undefined\n", prefix);
    break;

  } /* Solver */

  /* Preconditioner */

  switch (saddlep->precond) {

  case CS_PARAM_SADDLE_PRECOND_NONE:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: None\n", prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_DIAG_SCHUR:
    cs_log_printf(CS_LOG_SETUP,
                  "%s Precond: Diagonal blocks with Schur approx.\n",
                  prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_LOWER_SCHUR:
    cs_log_printf(CS_LOG_SETUP,
                  "%s Precond: Lower triangular blocks with Schur approx.\n",
                  prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_UPPER_SCHUR:
    cs_log_printf(CS_LOG_SETUP,
                  "%s Precond: Upper triangular blocks with Schur approx.\n",
                  prefix);
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Undefined\n", prefix);
    break;
  }

  /* Schur complement */
  /* ---------------- */

  /* Type of approximation */

  cs_param_get_schur_approx_name(saddlep->schur_approximation);

  /* SLES parameters */

  cs_param_sles_log(saddlep->schur_sles_param);

  BFT_FREE(prefix);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
