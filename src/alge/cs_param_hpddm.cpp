/*============================================================================
 * Routines and structure to handle the HPDDM settings
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_base.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_param_hpddm.h"

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
 * \brief Create and initialize with the default settings a new structure
 *        storing a set of parameters used when calling HPDDM.
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_param_hpddm_t *
cs_param_hpddm_create(void)
{
  cs_param_hpddm_t *hpddmp = nullptr;

  CS_MALLOC(hpddmp, 1, cs_param_hpddm_t);

  hpddmp->use_neumann = false;

  /* Basic options */

  hpddmp->harmonic_overlap   = 5;
  hpddmp->nb_eigenvector     = 100;
  hpddmp->relative_threshold = 5e-6;
  // Empirical choice advised by P. Jolivet
  hpddmp->p =
    cs::max(1, cs::min(cs_glob_n_ranks / 2 - 1, 1 + cs_glob_n_ranks / 8));

  /* Advanced options */

  hpddmp->adaptative = false;

  hpddmp->min_harmonic_overlap = hpddmp->harmonic_overlap;
  hpddmp->max_harmonic_overlap = hpddmp->harmonic_overlap;

  hpddmp->min_nb_eigenvector = hpddmp->nb_eigenvector;
  hpddmp->max_nb_eigenvector = hpddmp->nb_eigenvector;

  return hpddmp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy into a new structure the given set of parameters used when
 *        calling HPDDM
 *
 * \param[in] hpddmp   set of hpddm parameters
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_param_hpddm_t *
cs_param_hpddm_copy(const cs_param_hpddm_t *hpddmp)
{
  cs_param_hpddm_t *cpy = cs_param_hpddm_create();

  cpy->use_neumann = hpddmp->use_neumann;

  cpy->harmonic_overlap   = hpddmp->harmonic_overlap;
  cpy->nb_eigenvector     = hpddmp->nb_eigenvector;
  cpy->relative_threshold = hpddmp->relative_threshold;
  cpy->p                  = hpddmp->p;

  /* Advanced options */

  cpy->adaptative = hpddmp->adaptative;

  cpy->min_iter = hpddmp->min_iter;
  cpy->max_iter = hpddmp->max_iter;

  cpy->min_harmonic_overlap = hpddmp->min_harmonic_overlap;
  cpy->max_harmonic_overlap = hpddmp->max_harmonic_overlap;

  cpy->min_nb_eigenvector = hpddmp->min_nb_eigenvector;
  cpy->max_nb_eigenvector = hpddmp->max_nb_eigenvector;

  return cpy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the structure storing the set of parameters used with HPDDM
 *
 * \param[in] name     name related to the current SLES
 * \param[in] hpddmp   set of hpddm parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_param_hpddm_log(const char *name, const cs_param_hpddm_t *hpddmp)
{
  if (hpddmp == nullptr)
    return;

  if (hpddmp->use_neumann) {
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Use Neumann matrix for the coarse problem. \n",
                  name);
  }

  if (!hpddmp->use_neumann) {
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Harmonic overlap:   %d\n",
                  name,
                  hpddmp->harmonic_overlap);
  }

  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Number of eigenvectors:      %d\n",
                name,
                hpddmp->nb_eigenvector);

  if (!hpddmp->use_neumann)
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | SVD relative threshold:   %e\n",
                  name,
                  hpddmp->relative_threshold);

  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Number of mpi to solve coarse problem:      %d\n",
                name,
                hpddmp->p);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
