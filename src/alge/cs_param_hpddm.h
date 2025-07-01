#ifndef __CS_PARAM_HPDDM_H__
#define __CS_PARAM_HPDDM_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_param_hpddm.h

  \brief Routines and structure to handle the HPDDM setup. The structure is
         used as a context structure of a \ref cs_param_sles_t structure
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_param_hpddm_t
 *  \brief Set of parameters to specify additional options to HPDDM
 *  For more advanced settings, one has to use the \ref cs_user_sles_hpddm_hook
 *  function. Please also refer to the HPDDM user guide for more details.
 */

typedef struct {
  /* There is three important parameters */
  /* Left bounds: easy problem and low cost */
  /* Right bounds: hard problem and high cost */
  /* - 1 <= harmonic_overlap <= 10 */
  /* - 40 <= nb_eigenvector =< 500 */
  /* - 0.5 >= relative_threshold >= 1.e-6 */

  bool use_neumann; /*!< Use neumann matrix on each subdomains */

  int harmonic_overlap; /*!< Number of harmonic overlap */

  int nb_eigenvector; /*!< Number of eigenvector to compute
                       *    (= svd_nsv or eps_nev) */

  double relative_threshold; /*! < Thresold on eigenvalue to keep */

  /* Advanced options */
  bool adaptative; /* Compute parameters using adaptative algorithm */

  int min_harmonic_overlap,
    max_harmonic_overlap; /*!< Min/Max number of harmonic overlap */

  int min_nb_eigenvector,
    max_nb_eigenvector; /*!< Min/Max number of eigenvector to compute
                         *    (= svd_nsv or eps_nev) */

  int min_iter, max_iter; /*!< Adaptation of setup if nb_iter <= min_iter
                           *    or nb_iter >= max_iter */

  int nb_iter_prev; /* Number of iterations of previous solve */

} cs_param_hpddm_t;

/*============================================================================
 * Global variables
 *============================================================================*/

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
cs_param_hpddm_create(void);

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
cs_param_hpddm_copy(const cs_param_hpddm_t *hpddmp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the structure storing the set of parameters used with HPDDM
 *
 * \param[in] name     name related to the current SLES
 * \param[in] hpddmp   set of hpddm parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_param_hpddm_log(const char *name, const cs_param_hpddm_t *hpddmp);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_HPDDM_H__ */
