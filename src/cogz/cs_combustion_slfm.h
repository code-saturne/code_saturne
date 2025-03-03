#ifndef CS_COMBUSTION_SLFM_H
#define CS_COMBUSTION_SLFM_H

/*============================================================================
 * Steady laminar flamelet gas combustion model.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_field.h"
#include "pprt/cs_physical_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for slfm gas combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_fields_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute physical properties for steady laminar flamelet model.
 *
 * \param[in]     iterns     Navier-Stokes sub-iterations indicator:
 *                           - if strictly negative, indicate that this
 *                                function is called outside Navier-Stokes loop
 *                           - if positive, Navier-Stokes iteration number.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_physical_properties(int   iterns);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute scalar dissipation rate for steady laminar flamelet model.
 *
 * \param[in]     f          pointer to the scalar field used to compute the
 *                           dissipation rate
 * \param[in]     cpro_rho   pointer to the density
 * \param[in]     fp2m       the variance associated to the scalar field
 * \param[in,out] cpro_totki pointer to the scalar dissipation rate
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_scalar_dissipation_rate(const cs_field_t   *f,
                                      const cs_real_t    *cpro_rho,
                                      const cs_real_t    *fp2m,
                                      cs_real_t          *cpro_totki);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct scalar variance in case of transporting 2nd order moment.
 *
 * \param[in]     fm         pointer to the mixture fraction
 * \param[in]     fsqm       pointer to the 2nd order moment of mixture fraction
 * \param[in,out] recvr      pointer to the reconstructed variance
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_reconstruct_variance(const cs_real_t   *fm,
                                   const cs_real_t   *fsqm,
                                   cs_real_t         *recvr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve maximal, middle and minimal values of progress variable
 * respectively on stable/unstable/mixing branches with a given mixture fraction
 *   - Neither heat loss nor variance is considered
 *
 * \param[in]     zm         mixture fraction value
 * \param[in]     cmax       maximal value of progress variable at given zm
 * \param[in]     cmid       middle value of progress variable at given zm
 * \param[in]     cmin       minimal value of progress variable at given zm
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_max_mid_min_progvar(const cs_real_t   zm,
                                       cs_real_t        *cmax,
                                       cs_real_t        *cmid,
                                       cs_real_t        *cmin);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Defines the source terms for the soot mass fraction and the precursor
 *        number for soot model of Moss et al for one time step.
 *
 *  The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 *  \f$ rovsdt \f$ et \f$ smbrs \f$ could already contain source term
 *  and don't have to be erased but incremented.
 *
 *  For stability sake, only positive terms should be add in \f$ rovsdt \f$.
 *  There is no constrain for \f$ smbrs \f$.
 *
 *  For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
 *           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
 *           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
 *
 *  Here are set \f$ rovsdt \f$ and \f$ smbrs \f$ containning \f$ \rho \Omega \f$
 *   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
 *     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.C.s^{-1} \f$,
 *     for enthalpy: \f$ J.s^{-1} \f$)
 *   - \f$ rovsdt \f$ en \f$ kg.s^{-1} \f$
 *
 * \param[in]      f_sc          pointer to scalar field
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_source_terms(cs_field_t  *f_sc,
                                cs_real_t    smbrs[],
                                cs_real_t    rovsdt[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COMBUSTION_SLFM_H */
