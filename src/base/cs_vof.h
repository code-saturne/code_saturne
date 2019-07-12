#ifndef __CS_VOF_H__
#define __CS_VOF_H__

/*============================================================================
 * Functions associated to VOF model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_domain.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
   @addtogroup vof

   \defgroup vof_masks Masks used to specify Volume of Fluid models

   @addtogroup vof_masks
   @{

 */

/*! Volume of Fluid model */
#define CS_VOF_ENABLED (1 << 0)

/*! Free surface model */
#define CS_VOF_FREE_SURFACE (1 << 1)

/*! Mass transfer Merkle model for vaporization / condensation (cavitation) */
#define CS_VOF_MERKLE_MASS_TRANSFER (1 << 2)

/*!
    @}

    @}
*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* VOF model parameters */
/*----------------------*/

typedef struct {

  unsigned      vof_model;   /* VoF model (sum of masks defining Volume of
                                Fluid model and submodels */

  double        rho1;        /* density */

  double        rho2;

  double        mu1;         /* viscosity */

  double        mu2;

} cs_vof_parameters_t;

/* Cavitation parameters */
/*-----------------------*/

typedef struct {

  cs_real_t        presat;  /* reference saturation pressure */
  cs_real_t        uinf;    /* reference velocity */
  cs_real_t        linf;    /* reference length scale */
  cs_real_t        cdest;   /* constant of condensation model (Merkle) */
  cs_real_t        cprod;   /* constant of vaporization model (Merkle) */
  int              icvevm;  /* eddy-viscosity correction indicator */
  cs_real_t        mcav;    /* eddy-viscosity correction cstt (Reboud) */
  int              itscvi;  /* eddy-viscosity correction indicator */

} cs_cavitation_parameters_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/* pointer to VOF model parameters structure */

extern const cs_vof_parameters_t *cs_glob_vof_parameters;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to VOF structure.
 */
/*----------------------------------------------------------------------------*/

cs_vof_parameters_t *
cs_get_glob_vof_parameters(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mixture density, mixture dynamic viscosity given fluid
 *         volume fractions and the reference density and dynamic viscosity
 *         \f$ \rho_l, \mu_l \f$ (liquid), \f$ \rho_v, \mu_v \f$ (gas).
 *
 * Computation is done as follows on cells:
 * \f[
 * \rho_\celli = \alpha_\celli \rho_v + (1-\alpha_\celli) \rho_l,
 * \f]
 * \f[
 * \mu_\celli = \alpha_\celli \mu_v + (1-\alpha_\celli) \mu_l,
 * \f]
 *
 * A similar linear formula is followed on boundary using fluid volume fraction
 * value on the boundary.
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_compute_linear_rho_mu(const cs_domain_t *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the mixture density, mixture dynamic viscosity and mixture
 *        mass flux given the volumetric flux, the volume fraction and the
 *        reference density and dynamic viscosity \f$ \rho_l, \mu_l \f$
 *        (liquid), \f$ \rho_v, \mu_v \f$ (gas).
 *
 * For the computation of mixture density, mixture dynamic viscosity, see
 * \ref cs_vof_compute_linear_rho_mu.
 *
 * Computation of mass flux is as follows:
 * \f[
 * \left( \rho\vect{u}\cdot\vect{S} \right)_\ij = \\ \left\lbrace
 * \begin{array}{ll}
 *   \rho_\celli (\vect{u}\cdot\vect{S})_\ij
 *  &\text{ if } (\vect{u}\cdot\vect{S})_\ij>0, \\
 *   \rho_\cellj (\vect{u}\cdot\vect{S})_\ij
 *  &\text{ otherwise },
 * \end{array} \right.
 * \f]
 * \f[
 * \left( \rho\vect{u}\cdot\vect{S} \right)_\ib = \\ \left\lbrace
 * \begin{array}{ll}
 *   \rho_\celli (\vect{u}\cdot\vect{S})_\ib
 *  &\text{ if } (\vect{u}\cdot\vect{S})_\ib>0, \\
 *   \rho_b (\vect{u}\cdot\vect{S})_\ib
 *  &\text{ otherwise }.
 * \end{array} \right.
 * \f]
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_update_phys_prop(const cs_domain_t *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write in main log the global mixture mass budget:
 * \f[
 * \sum_i\left(
 * |\Omega_i|\dfrac{\alpha_i^n - \alpha_i^{n-1}}{\Delta t} +
 * \sum_{j\in\Face{\celli}}\left(\rho\vect{u}\vect{S}\right)_{ij}^n
 * \right).
 * \f]
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_log_mass_budget(const cs_domain_t *domain);

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cavitation parameters structure.
 */
/*----------------------------------------------------------------------------*/

cs_cavitation_parameters_t *
cs_get_glob_cavitation_parameters(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VOF_H__ */
