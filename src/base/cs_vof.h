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

/*============================================================================
 * Type definitions
 *============================================================================*/

/* VOF mixture properties */
/*------------------------*/

typedef struct {

  double        rho1;        /* density */
  double        rho2;
  double        mu1;         /* viscosity */
  double        mu2;

} cs_vof_parameters_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/* pointer to VOF model indicator */

extern int cs_glob_vof_model;

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
 * \brief  Compute the mixture density, mixture dynamic viscosity and mixture
 *        mass flux given the volumetric flux, the volume fraction and the
 *        reference density and dynamic viscosity \f$ \rho_l, \mu_l \f$
 *        (liquid), \f$ \rho_v, \mu_v \f$ (gas) as follows:
 *
 * \f[
 * \rho_\celli = \alpha_\celli \rho_v + (1-\alpha_\celli) \rho_l,
 * \f]
 * \f[
 * \mu_\celli = \alpha_\celli \mu_v + (1-\alpha_\celli) \mu_l,
 * \f]
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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VOF_H__ */
