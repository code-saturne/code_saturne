#ifndef __CS_CF_BOUNDARY_FLUX_H__
#define __CS_CF_BOUNDARY_FLUX_H__

/*============================================================================
 * Compute the flux at the boundary.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Compute the flux at the boundary.
 *
 * \param[in]  f_id        face id
 * \param[in]  val_ext_en  Dirichlet value for the total energy
 * \param[in]  val_ext_p   Dirichlet value for the pressure
 * \param[in]  val_ext_v   Dirichlet value for the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_analytical_flux(const cs_lnum_t    f_id,
                               const cs_real_t   *val_ext_en,
                               const cs_real_t   *val_ext_p,
                               const cs_real_3_t *val_ext_v);

/*----------------------------------------------------------------------------*/
/*
 * \param[in]       f_id        face id
 * \param[in]       val_ext_en  Dirichlet value for the total energy
 * \param[in, out]  val_ext_p   Dirichlet value for the pressure
 * \param[in, out]  val_ext_v   Dirichlet value for the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_rusanov(const cs_lnum_t  f_id,
                       const cs_real_t *val_ext_en,
                       cs_real_t       *val_ext_p,
                       cs_real_3_t     *val_ext_v);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CF_BOUNDARY_FLUX_H__ */
