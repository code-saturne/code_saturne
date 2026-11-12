#ifndef CS_VOLUME_MASS_INJECTION_H
#define CS_VOLUME_MASS_INJECTION_H

/*============================================================================
 * Mass source terms computation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flag volume zones with the appropriate
 *        CS_VOLUME_ZONE_MASS_SOURCE_TERM flag when at least one volume
 *        mass injection on that zone is present.
 *
 * This is necessary for the reverse zone indexing required by the legacy code
 * to function with defintions that are partially unrolled an not purely
 * zone-based.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_flag_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the list and zone ids of cells with volume mass injection.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_build_lists(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate contributions to volume mass injection.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_eval(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointers to the mass source term arrays.
 *
 * \param[in]   f         pointer to associated field
 * \param[out]  n_elts    number of cells with mass source terms
 * \param[out]  elt_ids   pointer to source mass cells list (1-based numbering)
 * \param[out]  mst_type  mass source types (0: ambient value, 1: s_val value)
 * \param[out]  mst_val   pointer to mass source values
 * \param[out]  gamma     pointer to flow mass value
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_get_arrays(const cs_field_t   *f,
                                    cs_lnum_t          *n_elts,
                                    const cs_lnum_t   **elt_ids,
                                    int               **mst_type,
                                    cs_real_t         **mst_val,
                                    cs_real_t         **gamma);

/*----------------------------------------------------------------------------*/

#endif /* CS_VOLUME_MASS_INJECTION_H */
