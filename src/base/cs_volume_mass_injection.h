#ifndef __CS_VOLUME_MASS_INJECTION_H__
#define __CS_VOLUME_MASS_INJECTION_H__

/*============================================================================
 * Mass source terms computation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 *
 * \param[in]   n_cells       number of cells in mass source term zones
 * \param[out]  cell_num      numbers (1-based) cells in mass source term zones
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_build_lists(cs_lnum_t   n_cells,
                                     cs_lnum_t   cell_num[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate contributions to volume mass injection.
 *
 * \param[in]     nvar          total number of variables
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     itypsm        mass source type for the working variable
 *                              size: [nvar][ncesmp]
 * \param[in]     smacel        values of the variables associated to the
 *                              mass source (for the pressure variable,
 *                              smacel is the mass flux)
 *                              size: [nvar][ncesmp]
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_eval(int        nvar,
                              cs_lnum_t  ncesmp,
                              int        itypsm[],
                              cs_real_t  smacel[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointers to the mass source term arrays.
 *
 * \param[in]   f         pointer to associated field
 * \param[out]  ncesmp    number of cells with mass source terms
 * \param[out]  icetsm    pointet to source mass cells list (1-based numbering)
 * \param[out]  itpsmp    mass source type for the working variable
 * \param[out]  s_type    mass source types (0: ambient value, 1: s_val value)
 * \param[out]  smcelp    pointer to mass source values
 * \param[out]  gamma     pointer to flow mass value
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_get_arrays(const cs_field_t  *f,
                                    cs_lnum_t         *ncesmp,
                                    cs_lnum_t         **icetsm,
                                    int               **itpsmp,
                                    cs_real_t         **smcelp,
                                    cs_real_t         **gamma);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VOLUME_MASS_INJECTION_H__ */
