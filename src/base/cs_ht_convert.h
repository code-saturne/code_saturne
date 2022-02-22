#ifndef __CS_HT_CONVERT_H__
#define __CS_HT_CONVERT_H__

/*============================================================================
 * Enthaly to and from temperature conversion.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at all cells.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   h   enthalpy values
 * \param[out]  t   temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_h_to_t_cells(const cs_real_t  h[],
                           cs_real_t        t[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at solid cells only.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_h_to_t_cells_solid(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at all boundary faces.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   h   enthalpy values
 * \param[out]  t   temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_h_to_t_faces(const cs_real_t  h[],
                           cs_real_t        t[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at selected boundary faces.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   n_faces   number of selected boundary faces
 * \param[in]   face_ids  list of associated face ids
 * \param[in]   t         temperature values
 * \param[out]  h         enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_t_to_h_faces_l(cs_lnum_t        n_faces,
                             const cs_lnum_t  face_ids[],
                             const cs_real_t  t[],
                             cs_real_t        h[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy for a given boundary zone,
 *        using dense storage for temperature and enthalpy arrays.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   z  pointer to selected zone.
 * \param[in]   t  temperature values
 * \param[out]  h  enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_t_to_h_faces_z(const cs_zone_t *z,
                             const cs_real_t  t[],
                             cs_real_t        h[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HT_CONVERT_H__ */
