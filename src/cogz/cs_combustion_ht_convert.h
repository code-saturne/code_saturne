#ifndef CS_COMBUSTION_HT_CONVERT_H
#define CS_COMBUSTION_HT_CONVERT_H

/*============================================================================
 * Gas combustion model: enthaly to and from temperature conversion.
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

#include "base/cs_defs.h"

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
 * \brief Convert an enthalpy to temperature value for gas combustion.
 *
 * \param[in]     x_sp    mass fraction of constituents
 * \param[in]     h       enthalpy
 *
 * \return  temperature
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_combustion_h_to_t(const cs_real_t   x_sp[],
                     cs_real_t         h);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert an enthalpy to temperature value for gas combustion.
 *
 * \deprecated Use cs_combustion_h_to_t instead.
 *
 * \param[in]     x_sp    mass fraction of constituents
 * \param[in]     h       enthalpy
 *
 * \return  temperature
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_gas_combustion_h_to_t(const cs_real_t   x_sp[],
                         cs_real_t         h)
{
  return cs_combustion_h_to_t(x_sp, h);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert a temperature to enthalpy value for gas combustion.
 *
 * \param[in]     x_sp    mass fraction of constituents
 * \param[in]     t       temperature at cells
 *
 * \return  enthalpy
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_combustion_t_to_h(const cs_real_t   x_sp[],
                     cs_real_t         t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert a temperature to enthalpy value for gas combustion.
 *
 * \deprecated Use cs_combustion_t_to_h instead.
 *
 * \param[in]     x_sp    mass fraction of constituents
 * \param[in]     t       temperature at cells
 *
 * \return  enthalpy
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_gas_combustion_t_to_h(const cs_real_t   x_sp[],
                         cs_real_t         t)
{
  return cs_combustion_t_to_h(x_sp, t);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Convert enthalpy to temperature at boundary for gas combustion.
 *
 * \param[in]   h    enthalpy at boundary
 * \param[out]  t    temperature at boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_ht_convert_h_to_t_faces(const cs_real_t  h[],
                                      cs_real_t        t[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy for gas combustion.
 *
 * \param[in]   xsp  masss fraction of constituents
 * \param[in]   t    temperature value
 *
 * \return   enthalpy value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_combustion_ht_convert_t_to_h(cs_real_t  xsp[],
                                cs_real_t  t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at selected boundary faces
 *        for gas combustion.
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
cs_combustion_ht_convert_t_to_h_faces_l(cs_lnum_t        n_faces,
                                        const cs_lnum_t  face_ids[],
                                        const cs_real_t  t[],
                                        cs_real_t        h[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy for a given boundary zone with
 *        a gas combustio model, using dense storage for temperature
 *        and enthalpy arrays.
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
cs_combustion_ht_convert_t_to_h_faces_z(const cs_zone_t *z,
                                        const cs_real_t  t[],
                                        cs_real_t        h[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COMBUSTION_HT_CONVERT_H */
