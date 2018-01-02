#ifndef __CS_INTERPOLATE_H__
#define __CS_INTERPOLATE_H__

/*============================================================================
 * Interpolation function handling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Field values interpolation type
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for interpolatation of values defined on
 *        a mesh location at a given set of points.
 *
 * \param[in, out]  input           pointer to optional (untyped) value
 *                                  or structure.
 * \param[in]       datatype        associated datatype
 * \param[in]       val_dim         dimension of data values
 * \param[in]       n_points        number of interpolation points
 * \param[in]       point_location  location of points in mesh elements
 * \param[in]       point_coords    point coordinates
 * \param[in]       location_vals   values at mesh location
 * \param[out]      point_vals      interpolated values at points
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_interpolate_from_location_t) (void                *input,
                                  cs_datatype_t        datatype,
                                  int                  val_dim,
                                  cs_lnum_t            n_points,
                                  const cs_lnum_t      point_location[],
                                  const cs_real_3_t    point_coords[],
                                  const void          *location_vals,
                                  void                *point_vals);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interpolate values defined on a mesh location at a given set of
 *        points using a P0 interpolation.
 *
 * This function allows unlocated points (with point_location < 0),
 * to which the value 0 is assigned.
 *
 * \param[in, out]  input           pointer to optional (untyped) value
 *                                  or structure.
 * \param[in]       datatype        associated datatype
 * \param[in]       val_dim         dimension of data values
 * \param[in]       n_points        number of interpolation points
 * \param[in]       point_location  location of points in mesh elements
 * \param[in]       point_coords    point coordinates
 * \param[in]       location_vals   values at mesh location
 * \param[out]      point_vals      interpolated values at points
 */
/*----------------------------------------------------------------------------*/

void
cs_interpolate_from_location_p0(void                *input,
                                cs_datatype_t        datatype,
                                int                  val_dim,
                                cs_lnum_t            n_points,
                                const cs_lnum_t      point_location[],
                                const cs_real_3_t    point_coords[],
                                const void          *location_vals,
                                void                *point_vals);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INTERPOLATE_H__ */
