#ifndef __CS_FIELD_OPERATOR_H__
#define __CS_FIELD_OPERATOR_H__

/*============================================================================
 * Field based algebraic operators.
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
#include "cs_field.h"
#include "cs_gradient.h"

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

typedef enum {

  CS_FIELD_INTERPOLATE_MEAN,      /* mean element value (P0 interpolation) */
  CS_FIELD_INTERPOLATE_GRADIENT   /* mean + gradient correction (pseudo-P1) */

} cs_field_interpolate_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_field_gradient_scalar(const cs_field_t          *f,
                         bool                       use_previous_t,
                         int                        inc,
                         bool                       recompute_cocg,
                         cs_real_3_t      *restrict grad);

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   f_ext          <-- exterior force generating the hydrostatic pressure
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_field_gradient_potential(const cs_field_t          *f,
                            bool                       use_previous_t,
                            int                        inc,
                            bool                       recompute_cocg,
                            int                        hyd_p_flag,
                            cs_real_3_t                f_ext[],
                            cs_real_3_t      *restrict grad);

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_field_gradient_vector(const cs_field_t          *f,
                         bool                       use_previous_t,
                         int                        inc,
                         cs_real_33_t     *restrict grad);

/*----------------------------------------------------------------------------
 * Compute cell gradient of tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_field_gradient_tensor(const cs_field_t          *f,
                         bool                       use_previous_t,
                         int                        inc,
                         cs_real_63_t     *restrict grad);

/*----------------------------------------------------------------------------
 * Interpolate field values at a given set of points.
 *
 * parameters:
 *   f                  <-- pointer to field
 *   interpolation_type <-- interpolation type
 *   n_points           <-- number of points at which interpolation
 *                          is required
 *   point_location     <-- location of points in mesh elements
 *                          (based on the field location)
 *   point_coords       <-- point coordinates
 *   val                --> interpolated values
 *----------------------------------------------------------------------------*/

void
cs_field_interpolate(cs_field_t              *f,
                     cs_field_interpolate_t   interpolation_type,
                     cs_lnum_t                n_points,
                     const cs_lnum_t          point_location[],
                     const cs_real_3_t        point_coords[],
                     cs_real_t               *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find local extrema of a given scalar field at each cell
 *
 * This assumes the field values have been synchronized.
 *
 * \param[in]       f id        scalar field id
 * \param[in]       halo_type   halo type
 * \param[in, out]  local_max   local maximum value
 * \param[in, out]  local_min   local minimum value
 */
/*----------------------------------------------------------------------------*/

void
cs_field_local_extrema_scalar(int              f_id,
                              cs_halo_type_t   halo_type,
                              cs_real_t       *local_max,
                              cs_real_t       *local_min);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize current parallel and periodic field values.
 *
 * This function currently only upates fields based on CS_MESH_LOCATION_CELLS.
 *
 * \param[in, out]   f           pointer to field
 * \param[in]        halo_type   halo type
 */
/*----------------------------------------------------------------------------*/

void
cs_field_synchronize(cs_field_t      *f,
                     cs_halo_type_t   halo_type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_OPERATOR_H__ */
