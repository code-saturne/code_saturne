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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_interpolate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_interpolate.c
        Interpolation function handling.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
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
                                void                *point_vals)
{
  CS_UNUSED(input);
  CS_UNUSED(point_coords);

  switch(datatype) {

  case CS_INT32:
    {
      const int32_t *l_vals = (const int32_t *)location_vals;
      int32_t *p_vals = point_vals;
      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t e_id = point_location[i];
        if (e_id > -1) {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = l_vals[e_id*val_dim + j];
        }
        else {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = 0;
        }
      }
    }
    break;

  case CS_INT64:
    {
      const int64_t *l_vals = (const int64_t *)location_vals;
      int64_t *p_vals = point_vals;
      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t e_id = point_location[i];
        if (e_id > -1) {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = l_vals[e_id*val_dim + j];
        }
        else {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = 0;
        }
      }
    }
    break;

  case CS_DOUBLE:
    {
      const cs_real_t *l_vals = (const cs_real_t *)location_vals;
      cs_real_t *p_vals = point_vals;
      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t e_id = point_location[i];
        if (e_id > -1) {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = l_vals[e_id*val_dim + j];
        }
        else {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = 0;
        }
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Function %s does not currently handle %s data type."),
              __func__, cs_datatype_name[datatype]);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
