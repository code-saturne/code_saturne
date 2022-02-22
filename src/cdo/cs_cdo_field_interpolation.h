#ifndef __CS_CDO_FIELD_INTERPOLATION_H__
#define __CS_CDO_FIELD_INTERPOLATION_H__

/*============================================================================
 * Functions to handle field interpolation with CDO schemes
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

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @defgroup cdo_field_interpolation Interpolation flags for CDO
 *   Flags specifying which kind of interpolation will be used
 *   This enables to activate/add equations at the setup stage
 * @{
 */

/*!  1: Perform an interpolation at vertices of a scalar-valued field defined
 *      at cells. Activating this flag yields the creation of the equation and
 *      its related field called "scalar_c2v_field_interpolation"
 */

#define CS_CDO_FIELD_INTERPOLATION_SCALAR_C2V       (1 << 0)

/*!  2: Perform an interpolation at faces of a scalar-valued field defined
 *      at cells. Activating this flag yields the creation of the equation and
 *      its related field called "scalar_c2f_field_interpolation"
 */

#define CS_CDO_FIELD_INTERPOLATION_SCALAR_C2F       (1 << 1)

/*! @} */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate equation(s) used for the field interpolation
 *         Two choices are available which can be combined:
 *         CS_CDO_FIELD_INTERPOLATION_SCALAR_C2V
 *         CS_CDO_FIELD_INTERPOLATION_SCALAR_C2F
 *
 * \param[in]      mode            kind of interpolation to perform
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_field_interpolation_activate(cs_flag_t     mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate an array defined at vertices from an array defined at
 *         cells
 *
 * \param[in]      mesh            pointer to a mesh structure
 * \param[in]      cell_values     values at cells
 * \param[in, out] vtx_values      interpolated values at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_field_interpolation_cell_to_vertices(const cs_mesh_t    *mesh,
                                            const cs_real_t    *cell_values,
                                            cs_real_t          *vtx_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate an array defined at faces from an array defined at
 *         cells
 *
 * \param[in]      mesh            pointer to a mesh structure
 * \param[in]      cell_values     values at cells
 * \param[in, out] face_values     interpolated values at faces
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_field_interpolation_cell_to_faces(const cs_mesh_t    *mesh,
                                         const cs_real_t    *cell_values,
                                         cs_real_t          *face_values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_FIELD_INTERPOLATION_H__ */
