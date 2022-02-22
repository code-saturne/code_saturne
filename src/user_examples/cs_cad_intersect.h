#ifndef __CS_CAD_INTERSECT_H__
#define __CS_CAD_INTERSECT_H__

/*============================================================================
 * Intersect cells with CAD object.
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/* File modes */

typedef enum {

  CS_CAD_INTERSECT_COMMON,   /* boolean common operation */
  CS_CAD_INTERSECT_CUT       /* boolean cut operation */

} cs_cad_intersect_op_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Intersect selected cells with CAD shape.
 *
 * \param[in]       m                pointer to a mesh structure
 * \param[in]       path             path to CAD file
 * \param[in]       op               common if CAD represents fluid domain
 *                                   cut if CAD represents solid complement
 * \param[in]       n_cells          number of selected cells
 * \param[in]       cell_ids         ids of selected cells
 * \param[in, out]  cell_porosity    cell porosity
 * \param[in, out]  cell_f_center    cell fluid center, or NULL
 * \param[in, out]  i_face_porosity  interior face porosity, or NULL
 * \param[in, out]  i_face_f_center  interior face fluid center, or NULL
 * \param[in, out]  b_face_porosity  boundary face porosity, or NULL
 * \param[in, out]  b_face_f_center  boundary face fluid center, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_cad_intersect(const cs_mesh_t        *m,
                 const char             *path,
                 cs_cad_intersect_op_t   op,
                 cs_lnum_t               n_cells,
                 const cs_lnum_t         cell_ids[],
                 cs_real_t               cell_porosity[],
                 cs_real_t               cell_f_center[][3],
                 cs_real_t               i_face_porosity[],
                 cs_real_t               i_face_f_center[][3],
                 cs_real_t               b_face_porosity[],
                 cs_real_t               b_face_f_center[][3]);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CAD_INTERSECT_H__ */
