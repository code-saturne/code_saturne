#ifndef __CS_MESH_EXTRUDE_H__
#define __CS_MESH_EXTRUDE_H__

/*============================================================================
 * Mesh extrusion.
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

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extrude mesh boundary faces in the normal direction.
 *
 * Extrusion is defined on selected boundary faces, and the number of layers
 * for each associated vertex may be (slightly) variable, to account for
 * cluttered areas where extrusion may be constrained, or more complex
 * extrusions.
 *
 * \param[in, out]  m             mesh
 * \param[in]       interior_gc   if true, maintain group classes of
 *                                interior faces previously on boundary
 * \param[in]       n_faces       number of selected boundary faces
 * \param[in]       n_vertices    number of selected vertices
 * \param[in]       faces         list of selected boundary faces (0 to n-1),
 *                                or NULL if no indirection is needed
 * \param[in]       vertices      ids of selected vertices (0 to n-1),
 *                                or NULL if no indirection is needed
 * \param[in]       n_layers      number of layers for each vertex
 * \param[in]       coord_shift   extrusion vector for each vertex
 * \param[in]       distribution  optional distribution of resulting vertices
 *                                along each extrusion vector
 *                                (size: n_vertices*n_layers) with values
 *                                in range ]0, 1].
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude(cs_mesh_t          *m,
                bool                interior_gc,
                cs_lnum_t           n_faces,
                cs_lnum_t           n_vertices,
                const cs_lnum_t     faces[],
                const cs_lnum_t     vertices[],
                const cs_lnum_t     n_layers[],
                const cs_coord_3_t  coord_shift[],
                const float         distribution[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extrude mesh boundary faces in the normal direction by a constant
 *        thickness.
 *
 * \param[in, out]  m                 mesh
 * \param[in]       interior_gc       if true, maintain group classes of
 *                                    interior faces previously on boundary
 * \param[in]       n_layers          number of layers
 * \param[in]       thickness         extrusion thickness
 * \param[in]       expansion_factor  geometric expansion factor for
 *                                    extrusion refinement
 * \param[in]       n_faces           number of selected boundary faces
 * \param[in]       faces             list of selected boundary faces (0 to n-1),
 *                                    or NULL if no indirection is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_constant(cs_mesh_t        *m,
                         bool              interior_gc,
                         cs_lnum_t         n_layers,
                         double            thickness,
                         double            expansion_factor,
                         cs_lnum_t         n_faces,
                         const cs_lnum_t   faces[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_EXTRUDE_H__ */
