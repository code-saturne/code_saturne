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

/*! Per face mesh extrusion settings;

 * This structure usually created or updated with utility functions, and may
 * be modified by the user in case fine control is needed */

typedef struct {

  cs_lnum_t   *n_layers;          /*!< number of layers for each boundary face;
                                   *   (< 0 for non-extruded faces) */

  cs_real_t   *distance;          /*!< total distance for each boundary face
                                    (if < 0, absolute value used as
                                    multiplier for boundary cell thickness) */

  float       *expansion_factor;  /*!< expansion factor for each boundary face */

  cs_real_t   *thickness_s;       /*!< optional start thickness for each boundary
                                   *   face; ignored if <= 0 */
  cs_real_t   *thickness_e;       /*!< optional end thickness for each boundary
                                   *   face; ignored if <= 0 */

} cs_mesh_extrude_face_info_t;

/*! Mesh extrusion vectors definition;

 * This structure defines local extrusion vectors; it is usually created
 * or updated with utility functions, and may be modified by the user
 * in case fine control is needed */

typedef struct {

  cs_lnum_t       n_faces;           /*!< number of associated faces */
  cs_lnum_t       n_vertices;        /*!< number of associated vertices */
  cs_lnum_t      *face_ids;          /*!< ids of associated faces, or NULL */
  cs_lnum_t      *vertex_ids;        /*!< ids of associated vertices, or NULL */
  cs_lnum_t      *n_layers;          /*!< number of layers for each vertex */
  cs_coord_3_t   *coord_shift;       /*!< extrusion vector for each vertex */
  cs_lnum_t      *distribution_idx;  /*!< index of optional distribution */
  float          *distribution;      /*!< optional distribution of resulting
                                      *   vertices along each extrusion vector,
                                      *   with values in range ]0, 1], or NULL
                                      *   (size: distribution_idx[n_vertices]) */

} cs_mesh_extrude_vectors_t;

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
 * \param[in]       e             extrusion vector definitions
 * \param[in]       interior_gc   if true, maintain group classes of
 *                                interior faces previously on boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude(cs_mesh_t                        *m,
                const cs_mesh_extrude_vectors_t  *e,
                bool                              interior_gc);

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
/*!
 * \brief Create a mesh extrusion face information structure.
 *
 * \param[in]  m  mesh
 *
 * \return pointer to new mesh extrusion face information structure.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_extrude_face_info_t *
cs_mesh_extrude_face_info_create(const cs_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a mesh extrusion face information structure.
 *
 * \param[in, out]  e  pointer to pointer to mesh extrusion face information.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_face_info_destroy(cs_mesh_extrude_face_info_t **efi);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set face extrusion information by zone.
 *
 * \param[in, out]  efi               mesh extrusion face information
 * \param[in]       n_layers          number of layers for selected faces
 * \param[in]       distance          extrusion distance for selected faces
 *                                    (if < 0, absolute value used as
 *                                    multiplier for boundary cell thickness)
 * \param[in]       expansion_factor  expansion factor for selected faces
 * \param[in]       n_faces           number of selected faces
 * \param[in]       face_ids          ids of selected faces, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_set_info_by_zone(cs_mesh_extrude_face_info_t  *efi,
                                 int                           n_layers,
                                 double                        distance,
                                 float                         expansion_factor,
                                 const cs_lnum_t               n_faces,
                                 const cs_lnum_t               face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and build a mesh extrusion vectors definition.
 *
 * Extrusion vectors will be computed based on the provided extrusion
 * face information structure. If no such structure is provided, an empty
 * structure is returned.
 *
 * \param[in]  efi  mesh extrusion face information, or NULL
 *
 * \return pointer to created mesh extrusion vectors definition.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_extrude_vectors_t *
cs_mesh_extrude_vectors_create(const cs_mesh_extrude_face_info_t  *efi);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a mesh extrusion vectors definition.
 *
 *
 * \param[in, out]  e  pointer to pointer to mesh extrusion vectors definition.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_vectors_destroy(cs_mesh_extrude_vectors_t  **e);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_EXTRUDE_H__ */
