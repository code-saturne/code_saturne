#ifndef CS_FACE_TO_VERTEX_H
#define CS_FACE_TO_VERTEX_H

/*============================================================================
 * Face to vertex interpolation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Face to vertex computation method
 *----------------------------------------------------------------------------*/

typedef enum {
  CS_FACE_TO_VERTEX_UNWEIGHTED, /*!< Uniform (constant) weights */
  CS_FACE_TO_VERTEX_SURFACE,    /*!< weights by surface */
  CS_FACE_TO_VERTEX_SHEPARD,    /*!< Shepard interpolation
                                  (weights by inverse distance) */

} cs_face_to_vertex_type_t;

/*----------------------------------------------------------------------------
 * Main structure for interpolation
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
class cs_face_to_vertex_t {
private:
  cs_lnum_t        _n_faces;    /*!< Number of faces */
  const cs_lnum_t *_list_faces; /*!< List of faces */

  cs_lnum_t        _n_vtx;    /*!< Number of vertices */
  const cs_lnum_t *_list_vtx; /*!< List of verticies */

  cs_real_t *_v_w;   /*!< Weight for verticies */
  cs_real_t *_v_var; /*!< Values for verticies */

public:
  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Initialize internal strucutres.
   *
   * \param[in]       n_faces         number of faces
   * \param[in]       list_faces      list of face's ids or nullptr
   *                                  if all faces are used
   * \param[in]       n_vtx           number of verticies
   * \param[in]       list_vtx        list of vertex's ids or nullptr
   *                                  if all verticies are used
   */
  /*----------------------------------------------------------------------------*/

  void
  initialize(cs_lnum_t        n_faces,
             const cs_lnum_t *list_faces,
             cs_lnum_t        n_vtx,
             const cs_lnum_t *list_vtx);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief  Interpolate boundary faces values to vertex values.
   *         The size of the vector-values have to be compatible with the
   *         initialization.
   *
   * \param[in]       method            interpolation method
   * \param[in]       verbosity         verbosity level
   * \param[in]       ignore_rot_perio  if true, ignore periodicity of rotation
   * \param[in]       b_var             boundary-face values
   * \param[out]      v_var             vertex-based values
   */
  /*----------------------------------------------------------------------------*/

  void
  compute_on_boundary(cs_face_to_vertex_type_t method,
                      bool                     ignore_rot_perio,
                      const cs_real_t         *b_var,
                      cs_real_t               *v_var);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief  Free internal structures.
   */
  /*----------------------------------------------------------------------------*/

  void
  free();
};

/*============================================================================
 * Global variables
 *============================================================================*/

/* Short names for face to vertex methods */

extern const char *cs_face_to_vertex_type_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#endif /* CS_FACE_TO_VERTEX */
