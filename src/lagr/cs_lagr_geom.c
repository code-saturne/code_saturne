/*============================================================================
 * Methods for particle parameters
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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_lagr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_geom.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build geometric information needed by the deposition model.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_geom(void)
{
  cs_mesh_t            *mesh = cs_glob_mesh;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->b_face_normal;
  const cs_real_3_t *vtx_coord
    = (const cs_real_3_t *)(cs_glob_mesh->vtx_coord);

  BFT_REALLOC(cs_glob_lagr_b_face_proj, mesh->n_b_faces, cs_real_33_t);

/* ==============================================================================*/
/* 1. The reference frame chang matrix is:
 *
 *  ( nx  ny  nz )
 *  ( t1x t1y t1z)
 *  ( t2x t2y t2z)
 *
 *  With n, t1, t2 three unit orthogonal vectors and n the outwrding normal to the
 *  boundary face.
 * */
/* ==============================================================================*/

  for (cs_lnum_t face_id = 0; face_id < mesh->n_b_faces; face_id++) {

    /* normal vector coordinates */
    cs_real_3_t normal;
    cs_math_3_normalise(b_face_normal[face_id], normal);

    /* Recover the first face nodes */
    cs_lnum_t v_id0  = mesh->b_face_vtx_lst[mesh->b_face_vtx_idx[face_id]];
    cs_lnum_t v_id1  = mesh->b_face_vtx_lst[mesh->b_face_vtx_idx[face_id] + 1];

    cs_real_3_t v0v1 = {
      vtx_coord[v_id1][0] - vtx_coord[v_id0][0],
      vtx_coord[v_id1][1] - vtx_coord[v_id0][1],
      vtx_coord[v_id1][2] - vtx_coord[v_id0][2]};

    /* tangential projection to the wall:
     * (Id -n (x) n) vect */
    cs_real_3_t t1, t1p, t2;
    cs_math_3_orthogonal_projection(normal, v0v1, t1p);

    cs_math_3_normalise(t1p, t1);

    /* t2 = n ^ t1 */
    cs_math_3_cross_product(normal, t1, t2);

    /* Matrix of Reference frame change    */
    cs_glob_lagr_b_face_proj[face_id][0][0] = normal[0];
    cs_glob_lagr_b_face_proj[face_id][0][1] = normal[1];
    cs_glob_lagr_b_face_proj[face_id][0][2] = normal[2];
    cs_glob_lagr_b_face_proj[face_id][1][0] = t1[0];
    cs_glob_lagr_b_face_proj[face_id][1][1] = t1[1];
    cs_glob_lagr_b_face_proj[face_id][1][2] = t1[2];
    cs_glob_lagr_b_face_proj[face_id][2][0] = t2[0];
    cs_glob_lagr_b_face_proj[face_id][2][1] = t2[1];
    cs_glob_lagr_b_face_proj[face_id][2][2] = t2[2];

  }

}
