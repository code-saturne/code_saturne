/*============================================================================
 * Filter symmetry faces whose effects cancel out
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_symmetry_faces_filter.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_symmetry_faces_filter.c
 * \brief Filter symmetry faces whose effects cancel out
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Filter selected faces whose effects should cancel out.
 *
 * This function simply checks if the sum of associated cell face normals
 * cancels out, and deselects faces for which this is not verified..
 *
 * \param[in]       m          pointer to mesh
 * \param[in]       mq         pointer to mesh quantities
 * \param[in, out]  n_faces    number of selected boundary faces
 * \param[in, out]  face_ids   ids of selected boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_symmetry_faces_filter_cancel(const cs_mesh_t             *m,
                                const cs_mesh_quantities_t  *mq,
                                cs_lnum_t                   *n_faces,
                                cs_lnum_t                    face_ids[])
{
  cs_lnum_t _n_faces = *n_faces;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;

  /* Flag selected faces */

  char *flag;
  BFT_MALLOC(flag, m->n_b_faces, char);

  for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
    flag[i] = 0;

  for (cs_lnum_t i = 0; i < _n_faces; i++)
    flag[face_ids[i]] = 1;

  char *c_flag;
  BFT_MALLOC(c_flag, m->n_cells, char);
  cs_real_4_t  *sum;
  BFT_MALLOC(sum, m->n_cells, cs_real_4_t);

  for (cs_lnum_t i = 0; i < m->n_cells; i++) {
    for (cs_lnum_t j = 0; j < 4; j++)
      sum[i][j] = 0;
  }

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    for (cs_lnum_t j = 0; j < 3; j++)
      sum[c_id][j] += b_face_normal[f_id][j];
    sum[c_id][3] += cs_math_3_norm(b_face_normal[f_id]);
  }

  /* Now filter flag */

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
    if (cs_math_3_norm(sum[c_id]) < 1e-2 * sum[c_id][3])
      c_flag[c_id] = 1;
    else
      c_flag[c_id] = 0;
  }

  BFT_FREE(sum);

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    if (c_flag[c_id] == 0)
      flag[f_id] = 0;
  }

  /* Rebuild list */

  _n_faces = 0;
  for (cs_lnum_t i = 0; i < m->n_b_faces; i++) {
    if (flag[i] != 0) {
      face_ids[_n_faces] = i;
      _n_faces++;
    }
  }

  *n_faces = _n_faces;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
