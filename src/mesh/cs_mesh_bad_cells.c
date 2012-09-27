/*============================================================================
 * \file Detect bad cells within meshes.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * BFT and FVM library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_bad_cells.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_MESH_QUALITY_N_SUBS  10

#undef _CROSS_PRODUCT_3D
#undef _DOT_PRODUCT_3D
#undef _MODULE_3D
#undef _COSINE_3D

#define _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2) ( \
 cross_v1_v2[0] = v1[1]*v2[2] - v1[2]*v2[1],   \
 cross_v1_v2[1] = v1[2]*v2[0] - v1[0]*v2[2],   \
 cross_v1_v2[2] = v1[0]*v2[1] - v1[1]*v2[0]  )

#define _DOT_PRODUCT_3D(v1, v2) ( \
 v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define _MODULE_3D(v) \
 sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

#define _COSINE_3D(v1, v2) (\
 _DOT_PRODUCT_3D(v1, v2) / (_MODULE_3D(v1) * _MODULE_3D(v2)) )

/*============================================================================
 * Private function definitions
 *============================================================================*/

static const int _type_flag_mask[] = {CS_BAD_CELL_ORTHO_NORM,
                                      CS_BAD_CELL_OFFSET,
                                      CS_BAD_CELL_LSQ_GRAD,
                                      CS_BAD_CELL_RATIO,
                                      CS_BAD_CELL_GUILT,
                                      CS_BAD_CELL_USER};

/*----------------------------------------------------------------------------
 * Evaluate cell's non-orthogonality.
 *
 * Compute orthogonal normal using the distance between two consecutive
 * cell centers and the surface vector orthogonal to the face.
 * Evaluates a level of non-orthogonality and tags identified bad cells.
 *
 * parameters:
 *   mesh                 <-- pointer to associated mesh structure.
 *   mesh_quantities      <-- pointer to associated mesh quantities structure
 *   bad_cell_flag        <-- array of bad cell flags for various uses
 *----------------------------------------------------------------------------*/

static void
_compute_ortho_norm(const cs_mesh_t             *mesh,
                    const cs_mesh_quantities_t  *mesh_quantities,
                    cs_real_t                    i_face_ortho[],
                    cs_real_t                    b_face_ortho[],
                    unsigned                     bad_cell_flag[])
{
  cs_lnum_t  i, face_id, cell1, cell2;

  cs_real_t  cell_center1[3], cell_center2[3];
  cs_real_t  face_center[3];
  cs_real_t  face_normal[3], vect[3];

  double  cos_alpha;

  const cs_lnum_t  dim = 3;

  /* Loop on interior faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[2 * face_id] - 1;
    cell2 = mesh->i_face_cells[2 * face_id + 1] - 1;

    /* Get information on mesh quantities */

    for (i = 0; i < dim; i++) {

      /* Center of gravity for each cell */
      cell_center1[i] = mesh_quantities->cell_cen[cell1*dim + i];
      cell_center2[i] = mesh_quantities->cell_cen[cell2*dim + i];

      /* Surface vector (orthogonal to the face) */
      face_normal[i] = mesh_quantities->i_face_normal[face_id*dim + i];

    }

    /* Evaluate the non-orthogonality. */

    for (i = 0; i < dim; i++)
      vect[i] = cell_center2[i] - cell_center1[i];

    cos_alpha = _COSINE_3D(vect, face_normal);

    i_face_ortho[face_id] = cos_alpha;

    if (i_face_ortho[face_id] < 0.1) {
      bad_cell_flag[cell1] = bad_cell_flag[cell1] | _type_flag_mask[0];
    }
  }

  /* Loop on boundary faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {

    cell1 = mesh->b_face_cells[face_id] - 1;

    /* Get information on mesh quantities */

    for (i = 0; i < dim; i++) {

      /* Center of gravity of the cell */
      cell_center1[i] = mesh_quantities->cell_cen[cell1*dim + i];

      /* Face center coordinates */
      face_center[i] = mesh_quantities->b_face_cog[face_id*dim + i];

      /* Surface vector (orthogonal to the face) */
      face_normal[i] = mesh_quantities->b_face_normal[face_id*dim + i];

    }

    /* Evaluate the non-orthogonality. */
    for (i = 0; i < dim; i++)
      vect[i] = face_center[i] - cell_center1[i];

    cos_alpha = _COSINE_3D(vect, face_normal);

    b_face_ortho[face_id] = cos_alpha;

    if (b_face_ortho[face_id] < 0.1) {
      bad_cell_flag[cell1] = bad_cell_flag[cell1] | _type_flag_mask[0];
    }
  }
}

/*----------------------------------------------------------------------------
 * Evaluate cell's center offsetting.
 *
 * Compute center offsetting coefficient for interior faces.
 * Evaluates a non-matching level (offset) and tags identified bad cells.
 *
 * parameters:
 *   mesh                   <-- pointer to associated mesh structure
 *   mesh_quantities        <-- pointer to associated mesh quantities
 *                              structure
 *   bad_cell_flag          <-- array of bad cell flags for various uses
 *----------------------------------------------------------------------------*/

static void
_compute_weighting_offsetting(const cs_mesh_t         *mesh,
                              const cs_mesh_quantities_t  *mesh_quantities,
                              cs_real_t                weighting[],
                              cs_real_t                offsetting[],
                              unsigned                 bad_cell_flag[])
{
  cs_lnum_t  i, face_id, cell1, cell2;
  cs_real_t  intersection;
  cs_real_t  cell_center1[3], cell_center2[3];
  cs_real_t  face_center[3], face_normal[3];
  cs_real_t  v0[3], v1[3], v2[3];

  double  coef0 = 0, coef1 = 0, coef2 = 0;

  const cs_lnum_t  dim = 3;

  /* Compute weighting coefficient */
  /*-------------------------------*/

  /* Loop on interior faces */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[2 * face_id] - 1;
    cell2 = mesh->i_face_cells[2 * face_id + 1] - 1;

    /* Get information on mesh quantities */

    for (i = 0; i < dim; i++) {

      /* Center of gravity for each cell */
      cell_center1[i] = mesh_quantities->cell_cen[cell1*dim + i];
      cell_center2[i] = mesh_quantities->cell_cen[cell2*dim + i];

      /* Face center coordinates */
      face_center[i] = mesh_quantities->i_face_cog[face_id*dim + i];

      /* Surface vector (orthogonal to the face) */
      face_normal[i] = mesh_quantities->i_face_normal[face_id*dim + i];
    }

    /* Compute weighting coefficient with two approaches. Keep the max value. */

    for (i = 0; i < dim; i++) {
      v0[i] = cell_center2[i] - cell_center1[i];
      v1[i] = face_center[i] - cell_center1[i];
      v2[i] = cell_center2[i] - face_center[i];
    }

    coef0 = _DOT_PRODUCT_3D(v0, face_normal);
    coef1 = _DOT_PRODUCT_3D(v1, face_normal)/coef0;
    coef2 = _DOT_PRODUCT_3D(v2, face_normal)/coef0;

    weighting[face_id] = CS_MAX(coef1, coef2);

    /* Compute center offsetting coefficient */
    /*---------------------------------------*/

    /* Compute intersection between face and segment defined by the two cell
       centers */

    for (i = 0; i < dim; i++) {
      intersection =  (1 - weighting[face_id]) * cell_center1[i]
                         + weighting[face_id]  * cell_center2[i];
      v1[i] = intersection - face_center[i];
      v2[i] = cell_center2[i] - cell_center1[i];
    }

    offsetting[face_id] = 1 - sqrt(  _DOT_PRODUCT_3D(v1, v1)
                                   / _MODULE_3D(face_normal));

    if (offsetting[face_id] < 0.1) {
      bad_cell_flag[cell1] = bad_cell_flag[cell1] | _type_flag_mask[1];
    }
  }
}

/*----------------------------------------------------------------------------
 * Evaluate cell's distorsion.
 *
 * Compute Least Squares Gradient coefficient for cells.
 * Evaluates a distorsion level (based on LSQ Gradient Method) and tags
 * identified bad cells.
 *
 * parameters:
 *   mesh               <-- pointer to associated mesh structure
 *   mesh_quantities    <-- pointer to associated mesh quantities structure
 *   bad_cell_flag      --> array of bad cell flags for various uses
 *----------------------------------------------------------------------------*/

static void
_compute_least_squares(const cs_mesh_t             *mesh,
                       const cs_mesh_quantities_t  *mesh_quantities,
                       cs_real_t                    lsq[],
                       unsigned                     bad_cell_flag[])
{
  const cs_lnum_t  dim = mesh->dim;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;

  const cs_real_t *surfbo  = mesh_quantities->b_face_normal;

  int           sweep_id, n_sweeps;
  cs_lnum_t     i, k, face_id, cell1, cell2, cell_id;
  cs_real_3_t   cell_center1, cell_center2, vect, dij, eigenvalues;
  cs_real_33_t  w2, rot;

  double unsdij, surfn, surf_n_inv, min_diag, max_diag;
  double pivot, y, t, s, c;

  cs_lnum_t  indexl = 0, indexc = 0;
  cs_real_t *w1 = NULL;
  BFT_MALLOC(w1, 6 * n_cells_wghosts, cs_real_t);

  for (i = 0; i < 6 * n_cells_wghosts; i++)
    w1[i] = 0.;

  /* Loop on interior faces */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[2 * face_id] - 1;
    cell2 = mesh->i_face_cells[2 * face_id + 1] - 1;

    /* Center of gravity for each cell */

    for (i = 0; i < dim; i++) {
      cell_center1[i] = mesh_quantities->cell_cen[cell1*dim + i];
      cell_center2[i] = mesh_quantities->cell_cen[cell2*dim + i];
      vect[i] = cell_center2[i] - cell_center1[i];
    }

    unsdij = 1.0 / sqrt(pow(vect[0], 2) + pow(vect[1], 2) + pow(vect[2], 2));

    for (i = 0; i < dim; i++)
      dij[i] = vect[i] * unsdij;

    w1[cell1] += dij[0] * dij[0];
    w1[cell1 + n_cells_wghosts] += dij[1] * dij[1];
    w1[cell1 + 2 * n_cells_wghosts] += dij[2] * dij[2];
    w1[cell1 + 3 * n_cells_wghosts] += dij[0] * dij[1];
    w1[cell1 + 4 * n_cells_wghosts] += dij[0] * dij[2];
    w1[cell1 + 5 * n_cells_wghosts] += dij[1] * dij[2];

    w1[cell2] += dij[0] * dij[0];
    w1[cell2 + n_cells_wghosts] += dij[1] * dij[1];
    w1[cell2 + 2 * n_cells_wghosts] += dij[2] * dij[2];
    w1[cell2 + 3 * n_cells_wghosts] += dij[0] * dij[1];
    w1[cell2 + 4 * n_cells_wghosts] += dij[0] * dij[2];
    w1[cell2 + 5 * n_cells_wghosts] += dij[1] * dij[2];

  }

  /* Loop on boundary faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
    cell1 = mesh->b_face_cells[face_id] - 1;

    surfn = sqrt(  pow(surfbo[face_id * 3], 2)
                 + pow(surfbo[face_id * 3 + 1], 2)
                 + pow(surfbo[face_id * 3 + 2], 2));

    surf_n_inv = 1.0 / surfn;

    for (i = 0; i < dim; i++)
      dij[i] = surfbo[face_id * 3 + i] * surf_n_inv;

    w1[cell1] += dij[0] * dij[0];
    w1[cell1 + n_cells_wghosts] += dij[1] * dij[1];
    w1[cell1 + 2 * n_cells_wghosts] += dij[2] * dij[2];
    w1[cell1 + 3 * n_cells_wghosts] += dij[0] * dij[1];
    w1[cell1 + 4 * n_cells_wghosts] += dij[0] * dij[2];
    w1[cell1 + 5 * n_cells_wghosts] += dij[1] * dij[2];
  }

  /* Use of a Jacobi transformation to estimate the min/max eigenvalues ratio */

  n_sweeps = 100;

  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    w2[0][0] = w1[cell_id];
    w2[1][1] = w1[cell_id + n_cells_wghosts];
    w2[2][2] = w1[cell_id + 2 * n_cells_wghosts];
    w2[0][1] = w1[cell_id + 3 * n_cells_wghosts];
    w2[0][2] = w1[cell_id + 4 * n_cells_wghosts];
    w2[1][2] = w1[cell_id + 5 * n_cells_wghosts];
    w2[1][0] = w1[cell_id + 3 * n_cells_wghosts];
    w2[2][0] = w1[cell_id + 4 * n_cells_wghosts];
    w2[2][1] = w1[cell_id + 5 * n_cells_wghosts];

    for (sweep_id = 0; sweep_id < n_sweeps; sweep_id++) {

      pivot  = 0.;

      for (i = 0; i < 3; i++)
        eigenvalues[i] = w2[i][i];

      for (i = 0; i < 3; i++) {
        for (k = 0; k < 3; k++) {
          if (fabs(w2[i][k]) >= pivot && i != k) {
            pivot = fabs(w2[i][k]);
            indexl = i;
            indexc = k;
          }
        }
      }

      if (pivot < 1.e-10)
        break;

      y = (w2[indexc][indexc] - w2[indexl][indexl]) / (2.*w2[indexl][indexc]);

      if (y >= 0.)
        t = 1. / (fabs(y) + sqrt(y*y + 1));
      else
        t = -1. / (fabs(y) + sqrt(y*y + 1));

      c = 1. / sqrt(t*t + 1);
      s = t*c;

      for (i = 0; i < 3; i++) {
        for (k = 0; k < 3; k++)
          rot[i][k] = 0.;
      }

      for (i = 0; i < 3; i++)
        rot[i][i] = 1.;

      rot[indexl][indexl] = c;
      rot[indexc][indexc] = c;
      rot[indexl][indexc] = -s;
      rot[indexc][indexl] = s;

      for (i = 0; i < 3; i++)
        for (k = 0; k < 3; k++)
          w2[i][k] = rot[i][0]*w2[0][k] +
                     rot[i][1]*w2[1][k] +
                     rot[i][2]*w2[2][k];

      rot[indexl][indexc] = s;
      rot[indexc][indexl] = -s;

      for (i = 0; i < 3; i++)
        for (k = 0; k < 3; k++)
          w2[i][k] = w2[i][0]*rot[0][k] +
                     w2[i][1]*rot[1][k] +
                     w2[i][2]*rot[2][k];
    }

    min_diag = 1.e15;
    max_diag = 0.0;

    for (i = 0; i < 3; i++) {
      min_diag = fmin(min_diag, fabs(eigenvalues[i]));
      max_diag = fmax(max_diag, fabs(eigenvalues[i]));
    }

    lsq[cell_id] = min_diag / max_diag;

    if (lsq[cell_id] < 0.1) {
      bad_cell_flag[cell_id] = _type_flag_mask[2];
    }
  }

  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------
 * Evaluate cell's volume ratio.
 *
 * Compute Volume Ratio coefficient for cells.
 * Evaluates a the cell's geometric continuity (based on volume ratio) and
 * tags identified bad cells.
 *
 * parameters:
 *   mesh                <-- pointer to associated mesh structure
 *   mesh_quantities     <-- pointer to associated mesh quantities structure
 *   bad_cell_flag       --> array of bad cell flags for various uses
 *----------------------------------------------------------------------------*/

static void
_compute_volume_ratio(const cs_mesh_t             *mesh,
                      const cs_mesh_quantities_t  *mesh_quantities,
                      cs_real_t                    vol_ratio[],
                      unsigned                     bad_cell_flag[])
{
  cs_real_t *volume = mesh_quantities->cell_vol;

  cs_lnum_t   face_id, cell1, cell2;

  /* Loop on interior faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[2 * face_id] - 1;
    cell2 = mesh->i_face_cells[2 * face_id + 1] - 1;

    vol_ratio[face_id] = fmin(volume[cell1] / volume[cell2],
                              volume[cell2] / volume[cell1 ]);

    if (vol_ratio[face_id] < 0.1*0.1) {
      bad_cell_flag[cell1] = _type_flag_mask[3];
      bad_cell_flag[cell2] = _type_flag_mask[3];
    }
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute bad cell quality indicators.
 *
 * \param [in]       mesh             pointer to associated mesh structure
 * \param [in, out]  mesh_quantities  pointer to associated mesh quantities
 *                                    structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_detect(const cs_mesh_t       *mesh,
                         cs_mesh_quantities_t  *mesh_quantities)
{
  double  *working_array = NULL;
  unsigned *bad_cell_flag = NULL;

  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;

  cs_lnum_t i;
  cs_gnum_t n_cells_tot, iwarning, ibad;

  /* Check input data */

  assert(mesh_quantities->i_face_normal != NULL);
  assert(mesh_quantities->i_face_cog != NULL);
  assert(mesh_quantities->cell_cen != NULL);
  assert(mesh_quantities->cell_vol != NULL);

  /* Global bad cells storing array initialization */

  if (mesh_quantities->bad_cell_flag == NULL)
    BFT_MALLOC(mesh_quantities->bad_cell_flag,
               mesh->n_cells_with_ghosts,
               unsigned);

  bad_cell_flag = mesh_quantities->bad_cell_flag;

  for (i = 0; i < n_cells_wghosts; i++)
    bad_cell_flag[i] = 0;

  /* Possible warning printed in the listing --> flag initialization */

  iwarning = 0;
  n_cells_tot = mesh->n_g_cells;


  /* Evaluate mesh quality criteria */
  /*--------------------------------*/

  /* Condition 1: Orthogonal Normal */
  /*--------------------------------*/

  cs_real_t  *i_face_ortho = NULL, *b_face_ortho = NULL;

  BFT_MALLOC(working_array, n_i_faces + n_b_faces, cs_real_t);

  for (i = 0; i < n_i_faces + n_b_faces; i++)
    working_array[i] = 0.;

  i_face_ortho = working_array;
  b_face_ortho = working_array + n_i_faces;

  _compute_ortho_norm(mesh,
                      mesh_quantities,
                      i_face_ortho,
                      b_face_ortho,
                      bad_cell_flag);

  ibad = 0;
  for (i = 0; i < n_cells; i++) {
    if (bad_cell_flag[i] & _type_flag_mask[0]) {
      ibad++;
      iwarning++;
    }
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_counter(&ibad, 1);
    cs_parall_counter(&iwarning, 1);
  }

  /* Display log output */
  bft_printf(_("\n  Criteria 1: Orthogonality:\n"));
  bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
             (unsigned long long)ibad,
             (double)ibad / (double)n_cells_tot * 100.0);

  BFT_FREE(working_array);

  /* Condition 2: Orthogonal A-Frame */
  /*---------------------------------*/

  cs_real_t  *weighting = NULL, *offsetting = NULL;

  /* Only defined on interior faces */
  BFT_MALLOC(working_array, 2*n_i_faces, cs_real_t);

  for (i = 0; i < 2*n_i_faces; i++)
    working_array[i] = 0.;

  weighting = working_array;
  offsetting = working_array + n_i_faces;

  _compute_weighting_offsetting(mesh,
                                mesh_quantities,
                                weighting,
                                offsetting,
                                bad_cell_flag);

  ibad = 0;
  for (i = 0; i < n_cells; i++) {
    if (bad_cell_flag[i] & _type_flag_mask[1]) {
      ibad++;
      iwarning++;
    }
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_counter(&ibad, 1);
    cs_parall_counter(&iwarning, 1);
  }

  /* Display listing output */
  bft_printf(_("\n  Criteria 2: Offset:\n"));
  bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
             (unsigned long long)ibad,
             (double)ibad / (double)n_cells_tot * 100.0);

  BFT_FREE(working_array);

  /* Condition 3: Least Squares Gradient */
  /*-------------------------------------*/

  cs_real_t  *lsq = NULL;

  BFT_MALLOC(working_array, n_cells_wghosts, cs_real_t);

  for (i = 0; i < n_cells_wghosts; i++)
    working_array[i] = 0.;

  lsq = working_array;

  _compute_least_squares(mesh,
                         mesh_quantities,
                         lsq,
                         bad_cell_flag);

  ibad = 0;
  for (i = 0; i < n_cells; i++) {
    if (bad_cell_flag[i] & _type_flag_mask[2]) {
      ibad++;
      iwarning++;
    }
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_counter(&ibad, 1);
    cs_parall_counter(&iwarning, 1);
  }

  /* Display log output */
  bft_printf(_("\n  Criteria 3: Least-Squares Gradient Quality:\n"));
  bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
             (unsigned long long)ibad,
             (double)ibad / (double)n_cells_tot * 100.0);

  BFT_FREE(working_array);

  /* Condition 4: Volume Ratio */
  /*---------------------------*/

  cs_real_t  *vol_ratio = NULL;

  BFT_MALLOC(working_array, n_i_faces, cs_real_t);

  for (i = 0; i < n_i_faces; i++)
    working_array[i] = 0.;

  vol_ratio = working_array;

  _compute_volume_ratio(mesh,
                        mesh_quantities,
                        vol_ratio,
                        bad_cell_flag);

  ibad = 0;
  for (i = 0; i < n_cells; i++) {
    if (bad_cell_flag[i] & _type_flag_mask[3]) {
      ibad++;
      iwarning++;
    }
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_counter(&ibad, 1);
    cs_parall_counter(&iwarning, 1);
  }

  /* Display listing output */
  bft_printf(_("\n  Criteria 4: Cells Volume Ratio:\n"));
  bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
             (unsigned long long)ibad,
             (double)ibad / (double)n_cells_tot * 100.0);

  BFT_FREE(working_array);

  /* Guilt by association */
  /*----------------------*/

  cs_lnum_t face_id, cell1, cell2;

  cs_lnum_t  *bad_guilt_cells = NULL;

  BFT_MALLOC(bad_guilt_cells, n_cells_wghosts, cs_lnum_t);

  for (i = 0; i < n_cells_wghosts; i++)
    bad_guilt_cells[i] = 0;

  /* Loop on interior faces */
  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[2 * face_id] - 1;
    cell2 = mesh->i_face_cells[2 * face_id + 1] - 1;

    if (bad_cell_flag[cell2] != 0)
      bad_guilt_cells[cell1]++;
  }

  ibad = 0;
  for (i = 0; i < n_cells; i++) {
    if (bad_guilt_cells[i] >= 5 && bad_cell_flag[i] == 0) {
      ibad++;
      iwarning++;
      bad_cell_flag[i] = _type_flag_mask[4];
    }
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_counter(&ibad, 1);
    cs_parall_counter(&iwarning, 1);
  }

  /* Display listing output */
  bft_printf(_("\n  Criteria 5: Guilt by Association:\n"));
  bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
             (unsigned long long)ibad,
             (double)ibad / (double)n_cells_tot * 100.0);

  /* Warning printed in the log file */
  /*---------------------------------*/

  if (iwarning > 0) {
    bft_printf
      (_("\n Warning:\n"
         " --------\n"
         "    Mesh quality issue has been detected\n\n"
         "    The mesh should be re-considered using the listed criteria.\n\n"
         "    The calculation will run but the solution quality may be"
         " degraded...\n"));
  }

  BFT_FREE(bad_guilt_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Post-process bad cell quality indicators.
 *
 * \param [in]  mesh             pointer to associated mesh structure
 * \param [in]  mesh_quantities  pointer to associated mesh quantities
 *                               structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_postprocess(const cs_mesh_t             *mesh,
                              const cs_mesh_quantities_t  *mesh_quantities)
{
  int i;

  cs_lnum_t  *bad_cells_v = NULL;

  const cs_lnum_t  n_cells       = mesh->n_cells;
  const unsigned  *bad_cell_flag = mesh_quantities->bad_cell_flag;

  const unsigned criterion[] = {CS_BAD_CELL_ORTHO_NORM,
                                CS_BAD_CELL_OFFSET,
                                CS_BAD_CELL_LSQ_GRAD,
                                CS_BAD_CELL_RATIO,
                                CS_BAD_CELL_GUILT,
                                CS_BAD_CELL_USER};

  const char *criterion_name[] = {N_("Bad Cell Ortho Norm"),
                                  N_("Bad Cell Offset"),
                                  N_("Bad Cell LSQ Gradient"),
                                  N_("Bad Cell Volume Ratio"),
                                  N_("Bad Cell Association"),
                                  N_("Bad Cell by User")};

  const int n_criteria = 6;

  /* Activate default post-processing writer */

  cs_post_activate_writer(-1, true);

  BFT_MALLOC(bad_cells_v, n_cells, int);

  /* Loop on criteria */
  /*------------------*/

  for (i = 0; i < n_criteria; i++) {

    cs_lnum_t j;
    cs_lnum_t crit_flag = 0;

    for (j = 0; j < n_cells; j++) {
      if (bad_cell_flag[j] & criterion[i]) {
        bad_cells_v[j] = 1;
        crit_flag = 1;
      }
      else
        bad_cells_v[j] = 0;
    }

    cs_parall_counter_max(&crit_flag, 1);

    if (crit_flag > 0)
      cs_post_write_var(-1,
                        _(criterion_name[i]),
                        1,
                        false,
                        true,
                        CS_POST_TYPE_int,
                        -1,
                        0.0,
                        bad_cells_v,
                        NULL,
                        NULL);

  }

  BFT_FREE(bad_cells_v);
}

/*----------------------------------------------------------------------------*/

/* Delete local macros */

#undef _CROSS_PRODUCT_3D
#undef _DOT_PRODUCT_3D
#undef _MODULE_3D
#undef _COSINE_3D

/*----------------------------------------------------------------------------*/

END_C_DECLS

