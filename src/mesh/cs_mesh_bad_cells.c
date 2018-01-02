/*============================================================================
 * \file Detect bad cells within meshes.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_halo.h"
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

/*----------------------------------------------------------------------------*/
/*! \file cs_mesh_bad_cells.c
 *
 * \brief Detect bad cells within meshes.
 *
 * Please refer to the
 * <a href="../../theory.pdf#badcells"><b>flagging of bad cells</b></a>
 * section of the theory guide for more informations.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
 * Static global variables
 *============================================================================*/

/* compute and visualize flags (-1 initially, mask afterwards;
   first value: at initialization: second value: at each time step) */

static int  _type_flag_compute[] = {-1, 0};
static int  _type_flag_visualize[] = {0, 0};
static int  _call_type_compute = 0;
static int  _call_type_visualize = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
                    unsigned                     bad_cell_flag[])
{
  cs_lnum_t  i, face_id, cell1, cell2;

  cs_real_t  cell_center1[3], cell_center2[3];
  cs_real_t  face_center[3];
  cs_real_t  face_normal[3], vect[3];

  double  cos_alpha, i_face_ortho, b_face_ortho;

  const cs_lnum_t  dim = 3;

  /* Loop on interior faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[face_id][0];
    cell2 = mesh->i_face_cells[face_id][1];

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

    i_face_ortho = cos_alpha;

    if (i_face_ortho < 0.1) {
      bad_cell_flag[cell1] |= CS_BAD_CELL_ORTHO_NORM;
      bad_cell_flag[cell2] |= CS_BAD_CELL_ORTHO_NORM;
    }
  }

  /* Loop on boundary faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {

    cell1 = mesh->b_face_cells[face_id];

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

    b_face_ortho = cos_alpha;

    if (b_face_ortho < 0.1)
      bad_cell_flag[cell1] |= CS_BAD_CELL_ORTHO_NORM;
  }

  if (mesh->halo != NULL)
    cs_halo_sync_untyped(mesh->halo,
                         CS_HALO_EXTENDED,
                         sizeof(unsigned),
                         bad_cell_flag);
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
_compute_offsetting(const cs_mesh_t             *mesh,
                    const cs_mesh_quantities_t  *mesh_quantities,
                    unsigned                     bad_cell_flag[])
{
  cs_lnum_t  face_id, cell1, cell2;
  double  of_n, off_1, off_2;

  const cs_real_t *v_of, *v_n;

  /* Loop on interior faces */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[face_id][0];
    cell2 = mesh->i_face_cells[face_id][1];

    /* Compute center offsetting coefficient,
       in a manner consistent with iterative gradient reconstruction */

    v_of = &(mesh_quantities->dofij[face_id*3]);
    v_n = &(mesh_quantities->i_face_normal[face_id*3]);
    of_n = _MODULE_3D(v_of) * _MODULE_3D(v_n);

    off_1 = 1 - pow(of_n / mesh_quantities->cell_vol[cell1], 1/3.);
    off_2 = 1 - pow(of_n / mesh_quantities->cell_vol[cell2], 1/3.);

    if (off_1 < 0.1)
      bad_cell_flag[cell1] = bad_cell_flag[cell1] | CS_BAD_CELL_OFFSET;
    if (off_2 < 0.1)
      bad_cell_flag[cell2] = bad_cell_flag[cell2] | CS_BAD_CELL_OFFSET;

  }

  if (mesh->halo != NULL)
    cs_halo_sync_untyped(mesh->halo,
                         CS_HALO_EXTENDED,
                         sizeof(unsigned),
                         bad_cell_flag);
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
                       unsigned                     bad_cell_flag[])
{
  cs_real_t lsq;
  const cs_lnum_t  dim = mesh->dim;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;

  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mesh_quantities->b_face_normal;

  cs_lnum_t     i, k, face_id, cell1, cell2, cell_id;
  cs_real_3_t   cell_center1, cell_center2, vect, dij, eigenvalues;
  cs_real_33_t  w2;

  double unsdij, surfn, surf_n_inv, min_diag, max_diag;
  double xam, q, p, r, phi;

  cs_real_t *w1 = NULL;

  const double pi = 4 * atan(1);

  BFT_MALLOC(w1, 6 * n_cells_wghosts, cs_real_t);

  for (i = 0; i < 6 * n_cells_wghosts; i++)
    w1[i] = 0.;

  /* Loop on interior faces */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[face_id][0];
    cell2 = mesh->i_face_cells[face_id][1];

    /* Center of gravity for each cell */

    for (i = 0; i < dim; i++) {
      cell_center1[i] = mesh_quantities->cell_cen[cell1*dim + i];
      cell_center2[i] = mesh_quantities->cell_cen[cell2*dim + i];
      vect[i] = cell_center2[i] - cell_center1[i];
    }

    unsdij = 1.0 / _MODULE_3D(vect);

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
    cell1 = mesh->b_face_cells[face_id];

    surfn = _MODULE_3D(b_face_normal[face_id]);

    surf_n_inv = 1.0 / surfn;

    for (i = 0; i < dim; i++)
      dij[i] = b_face_normal[face_id][i] * surf_n_inv;

    w1[cell1] += dij[0] * dij[0];
    w1[cell1 + n_cells_wghosts] += dij[1] * dij[1];
    w1[cell1 + 2 * n_cells_wghosts] += dij[2] * dij[2];
    w1[cell1 + 3 * n_cells_wghosts] += dij[0] * dij[1];
    w1[cell1 + 4 * n_cells_wghosts] += dij[0] * dij[2];
    w1[cell1 + 5 * n_cells_wghosts] += dij[1] * dij[2];
  }

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

    /* Compute the eigenvalues for a given real symmetric 3x3 matrix */

    xam = w2[0][1] * w2[0][1] + w2[0][2] * w2[0][2] + w2[1][2] * w2[1][2];

    /* First check if the matrix is diagonal */
    if (xam <= 0.) {
      for (i = 0; i < 3; i++)
        eigenvalues[i] = w2[i][i];
    }

    /* If the matrix is not diagonal, we get the eigenvalues from a
       trigonometric solution                                       */
    else {
      q = (w2[0][0] + w2[1][1] + w2[2][2]) / 3.;

      p = (w2[0][0] - q) * (w2[0][0] - q) +
          (w2[1][1] - q) * (w2[1][1] - q) +
          (w2[2][2] - q) * (w2[2][2] - q) + 2. * xam;

      p = sqrt(p / 6.);

      for (i = 0; i < 3; i++) {
        for (k = 0; k < 3; k++) {
          if (i == k)
            w2[i][k] = (1. / p) * (w2[i][k] - q);
          else
            w2[i][k] = (1. / p) * (w2[i][k]);
        }
      }

      r =   w2[0][0] * w2[1][1] * w2[2][2]
          + w2[0][1] * w2[1][2] * w2[2][0]
          + w2[0][2] * w2[1][0] * w2[2][1]
          - w2[0][2] * w2[1][1] * w2[2][0]
          - w2[0][1] * w2[1][0] * w2[2][2]
          - w2[0][0] * w2[1][2] * w2[2][1];

      r *= 0.5;

      /* In exact arithmetic for a symmetric matrix  -1 <= r <= 1
         but computation error can leave it slightly outside this range */
      if (r <= -1.)
        phi = pi / 3.;
      else if (r >= 1.)
        phi = 0.;
      else
        phi = acos(r) / 3.;

      /* The eigenvalues satisfy eig3 <= eig2 <= eig1
         with tr(w2) = eig1 + eig2 + eig3             */
      eigenvalues[0] = q + 2. * p * cos(phi);
      eigenvalues[2] = q + 2. * p * cos(phi + (2. * pi / 3.));
      eigenvalues[1] = 3. * q - eigenvalues[0] - eigenvalues[2];
    }

    min_diag = 1.e15;
    max_diag = 0.;

    for (i = 0; i < 3; i++) {
      min_diag = fmin(min_diag, fabs(eigenvalues[i]));
      max_diag = fmax(max_diag, fabs(eigenvalues[i]));
    }

    lsq = min_diag / max_diag;

    if (lsq < 0.1)
      bad_cell_flag[cell_id] = bad_cell_flag[cell_id] | CS_BAD_CELL_LSQ_GRAD;
  }

  BFT_FREE(w1);

  if (mesh->halo != NULL)
    cs_halo_sync_untyped(mesh->halo,
                         CS_HALO_EXTENDED,
                         sizeof(unsigned),
                         bad_cell_flag);
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
                      unsigned                     bad_cell_flag[])
{
  double vol_ratio;
  cs_real_t *volume = mesh_quantities->cell_vol;

  cs_lnum_t   face_id, cell1, cell2;

  /* Loop on interior faces */
  /*------------------------*/

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell1 = mesh->i_face_cells[face_id][0];
    cell2 = mesh->i_face_cells[face_id][1];

    vol_ratio = fmin(volume[cell1] / volume[cell2],
                     volume[cell2] / volume[cell1]);

    if (vol_ratio < 0.1*0.1) {
      bad_cell_flag[cell1] = bad_cell_flag[cell1] | CS_BAD_CELL_RATIO;
      bad_cell_flag[cell2] = bad_cell_flag[cell2] | CS_BAD_CELL_RATIO;
    }
  }

  if (mesh->halo != NULL)
    cs_halo_sync_untyped(mesh->halo,
                         CS_HALO_EXTENDED,
                         sizeof(unsigned),
                         bad_cell_flag);
}

/*----------------------------------------------------------------------------
 * Post-process bad cell quality indicators.
 *
 * parameters:
 *   mesh            <-- pointer to a mesh structure.
 *   mesh_quantities <-- pointer to a mesh quantities structures.
 *   call_type       <-- visualization type id (0: fixed; 1: time varying)
 *   ts              <-- time step structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_bad_cells_post(const cs_mesh_t             *mesh,
                const cs_mesh_quantities_t  *mesh_quantities,
                int                          call_type,
                const cs_time_step_t        *ts)
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

  if (_type_flag_visualize[call_type] == 0)
    return;

  BFT_MALLOC(bad_cells_v, n_cells, int);

  /* Loop on criteria */
  /*------------------*/

  for (i = 0; i < n_criteria; i++) {

    if (_type_flag_visualize[call_type] & criterion[i]) {

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
        cs_post_write_var(CS_POST_MESH_VOLUME,
                          CS_POST_WRITER_ALL_ASSOCIATED,
                          _(criterion_name[i]),
                          1,
                          false,
                          true,
                          CS_POST_TYPE_int,
                          bad_cells_v,
                          NULL,
                          NULL,
                          ts);

    }

  }

  BFT_FREE(bad_cells_v);
}

/*----------------------------------------------------------------------------
 * Post-process bad cell quality indicators.
 *
 * parameters:
 *   mesh  <--  Void pointer to associated mesh structure
 *   ts    <-- time step structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_bad_cells_post_function(void                  *mesh,
                         const cs_time_step_t  *ts)
{
  /* TODO: enable this function with other meshes thant the
     global mesh (will be easier when mesh_quantities becomes a member
     of mesh). */

  if (mesh != cs_glob_mesh)
    return;

  _bad_cells_post(mesh,
                  cs_glob_mesh_quantities,
                  1,
                  ts);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define which cell quality indicators are used and when.
 *
 * \note
 * We assume that if a given criterion is computed at each time
 * step, it is also computed at initialization, but for visualization,
 * it is either one or the other, as visualization formats and tools
 * may not always accept both a fixed and time-varying instance of a
 * given variable.
 *
 * \param[in]   type_flag_mask   criterion type mask (0 for all)
 * \param[in]   compute          0: never compute;
 *                               1: compute at initialization;
 *                               2: compute at each time step
 * \param[in]   visualize        0: never visualize
 *                               1: visualize at initialization;
 *                               2: visualize at each time step
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_set_options(int  type_flag_mask,
                              int  compute,
                              int  visualize)
{
  int i;

  for (i = 0; i < 2; i++) {
    _type_flag_compute[i] = 0;
    _type_flag_visualize[i] = 0;
  }

  for (i = 0; i < 6; i++) {
    int mask = (1 << i);
    if (type_flag_mask == 0 || (type_flag_mask & mask)) {
      if (compute > 0) {
        _type_flag_compute[0] = _type_flag_compute[0] | mask;
        if (visualize == 1)
          _type_flag_visualize[0] = _type_flag_visualize[0] | mask;
        if (compute > 1) {
          _type_flag_compute[1] = _type_flag_compute[1] | mask;
          if (visualize > 1)
            _type_flag_visualize[1] = _type_flag_visualize[1] | mask;
        }
      }
    }
  }

  /* Register post processing function if required */

  if (_type_flag_visualize[1] != 0)
    cs_post_add_time_dep_output(_bad_cells_post_function,
                                (void *)cs_glob_mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate which cell quality indicators are used and when.
 *
 * Each array is optional, and returns 2 flags; the first flag is used at
 * initialization, the second one at each time step.
 *
 * A flag is a mask to be compared using an "and" (&) operation with a given
 * criteria type mask (CS_BAD_CELL_ORTHO_NORM, CS_BAD_CELL_OFFSET, ...).
 *
 * \param [out]  compute    computation mask (initialization, per time step),
 *                          or NULL
 * \param [out]  visualize  visualization mask (initialization, per time step),
                            or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_get_options(int  compute[2],
                              int  visualize[2])
{
  if (_type_flag_compute[0] < 0)  /* Set default if not done yet */
    cs_mesh_bad_cells_set_options(0, 1, 1);

  if (compute != NULL) {
    compute[0] = _type_flag_compute[0] = 0;
    compute[1] = _type_flag_compute[1] = 0;
  }

  if (visualize != NULL) {
    visualize[0] = _type_flag_visualize[0] = 0;
    visualize[1] = _type_flag_visualize[1] = 0;
  }
}

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
  unsigned *bad_cell_flag = NULL;

  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;

  cs_lnum_t i;
  cs_gnum_t n_cells_tot, iwarning, ibad;

  const int call_type_log = _call_type_compute;

  /* If bad cell data has been destroyed and this function is
     called, we have a call type 0, even if it was called before */

  if (mesh_quantities->bad_cell_flag == NULL)
    _call_type_compute = 0;

  /* Initialization or per time step ? */

  const int call_type = _call_type_compute;

  /* Set defaults if not done yet */

  if (_type_flag_compute[0] < 0)
    cs_mesh_bad_cells_set_options(0, 1, 1);

  if (_type_flag_compute[call_type] == 0)
    return;

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

  if (_type_flag_compute[call_type] & CS_BAD_CELL_ORTHO_NORM)
    _compute_ortho_norm(mesh,
                        mesh_quantities,
                        bad_cell_flag);


  if (_type_flag_compute[call_type_log] & CS_BAD_CELL_ORTHO_NORM) {

    ibad = 0;
    for (i = 0; i < n_cells; i++) {
      if (bad_cell_flag[i] & CS_BAD_CELL_ORTHO_NORM) {
        ibad++;
        iwarning++;
      }
    }

    if (cs_glob_rank_id >= 0) {
      cs_parall_counter(&ibad, 1);
      cs_parall_counter(&iwarning, 1);
    }

    /* Display log output */
    bft_printf(_("\n  Criterion 1: Orthogonality:\n"));
    bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
               (unsigned long long)ibad,
               (double)ibad / (double)n_cells_tot * 100.0);

  }

  /* Condition 2: Orthogonal A-Frame */
  /*---------------------------------*/

  if (   _type_flag_compute[call_type] & CS_BAD_CELL_OFFSET
      && cs_glob_mesh_quantities->min_vol >= 0.)
    _compute_offsetting(mesh,
                        mesh_quantities,
                        bad_cell_flag);

  if (   _type_flag_compute[call_type_log] & CS_BAD_CELL_OFFSET
      && cs_glob_mesh_quantities->min_vol >= 0.) {

    ibad = 0;
    for (i = 0; i < n_cells; i++) {
      if (bad_cell_flag[i] & CS_BAD_CELL_OFFSET) {
        ibad++;
        iwarning++;
      }
    }

    if (cs_glob_rank_id >= 0) {
      cs_parall_counter(&ibad, 1);
      cs_parall_counter(&iwarning, 1);
    }

    /* Display listing output */
    bft_printf(_("\n  Criterion 2: Offset:\n"));
    bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
               (unsigned long long)ibad,
               (double)ibad / (double)n_cells_tot * 100.0);

  }

  /* Condition 3: Least Squares Gradient */
  /*-------------------------------------*/

  if (_type_flag_compute[call_type] & CS_BAD_CELL_LSQ_GRAD)
    _compute_least_squares(mesh,
                           mesh_quantities,
                           bad_cell_flag);

  if (_type_flag_compute[call_type_log] & CS_BAD_CELL_LSQ_GRAD) {

    ibad = 0;
    for (i = 0; i < n_cells; i++) {
      if (bad_cell_flag[i] & CS_BAD_CELL_LSQ_GRAD) {
        ibad++;
        iwarning++;
      }
    }

    if (cs_glob_rank_id >= 0) {
      cs_parall_counter(&ibad, 1);
      cs_parall_counter(&iwarning, 1);
    }

    /* Display log output */
    bft_printf(_("\n  Criterion 3: Least-Squares Gradient Quality:\n"));
    bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
               (unsigned long long)ibad,
               (double)ibad / (double)n_cells_tot * 100.0);

  }

  /* Condition 4: Volume Ratio */
  /*---------------------------*/

  if (_type_flag_compute[call_type] & CS_BAD_CELL_RATIO)
    _compute_volume_ratio(mesh,
                          mesh_quantities,
                          bad_cell_flag);

  if (_type_flag_compute[call_type_log] & CS_BAD_CELL_RATIO) {

    ibad = 0;
    for (i = 0; i < n_cells; i++) {
      if (bad_cell_flag[i] & CS_BAD_CELL_RATIO) {
        ibad++;
        iwarning++;
      }
    }

    if (cs_glob_rank_id >= 0) {
      cs_parall_counter(&ibad, 1);
      cs_parall_counter(&iwarning, 1);
    }

    /* Display listing output */
    bft_printf(_("\n  Criterion 4: Cells Volume Ratio:\n"));
    bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
               (unsigned long long)ibad,
               (double)ibad / (double)n_cells_tot * 100.0);

  }

  /* Guilt by association */
  /*----------------------*/

  if (_type_flag_compute[call_type] & CS_BAD_CELL_GUILT) {

    cs_lnum_t face_id, cell1, cell2;

    cs_lnum_t  *bad_guilt_cells = NULL;

    BFT_MALLOC(bad_guilt_cells, n_cells_wghosts, cs_lnum_t);

    for (i = 0; i < n_cells_wghosts; i++)
      bad_guilt_cells[i] = 0;

    /* Loop on interior faces */
    for (face_id = 0; face_id < n_i_faces; face_id++) {

      cell1 = mesh->i_face_cells[face_id][0];
      cell2 = mesh->i_face_cells[face_id][1];

      if (bad_cell_flag[cell2] != 0)
        bad_guilt_cells[cell1]++;
    }

    ibad = 0;
    for (i = 0; i < n_cells; i++) {
      if (bad_guilt_cells[i] >= 5 && bad_cell_flag[i] == 0) {
        ibad++;
        iwarning++;
        bad_cell_flag[i] = CS_BAD_CELL_GUILT;
      }
    }

    BFT_FREE(bad_guilt_cells);

    if (_type_flag_compute[call_type_log] & CS_BAD_CELL_GUILT) {

      if (cs_glob_rank_id >= 0) {
        cs_parall_counter(&ibad, 1);
        cs_parall_counter(&iwarning, 1);
      }

      /* Display listing output */
      bft_printf(_("\n  Criterion 5: Guilt by Association:\n"));
      bft_printf(_("    Number of bad cells detected: %llu --> %3.0f %%\n"),
                 (unsigned long long)ibad,
                 (double)ibad / (double)n_cells_tot * 100.0);

    }

  }

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

  /* After first call, we assume others are done at each time step */

  _call_type_compute = 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Post-process time-independent bad cell quality indicators.
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
  /* Initialization or per time step ? */

  const int call_type = _call_type_visualize;

  /* Set defaults if not done yet */

  if (_type_flag_visualize[0] < 0)
    cs_mesh_bad_cells_set_options(0, 1, 1);

  if (_type_flag_visualize[call_type] == 0)
    return;

  cs_post_activate_writer(-1, true);

  _bad_cells_post(mesh,
                  mesh_quantities,
                  0,
                  NULL);

  _call_type_visualize = 1; /* Prevent future calls from doing anything */
}

/*----------------------------------------------------------------------------*/

/* Delete local macros */

#undef _CROSS_PRODUCT_3D
#undef _DOT_PRODUCT_3D
#undef _MODULE_3D
#undef _COSINE_3D

/*----------------------------------------------------------------------------*/

END_C_DECLS
