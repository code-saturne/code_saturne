/*============================================================================
 * Management of mesh quantities
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
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_bad_cells_regularisation.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_mesh_quantities.c
 *
 * \brief Management of mesh quantities.
 *
 * Please refer to the
 * <a href="../../theory.pdf#meshquantities"><b>geometric quantities</b></a>
 * section of the theory guide for more informations.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to cs_mesh_quantities_t structure for the main mesh */

cs_mesh_quantities_t  *cs_glob_mesh_quantities = NULL;

/* Choice of the algorithm for computing gravity centers of the cells */

static int _cell_cen_algorithm = 0;
static int _ajust_face_cog_compat_v11_v52 = 0;

/* Choice of the option for computing cocg
   (iterative or least squares method) or not */

static bool _compute_cocg_s_it = false;
static bool _compute_cocg_it = false;
static bool _compute_cocg_lsq = false;

/* Flag (mask) to activate bad cells correction
 * CS_BAD_CELLS_WARPED_CORRECTION
 * CS_FACE_DISTANCE_CLIP
 * CS_FACE_RECONSTRUCTION_CLIP
 * are set as default options
 * */
unsigned cs_glob_mesh_quantities_flag =
  CS_BAD_CELLS_WARPED_CORRECTION + CS_FACE_DISTANCE_CLIP
  + CS_FACE_RECONSTRUCTION_CLIP;

/* Choice of the porous model */
int cs_glob_porous_model = 0;

/* Number of computation updates */

static int _n_computations = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient iterative algorithm
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *----------------------------------------------------------------------------*/

static void
_compute_cell_cocg_s_it(const cs_mesh_t         *m,
                        cs_mesh_quantities_t    *fvq)
{
  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;

  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;

  cs_real_33_t   *restrict cocgb;
  cs_real_33_t   *restrict cocg;
  cocg = fvq->cocg_s_it;
  cocgb = fvq->cocgb_s_it;

  if (cocg == NULL) {
    BFT_MALLOC(cocg, n_cells_ext, cs_real_33_t);
    BFT_MALLOC(cocgb, m->n_b_cells, cs_real_33_t);
    fvq->cocgb_s_it = cocgb;
    fvq->cocg_s_it = cocg;
  }

  /* Compute cocg */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    cocg[cell_id][0][0] = cell_vol[cell_id];
    cocg[cell_id][0][1] = 0.0;
    cocg[cell_id][0][2] = 0.0;
    cocg[cell_id][1][0] = 0.0;
    cocg[cell_id][1][1] = cell_vol[cell_id];
    cocg[cell_id][1][2] = 0.0;
    cocg[cell_id][2][0] = 0.0;
    cocg[cell_id][2][1] = 0.0;
    cocg[cell_id][2][2] = cell_vol[cell_id];
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_real_t  fctb[4];

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++) {
            fctb[mm] = -dofij[face_id][mm] * 0.5 * i_face_normal[face_id][ll];
            cocg[ii][ll][mm] += fctb[mm];
            cocg[jj][ll][mm] -= fctb[mm];
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Save partial cocg at interior faces of boundary cells */

# pragma omp parallel for if (m->n_b_cells > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < m->n_b_cells; ii++) {
    cs_lnum_t cell_id = m->b_cells[ii];
    for (cs_lnum_t ll = 0; ll < 3; ll++) {
      for (cs_lnum_t mm = 0; mm < 3; mm++)
        cocgb[ii][ll][mm] = cocg[cell_id][ll][mm];
    }
  }

  /* Invert for all cells. */
  /*-----------------------*/

  /* The cocg term for interior cells only changes if the mesh does */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    cs_math_33_inv_cramer_in_place(cocg[cell_id]);
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient least squares algorithm
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *   ce   <->  coupling entity
 *----------------------------------------------------------------------------*/

static void
_compute_cell_cocg_lsq(const cs_mesh_t        *m,
                       cs_mesh_quantities_t   *fvq,
                       cs_internal_coupling_t *ce)
{
  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;

  cs_real_33_t   *restrict cocgb;
  if (ce == NULL) {
    cocgb = fvq->cocgb_s_lsq;
  } else {
    cocgb = ce->cocgb_s_lsq;
  }
  cs_real_33_t   *restrict cocg = fvq->cocg_lsq;

  const bool *coupled_faces;
  if (ce != NULL) {
    coupled_faces = ce->coupled_faces;
  }

  if (ce == NULL) {
    if (cocg == NULL) {
      BFT_MALLOC(cocg, n_cells_ext, cs_real_33_t);
      fvq->cocg_lsq = cocg;
    }
    if (cocgb == NULL) {
      BFT_MALLOC(cocgb, m->n_b_cells, cs_real_33_t);
      fvq->cocgb_s_lsq = cocgb;
    }
  } else if (ce != NULL) {
    if (cocgb == NULL) {
      BFT_MALLOC(cocgb, m->n_b_cells, cs_real_33_t);
      ce->cocgb_s_lsq = cocgb;
    }
  }

  /* Initialization */

# pragma omp parallel
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (cs_lnum_t ll = 0; ll < 3; ll++) {
      for (cs_lnum_t mm = 0; mm < 3; mm++)
        cocg[cell_id][ll][mm] = 0.0;
    }
  }

  /* Contribution from interior faces */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_real_t dc[3];
        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];
        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++)
            cocg[ii][ll][mm] += dc[ll] * dc[mm] * ddc;
        }
        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++)
            cocg[jj][ll][mm] += dc[ll] * dc[mm] * ddc;
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution for internal coupling */
  if (ce != NULL) {
    cs_internal_coupling_lsq_cocg_contribution(ce, cocg);
  }

  /* Contribution from extended neighborhood */

  if (m->halo_type == CS_HALO_EXTENDED) {

    /* Not compatible with internal coupling */
    if (ce != NULL) {
      bft_error(__FILE__, __LINE__, 0,
                "Extended least-square gradient reconstruction \
                 is not supported with internal coupling");
    }

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
      for (cs_lnum_t cidx = cell_cells_idx[ii];
           cidx < cell_cells_idx[ii+1];
           cidx++) {

        cs_real_t dc[3];
        cs_lnum_t jj = cell_cells_lst[cidx];

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dc[ll] = cell_cen[jj][ll] - cell_cen[ii][ll];
        cs_real_t ddc = 1. / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (cs_lnum_t ll = 0; ll < 3; ll++) {
          for (cs_lnum_t mm = 0; mm < 3; mm++)
            cocg[ii][ll][mm] += dc[ll] * dc[mm] * ddc;
        }

      }
    }

  } /* End for extended neighborhood */

  /* Save partial cocg at interior faces of boundary cells */

# pragma omp parallel for
  for (cs_lnum_t ii = 0; ii < m->n_b_cells; ii++) {
    cs_lnum_t cell_id = m->b_cells[ii];
    for (cs_lnum_t ll = 0; ll < 3; ll++) {
      for (cs_lnum_t mm = 0; mm < 3; mm++)
        cocgb[ii][ll][mm] = cocg[cell_id][ll][mm];
    }
  }

  /* Contribution from boundary faces, assuming symmetry everywhere
     so as to avoid obtaining a non-invertible matrix in 2D cases. */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        if (ce == NULL || !coupled_faces[face_id]) {

          cs_lnum_t ii = b_face_cells[face_id];

          cs_real_t udbfs = 1. / b_face_surf[face_id];

          cs_real_t dddij[3];
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            dddij[ll] =   udbfs * b_face_normal[face_id][ll];

          for (cs_lnum_t ll = 0; ll < 3; ll++) {
            for (cs_lnum_t mm = 0; mm < 3; mm++)
              cocg[ii][ll][mm] += dddij[ll]*dddij[mm];
          }

        } /* face without internal coupling */

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Invert for all cells. */
  /*-----------------------*/

  /* The cocg term for interior cells only changes if the mesh does */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    cs_math_33_inv_cramer_in_place(cocg[cell_id]);
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the iterative algorithm
 *
 * parameters:
 *   m               <--  mesh
 *   fvq             <->  mesh quantities
 *   ce              <->  coupling entity
 *----------------------------------------------------------------------------*/

static void
_compute_cell_cocg_it(const cs_mesh_t        *m,
                      cs_mesh_quantities_t   *fvq,
                      cs_internal_coupling_t *ce)
{
  /* Local variables */

  const int n_cells = m->n_cells;
  const int n_cells_with_ghosts = m->n_cells_with_ghosts;
  const int n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;

  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;
  cs_real_33_t *restrict cocg
    = fvq->cocg_it;

  if (ce == NULL) {
    cocg = fvq->cocg_it;
  } else {
    cocg = ce->cocg_it;
  }

  cs_lnum_t  cell_id, face_id, i, j;
  cs_real_t  pfac, vecfac;
  cs_real_t  dvol1, dvol2;

  if (cocg == NULL) {
    BFT_MALLOC(cocg, n_cells_with_ghosts, cs_real_33_t);
    if (ce == NULL)
      fvq->cocg_it = cocg;
    else
      ce->cocg_it = cocg;
  }

  /* compute the dimensionless matrix COCG for each cell*/

  for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    cocg[cell_id][0][0]= 1.0;
    cocg[cell_id][0][1]= 0.0;
    cocg[cell_id][0][2]= 0.0;
    cocg[cell_id][1][0]= 0.0;
    cocg[cell_id][1][1]= 1.0;
    cocg[cell_id][1][2]= 0.0;
    cocg[cell_id][2][0]= 0.0;
    cocg[cell_id][2][1]= 0.0;
    cocg[cell_id][2][2]= 1.0;
  }

  /* Interior face treatment */

  for (face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    dvol1 = 1./cell_vol[cell_id1];
    dvol2 = 1./cell_vol[cell_id2];

    for (i = 0; i < 3; i++) {

      pfac = -0.5*dofij[face_id][i];

      for (j = 0; j < 3; j++) {
        vecfac = pfac*i_face_normal[face_id][j];
        cocg[cell_id1][i][j] += vecfac * dvol1;
        cocg[cell_id2][i][j] -= vecfac * dvol2;
      }
    }
  }

  /* Contribution for internal coupling */
  if (ce != NULL) {
    cs_internal_coupling_it_cocg_contribution(ce, cocg);
  }

  /* 3x3 Matrix inversion */

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells; cell_id++)
    cs_math_33_inv_cramer_in_place(cocg[cell_id]);
}

/*----------------------------------------------------------------------------
 * Build the geometrical matrix linear gradient correction
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *----------------------------------------------------------------------------*/

static void
_compute_corr_grad_lin(const cs_mesh_t       *m,
                       cs_mesh_quantities_t  *fvq)
{
  /* Local variables */

  const int n_cells = m->n_cells;
  const int n_cells_with_ghosts = m->n_cells_with_ghosts;
  const int n_i_faces = m->n_i_faces;
  const int n_b_faces = m->n_b_faces;

  const cs_lnum_t  *b_face_cells = m->b_face_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;

  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;

  if (fvq->corr_grad_lin_det == NULL)
    BFT_MALLOC(fvq->corr_grad_lin_det, n_cells_with_ghosts, cs_real_t);
  if (fvq->corr_grad_lin == NULL)
    BFT_MALLOC(fvq->corr_grad_lin, n_cells_with_ghosts, cs_real_33_t);

  cs_real_t    *restrict corr_grad_lin_det = fvq->corr_grad_lin_det;
  cs_real_33_t *restrict corr_grad_lin     = fvq->corr_grad_lin;

  /* Initialization */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        corr_grad_lin[cell_id][i][j] = 0.;
    }
  }

  /* Internal faces contribution */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_real_t flux = i_face_cog[face_id][i] * i_face_normal[face_id][j];
        corr_grad_lin[cell_id1][i][j] += flux;
        corr_grad_lin[cell_id2][i][j] -= flux;
      }
    }
  }

  /* Boundary faces contribution */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t cell_id = b_face_cells[face_id];
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_real_t flux = b_face_cog[face_id][i] * b_face_normal[face_id][j];
        corr_grad_lin[cell_id][i][j] += flux;
      }
    }
  }

  /* Matrix inversion */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    double cocg11 = corr_grad_lin[cell_id][0][0] / cell_vol[cell_id];
    double cocg12 = corr_grad_lin[cell_id][1][0] / cell_vol[cell_id];
    double cocg13 = corr_grad_lin[cell_id][2][0] / cell_vol[cell_id];
    double cocg21 = corr_grad_lin[cell_id][0][1] / cell_vol[cell_id];
    double cocg22 = corr_grad_lin[cell_id][1][1] / cell_vol[cell_id];
    double cocg23 = corr_grad_lin[cell_id][2][1] / cell_vol[cell_id];
    double cocg31 = corr_grad_lin[cell_id][0][2] / cell_vol[cell_id];
    double cocg32 = corr_grad_lin[cell_id][1][2] / cell_vol[cell_id];
    double cocg33 = corr_grad_lin[cell_id][2][2] / cell_vol[cell_id];

    double a11 = cocg22 * cocg33 - cocg32 * cocg23;
    double a12 = cocg32 * cocg13 - cocg12 * cocg33;
    double a13 = cocg12 * cocg23 - cocg22 * cocg13;
    double a21 = cocg31 * cocg23 - cocg21 * cocg33;
    double a22 = cocg11 * cocg33 - cocg31 * cocg13;
    double a23 = cocg21 * cocg13 - cocg11 * cocg23;
    double a31 = cocg21 * cocg32 - cocg31 * cocg22;
    double a32 = cocg31 * cocg12 - cocg11 * cocg32;
    double a33 = cocg11 * cocg22 - cocg21 * cocg12;

    double det_inv = cocg11 * a11 + cocg21 * a12 + cocg31 * a13;

    if (fabs(det_inv) >= 1.e-15) {
      det_inv = 1. / det_inv;

      corr_grad_lin[cell_id][0][0] = a11 * det_inv;
      corr_grad_lin[cell_id][0][1] = a12 * det_inv;
      corr_grad_lin[cell_id][0][2] = a13 * det_inv;
      corr_grad_lin[cell_id][1][0] = a21 * det_inv;
      corr_grad_lin[cell_id][1][1] = a22 * det_inv;
      corr_grad_lin[cell_id][1][2] = a23 * det_inv;
      corr_grad_lin[cell_id][2][0] = a31 * det_inv;
      corr_grad_lin[cell_id][2][1] = a32 * det_inv;
      corr_grad_lin[cell_id][2][2] = a33 * det_inv;

      double a1 = corr_grad_lin[cell_id][0][0];
      double a2 = corr_grad_lin[cell_id][0][1];
      double a3 = corr_grad_lin[cell_id][0][2];
      double a4 = corr_grad_lin[cell_id][1][0];
      double a5 = corr_grad_lin[cell_id][1][1];
      double a6 = corr_grad_lin[cell_id][1][2];
      double a7 = corr_grad_lin[cell_id][2][0];
      double a8 = corr_grad_lin[cell_id][2][1];
      double a9 = corr_grad_lin[cell_id][2][2];

      double determinant =  a1 * (a5*a9 - a8*a6)
                          - a2 * (a4*a9 - a7*a6)
                          + a3 * (a4*a8 - a7*a5);

      corr_grad_lin_det[cell_id] = determinant;
    }
    else {
      corr_grad_lin[cell_id][0][0] = 0.;
      corr_grad_lin[cell_id][0][1] = 0.;
      corr_grad_lin[cell_id][0][2] = 0.;
      corr_grad_lin[cell_id][1][0] = 0.;
      corr_grad_lin[cell_id][1][1] = 0.;
      corr_grad_lin[cell_id][1][2] = 0.;
      corr_grad_lin[cell_id][2][0] = 0.;
      corr_grad_lin[cell_id][2][1] = 0.;
      corr_grad_lin[cell_id][2][2] = 0.;

      corr_grad_lin_det[cell_id] = 1.;
    }
  }

  if (m->halo != NULL) {
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, corr_grad_lin_det);
    cs_halo_sync_var_strided(m->halo, CS_HALO_STANDARD,
                             (cs_real_t *)corr_grad_lin, 9);
    /* TODO handle rotational periodicity */
  }
}

/*----------------------------------------------------------------------------
 * Compute quantities associated to faces (border or internal)
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx_lst    <--  "face -> vertices" connectivity list
 *   face_normal     -->  surface normal of the face
 *
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_compute_face_normal(cs_lnum_t         n_faces,
                     const cs_real_t   vtx_coord[][3],
                     const cs_lnum_t   face_vtx_idx[],
                     const cs_lnum_t   face_vtx[],
                     cs_real_t         face_normal[][3])
{
  /* Checking */

  assert(face_normal != NULL || n_faces == 0);

  /* Loop on faces */

# pragma omp parallel for  if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    if (n_face_vertices == 3) {
      const cs_lnum_t v0 = face_vtx[s_id];
      const cs_lnum_t v1 = face_vtx[s_id+1];
      const cs_lnum_t v2 = face_vtx[s_id+2];
      cs_real_t v01[3], v02[3], vn[3];
      for (cs_lnum_t i = 0; i < 3; i++)
        v01[i] = vtx_coord[v1][i] - vtx_coord[v0][i];
      for (cs_lnum_t i = 0; i < 3; i++)
        v02[i] = vtx_coord[v2][i] - vtx_coord[v0][i];
      cs_math_3_cross_product(v01, v02, vn);
      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = 0.5*vn[i];
    }

    else {

      /* Compute approximate face center coordinates for the polygon */

      cs_real_t a_center[3] = {0, 0, 0};
      cs_real_t f_norm[3] = {0, 0, 0};

      for (cs_lnum_t j = s_id; j < e_id; j++) {
        const cs_lnum_t v0 = face_vtx[j];
        for (cs_lnum_t i = 0; i < 3; i++)
          a_center[i] += vtx_coord[v0][i];
      }

      for (cs_lnum_t i = 0; i < 3; i++)
        a_center[i] /= n_face_vertices;

      /* loop on edges, with implied subdivision into triangles
         defined by edge and cell center */

      cs_real_t vc0[3], vc1[3], vn[3];

      for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

        const cs_lnum_t v0 = face_vtx[s_id + tri_id];
        const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

        for (cs_lnum_t i = 0; i < 3; i++) {
          vc0[i] = vtx_coord[v0][i] - a_center[i];
          vc1[i] = vtx_coord[v1][i] - a_center[i];
        }

        cs_math_3_cross_product(vc0, vc1, vn);

        for (cs_lnum_t i = 0; i < 3; i++)
          f_norm[i] += vn[i];

      }

      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = 0.5*f_norm[i];

    } /* end of test on triangle */

  } /* end of loop on faces */
}

/*----------------------------------------------------------------------------
 * Compute quantities associated to faces (border or internal)
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity
 *   face_cog        -->  coordinates of the center of gravity of the faces
 *   face_normal     -->  face surface normals
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_compute_face_quantities(const cs_lnum_t    n_faces,
                         const cs_real_3_t  vtx_coord[],
                         const cs_lnum_t    face_vtx_idx[],
                         const cs_lnum_t    face_vtx[],
                         cs_real_3_t        face_cog[],
                         cs_real_3_t        face_normal[])
{
  /* Checking */

  assert(face_cog != NULL || n_faces == 0);
  assert(face_normal != NULL || n_faces == 0);

  const cs_real_t one_third = 1./3.;
  const cs_real_t s_epsilon = 1.e-32; /* TODO define better "zero" threshold */

  /* Loop on faces */

# pragma omp parallel for  if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    if (n_face_vertices == 3) {
      const cs_lnum_t v0 = face_vtx[s_id];
      const cs_lnum_t v1 = face_vtx[s_id+1];
      const cs_lnum_t v2 = face_vtx[s_id+2];
      cs_real_t v01[3], v02[3], vn[3];
      for (cs_lnum_t i = 0; i < 3; i++)
        face_cog[f_id][i] = one_third * (  vtx_coord[v0][i]
                                         + vtx_coord[v1][i]
                                         + vtx_coord[v2][i]);
      for (cs_lnum_t i = 0; i < 3; i++)
        v01[i] = vtx_coord[v1][i] - vtx_coord[v0][i];
      for (cs_lnum_t i = 0; i < 3; i++)
        v02[i] = vtx_coord[v2][i] - vtx_coord[v0][i];
      cs_math_3_cross_product(v01, v02, vn);
      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = 0.5*vn[i];
    }

    else { /* For non-triangle faces, assume a division into triangles
              joining edges and an approximate face center */

      /* Compute approximate face center coordinates for the polygon */

      cs_real_t a_center[3] = {0, 0, 0};
      cs_real_t f_center[3] = {0, 0, 0}, f_norm[3] = {0, 0, 0};

      for (cs_lnum_t j = s_id; j < e_id; j++) {
        const cs_lnum_t v0 = face_vtx[j];
        for (cs_lnum_t i = 0; i < 3; i++)
          a_center[i] += vtx_coord[v0][i];
      }

      for (cs_lnum_t i = 0; i < 3; i++)
        a_center[i] /= n_face_vertices;

      cs_real_t sum_w = 0;

      /* In most cases, the following 2 loops can be merged into a single loop,
         but for very bad quality faces, some sub-triangles could be oriented
         differently, so we use 2 passes for safety. */

      if (n_face_vertices < 8) { /* version with local caching for most cases */

        cs_real_t vc0[3], vc1[3], vn[8][3], vtc[8][3];

        /* First pass (face normal) */

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
            vtc[tri_id][i] = vtx_coord[v1][i] + vtx_coord[v0][i] + a_center[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn[tri_id]);

          for (cs_lnum_t i = 0; i < 3; i++)
            f_norm[i] += vn[tri_id][i];

        }

        for (cs_lnum_t i = 0; i < 3; i++)
          f_norm[i] = 0.5*f_norm[i];

        /* Second pass (face center) */

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          cs_real_t w = cs_math_3_norm(vn[tri_id]);

          if (cs_math_3_dot_product(vn[tri_id], f_norm) < 0.0)
            w *= -1.0;

          for (cs_lnum_t i = 0; i < 3; i++)
            f_center[i] += w*vtc[tri_id][i];

          sum_w += w;

        }

      }
      else  { /* generic version */

        cs_real_t vc0[3], vc1[3], vn[3], vtc[3];

        /* First pass (face normal) */

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn);

          for (cs_lnum_t i = 0; i < 3; i++)
            f_norm[i] += vn[i];

        }

        for (cs_lnum_t i = 0; i < 3; i++)
          f_norm[i] = 0.5*f_norm[i];

        /* Second pass (face center) */

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
            vtc[i] = vtx_coord[v1][i] + vtx_coord[v0][i] + a_center[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn);

          cs_real_t w = cs_math_3_norm(vn);

          if (cs_math_3_dot_product(vn, f_norm) < 0.0)
            w *= -1.0;

          for (cs_lnum_t i = 0; i < 3; i++)
            f_center[i] += w*vtc[i];

          sum_w += w;

        }

      }

      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = f_norm[i];

      if (sum_w > s_epsilon) {
        for (cs_lnum_t i = 0; i < 3; i++)
          face_cog[f_id][i] = one_third * f_center[i]/sum_w;
      }
      else {
        for (cs_lnum_t i = 0; i < 3; i++)
          face_cog[f_id][i] = a_center[i];
      }

    } /* end of test on triangle */

  } /* end of loop on faces */
}

/*----------------------------------------------------------------------------
 * Compute face surfaces based on face norms.
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   face_norm       <--  face surface normals
 *   face_surf       -->  face surfaces
 *----------------------------------------------------------------------------*/

static void
_compute_face_surface(cs_lnum_t        n_faces,
                      const cs_real_t  face_norm[],
                      cs_real_t        face_surf[])
{
# pragma omp parallel for  if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
    face_surf[f_id] = cs_math_3_norm(face_norm + f_id*3);
}

/*----------------------------------------------------------------------------
 * Adjust the position of face centers to obtain a contribution to cell
 * volumes consistent with that obtained using a splitting of the face
 * into triangles using face edges and an approxmate center.
 *
 * This adjustment was done by default from Code_Saturne 1.1 to 5.2.
 * To handle warped faces, a correction of the cell center is used so
 * that the contribution to the cell volume using Green's theorem is
 * consistent with the volume computed for sub-faces.
 *
 * It is based on using the coordinates origin instead of the cell centers
 * to compute the correction, whereas the volume computation uses the cell
 * center for better rounding behavior, and uses the cell vertices center
 * instead of the computed face center, so the consistency gain with
 * warped faces (where this is supposed to be useful) would need to be
 * checked.
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity list
 *   face_cog        <->  coordinates of the center of gravity of the faces
 *   face_norm       <--  face surface normals
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_adjust_face_cog_v11_v52(cs_lnum_t         n_faces,
                         const cs_real_t   vtx_coord[][3],
                         const cs_lnum_t   face_vtx_idx[],
                         const cs_lnum_t   face_vtx[],
                         cs_real_t         face_cog[][3],
                         const cs_real_t   face_norm[][3])
{
  const cs_real_t one_third = 1./3.;

  /* Loop on faces
   --------------- */

  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    cs_real_t tri_vol_part = 0.;

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    /* No warping - related correction required for triangles*/
    if (n_face_vertices < 4)
      continue;

    /* Compute approximate face center coordinates for the polygon */

    cs_real_t a_center[3] = {0, 0, 0}, f_center[3] = {0, 0, 0};

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      const cs_lnum_t v0 = face_vtx[j];
      for (cs_lnum_t i = 0; i < 3; i++)
        a_center[i] += vtx_coord[v0][i];
    }

    for (cs_lnum_t i = 0; i < 3; i++)
      a_center[i] /= n_face_vertices;

    /* loop on edges, with implied subdivision into triangles
       defined by edge and cell center */

    cs_real_t vc0[3], vc1[3], vn[3], vtc[3];
    cs_real_t sum_w = 0;

    for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

      const cs_lnum_t v0 = face_vtx[s_id + tri_id];
      const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

      for (cs_lnum_t i = 0; i < 3; i++) {
        vc0[i] = vtx_coord[v0][i] - a_center[i];
        vc1[i] = vtx_coord[v1][i] - a_center[i];
        vtc[i] = vtx_coord[v1][i] + vtx_coord[v0][i] + a_center[i];
      }

      cs_math_3_cross_product(vc0, vc1, vn);

      tri_vol_part += one_third * 0.5 * cs_math_3_dot_product(vtc, vn);

      cs_real_t w = 0.5 * cs_math_3_norm(vn);

      if (cs_math_3_dot_product(vn, face_norm[f_id]) < 0.0)
        w *= -1.0;

      sum_w += w;

      for (cs_lnum_t i = 0; i < 3; i++)
        f_center[i] += w * vtc[i];

    } /* End of loop on triangles of the face */

    /*
     * Compute the center of gravity G(P) of the polygon P, and
     * compute the part of volume of the polygon (before rectification)
     *
     *  -->    ->
     *  OG(P).N(P)
     *------------ */

    cs_real_t face_vol_part = 0.0;

    for (cs_lnum_t i = 0; i < 3; i++) {
      f_center[i] *= one_third / sum_w;
      face_vol_part += (f_center[i] * face_norm[f_id][i]);
    }

    cs_real_t rectif_cog =    (tri_vol_part - face_vol_part)
                            / (sum_w * sum_w);

    for (cs_lnum_t i = 0; i < 3; i++)
      face_cog[f_id][i] = f_center[i] + rectif_cog * face_norm[f_id][i];

  } /* End of loop on faces */
}

/*----------------------------------------------------------------------------
 * Refine face center computation based on initial value.
 *
 * This may be useful for warped faces; for plance faces, the face center
 * should remain identical.
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity
 *   face_cog        <->  coordinates of the center of gravity of the faces
 *   face_norm       <--  face surface normals
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_refine_warped_face_centers(cs_lnum_t          n_faces,
                            const cs_real_3_t  vtx_coord[],
                            const cs_lnum_t    face_vtx_idx[],
                            const cs_lnum_t    face_vtx[],
                            cs_real_t         face_cog[][3],
                            const cs_real_t   face_norm[][3])
{
  /* Checking */

  assert(face_cog != NULL || n_faces == 0);

  const cs_real_t one_third = 1./3.;
  const cs_real_t s_epsilon = 1.e-32; /* TODO define better "zero" threshold */

  /* Loop on faces */

  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    if (n_face_vertices > 3) {

      cs_real_t ref_size = sqrt(cs_math_3_norm(face_norm[f_id]));

      for (int ite = 0; ite < 5; ite++) {

        /* Use previous face center coordinates for the polygon */

        cs_real_t a_center[3];
        cs_real_t f_center[3] = {0, 0, 0};

        for (cs_lnum_t i = 0; i < 3; i++)
          a_center[i] = face_cog[f_id][i];

        /* loop on edges, with implied subdivision into triangles
           defined by edge and cell center */

        cs_real_t vc0[3], vc1[3], vn[3], vtc[3];
        cs_real_t sum_w = 0;

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
            vtc[i] = vtx_coord[v1][i] + vtx_coord[v0][i] + a_center[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn);

          cs_real_t w = cs_math_3_norm(vn);

          if (cs_math_3_dot_product(vn, face_norm[f_id]) < 0.0)
            w *= -1.0;

          sum_w += w;

          for (cs_lnum_t i = 0; i < 3; i++)
            f_center[i] += w*vtc[i];

        }

        if (sum_w > s_epsilon) {
          for (cs_lnum_t i = 0; i < 3; i++)
            face_cog[f_id][i] = one_third * f_center[i]/sum_w;
          if (  cs_math_3_distance(face_cog[f_id], a_center)/ref_size < 1e-8)
            break;
        }
        else
          break;

      } /* end of face iterations */

    } /* end of test on triangle */

  } /* end of loop on faces */
}

/*----------------------------------------------------------------------------
 * Recompute quantities associated to faces (border or internal) when
 * the quality of the mesh is not good enough
 *----------------------------------------------------------------------------*/

static void
_correct_cell_face_center(const cs_mesh_t  *mesh,
                          const cs_lnum_t   n_cells_with_ghosts,
                          const cs_lnum_t   n_i_faces,
                          const cs_lnum_t   n_b_faces,
                          const cs_lnum_2_t i_face_cells[],
                          const cs_lnum_t   b_face_cells[],
                          cs_real_3_t       cell_cen[],
                          cs_real_3_t       i_face_cog[],
                          cs_real_3_t       b_face_cog[],
                          cs_real_3_t       i_face_normal[],
                          cs_real_3_t       b_face_normal[])
{
  int nitmax = 500;
  cs_real_3_t *i_face_cog0, *b_face_cog0;
  cs_real_3_t *i_face_cen, *b_face_cen;

  cs_real_t *relaxf;
  cs_real_t *relaxb;

  cs_real_33_t *dxidxj;
  cs_real_t *determinant;

  BFT_MALLOC(i_face_cog0, n_i_faces, cs_real_3_t);
  BFT_MALLOC(b_face_cog0, n_b_faces, cs_real_3_t);
  BFT_MALLOC(i_face_cen, n_i_faces, cs_real_3_t);
  BFT_MALLOC(b_face_cen, n_b_faces, cs_real_3_t);
  BFT_MALLOC(relaxf, n_i_faces, cs_real_t);
  BFT_MALLOC(relaxb, n_b_faces, cs_real_t);
  BFT_MALLOC(dxidxj, n_cells_with_ghosts, cs_real_33_t);
  BFT_MALLOC(determinant, n_cells_with_ghosts, cs_real_t);

  /* Iterative process */
  for (int sweep = 0; sweep < 1; sweep++) {

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      for (int i = 0; i < 3; i++)
        i_face_cog0[face_id][i] = i_face_cog[face_id][i];

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      for (int i = 0; i < 3; i++)
        b_face_cog0[face_id][i] = b_face_cog[face_id][i];

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      cs_lnum_t cell_id2 = i_face_cells[face_id][1];

      /* IF0 . S */
      double ps1 = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                  i_face_cog0[face_id],
                                                  i_face_normal[face_id]);

      /* IJ . S */
      double ps2 = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                  cell_cen[cell_id2],
                                                  i_face_normal[face_id]);

      double lambda = 0.5;
      if (CS_ABS(ps2) > 1.e-20)
        lambda = ps1 / ps2;

      lambda = CS_MAX(lambda, 1./3.);
      lambda = CS_MIN(lambda, 2./3.);

      /* F = I + lambda * IJ */
      for (int i = 0; i < 3; i++)
        i_face_cen[face_id][i] = cell_cen[cell_id1][i]
                               + lambda * (   cell_cen[cell_id2][i]
                                            - cell_cen[cell_id1][i]);
    }

    /* Compute the projection of I on the boundary face: b_face_cen */
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_3_t normal;
      /* Normal is vector 0 if the b_face_normal norm is too small */
      cs_math_3_normalise(b_face_normal[face_id], normal);

      if (cs_math_3_norm(normal) > 0.) {
        double lambda = cs_math_3_distance_dot_product(cell_cen[cell_id],
                                                       b_face_cog[face_id],
                                                       normal);

        for (int i = 0; i < 3; i++)
          b_face_cen[face_id][i] = cell_cen[cell_id][i] + lambda * normal[i];
      } else {
        for (int i = 0; i < 3; i++)
          b_face_cen[face_id][i] =  b_face_cog[face_id][i];
      }
    }

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      relaxf[face_id] = 1.;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      relaxb[face_id] = 1.;

    int iiter = 0;
    cs_gnum_t irelax = 0;

    do
    {
      iiter +=1;

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        for (int i = 0; i < 3; i++)
          i_face_cog[face_id][i] = (1. - relaxf[face_id]) * i_face_cog0[face_id][i]
                                       + relaxf[face_id]  * i_face_cen[face_id][i];
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        for (int i = 0; i < 3; i++)
          b_face_cog[face_id][i] = (1. - relaxb[face_id]) * b_face_cog0[face_id][i]
                                       + relaxb[face_id]  * b_face_cen[face_id][i];
      }

      for (cs_lnum_t cell_id = 0; cell_id <  n_cells_with_ghosts; cell_id++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            dxidxj[cell_id][i][j] = 0.;

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            double fluxij = i_face_cog[face_id][i]//TODO minus celli
              * i_face_normal[face_id][j];
            dxidxj[cell_id1][i][j] += fluxij;
            dxidxj[cell_id2][i][j] -= fluxij;
          }
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            double fluxij =  b_face_cog[face_id][i]
                           * b_face_normal[face_id][j];
            dxidxj[cell_id][i][j] += fluxij;
          }
      }

      for (cs_lnum_t cell_id = 0; cell_id <  mesh->n_cells; cell_id++) {
        double vol = (   dxidxj[cell_id][0][0]
                       + dxidxj[cell_id][1][1]
                       + dxidxj[cell_id][2][2] ) / 3.;

        //FIXME
        if (vol >= 0)
          vol = CS_MAX(vol,1.e-20);

        double a1 = dxidxj[cell_id][0][0] / vol;
        double a2 = dxidxj[cell_id][0][1] / vol;
        double a3 = dxidxj[cell_id][0][2] / vol;
        double a4 = dxidxj[cell_id][1][0] / vol;
        double a5 = dxidxj[cell_id][1][1] / vol;
        double a6 = dxidxj[cell_id][1][2] / vol;
        double a7 = dxidxj[cell_id][2][0] / vol;
        double a8 = dxidxj[cell_id][2][1] / vol;
        double a9 = dxidxj[cell_id][2][2] / vol;

        determinant[cell_id] = fabs(  a1 * (a5*a9 - a8*a6)
                                    - a2 * (a4*a9 - a7*a6)
                                    + a3 * (a4*a8 - a7*a5) );

        //FIXME
        determinant[cell_id] = CS_MAX(determinant[cell_id],1.e-20);
        determinant[cell_id] = CS_MIN(determinant[cell_id],1./determinant[cell_id]);

        /* M-matrix structure control */
        double mmatrice1a = a1 - CS_ABS(a2) - CS_ABS(a3);
        double mmatrice1b = a1 - CS_ABS(a4) - CS_ABS(a7);
        double mmatrice2a = a5 - CS_ABS(a4) - CS_ABS(a6);
        double mmatrice2b = a5 - CS_ABS(a2) - CS_ABS(a8);
        double mmatrice3a = a9 - CS_ABS(a7) - CS_ABS(a8);
        double mmatrice3b = a9 - CS_ABS(a6) - CS_ABS(a3);

        if (   mmatrice1a <= 0. || mmatrice1b <= 0.
            || mmatrice2a <= 0. || mmatrice2b <= 0.
            || mmatrice3a <= 0. || mmatrice3b <= 0.)
          determinant[cell_id] = 0.;
      }

      if (mesh->halo != NULL)
        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, determinant);

      irelax = 0;

      //FIXME test was 0.001
      cs_real_t threshold = 0.1;
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        if (determinant[cell_id1] < threshold || determinant[cell_id2] < threshold) {
          irelax +=1;
          relaxf[face_id] *= 0.95;
        }
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        if (determinant[cell_id] < threshold) {
          irelax += 1;
          relaxb[face_id] *= 0.95;
        }
      }

      cs_parall_counter(&irelax, 1);

    } while (iiter < nitmax && irelax > 0);

  }
  BFT_FREE(i_face_cog0);
  BFT_FREE(b_face_cog0);
  BFT_FREE(i_face_cen);
  BFT_FREE(b_face_cen);
  BFT_FREE(relaxf);
  BFT_FREE(relaxb);
  BFT_FREE(dxidxj);
  BFT_FREE(determinant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell centers and volumes.
 *
 * \param[in]   mesh         pointer to mesh structure
 * \param[in]   i_face_norm  surface normal of internal faces
 * \param[in]   i_face_cog   center of gravity of internal faces
 * \param[in]   b_face_norm  surface normal of border faces
 * \param[in]   b_face_cog   center of gravity of border faces
 * \param[out]  cell_cen     cell centers
 * \param[out]  cell_vol     cell volumes
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_quantities(const cs_mesh_t    *mesh,
                         const cs_real_3_t   i_face_norm[],
                         const cs_real_3_t   i_face_cog[],
                         const cs_real_3_t   b_face_norm[],
                         const cs_real_3_t   b_face_cog[],
                         cs_real_3_t         cell_cen[restrict],
                         cs_real_t           cell_vol[restrict])
{
  /* Mesh connectivity */

  const  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const  cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const  cs_lnum_t  n_cells = mesh->n_cells;
  const  cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const  cs_lnum_2_t  *i_face_cells
    = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const  cs_lnum_t  *b_face_cells = mesh->b_face_cells;

  /* Checking */

  assert(cell_cen != NULL);
  assert(cell_vol != NULL);

  /* Compute approximate cell center using face centers */

  cs_real_3_t *a_cell_cen;
  BFT_MALLOC(a_cell_cen, n_cells_ext, cs_real_3_t);

  cs_mesh_quantities_cell_faces_cog(mesh,
                                    (const cs_real_t *)i_face_norm,
                                    (const cs_real_t *)i_face_cog,
                                    (const cs_real_t *)b_face_norm,
                                    (const cs_real_t *)b_face_cog,
                                    (cs_real_t *)a_cell_cen);

  /* Initialization */

  for (cs_lnum_t j = 0; j < n_cells_ext; j++)
    cell_vol[j] = 0.;

  for (cs_lnum_t j = 0; j < n_cells_ext; j++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[j][i] = 0.;
  }

  /* Loop on interior faces
     ---------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    /* For each cell sharing the internal face, we update
     * cell_cen and cell_area */

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    /* Implicit subdivision of cell into face vertices-cell-center pyramids */

    if (c_id1 > -1) {
      cs_real_t pyra_vol_3 = cs_math_3_distance_dot_product(a_cell_cen[c_id1],
                                                            i_face_cog[f_id],
                                                            i_face_norm[f_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[c_id1][i] += pyra_vol_3 *(  0.75*i_face_cog[f_id][i]
                                           + 0.25*a_cell_cen[c_id1][i]);
      cell_vol[c_id1] += pyra_vol_3;
    }
    if (c_id2 > -1) {
      cs_real_t pyra_vol_3 = cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                            a_cell_cen[c_id2],
                                                            i_face_norm[f_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[c_id2][i] += pyra_vol_3 *(  0.75*i_face_cog[f_id][i]
                                           + 0.25*a_cell_cen[c_id2][i]);
      cell_vol[c_id2] += pyra_vol_3;
    }

  } /* End of loop on interior faces */

  /* Loop on boundary faces
     --------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* For each cell sharing a border face, we update the numerator
     * of cell_cen and cell_area */

    cs_lnum_t c_id1 = b_face_cells[f_id];

    /* Computation of the area of the face
       (note that c_id1 == -1 may happen for isolated faces,
       which are cleaned afterwards) */

    if (c_id1 > -1) {
      cs_real_t pyra_vol_3 = cs_math_3_distance_dot_product(a_cell_cen[c_id1],
                                                            b_face_cog[f_id],
                                                            b_face_norm[f_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[c_id1][i] += pyra_vol_3 *(  0.75*b_face_cog[f_id][i]
                                           + 0.25*a_cell_cen[c_id1][i]);
      cell_vol[c_id1] += pyra_vol_3;
    }

  } /* End of loop on boundary faces */

  BFT_FREE(a_cell_cen);

  /* Loop on cells to finalize the computation
     ----------------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[c_id][i] /= cell_vol[c_id];

    cell_vol[c_id] /= 3.0;

  }
}

/*----------------------------------------------------------------------------*
 * Compute new cell centers by minimizing the distance to faces
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   i_face_normal  <--  surface normal of internal faces
 *   i_face_cog     <--  center of gravity of internal faces
 *   b_face_normal  <--  surface normal of border faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       -->  center of gravity of cells
 *----------------------------------------------------------------------------*/

static void
_recompute_cell_cen_face(const cs_mesh_t     *mesh,
                         const cs_real_3_t   i_face_normal[],
                         const cs_real_3_t   i_face_cog[],
                         const cs_real_3_t   b_face_normal[],
                         const cs_real_3_t   b_face_cog[],
                         cs_real_3_t         cell_cen[])
{
  const  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const  cs_lnum_t  n_b_faces = mesh->n_b_faces;

  const  cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  const  cs_lnum_2_t  *i_face_cells
    = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const  cs_lnum_t  *b_face_cells = mesh->b_face_cells;

  /* First pass of verification */
  int *pb1;
  BFT_MALLOC(pb1, n_cells_with_ghosts, int);

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
    pb1[cell_id] = 0;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* IF . S */
    double psi1 = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                 i_face_cog[face_id],
                                                 i_face_normal[face_id]);
    /* JF . S */
    double psj1 = cs_math_3_distance_dot_product(cell_cen[cell_id2],
                                                 i_face_cog[face_id],
                                                 i_face_normal[face_id]);
    if (psi1 < 0.)
      pb1[cell_id1]++;
    if (psj1 > 0.)
      pb1[cell_id2]++;
  }

  cs_gnum_t cpt1 = 0;
  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
    if (pb1[cell_id] > 0) cpt1++;
  cs_parall_counter(&cpt1, 1);

  if (cpt1 > 0) {
    bft_printf("Total number of cell centers on the other side of a face "
               "(before correction) = %lu / %d\n", cpt1, mesh->n_cells);

    /* Second pass */
    cs_real_33_t *a;
    cs_real_3_t  *b;
    cs_real_3_t  *cdgbis;

    BFT_MALLOC(a, n_cells_with_ghosts, cs_real_33_t);
    BFT_MALLOC(b, n_cells_with_ghosts, cs_real_3_t);
    BFT_MALLOC(cdgbis, n_cells_with_ghosts, cs_real_3_t);

    /* init matrice et second membre */
    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++) {
      for (int i = 0; i < 3; i++) {
        b[cell_id][i] = 0.;
        for (int j = 0; j < 3; j++)
          a[cell_id][i][j] = 0.;
      }
    }

    /* Contribution from interior faces */
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      cs_lnum_t cell_id2 = i_face_cells[face_id][1];
      double surfn = cs_math_3_norm(i_face_normal[face_id]);

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          a[cell_id1][i][j] +=   i_face_normal[face_id][i]
                               * i_face_normal[face_id][j] / surfn;
          a[cell_id2][i][j] +=   i_face_normal[face_id][i]
                               * i_face_normal[face_id][j] / surfn;
        }

      double ps = cs_math_3_dot_product(i_face_normal[face_id],
                                        i_face_cog[face_id]);

      for (int i = 0; i < 3; i++) {
        b[cell_id1][i] += ps * i_face_normal[face_id][i] / surfn;
        b[cell_id2][i] += ps * i_face_normal[face_id][i] / surfn;
      }

    }

    /* Contribution from boundary faces */
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      cs_lnum_t cell_id = b_face_cells[face_id];
      double surfn = cs_math_3_norm(b_face_normal[face_id]);

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          a[cell_id][i][j] +=   b_face_normal[face_id][i]
                              * b_face_normal[face_id][j] / surfn;
        }

      double ps = cs_math_3_dot_product(b_face_normal[face_id],
                                        b_face_cog[face_id]);

      for (int i = 0; i < 3; i++) {
        b[cell_id][i] += ps * b_face_normal[face_id][i] / surfn;
      }

    }

    /* invert system */
    double aainv[3][3];
    double bb[3];
    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {

      cdgbis[cell_id][0] = cell_cen[cell_id][0];
      cdgbis[cell_id][1] = cell_cen[cell_id][1];
      cdgbis[cell_id][2] = cell_cen[cell_id][2];

      if (pb1[cell_id] > 0) {

        double adim = a[cell_id][0][0] + a[cell_id][1][1] + a[cell_id][2][2];

        if (adim > 0.) {
          bb[0] = b[cell_id][0] / adim;
          bb[1] = b[cell_id][1] / adim;
          bb[2] = b[cell_id][2] / adim;

          /* Matrix inversion */
          double cocg11 = a[cell_id][0][0] / adim;
          double cocg12 = a[cell_id][0][1] / adim;
          double cocg13 = a[cell_id][0][2] / adim;
          double cocg21 = a[cell_id][1][0] / adim;
          double cocg22 = a[cell_id][1][1] / adim;
          double cocg23 = a[cell_id][1][2] / adim;
          double cocg31 = a[cell_id][2][0] / adim;
          double cocg32 = a[cell_id][2][1] / adim;
          double cocg33 = a[cell_id][2][2] / adim;

          double a11 = cocg22 * cocg33 - cocg32 * cocg23;
          double a12 = cocg32 * cocg13 - cocg12 * cocg33;
          double a13 = cocg12 * cocg23 - cocg22 * cocg13;
          double a21 = cocg31 * cocg23 - cocg21 * cocg33;
          double a22 = cocg11 * cocg33 - cocg31 * cocg13;
          double a23 = cocg21 * cocg13 - cocg11 * cocg23;
          double a31 = cocg21 * cocg32 - cocg31 * cocg22;
          double a32 = cocg31 * cocg12 - cocg11 * cocg32;
          double a33 = cocg11 * cocg22 - cocg21 * cocg12;

          double det_inv = cocg11 * a11 + cocg21 * a12 + cocg31 * a13;

          if (CS_ABS(det_inv) >= 1.e-15) {
            det_inv = 1. / det_inv;

            aainv[0][0] = a11 * det_inv;
            aainv[0][1] = a12 * det_inv;
            aainv[0][2] = a13 * det_inv;
            aainv[1][0] = a21 * det_inv;
            aainv[1][1] = a22 * det_inv;
            aainv[1][2] = a23 * det_inv;
            aainv[2][0] = a31 * det_inv;
            aainv[2][1] = a32 * det_inv;
            aainv[2][2] = a33 * det_inv;

            for (int i = 0; i < 3; i++)
              cdgbis[cell_id][i] =   aainv[i][0] * bb[0]
                                   + aainv[i][1] * bb[1]
                                   + aainv[i][2] * bb[2];
          }
        }
      }
    }

    if (mesh->halo != NULL) {
      cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                               (cs_real_t *)cdgbis, 3);
      if (mesh->n_init_perio > 0)
        cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                                  (cs_real_t *)cdgbis);
    }

    /* Second verification */

    int *pb2;
    BFT_MALLOC(pb2, n_cells_with_ghosts, int);

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
      pb2[cell_id] = 0;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      cs_lnum_t cell_id2 = i_face_cells[face_id][1];

      /* IF . S */
      double psi1 = cs_math_3_distance_dot_product(cdgbis[cell_id1],
                                                   i_face_cog[face_id],
                                                   i_face_normal[face_id]);
      /* JF . S */
      double psj1 = cs_math_3_distance_dot_product(cdgbis[cell_id2],
                                                   i_face_cog[face_id],
                                                   i_face_normal[face_id]);
      if (psi1 < 0.)
        pb2[cell_id1]++;
      if (psj1 > 0.)
        pb2[cell_id2]++;
    }

    cs_gnum_t cpt2 = 0;
    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
      if (pb2[cell_id] > 0)
        cpt2++;
    cs_parall_counter(&cpt2, 1);

    bft_printf("Total number of cell centers on the other side of a face "
               "(after correction) = %lu / %d\n", cpt2, mesh->n_cells);

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
      if (pb1[cell_id] > 0 && pb2[cell_id] == 0) {

          cell_cen[cell_id][0] = cdgbis[cell_id][0];
          cell_cen[cell_id][1] = cdgbis[cell_id][1];
          cell_cen[cell_id][2] = cdgbis[cell_id][2];
      }
    }

    if (mesh->halo != NULL) {
      cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                               (cs_real_t *)cell_cen, 3);
      if (mesh->n_init_perio > 0)
        cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                                  (cs_real_t *)cell_cen);
    }

    BFT_FREE(a);
    BFT_FREE(b);
    BFT_FREE(cdgbis);
    BFT_FREE(pb2);
  }

  /* Free memory */
  BFT_FREE(pb1);
}

/*----------------------------------------------------------------------------
 * Compute the volume of cells C from their n faces F(i) and their center of
 * gravity G(Fi) where i=0, n-1
 *
 *         1    n-1
 *  G(C) = - .  Sum  Surf(Fi) G(Fi)
 *         3    i=0
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   i_face_norm    <--  surface normal of internal faces
 *   i_face_cog     <--  center of gravity of internal faces
 *   b_face_norm    <--  surface normal of border faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       <--  center of gravity of cells
 *   cell_vol       -->  cells volume
 *----------------------------------------------------------------------------*/

static void
_compute_cell_volume(const cs_mesh_t   *mesh,
                     const cs_real_3_t  i_face_norm[],
                     const cs_real_3_t  i_face_cog[],
                     const cs_real_3_t  b_face_norm[],
                     const cs_real_3_t  b_face_cog[],
                     const cs_real_3_t  cell_cen[],
                     cs_real_t          cell_vol[])
{
  const cs_real_t  a_third = 1.0/3.0;

  /* Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
    cell_vol[cell_id] = 0;

  /* Loop on internal faces */

  for (cs_lnum_t fac_id = 0; fac_id < mesh->n_i_faces; fac_id++) {

    cs_lnum_t cell_id1 = mesh->i_face_cells[fac_id][0];
    cs_lnum_t cell_id2 = mesh->i_face_cells[fac_id][1];

    cell_vol[cell_id1] += cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                         i_face_cog[fac_id],
                                                         i_face_norm[fac_id]);
    cell_vol[cell_id2] -= cs_math_3_distance_dot_product(cell_cen[cell_id2],
                                                         i_face_cog[fac_id],
                                                         i_face_norm[fac_id]);
  }

  /* Loop on border faces */

  for (cs_lnum_t fac_id = 0; fac_id < mesh->n_b_faces; fac_id++) {

    cs_lnum_t cell_id1 = mesh->b_face_cells[fac_id];

    cell_vol[cell_id1] += cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                         b_face_cog[fac_id],
                                                         b_face_norm[fac_id]);
  }

  /* First Computation of the volume */

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
    cell_vol[cell_id] *= a_third;
}

/*----------------------------------------------------------------------------
 * Correction of small or negative volumes.
 *
 * Based on smoothing with neighbors; does not conserve the total volume.
 *----------------------------------------------------------------------------*/

static void
_cell_bad_volume_correction(const cs_mesh_t   *mesh,
                            cs_real_t          cell_vol[])
{
  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, cell_vol);

  /* Iterations in order to get vol_I / max(vol_J) > critmin */

  double *vol_neib_max;
  BFT_MALLOC(vol_neib_max, mesh->n_cells_with_ghosts, double);

  for (int iter = 0; iter < 10; iter++) {

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
      vol_neib_max[cell_id] = 0.;

    for (cs_lnum_t fac_id = 0; fac_id < mesh->n_i_faces; fac_id++) {

      cs_lnum_t cell_id1 = mesh->i_face_cells[fac_id][0];
      cs_lnum_t cell_id2 = mesh->i_face_cells[fac_id][1];
      double vol1 = cell_vol[cell_id1];
      double vol2 = cell_vol[cell_id2];

      if (vol2 > 0.)
        vol_neib_max[cell_id1] = CS_MAX(vol_neib_max[cell_id1], vol2);

      if (vol1 > 0.)
        vol_neib_max[cell_id2] = CS_MAX(vol_neib_max[cell_id2], vol1);
    }

    /* Previous value of 0.2 sometimes leads to computation divergence */
    /* 0.01 seems better and safer for the moment */
    double critmin = 0.01;

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
      cell_vol[cell_id] = CS_MAX(cell_vol[cell_id],
                                 critmin * vol_neib_max[cell_id]);

    if (mesh->halo != NULL)
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, cell_vol);
  }

  BFT_FREE(vol_neib_max);
}

/*----------------------------------------------------------------------------
 * Compute the total, min, and max volumes of cells.
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   cell_vol       <--  cells volume
 *   min_vol        -->  minimum control volume
 *   max_vol        -->  maximum control volume
 *   tot_vol        -->  total   control volume
 *----------------------------------------------------------------------------*/

static void
_cell_volume_reductions(const cs_mesh_t  *mesh,
                        cs_real_t         cell_vol[],
                        cs_real_t        *min_vol,
                        cs_real_t        *max_vol,
                        cs_real_t        *tot_vol)
{
  *min_vol =  cs_math_infinite_r;
  *max_vol = -cs_math_infinite_r;
  *tot_vol = 0.;

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
    *min_vol = CS_MIN(*min_vol, cell_vol[cell_id]);
    *max_vol = CS_MAX(*max_vol, cell_vol[cell_id]);
    *tot_vol = *tot_vol + cell_vol[cell_id];
  }
}

/*----------------------------------------------------------------------------
 * Compute some distances relative to faces and associated weighting.
 *
 * parameters:
 *   n_i_faces      <--  number of interior faces
 *   n_b_faces      <--  number of border  faces
 *   i_face_cells   <--  interior "faces -> cells" connectivity
 *   b_face_cells   <--  border "faces -> cells" connectivity
 *   i_face_norm    <--  surface normal of interior faces
 *   b_face_norm    <--  surface normal of border faces
 *   i_face_cog     <--  center of gravity of interior faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       <--  cell center
 *   i_dist         -->  distance IJ.Nij for interior faces
 *   b_dist         -->  likewise for border faces
 *   weight         -->  weighting factor (Aij=pond Ai+(1-pond)Aj)
 *----------------------------------------------------------------------------*/

static void
_compute_face_distances(cs_lnum_t          n_i_faces,
                        cs_lnum_t          n_b_faces,
                        const cs_lnum_2_t  i_face_cells[],
                        const cs_lnum_t    b_face_cells[],
                        const cs_real_t    i_face_normal[][3],
                        const cs_real_t    b_face_normal[][3],
                        const cs_real_t    i_face_cog[][3],
                        const cs_real_t    b_face_cog[][3],
                        const cs_real_t    cell_cen[][3],
                        cs_real_t          i_dist[],
                        cs_real_t          b_dist[],
                        cs_real_t          weight[])
{
  cs_gnum_t w_count = 0;

  /* Interior faces */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    const cs_real_t *face_nomal = i_face_normal[face_id];
    cs_real_t normal[3];
    cs_math_3_normalise(face_nomal, normal);

    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* Distance between the neighbor cell centers
     * and dot-product with the normal */
    i_dist[face_id] = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                     cell_cen[cell_id2],
                                                     normal);

    if (CS_ABS(i_dist[face_id]) > 1e-20) {
      /* Distance between the face center of gravity
         and the neighbor cell center
         and dot-product with the normal */
      cs_real_t dist2f = cs_math_3_distance_dot_product(i_face_cog[face_id],
                                                        cell_cen[cell_id2],
                                                        normal);
      weight[face_id] = dist2f / i_dist[face_id];
    }
    else {
      weight[face_id] = 0.5;
    }
    double distmax = cs_math_3_distance(cell_cen[cell_id1],
                                        cell_cen[cell_id2]);

    /* Clipping of cell cell distances */
    if (cs_glob_mesh_quantities_flag & CS_FACE_DISTANCE_CLIP) {
      if (i_dist[face_id] < 0.2 * distmax) {
        w_count++;
        i_dist[face_id] = CS_MAX(i_dist[face_id], 0.2 * distmax);
      }

      /* Clipping of weighting */
      weight[face_id] = CS_MAX(weight[face_id], 0.001);
      weight[face_id] = CS_MIN(weight[face_id], 0.999);
    }
  }

  cs_parall_counter(&w_count, 1);

  if (w_count > 0)
    bft_printf(_("\n"
                 "%llu faces have a too small distance between centers.\n"
                 "For these faces, the weight may be clipped.\n"),
               (unsigned long long)w_count);

  /* Border faces */

  w_count = 0;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    const cs_real_t *face_nomal = b_face_normal[face_id];
    cs_real_t normal[3];
    cs_math_3_normalise(face_nomal, normal);

    cs_lnum_t cell_id = b_face_cells[face_id];

    /* Distance between the face center of gravity
       and the neighbor cell center */
    b_dist[face_id] = cs_math_3_distance_dot_product(cell_cen[cell_id],
                                                     b_face_cog[face_id],
                                                     normal);
    /* Clipping of cell boundary distances */
    if (cs_glob_mesh_quantities_flag & CS_FACE_DISTANCE_CLIP) {
      cs_real_t distmax = cs_math_3_distance(cell_cen[cell_id],
                                             b_face_cog[face_id]);
      if (b_dist[face_id] < 0.2 * distmax) {
        w_count++;
        b_dist[face_id] = CS_MAX(b_dist[face_id], 0.2 * distmax);
      }

    }
  }

  cs_parall_counter(&w_count, 1);

  if (w_count > 0)
    bft_printf(_("\n"
                 "%llu boundary faces have a too small distance between\n"
                 "cell centre and face centre.\n"),
               (unsigned long long)w_count);
}

/*----------------------------------------------------------------------------
 * Compute some vectors to handle non-orthogonalities.
 *
 * Let a face and I, J the centers of neighboring cells
 *   (only I is defined for a border face)
 *
 * The face is oriented from I to J, with Nij its normal.
 *   (border faces are oriented towards the exterior)
 * The norm of Nij is 1.
 * The face surface is Sij.
 *
 * I' and J' are defined as the orthogonal projection of I and J on the line
 * orthogonal to the face passing through the center of gravity F of the face.
 *   (only I' is defined for a border face)
 *
 * We compute here the vector I'J' for interior faces (dijpf)
 *                 the vector II'  for border faces   (diipb)
 *                 the vector OF   for interior faces (dofij)
 *
 * We also have the following formulae
 *   II' = IG - (IG.Nij)Nij
 *   JJ' = JG - (JG.Nij)Nij
 *
 * parameters:
 *   dim            <--  dimension
 *   n_i_faces      <--  number of interior faces
 *   n_b_faces      <--  number of border  faces
 *   i_face_cells   <--  interior "faces -> cells" connectivity
 *   b_face_cells   <--  border "faces -> cells" connectivity
 *   i_face_norm    <--  surface normal of interior faces
 *   b_face_norm    <--  surface normal of border faces
 *   i_face_cog     <--  center of gravity of interior faces
 *   b_face_cog     <--  center of gravity of border faces
 *   i_face_surf    <--  interior faces surface
 *   b_face_surf    <--  border faces surface
 *   cell_cen       <--  cell center
 *   weight         <--  weighting factor (Aij=pond Ai+(1-pond)Aj)
 *   dijpf          -->  vector i'j' for interior faces
 *   diipb          -->  vector ii'  for border faces
 *   dofij          -->  vector OF   for interior faces
 *----------------------------------------------------------------------------*/

static void
_compute_face_vectors(int                dim,
                      const cs_lnum_t    n_i_faces,
                      const cs_lnum_t    n_b_faces,
                      const cs_lnum_2_t  i_face_cells[],
                      const cs_lnum_t    b_face_cells[],
                      const cs_real_t    i_face_normal[],
                      const cs_real_t    b_face_normal[],
                      const cs_real_t    i_face_cog[],
                      const cs_real_t    b_face_cog[],
                      const cs_real_t    i_face_surf[],
                      const cs_real_t    b_face_surf[],
                      const cs_real_t    cell_cen[],
                      const cs_real_t    weight[],
                      const cs_real_t    b_dist[],
                      cs_real_t          dijpf[],
                      cs_real_t          diipb[],
                      cs_real_t          dofij[])
{
  cs_lnum_t face_id;
  cs_lnum_t cell_id;

  cs_real_t dipjp, psi, pond;
  cs_real_t surfnx, surfny, surfnz;
  cs_real_t vecigx, vecigy, vecigz, vecijx, vecijy, vecijz;

  /* Interior faces */

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* Normalized normal */
    surfnx = i_face_normal[face_id*dim]     / i_face_surf[face_id];
    surfny = i_face_normal[face_id*dim + 1] / i_face_surf[face_id];
    surfnz = i_face_normal[face_id*dim + 2] / i_face_surf[face_id];

    /* ---> IJ */
    vecijx = cell_cen[cell_id2*dim]     - cell_cen[cell_id1*dim];
    vecijy = cell_cen[cell_id2*dim + 1] - cell_cen[cell_id1*dim + 1];
    vecijz = cell_cen[cell_id2*dim + 2] - cell_cen[cell_id1*dim + 2];

    /* ---> DIJPP = IJ.NIJ */
    dipjp = vecijx*surfnx + vecijy*surfny + vecijz*surfnz;

    /* ---> DIJPF = (IJ.NIJ).NIJ */
    dijpf[face_id*dim]     = dipjp*surfnx;
    dijpf[face_id*dim + 1] = dipjp*surfny;
    dijpf[face_id*dim + 2] = dipjp*surfnz;

    pond = weight[face_id];

    /* ---> DOFIJ = OF */
    dofij[face_id*dim]     = i_face_cog[face_id*dim]
      - (        pond *cell_cen[cell_id1*dim]
         + (1. - pond)*cell_cen[cell_id2*dim]);

    dofij[face_id*dim + 1] = i_face_cog[face_id*dim + 1]
      - (        pond *cell_cen[cell_id1*dim + 1]
         + (1. - pond)*cell_cen[cell_id2*dim + 1]);

    dofij[face_id*dim + 2] = i_face_cog[face_id*dim + 2]
      - (        pond *cell_cen[cell_id1*dim + 2]
         + (1. - pond)*cell_cen[cell_id2*dim + 2]);
  }

  /* Border faces */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    cell_id = b_face_cells[face_id];

    /* Normalized normal */
    surfnx = b_face_normal[face_id*dim]     / b_face_surf[face_id];
    surfny = b_face_normal[face_id*dim + 1] / b_face_surf[face_id];
    surfnz = b_face_normal[face_id*dim + 2] / b_face_surf[face_id];

    /* ---> IG */
    vecigx = b_face_cog[face_id*dim]     - cell_cen[cell_id*dim];
    vecigy = b_face_cog[face_id*dim + 1] - cell_cen[cell_id*dim + 1];
    vecigz = b_face_cog[face_id*dim + 2] - cell_cen[cell_id*dim + 2];

    /* ---> PSI = IG.NIJ */
    psi = vecigx*surfnx + vecigy*surfny + vecigz*surfnz;

    /* ---> DIIPB = IG - (IG.NIJ)NIJ */
    diipb[face_id*dim]     = vecigx - psi*surfnx;
    diipb[face_id*dim + 1] = vecigy - psi*surfny;
    diipb[face_id*dim + 2] = vecigz - psi*surfnz;

    /* Limiter on boundary face reconstruction */
    if (cs_glob_mesh_quantities_flag & CS_FACE_RECONSTRUCTION_CLIP) {
      double iip = sqrt(   diipb[face_id*dim]     * diipb[face_id*dim]
                         + diipb[face_id*dim + 1] * diipb[face_id*dim + 1]
                         + diipb[face_id*dim + 2] * diipb[face_id*dim + 2]);

      cs_real_t corri = 1.;
      if (iip > 0.5 * b_dist[face_id])
        corri = 0.5 * b_dist[face_id] / iip;

      diipb[face_id*dim]    *= corri;
      diipb[face_id*dim +1] *= corri;
      diipb[face_id*dim +2] *= corri;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute some vectors to handle non-orthogonalities.
 *
 * Let a face and I, J the centers of neighboring cells
 *   (only I is defined for a border face)
 *
 * The face is oriented from I to J, with Nij its normal.
 *   (border faces are oriented towards the exterior)
 * The norm of Nij is 1.
 * The face surface is Sij.
 *
 * I' and J' are defined as the orthogonal projection of I and J on the line
 * orthogonal to the face passing through the center of gravity F of the face.
 *   (only I' is defined for a border face)
 *
 * We compute here the vector II' for interior faces (diipf)
 *                 the vector JJ' for interior faces (djjpf)
 *
 * We also have the following formulae
 *   II' = IF - (IF.Nij)Nij
 *   JJ' = JF - (JF.Nij)Nij
 *
 * parameters:
 *   dim            <--  dimension
 *   n_i_faces      <--  number of interior faces
 *   i_face_cells   <--  interior "faces -> cells" connectivity
 *   i_face_norm    <--  surface normal of interior faces
 *   i_face_cog     <--  center of gravity of interior faces
 *   i_face_surf    <--  interior faces surface
 *   cell_cen       <--  cell center
 *   cell_vol       <--  cell volume
 *   dist           <--  interior distance
 *   diipf          -->  vector ii' for interior faces
 *   djjpf          -->  vector jj' for interior faces
 *----------------------------------------------------------------------------*/

static void
_compute_face_sup_vectors(const cs_lnum_t    n_i_faces,
                          const cs_lnum_2_t  i_face_cells[],
                          const cs_real_t    i_face_normal[][3],
                          const cs_real_t    i_face_cog[][3],
                          const cs_real_t    cell_cen[][3],
                          const cs_real_t    cell_vol[],
                          const cs_real_t    dist[],
                          cs_real_t          diipf[][3],
                          cs_real_t          djjpf[][3])
{

  /* Interior faces */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* Normalized normal */
    cs_real_3_t normal;
    cs_math_3_normalise(i_face_normal[face_id], normal);

    /* ---> IF and JF */
    cs_real_t vec_if[3] = {
      i_face_cog[face_id][0] - cell_cen[cell_id1][0],
      i_face_cog[face_id][1] - cell_cen[cell_id1][1],
      i_face_cog[face_id][2] - cell_cen[cell_id1][2]};

    cs_real_t vec_jf[3] = {
      i_face_cog[face_id][0] - cell_cen[cell_id2][0],
      i_face_cog[face_id][1] - cell_cen[cell_id2][1],
      i_face_cog[face_id][2] - cell_cen[cell_id2][2]};

    /* ---> diipf = IF - (IF.Nij)Nij */
    cs_math_3_orthogonal_projection(normal, vec_if, diipf[face_id]);

    /* ---> djjpf = JF - (JF.Nij)Nij */
    cs_math_3_orthogonal_projection(normal, vec_jf, djjpf[face_id]);

    /* Limiter on interior face reconstruction */
    if (cs_glob_mesh_quantities_flag & CS_FACE_RECONSTRUCTION_CLIP) {

      cs_real_t surfn = cs_math_3_norm(i_face_normal[face_id]);
      cs_real_t iip   = cs_math_3_norm(diipf[face_id]);

      cs_real_t corri = 1.;

      if (0.5 * dist[face_id] < iip)
        corri = 0.5 * dist[face_id] / iip;

      diipf[face_id][0] *= corri;
      diipf[face_id][1] *= corri;
      diipf[face_id][2] *= corri;

      iip = cs_math_3_norm(diipf[face_id]);

      corri = 1.;

      if (0.9 * cell_vol[cell_id1] < surfn * iip)
        corri = 0.9 * cell_vol[cell_id1] / (surfn * iip);

      diipf[face_id][0] *= corri;
      diipf[face_id][1] *= corri;
      diipf[face_id][2] *= corri;

      cs_real_t jjp = cs_math_3_norm(djjpf[face_id]);

      cs_real_t corrj = 1.;

      if (0.5 * dist[face_id] < jjp)
        corrj = 0.5 * dist[face_id] / jjp;

      djjpf[face_id][0] *= corrj;
      djjpf[face_id][1] *= corrj;
      djjpf[face_id][2] *= corrj;

      jjp = cs_math_3_norm(djjpf[face_id]);

      corrj = 1.;
      if (0.9 * cell_vol[cell_id2] < surfn * jjp)
        corrj = 0.9 * cell_vol[cell_id2] / (surfn * jjp);

      djjpf[face_id][0] *= corrj;
      djjpf[face_id][1] *= corrj;
      djjpf[face_id][2] *= corrj;

    }
  }
}

/*----------------------------------------------------------------------------
 * Evaluate boundary thickness.
 *
 * parameters:
 *   m             <-- pointer to mesh structure
 *   m_quantities  <-- pointer to mesh quantities structures.
 *   b_thickness   --> boundary thickness
 *----------------------------------------------------------------------------*/

static void
_b_thickness(const cs_mesh_t             *m,
             const cs_mesh_quantities_t  *mq,
             cs_real_t                    b_thickness[])
{
  const cs_real_3_t  *cell_cen
    = (const cs_real_3_t  *)(mq->cell_cen);
  const cs_real_3_t  *b_face_cog
    = (const cs_real_3_t  *)(mq->b_face_cog);
  const cs_real_3_t  *b_face_normal
    = (const cs_real_3_t  *)(mq->b_face_normal);
  const cs_real_t  *b_face_surf
    = (const cs_real_t *)(mq->b_face_surf);

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t c_id = m->b_face_cells[f_id];
    b_thickness[f_id]
      = (  (b_face_cog[f_id][0] - cell_cen[c_id][0])*b_face_normal[f_id][0]
         + (b_face_cog[f_id][1] - cell_cen[c_id][1])*b_face_normal[f_id][1]
         + (b_face_cog[f_id][2] - cell_cen[c_id][2])*b_face_normal[f_id][2])
        * 2.0 / b_face_surf[f_id];
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set behavior for computing the cocg matrixes for the iterative algo
 * and for the least squares method for scalar and vector gradients.
 *
 * Fortran interface :
 *
 * subroutine comcoc (imrgra)
 * *****************
 *
 * integer          imrgra        : <-- : gradient reconstruction option
 *----------------------------------------------------------------------------*/

void
CS_PROCF (comcoc, COMCOC) (const cs_int_t  *const imrgra)
{
  cs_mesh_quantities_set_cocg_options(*imrgra);
}

/*----------------------------------------------------------------------------
 * Set porous model
 *
 * Fortran interface :
 *
 * subroutine compor (iporos)
 * *****************
 *
 * integer          iporos        : <-- : porous model
 *----------------------------------------------------------------------------*/

void
CS_PROCF (compor, COMPOR) (const cs_int_t  *const iporos)
{
  cs_mesh_quantities_set_porous_model(*iporos);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query or modification of the option for computing cell centers.
 *
 * \param[in]  algo_choice  < 0 : query
 *                            0 : computation based on face centers (default)
 *                            1 : computation by cell sub-volumes
 *
 * \return  0 or 1 according to the selected algorithm
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_quantities_cell_cen_choice(int  algo_choice)
{
  if (algo_choice > -1 && algo_choice < 2)
    _cell_cen_algorithm = algo_choice;

  return _cell_cen_algorithm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query or modification of the option for computing face centers.
 *
 * \param[in]  algo_choice  < 0 : query
 *                            0 : standard computation
 *                            1 : use adjustment for volume
 *                                from versions 1.1 to 5.3
 *
 * \return  0 or 1 according to the selected algorithm
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_quantities_face_cog_choice(int  algo_choice)
{
  if (algo_choice > -1 && algo_choice < 2)
    _ajust_face_cog_compat_v11_v52 = algo_choice;

  return _ajust_face_cog_compat_v11_v52;
}

/*----------------------------------------------------------------------------
 * Compute cocg for iterative gradient reconstruction for scalars.
 *
 * parameters:
 *   gradient_option <-- gradient option (Fortran IMRGRA)
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_set_cocg_options(int  gradient_option)
{
  int _gradient_option = CS_ABS(gradient_option);

  assert(_gradient_option <= 16);

  switch (_gradient_option) {
  case 0:
    _compute_cocg_s_it = true;
    break;
  case 1:
  case 2:
  case 3:
    _compute_cocg_lsq = true;
    break;
  case 4:
  case 5:
  case 6:
    _compute_cocg_lsq = true;
    break;
  /* deprecated options */
  case 10:
    _compute_cocg_s_it = true;
    break;
  case 11:
  case 12:
  case 13:
    _compute_cocg_lsq = true;
    break;
  case 14:
  case 15:
  case 16:
    _compute_cocg_s_it = true;
    _compute_cocg_lsq = true;
    break;

  default:
    break;
  }

  if (gradient_option < 0)
    _compute_cocg_s_it = true;

  _compute_cocg_it = _compute_cocg_s_it;
}

/*----------------------------------------------------------------------------
 * Compute Fluid volumes and fluid surface in addition to
 * cell volumes and surfaces.
 *
 * parameters:
 *   porous_model <-- gradient option (Fortran iporos)
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_set_porous_model(int  porous_model)
{
  cs_glob_porous_model = porous_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a mesh quantities structure.
 *
 * \return  pointer to created cs_mesh_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_quantities_t *
cs_mesh_quantities_create(void)
{
  cs_mesh_quantities_t  *mesh_quantities = NULL;

  BFT_MALLOC(mesh_quantities, 1, cs_mesh_quantities_t);

  mesh_quantities->cell_cen = NULL;
  mesh_quantities->cell_vol = NULL;
  mesh_quantities->cell_f_vol = NULL;
  mesh_quantities->i_face_normal = NULL;
  mesh_quantities->b_face_normal = NULL;
  mesh_quantities->i_f_face_normal = NULL;
  mesh_quantities->b_f_face_normal = NULL;
  mesh_quantities->i_face_cog = NULL;
  mesh_quantities->b_face_cog = NULL;
  mesh_quantities->i_face_surf = NULL;
  mesh_quantities->b_face_surf = NULL;
  mesh_quantities->i_f_face_surf = NULL;
  mesh_quantities->b_f_face_surf = NULL;
  mesh_quantities->i_f_face_factor = NULL;
  mesh_quantities->b_f_face_factor = NULL;
  mesh_quantities->i_dist = NULL;
  mesh_quantities->b_dist = NULL;
  mesh_quantities->weight = NULL;
  mesh_quantities->dijpf = NULL;
  mesh_quantities->diipb = NULL;
  mesh_quantities->dofij = NULL;
  mesh_quantities->diipf = NULL;
  mesh_quantities->djjpf = NULL;
  mesh_quantities->cocgb_s_it = NULL;
  mesh_quantities->cocg_s_it = NULL;
  mesh_quantities->cocgb_s_lsq = NULL;
  mesh_quantities->cocg_it = NULL;
  mesh_quantities->cocg_lsq = NULL;
  mesh_quantities->corr_grad_lin_det = NULL;
  mesh_quantities->corr_grad_lin = NULL;
  mesh_quantities->b_sym_flag = NULL;
  mesh_quantities->c_solid_flag = NULL;
  mesh_quantities->bad_cell_flag = NULL;

  return (mesh_quantities);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a mesh quantities structure.
 *
 * \param[in]  mq  pointer to mesh quantities structures
 *
 * \return  NULL
 */
/*----------------------------------------------------------------------------*/

cs_mesh_quantities_t *
cs_mesh_quantities_destroy(cs_mesh_quantities_t  *mesh_quantities)
{
  cs_mesh_quantities_free_all(mesh_quantities);

  BFT_FREE(mesh_quantities);

  return (mesh_quantities);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a mesh quantities structure to its empty initial state.
 *
 * \param[in]  mq  pointer to mesh quantities structures
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_free_all(cs_mesh_quantities_t  *mq)
{
  BFT_FREE(mq->cell_cen);
  BFT_FREE(mq->cell_vol);
  if (cs_glob_porous_model > 0)
    BFT_FREE(mq->cell_f_vol);
  BFT_FREE(mq->i_face_normal);
  BFT_FREE(mq->b_face_normal);
  if (cs_glob_porous_model == 3) {
    BFT_FREE(mq->i_f_face_normal);
    BFT_FREE(mq->b_f_face_normal);
  }
  BFT_FREE(mq->i_face_cog);
  BFT_FREE(mq->b_face_cog);
  BFT_FREE(mq->i_face_surf);
  BFT_FREE(mq->b_face_surf);
  if (cs_glob_porous_model == 3) {
    BFT_FREE(mq->i_f_face_surf);
    BFT_FREE(mq->b_f_face_surf);
    BFT_FREE(mq->i_f_face_factor);
    BFT_FREE(mq->b_f_face_factor);
  }
  BFT_FREE(mq->i_dist);
  BFT_FREE(mq->b_dist);
  BFT_FREE(mq->weight);
  BFT_FREE(mq->dijpf);
  BFT_FREE(mq->diipb);
  BFT_FREE(mq->dofij);
  BFT_FREE(mq->diipf);
  BFT_FREE(mq->djjpf);
  BFT_FREE(mq->cocgb_s_it);
  BFT_FREE(mq->cocg_s_it);
  BFT_FREE(mq->cocgb_s_lsq);
  BFT_FREE(mq->cocg_it);
  BFT_FREE(mq->cocg_lsq);
  BFT_FREE(mq->corr_grad_lin_det);
  BFT_FREE(mq->corr_grad_lin);
  BFT_FREE(mq->b_sym_flag);
  BFT_FREE(mq->c_solid_flag);
  BFT_FREE(mq->bad_cell_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute mesh quantities needed for preprocessing.
 *
 * \param[in]       m   pointer to mesh structure
 * \param[in, out]  mq  pointer to mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute_preprocess(const cs_mesh_t       *mesh,
                                      cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  cs_lnum_t  n_b_faces = mesh->n_b_faces;
  cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  /* If this is not an update, allocate members of the structure */

  if (mesh_quantities->i_face_normal == NULL)
    BFT_MALLOC(mesh_quantities->i_face_normal, n_i_faces*3, cs_real_t);

  if (mesh_quantities->i_face_cog == NULL)
    BFT_MALLOC(mesh_quantities->i_face_cog, n_i_faces*3, cs_real_t);

  if (mesh_quantities->b_face_normal == NULL)
    BFT_MALLOC(mesh_quantities->b_face_normal, n_b_faces*3, cs_real_t);

  if (mesh_quantities->b_face_cog == NULL)
    BFT_MALLOC(mesh_quantities->b_face_cog, n_b_faces*3, cs_real_t);

  if (mesh_quantities->cell_cen == NULL)
    BFT_MALLOC(mesh_quantities->cell_cen, n_cells_with_ghosts*3, cs_real_t);

  if (mesh_quantities->cell_vol == NULL)
    BFT_MALLOC(mesh_quantities->cell_vol, n_cells_with_ghosts, cs_real_t);

  if (mesh_quantities->i_face_surf == NULL)
    BFT_MALLOC(mesh_quantities->i_face_surf, n_i_faces, cs_real_t);

  if (mesh_quantities->b_face_surf == NULL)
    BFT_MALLOC(mesh_quantities->b_face_surf, n_b_faces, cs_real_t);

  /* Compute face centers of gravity, normals, and surfaces */

  _compute_face_quantities(n_i_faces,
                           (const cs_real_3_t *)mesh->vtx_coord,
                           mesh->i_face_vtx_idx,
                           mesh->i_face_vtx_lst,
                           (cs_real_3_t *)mesh_quantities->i_face_cog,
                           (cs_real_3_t *)mesh_quantities->i_face_normal);

  _compute_face_surface(n_i_faces,
                        mesh_quantities->i_face_normal,
                        mesh_quantities->i_face_surf);

  _compute_face_quantities(n_b_faces,
                           (const cs_real_3_t *)mesh->vtx_coord,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           (cs_real_3_t *)mesh_quantities->b_face_cog,
                           (cs_real_3_t *)mesh_quantities->b_face_normal);

  _compute_face_surface(n_b_faces,
                        mesh_quantities->b_face_normal,
                        mesh_quantities->b_face_surf);

  if (cs_glob_mesh_quantities_flag & CS_FACE_CENTER_REFINE) {
    _refine_warped_face_centers
      (n_i_faces,
       (const cs_real_3_t *)mesh->vtx_coord,
       mesh->i_face_vtx_idx,
       mesh->i_face_vtx_lst,
       (cs_real_3_t *)mesh_quantities->i_face_cog,
       (const cs_real_3_t *)mesh_quantities->i_face_normal);

    _refine_warped_face_centers
      (n_b_faces,
       (const cs_real_3_t *)mesh->vtx_coord,
       mesh->b_face_vtx_idx,
       mesh->b_face_vtx_lst,
       (cs_real_3_t *)mesh_quantities->b_face_cog,
       (const cs_real_3_t *)mesh_quantities->b_face_normal);
  }

  if (_ajust_face_cog_compat_v11_v52) {
    _adjust_face_cog_v11_v52
      (n_i_faces,
       (const cs_real_3_t *)mesh->vtx_coord,
       mesh->i_face_vtx_idx,
       mesh->i_face_vtx_lst,
       (cs_real_3_t *)mesh_quantities->i_face_cog,
       (const cs_real_3_t *)mesh_quantities->i_face_normal);

    _adjust_face_cog_v11_v52
      (n_b_faces,
       (const cs_real_3_t *)mesh->vtx_coord,
       mesh->b_face_vtx_idx,
       mesh->b_face_vtx_lst,
       (cs_real_3_t *)mesh_quantities->b_face_cog,
       (const cs_real_3_t *)mesh_quantities->b_face_normal);
  }

  /* Compute cell centers from face barycenters or vertices */

  bool volume_computed = false;

  switch (_cell_cen_algorithm) {

  case 0:
    cs_mesh_quantities_cell_faces_cog(mesh,
                                      mesh_quantities->i_face_normal,
                                      mesh_quantities->i_face_cog,
                                      mesh_quantities->b_face_normal,
                                      mesh_quantities->b_face_cog,
                                      mesh_quantities->cell_cen);

    break;
  case 1:
    _compute_cell_quantities(mesh,
                             (const cs_real_3_t *)mesh_quantities->i_face_normal,
                             (const cs_real_3_t *)mesh_quantities->i_face_cog,
                             (const cs_real_3_t *)mesh_quantities->b_face_normal,
                             (const cs_real_3_t *)mesh_quantities->b_face_cog,
                             (cs_real_3_t *)mesh_quantities->cell_cen,
                             mesh_quantities->cell_vol);
    volume_computed = true;
    break;

  default:
    assert(0);

  }

  if (cs_glob_mesh_quantities_flag & CS_CELL_CENTER_CORRECTION) {
    _recompute_cell_cen_face(mesh,
                             (const cs_real_3_t *)(mesh_quantities->i_face_normal),
                             (const cs_real_3_t *)(mesh_quantities->i_face_cog),
                             (const cs_real_3_t *)(mesh_quantities->b_face_normal),
                             (const cs_real_3_t *)(mesh_quantities->b_face_cog),
                             (cs_real_3_t *)(mesh_quantities->cell_cen));
    volume_computed = false; /* should not be different with plane faces,
                                not sure with warped faces */
  }

  /* Recompute face centers as the middle of two cell centers if possible */
  if (cs_glob_mesh_quantities_flag & CS_CELL_FACE_CENTER_CORRECTION) {

     if (mesh->halo != NULL) {
       cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                                mesh_quantities->cell_cen, 3);
       if (mesh->n_init_perio > 0)
         cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                                   mesh_quantities->cell_cen);
     }

    _correct_cell_face_center(mesh,
                              n_cells_with_ghosts,
                              n_i_faces,
                              n_b_faces,
                              (const cs_lnum_2_t *)(mesh->i_face_cells),
                              mesh->b_face_cells,
                              (cs_real_3_t *)(mesh_quantities->cell_cen),
                              (cs_real_3_t *)(mesh_quantities->i_face_cog),
                              (cs_real_3_t *)(mesh_quantities->b_face_cog),
                              (cs_real_3_t *)(mesh_quantities->i_face_normal),
                              (cs_real_3_t *)(mesh_quantities->b_face_normal));

    volume_computed = false;

  }

  /* Compute the volume of cells */

  if (volume_computed == false)
    _compute_cell_volume
      (mesh,
       (const cs_real_3_t *)(mesh_quantities->i_face_normal),
       (const cs_real_3_t *)(mesh_quantities->i_face_cog),
       (const cs_real_3_t *)(mesh_quantities->b_face_normal),
       (const cs_real_3_t *)(mesh_quantities->b_face_cog),
       (const cs_real_3_t *)(mesh_quantities->cell_cen),
       mesh_quantities->cell_vol);


  /* Correction of small or negative volumes
     (doesn't conserve the total volume) */

  if (cs_glob_mesh_quantities_flag & CS_CELL_VOLUME_RATIO_CORRECTION)
    _cell_bad_volume_correction(mesh,
                                mesh_quantities->cell_vol);

  /* Synchronize geometric quantities */

  if (mesh->halo != NULL) {

    cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                             mesh_quantities->cell_cen, 3);
    if (mesh->n_init_perio > 0)
      cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                                mesh_quantities->cell_cen);

    cs_halo_sync_var(mesh->halo, CS_HALO_EXTENDED, mesh_quantities->cell_vol);

  }

  _cell_volume_reductions(mesh,
                          mesh_quantities->cell_vol,
                          &(mesh_quantities->min_vol),
                          &(mesh_quantities->max_vol),
                          &(mesh_quantities->tot_vol));

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    cs_real_t  _min_vol, _max_vol, _tot_vol;

    MPI_Allreduce(&(mesh_quantities->min_vol), &_min_vol, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    MPI_Allreduce(&(mesh_quantities->max_vol), &_max_vol, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    MPI_Allreduce(&(mesh_quantities->tot_vol), &_tot_vol, 1, CS_MPI_REAL,
                  MPI_SUM, cs_glob_mpi_comm);

    mesh_quantities->min_vol = _min_vol;
    mesh_quantities->max_vol = _max_vol;
    mesh_quantities->tot_vol = _tot_vol;

  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute mesh quantities.
 *
 * \param[in]       m   pointer to mesh structure
 * \param[in, out]  mq  pointer to mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute(const cs_mesh_t       *mesh,
                           cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  dim = mesh->dim;
  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  cs_lnum_t  n_b_faces = mesh->n_b_faces;
  cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  /* Update the number of passes */

  _n_computations++;

  cs_mesh_quantities_compute_preprocess(mesh, mesh_quantities);

  /* Balance porous model */
  if (cs_glob_porous_model == 3) {
    if (mesh_quantities->i_f_face_normal == NULL)
      BFT_MALLOC(mesh_quantities->i_f_face_normal, n_i_faces*dim, cs_real_t);

    if (mesh_quantities->b_f_face_normal == NULL)
      BFT_MALLOC(mesh_quantities->b_f_face_normal, n_b_faces*dim, cs_real_t);

    if (mesh_quantities->i_f_face_surf == NULL)
      BFT_MALLOC(mesh_quantities->i_f_face_surf, n_i_faces, cs_real_t);

    if (mesh_quantities->b_f_face_surf == NULL)
      BFT_MALLOC(mesh_quantities->b_f_face_surf, n_b_faces, cs_real_t);

    if (mesh_quantities->i_f_face_factor == NULL)
      BFT_MALLOC(mesh_quantities->i_f_face_factor, n_i_faces, cs_real_2_t);

    if (mesh_quantities->b_f_face_factor == NULL)
      BFT_MALLOC(mesh_quantities->b_f_face_factor, n_b_faces, cs_real_t);

  }
  else {
    mesh_quantities->i_f_face_normal = mesh_quantities->i_face_normal;
    mesh_quantities->b_f_face_normal = mesh_quantities->b_face_normal;
    mesh_quantities->i_f_face_surf = mesh_quantities->i_face_surf;
    mesh_quantities->b_f_face_surf = mesh_quantities->b_face_surf;
  }

  /* Porous models */
  if (cs_glob_porous_model > 0) {
    if (mesh_quantities->cell_f_vol == NULL)
      BFT_MALLOC(mesh_quantities->cell_f_vol, n_cells_with_ghosts, cs_real_t);

    if (mesh_quantities->c_solid_flag == NULL) {
      BFT_MALLOC(mesh_quantities->c_solid_flag, n_cells_with_ghosts, cs_int_t);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++)
        mesh_quantities->c_solid_flag[cell_id] = 0;
    }

  }
  else {
    mesh_quantities->cell_f_vol = mesh_quantities->cell_vol;

    if (mesh_quantities->c_solid_flag == NULL) {
      BFT_MALLOC(mesh_quantities->c_solid_flag, 1, cs_int_t);
      mesh_quantities->c_solid_flag[0] = 0;
    }
  }

  if (mesh_quantities->i_dist == NULL)
    BFT_MALLOC(mesh_quantities->i_dist, n_i_faces, cs_real_t);

  if (mesh_quantities->b_dist == NULL)
    BFT_MALLOC(mesh_quantities->b_dist, n_b_faces, cs_real_t);

  if (mesh_quantities->weight == NULL)
    BFT_MALLOC(mesh_quantities->weight, n_i_faces, cs_real_t);

  if (mesh_quantities->dijpf == NULL)
    BFT_MALLOC(mesh_quantities->dijpf, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->diipb == NULL)
    BFT_MALLOC(mesh_quantities->diipb, n_b_faces*dim, cs_real_t);

  if (mesh_quantities->dofij == NULL)
    BFT_MALLOC(mesh_quantities->dofij, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->diipf == NULL)
    BFT_MALLOC(mesh_quantities->diipf, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->djjpf == NULL)
    BFT_MALLOC(mesh_quantities->djjpf, n_i_faces*dim, cs_real_t);


  /* Compute 3x3 cocg dimensionless matrix */

  if (_compute_cocg_it == 1) {
    if (mesh_quantities->cocg_it == NULL)
      BFT_MALLOC(mesh_quantities->cocg_it, n_cells_with_ghosts, cs_real_33_t);
  }

  if (_compute_cocg_lsq == 1) {
    if (mesh_quantities->cocg_lsq == NULL) {
      BFT_MALLOC(mesh_quantities->cocg_lsq, n_cells_with_ghosts, cs_real_33_t);
    }
  }

  if (mesh_quantities->b_sym_flag == NULL)
    BFT_MALLOC(mesh_quantities->b_sym_flag, n_b_faces, cs_int_t);

  /* Compute some distances relative to faces and associated weighting */

  _compute_face_distances(mesh->n_i_faces,
                          mesh->n_b_faces,
                          (const cs_lnum_2_t *)(mesh->i_face_cells),
                          mesh->b_face_cells,
                          (const cs_real_3_t *)(mesh_quantities->i_face_normal),
                          (const cs_real_3_t *)(mesh_quantities->b_face_normal),
                          (const cs_real_3_t *)(mesh_quantities->i_face_cog),
                          (const cs_real_3_t *)(mesh_quantities->b_face_cog),
                          (const cs_real_3_t *)(mesh_quantities->cell_cen),
                          mesh_quantities->i_dist,
                          mesh_quantities->b_dist,
                          mesh_quantities->weight);

  /* Compute some vectors relative to faces to handle non-orthogonalities */

  _compute_face_vectors(dim,
                        mesh->n_i_faces,
                        mesh->n_b_faces,
                        (const cs_lnum_2_t *)(mesh->i_face_cells),
                        mesh->b_face_cells,
                        mesh_quantities->i_face_normal,
                        mesh_quantities->b_face_normal,
                        mesh_quantities->i_face_cog,
                        mesh_quantities->b_face_cog,
                        mesh_quantities->i_face_surf,
                        mesh_quantities->b_face_surf,
                        mesh_quantities->cell_cen,
                        mesh_quantities->weight,
                        mesh_quantities->b_dist,
                        mesh_quantities->dijpf,
                        mesh_quantities->diipb,
                        mesh_quantities->dofij);

  /* Compute additional vectors relative to faces to handle non-orthogonalities */

  _compute_face_sup_vectors
    (mesh->n_i_faces,
     (const cs_lnum_2_t *)(mesh->i_face_cells),
     (const cs_real_3_t *)(mesh_quantities->i_face_normal),
     (const cs_real_3_t *)(mesh_quantities->i_face_cog),
     (const cs_real_3_t *)(mesh_quantities->cell_cen),
     mesh_quantities->cell_vol,
     mesh_quantities->i_dist,
     (cs_real_3_t *)(mesh_quantities->diipf),
     (cs_real_3_t *)(mesh_quantities->djjpf));

  /* Compute 3x3 cocg matrixes */

  if (_compute_cocg_s_it == 1)
    _compute_cell_cocg_s_it(mesh, mesh_quantities);

  if (_compute_cocg_lsq == 1)
    _compute_cell_cocg_lsq(mesh, mesh_quantities, NULL);

  if (_compute_cocg_it == 1)
    _compute_cell_cocg_it(mesh, mesh_quantities, NULL);

  /* Build the geometrical matrix linear gradient correction */
  if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_WARPED_CORRECTION)
    _compute_corr_grad_lin(mesh, mesh_quantities);

  /* Print some information on the control volumes, and check min volume */

  if (_n_computations == 1)
    bft_printf(_(" --- Information on the volumes\n"
                 "       Minimum control volume      = %14.7e\n"
                 "       Maximum control volume      = %14.7e\n"
                 "       Total volume for the domain = %14.7e\n"),
               mesh_quantities->min_vol, mesh_quantities->max_vol,
               mesh_quantities->tot_vol);
  else {
    if (mesh_quantities->min_vol <= 0.) {
      bft_printf(_(" --- Information on the volumes\n"
                   "       Minimum control volume      = %14.7e\n"
                   "       Maximum control volume      = %14.7e\n"
                   "       Total volume for the domain = %14.7e\n"),
                 mesh_quantities->min_vol, mesh_quantities->max_vol,
                 mesh_quantities->tot_vol);
      bft_printf(_("\nAbort due to the detection of a negative control "
                   "volume.\n"));
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute fluid mesh quantities
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_fluid_compute(const cs_mesh_t       *mesh,
                                 cs_mesh_quantities_t  *mesh_quantities)
{
  CS_UNUSED(mesh);
  CS_UNUSED(mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Compute fluid section mesh quantities at the initial step
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_fluid_sections(const cs_mesh_t       *mesh,
                            cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  cs_lnum_t  n_b_faces = mesh->n_b_faces;

  cs_real_3_t *restrict i_face_normal
    = (cs_real_3_t *restrict)mesh_quantities->i_face_normal;
  cs_real_3_t *restrict b_face_normal
    = (cs_real_3_t *restrict)mesh_quantities->b_face_normal;
  cs_real_3_t *restrict i_f_face_normal
    = (cs_real_3_t *restrict)mesh_quantities->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal
    = (cs_real_3_t *restrict)mesh_quantities->b_f_face_normal;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    mesh_quantities->i_f_face_surf[face_id] =
      mesh_quantities->i_face_surf[face_id];

    for (int i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = i_face_normal[face_id][i];

    mesh_quantities->i_f_face_factor[face_id][0] = 1.;
    mesh_quantities->i_f_face_factor[face_id][1] = 1.;
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    mesh_quantities->b_f_face_surf[face_id]
      = mesh_quantities->b_face_surf[face_id];

    for (int i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = b_face_normal[face_id][i];

    mesh_quantities->b_f_face_factor[face_id] = 1.;
  }
}

/*----------------------------------------------------------------------------
 * Compute mesh quantities -> vectors II' and JJ'
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_sup_vectors(const cs_mesh_t       *mesh,
                               cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  dim = mesh->dim;
  cs_lnum_t  n_i_faces = mesh->n_i_faces;

  if (mesh_quantities->diipf == NULL)
    BFT_MALLOC(mesh_quantities->diipf, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->djjpf == NULL)
    BFT_MALLOC(mesh_quantities->djjpf, n_i_faces*dim, cs_real_t);

  _compute_face_sup_vectors
    (mesh->n_i_faces,
     (const cs_lnum_2_t *)(mesh->i_face_cells),
     (const cs_real_3_t *)(mesh_quantities->i_face_normal),
     (const cs_real_3_t *)(mesh_quantities->i_face_cog),
     (const cs_real_3_t *)(mesh_quantities->cell_cen),
     mesh_quantities->cell_vol,
     mesh_quantities->i_dist,
     (cs_real_3_t *)(mesh_quantities->diipf),
     (cs_real_3_t *)(mesh_quantities->djjpf));
}

/*----------------------------------------------------------------------------
 * Compute internal and border face normal.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_normal <-> pointer to the internal face normal array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_face_normal(const cs_mesh_t   *mesh,
                               cs_real_t         *p_i_face_normal[],
                               cs_real_t         *p_b_face_normal[])
{
  cs_real_t  *i_face_normal = NULL, *b_face_normal = NULL;

  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;

  /* Internal face treatment */

  BFT_MALLOC(i_face_normal, n_i_faces * 3, cs_real_t);

  _compute_face_normal(mesh->n_i_faces,
                       (const cs_real_3_t *)mesh->vtx_coord,
                       mesh->i_face_vtx_idx,
                       mesh->i_face_vtx_lst,
                       (cs_real_3_t *)i_face_normal);

  *p_i_face_normal = i_face_normal;

  /* Boundary face treatment */

  BFT_MALLOC(b_face_normal, n_b_faces * 3, cs_real_t);

  _compute_face_normal(mesh->n_b_faces,
                       (const cs_real_3_t *)mesh->vtx_coord,
                       mesh->b_face_vtx_idx,
                       mesh->b_face_vtx_lst,
                       (cs_real_3_t *)b_face_normal);

  *p_b_face_normal = b_face_normal;
}

/*----------------------------------------------------------------------------
 * Compute interior face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_cog    <-> pointer to the interior face center array
 *   p_i_face_normal <-> pointer to the interior face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_i_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_i_face_cog[],
                           cs_real_t         *p_i_face_normal[])
{
  cs_real_t  *i_face_cog = NULL, *i_face_normal = NULL;

  BFT_MALLOC(i_face_cog, mesh->n_i_faces * mesh->dim, cs_real_t);
  BFT_MALLOC(i_face_normal, mesh->n_i_faces * mesh->dim, cs_real_t);

  _compute_face_quantities(mesh->n_i_faces,
                           (const cs_real_3_t *)mesh->vtx_coord,
                           mesh->i_face_vtx_idx,
                           mesh->i_face_vtx_lst,
                           (cs_real_3_t *)i_face_cog,
                           (cs_real_3_t *)i_face_normal);

  *p_i_face_cog = i_face_cog;
  *p_i_face_normal = i_face_normal;
}

/*----------------------------------------------------------------------------
 * Compute border face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_b_face_cog    <-> pointer to the border face center array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_b_face_cog[],
                           cs_real_t         *p_b_face_normal[])
{
  cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;

  BFT_MALLOC(b_face_cog, mesh->n_b_faces * mesh->dim, cs_real_t);
  BFT_MALLOC(b_face_normal, mesh->n_b_faces * mesh->dim, cs_real_t);

  _compute_face_quantities(mesh->n_b_faces,
                           (const cs_real_3_t *)mesh->vtx_coord,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           (cs_real_3_t *)b_face_cog,
                           (cs_real_3_t *)b_face_normal);

  *p_b_face_cog = b_face_cog;
  *p_b_face_normal = b_face_normal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute approximate cells centers as the mean of the given face
 *         centers weighted by the associated surfaces.
 *
 *           n-1
 *           Sum  Surf(Fi) G(Fi)
 *           i=0
 *  G(C) = -----------------------
 *           n-1
 *           Sum  Surf(Fi)
 *           i=0
 *
 * \param[in]   mesh         pointer to mesh structure
 * \param[in]   i_face_norm  surface normal of internal faces
 * \param[in]   i_face_cog   center of gravity of internal faces
 * \param[in]   b_face_norm  surface normal of border faces
 * \param[in]   b_face_cog   center of gravity of border faces
 * \param[out]  cell_cen     cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_cell_faces_cog(const cs_mesh_t  *mesh,
                                  const cs_real_t   i_face_norm[],
                                  const cs_real_t   i_face_cog[],
                                  const cs_real_t   b_face_norm[],
                                  const cs_real_t   b_face_cog[],
                                  cs_real_t         cell_cen[])
{
  cs_real_t  *cell_area = NULL;

  /* Mesh connectivity */

  const  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const  cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const  cs_lnum_t  n_cells = mesh->n_cells;
  const  cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;
  const  cs_lnum_2_t  *i_face_cells
    = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const  cs_lnum_t  *b_face_cells = mesh->b_face_cells;

  /* Return if ther is not enough data (Solcom case except rediative module
     or Pre-processor 1.2.d without option "-n") */

  if (mesh->i_face_vtx_lst == NULL && mesh->b_face_vtx_lst == NULL)
    return;

  /* Checking */

  assert(cell_cen != NULL);

  /* Initialization */

  BFT_MALLOC(cell_area, n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t j = 0; j < n_cells_with_ghosts; j++) {

    cell_area[j] = 0.;

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[3*j + i] = 0.;

  }

  /* Loop on interior faces
     ---------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    /* For each cell sharing the internal face, we update
     * cell_cen and cell_area */

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    /* Computation of the area of the face */

    cs_real_t area = cs_math_3_norm(i_face_norm + 3*f_id);

    if (c_id1 > -1) {
      cell_area[c_id1] += area;
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[3*c_id1 + i] += i_face_cog[3*f_id + i]*area;
    }
    if (c_id2 > -1) {
      cell_area[c_id2] += area;
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[3*c_id2 + i] += i_face_cog[3*f_id + i]*area;
    }

  } /* End of loop on interior faces */

  /* Loop on boundary faces
     --------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* For each cell sharing a border face, we update the numerator
     * of cell_cen and cell_area */

    cs_lnum_t c_id1 = b_face_cells[f_id];

    /* Computation of the area of the face
       (note that c_id1 == -1 may happen for isolated faces,
       which are cleaned afterwards) */

    if (c_id1 > -1) {

      cs_real_t area = cs_math_3_norm(b_face_norm + 3*f_id);

      cell_area[c_id1] += area;

      /* Computation of the numerator */

      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[3*c_id1 + i] += b_face_cog[3*f_id + i]*area;

    }

  } /* End of loop on boundary faces */

  /* Loop on cells to finalize the computation of center of gravity
     -------------------------------------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[c_id*3 + i] /= cell_area[c_id];

  }

  /* Free memory */

  BFT_FREE(cell_area);
}

/*----------------------------------------------------------------------------
 * Compute cell volumes.
 *
 * The corresponding array is allocated by this function, and it is the
 * caller's responsability to free it when they are no longer needed.
 *
 * parameters:
 *   mesh     <-- pointer to a cs_mesh_t structure
 *
 * return:
 *   pointer to newly allocated cell volumes array
 *----------------------------------------------------------------------------*/

cs_real_t *
cs_mesh_quantities_cell_volume(const cs_mesh_t  *mesh)
{
  cs_real_t *cell_vol;
  BFT_MALLOC(cell_vol, mesh->n_cells_with_ghosts, cs_real_t);

  cs_real_3_t *cell_cen;
  BFT_MALLOC(cell_cen, mesh->n_cells_with_ghosts, cs_real_3_t);

  cs_real_t  *i_face_cog = NULL, *i_face_normal = NULL;
  cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;

  cs_mesh_quantities_i_faces(mesh, &i_face_cog, &i_face_normal);
  cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

  _compute_cell_quantities(mesh,
                           (const cs_real_3_t *)i_face_normal,
                           (const cs_real_3_t *)i_face_cog,
                           (const cs_real_3_t *)b_face_normal,
                           (const cs_real_3_t *)b_face_cog,
                           (cs_real_3_t *)cell_cen,
                           cell_vol);

  BFT_FREE(cell_cen);
  BFT_FREE(b_face_normal);
  BFT_FREE(b_face_cog);
  BFT_FREE(i_face_normal);
  BFT_FREE(i_face_cog);

  return cell_vol;
}

/*----------------------------------------------------------------------------
 * Check that no negative volumes are present, and exit on error otherwise.
 *
 * parameters:
 *   mesh            <-- pointer to mesh structure
 *   mesh_quantities <-- pointer to mesh quantities structure
 *   allow_error     <-- 1 if errors are allowed, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_check_vol(const cs_mesh_t             *mesh,
                             const cs_mesh_quantities_t  *mesh_quantities,
                             int                          allow_error)
{
  cs_lnum_t  cell_id;

  cs_gnum_t  error_count = 0;

  for (cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
    if (mesh_quantities->cell_vol[cell_id] < 0.0)
      error_count += 1;
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t tot_error_count = 0;
    MPI_Allreduce(&error_count, &tot_error_count, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    error_count = tot_error_count;
  }
#endif

  /* Exit with error */

  if (error_count > 0) {
    const char fmt[]
      = N_("  %llu cells have a Negative volume.\n"
           " Run mesh quality check for post-processing output.\n"
           " In case of mesh joining, this may be due to overly "
           " agressive joining parameters.");

    if (allow_error) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_(fmt), (unsigned long long)error_count);
      bft_printf("\n\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(fmt), (unsigned long long)error_count);
  }
}

/*----------------------------------------------------------------------------
 * Update mesh quantities relative to extended ghost cells when the
 * neighborhood is reduced.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_reduce_extended(const cs_mesh_t       *mesh,
                                   cs_mesh_quantities_t  *mesh_quantities)
{
  if (_compute_cocg_lsq == 1)
    _compute_cell_cocg_lsq(mesh, mesh_quantities, NULL);
}

/*----------------------------------------------------------------------------
 * Return the number of times mesh quantities have been computed.
 *
 * returns:
 *   number of times mesh quantities have been computed
 *----------------------------------------------------------------------------*/

int
cs_mesh_quantities_compute_count(void)
{
  return _n_computations;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine local boundary thickness around each vertex.
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   mq           pointer to mesh quantities structures.
 * \param[in]   n_passes     number of smoothing passes
 * \param[out]  b_thickness  thickness for each mesh vertex
 *                           (0 at non-boundary vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_thickness_v(const cs_mesh_t             *m,
                                 const cs_mesh_quantities_t  *mq,
                                 int                          n_passes,
                                 cs_real_t                    b_thickness[])
{
  cs_real_t *v_sum = NULL;
  cs_real_t *f_b_thickness = NULL;

  BFT_MALLOC(v_sum, m->n_vertices*2, cs_real_t);

  BFT_MALLOC(f_b_thickness, m->n_b_faces*2, cs_real_t);
  _b_thickness(m, mq, f_b_thickness);

  if (n_passes < 1)
    n_passes = 1;

  for (int i = 0; i < n_passes; i++) {

    for (cs_lnum_t j = 0; j < m->n_vertices*2; j++)
      v_sum[j] = 0.;

    for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      const cs_real_t f_s = mq->b_face_surf[f_id];
      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t v_id = m->b_face_vtx_lst[k];
        v_sum[v_id*2]   += f_s * f_b_thickness[f_id];
        v_sum[v_id*2+1] += f_s;
      }
    }

    if (m->vtx_interfaces != NULL)
      cs_interface_set_sum(m->vtx_interfaces,
                           m->n_vertices,
                           2,
                           true,
                           CS_REAL_TYPE,
                           v_sum);

    /* Prepare for next smoothing */

    if (i < n_passes-1) {

      for (cs_lnum_t j = 0; j < m->n_b_faces*2; j++)
        f_b_thickness[j] = 0.;

      for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
        cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
        cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
        for (cs_lnum_t k = s_id; k < e_id; k++) {
          cs_lnum_t v_id = m->b_face_vtx_lst[k];
          f_b_thickness[f_id] += v_sum[v_id*2];
          f_b_thickness[f_id + m->n_b_faces] += v_sum[v_id*2 + 1];
        }
      }

      for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
        if (f_b_thickness[f_id + m->n_b_faces] > 0)
          f_b_thickness[f_id] /= f_b_thickness[f_id + m->n_b_faces];
      }

    }

  }

  BFT_FREE(f_b_thickness);

  for (cs_lnum_t j = 0; j < m->n_vertices; j++) {
    if (v_sum[j*2+1] > 0)
      b_thickness[j] = v_sum[j*2] / v_sum[j*2+1];
    else
      b_thickness[j] = 0;
  }

  BFT_FREE(v_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine local boundary thickness around each boundary face.
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   mq           pointer to mesh quantities structures.
 * \param[in]   n_passes     number of optional smoothing passes
 * \param[out]  b_thickness  thickness for each mesh boundary face
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_thickness_f(const cs_mesh_t             *m,
                                 const cs_mesh_quantities_t  *mq,
                                 int                          n_passes,
                                 cs_real_t                    b_thickness[])
{
  if (n_passes < 1)
    _b_thickness(m, mq, b_thickness);

  else {

    cs_real_t *v_b_thickness = NULL;

    BFT_MALLOC(v_b_thickness, m->n_vertices, cs_real_t);

    cs_mesh_quantities_b_thickness_v(m,
                                     mq,
                                     n_passes,
                                     v_b_thickness);

    for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
      b_thickness[f_id] = 0;
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t v_id = m->b_face_vtx_lst[k];
        b_thickness[f_id] += v_b_thickness[v_id];
      }
      b_thickness[f_id] /= (e_id - s_id);
    }

    BFT_FREE(v_b_thickness);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log mesh quantities options to setup file.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_log_setup(void)
{
  if (cs_glob_mesh_quantities_flag != 0 || _cell_cen_algorithm != 1)
    cs_log_printf(CS_LOG_SETUP, _("\n"
                                  "Mesh quantity computation options\n"
                                  "---------------------------------\n\n"));

  const char *cen_type_name[] = {N_("weighted center of face centers"),
                                 N_("center of mass")};
  cs_log_printf(CS_LOG_SETUP,
                _("  Cell centers: %s\n"),
                _(cen_type_name[_cell_cen_algorithm]));

  if (cs_glob_mesh_quantities_flag != 0) {

    const char *correction_name[] = {"CS_BAD_CELLS_WARPED_CORRECTION",
                                     "CS_BAD_CELLS_WARPED_REGULARISATION",
                                     "CS_CELL_FACE_CENTER_CORRECTION",
                                     "CS_CELL_FACE_CENTER_CORRECTION",
                                     "CS_FACE_DISTANCE_CLIP",
                                     "CS_FACE_RECONSTRUCTION_CLIP",
                                     "CS_CELL_VOLUME_RATIO_CORRECTION",
                                     "CS_FACE_CENTER_REFINE"};

    cs_log_printf(CS_LOG_SETUP,
       ("\n"
        "   Mesh quantity corrections:\n"));

    for (int i = 0; i < 8; i++) {
      if (cs_glob_mesh_quantities_flag & (1 << i))
        cs_log_printf(CS_LOG_SETUP, "      %s\n", correction_name[i]);
    }

  }
}

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_quantities_t structure
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-- pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_dump(const cs_mesh_t             *mesh,
                        const cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  i;

  const cs_lnum_t  n_cells = mesh->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;

  const cs_real_t  *cell_cen = mesh_quantities->cell_cen;
  const cs_real_t  *cell_vol = mesh_quantities->cell_vol;
  const cs_real_t  *i_fac_norm = mesh_quantities->i_face_normal;
  const cs_real_t  *b_fac_norm = mesh_quantities->b_face_normal;
  const cs_real_t  *i_fac_cog = mesh_quantities->i_face_cog;
  const cs_real_t  *b_fac_cog = mesh_quantities->b_face_cog;
  const cs_real_t  *i_fac_surf = mesh_quantities->i_face_surf;
  const cs_real_t  *b_fac_surf = mesh_quantities->b_face_surf;

  bft_printf("\n\nDUMP OF A MESH QUANTITIES STRUCTURE: %p\n\n",
             (const void *)mesh_quantities);

  if (mesh_quantities == NULL)
    return;

  /* Cell data */

  bft_printf("\n\n"
             "    ---------------"
             "    Cell quantities"
             "    ---------------\n\n");

  bft_printf("Cell center coordinates:\n");
  for (i = 0; i < n_cells; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               cell_cen[3*i], cell_cen[3*i+1], cell_cen[3*i+2]);

  bft_printf("\nCell volume:\n");
  for (i = 0; i < n_cells; i++)
    bft_printf("    < %d >    %.3f\n", i+1, cell_vol[i]);

  /* Internal faces data */

  bft_printf("\n\n"
             "    ------------------------"
             "    Interior face quantities"
             "    ------------------------\n\n");

  bft_printf("\nInterior face normals\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               i_fac_norm[3*i], i_fac_norm[3*i+1], i_fac_norm[3*i+2]);

  bft_printf("\nInterior face centers\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               i_fac_cog[3*i], i_fac_cog[3*i+1], i_fac_cog[3*i+2]);

  bft_printf("\nInterior face surfaces\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %d >    %.3f\n", i+1, i_fac_surf[i]);

  /* Border faces data */

  bft_printf("\n\n"
             "    ------------------------"
             "    Boundary face quantities"
             "    ------------------------\n\n");

  bft_printf("\nBoundary face normals\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               b_fac_norm[3*i], b_fac_norm[3*i+1], b_fac_norm[3*i+2]);

  bft_printf("\nBoundary faces centers\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               b_fac_cog[3*i], b_fac_cog[3*i+1], b_fac_cog[3*i+2]);

  bft_printf("\nBoundary face surfaces\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %d >    %.3f\n", i+1, b_fac_surf[i]);

  bft_printf("\n\nEND OF DUMP OF MESH QUANTITIES STRUCTURE\n\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient least squares algorithm
 * adapted for internal coupling.
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *   ce   <->  coupling
 *----------------------------------------------------------------------------*/

void
cs_compute_cell_cocg_lsq_coupling(const cs_mesh_t         *m,
                                  cs_mesh_quantities_t    *fvq,
                                  cs_internal_coupling_t  *ce)
{
  _compute_cell_cocg_lsq(m, fvq, ce);
}

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient iterative algorithm
 * adapted for internal coupling.
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *   ce   <->  coupling
 *----------------------------------------------------------------------------*/

void
cs_compute_cell_cocg_it_coupling(const cs_mesh_t         *m,
                                 cs_mesh_quantities_t    *fvq,
                                 cs_internal_coupling_t  *ce)
{
  _compute_cell_cocg_it(m, fvq, ce);
}

/*----------------------------------------------------------------------------*/

#if 0 /* Test if face orientation is OK */

  cs_lnum_t  i, fac_id, cell_id;
  cs_real_t  cogfac[3];
  cs_real_t  cogcel[3];
  cs_real_t  normal[3];
  cs_real_t  pscal;

  for (fac_id = 0 ; fac_id < mesh->n_b_faces ; fac_id++) {

    cell_id = mesh->b_face_cells[fac_id];
    pscal = 0;

    for (i = 0 ; i < 3 ; i++) {
      cogcel[i]  = cs_glob_mesh_quantities->cell_cen[cell_id*3 + i];
      cogfac[i]  = cs_glob_mesh_quantities->b_face_cog[fac_id*3 + i];
      normal[i] = cs_glob_mesh_quantities->b_face_normal[fac_id*3 + i];
      pscal += normal[i] * (cogfac[i] - cogcel[i]);
    }

    if (pscal < 0.0)
      printf("num_fac_brd = %d, num_cel = %d, pscal = %f\n",
             fac_id + 1, cell_id + 1, pscal);
  }

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS
