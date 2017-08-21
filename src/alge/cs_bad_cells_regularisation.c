/*============================================================================
 * Regulation on bad cells.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 * Code_Saturne library headers
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_bad_cells.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_bad_cells_regularisation.h"

/*----------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

BEGIN_C_DECLS

void
_tag_bad_cells(const cs_mesh_t      *mesh,
               cs_mesh_quantities_t *mq)
{
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = mesh->n_cells;
  cs_lnum_t n_i_faces = mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const int *b_face_cells = mesh->b_face_cells;

  unsigned *bad = mq->bad_cell_flag;
  const cs_real_t *surfn = mq->i_face_surf;
  const cs_real_t *surfbn = mq->b_face_surf;
  double *dist = mq->i_dist;
  double *distbr = mq->b_dist;
  double *volume  = mq->cell_vol;

  const cs_real_3_t *xyzcen = (const cs_real_3_t *) mq->cell_cen;
  const cs_real_3_t *cdgfac = (const cs_real_3_t *) mq->i_face_cog;
  const cs_real_3_t *cdgfbo = (const cs_real_3_t *) mq->b_face_cog;
  const cs_real_3_t *surfac = (const cs_real_3_t *) mq->i_face_normal;
  const cs_real_3_t *surfbo = (const cs_real_3_t *) mq->b_face_normal;

  static cs_gnum_t nb_bad_cells = 0;

  /* If cells are not yet tagged, tag them */
  if (mq->bad_cell_indic == NULL) {
    BFT_MALLOC(mq->bad_cell_indic, n_cells_ext, int);

    cs_field_t *f = cs_field_by_name("regul");

    cs_real_3_t *vol;
    BFT_MALLOC(vol, n_cells_ext, cs_real_3_t);

    //FIXME tensor ?
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      vol[cell_id][0] = 0.;
      vol[cell_id][1] = 0.;
      vol[cell_id][2] = 0.;
    }
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      cs_lnum_t cell_id2 = i_face_cells[face_id][1];
      vol[cell_id1][0] += cdgfac[face_id][0] * surfac[face_id][0];
      vol[cell_id1][1] += cdgfac[face_id][1] * surfac[face_id][1];
      vol[cell_id1][2] += cdgfac[face_id][2] * surfac[face_id][2];
      vol[cell_id2][0] -= cdgfac[face_id][0] * surfac[face_id][0];
      vol[cell_id2][1] -= cdgfac[face_id][1] * surfac[face_id][1];
      vol[cell_id2][2] -= cdgfac[face_id][2] * surfac[face_id][2];
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_lnum_t cell_id = b_face_cells[face_id];
      vol[cell_id][0] += cdgfbo[face_id][0] * surfbo[face_id][0];
      vol[cell_id][1] += cdgfbo[face_id][1] * surfbo[face_id][1];
      vol[cell_id][2] += cdgfbo[face_id][2] * surfbo[face_id][2];
    }

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      mq->bad_cell_indic[cell_id] = 0;

      double determinant = mq->corr_grad_lin_det[cell_id];

      int probleme = 0;
      if (determinant <= 0.) {
        probleme = 1;
      }
      else {
        determinant = CS_MAX(determinant, 1.e-10);
        determinant = CS_MAX(determinant, 1./determinant);
        if (determinant > 2)
          probleme = 1;
      }

      // FIXME not invariant by rotation
      if (probleme > 0 || (bad[cell_id] & CS_BAD_CELL_RATIO)
          || vol[cell_id][0] < 0. || vol[cell_id][1] < 0. || vol[cell_id][2] < 0.
          || (bad[cell_id] & CS_BAD_CELL_ORTHO_NORM)) {
        mq->bad_cell_indic[cell_id] = 1;
        nb_bad_cells ++;
      }

      f->val[cell_id] = mq->bad_cell_indic[cell_id];//TODO post
    }

    if (mesh->halo != NULL)
      cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, mq->bad_cell_indic);

    cs_parall_counter(&nb_bad_cells, 1);

    bft_printf("Number of bad cells with regularisation = %d\n", nb_bad_cells);

    BFT_FREE(vol);

    /* If no bad cells, no need of regularisation */
    if (nb_bad_cells <= 0) {
      bft_printf("No need of regularisation\n");
      cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_REGULARISATION;
    }
  }
  return;

}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for scalars
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_scalar(cs_real_t *var)
{


  const cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = mesh->n_cells;
  cs_lnum_t n_i_faces = mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const int *b_face_cells = mesh->b_face_cells;

  unsigned *bad = mq->bad_cell_flag;
  const cs_real_t *surfn = mq->i_face_surf;
  const cs_real_t *surfbn = mq->b_face_surf;
  double *dist = mq->i_dist;
  double *volume  = mq->cell_vol;


  /* Tag bad cells, if no bad cells, CS_BAD_CELLS_REGULARISATION is
   * switch off */
  _tag_bad_cells(mesh,
                 mq);

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_real_t *xam, *dam, *rhs;

  double varmin = 1.e20;
  double varmax =-1.e20;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    if (mq->bad_cell_indic[cell_id] == 0) {
      varmin = CS_MIN(varmin, var[cell_id]);
      varmax = CS_MAX(varmax, var[cell_id]);
    }

  cs_parall_min(1, CS_DOUBLE, &varmin);
  cs_parall_max(1, CS_DOUBLE, &varmax);

  BFT_MALLOC(xam, n_i_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_t);
  BFT_MALLOC(rhs, n_cells_ext, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    dam[cell_id] = 0.;
    rhs[cell_id] = 0.;
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];
    xam[face_id] = 0.;

    double surf = surfn[face_id];
    double vol = 0.5 * (volume[cell_id1] + volume[cell_id2]);
    surf = CS_MAX(surf, 0.1*vol/dist[face_id]);
    double ssd = surf / dist[face_id];

    dam[cell_id1] += ssd;
    dam[cell_id2] += ssd;

    if (mq->bad_cell_indic[cell_id1] > 0 && mq->bad_cell_indic[cell_id2] > 0) {
      xam[face_id] = -ssd;
    }
    else if (mq->bad_cell_indic[cell_id1] > 0) {
      rhs[cell_id1] += ssd * var[cell_id2];
      rhs[cell_id2] += ssd * var[cell_id2];
    }
    else if (mq->bad_cell_indic[cell_id2] > 0) {
      rhs[cell_id2] += ssd * var[cell_id1];
      rhs[cell_id1] += ssd * var[cell_id1];
    }
    else {
      rhs[cell_id1] += ssd * var[cell_id1];
      rhs[cell_id2] += ssd * var[cell_id2];
    }
  }

  cs_real_t rnorm = sqrt(cs_gdot(n_cells, rhs, rhs));

  /*  Solver residual */
  double ressol = 0.;
  int niterf = 0;

  int nitmax = 10000;
  cs_real_t epsilp = 1.e-12;

  /* Matrix block size */
  int ibsize = 1;
  int db_size[4] = {ibsize, ibsize, ibsize, ibsize*ibsize};

  cs_sles_solve_native(-1, /* f_id */
                       "potential_regularisation_scalar",
                       true, /* symmetric */
                       db_size,
                       NULL, /* eb_size */
                       (cs_real_t *)dam,
                       xam,
                       CS_HALO_ROTATION_COPY,
                       epsilp,
                       rnorm,
                       &niterf,
                       &ressol,
                       (cs_real_t *)rhs,
                       (cs_real_t *)var);

  bft_printf("Solving %s: N iter: %d, Res: %12.5e, Norm: %12.5e\n",
                 "potential_regularisation_scalar", niterf, ressol, rnorm);
  //FIXME use less!
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    var[cell_id] = CS_MIN(var[cell_id], varmax);
    var[cell_id] = CS_MAX(var[cell_id], varmin);
  }

  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, var);

  /* Free solver setup */
  cs_sles_free_native(-1, /* f_id*/
                      "potential_regularisation_scalar");

  /* Free memory */
  BFT_FREE(xam);
  BFT_FREE(dam);
  BFT_FREE(rhs);

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for vectors
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_vector(cs_real_3_t *var,
                                   int         boundary_projection)
{


  const cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = mesh->n_cells;
  cs_lnum_t n_i_faces = mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const int *b_face_cells = mesh->b_face_cells;

  unsigned *bad = mq->bad_cell_flag;//FIXME use it
  const cs_real_t *surfn = mq->i_face_surf;
  const cs_real_t *surfbn = mq->b_face_surf;
  double *dist = mq->i_dist;
  double *distbr = mq->b_dist;
  double *volume  = mq->cell_vol;

  const cs_real_3_t *surfbo = (const cs_real_3_t *) mq->b_face_normal;

  /* Tag bad cells, if no bad cells, CS_BAD_CELLS_REGULARISATION is
   * switch off */
  _tag_bad_cells(mesh,
                 mq);

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_real_33_t *dam;
  cs_real_3_t *rhs;
  cs_real_t *xam;
#if 1
  double varmin[3] = {1.e20, 1.e20, 1.e20};
  double varmax[3] = {-1.e20, -1.e20,-1.e20};

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    if (mq->bad_cell_indic[cell_id] == 0) {
      for (int i = 0; i < 3; i++) {
        varmin[i] = CS_MIN(varmin[i], var[cell_id][i]);
        varmax[i] = CS_MAX(varmax[i], var[cell_id][i]);
      }
    }

  for (int i = 0; i < 3; i++) {
    cs_parall_min(1, CS_DOUBLE, &varmin[i]);
    cs_parall_max(1, CS_DOUBLE, &varmax[i]);
  }
#endif

  BFT_MALLOC(xam, n_i_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_33_t);
  BFT_MALLOC(rhs, n_cells_ext, cs_real_3_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        dam[cell_id][i][j] = 0.;
      }
      rhs[cell_id][i] = 0.;
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];
    xam[face_id] = 0.;

    //FIXME usefull?
    double surf = surfn[face_id];
    double vol = 0.5 * (volume[cell_id1] + volume[cell_id2]);
    surf = CS_MAX(surf, 0.1*vol/dist[face_id]);
    double ssd = surf / dist[face_id];

    for (int i = 0; i < 3; i++) {
      dam[cell_id1][i][i] += ssd;
      dam[cell_id2][i][i] += ssd;
    }

    if (mq->bad_cell_indic[cell_id1] > 0 && mq->bad_cell_indic[cell_id2] > 0) {
      xam[face_id] = -ssd;
    }
    else if (mq->bad_cell_indic[cell_id1] > 0) {
      for (int i = 0; i < 3; i++) {
        rhs[cell_id1][i] += ssd * var[cell_id2][i];
        rhs[cell_id2][i] += ssd * var[cell_id2][i];
      }
    }
    else if (mq->bad_cell_indic[cell_id2] > 0) {
      for (int i = 0; i < 3; i++) {
        rhs[cell_id2][i] += ssd * var[cell_id1][i];
        rhs[cell_id1][i] += ssd * var[cell_id1][i];
      }
    }
    else {
      for (int i = 0; i < 3; i++) {
        rhs[cell_id1][i] += ssd * var[cell_id1][i];
        rhs[cell_id2][i] += ssd * var[cell_id2][i];
      }
    }
  }

  /* Boudanry projection... should be consistent with BCs... */
  if (boundary_projection == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (cs_glob_bc_type[face_id] == CS_SMOOTHWALL ||
          cs_glob_bc_type[face_id] == CS_ROUGHWALL  ||
          cs_glob_bc_type[face_id] == CS_SYMMETRY     ) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        if (mq->bad_cell_indic[cell_id] > 0) {
          double ssd = surfbn[face_id] / distbr[face_id];
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              double nn = surfbo[face_id][i]/surfbn[face_id] * surfbo[face_id][j]/surfbn[face_id];
              dam[cell_id][i][j] += ssd * nn;
            }
          }
        }
      }
    }
  }

  cs_real_t rnorm = sqrt(cs_gdot(3*n_cells, rhs, rhs ));

  /*  Solver residual */
  double ressol = 0.;
  int niterf = 0;

  int nitmax = 10000;
  cs_real_t epsilp = 1.e-12;

  /* Matrix block size */
  int ibsize = 3;
  int db_size[4] = {ibsize, ibsize, ibsize, ibsize*ibsize};

  cs_sles_solve_native(-1, /* f_id */
                       "potential_regularisation_vector",
                       true, /* symmetric */
                       db_size,
                       NULL, /* eb_size */
                       (cs_real_t *)dam,
                       xam,
                       CS_HALO_ROTATION_COPY,
                       epsilp,
                       rnorm,
                       &niterf,
                       &ressol,
                       (cs_real_t *)rhs,
                       (cs_real_t *)var);

  bft_printf("Solving %s: N iter: %d, Res: %12.5e, Norm: %12.5e\n",
                 "potential_regularisation_vector", niterf, ressol, rnorm);

  //FIXME useless clipping? the matrix is min/max preserving..
#if 1
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 3; i++) {
      var[cell_id][i] = CS_MIN(var[cell_id][i], varmax[i]);
      var[cell_id][i] = CS_MAX(var[cell_id][i], varmin[i]);
    }
  }
#endif

  //FIXME periodicity of rotation
  if (mesh->halo != NULL)
    cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD, (cs_real_t *)var, 3);

  /* Free solver setup */
  cs_sles_free_native(-1, /* f_id*/
                      "potential_regularisation_vector");

  /* Free memory */
  BFT_FREE(xam);
  BFT_FREE(dam);
  BFT_FREE(rhs);

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for symmetric tensors
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_sym_tensor(cs_real_6_t *var,
                                       int         boundary_projection)
{

  const cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = mesh->n_cells;
  cs_lnum_t n_i_faces = mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const int *b_face_cells = mesh->b_face_cells;

  unsigned *bad = mq->bad_cell_flag;
  const cs_real_t *surfn = mq->i_face_surf;
  const cs_real_t *surfbn = mq->b_face_surf;
  double *dist = mq->i_dist;
  double *distbr = mq->b_dist;
  double *volume  = mq->cell_vol;

  const cs_real_3_t *surfbo = (const cs_real_3_t *) mq->b_face_normal;

  /* Tag bad cells, if no bad cells, CS_BAD_CELLS_REGULARISATION is
   * switch off */
  _tag_bad_cells(mesh,
                 mq);

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_real_66_t *dam;
  cs_real_6_t *rhs;
  cs_real_t *xam;
#if 1
  double varmin[6] = { 1.e20,  1.e20, 1.e20,  1.e20,  1.e20, 1.e20};
  double varmax[6] = {-1.e20, -1.e20,-1.e20, -1.e20, -1.e20,-1.e20};

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    if (mq->bad_cell_indic[cell_id] == 0) {
      for (int i = 0; i < 6; i++) {
        varmin[i] = CS_MIN(varmin[i], var[cell_id][i]);
        varmax[i] = CS_MAX(varmax[i], var[cell_id][i]);
      }
    }

  for (int i = 0; i < 6; i++) {
    cs_parall_min(1, CS_DOUBLE, &varmin[i]);
    cs_parall_max(1, CS_DOUBLE, &varmax[i]);
  }
#endif

  BFT_MALLOC(xam, n_i_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_66_t);
  BFT_MALLOC(rhs, n_cells_ext, cs_real_6_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        dam[cell_id][i][j] = 0.;
      }
      rhs[cell_id][i] = 0.;
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];
    xam[face_id] = 0.;

    //FIXME usefull?
    double surf = surfn[face_id];
    double vol = 0.5 * (volume[cell_id1] + volume[cell_id2]);
    surf = CS_MAX(surf, 0.1*vol/dist[face_id]);
    double ssd = surf / dist[face_id];

    for (int i = 0; i < 6; i++) {
      dam[cell_id1][i][i] += ssd;
      dam[cell_id2][i][i] += ssd;
    }

    if (mq->bad_cell_indic[cell_id1] > 0 && mq->bad_cell_indic[cell_id2] > 0) {
      xam[face_id] = -ssd;
    }
    else if (mq->bad_cell_indic[cell_id1] > 0) {
      for (int i = 0; i < 6; i++) {
        rhs[cell_id1][i] += ssd * var[cell_id2][i];
        rhs[cell_id2][i] += ssd * var[cell_id2][i];
      }
    }
    else if (mq->bad_cell_indic[cell_id2] > 0) {
      for (int i = 0; i < 6; i++) {
        rhs[cell_id2][i] += ssd * var[cell_id1][i];
        rhs[cell_id1][i] += ssd * var[cell_id1][i];
      }
    }
    else {
      for (int i = 0; i < 6; i++) {
        rhs[cell_id1][i] += ssd * var[cell_id1][i];
        rhs[cell_id2][i] += ssd * var[cell_id2][i];
      }
    }
  }

  /* Boudanry projection... should be consistent with BCs... */
  if (boundary_projection == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (cs_glob_bc_type[face_id] == CS_SMOOTHWALL ||
          cs_glob_bc_type[face_id] == CS_ROUGHWALL  ||
          cs_glob_bc_type[face_id] == CS_SYMMETRY     ) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        if (mq->bad_cell_indic[cell_id] > 0) {
          double ssd = surfbn[face_id] / distbr[face_id];
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              double nn = surfbo[face_id][i]/surfbn[face_id] * surfbo[face_id][j]/surfbn[face_id];
//TODO ???              dam[cell_id][i][j] += ssd * nn;
            }
          }
        }
      }
    }
  }

  cs_real_t rnorm = sqrt(cs_gdot(6*n_cells, rhs, rhs));

  /*  Solver residual */
  double ressol = 0.;
  int niterf = 0;

  int nitmax = 10000;
  cs_real_t epsilp = 1.e-12;

  /* Matrix block size */
  int ibsize = 6;
  int db_size[4] = {ibsize, ibsize, ibsize, ibsize*ibsize};

  cs_sles_solve_native(-1, /* f_id */
                       "potential_regularisation_sym_tensor",
                       true, /* symmetric */
                       db_size,
                       NULL, /* eb_size */
                       (cs_real_t *)dam,
                       xam,
                       CS_HALO_ROTATION_COPY,
                       epsilp,
                       rnorm,
                       &niterf,
                       &ressol,
                       (cs_real_t *)rhs,
                       (cs_real_t *)var);

  bft_printf("Solving %s: N iter: %d, Res: %12.5e, Norm: %12.5e\n",
                 "potential_regularisation_sym_tensor", niterf, ressol, rnorm);
  //FIXME useless clipping? the matrix is min/max preserving..
#if 1
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 6; i++) {
      var[cell_id][i] = CS_MIN(var[cell_id][i], varmax[i]);
      var[cell_id][i] = CS_MAX(var[cell_id][i], varmin[i]);
    }
  }
#endif

  //FIXME periodicity of rotation
  if (mesh->halo != NULL)
    cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD, (cs_real_t *)var, 6);

  /* Free solver setup */
  cs_sles_free_native(-1, /* f_id*/
                      "potential_regularisation_sym_tensor");

  /* Free memory */
  BFT_FREE(xam);
  BFT_FREE(dam);
  BFT_FREE(rhs);

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Regularisation on bad cells for tensors
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_bad_cells_regularisation_tensor(cs_real_9_t *var,
                                   int         boundary_projection)
{

  const cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = mesh->n_cells;
  cs_lnum_t n_i_faces = mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const int *b_face_cells = mesh->b_face_cells;

  unsigned *bad = mq->bad_cell_flag;
  const cs_real_t *surfn = mq->i_face_surf;
  const cs_real_t *surfbn = mq->b_face_surf;
  double *dist = mq->i_dist;
  double *distbr = mq->b_dist;
  double *volume  = mq->cell_vol;

  const cs_real_3_t *surfbo = (const cs_real_3_t *) mq->b_face_normal;

  /* Tag bad cells, if no bad cells, CS_BAD_CELLS_REGULARISATION is
   * switch off */
  _tag_bad_cells(mesh,
                 mq);

  if (!(cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION))
    return;

  cs_real_99_t *dam;
  cs_real_9_t *rhs;
  cs_real_t *xam;
#if 1
  double varmin[9] = { 1.e20,  1.e20, 1.e20,  1.e20,  1.e20, 1.e20,  1.e20,  1.e20, 1.e20};
  double varmax[9] = {-1.e20, -1.e20,-1.e20, -1.e20, -1.e20,-1.e20, -1.e20, -1.e20,-1.e20};

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    if (mq->bad_cell_indic[cell_id] == 0) {
      for (int i = 0; i < 9; i++) {
        varmin[i] = CS_MIN(varmin[i], var[cell_id][i]);
        varmax[i] = CS_MAX(varmax[i], var[cell_id][i]);
      }
    }

  for (int i = 0; i < 9; i++) {
    cs_parall_min(1, CS_DOUBLE, &varmin[i]);
    cs_parall_max(1, CS_DOUBLE, &varmax[i]);
  }
#endif

  BFT_MALLOC(xam, n_i_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_99_t);
  BFT_MALLOC(rhs, n_cells_ext, cs_real_9_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        dam[cell_id][i][j] = 0.;
      }
      rhs[cell_id][i] = 0.;
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];
    xam[face_id] = 0.;

    //FIXME usefull?
    double surf = surfn[face_id];
    double vol = 0.5 * (volume[cell_id1] + volume[cell_id2]);
    surf = CS_MAX(surf, 0.1*vol/dist[face_id]);
    double ssd = surf / dist[face_id];

    for (int i = 0; i < 9; i++) {
      dam[cell_id1][i][i] += ssd;
      dam[cell_id2][i][i] += ssd;
    }

    if (mq->bad_cell_indic[cell_id1] > 0 && mq->bad_cell_indic[cell_id2] > 0) {
      xam[face_id] = -ssd;
    }
    else if (mq->bad_cell_indic[cell_id1] > 0) {
      for (int i = 0; i < 9; i++) {
        rhs[cell_id1][i] += ssd * var[cell_id2][i];
        rhs[cell_id2][i] += ssd * var[cell_id2][i];
      }
    }
    else if (mq->bad_cell_indic[cell_id2] > 0) {
      for (int i = 0; i < 9; i++) {
        rhs[cell_id2][i] += ssd * var[cell_id1][i];
        rhs[cell_id1][i] += ssd * var[cell_id1][i];
      }
    }
    else {
      for (int i = 0; i < 9; i++) {
        rhs[cell_id1][i] += ssd * var[cell_id1][i];
        rhs[cell_id2][i] += ssd * var[cell_id2][i];
      }
    }
  }

  /* Boudanry projection... should be consistent with BCs... */
  if (boundary_projection == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (cs_glob_bc_type[face_id] == CS_SMOOTHWALL ||
          cs_glob_bc_type[face_id] == CS_ROUGHWALL  ||
          cs_glob_bc_type[face_id] == CS_SYMMETRY     ) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        if (mq->bad_cell_indic[cell_id] > 0) {
          double ssd = surfbn[face_id] / distbr[face_id];
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              double nn = surfbo[face_id][i]/surfbn[face_id] * surfbo[face_id][j]/surfbn[face_id];
//TODO ???              dam[cell_id][i][j] += ssd * nn;
            }
          }
        }
      }
    }
  }

  cs_real_t rnorm = sqrt(cs_gdot(9*n_cells, rhs, rhs));

  /*  Solver residual */
  double ressol = 0.;
  int niterf = 0;

  int nitmax = 10000;
  cs_real_t epsilp = 1.e-12;

  /* Matrix block size */
  int ibsize = 9;
  int db_size[4] = {ibsize, ibsize, ibsize, ibsize*ibsize};

  cs_sles_solve_native(-1, /* f_id */
                       "potential_regularisation_tensor",
                       true, /* symmetric */
                       db_size,
                       NULL, /* eb_size */
                       (cs_real_t *)dam,
                       xam,
                       CS_HALO_ROTATION_COPY,
                       epsilp,
                       rnorm,
                       &niterf,
                       &ressol,
                       (cs_real_t *)rhs,
                       (cs_real_t *)var);

  bft_printf("Solving %s: N iter: %d, Res: %12.5e, Norm: %12.5e\n",
                 "potential_regularisation_tensor", niterf, ressol, rnorm);
  //FIXME useless clipping? the matrix is min/max preserving..
#if 1
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 9; i++) {
      var[cell_id][i] = CS_MIN(var[cell_id][i], varmax[i]);
      var[cell_id][i] = CS_MAX(var[cell_id][i], varmin[i]);
    }
  }
#endif

  //FIXME periodicity of rotation
  if (mesh->halo != NULL)
    cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD, (cs_real_t *)var, 9);

  /* Free solver setup */
  cs_sles_free_native(-1, /* f_id*/
                      "potential_regularisation_tensor");

  /* Free memory */
  BFT_FREE(xam);
  BFT_FREE(dam);
  BFT_FREE(rhs);

  return;
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
