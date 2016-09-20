/*============================================================================
 * Internal coupling: coupling for one instance of Code_Saturne
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"
#include "bft_error.h"

#include "fvm_defs.h"
#include "fvm_selector.h"

#include "cs_defs.h"
#include "cs_math.h"
#include "cs_sort.h"
#include "cs_search.h"
#include "cs_mesh_connect.h"
#include "cs_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_convection_diffusion.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_internal_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_internal_coupling_t  *_internal_coupling = NULL;
static int                      _n_internal_couplings = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the locator associated with given coupling entity and group number
 *
 * parameters:
 *   cpl  <-- pointer to coupling structure
 *   pov  <-- group number
 *----------------------------------------------------------------------------*/

static ple_locator_t*
_cs_internal_coupling_create_locator(cs_internal_coupling_t  *cpl,
                                     int                      pov)
{
  cs_lnum_t ifac;
  int i, j, n_local = 0, n_distant = 0;

  cs_lnum_t *faces_local = NULL;
  cs_lnum_t *faces_distant = NULL;
  fvm_nodal_t* nm = NULL;

  ple_coord_t* point_coords = NULL;

  char mesh_name[16];
  sprintf(mesh_name, "locator_%i", pov);


#if defined(PLE_HAVE_MPI)
  ple_locator_t *locator = ple_locator_create(cs_glob_mpi_comm,
                                              cs_glob_n_ranks,
                                              0);
#else
  ple_locator_t *locator = ple_locator_create();
#endif

  if (pov == 2) {
    n_local = cpl->n_1;
    faces_local = cpl->faces_1;

    n_distant = cpl->n_2;
    faces_distant = cpl->faces_2;
  }
  else if (pov == 1) {
    n_local = cpl->n_2;
    faces_local = cpl->faces_2;

    n_distant = cpl->n_1;
    faces_distant = cpl->faces_1;
  }
  else {
    bft_error(__FILE__, __LINE__, 0, "Wrong selector (pov)\n");
  }

  nm = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                      mesh_name,
                                      false,
                                      0,
                                      n_local,
                                      NULL,
                                      faces_local); /* 1..n */

  /* Creation of distant group cell centers */

  BFT_MALLOC(point_coords, 3*n_distant, cs_real_t);


  for (i = 0; i < n_distant; i++) {
    ifac = faces_distant[i] - 1; /* 0..n-1 */
    for (j = 0; j < 3; j++) {
      point_coords[3*i+j] = cs_glob_mesh_quantities->b_face_cog[3*ifac+j];
    }
  }

  /* Locator initialization */
  ple_locator_set_mesh(locator,
                       nm,
                       NULL,
                       0,
                       1.1, /* TODO */
                       cpl->dim,
                       n_distant,
                       NULL,
                       NULL,
                       point_coords,
                       NULL,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  /* Free memory */
  nm = fvm_nodal_destroy(nm);
  BFT_FREE(point_coords);
  return locator;
}

/*----------------------------------------------------------------------------
 * Destruction of given internal coupling structure.
 *
 * parameters:
 *   cpl <-> pointer to coupling structure to destroy
 *----------------------------------------------------------------------------*/

static void
_cs_internal_coupling_destroy_entity(cs_internal_coupling_t  *cpl)
{
  BFT_FREE(cpl->faces_1);
  BFT_FREE(cpl->faces_2);

  BFT_FREE(cpl->hint_1);
  BFT_FREE(cpl->hint_2);

  BFT_FREE(cpl->hext_1);
  BFT_FREE(cpl->hext_2);

  BFT_FREE(cpl->gweight_1);
  BFT_FREE(cpl->gweight_2);

  BFT_FREE(cpl->ij_1);
  BFT_FREE(cpl->ij_2);

  BFT_FREE(cpl->ofij_1);
  BFT_FREE(cpl->ofij_2);

  BFT_FREE(cpl->coupled_faces);

  BFT_FREE(cpl->cocgb_s_lsq);
  BFT_FREE(cpl->cocgb_s_it);
  BFT_FREE(cpl->cocg_s_it);

  BFT_FREE(cpl->namesca);

  ple_locator_destroy(cpl->locator_1);
  ple_locator_destroy(cpl->locator_2);
}

/*----------------------------------------------------------------------------
 * Compare local numbers (qsort function).
 *
 * parameters:
 *   a <-- pointer to first value
 *   b <-- pointer to second value
 *
 * returns:
 *   a-b
 *----------------------------------------------------------------------------*/

static int
_compare_int(const void* a,
             const void* b)
{
  return (*(const int*)a - *(const int*)b);
}

/*----------------------------------------------------------------------------
 * Compute weights around coupling interface.
 *
 * parameters:
 *   cpl <-- pointer to coupling structure
 *----------------------------------------------------------------------------*/

static void
_cs_internal_coupling_exchange_gweight(const cs_internal_coupling_t  *cpl)
{
  int ii;
  cs_lnum_t face_id, cell_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_lnum_t n_dist_1, n_dist_2;
  cs_lnum_t *dist_loc_1, *dist_loc_2;

  cs_real_3_t *ij_1, *ij_2;

  cs_real_t* gweight_1 = cpl->gweight_1;
  cs_real_t* gweight_2 = cpl->gweight_2;
  cs_real_t *gweight_distant_1, *gweight_distant_2;

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_mesh_t             *m   = cs_glob_mesh;

  const cs_real_t* cell_cen = fvq->cell_cen;
  const cs_real_t* b_face_cog = fvq->b_face_cog;
  const cs_real_t* diipb = fvq->diipb;
  const cs_real_t* b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_t* b_face_normal = cs_glob_mesh_quantities->b_face_normal;

  cs_internal_coupling_coupled_faces(cpl,
             &n_1,
             &n_2,
             &faces_1,
             &faces_2,
             &n_dist_1,
             &n_dist_2,
             &dist_loc_1,
             &dist_loc_2);

  ij_1 = cpl->ij_1;
  ij_2 = cpl->ij_2;

  /* Store local FI' distance in gweight_distant_* */
  BFT_MALLOC(gweight_distant_1, n_dist_1, cs_real_t);
  for (ii = 0; ii < n_dist_1; ii++) {
    cs_real_t dv[3];
    face_id = dist_loc_1[ii] - 1;
    cell_id = m->b_face_cells[face_id];

    for (int jj = 0; jj < 3; jj++)
      dv[jj] =  - diipb[3*face_id + jj] - cell_cen[3*cell_id +jj]
                + b_face_cog[3*face_id +jj];

    gweight_distant_1[ii] = cs_math_3_norm(dv);
  }

  BFT_MALLOC(gweight_distant_2, n_dist_2, cs_real_t);
  for (ii = 0; ii < n_dist_2; ii++) {
    cs_real_t dv[3];
    face_id = dist_loc_2[ii] - 1;
    cell_id = m->b_face_cells[face_id];

    for (int jj = 0; jj < 3; jj++)
      dv[jj] =  - diipb[3*face_id + jj] - cell_cen[3*cell_id +jj]
                + b_face_cog[3*face_id +jj];

    gweight_distant_2[ii] = cs_math_3_norm(dv);
  }

  /* Groups 1 and 2 exchange FI' distances */
  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    gweight_distant_1,
                                    gweight_distant_2,
                                    gweight_1,
                                    gweight_2);

  BFT_FREE(gweight_distant_1);
  BFT_FREE(gweight_distant_2);

  /* Normalise the distance to obtain weights */
  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    gweight_1[ii] /= ( ij_1[ii][0]*b_face_normal[3*face_id+0]
                     + ij_1[ii][1]*b_face_normal[3*face_id+1]
                     + ij_1[ii][2]*b_face_normal[3*face_id+2])
                     / b_face_surf[face_id];
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    gweight_2[ii] /= ( ij_2[ii][0]*b_face_normal[3*face_id+0]
                     + ij_2[ii][1]*b_face_normal[3*face_id+1]
                     + ij_2[ii][2]*b_face_normal[3*face_id+2])
                     / b_face_surf[face_id];
  }
}

/*----------------------------------------------------------------------------
 * Compute rweight_* around coupling interface based on diffusivity c_weight.
 *
 * parameters:
 *   cpl          <-- pointer to coupling structure
 *   c_weight[]   <-- diffusivity
 *   rweight_1[]  -> rhs weight for group 1
 *   rweight_2[]  -> rhs weight for group 2
 *----------------------------------------------------------------------------*/

static void
_cs_internal_coupling_exchange_rhs_weight(const cs_internal_coupling_t  *cpl,
                                          const cs_real_t        c_weight[],
                                          cs_real_t              rweight_1[],
                                          cs_real_t              rweight_2[])
{
  cs_real_t ki, kj, pond;

  const cs_real_t* gweight_1 = cpl->gweight_1;
  const cs_real_t* gweight_2 = cpl->gweight_2;

  cs_real_t *kj_1, *kj_2;

  int ii;
  cs_lnum_t face_id, cell_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  BFT_MALLOC(kj_1, n_1, cs_real_t);
  BFT_MALLOC(kj_2, n_2, cs_real_t);

  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           1,
                                           c_weight,
                                           kj_1,
                                           kj_2);


  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    ki = c_weight[cell_id];
    kj = kj_1[ii];
    pond = gweight_1[ii];

    rweight_1[ii] = kj / ( pond * ki + (1. - pond) * kj);
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];
    ki = c_weight[cell_id];
    kj = kj_2[ii];
    pond = gweight_2[ii];

    rweight_2[ii] = kj / ( pond * ki + (1. - pond) * kj);
  }

  BFT_FREE(kj_1);
  BFT_FREE(kj_2);
}

/*----------------------------------------------------------------------------
 * Compute vector OFIJ on coupled faces
 *
 * parameters:
 *   cpl <-> pointer to coupling entity
 *----------------------------------------------------------------------------*/

static void
_cs_internal_coupling_ofij(const cs_internal_coupling_t  *cpl)
{
  int ii, jj;
  cs_lnum_t face_id, cell_id;

  int *faces_1, *faces_2;
  cs_lnum_t n_1, n_2;

  cs_real_3_t* ofij_1 = cpl->ofij_1;
  cs_real_3_t* ofij_2 = cpl->ofij_2;

  const cs_real_t* gweight_1 = cpl->gweight_1;
  const cs_real_t* gweight_2 = cpl->gweight_2;

  const cs_mesh_quantities_t* mq = cs_glob_mesh_quantities;
  const cs_mesh_t* m = cs_glob_mesh;
  const cs_real_t* cell_cen = mq->cell_cen;
  const cs_real_t* b_face_cog = mq->b_face_cog;
  const cs_lnum_t* b_face_cells = m->b_face_cells;

  cs_real_t *cell_cen_1, *cell_cen_2;
  cs_real_t xxd, xxl;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  BFT_MALLOC(cell_cen_1, 3 * n_1, cs_real_t);
  BFT_MALLOC(cell_cen_2, 3 * n_2, cs_real_t);

  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           3,
                                           mq->cell_cen,
                                           cell_cen_1,
                                           cell_cen_2);

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];

    for (jj = 0; jj < 3; jj++) {
      xxd = cell_cen_1[3*ii + jj];
      xxl = cell_cen[3*cell_id + jj];
      ofij_1[ii][jj] = b_face_cog[3*face_id + jj]
        - (        gweight_1[ii]*xxl
           + (1. - gweight_1[ii])*xxd);
    }
  }
  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    for (jj = 0; jj < 3; jj++) {
      xxd = cell_cen_2[3*ii + jj];
      xxl = cell_cen[3*cell_id + jj];
      ofij_2[ii][jj] = b_face_cog[3*face_id + jj]
        - (        gweight_2[ii]*xxl
           + (1. - gweight_2[ii])*xxd);
    }
  }

  BFT_FREE(cell_cen_1);
  BFT_FREE(cell_cen_2);
}

/*----------------------------------------------------------------------------
 * Update local variable using distant one.
 *
 * parameters:
 *   cpl <-- pointer to coupling structure
 *   stride          <-- number of values (non interlaced) by entity
 *   pov             <-- number of values (non interlaced) by entity
 *   local[]         --> local variable updated
 *   distant_var[]   <-- distant variable
 *----------------------------------------------------------------------------*/

static void
_cs_internal_coupling_send_var(const cs_internal_coupling_t  *cpl,
                               int                            stride,
                               int                            pov,
                               cs_real_t                      local[],
                               cs_real_t                      distant_var[])
{
  ple_locator_t *locator = NULL;

  switch(pov) {
  case 1:
    locator = cpl->locator_1;
    break;
  case 2:
    locator = cpl->locator_2;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0, "Wrong selector (pov=%d)", pov);
  }

  ple_locator_exchange_point_var(locator,
                                 distant_var,
                                 local,
                                 NULL,
                                 sizeof(cs_real_t),
                                 stride,
                                 0);
}

/*----------------------------------------------------------------------------
 * Define component coupled_faces[] of given coupling entity.
 *
 * parameters:
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

static void
_cs_internal_coupling_select_coupled_faces(
    cs_internal_coupling_t *cpl)

{
  int ii;
  cs_lnum_t face_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  const cs_mesh_t* m = cs_glob_mesh;
  bool *facoup = cpl->coupled_faces;
  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  for (face_id = 0; face_id < m-> n_b_faces; face_id++) {
    facoup[face_id] = false;
  }

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    facoup[face_id] = true;
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    facoup[face_id] = true;
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add contribution from coupled faces (internal coupling) to initialisation
 * for iterative scalar gradient calculation
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
 *   c_weight <-- weighted gradient coefficient variable, or NULL
 *   pvar     <-- variable
 *   grad     <-> gradient
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_initial_contribution(const cs_internal_coupling_t  *cpl,
                                          const cs_real_t                c_weight[],
                                          const cs_real_t                pvar[],
                                          cs_real_3_t          *restrict grad)
{
  int ii;
  cs_lnum_t face_id, cell_id;

  const cs_real_t* gweight_1 = cpl->gweight_1;
  const cs_real_t* gweight_2 = cpl->gweight_2;
  cs_real_t *rweight_1 = NULL, *rweight_2 = NULL; /* Weight with c_weight */

  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_lnum_t n_dist_1, n_dist_2;
  cs_lnum_t *dist_loc_1, *dist_loc_2;

  cs_real_t *pvar_local_1, *pvar_local_2;
  cs_real_t *pvar_distant_1, *pvar_distant_2;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_mesh_quantities_t* fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     &n_dist_1,
                                     &n_dist_2,
                                     &dist_loc_1,
                                     &dist_loc_2);

  /* 1 & 2 exchange pvar */
  BFT_MALLOC(pvar_distant_1, n_dist_1, cs_real_t);
  BFT_MALLOC(pvar_distant_2, n_dist_2, cs_real_t);
  for (ii = 0; ii < n_dist_1; ii++) {
    face_id = dist_loc_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    pvar_distant_1[ii] = pvar[cell_id];
  }
  for (ii = 0; ii < n_dist_2; ii++) {
    face_id = dist_loc_2[ii] - 1;
    cell_id = b_face_cells[face_id];
    pvar_distant_2[ii] = pvar[cell_id];
  }
  BFT_MALLOC(pvar_local_1, n_1, cs_real_t);
  BFT_MALLOC(pvar_local_2, n_2, cs_real_t);
  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    pvar_distant_1,
                                    pvar_distant_2,
                                    pvar_local_1,
                                    pvar_local_2);
  BFT_FREE(pvar_distant_1);
  BFT_FREE(pvar_distant_2);

  /* Add contribution */
  if (c_weight != NULL) { /* Heterogenous diffusivity */
    BFT_MALLOC(rweight_1, n_1, cs_real_t);
    BFT_MALLOC(rweight_2, n_2, cs_real_t);
    _cs_internal_coupling_exchange_rhs_weight(cpl,
                                              c_weight, /* diffusivity */
                                              rweight_1, /* rhs weight */
                                              rweight_2); /* rhs weight */
    /* Redefinition of rweight_* :
         Before : (1-gweight_1)*rweight_1 <==> 1 - ktpond
                  (1-gweight_1)*rweight_1 + (1-gweight_2)*rweight_2 = 1
         Modif : rweight_1 = ktpond
                 rweight_1 + rweight_2 = 1
       Scope of this modification is local */
    for (ii = 0; ii < n_1; ii++)
      rweight_1[ii] = 1.0 - (1.0-gweight_1[ii]) * rweight_1[ii];
    for (ii = 0; ii < n_2; ii++)
      rweight_2[ii] = 1.0 - (1.0-gweight_2[ii]) * rweight_2[ii];
  }

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];

    /* compared to _initialize_scalar_gradient :
       1 - rweight_1 <==> 1 - ktpond
       gweight_1 <==> weight = alpha_ij
       pvar[cell_id] <==> pvar[ii]
       pvar_local_1[ii] <==> pvar[jj]
       b_f_face_normal <==> i_f_face_normal */
    cs_real_t pfaci = (c_weight == NULL) ?
      (1.0-gweight_1[ii]) * (pvar_local_1[ii] - pvar[cell_id]) :
      (1.0-rweight_1[ii]) * (pvar_local_1[ii] - pvar[cell_id]);

    for (int j = 0; j < 3; j++) {
      grad[cell_id][j] += pfaci * b_f_face_normal[face_id][j];
    }

  }
  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    /* compared to _initialize_scalar_gradient :
       1 - rweight_2  = rweight_1 <==> ktpond
       gweight_1 <==> weight = alpha_ij
       pvar[cell_id] <==> pvar[jj]
       pvar_local_2[ii] <==> pvar[ii]
       b_f_face_normal <==> -i_f_face_normal */
    cs_real_t pfacj = (c_weight == NULL) ?
      (1.0-gweight_2[ii]) * (pvar_local_2[ii] - pvar[cell_id]) :
      (1.0-rweight_2[ii]) * (pvar_local_2[ii] - pvar[cell_id]);

    for (int j = 0; j < 3; j++) {
      grad[cell_id][j] += pfacj * b_f_face_normal[face_id][j];
    }

  }

  if (c_weight != NULL) {
    BFT_FREE(rweight_1);
    BFT_FREE(rweight_2);
  }
  BFT_FREE(pvar_local_1);
  BFT_FREE(pvar_local_2);
}

/*----------------------------------------------------------------------------
 * Add internal coupling rhs contribution for iterative gradient calculation
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
 *   c_weight <-- weighted gradient coefficient variable, or NULL
 *   grad     <-- pointer to gradient
 *   pvar     <-- pointer to variable
 *   rhs      <-> pointer to rhs contribution
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_iter_rhs(const cs_internal_coupling_t  *cpl,
                              const cs_real_t                c_weight[],
                              cs_real_3_t          *restrict grad,
                              const cs_real_t                pvar[],
                              cs_real_3_t                    rhs[])
{
  int ii, ll;
  cs_lnum_t face_id, cell_id;

  const cs_real_t* gweight_1 = cpl->gweight_1;
  const cs_real_t* gweight_2 = cpl->gweight_2;
  cs_real_t *rweight_1 = NULL, *rweight_2 = NULL;

  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_lnum_t n_dist_1, n_dist_2;
  cs_lnum_t *dist_loc_1, *dist_loc_2;
  cs_real_3_t *ofij_1, *ofij_2;

  cs_real_t *grad_local_1, *grad_local_2;
  cs_real_t *grad_distant_1, *grad_distant_2;
  cs_real_t *pvar_local_1, *pvar_local_2;
  cs_real_t *pvar_distant_1, *pvar_distant_2;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_mesh_quantities_t* fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     &n_dist_1,
                                     &n_dist_2,
                                     &dist_loc_1,
                                     &dist_loc_2);

  ofij_1 = cpl->ofij_1;
  ofij_2 = cpl->ofij_2;

  /* 1 & 2 exchange grad and pvar */
  BFT_MALLOC(grad_distant_1, 3*n_dist_1, cs_real_t);
  BFT_MALLOC(grad_distant_2, 3*n_dist_2, cs_real_t);
  BFT_MALLOC(pvar_distant_1, n_dist_1, cs_real_t);
  BFT_MALLOC(pvar_distant_2, n_dist_2, cs_real_t);

  for (ii = 0; ii < n_dist_1; ii++) {
    face_id = dist_loc_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    pvar_distant_1[ii] = pvar[cell_id];
    for (ll = 0; ll < 3; ll++) {
      grad_distant_1[3*ii+ll] = grad[cell_id][ll];
    }
  }

  for (ii = 0; ii < n_dist_2; ii++) {
    face_id = dist_loc_2[ii] - 1;
    cell_id = b_face_cells[face_id];
    pvar_distant_2[ii] = pvar[cell_id];
    for (ll = 0; ll < 3; ll++) {
      grad_distant_2[3*ii+ll] = grad[cell_id][ll];
    }
  }

  BFT_MALLOC(grad_local_1, 3*n_1, cs_real_t);
  BFT_MALLOC(grad_local_2, 3*n_2, cs_real_t);
  BFT_MALLOC(pvar_local_1, n_1, cs_real_t);
  BFT_MALLOC(pvar_local_2, n_2, cs_real_t);

  cs_internal_coupling_exchange_var(cpl,
                                    3,
                                    grad_distant_1,
                                    grad_distant_2,
                                    grad_local_1,
                                    grad_local_2);
  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    pvar_distant_1,
                                    pvar_distant_2,
                                    pvar_local_1,
                                    pvar_local_2);

  BFT_FREE(grad_distant_1);
  BFT_FREE(grad_distant_2);
  BFT_FREE(pvar_distant_1);
  BFT_FREE(pvar_distant_2);

  /* Compute rhs */
  if (c_weight != NULL) { /* Heterogenous diffusivity */
    BFT_MALLOC(rweight_1, n_1, cs_real_t);
    BFT_MALLOC(rweight_2, n_2, cs_real_t);
    _cs_internal_coupling_exchange_rhs_weight(cpl,
                                              c_weight, /* diffusivity */
                                              rweight_1, /* rhs weight */
                                              rweight_2); /* rhs weight */
    /* Redefinition of rweight_* :
         Before : (1-gweight_1)*rweight_1 <==> 1 - ktpond
                  (1-gweight_1)*rweight_1 + (1-gweight_2)*rweight_2 = 1
         Modif : rweight_1 = ktpond
                 rweight_1 + rweight_2 = 1
       Scope of this modification is local */
    for (ii = 0; ii < n_1; ii++)
      rweight_1[ii] = 1.0 - (1.0-gweight_1[ii]) * rweight_1[ii];
    for (ii = 0; ii < n_2; ii++)
      rweight_2[ii] = 1.0 - (1.0-gweight_2[ii]) * rweight_2[ii];
  }

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];

    /*
       Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                + (1-\alpha_\ij) \varia_\cellj\f$
               but for the cell \f$ \celli \f$ we remove
               \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
               and for the cell \f$ \cellj \f$ we remove
               \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
    */

    /* Reconstruction part
         compared to _iterative_scalar_gradient :
         gweight_1 <==> weight = alpha_ij
         pvar[cell_id] <==> pvar[cell_id1]
         pvar_local_1[ii] <==> pvar[cell_id2]
         b_f_face_normal <==> i_f_face_normal */
    cs_real_t pfaci = 0.5;
    pfaci *= ofij_1[ii][0]*(grad_local_1[3*ii  ]+grad[cell_id][0])
            +ofij_1[ii][1]*(grad_local_1[3*ii+1]+grad[cell_id][1])
            +ofij_1[ii][2]*(grad_local_1[3*ii+2]+grad[cell_id][2]);
    if (c_weight != NULL) {
      pfaci += (1.0-rweight_1[ii]) * (pvar_local_1[ii] - pvar[cell_id]);
    } else {
      pfaci += (1.0-gweight_1[ii]) * (pvar_local_1[ii] - pvar[cell_id]);
    }

    for (int j = 0; j < 3; j++) {
      rhs[cell_id][j] += pfaci * b_f_face_normal[face_id][j];
    }

  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    /*
       Remark: \f$ \varia_\face = \alpha_\ij \varia_\celli
                                + (1-\alpha_\ij) \varia_\cellj\f$
               but for the cell \f$ \celli \f$ we remove
               \f$ \varia_\celli \sum_\face \vect{S}_\face = \vect{0} \f$
               and for the cell \f$ \cellj \f$ we remove
               \f$ \varia_\cellj \sum_\face \vect{S}_\face = \vect{0} \f$
    */

    /* Reconstruction part
         compared to _iterative_scalar_gradient :
         1-gweight_2 = gweight_1 <==> weight = alpha_ij
         pvar[cell_id] <==> pvar[cell_id2]
         pvar_local_2[ii] <==> pvar[cell_id1]
         b_f_face_normal <==> -i_f_face_normal */
    cs_real_t pfacj = 0.5;
    pfacj *= ofij_2[ii][0]*(grad_local_2[3*ii  ]+grad[cell_id][0])
            +ofij_2[ii][1]*(grad_local_2[3*ii+1]+grad[cell_id][1])
            +ofij_2[ii][2]*(grad_local_2[3*ii+2]+grad[cell_id][2]);
    if (c_weight != NULL) {
      pfacj += (1.0-rweight_2[ii]) * (pvar_local_2[ii] - pvar[cell_id]);
    } else {
      pfacj += (1.0-gweight_2[ii]) * (pvar_local_2[ii] - pvar[cell_id]);
    }

    for (int j = 0; j < 3; j++) {
      rhs[cell_id][j] += pfacj * b_f_face_normal[face_id][j];
    }

  }

  if (c_weight != NULL) {
    BFT_FREE(rweight_1);
    BFT_FREE(rweight_2);
  }
  BFT_FREE(grad_local_1);
  BFT_FREE(grad_local_2);
  BFT_FREE(pvar_local_1);
  BFT_FREE(pvar_local_2);
}

/*----------------------------------------------------------------------------
 * Add internal coupling rhs contribution for LSQ gradient calculation
 *
 * parameters:
 *   cpl       <-- pointer to coupling entity
 *   c_weight <-- weighted gradient coefficient variable, or NULL
 *   rhsv     <-> rhs contribution modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_rhs(const cs_internal_coupling_t  *cpl,
                             const cs_real_t                c_weight[],
                             cs_real_4_t                    rhsv[])
{
  int ii, ll;
  cs_lnum_t face_id, cell_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_lnum_t n_dist_1, n_dist_2;
  cs_lnum_t *dist_loc_1, *dist_loc_2;
  cs_real_t pfac;
  cs_real_3_t dc, fctb;
  cs_real_t *pvar_local_1, *pvar_local_2;
  cs_real_t *pvar_distant_1, *pvar_distant_2;
  cs_real_t *rweight_1, *rweight_2;
  cs_real_3_t *ij_1, *ij_2;
  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     &n_dist_1,
                                     &n_dist_2,
                                     &dist_loc_1,
                                     &dist_loc_2);

  ij_1 = cpl->ij_1;
  ij_2 = cpl->ij_2;

  /* 1 & 2 exchange pvar stored in rhsv[][3] */
  BFT_MALLOC(pvar_distant_1, n_dist_1, cs_real_t);
  BFT_MALLOC(pvar_distant_2, n_dist_2, cs_real_t);

  for (ii = 0; ii < n_dist_1; ii++) {
    face_id = dist_loc_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    pvar_distant_1[ii] = rhsv[cell_id][3];
  }

  for (ii = 0; ii < n_dist_2; ii++) {
    face_id = dist_loc_2[ii] - 1;
    cell_id = b_face_cells[face_id];
    pvar_distant_2[ii] = rhsv[cell_id][3];
  }

  BFT_MALLOC(pvar_local_1, n_1, cs_real_t);
  BFT_MALLOC(pvar_local_2, n_2, cs_real_t);

  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    pvar_distant_1,
                                    pvar_distant_2,
                                    pvar_local_1,
                                    pvar_local_2);

  BFT_FREE(pvar_distant_1);
  BFT_FREE(pvar_distant_2);

  /* Compute rhs */
  if (c_weight != NULL) { /* Heterogenous diffusivity */
    BFT_MALLOC(rweight_1, n_1, cs_real_t);
    BFT_MALLOC(rweight_2, n_2, cs_real_t);
    _cs_internal_coupling_exchange_rhs_weight(cpl,
                                              c_weight, /* diffusivity */
                                              rweight_1, /* rhs weight */
                                              rweight_2); /* rhs weight */
  }

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    for (ll = 0; ll < 3; ll++)
      dc[ll] = ij_1[ii][ll];

    pfac = (pvar_local_1[ii] - rhsv[cell_id][3])
         / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
    if (c_weight != NULL) pfac *= rweight_1[ii];

    for (ll = 0; ll < 3; ll++)
      fctb[ll] = dc[ll] * pfac;

    for (ll = 0; ll < 3; ll++)
      rhsv[cell_id][ll] += fctb[ll];
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    for (ll = 0; ll < 3; ll++)
      dc[ll] = ij_2[ii][ll];

    pfac = (pvar_local_2[ii] - rhsv[cell_id][3])
         / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
    if (c_weight != NULL) pfac *= rweight_2[ii];

    for (ll = 0; ll < 3; ll++)
      fctb[ll] = dc[ll] * pfac;

    for (ll = 0; ll < 3; ll++)
      rhsv[cell_id][ll] += fctb[ll];
  }

  if (c_weight != NULL) {
    BFT_FREE(rweight_1);
    BFT_FREE(rweight_2);
  }
  BFT_FREE(pvar_local_1);
  BFT_FREE(pvar_local_2);
}

/*----------------------------------------------------------------------------
 * Modify LSQ COCG matrix to include internal coupling
 *
 * parameters:
 *   coupling <-- pointer to coupling entity
 *   cocg            <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                           cs_real_33_t                   cocg[])
{
  int ii, ll, mm;
  cs_lnum_t face_id, cell_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_real_t umdddij;
  cs_real_3_t dddij;
  cs_real_3_t *ij_1, *ij_2;
  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  ij_1 = cpl->ij_1;
  ij_2 = cpl->ij_2;

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    for (ll = 0; ll < 3; ll++)
      dddij[ll] = ij_1[ii][ll];

    umdddij = 1./ cs_math_3_norm(dddij);
    for (ll = 0; ll < 3; ll++)
      dddij[ll] *= umdddij;

    for (ll = 0; ll < 3; ll++) {
      for (mm = 0; mm < 3; mm++) {
        cocg[cell_id][ll][mm] += dddij[ll]*dddij[mm];
      }
    }
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];
    for (ll = 0; ll < 3; ll++)
      dddij[ll] = ij_2[ii][ll];

    umdddij = 1./ cs_math_3_norm(dddij);
    for (ll = 0; ll < 3; ll++)
      dddij[ll] *= umdddij;

    for (ll = 0; ll < 3; ll++) {
      for (mm = 0; mm < 3; mm++) {
        cocg[cell_id][ll][mm] += dddij[ll]*dddij[mm];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Modify iterative COCG matrix to include internal coupling
 *
 * parameters:
 *   coupling <-- pointer to coupling entity
 *   cocg            <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_it_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                          cs_real_33_t                   cocg[])
{
  int ii, ll, mm;
  cs_lnum_t cell_id, face_id;
  cs_real_t dvol;

  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_real_3_t *ofij_1, *ofij_2;

  const cs_mesh_t* m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_mesh_quantities_t* fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  ofij_1 = cpl->ofij_1;
  ofij_2 = cpl->ofij_2;

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];

    for (ll = 0; ll < 3; ll++) {
      for (mm = 0; mm < 3; mm++) {
        cocg[cell_id][ll][mm] -= 0.5 * ofij_1[ii][ll] * b_f_face_normal[face_id][mm];
      }
    }
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    for (ll = 0; ll < 3; ll++) {
      for (mm = 0; mm < 3; mm++) {
        cocg[cell_id][ll][mm] -= 0.5 * ofij_2[ii][ll] * b_f_face_normal[face_id][mm];
      }
    }
  }

# pragma omp parallel for private(dvol, ll, mm)
  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    dvol = 1. / cell_vol[cell_id];
    for (ll = 0; ll < 3; ll++) {
      for (mm = 0; mm < 3; mm++) {
         cocg[cell_id][ll][mm] *= dvol;
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Initialize coupling criteria from strings.
 *
 * parameters:
 *   criteria_cells_1  <-- string criteria for the first group of cells
 *   criteria_cells_2  <-- string criteria for the second group of cells
 *   criteria_juncture <-- string criteria for the juncture, which is a
 *                         group of faces
 *   cpl               --> pointer to coupling structure to initialize
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_criteria_initialize(const char   criteria_cells_1[],
                                         const char   criteria_cells_2[],
                                         const char   criteria_juncture[],
                                         cs_internal_coupling_t  *cpl)
{
  BFT_MALLOC(cpl->criteria_cells_1,
             strlen(criteria_cells_1)+1,
             char);
  strcpy(cpl->criteria_cells_1, criteria_cells_1);

  BFT_MALLOC(cpl->criteria_cells_2,
             strlen(criteria_cells_2)+1,
             char);
  strcpy(cpl->criteria_cells_2, criteria_cells_2);

  BFT_MALLOC(cpl->criteria_juncture,
             strlen(criteria_juncture)+1,
             char);
  strcpy(cpl->criteria_juncture, criteria_juncture);
}

/*----------------------------------------------------------------------------
 * Initialize locators using selection criteria.
 *
 * parameters:
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_locators_initialize(cs_internal_coupling_t  *cpl)
{
  int k1, k2;
  cs_lnum_t   face_id, cell_id;
  cs_lnum_t   n_selected_faces = 0, n_selected_cells_1 = 0,
    n_selected_cells_2 = 0;
  cs_lnum_t  *selected_faces = NULL,
    *selected_cells_1 = NULL, *selected_cells_2 = NULL;
  cs_lnum_t *faces_1 = NULL, *faces_2 = NULL;

  char* tag_cells;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t ncel = m->n_cells_with_ghosts;
  const cs_lnum_t* b_face_cells = m->b_face_cells;

  /* Selection of the juncture */

  BFT_MALLOC(selected_faces, m->n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list(cpl->criteria_juncture,
                              &n_selected_faces,
                              selected_faces); /* 0..n-1 */

  qsort(selected_faces, n_selected_faces, sizeof(int), _compare_int);

  /* Selection of cells_1 using volumic selection criteria*/

  BFT_MALLOC(selected_cells_1, ncel, cs_lnum_t);

  cs_selector_get_cell_num_list(cpl->criteria_cells_1,
                                &n_selected_cells_1,
                                selected_cells_1);
  /* convert to 0..n-1 */
  for (int ii = 0; ii < n_selected_cells_1; ii++) {
    selected_cells_1[ii]--;
  }

  /* Selection of cells_2 using volumic selection criteria */

  BFT_MALLOC(selected_cells_2, ncel, cs_lnum_t);

  cs_selector_get_cell_num_list(cpl->criteria_cells_2,
                                &n_selected_cells_2,
                                selected_cells_2);
  /* convert to 0..n-1 */
  for (int ii = 0; ii < n_selected_cells_2; ii++) {
    selected_cells_2[ii]--;
  }

  /* Sort arrays of selected cells */
  qsort(selected_cells_1, n_selected_cells_1, sizeof(int), _compare_int);
  qsort(selected_cells_2, n_selected_cells_2, sizeof(int), _compare_int);

  BFT_MALLOC(tag_cells, ncel, char);

  for (int ii = 0; ii < n_selected_cells_1; ii++) {
    cell_id = selected_cells_1[ii];
    tag_cells[cell_id] = 1;
  }

  for (int ii = 0; ii < n_selected_cells_2; ii++) {
    cell_id = selected_cells_2[ii];
    tag_cells[cell_id] = 2;
  }

  BFT_FREE(selected_cells_1);
  BFT_FREE(selected_cells_2);

  BFT_MALLOC(faces_1, n_selected_faces, cs_lnum_t);
  BFT_MALLOC(faces_2, n_selected_faces, cs_lnum_t);

  k1 = k2 = 0;
  for (int ii = 0; ii < n_selected_faces; ii++) {
    face_id = selected_faces[ii];
    cell_id = b_face_cells[face_id];
    if (tag_cells[cell_id] == 1) {
      faces_1[k1++] = face_id + 1;
    } else if (tag_cells[cell_id] == 2) {
      faces_2[k2++] = face_id + 1;
    } else {
      bft_error(__FILE__, __LINE__, 0, "Cell tag not found");
    }
  }

  BFT_REALLOC(faces_1, k1, cs_lnum_t);
  BFT_REALLOC(faces_2, k2, cs_lnum_t);

  BFT_FREE(tag_cells);

  /* Initialize locators */
  cpl->faces_1 = faces_1;
  cpl->n_1 = k1;
  cpl->faces_2 = faces_2;
  cpl->n_2 = k2;

  cpl->locator_1 =
    _cs_internal_coupling_create_locator(cpl, 1);

  cpl->n_dist_1 =
    ple_locator_get_n_dist_points(cpl->locator_1);

  BFT_MALLOC(cpl->dist_loc_1,
             cpl->n_dist_1,
             cs_lnum_t);
  memcpy(cpl->dist_loc_1,
         ple_locator_get_dist_locations(cpl->locator_1),
         (cpl->n_dist_1)*sizeof(cs_lnum_t));

  cpl->locator_2 =
    _cs_internal_coupling_create_locator(cpl, 2);

  cpl->n_dist_2 =
    ple_locator_get_n_dist_points(cpl->locator_2);

  BFT_MALLOC(cpl->dist_loc_2,
             cpl->n_dist_2,
             cs_lnum_t);
  memcpy(cpl->dist_loc_2,
         ple_locator_get_dist_locations(cpl->locator_2),
         (cpl->n_dist_2)*sizeof(cs_lnum_t));

  /* Geometric quantities */
  BFT_MALLOC(cpl->gweight_1, cpl->n_1, cs_real_t);
  BFT_MALLOC(cpl->gweight_2, cpl->n_2, cs_real_t);
  BFT_MALLOC(cpl->ij_1, cpl->n_1, cs_real_3_t);
  BFT_MALLOC(cpl->ij_2, cpl->n_2, cs_real_3_t);
  BFT_MALLOC(cpl->ofij_1, cpl->n_1, cs_real_3_t);
  BFT_MALLOC(cpl->ofij_2, cpl->n_2, cs_real_3_t);

  cs_internal_coupling_exchange_ij(cpl);

  /* Allocate coupling exchange coefficients */

  BFT_MALLOC(cpl->hint_1, cpl->n_1, cs_real_t);
  BFT_MALLOC(cpl->hext_1, cpl->n_1, cs_real_t);

  BFT_MALLOC(cpl->hint_2, cpl->n_2, cs_real_t);
  BFT_MALLOC(cpl->hext_2, cpl->n_2, cs_real_t);

  BFT_MALLOC(cpl->coupled_faces, m->n_b_faces, bool);

  cpl->cocgb_s_lsq = NULL;
  cpl->cocgb_s_it = NULL;
  cpl->cocg_s_it = NULL;

  /* Release allocated memory */
  BFT_FREE(selected_faces);
}

/*----------------------------------------------------------------------------
 * Destruction of all internal coupling related structures.
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_finalize(void)
{
  cs_internal_coupling_t* cpl;
  for (int ii = 0; ii < _n_internal_couplings; ii++) {
    cpl = _internal_coupling + ii;
    _cs_internal_coupling_destroy_entity(cpl);
  }
  BFT_FREE(_internal_coupling);
  _n_internal_couplings = 0;
}

/*----------------------------------------------------------------------------
 * Exchange variable between groups using face id
 *
 * parameters:
 *   cpl     <-- pointer to coupling entity
 *   stride  <-- number of values (non interlaced) by entity
 *   tab     <-- variable exchanged
 *   local_1 --> local data for group 1
 *   local_2 --> local data for group 2
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_face_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local_1[],
                                         cs_real_t                      local_2[])
{
  int ii, jj;
  cs_lnum_t face_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_lnum_t n_dist_1, n_dist_2;
  cs_lnum_t *dist_loc_1, *dist_loc_2;
  cs_real_t *distant_1, *distant_2;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     &n_dist_1,
                                     &n_dist_2,
                                     &dist_loc_1,
                                     &dist_loc_2);


  BFT_MALLOC(distant_1, n_dist_1*stride, cs_real_t);
  BFT_MALLOC(distant_2, n_dist_2*stride, cs_real_t);

  for (ii = 0; ii < n_dist_1; ii++) {
    face_id = dist_loc_1[ii] - 1;
    for (jj = 0; jj < stride; jj++) {
      distant_1[stride * ii + jj] = tab[stride * face_id + jj];
    }
  }

  for (ii = 0; ii < n_dist_2; ii++) {
    face_id = dist_loc_2[ii] - 1;
    for (jj = 0; jj < stride; jj++) {
      distant_2[stride * ii + jj] = tab[stride * face_id + jj];
    }
  }

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    distant_1,
                                    distant_2,
                                    local_1,
                                    local_2);

  BFT_FREE(distant_1);
  BFT_FREE(distant_2);
}

/*----------------------------------------------------------------------------
 * Exchange variable between groups using cell id
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
 *   stride   <-- number of values (non interlaced) by entity
 *   tab      <-- variable exchanged
 *   local_1  --> local data for group 1
 *   local_2  --> local data for group 2
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_cell_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local_1[],
                                         cs_real_t                      local_2[])
{
  int ii, jj;
  cs_lnum_t face_id, cell_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_lnum_t n_dist_1, n_dist_2;
  cs_lnum_t *dist_loc_1, *dist_loc_2;

  cs_real_t *distant_1, *distant_2;

  const cs_mesh_t* m = cs_glob_mesh;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     &n_dist_1,
                                     &n_dist_2,
                                     &dist_loc_1,
                                     &dist_loc_2);


  BFT_MALLOC(distant_1, n_dist_1*stride, cs_real_t);
  BFT_MALLOC(distant_2, n_dist_2*stride, cs_real_t);

  for (ii = 0; ii < n_dist_1; ii++) {
    face_id = dist_loc_1[ii] - 1;
    cell_id = b_face_cells[face_id];
    for (jj = 0; jj < stride; jj++) {
      distant_1[stride * ii + jj] = tab[stride * cell_id + jj];
    }
  }

  for (ii = 0; ii < n_dist_2; ii++) {
    face_id = dist_loc_2[ii] - 1;
    cell_id = b_face_cells[face_id];
    for (jj = 0; jj < stride; jj++) {
      distant_2[stride * ii + jj] = tab[stride * cell_id + jj];
    }
  }

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    distant_1,
                                    distant_2,
                                    local_1,
                                    local_2);

  BFT_FREE(distant_1);
  BFT_FREE(distant_2);
}

/*----------------------------------------------------------------------------
 * Exchange quantities from distant to local
 *
 * parameters:
 *   cpl       <-- pointer to coupling entity
 *   stride    <-- Stride (e.g. 1 for double, 3 for interleaved coordinates)
 *   distant_1 <-- Distant values 1, size coupling->n_dist_1
 *   distant_2 <-- Distant values 2, size coupling->n_dist_2
 *   local_1   --> Local values 1, size coupling->n_1
 *   local_2   --> Local values 2, size coupling->n_2
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_var(const cs_internal_coupling_t  *cpl,
                                  int                            stride,
                                  cs_real_t                      distant_1[],
                                  cs_real_t                      distant_2[],
                                  cs_real_t                      local_1[],
                                  cs_real_t                      local_2[])
{
  _cs_internal_coupling_send_var(cpl,
                                 stride,
                                 1,
                                 local_1,
                                 distant_1);

  _cs_internal_coupling_send_var(cpl,
                                 stride,
                                 2,
                                 local_2,
                                 distant_2);
}

/*----------------------------------------------------------------------------
 * Compute and exchange ij vectors
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_ij(const cs_internal_coupling_t  *cpl)
{
  int ii;
  cs_lnum_t face_id, cell_id;

  int *faces_1, *faces_2;
  cs_lnum_t n_1, n_2;

  cs_real_3_t* ij_1 = cpl->ij_1;
  cs_real_3_t* ij_2 = cpl->ij_2;

  cs_real_t xxd, yyd, zzd, xxl, yyl, zzl;

  const cs_mesh_quantities_t* mq = cs_glob_mesh_quantities;
  const cs_mesh_t* m = cs_glob_mesh;

  const cs_real_t* cell_cen = mq->cell_cen;
  const cs_lnum_t* b_face_cells = m->b_face_cells;

  cs_real_t *cell_cen_1, *cell_cen_2;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  BFT_MALLOC(cell_cen_1, 3 * n_1, cs_real_t);
  BFT_MALLOC(cell_cen_2, 3 * n_2, cs_real_t);

  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           3,
                                           mq->cell_cen,
                                           cell_cen_1,
                                           cell_cen_2);

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];

    xxd = cell_cen_1[3*ii];
    yyd = cell_cen_1[3*ii + 1];
    zzd = cell_cen_1[3*ii + 2];

    xxl = cell_cen[3*cell_id];
    yyl = cell_cen[3*cell_id + 1];
    zzl = cell_cen[3*cell_id + 2];

    ij_1[ii][0] = xxd - xxl;
    ij_1[ii][1] = yyd - yyl;
    ij_1[ii][2] = zzd - zzl;
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    xxd = cell_cen_2[3*ii];
    yyd = cell_cen_2[3*ii + 1];
    zzd = cell_cen_2[3*ii + 2];

    xxl = cell_cen[3*cell_id];
    yyl = cell_cen[3*cell_id + 1];
    zzl = cell_cen[3*cell_id + 2];

    ij_2[ii][0] = xxd - xxl;
    ij_2[ii][1] = yyd - yyl;
    ij_2[ii][2] = zzd - zzl;
  }

  BFT_FREE(cell_cen_1);
  BFT_FREE(cell_cen_2);

  /* Compute geometric weights and iterative reconstruction vector */
  _cs_internal_coupling_exchange_gweight(cpl);
  _cs_internal_coupling_ofij(cpl);
}

/*----------------------------------------------------------------------------
 * Return pointers to coupling components
 *
 * parameters:
 *   cpl             <-- pointer to coupling entity
 *   n1              --> NULL or pointer to component n1
 *   n2              --> NULL or pointer to component n2
 *   fac_1[]         --> NULL or pointer to component fac_1[]
 *   fac_2[]         --> NULL or pointer to component fac_2[]
 *   n_dist_1        --> NULL or pointer to component n_dist_1
 *   n_dist_2        --> NULL or pointer to component n_dist_2
 *   dist_loc_1[]    --> NULL or pointer to component dist_loc_1[]
 *   dist_loc_2[]    --> NULL or pointer to component dist_loc_2[]
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_coupled_faces(const cs_internal_coupling_t  *cpl,
                                   cs_lnum_t                     *n1,
                                   cs_lnum_t                     *n2,
                                   cs_lnum_t                     *fac_1[],
                                   cs_lnum_t                     *fac_2[],
                                   cs_lnum_t                     *n_dist_1,
                                   cs_lnum_t                     *n_dist_2,
                                   cs_lnum_t                     *dist_loc_1[],
                                   cs_lnum_t                     *dist_loc_2[])
{
  if (n1 != NULL) {
    *n1 = cpl->n_1;
  }
  if (n2 != NULL) {
    *n2 = cpl->n_2;
  }
  if (fac_1 != NULL) {
    *fac_1 = cpl->faces_1;
  }
  if (fac_2 != NULL) {
    *fac_2 = cpl->faces_2;
  }
  if (n_dist_1 != NULL) {
    *n_dist_1 = cpl->n_dist_1;
  }
  if (n_dist_2 != NULL) {
    *n_dist_2 = cpl->n_dist_2;
  }
  if (dist_loc_1 != NULL) {
    *dist_loc_1 = cpl->dist_loc_1;
  }
  if (dist_loc_2 != NULL) {
    *dist_loc_2 = cpl->dist_loc_2;
  }
}

/*----------------------------------------------------------------------------
 * Return the coupling associated with a given coupling_id.
 *
 * parameters:
 *   coupling_id <-> id associated with a coupling entity
 *----------------------------------------------------------------------------*/

cs_internal_coupling_t *
cs_internal_coupling_by_id(int coupling_id)
{
  if (coupling_id > -1 && coupling_id < _n_internal_couplings) {
    return _internal_coupling + coupling_id;
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              "coupling_id = %d provided is invalid", coupling_id);
  }
  return (cs_internal_coupling_t*)NULL;
}

/*----------------------------------------------------------------------------
 * Modify matrix-vector product in case of internal coupling
 *
 * parameters:
 *   exclude_diag <-- extra diagonal flag
 *   matrix       <-- matrix m in m * x = y
 *   x            <-- vector x in m * x = y
 *   y            <-> vector y in m * x = y
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_spmv_contribution(bool              exclude_diag,
                                       void             *input,
                                       const cs_real_t  *restrict x,
                                       cs_real_t        *restrict y)
{
  cs_lnum_t face_id, cell_id;
  int ii;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;
  cs_real_t *x_j_1, *x_j_2;
  cs_real_t pi, pj;
  cs_real_t hint, hext, heq;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_glob_mesh->b_face_cells;
  const cs_internal_coupling_t* cpl = (const cs_internal_coupling_t *)input;

  const cs_real_t thetap = cpl->thetav;
  const int       idiffp = cpl->idiff;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  BFT_MALLOC(x_j_1, n_1, cs_real_t);
  BFT_MALLOC(x_j_2, n_2, cs_real_t);

  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           1,
                                           x,
                                           x_j_1,
                                           x_j_2);

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    cell_id = b_face_cells[face_id];

    cs_real_t fluxi = 0.;
    pi = exclude_diag ? 0. : x[cell_id];/* If exclude_diag, no diagonal term */
    pj = x_j_1[ii];

    hint = cpl->hint_1[ii];
    hext = cpl->hext_1[ii];
    heq = hint * hext / (hint + hext);

    cs_b_diff_flux_coupling(idiffp,
                            pi,
                            pj,
                            heq,
                            &fluxi);

    y[cell_id] += thetap * fluxi;
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    cell_id = b_face_cells[face_id];

    cs_real_t fluxi = 0.;
    pi = exclude_diag ? 0. : x[cell_id];
    pj = x_j_2[ii];

    hint = cpl->hint_2[ii];
    hext = cpl->hext_2[ii];
    heq = hint * hext / (hint + hext);

    cs_b_diff_flux_coupling(idiffp,
                            pi,
                            pj,
                            heq,
                            &fluxi);

    y[cell_id] += thetap * fluxi;
  }
  BFT_FREE(x_j_1);
  BFT_FREE(x_j_2);
}

/*----------------------------------------------------------------------------
 * Initialize internal coupling related structures.
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize(void)
{
  int field_id;
  cs_field_t* f, *f_diff;
  int coupling_key_id = cs_field_key_id("coupling_entity");
  int coupling_id = 0;

  const int diffusivity_key_id = cs_field_key_id("scalar_diffusivity_id");
  int diffusivity_id;

  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  const int n_fields = cs_field_n_fields();

  /* Definition of coupling_ids as keys of variable fields */
  for (field_id = 0; field_id < n_fields; field_id++) {
    f = cs_field_by_id(field_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      if (var_cal_opt.icoupl > 0) {
        /* Set same coupling_id for associated diffusivity field */
        cs_field_set_key_int(f, coupling_key_id, coupling_id);
        diffusivity_id = cs_field_get_key_int(f, diffusivity_key_id);
        if (diffusivity_id > -1) {
          f_diff = cs_field_by_id(diffusivity_id);
          cs_field_set_key_int(f_diff, coupling_key_id, coupling_id);
        }
        coupling_id++;
      }
    }
  }

  /* Initialization of coupling entities */
  coupling_id = 0;
  for (field_id = 0; field_id < n_fields; field_id++) {
    f = cs_field_by_id(field_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      if (var_cal_opt.icoupl > 0) {
        cs_internal_coupling_t* cpl =
          _internal_coupling + coupling_id;

        /* Definition of var_cal_opt options
         * (needed for matrix.vector multiply) */
        cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
        cpl->thetav = var_cal_opt.thetav;
        cpl->idiff = var_cal_opt.idiff;
        cpl->dim = cs_glob_mesh->dim;

        /* Initialize locators */
        cs_internal_coupling_locators_initialize(cpl);

        /* Initialize coupled_faces */
        _cs_internal_coupling_select_coupled_faces(cpl);

        /* Initialize cocg & cocgb */
        if (var_cal_opt.imrgra == 0) {
          cs_compute_cell_cocg_s_it_coupling(cs_glob_mesh,
                                              cs_glob_mesh_quantities,
                                              cpl);
        } else if (var_cal_opt.imrgra == 1) {
          cs_compute_cell_cocg_s_lsq_coupling(cs_glob_mesh,
                                              cs_glob_mesh_quantities,
                                              cpl);
        }

        /* Update user information */
        BFT_MALLOC(cpl->namesca, strlen(f->name) + 1, char);
        strcpy(cpl->namesca, f->name);

        coupling_id++;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Print informations about the given coupling entity
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_print(const cs_internal_coupling_t  *cpl)
{

  if (cpl == NULL)
    return;

  int ii;
  cs_lnum_t n_1, n_2, face_id;
  cs_lnum_t *faces_1, *faces_2;
  const cs_real_t *hint_1 = cpl->hint_1,
    *hext_1 = cpl->hext_1;
  const cs_real_t *hint_2 = cpl->hint_2,
    *hext_2 = cpl->hext_2;


  const cs_real_t* cog = cs_glob_mesh_quantities->b_face_cog;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  const cs_real_t* gweight_1 = cpl->gweight_1;
  const cs_real_t* gweight_2 = cpl->gweight_2;

  cs_real_t *rweight_1 = NULL, *rweight_2 = NULL;
  BFT_MALLOC(rweight_1, n_1, cs_real_t);
  BFT_MALLOC(rweight_2, n_2, cs_real_t);
  cs_field_t *f = cs_field_by_name(cpl->namesca);
  int key_id = cs_field_key_id("gradient_weighting_id");
  int diff_id = cs_field_get_key_int(f, key_id);
  cs_field_t *weight_f = cs_field_by_id(diff_id);
  cs_real_t *c_weight = weight_f->val;
  if (c_weight != NULL) {
    _cs_internal_coupling_exchange_rhs_weight(cpl,
                                              c_weight, /* diffusivity */
                                              rweight_1, /* rhs weight */
                                              rweight_2); /* rhs weight */
  } else {
    for (ii = 0; ii < n_1; ii++) {
      rweight_1[ii] = 0;
    }
    for (ii = 0; ii < n_2; ii++) {
      rweight_2[ii] = 0;
    }
  }

  bft_printf("Coupled scalar %s\n"
       "Group 1 selection criterion : %s\n"
             "Group 2 selection criterion : %s\n"
       "Juncture selection criterion : %s\n",
       cpl->namesca,
       cpl->criteria_cells_1,
       cpl->criteria_cells_2,
       cpl->criteria_juncture);

  if (n_1 > 0) {
  bft_printf("\nFaces in group 1 (%d face%s) :\n"
             "--------------------------------------------------------------\n"
             "face\tx\t\ty\t\tz\t\thint\t\thext\t\tgweight\t\trweight\t\t\n"
             "--------------------------------------------------------------\n",
             n_1,
             n_1 <=1 ? "" : "s");
  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    bft_printf("%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t\n",
               face_id,
               cog[3*face_id],
               cog[3*face_id+1],
               cog[3*face_id+2],
               hint_1[ii],
               hext_1[ii],
               gweight_1[ii],
               rweight_1[ii]);
  }
  bft_printf("-------------------------------------------------------------\n");
  }
  else {
    bft_printf("\nNo faces in group 1\n");
  }


  if (n_2 > 0) {
    bft_printf("\nFaces in group 2 (%d face%s) :\n"
               "--------------------------------------------------------------\n"
               "face\tx\t\ty\t\tz\t\thint\t\thext\t\tgweight\t\trweight\t\t\n"
               "--------------------------------------------------------------\n",
               n_2,
               n_2<=1 ? "" : "s");
  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    bft_printf("%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t\n",
               face_id,
               cog[3*face_id],
               cog[3*face_id+1],
               cog[3*face_id+2],
               hint_2[ii],
               hext_2[ii],
               gweight_2[ii],
               rweight_2[ii]);
  }
  bft_printf("---------------------------------------------------------\n");
  } else {
    bft_printf("\nNo faces in group 2\n");
  }

  BFT_FREE(rweight_1);
  BFT_FREE(rweight_2);
}

/*----------------------------------------------------------------------------
 * Print informations about all coupling entities
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_dump(void)
{
  int ii;
  cs_internal_coupling_t* cpl;
  bft_printf("\ncs_internal_coupling_dump\n");
  for (ii = 0; ii < _n_internal_couplings; ii++) {
    cpl = _internal_coupling + ii;
    bft_printf("coupling_id = %d\n", ii);
    cs_internal_coupling_print(cpl);
  }
}

/*----------------------------------------------------------------------------
 * Define coupling entity using given criterias
 *
 * parameters:
 *   field_id   <-- id of the field
 *   volume_1[] <-- string criteria for the first group of cells
 *   volume_2[] <-- string criteria for the second group of cells
 *   juncture[] <-- string criteria for the juncture, which is a
 *                  group of faces
 *----------------------------------------------------------------------------*/

int
cs_internal_coupling_add_entity(int        f_id,
                                const char volume_1[],
                                const char volume_2[],
                                const char juncture[])
{
  int coupling_id;
  cs_internal_coupling_t* cpl;

  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  cs_field_t* f = cs_field_by_id(f_id);

  if (f->type & CS_FIELD_VARIABLE) {
    BFT_REALLOC(_internal_coupling, _n_internal_couplings + 1,
                cs_internal_coupling_t);
    cpl = _internal_coupling + _n_internal_couplings;

    cs_internal_coupling_criteria_initialize(volume_1,
                                             volume_2,
                                             juncture,
                                             cpl);

    cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.icoupl = 1;
    cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);

    coupling_id = _n_internal_couplings;
    _n_internal_couplings++;
    return coupling_id;
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              "field id = %d provided is invalid."
              " The field must be a variable.",
              f_id);
    return -1;
  }
}

/*----------------------------------------------------------------------------
 * Update components hint_* and hext_* using hbord
 *   in the coupling entity associated with given field_id
 *
 * parameters:
 *   field_id <-- id of the field
 *   hbord    <-- array used to update hint_* and hext_*
 *----------------------------------------------------------------------------*/

void
cs_ic_set_exchcoeff(const int         field_id,
                    const cs_real_t  *hbord)
{
  int ii;
  cs_lnum_t face_id;
  cs_lnum_t n_1, n_2;
  cs_lnum_t *faces_1, *faces_2;

  cs_real_t surf;
  const cs_real_t* b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  const cs_field_t* f = cs_field_by_id(field_id);
  const cs_int_t coupling_key_id = cs_field_key_id("coupling_entity");
  int coupling_id = cs_field_get_key_int(f, coupling_key_id);
  const cs_internal_coupling_t  *cpl
    = cs_internal_coupling_by_id(coupling_id);

  cs_real_t *hint_1 = cpl->hint_1;
  cs_real_t *hext_1 = cpl->hext_1;

  cs_real_t *hint_2 = cpl->hint_2;
  cs_real_t *hext_2 = cpl->hext_2;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_1,
                                     &n_2,
                                     &faces_1,
                                     &faces_2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL);

  cs_internal_coupling_exchange_by_face_id(cpl,
                                           1,
                                           hbord,
                                           hext_1,
                                           hext_2);

  for (ii = 0; ii < n_1; ii++) {
    face_id = faces_1[ii] - 1;
    surf = b_face_surf[face_id];
    hint_1[ii] = hbord[face_id] * surf;
    hext_1[ii] *= surf;
  }

  for (ii = 0; ii < n_2; ii++) {
    face_id = faces_2[ii] - 1;
    surf = b_face_surf[face_id];
    hint_2[ii] = hbord[face_id] * surf;
    hext_2[ii] *= surf;
  }
}

/*----------------------------------------------------------------------------
 * Add contribution from coupled faces (internal coupling) to polynomial
 * preconditionning.
 *
 * This function is common to most solvers
 *
 * parameters:
 *   input  <-- input
 *   ad     <-> diagonal part of linear equation matrix
 *----------------------------------------------------------------------------*/

void
cs_matrix_preconditionning_add_coupling_contribution(void       *input,
                                                     cs_real_t  *ad)
{
  const cs_internal_coupling_t* cpl = (const cs_internal_coupling_t *)input;
  if (cpl != NULL) {
    int ii;
    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *restrict)cs_glob_mesh->b_face_cells;
    cs_lnum_t *faces_1, *faces_2, face_id, cell_id;
    int n_1, n_2;
    cs_real_t hint, hext, heq;

    cs_internal_coupling_coupled_faces(cpl,
                                       &n_1,
                                       &n_2,
                                       &faces_1,
                                       &faces_2,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL);

    for (ii = 0; ii < n_1; ii++) {
      face_id = faces_1[ii] - 1;
      cell_id = b_face_cells[face_id];

      hint = cpl->hint_1[ii];
      hext = cpl->hext_1[ii];
      heq = hint * hext / (hint + hext);

      ad[cell_id] += heq;
    }

    for (ii = 0; ii < n_2; ii++) {
      face_id = faces_2[ii] - 1;
      cell_id = b_face_cells[face_id];

      hint = cpl->hint_2[ii];
      hext = cpl->hext_2[ii];
      heq = hint * hext / (hint + hext);

      ad[cell_id] += heq;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

