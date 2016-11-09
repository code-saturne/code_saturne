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
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_boundary.h"
#include "cs_mesh_quantities.h"
#include "cs_convection_diffusion.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_selector.h"
#include "cs_prototypes.h"

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
 *----------------------------------------------------------------------------*/

static ple_locator_t*
_create_locator(cs_internal_coupling_t  *cpl)
{
  cs_lnum_t i, j;
  cs_lnum_t ifac;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t n_distant = n_local;//FIXME ???
  const int *tag = cpl->c_tag;
  cs_lnum_t *faces_local_num = NULL;
  const cs_lnum_t *faces_distant = cpl->faces_local;//FIXME

  fvm_nodal_t* nm = NULL;
  int *tag_nm = NULL;
  cs_lnum_t *faces_in_nm = NULL;

  ple_coord_t* point_coords = NULL;

  char mesh_name[16] = "locator";

  /* Create locator */

#if defined(PLE_HAVE_MPI)
  ple_locator_t *locator = ple_locator_create(cs_glob_mpi_comm,
                                              cs_glob_n_ranks,
                                              0);
#else
  ple_locator_t *locator = ple_locator_create();
#endif

  /* Create fvm_nodal_t structure */

  BFT_MALLOC(faces_local_num, n_local, cs_lnum_t);
  for (cs_lnum_t face_id = 0; face_id < n_local; face_id++)
    faces_local_num[face_id] = cpl->faces_local[face_id] + 1;

  nm = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                      mesh_name,
                                      false,
                                      0,
                                      n_local,
                                      NULL,
                                      faces_local_num); /* 1..n */

  /* Tag fvm_nodal_t structure */

  /* Number of faces to tag */
  const int nfac_in_nm = fvm_nodal_get_n_entities(nm, 2);
  /* Memory allocation */
  BFT_MALLOC(faces_in_nm, nfac_in_nm, cs_lnum_t);
  BFT_MALLOC(tag_nm, nfac_in_nm, int);
  /* Get id of faces to tag in parent */
  fvm_nodal_get_parent_num(nm, 2, faces_in_nm);
  /* Tag faces */
  for (cs_lnum_t ii = 0; ii < nfac_in_nm; ii++) {
    /* Default tag is 0 */
    tag_nm[ii] = 0;
    for (cs_lnum_t jj = 0; jj < n_local; jj++) {
      if (faces_in_nm[ii] == faces_local_num[jj]){
        tag_nm[ii] = tag[jj];
        break;
      }
    }
  }
  fvm_nodal_set_tag(nm, tag_nm, 2);
  /* Free memory */
  BFT_FREE(faces_in_nm);
  BFT_FREE(tag_nm);
  BFT_FREE(faces_local_num);

  /* Creation of distant group cell centers */

  BFT_MALLOC(point_coords, 3*n_distant, cs_real_t);

  for (i = 0; i < n_distant; i++) {
    ifac = faces_distant[i]; /* 0..n-1 */
    for (j = 0; j < 3; j++)
      point_coords[3*i+j] = cs_glob_mesh_quantities->b_face_cog[3*ifac+j];
  }

  /* Locator initialization */

  ple_locator_set_mesh(locator,
                       nm,
                       NULL,
                       0,
                       1.1, /* TODO */
                       cs_glob_mesh->dim,
                       n_distant,
                       NULL,
                       tag,
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
_destroy_entity(cs_internal_coupling_t  *cpl)
{
  BFT_FREE(cpl->c_tag);
  BFT_FREE(cpl->faces_local);
  BFT_FREE(cpl->faces_distant);
  BFT_FREE(cpl->h_int);
  BFT_FREE(cpl->h_ext);
  BFT_FREE(cpl->g_weight);
  BFT_FREE(cpl->ci_cj_vect);
  BFT_FREE(cpl->offset_vect);
  BFT_FREE(cpl->coupled_faces);
  BFT_FREE(cpl->cocgb_s_lsq);
  BFT_FREE(cpl->cocgb_s_it);
  BFT_FREE(cpl->cocg_s_it);
  BFT_FREE(cpl->cells_criteria);
  BFT_FREE(cpl->namesca);
  ple_locator_destroy(cpl->locator);
}

/*----------------------------------------------------------------------------
 * Compute geometrical face weights around coupling interface.
 *
 * parameters:
 *   cpl <-- pointer to coupling structure
 *----------------------------------------------------------------------------*/

static void
_compute_geometrical_face_weight(const cs_internal_coupling_t  *cpl)
{
  int ii;
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;
  const cs_real_3_t *ci_cj_vect = (const cs_real_3_t *)cpl->ci_cj_vect;
  cs_real_t* g_weight = cpl->g_weight;

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_mesh_t             *m   = cs_glob_mesh;

  const cs_real_t* cell_cen = fvq->cell_cen;
  const cs_real_t* b_face_cog = fvq->b_face_cog;
  const cs_real_t* diipb = fvq->diipb;
  const cs_real_t* b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_t* b_face_normal = cs_glob_mesh_quantities->b_face_normal;

  /* Store local FI' distances in gweight_distant */

  cs_real_t *g_weight_distant = NULL;
  BFT_MALLOC(g_weight_distant, n_distant, cs_real_t);
  for (ii = 0; ii < n_distant; ii++) {
    cs_real_t dv[3];
    face_id = faces_distant[ii];
    cell_id = m->b_face_cells[face_id];

    for (int jj = 0; jj < 3; jj++)
      dv[jj] =  - diipb[3*face_id + jj] - cell_cen[3*cell_id +jj]
                + b_face_cog[3*face_id +jj];

    g_weight_distant[ii] = cs_math_3_norm(dv);
  }

  /* Exchange FI' distances */

  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    g_weight_distant,
                                    g_weight);
  /* Free memory */
  BFT_FREE(g_weight_distant);

  /* Normalise the distance to obtain geometrical weights */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    g_weight[ii] /= ( ci_cj_vect[ii][0]*b_face_normal[3*face_id+0]
                     + ci_cj_vect[ii][1]*b_face_normal[3*face_id+1]
                     + ci_cj_vect[ii][2]*b_face_normal[3*face_id+2])
                     / b_face_surf[face_id];
  }
}

/*----------------------------------------------------------------------------
 * Compute r_weight around coupling interface based on diffusivity c_weight.
 *
 * parameters:
 *   cpl        <-- pointer to coupling structure
 *   c_weight   <-- diffusivity
 *   r_weight   --> physical face weight
 *----------------------------------------------------------------------------*/

static void
_compute_physical_face_weight(const cs_internal_coupling_t  *cpl,
                              const cs_real_t        c_weight[],
                              cs_real_t              rweight[])
{
  int ii;
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_real_t* g_weight = cpl->g_weight;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Exchange c_weight */

  cs_real_t *c_weight_local = NULL;
  BFT_MALLOC(c_weight_local, n_local, cs_real_t);
  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           1,
                                           c_weight,
                                           c_weight_local);

  /* Compute rweight */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];
    cs_real_t ki = c_weight[cell_id];
    cs_real_t kj = c_weight_local[ii];
    cs_real_t pond = g_weight[ii];
    rweight[ii] = kj / ( pond * ki + (1. - pond) * kj);
  }

  /* Free memory */
  BFT_FREE(c_weight_local);
}

/*----------------------------------------------------------------------------
 * Compute offset vector on coupled faces
 *
 * parameters:
 *   cpl <-> pointer to coupling entity
 *----------------------------------------------------------------------------*/

static void
_compute_offset(const cs_internal_coupling_t  *cpl)
{
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_real_t *g_weight = cpl->g_weight;
  cs_real_3_t *offset_vect = cpl->offset_vect;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_real_t *cell_cen = mq->cell_cen;
  const cs_real_t *b_face_cog = mq->b_face_cog;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  /* Exchange cell center location */

  cs_real_t *cell_cen_local = NULL;
  BFT_MALLOC(cell_cen_local, 3 * n_local, cs_real_t);
  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           3,
                                           mq->cell_cen,
                                           cell_cen_local);

  /* Compute OF vectors */

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];

    for (int jj = 0; jj < 3; jj++) {
      cs_real_t xxd = cell_cen_local[3*ii + jj];
      cs_real_t xxl = cell_cen[3*cell_id + jj];
      offset_vect[ii][jj] = b_face_cog[3*face_id + jj]
        - (        g_weight[ii]*xxl
           + (1. - g_weight[ii])*xxd);
    }
  }
  /* Free memory */
  BFT_FREE(cell_cen_local);
}

/*----------------------------------------------------------------------------
 * Compute cell centers vectors IJ on coupled faces
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

static void
_compute_ci_cj_vect(const cs_internal_coupling_t  *cpl)
{
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  cs_real_3_t *ci_cj_vect = cpl->ci_cj_vect;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_mesh_t *m = cs_glob_mesh;

  const cs_real_t *cell_cen = mq->cell_cen;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  /* Exchange cell center location */

  cs_real_t *cell_cen_local = NULL;
  BFT_MALLOC(cell_cen_local, 3 * n_local, cs_real_t);
  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           3, /* dimension */
                                           mq->cell_cen,
                                           cell_cen_local);

  /* Compute IJ vectors */

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];

    for (int jj = 0; jj < 3; jj++) {
      cs_real_t xxd = cell_cen_local[3*ii + jj];
      cs_real_t xxl = cell_cen[3*cell_id + jj];
      ci_cj_vect[ii][jj] = xxd - xxl;
    }
  }
  /* Free memory */
  BFT_FREE(cell_cen_local);

  /* Compute geometric weights and iterative reconstruction vector */

  _compute_geometrical_face_weight(cpl);
  _compute_offset(cpl);
}


/*----------------------------------------------------------------------------
 * Define component coupled_faces[] of given coupling entity.
 *
 * parameters:
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

static void
_initialize_coupled_faces(cs_internal_coupling_t *cpl)
{
  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_mesh_t *m = cs_glob_mesh;
  bool *coupled_faces = cpl->coupled_faces;

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
    coupled_faces[face_id] = false;
  }

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = faces_local[ii];
    coupled_faces[face_id] = true;
  }
}

/*----------------------------------------------------------------------------
 * Initialize coupling criteria from strings.
 *
 * parameters:
 *   criteria_cells    <-- string criteria for the first group of cells
 *   criteria_faces    <-- string criteria for faces
 *   cpl               --> pointer to coupling structure to initialize
 *----------------------------------------------------------------------------*/

static void
_criteria_initialize(const char   criteria_cells[],
                     const char   criteria_faces[],
                     cs_internal_coupling_t  *cpl)
{
  BFT_MALLOC(cpl->cells_criteria,
             strlen(criteria_cells)+1,
             char);

  strcpy(cpl->cells_criteria, criteria_cells);

  if (criteria_faces == NULL)
    return;

  BFT_MALLOC(cpl->faces_criteria,
             strlen(criteria_faces)+1,
             char);

  strcpy(cpl->faces_criteria, criteria_faces);

}

/*----------------------------------------------------------------------------
 * Initialize internal coupling locators using cell selection criteria ONLY.
 *
 * parameters:
 *   m   <->  pointer to mesh structure to modify
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

static void
_volume_initialize(cs_mesh_t               *m,
                   cs_internal_coupling_t  *cpl)
{
  cs_lnum_t  n_selected_cells;
  cs_lnum_t *selected_cells = NULL;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  /* Selection of Volume zone using volumic selection criteria*/

  BFT_MALLOC(selected_cells, n_cells_ext, cs_lnum_t);
  cs_selector_get_cell_list(cpl->cells_criteria,
                            &n_selected_cells,
                            selected_cells);

  /* Store number of boundary faces before thin wall, the additional
   * boundary faces will be added at the end of the list */
  cs_lnum_t n_b_faces_prev = m->n_b_faces;

  cs_mesh_boundary_insert_separating_cells(m,
                                           NULL, /* group_name */
                                           n_selected_cells,
                                           selected_cells);

  /* Prepare locator */

  cpl->n_local = m->n_b_faces - n_b_faces_prev; /* boundary faces added at end */

  BFT_MALLOC(cpl->faces_local, cpl->n_local, cs_lnum_t);
  BFT_MALLOC(cpl->c_tag, cpl->n_local, int);

  /* Assign separate cell tags on each side for location */

  cs_lnum_t *cell_tag = NULL;

  BFT_MALLOC(cell_tag, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
    cell_tag[cell_id] = 0;

  for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++)
    cell_tag[selected_cells[ii]] = 1;

  BFT_FREE(selected_cells);

  /* Loop over new boundary faces */

  cs_lnum_t ii = 0;
  for (cs_lnum_t face_id = n_b_faces_prev; face_id < m->n_b_faces ; face_id++) {
    cs_lnum_t cell_id = m->b_face_cells[face_id];
    cpl->c_tag[ii] = cell_tag[cell_id];
    cpl->faces_local[ii] = face_id;
    ii++;
  }

  BFT_FREE(cell_tag);
}

/*----------------------------------------------------------------------------
 * Initialize internal coupling locators using cell and face selection criteria.
 *
 * parameters:
 *   m   <->  pointer to mesh structure to modify
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

static void
_volume_face_initialize(cs_mesh_t               *m,
                        cs_internal_coupling_t  *cpl)
{
  cs_lnum_t  n_selected_cells;
  cs_lnum_t *selected_faces = NULL;
  cs_lnum_t *selected_cells = NULL;

  cs_lnum_t *cell_tag = NULL;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  /* Selection of Volme zone using volumic selection criteria*/

  BFT_MALLOC(selected_cells, n_cells_ext, cs_lnum_t);
  cs_selector_get_cell_list(cpl->cells_criteria,
                            &n_selected_cells,
                            selected_cells);

  /* Initialization */
  BFT_MALLOC(cell_tag, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    cell_tag[cell_id] = 2;
  }

  /* Tag cells */
  for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++) {
    cs_lnum_t cell_id = selected_cells[ii];
    cell_tag[cell_id] = 1;
  }
  if (cs_glob_mesh->halo != NULL)
    cs_halo_sync_num(cs_glob_mesh->halo, CS_HALO_STANDARD, cell_tag);

  /* Free memory */
  BFT_FREE(selected_cells);

  /* Selection of the juncture */
  cs_lnum_t  n_selected_faces = 0;

  BFT_MALLOC(selected_faces, m->n_b_faces, cs_lnum_t);
  cs_selector_get_b_face_list(cpl->faces_criteria,
                              &n_selected_faces,
                              selected_faces);

  /* Prepare locator */

  cpl->n_local = n_selected_faces; /* WARNING: only valid for conformal meshes */

  BFT_MALLOC(cpl->faces_local, cpl->n_local, cs_lnum_t);
  BFT_MALLOC(cpl->c_tag, cpl->n_local, int);

  for (cs_lnum_t ii = 0; ii < cpl->n_local; ii++) {
    cs_lnum_t face_id = selected_faces[ii];
    cpl->faces_local[ii] = face_id;
    cs_lnum_t cell_id = m->b_face_cells[face_id];
    cpl->c_tag[ii] = cell_tag[cell_id];
  }

  /* Free memory */
  BFT_FREE(selected_faces);
  BFT_FREE(cell_tag);
}

/*----------------------------------------------------------------------------
 * Initialize locators
 *
 * parameters:
 *   m              <->  pointer to mesh structure to modify
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

static void
_locator_initialize(cs_mesh_t               *m,
                    cs_internal_coupling_t  *cpl)
{
  /* Initialize locator */

  cpl->locator = _create_locator(cpl);
  cpl->n_distant = ple_locator_get_n_dist_points(cpl->locator);
  BFT_MALLOC(cpl->faces_distant,
             cpl->n_distant,
             cs_lnum_t);
  const cs_lnum_t *faces_distant_num
    = ple_locator_get_dist_locations(cpl->locator);

  /* From 1..n to 0..n-1 */
  for (cs_lnum_t i = 0; i < cpl->n_distant; i++)
    cpl->faces_distant[i] = faces_distant_num[i] - 1;


  /* Geometric quantities */

  BFT_MALLOC(cpl->g_weight, cpl->n_local, cs_real_t);
  BFT_MALLOC(cpl->ci_cj_vect, cpl->n_local, cs_real_3_t);
  BFT_MALLOC(cpl->offset_vect, cpl->n_local, cs_real_3_t);

  _compute_ci_cj_vect(cpl);

  /* Allocate coupling exchange coefficients */

  BFT_MALLOC(cpl->h_int, cpl->n_local, cs_real_t);
  BFT_MALLOC(cpl->h_ext, cpl->n_local, cs_real_t);

  BFT_MALLOC(cpl->coupled_faces, m->n_b_faces, bool);

  cpl->cocgb_s_lsq = NULL;
  cpl->cocgb_s_it = NULL;
  cpl->cocg_s_it = NULL;
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

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;
  const cs_real_t* g_weight = cpl->g_weight;
  cs_real_t *r_weight = NULL;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_mesh_quantities_t* fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;

  /* Exchange pvar */

  cs_real_t *pvar_distant = NULL;
  BFT_MALLOC(pvar_distant, n_distant, cs_real_t);
  for (ii = 0; ii < n_distant; ii++) {
    face_id = faces_distant[ii];
    cell_id = b_face_cells[face_id];
    pvar_distant[ii] = pvar[cell_id];
  }
  cs_real_t *pvar_local = NULL;
  BFT_MALLOC(pvar_local, n_local, cs_real_t);
  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    pvar_distant,
                                    pvar_local);
  /* Free memory */
  BFT_FREE(pvar_distant);

  /* Preliminary step in case of heterogenous diffusivity */

  if (c_weight != NULL) {
    BFT_MALLOC(r_weight, n_local, cs_real_t);
    _compute_physical_face_weight(cpl,
                                  c_weight, /* diffusivity */
                                  r_weight); /* physical face weight */
    /* Redefinition of rweight :
         Before : (1-g_weight)*rweight <==> 1 - ktpond
         Modif : rweight = ktpond
       Scope of this modification is local */
    for (ii = 0; ii < n_local; ii++)
      r_weight[ii] = 1.0 - (1.0-g_weight[ii]) * r_weight[ii];
  }

  /* Add contribution */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];

    /* compared to _initialize_scalar_gradient :
       1 - rweight <==> 1 - ktpond
       g_weight <==> weight = alpha_ij
       pvar[cell_id] <==> pvar[ii]
       pvar_local[ii] <==> pvar[jj]
       b_f_face_normal <==> i_f_face_normal */
    cs_real_t pfaci = (c_weight == NULL) ?
      (1.0-g_weight[ii]) * (pvar_local[ii] - pvar[cell_id]) :
      (1.0-r_weight[ii]) * (pvar_local[ii] - pvar[cell_id]);

    for (int j = 0; j < 3; j++) {
      grad[cell_id][j] += pfaci * b_f_face_normal[face_id][j];
    }

  }
  /* Free memory */
  BFT_FREE(r_weight);
  BFT_FREE(pvar_local);
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

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;
  const cs_real_3_t *offset_vect = (const cs_real_3_t *)cpl->offset_vect;
  const cs_real_t* g_weight = cpl->g_weight;
  cs_real_t *r_weight = NULL;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_mesh_quantities_t* fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;

  /* Exchange grad and pvar */

  cs_real_t *grad_distant = NULL;
  BFT_MALLOC(grad_distant, 3*n_distant, cs_real_t);
  cs_real_t *pvar_distant = NULL;
  BFT_MALLOC(pvar_distant, n_distant, cs_real_t);
  for (ii = 0; ii < n_distant; ii++) {
    face_id = faces_distant[ii];
    cell_id = b_face_cells[face_id];
    pvar_distant[ii] = pvar[cell_id];
    for (ll = 0; ll < 3; ll++) {
      grad_distant[3*ii+ll] = grad[cell_id][ll];
    }
  }
  cs_real_t *grad_local = NULL;
  BFT_MALLOC(grad_local, 3*n_local, cs_real_t);
  cs_real_t *pvar_local = NULL;
  BFT_MALLOC(pvar_local, n_local, cs_real_t);
  cs_internal_coupling_exchange_var(cpl,
                                    3,
                                    grad_distant,
                                    grad_local);
  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    pvar_distant,
                                    pvar_local);
  /* Free memory */
  BFT_FREE(grad_distant);
  BFT_FREE(pvar_distant);

  /* Preliminary step in case of heterogenous diffusivity */

  if (c_weight != NULL) { /* Heterogenous diffusivity */
    BFT_MALLOC(r_weight, n_local, cs_real_t);
    _compute_physical_face_weight(cpl,
                                  c_weight, /* diffusivity */
                                  r_weight); /* physical face weight */
    /* Redefinition of rweight_* :
         Before : (1-g_weight)*rweight <==> 1 - ktpond
         Modif : rweight = ktpond
       Scope of this modification is local */
    for (ii = 0; ii < n_local; ii++)
      r_weight[ii] = 1.0 - (1.0-g_weight[ii]) * r_weight[ii];
  }

  /* Compute rhs */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
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
         g_weight <==> weight = alpha_ij
         pvar[cell_id] <==> pvar[cell_id1]
         pvar_local[ii] <==> pvar[cell_id2]
         b_f_face_normal <==> i_f_face_normal */
    cs_real_t pfaci = 0.5;
    pfaci *= offset_vect[ii][0]*(grad_local[3*ii  ]+grad[cell_id][0])
            +offset_vect[ii][1]*(grad_local[3*ii+1]+grad[cell_id][1])
            +offset_vect[ii][2]*(grad_local[3*ii+2]+grad[cell_id][2]);
    if (c_weight != NULL) {
      pfaci += (1.0-r_weight[ii]) * (pvar_local[ii] - pvar[cell_id]);
    } else {
      pfaci += (1.0-g_weight[ii]) * (pvar_local[ii] - pvar[cell_id]);
    }

    for (int j = 0; j < 3; j++) {
      rhs[cell_id][j] += pfaci * b_f_face_normal[face_id][j];
    }

  }

  BFT_FREE(r_weight);
  BFT_FREE(grad_local);
  BFT_FREE(pvar_local);
}

/*----------------------------------------------------------------------------
 * Add internal coupling rhs contribution for LSQ gradient calculation
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
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
  cs_real_t pfac;
  cs_real_3_t dc, fctb;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;
  const cs_real_3_t *ci_cj_vect = (const cs_real_3_t *)cpl->ci_cj_vect;
  cs_real_t *r_weight;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Exchange pvar stored in rhsv[][3] */

  cs_real_t *pvar_distant = NULL;
  BFT_MALLOC(pvar_distant, n_distant, cs_real_t);
  for (ii = 0; ii < n_distant; ii++) {
    face_id = faces_distant[ii];
    cell_id = b_face_cells[face_id];
    pvar_distant[ii] = rhsv[cell_id][3];
  }
  cs_real_t *pvar_local = NULL;
  BFT_MALLOC(pvar_local, n_local, cs_real_t);
  cs_internal_coupling_exchange_var(cpl,
                                    1,
                                    pvar_distant,
                                    pvar_local);
  /* Free memory */
  BFT_FREE(pvar_distant);

  /* Preliminary step in case of heterogenous diffusivity */

  if (c_weight != NULL) { /* Heterogenous diffusivity */
    BFT_MALLOC(r_weight, n_local, cs_real_t);
    _compute_physical_face_weight(cpl,
                                  c_weight, /* diffusivity */
                                  r_weight); /* physical face weight */
  }

  /* Compute rhs */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];
    for (ll = 0; ll < 3; ll++)
      dc[ll] = ci_cj_vect[ii][ll];

    pfac = (pvar_local[ii] - rhsv[cell_id][3])
         / (dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
    if (c_weight != NULL) pfac *= r_weight[ii];

    for (ll = 0; ll < 3; ll++)
      fctb[ll] = dc[ll] * pfac;

    for (ll = 0; ll < 3; ll++)
      rhsv[cell_id][ll] += fctb[ll];
  }
  /* Free memory */
  BFT_FREE(r_weight);
  BFT_FREE(pvar_local);
}

/*----------------------------------------------------------------------------
 * Modify LSQ COCG matrix to include internal coupling
 *
 * parameters:
 *   cpl  <-- pointer to coupling entity
 *   cocg <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                           cs_real_33_t                   cocg[])
{
  int ii, ll, mm;
  cs_lnum_t face_id, cell_id;
  cs_real_3_t dddij;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_real_3_t *ci_cj_vect = (const cs_real_3_t *)cpl->ci_cj_vect;

  const cs_mesh_t* m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];
    for (ll = 0; ll < 3; ll++)
      dddij[ll] = ci_cj_vect[ii][ll];

    cs_real_t umdddij = 1./ cs_math_3_norm(dddij);
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
 *   cpl  <-- pointer to coupling entity
 *   cocg <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_it_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                          cs_real_33_t                   cocg[])
{
  int ii, ll, mm;
  cs_lnum_t cell_id, face_id;
  cs_real_t dvol;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_real_3_t *offset_vect = (const cs_real_3_t *)cpl->offset_vect;

  const cs_mesh_t* m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_mesh_quantities_t* fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];

    for (ll = 0; ll < 3; ll++) {
      for (mm = 0; mm < 3; mm++) {
        cocg[cell_id][ll][mm] -= 0.5 * offset_vect[ii][ll] * b_f_face_normal[face_id][mm];
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
 * Destruction of all internal coupling related structures.
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_finalize(void)
{
  cs_internal_coupling_t* cpl;
  for (int ii = 0; ii < _n_internal_couplings; ii++) {
    cpl = _internal_coupling + ii;
    _destroy_entity(cpl);
  }
  BFT_FREE(_internal_coupling);
  _n_internal_couplings = 0;
}

/*----------------------------------------------------------------------------
 * Exchange variable between groups using face id
 *
 * parameters:
 *   cpl    <-- pointer to coupling entity
 *   stride <-- number of values (non interlaced) by entity
 *   tab    <-- variable exchanged
 *   local  --> local data
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_face_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local[])
{
  cs_lnum_t face_id;

  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;

  /* Initialize distant array */

  cs_real_t *distant = NULL;
  BFT_MALLOC(distant, n_distant*stride, cs_real_t);
  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    face_id = faces_distant[ii];
    for (int jj = 0; jj < stride; jj++) {
      distant[stride * ii + jj] = tab[stride * face_id + jj];
    }
  }

  /* Exchange variable */

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    distant,
                                    local);
  /* Free memory */
  BFT_FREE(distant);
}

/*----------------------------------------------------------------------------
 * Exchange variable between groups using cell id
 *
 * parameters:
 *   cpl    <-- pointer to coupling entity
 *   stride <-- number of values (non interlaced) by entity
 *   tab    <-- variable exchanged
 *   local  --> local data
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_cell_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local[])
{
  int ii, jj;
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;

  const cs_mesh_t* m = cs_glob_mesh;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Initialize distant array */

  cs_real_t *distant = NULL;
  BFT_MALLOC(distant, n_distant*stride, cs_real_t);
  for (ii = 0; ii < n_distant; ii++) {
    face_id = faces_distant[ii];
    cell_id = b_face_cells[face_id];
    for (jj = 0; jj < stride; jj++) {
      distant[stride * ii + jj] = tab[stride * cell_id + jj];
    }
  }

  /* Exchange variable */

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    distant,
                                    local);
  /* Free memory */
  BFT_FREE(distant);
}

/*----------------------------------------------------------------------------
 * Exchange quantities from distant to local (update local using distant)
 *
 * parameters:
 *   cpl     <-- pointer to coupling entity
 *   stride  <-- Stride (e.g. 1 for double, 3 for interleaved coordinates)
 *   distant <-- Distant values, size coupling->n_distant
 *   local   --> Local values, size coupling->n_local
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_var(const cs_internal_coupling_t  *cpl,
                                  int                            stride,
                                  cs_real_t                      distant[],
                                  cs_real_t                      local[])
{
  ple_locator_exchange_point_var(cpl->locator,
                                 distant,
                                 local,
                                 NULL,
                                 sizeof(cs_real_t),
                                 stride,
                                 0);
}

/*----------------------------------------------------------------------------
 * Return pointers to coupling components
 *
 * parameters:
 *   cpl             <-- pointer to coupling entity
 *   n_local         --> NULL or pointer to component n_local
 *   faces_local     --> NULL or pointer to component faces_local
 *   n_distant       --> NULL or pointer to component n_distant
 *   faces_distant   --> NULL or pointer to component faces_distant
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_coupled_faces(const cs_internal_coupling_t  *cpl,
                                   cs_lnum_t                     *n_local,
                                   cs_lnum_t                     *faces_local[],
                                   cs_lnum_t                     *n_distant,
                                   cs_lnum_t                     *faces_distant[])
{
  if (n_local != NULL) {
    *n_local = cpl->n_local;
  }
  if (faces_local != NULL) {
    *faces_local = cpl->faces_local;
  }
  if (n_distant != NULL) {
    *n_distant = cpl->n_distant;
  }
  if (faces_distant != NULL) {
    *faces_distant = cpl->faces_distant;
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
  int ii;
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_glob_mesh->b_face_cells;
  const cs_internal_coupling_t* cpl = (const cs_internal_coupling_t *)input;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  const cs_real_t thetap = cpl->thetav;
  const int       idiffp = cpl->idiff;

  /* Exchange x */

  cs_real_t *x_j = NULL;
  BFT_MALLOC(x_j, n_local, cs_real_t);
  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           1,
                                           x,
                                           x_j);

  /* Compute heq and update y */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];

    cs_real_t fluxi = 0.;
    cs_real_t pi = exclude_diag ?
      0. : x[cell_id]; /* If exclude_diag, no diagonal term */
    cs_real_t pj = x_j[ii];

    cs_real_t hint = cpl->h_int[ii];
    cs_real_t hext = cpl->h_ext[ii];
    cs_real_t heq = hint * hext / (hint + hext);

    cs_b_diff_flux_coupling(idiffp,
                            pi,
                            pj,
                            heq,
                            &fluxi);

    y[cell_id] += thetap * fluxi;
  }
  /* Free memory */
  BFT_FREE(x_j);
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

  for (int i = 0; i < _n_internal_couplings; i++) {
    cs_internal_coupling_t *cpl = _internal_coupling + i;
    _locator_initialize(cs_glob_mesh, cpl);
  }

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
        cs_internal_coupling_t *cpl
          = _internal_coupling + coupling_id;

        /* Definition of var_cal_opt options
         * (needed for matrix.vector multiply) */
        cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
        cpl->thetav = var_cal_opt.thetav;
        cpl->idiff = var_cal_opt.idiff;

        /* Initialize coupled_faces */
        _initialize_coupled_faces(cpl);

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
 * Log information about a given internal coupling entity
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_log(const cs_internal_coupling_t  *cpl)
{
  if (cpl == NULL)
    return;

  const cs_lnum_t n_local = cpl->n_local;

  if (cpl->faces_criteria == NULL)
    bft_printf("  Coupled scalar: %s\n"
               "  Cell group selection criterion: %s\n"
               "  Locator: - n dist points = %d\n"
               "           - n interior = %d\n"
               "           - n exterior = %d\n",
               cpl->namesca,
               cpl->cells_criteria,
               ple_locator_get_n_dist_points(cpl->locator),
               ple_locator_get_n_interior(cpl->locator),
               ple_locator_get_n_exterior(cpl->locator));
  else
    bft_printf("  Coupled scalar: %s\n"
               "  Cell group selection criterion: %s\n"
               "  Face group selection criterion: %s\n"
               "  Locator: - n dist points = %d\n"
               "           - n interior = %d\n"
               "           - n exterior = %d\n",
               cpl->namesca,
               cpl->cells_criteria,
               cpl->faces_criteria,
               ple_locator_get_n_dist_points(cpl->locator),
               ple_locator_get_n_interior(cpl->locator),
               ple_locator_get_n_exterior(cpl->locator));

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

  if (_n_internal_couplings == 0)
    return;

  bft_printf("\n Internal coupling\n");
  for (ii = 0; ii < _n_internal_couplings; ii++) {
    cpl = _internal_coupling + ii;
    bft_printf("coupling_id = %d\n", ii);
    cs_internal_coupling_log(cpl);
  }
}

/*----------------------------------------------------------------------------
 * Define coupling volume using given criterias. Then, this volume will be
 * seperated from the rest of the domain with thin walls.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure to modify
 *   criteria_cells <-- selection criteria for the first group of cells
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_volume(cs_mesh_t   *mesh,
                                const char   criteria_cells[])
{
  BFT_REALLOC(_internal_coupling,
              _n_internal_couplings + 1,
              cs_internal_coupling_t);

  cs_internal_coupling_t *cpl = _internal_coupling + _n_internal_couplings;

  cpl->locator = NULL;
  cpl->c_tag = NULL;
  cpl->cells_criteria = NULL;
  cpl->faces_criteria = NULL; /* optional here */

  cpl->n_local = 0;
  cpl->faces_local = NULL; /* Coupling boundary faces, numbered 0..n-1 */

  cpl->n_distant = 0;
  cpl->faces_distant = NULL;

  cpl->coupled_faces = NULL;

  cpl->h_int = NULL;
  cpl->h_ext = NULL;

  cpl->g_weight = NULL;
  cpl->ci_cj_vect = NULL;
  cpl->offset_vect = NULL;

  cpl->thetav = 0;
  cpl->idiff = 0;

  cpl->cocgb_s_lsq = NULL;
  cpl->cocgb_s_it = NULL;
  cpl->cocg_s_it = NULL;

  cpl->namesca = NULL;

  _criteria_initialize(criteria_cells, NULL, cpl);

  /* Initialization of locators */
  _volume_initialize(mesh, cpl);

  _n_internal_couplings++;
}

/*----------------------------------------------------------------------------
 * Define coupling volume using given criterias. Then, this volume has tobe
 * seperated from the rest of the domain with a wall.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure to modify
 *   criteria_cells <-- selection criteria for the first group of cells
 *   criteria_faces <-- selection criteria for faces to be joint
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_add(cs_mesh_t   *mesh,
                         const char   criteria_cells[],
                         const char   criteria_faces[])
{
  BFT_REALLOC(_internal_coupling,
              _n_internal_couplings + 1,
              cs_internal_coupling_t);

  cs_internal_coupling_t *cpl = _internal_coupling + _n_internal_couplings;

  cpl->locator = NULL;
  cpl->c_tag = NULL;
  cpl->cells_criteria = NULL;
  cpl->faces_criteria = NULL;

  cpl->n_local = 0;
  cpl->faces_local = NULL; /* Coupling boundary faces, numbered 0..n-1 */

  cpl->n_distant = 0;
  cpl->faces_distant = NULL;

  cpl->coupled_faces = NULL;

  cpl->h_int = NULL;
  cpl->h_ext = NULL;

  cpl->g_weight = NULL;
  cpl->ci_cj_vect = NULL;
  cpl->offset_vect = NULL;

  cpl->thetav = 0;
  cpl->idiff = 0;

  cpl->cocgb_s_lsq = NULL;
  cpl->cocgb_s_it = NULL;
  cpl->cocg_s_it = NULL;

  cpl->namesca = NULL;

  _criteria_initialize(criteria_cells, criteria_faces, cpl);

  /* Initialization of locators */
  _volume_face_initialize(mesh, cpl);

  _n_internal_couplings++;
}


/*----------------------------------------------------------------------------
 * Define coupling entity using given criterias
 *
 * parameters:
 *   f_id       <-- id of the field
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_entity(int        f_id)
{
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  cs_field_t* f = cs_field_by_id(f_id);

  if (f->type & CS_FIELD_VARIABLE) {
    cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.icoupl = 1;
    cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              "field id = %d provided is invalid."
              " The field must be a variable.",
              f_id);
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

  cs_real_t surf;
  const cs_real_t* b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  const cs_field_t* f = cs_field_by_id(field_id);
  const cs_int_t coupling_key_id = cs_field_key_id("coupling_entity");
  int coupling_id = cs_field_get_key_int(f, coupling_key_id);
  const cs_internal_coupling_t  *cpl
    = cs_internal_coupling_by_id(coupling_id);

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  cs_real_t *h_int = cpl->h_int;
  cs_real_t *h_ext = cpl->h_ext;

  /* Exchange hbord */

  cs_internal_coupling_exchange_by_face_id(cpl,
                                           1,
                                           hbord,
                                           h_ext);

  /* Compute hint and hext */

  for (ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    surf = b_face_surf[face_id];
    h_int[ii] = hbord[face_id] * surf;
    h_ext[ii] *= surf;
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
    cs_lnum_t face_id, cell_id;

    const cs_lnum_t n_local = cpl->n_local;
    const cs_lnum_t *faces_local = cpl->faces_local;

    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *restrict)cs_glob_mesh->b_face_cells;

    for (ii = 0; ii < n_local; ii++) {
      face_id = faces_local[ii];
      cell_id = b_face_cells[face_id];

      cs_real_t hint = cpl->h_int[ii];
      cs_real_t hext = cpl->h_ext[ii];
      cs_real_t heq = hint * hext / (hint + hext);

      ad[cell_id] += heq;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

