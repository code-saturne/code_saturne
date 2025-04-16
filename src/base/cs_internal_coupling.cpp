/*============================================================================
 * Internal coupling: coupling for one instance of code_saturne
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "bft/bft_printf.h"
#include "bft/bft_error.h"

#include "fvm/fvm_defs.h"
#include "fvm/fvm_selector.h"

#include "base/cs_dispatch.h"
#include "base/cs_defs.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_sort.h"
#include "base/cs_search.h"
#include "mesh/cs_mesh_connect.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_coupling.h"
#include "alge/cs_gradient_boundary.h"
#include "base/cs_halo.h"
#include "alge/cs_matrix.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_boundary.h"
#include "mesh/cs_mesh_group.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_selector.h"
#include "base/cs_parall.h"
#include "base/cs_prototypes.h"
#include "base/cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_internal_coupling.h"

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

static cs_internal_coupling_t  *_internal_coupling = nullptr;
static int                      _n_internal_couplings = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the equivalent heat transfer coefficient. If both terms are
 * below a given tolerance, 0. is returned.
 *
 * parameters:
 *   h1     <-- first exchange coefficient
 *   h2     <-- second exchange coefficient
 *
 * return:
 *   value of equivalent exchange coefficient
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_calc_heq(cs_real_t h1,
          cs_real_t h2)
{
  const cs_real_t h_eps = 1.e-12;

  cs_real_t heq = 0.;
  if (h1 + h2 > h_eps)
    heq = h1 * h2 / (h1 + h2);

  return heq;
}

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
  const int *c_tag = cpl->c_tag;

  fvm_nodal_t* nm = nullptr;
  int *tag_nm = nullptr;
  cs_lnum_t *faces_in_nm = nullptr;

  ple_coord_t* point_coords = nullptr;

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

  nm = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                      mesh_name,
                                      false,
                                      0,
                                      n_local,
                                      nullptr,
                                      cpl->faces_local);

  /* Tag fvm_nodal_t structure */

  /* Number of faces to tag */
  const int nfac_in_nm = fvm_nodal_get_n_entities(nm, 2);
  /* Memory allocation */
  CS_MALLOC(faces_in_nm, nfac_in_nm, cs_lnum_t);
  CS_MALLOC(tag_nm, nfac_in_nm, int);
  /* Get id of faces to tag in parent */
  fvm_nodal_get_parent_num(nm, 2, faces_in_nm);
  /* Tag faces */
  for (cs_lnum_t ii = 0; ii < nfac_in_nm; ii++) {
    /* Default tag is 0 */
    tag_nm[ii] = 0;
    for (cs_lnum_t jj = 0; jj < n_local; jj++) {
      if (faces_in_nm[ii] == cpl->faces_local[jj] + 1) {
        tag_nm[ii] = c_tag[jj];
        break;
      }
    }
  }
  fvm_nodal_set_tag(nm, tag_nm, 2);
  /* Free memory */
  CS_FREE(faces_in_nm);
  CS_FREE(tag_nm);

  /* Creation of distant group cell centers */

  CS_MALLOC(point_coords, 3*n_local, cs_real_t);

  for (i = 0; i < n_local; i++) {
    ifac = cpl->faces_local[i]; /* 0..n-1 */
    for (j = 0; j < 3; j++)
      point_coords[3*i+j] = cs_glob_mesh_quantities->b_face_cog[ifac][j];
  }

  /* Locator initialization */

  ple_locator_set_mesh(locator,
                       nm,
                       nullptr,
                       0,
                       1.1, /* TODO */
                       cs_glob_mesh->dim,
                       n_local,
                       nullptr,
                       c_tag,
                       point_coords,
                       nullptr,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);
  /* Free memory */
  nm = fvm_nodal_destroy(nm);
  CS_FREE(point_coords);
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
  CS_FREE(cpl->c_tag);
  CS_FREE(cpl->faces_local);
  CS_FREE(cpl->faces_distant);
  CS_FREE(cpl->ci_cj_vect);
  CS_FREE(cpl->coupled_faces);
  CS_FREE(cpl->cells_criteria);
  CS_FREE(cpl->faces_criteria);
  CS_FREE(cpl->interior_faces_group_name);
  CS_FREE(cpl->exterior_faces_group_name);
  CS_FREE(cpl->volume_zone_ids);
  ple_locator_destroy(cpl->locator);
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

  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  /* Exchange cell center location */

  cs_real_t *cell_cen_local = nullptr;
  CS_MALLOC(cell_cen_local, 3 * n_local, cs_real_t);
  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           3, /* dimension */
                                           (const cs_real_t *)mq->cell_cen,
                                           cell_cen_local);

  /* Compute IJ vectors */

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    face_id = faces_local[ii];
    cell_id = b_face_cells[face_id];

    for (cs_lnum_t jj = 0; jj < 3; jj++) {
      cs_real_t xxd = cell_cen_local[3*ii + jj];
      cs_real_t xxl = cell_cen[cell_id][jj];
      ci_cj_vect[ii][jj] = xxd - xxl;
    }
  }

  /* Free memory */
  CS_FREE(cell_cen_local);
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
 * Initialize to 0 or nullptr most of the fields in given coupling structure.
 *
 * parameters:
 *   cpl               --> pointer to coupling structure to initialize
 *----------------------------------------------------------------------------*/

static void
_cpl_initialize(cs_internal_coupling_t *cpl)
{
  cpl->locator = nullptr;
  cpl->c_tag = nullptr;
  cpl->cells_criteria = nullptr;
  cpl->faces_criteria = nullptr;
  cpl->interior_faces_group_name = nullptr;
  cpl->exterior_faces_group_name = nullptr;

  cpl->n_volume_zones = 0;
  cpl->volume_zone_ids = nullptr;

  cpl->n_local = 0;
  cpl->faces_local = nullptr; /* Coupling boundary faces, numbered 0..n-1 */

  cpl->n_distant = 0;
  cpl->faces_distant = nullptr;

  cpl->coupled_faces = nullptr;

  cpl->ci_cj_vect = nullptr;
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
_criteria_initialize(const char               criteria_cells[],
                     const char               criteria_faces[],
                     cs_internal_coupling_t  *cpl)
{
  CS_MALLOC(cpl->cells_criteria, strlen(criteria_cells)+1, char);
  strcpy(cpl->cells_criteria, criteria_cells);

  if (criteria_faces != nullptr) {
    CS_MALLOC(cpl->faces_criteria, strlen(criteria_faces)+1, char);
    strcpy(cpl->faces_criteria, criteria_faces);
  }
}

/*----------------------------------------------------------------------------
 * Define face to face mappings for internal couplings.
 *
 * The caller is responsible for freeing the list.
 *
 * parameters:
 *   cpl          <->  pointer to internal coupling structure
 *   coupling_id  <--  associated coupling id
 *----------------------------------------------------------------------------*/

static void
_auto_group_name(cs_internal_coupling_t  *cpl,
                 int                      coupling_id)
{
  char group_name[64];
  snprintf(group_name, 63, "auto:internal_coupling_%d", coupling_id);
  group_name[63] = '\0';
  CS_REALLOC(cpl->faces_criteria,
             strlen(group_name)+1,
             char);
  strcpy(cpl->faces_criteria, group_name);
}

/*----------------------------------------------------------------------------
 * Get selected cells list.
 *
 * parameters:
 *   m         <--  pointer to mesh structure to modify
 *   cpl       <-- pointer to coupling structure to modify
 *   n_cells   --> number of associated cells
 *   cell_list --> associated cells list
 *----------------------------------------------------------------------------*/

static void
_get_cell_list(cs_mesh_t               *m,
               cs_internal_coupling_t  *cpl,
               cs_lnum_t               *n_cells,
               cs_lnum_t              **cell_list)
{
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells = nullptr;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  CS_MALLOC(selected_cells, n_cells_ext, cs_lnum_t);

  if (cpl->cells_criteria != nullptr) {
    cs_selector_get_cell_list(cpl->cells_criteria,
                              &n_selected_cells,
                              selected_cells);
  }

  /* For zone selection, zones may not be built yet, so use selection
     mechanism directly. */

  else if (cpl->n_volume_zones > 0) {

    int *cell_flag;
    CS_MALLOC(cell_flag, n_cells_ext, int);
    for (cs_lnum_t i = 0; i < n_cells_ext; i++)
      cell_flag[i] = 0;

    for (int i = 0; i < cpl->n_volume_zones; i++) {
      const cs_zone_t *z = cs_volume_zone_by_id(cpl->volume_zone_ids[i]);
      const char *criteria
        = cs_mesh_location_get_selection_string(z->location_id);

      if (criteria == nullptr)
        bft_error
          (__FILE__, __LINE__, 0,
           _("Only zones based on selection criteria strings "
             "(not functions) are currently\n"
             "supperted for the selection of internal coupling volumes.\n\n"
             "This is not the case for zone: \"%s\"."), z->name);

      cs_selector_get_cell_list(criteria, &n_selected_cells, selected_cells);

      for (cs_lnum_t j = 0; j < n_selected_cells; j++)
        cell_flag[selected_cells[j]] = 1;
    }

    n_selected_cells = 0;
    for (cs_lnum_t i = 0; i < m->n_cells; i++) {
      if (cell_flag[i] == 1) {
        selected_cells[n_selected_cells] = i;
        n_selected_cells++;
      }
    }

    CS_FREE(cell_flag);
  }

  CS_REALLOC(selected_cells, n_selected_cells, cs_lnum_t);

  *n_cells = n_selected_cells;
  *cell_list = selected_cells;
}

/*----------------------------------------------------------------------------
 * Initialize internal coupling and insert boundaries
 * using cell selection criteria ONLY.
 *
 * parameters:
 *   m   <->  pointer to mesh structure to modify
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

static void
_volume_initialize_insert_boundary(cs_mesh_t               *m,
                                   cs_internal_coupling_t  *cpl)
{
  cs_lnum_t  n_sel_cells = 0;
  cs_lnum_t *sel_cells = nullptr;

  /* Selection of Volume zone using volumic selection criteria*/

  _get_cell_list(m,
                 cpl,
                 &n_sel_cells,
                 &sel_cells);

  int coupling_id = _n_internal_couplings - 1;

  _auto_group_name(cpl, coupling_id);

  cs_mesh_boundary_insert_separating_cells(m,
                                           cpl->faces_criteria,
                                           n_sel_cells,
                                           sel_cells);

  /* Select faces adjacent to volume zone and add appropriate group
     so as to be able to easily extract separate sides */

  {
    cs_lnum_t  n_sel_faces = 0;
    cs_lnum_t *sel_faces_ext = nullptr, *sel_faces_int = nullptr;
    int *cell_flag;

    CS_MALLOC(cell_flag, m->n_cells, int);
    for (cs_lnum_t i = 0; i < m->n_cells; i++)
      cell_flag[i] = 0;

    for (cs_lnum_t i = 0; i < n_sel_cells; i++)
      cell_flag[sel_cells[i]] = 1;

    CS_MALLOC(sel_faces_ext, m->n_b_faces, cs_lnum_t);
    cs_selector_get_b_face_list(cpl->faces_criteria,
                                &n_sel_faces,
                                sel_faces_ext);

    cs_lnum_t n_sel_int = 0, n_sel_ext = 0;
    CS_MALLOC(sel_faces_int, n_sel_faces, cs_lnum_t);

    for (cs_lnum_t i = 0; i < n_sel_faces; i++) {
      cs_lnum_t face_id = sel_faces_ext[i];
      if (cell_flag[m->b_face_cells[face_id]]) {
        sel_faces_ext[n_sel_ext] = face_id;
        n_sel_ext++;
      }
      else {
        sel_faces_int[n_sel_int] = face_id;
        n_sel_int++;
      }
    }

    CS_FREE(cell_flag);

    if (cpl->exterior_faces_group_name != nullptr) {
      cs_mesh_group_b_faces_add(m,
                                cpl->exterior_faces_group_name,
                                n_sel_ext,
                                sel_faces_ext);
    }

    if (cpl->interior_faces_group_name != nullptr) {
      cs_mesh_group_b_faces_add(m,
                                cpl->interior_faces_group_name,
                                n_sel_int,
                                sel_faces_int);
    }

    CS_FREE(sel_faces_int);
    CS_FREE(sel_faces_ext);
  }

  CS_FREE(sel_cells);
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
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_faces = nullptr;
  cs_lnum_t *selected_cells = nullptr;

  cs_lnum_t *cell_tag = nullptr;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  /* Selection of Volume zone using selection criteria */

  _get_cell_list(m,
                 cpl,
                 &n_selected_cells,
                 &selected_cells);

  /* Initialization */

  CS_MALLOC(cell_tag, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
    cell_tag[cell_id] = 2;

  /* Tag cells */

  for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++) {
    cs_lnum_t cell_id = selected_cells[ii];
    cell_tag[cell_id] = 1;
  }
  if (cs_glob_mesh->halo != nullptr)
    cs_halo_sync_num(cs_glob_mesh->halo, CS_HALO_STANDARD, cell_tag);

  /* Free memory */

  CS_FREE(selected_cells);

  /* Selection of the interface */

  cs_lnum_t  n_selected_faces = 0;

  CS_MALLOC(selected_faces, m->n_b_faces, cs_lnum_t);
  cs_selector_get_b_face_list(cpl->faces_criteria,
                              &n_selected_faces,
                              selected_faces);

  /* reorder selected faces */
  {
    cs_lnum_t n = 0;
    char  *b_face_flag;
    CS_MALLOC(b_face_flag, m->n_b_faces, char);
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      b_face_flag[i] = 0;
    for (cs_lnum_t i = 0; i < n_selected_faces; i++)
      b_face_flag[selected_faces[i]] = 1;
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++) {
      if (b_face_flag[i] == 1)
        selected_faces[n++] = i;
    }
    assert(n == n_selected_faces);
    CS_FREE(b_face_flag);
  }

  /* Prepare locator */

  cpl->n_local = n_selected_faces; /* WARNING: only numerically
                                      valid for conformal meshes */

  CS_MALLOC(cpl->faces_local, cpl->n_local, cs_lnum_t);
  CS_MALLOC(cpl->c_tag, cpl->n_local, int);

  for (cs_lnum_t ii = 0; ii < cpl->n_local; ii++) {
    cs_lnum_t face_id = selected_faces[ii];
    cpl->faces_local[ii] = face_id;
    cs_lnum_t cell_id = m->b_face_cells[face_id];
    cpl->c_tag[ii] = cell_tag[cell_id];
  }

  /* Free memory */

  CS_FREE(selected_faces);
  CS_FREE(cell_tag);
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
  CS_MALLOC(cpl->faces_distant,
            cpl->n_distant,
            cs_lnum_t);
  const cs_lnum_t *faces_distant_num
    = ple_locator_get_dist_locations(cpl->locator);

  /* From 1..n to 0..n-1 */
  for (cs_lnum_t i = 0; i < cpl->n_distant; i++)
    cpl->faces_distant[i] = faces_distant_num[i] - 1;

  /* Geometric quantities */

  CS_MALLOC(cpl->ci_cj_vect, cpl->n_local, cs_real_3_t);

  _compute_ci_cj_vect(cpl);

  CS_MALLOC(cpl->coupled_faces, m->n_b_faces, bool);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of defined internal couplings.
 *
 * \return  number of internal couplings
 */
/*----------------------------------------------------------------------------*/

int
cs_internal_coupling_n_couplings(void)
{
  return _n_internal_couplings;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupling volume using given selection criteria.
 *
 * Then, this volume must be separated from the rest of the domain with a wall.
 *
 * \param[in]  criteria_cells  criteria for the first group of cells
 * \param[in]  criteria_faces  criteria for faces to be joined
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add(const char  criteria_cells[],
                         const char  criteria_faces[])
{
  CS_REALLOC(_internal_coupling,
             _n_internal_couplings + 1,
             cs_internal_coupling_t);

  cs_internal_coupling_t *cpl = _internal_coupling + _n_internal_couplings;

  cpl->id = _n_internal_couplings;

  _cpl_initialize(cpl);

  _criteria_initialize(criteria_cells, criteria_faces, cpl);

  _n_internal_couplings++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupling volume using given criteria. Then, this volume will
 * be separated from the rest of the domain with thin walls.
 *
 * \param[in]  criteria_cells  criteria for the first group of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_volume(const char  criteria_cells[])
{
  if (_n_internal_couplings > 0)
    bft_error(__FILE__, __LINE__, 0,
              "Only one volume can be added in this version.");

  CS_REALLOC(_internal_coupling,
             _n_internal_couplings + 1,
             cs_internal_coupling_t);

  cs_internal_coupling_t *cpl = _internal_coupling + _n_internal_couplings;

  cpl->id = _n_internal_couplings;

  _cpl_initialize(cpl);

  _criteria_initialize(criteria_cells, nullptr, cpl);

  _n_internal_couplings++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupling volume using given cs_zone_t. Then, this volume will
 * be separated from the rest of the domain with thin walls.
 *
 * \param[in]  z  pointer to cs_volume_zone_t
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_volume_zone(const cs_zone_t *z)
{
  int z_ids[] = {z->id};

  cs_internal_coupling_add_volume_zones(1, z_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupling volume using given cs_zone_t. Then, this volume will
 * be separated from the rest of the domain with thin walls.
 *
 * \param[in]  n_zones   number of associated volume zones
 * \param[in]  zone_ids  ids of associated volume zones
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_volume_zones(int        n_zones,
                                      const int  zone_ids[])
{
  if (_n_internal_couplings > 0)
    bft_error(__FILE__, __LINE__, 0,
              "Only one volume can be added in this version.");

  CS_REALLOC(_internal_coupling,
             _n_internal_couplings + 1,
             cs_internal_coupling_t);

  cs_internal_coupling_t *cpl = _internal_coupling + _n_internal_couplings;

  _cpl_initialize(cpl);

  cpl->id = _n_internal_couplings;

  cpl->n_volume_zones = n_zones;
  CS_MALLOC(cpl->volume_zone_ids, n_zones, int);

  for (int i = 0; i < n_zones; i++)
    cpl->volume_zone_ids[i] = zone_ids[i];

  _n_internal_couplings++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling volume boundary group names.
 *
 * This is used only for internal couplings based on a separation of volumes
 * (cs_internal_coupling_add_volume, cs_internal_coupling_add_volume_zone,
 * cs_internal_coupling_add_volume_zones).
 *
 * The interior name is used for faces adjacent to the main volume, and
 * the exterio name for faces adjacent to the selected (exterior) volume.
 *
 * This allows filtering faces on each side of the boundary in a simpler manner.
 *
 * \param[in, out] cpl             pointer to mesh structure to modify
 * \param[in]      criteria_cells  criteria for the first group of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_boundary_groups(cs_internal_coupling_t  *cpl,
                                         const char              *interior_name,
                                         const char              *exterior_name)
{
  if (cpl != nullptr && interior_name != nullptr) {
    CS_REALLOC(cpl->interior_faces_group_name, strlen(interior_name) + 1, char);
    strcpy(cpl->interior_faces_group_name, interior_name);
  }

  if (cpl != nullptr && exterior_name != nullptr) {
    CS_REALLOC(cpl->exterior_faces_group_name, strlen(exterior_name) + 1, char);
    strcpy(cpl->exterior_faces_group_name, exterior_name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Impose wall BCs to internal coupled faces if not yet defined.
 *
 * \param[in, out] bc_type       face boundary condition type
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_bcs(int         bc_type[])
{
  cs_internal_coupling_t *cpl;

  for (int cpl_id = 0; cpl_id < _n_internal_couplings; cpl_id++) {
    cpl = _internal_coupling + cpl_id;

    const cs_lnum_t n_local = cpl->n_local;
    const cs_lnum_t *faces_local = cpl->faces_local;

    for (cs_lnum_t ii = 0; ii < n_local; ii++) {
      cs_lnum_t face_id = faces_local[ii];
      if (bc_type[face_id] == 0)
        bc_type[face_id] = CS_SMOOTHWALL;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destruction of all internal coupling related structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_finalize(void)
{
  cs_internal_coupling_t* cpl;
  for (int cpl_id = 0; cpl_id < _n_internal_couplings; cpl_id++) {
    cpl = _internal_coupling + cpl_id;
    _destroy_entity(cpl);
  }
  CS_FREE(_internal_coupling);
  _n_internal_couplings = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the coupling associated with a given coupling_id.
 *
 * \param[in]  coupling_id  associated with a coupling entity
 *
 * \return pointer to associated coupling structure
 */
/*----------------------------------------------------------------------------*/

cs_internal_coupling_t *
cs_internal_coupling_by_id(int coupling_id)
{
  if (coupling_id > -1 && coupling_id < _n_internal_couplings)
    return _internal_coupling + coupling_id;
  else
    bft_error(__FILE__, __LINE__, 0,
              "coupling_id = %d provided is invalid", coupling_id);
  return (cs_internal_coupling_t*)nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange quantities from distant to local
 * (update local using distant).
 *
 * \param[in]  cpl     pointer to coupling entity
 * \param[in]  stride  stride (e.g. 1 for double, 3 for interleaved coordinates)
 * \param[in]  distant distant values, size coupling->n_distant
 * \param[out] local   local values, size coupling->n_local
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_var(const cs_internal_coupling_t  *cpl,
                                  int                            stride,
                                  cs_real_t                      distant[],
                                  cs_real_t                      local[])
{
  ple_locator_exchange_point_var(cpl->locator,
                                 distant,
                                 local,
                                 nullptr,
                                 sizeof(cs_real_t),
                                 stride,
                                 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange variable between groups using cell id.
 *
 * \param[in]  cpl     pointer to coupling entity
 * \param[in]  stride  number of values (non interlaced) by entity
 * \param[in]  tab     variable exchanged
 * \param[out] local   local data
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_cell_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local[])
{
  int jj;
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;

  const cs_mesh_t* m = cs_glob_mesh;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  /* Initialize distant array */

  cs_real_t *distant = nullptr;
  CS_MALLOC(distant, n_distant*stride, cs_real_t);
  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    face_id = faces_distant[ii];
    cell_id = b_face_cells[face_id];
    for (jj = 0; jj < stride; jj++)
      distant[stride * ii + jj] = tab[stride * cell_id + jj];
  }

  /* Exchange variable */

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    distant,
                                    local);
  /* Free memory */
  CS_FREE(distant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange variable between groups using face id.
 *
 * \param[in]  cpl     pointer to coupling entity
 * \param[in]  stride  number of values (non interlaced) by entity
 * \param[in]  tab     variable exchanged
 * \param[out] local   local data
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_face_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local[])
{
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;

  /* Initialize distant array */

  cs_real_t *distant = nullptr;
  CS_MALLOC(distant, n_distant*stride, cs_real_t);
  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    cs_lnum_t face_id = faces_distant[ii];
    for (cs_lnum_t jj = 0; jj < stride; jj++)
      distant[stride * ii + jj] = tab[stride * face_id + jj];
  }

  /* Exchange variable */

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    distant,
                                    local);
  /* Free memory */
  CS_FREE(distant);
}

/*----------------------------------------------------------------------------
 * Return pointers to coupling components
 *
 * parameters:
 *   cpl             <-- pointer to coupling entity
 *   n_local         --> null or pointer to component n_local
 *   faces_local     --> null or pointer to component faces_local
 *   n_distant       --> null or pointer to component n_distant
 *   faces_distant   --> null or pointer to component faces_distant
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_coupled_faces(const cs_internal_coupling_t  *cpl,
                                   cs_lnum_t                     *n_local,
                                   const cs_lnum_t               *faces_local[],
                                   cs_lnum_t                     *n_distant,
                                   const cs_lnum_t               *faces_distant[])
{
  if (n_local != nullptr)
    *n_local = cpl->n_local;
  if (faces_local != nullptr)
    *faces_local = cpl->faces_local;
  if (n_distant != nullptr)
    *n_distant = cpl->n_distant;
  if (faces_distant != nullptr)
    *faces_distant = cpl->faces_distant;
}

/*----------------------------------------------------------------------------
 * Addition to matrix-vector product in case of internal coupling.
 *
 * parameters:
 *   exclude_diag <-- extra diagonal flag
 *   f            <-- associated field pointer
 *   x            <-- vector x in m * x = y
 *   y            <-> vector y in m * x = y
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_spmv_contribution(bool               exclude_diag,
                                       const cs_field_t  *f,
                                       const cs_real_t   *restrict x,
                                       cs_real_t         *restrict y)
{
  const cs_real_t* b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  cs_lnum_t face_id, cell_id;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)cs_glob_mesh->b_face_cells;

  int coupling_id = cs_field_get_key_int(f,
                                         cs_field_key_id("coupling_entity"));
  const cs_internal_coupling_t *cpl
    = cs_internal_coupling_by_id(coupling_id);

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
  cs_real_t thetap = 0.0;
  int idiffp = 0;

  if (eqp->icoupl > 0) {
    thetap = eqp->theta;
    idiffp = eqp->idiff;
  }

  /* Exchange x */

  cs_real_t *x_j = nullptr;
  CS_MALLOC(x_j, f->dim * n_local, cs_real_t);
  cs_internal_coupling_exchange_by_cell_id(cpl,
                                           f->dim,
                                           x,
                                           x_j);

  /* Compute heq and update y */

  cs_real_t *hintp = f->bc_coeffs->hint;
  cs_real_t *rcodcl2p = f->bc_coeffs->rcodcl2;

  if (f->dim == 1) {
    for (cs_lnum_t ii = 0; ii < n_local; ii++) {
      face_id = faces_local[ii];
      cell_id = b_face_cells[face_id];
      cs_real_t surf = b_face_surf[face_id];

      cs_real_t pi = exclude_diag ?
        0. : x[cell_id]; /* If exclude_diag, no diagonal term */
      cs_real_t pj = x_j[ii];

      cs_real_t hint = hintp[face_id];
      cs_real_t rcodcl2 = rcodcl2p[face_id];
      cs_real_t heq = _calc_heq(hint, rcodcl2)*surf;

      y[cell_id] += thetap * idiffp * heq * (pi - pj);
    }

  } else if (f->dim == 3) {

    cs_real_3_t *_y = (cs_real_3_t *)y;
    const cs_real_3_t *_x = (const cs_real_3_t *)x;
    for (cs_lnum_t ii = 0; ii < n_local; ii++) {
      face_id = faces_local[ii];
      cell_id = b_face_cells[face_id];
      cs_real_t surf = b_face_surf[face_id];
      cs_real_t pi[3];
      /* If exclude_diag, no diagonal term */
      if (exclude_diag) {
        for (cs_lnum_t k = 0; k < 3; k++)
          pi[k] = 0.;
      } else {
        for (cs_lnum_t k = 0; k < 3; k++)
          pi[k] = _x[cell_id][k];
      }
      cs_real_t pj[3] = {x_j[ii], x_j[ii+1], x_j[ii+2]};

      cs_real_t hint = hintp[face_id];
      cs_real_t rcodcl2 = rcodcl2p[face_id];
      cs_real_t heq = _calc_heq(hint, rcodcl2) * surf;

      for (cs_lnum_t k = 0; k < 3; k++)
        _y[cell_id][k] += thetap * idiffp * heq * (pi[k] - pj[k]);
    }

  }

  /* Free memory */
  CS_FREE(x_j);
}

/*----------------------------------------------------------------------------
 * Add coupling term coordinates to matrix assembler.
 *
 * parameters:
 *   coupling_id
 *   r_g_id   <-- global row ids (per cell)
 *   ma       <-> matrix assembler
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_matrix_add_ids(int                     coupling_id,
                                    const cs_gnum_t        *r_g_id,
                                    cs_matrix_assembler_t  *ma)
{
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)cs_glob_mesh->b_face_cells;
  const cs_internal_coupling_t *cpl
    = cs_internal_coupling_by_id(coupling_id);

  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t n_local = cpl->n_local;

  const cs_lnum_t block_size = 800;
  cs_gnum_t g_row_id[800];
  cs_gnum_t g_col_id[800];

  cs_gnum_t *g_id_l, *g_id_d;
  CS_MALLOC(g_id_l, cs::max(n_local, n_distant), cs_gnum_t);
  CS_MALLOC(g_id_d, n_local, cs_gnum_t);

  /* local to global preparation and exchange */

  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    cs_lnum_t face_id = cpl->faces_distant[ii];
    cs_lnum_t cell_id = b_face_cells[face_id]; /* boundary to cell */
    g_id_l[ii] = r_g_id[cell_id];
  }

  ple_locator_exchange_point_var(cpl->locator,
                                 g_id_l,
                                 g_id_d,
                                 nullptr,
                                 sizeof(cs_gnum_t),
                                 1,
                                 0);

  /* local side */

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = cpl->faces_local[ii];
    cs_lnum_t cell_id = b_face_cells[face_id]; /* boundary to cell */
    g_id_l[ii] = r_g_id[cell_id];
  }

  cs_lnum_t jj = 0;
  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    g_row_id[jj] = g_id_l[ii];
    g_col_id[jj] = g_id_d[ii];
    jj++;
    if (jj >= block_size - 1) {
      cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
      jj = 0;
    }
  }
  if (jj > 0)
    cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);

  CS_FREE(g_id_l);
  CS_FREE(g_id_d);
}

/*----------------------------------------------------------------------------
 * Add coupling terms to matrix values assembly.
 *
 * parameters:
 *   f        <-- associated field
 *   db_size  <-- diagonal block size
 *   eb_size  <-- extra-diagonal block size
 *   r_g_id   <-- global row ids (per cell)
 *   mav      <-> matrix values assembler
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_matrix_add_values(const cs_field_t              *f,
                                       cs_lnum_t                      db_size,
                                       cs_lnum_t                      eb_size,
                                       const cs_gnum_t                r_g_id[],
                                       cs_matrix_assembler_values_t  *mav)
{
  const cs_real_t* b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)cs_glob_mesh->b_face_cells;

  int coupling_id = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));
  const cs_internal_coupling_t *cpl
    = cs_internal_coupling_by_id(coupling_id);

  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t n_local = cpl->n_local;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
  cs_real_t thetap = 0.0;
  int idiffp = 0;

  if (eqp->icoupl > 0) {
    thetap = eqp->theta;
    idiffp = eqp->idiff;
  }

  /* Compute global ids and exchange coefficient */

  cs_real_t *hintp = f->bc_coeffs->hint;
  cs_real_t *rcodcl2p = f->bc_coeffs->rcodcl2;

  /* local to global preparation and exchange */

  cs_gnum_t *g_id_l, *g_id_d;
  CS_MALLOC(g_id_l, cs::max(n_local, n_distant), cs_gnum_t);
  CS_MALLOC(g_id_d, n_local, cs_gnum_t);

  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    cs_lnum_t face_id = cpl->faces_distant[ii];
    cs_lnum_t cell_id = b_face_cells[face_id]; /* boundary to cell */
    g_id_l[ii] = r_g_id[cell_id];
  }

  ple_locator_exchange_point_var(cpl->locator,
                                 g_id_l,
                                 g_id_d,
                                 nullptr,
                                 sizeof(cs_gnum_t),
                                 1,
                                 0);

  /* local side */

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = cpl->faces_local[ii];
    cs_lnum_t cell_id = b_face_cells[face_id]; /* boundary to cell */
    g_id_l[ii] = r_g_id[cell_id];
  }

  /* Compute contributions to matrix */

  const cs_lnum_t block_size = 514;
  cs_gnum_t d_g_row_id[514];
  cs_real_t d_aij[514];
  cs_gnum_t e_g_row_id[514];
  cs_gnum_t e_g_col_id[514];
  cs_real_t e_aij[514];

  cs_lnum_t db_stride = db_size*db_size;
  cs_lnum_t eb_stride = eb_size*eb_size;

  assert(block_size > db_stride && block_size > eb_size);

  cs_lnum_t jj = 0, kk = 0, db_fill = 0, eb_fill = 0;
  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = cpl->faces_local[ii];

    cs_real_t surf = b_face_surf[face_id];
    cs_real_t hint = hintp[face_id];
    cs_real_t rcodcl2 = rcodcl2p[face_id];
    cs_real_t heq = _calc_heq(hint, rcodcl2) * surf;
    cs_real_t c = thetap * idiffp * heq;

    d_g_row_id[jj] = g_id_l[ii];
    e_g_row_id[kk] = g_id_l[ii]; e_g_col_id[kk] = g_id_d[ii];

    for (cs_lnum_t ib = 0; ib < db_stride; ib++)
      d_aij[db_fill + ib] = 0;
    for (cs_lnum_t ib = 0; ib < db_size; ib++)
      d_aij[db_fill + ib*(db_size + 1)] = c;

    for (cs_lnum_t ib = 0; ib < eb_stride; ib++)
      e_aij[eb_fill + ib] = 0;
    for (cs_lnum_t ib = 0; ib < eb_size; ib++)
      e_aij[eb_fill + ib*(eb_size + 1)] = -c;

    jj += 1;
    kk += 1;
    db_fill += db_stride;
    eb_fill += eb_stride;

    if (db_fill >= block_size - 1) {
      cs_matrix_assembler_values_add_g(mav, jj, d_g_row_id, d_g_row_id, d_aij);
      jj = 0;
      db_fill = 0;
    }

    if (eb_fill >= block_size - 1) {
      cs_matrix_assembler_values_add_g(mav, kk, e_g_row_id, e_g_col_id, e_aij);
      kk = 0;
      eb_fill = 0;
    }
  }

  cs_matrix_assembler_values_add_g(mav, jj, d_g_row_id, d_g_row_id, d_aij);
  cs_matrix_assembler_values_add_g(mav, kk, e_g_row_id, e_g_col_id, e_aij);

  /* Free memory */

  CS_FREE(g_id_l);
  CS_FREE(g_id_d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup internal coupling related parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_setup(void)
{
  if (_n_internal_couplings < 1)
    return;

  /* Call deprecated functions first */
  cs_user_internal_coupling_add_volumes(cs_glob_mesh);
  cs_user_internal_coupling_from_disjoint_meshes(cs_glob_mesh);

  /* Now do setup proper */

  int field_id;
  cs_field_t *f;
  const int coupling_key_id = cs_field_key_id("coupling_entity");
  int coupling_id = 0;

  const int n_fields = cs_field_n_fields();

  /* Definition of coupling_ids as keys of variable fields */

  coupling_id = 0;
  for (field_id = 0; field_id < n_fields; field_id++) {
    f = cs_field_by_id(field_id);
    if (f->type & CS_FIELD_VARIABLE) {
      const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
      if (eqp->icoupl > 0) {
        cs_field_set_key_int(f, coupling_key_id, coupling_id);
        // coupling_id++;
      }
    }
  }

  /* Initialization of coupling entities */

  coupling_id = 0;
  for (field_id = 0; field_id < n_fields; field_id++) {
    f = cs_field_by_id(field_id);
    if (f->type & CS_FIELD_VARIABLE) {
      const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
      if (eqp->icoupl > 0) {
        coupling_id++;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize internal coupling related structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize(void)
{
  for (int i = 0; i < _n_internal_couplings; i++) {
    cs_internal_coupling_t *cpl = _internal_coupling + i;
    _locator_initialize(cs_glob_mesh, cpl);
    _initialize_coupled_faces(cpl);
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
  if (cpl == nullptr)
    return;

  cs_gnum_t n_local = cpl->n_local;

  cs_parall_counter(&n_local, 1);

  if (cpl->cells_criteria != nullptr)
    bft_printf("   Cell group selection criterion: %s\n",
               cpl->cells_criteria);

  if (cpl->faces_criteria != nullptr)
    bft_printf("   Face group selection criterion: %s\n",
               cpl->faces_criteria);

  if (cpl->interior_faces_group_name != nullptr)
    bft_printf("   Assign interior faces group name: %s\n",
               cpl->interior_faces_group_name);

  if (cpl->exterior_faces_group_name != nullptr)
    bft_printf("   Assign interior faces group name: %s\n",
               cpl->exterior_faces_group_name);

  if (cpl->n_volume_zones > 0) {
    bft_printf("   Volume zones:\n");
    for (int i = 0; i < cpl->n_volume_zones; i++) {
      const cs_zone_t *z = cs_volume_zone_by_id(cpl->volume_zone_ids[i]);
      bft_printf("      %s\n", z->name);
    }
  }

  bft_printf("\n"
             "   Locator: n dist points (total coupled boundary faces) = %llu\n",
             (unsigned long long)n_local);
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
  cs_internal_coupling_t* cpl;

  if (_n_internal_couplings == 0)
    return;

  bft_printf("\n Internal coupling\n");
  for (int cpl_id = 0; cpl_id < _n_internal_couplings; cpl_id++) {
    cpl = _internal_coupling + cpl_id;
    bft_printf("   coupling_id = %d\n", cpl_id);
    cs_internal_coupling_log(cpl);
  }
}

/*----------------------------------------------------------------------------
 * Add preprocessing operations required by coupling volume using given
 * criteria.
 *
 * The volume is separated from the rest of the domain with inserted
 * boundaries.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure to modify
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_preprocess(cs_mesh_t  *mesh)
{
  for (int i = 0; i < _n_internal_couplings; i++) {
    cs_internal_coupling_t *cpl = _internal_coupling + i;
    if (   (cpl->cells_criteria != nullptr || cpl->n_volume_zones > 0)
        && cpl->faces_criteria == nullptr) {
      /* Insert wall boundaries,
       * locators are initialized afterwards */
      _volume_initialize_insert_boundary(mesh, cpl);
    }
  }
}

/*----------------------------------------------------------------------------
 * Define face to face mappings for internal couplings.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure to modify
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_map(cs_mesh_t   *mesh)
{
  /* Initialization of locators  for all coupling entities */

  for (int cpl_id = 0; cpl_id < _n_internal_couplings; cpl_id++) {
    cs_internal_coupling_t  *cpl = _internal_coupling + cpl_id;
    if (cpl->faces_criteria == nullptr)
      _auto_group_name(cpl, cpl_id);
    _volume_face_initialize(mesh, cpl);
  }
}

/*----------------------------------------------------------------------------
 * Define coupling entity using given criteria.
 *
 * parameters:
 *   f_id       <-- id of the field
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_entity(int        f_id)
{
  cs_field_t* f = cs_field_by_id(f_id);

  if (f->type & CS_FIELD_VARIABLE) {
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    eqp->icoupl = 1;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "field id = %d provided is invalid."
              " The field must be a variable.",
              f_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update internal coupling coefficients of the field of the
 * given id using given boundary exchange coefficients passed by face id.
 *
 * \param[in] f     pointer to field
 * \param[in] hbnd  boundary exchange coefficients passed by face id
 */
/*----------------------------------------------------------------------------*/

void
cs_ic_field_set_exchcoeff(const cs_field_t *f,
                          const cs_real_t  *hbnd)
{
  const int coupling_key_id = cs_field_key_id("coupling_entity");
  int coupling_id = cs_field_get_key_int(f,
                                         coupling_key_id);
  const cs_internal_coupling_t  *cpl
    = cs_internal_coupling_by_id(coupling_id);

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;
  cs_real_t *hint = f->bc_coeffs->hint;
  cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;

  cs_real_t *hextloc = nullptr;
  CS_MALLOC(hextloc, n_local, cs_real_t);

  /* Exchange hbnd */
  cs_internal_coupling_exchange_by_face_id(cpl,
                                           1, /* Dimension */
                                           hbnd,
                                           hextloc);

  /* Compute hint and hext */
  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = faces_local[ii];
    hint[face_id] = hbnd[face_id];
    rcodcl2[face_id] = hextloc[ii];
  }

  CS_FREE(hextloc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get distant data using face id at all coupling faces for a given
 * field id.
 *
 * \param[in]  field_id    field id
 * \param[in]  stride      number of values (interlaced) by entity
 * \param[in]  tab_distant exchanged data by face id
 * \param[out] tab_local   local data by face id
 */
/*----------------------------------------------------------------------------*/

void
cs_ic_field_dist_data_by_face_id(const int         field_id,
                                 int               stride,
                                 const cs_real_t   tab_distant[],
                                 cs_real_t         tab_local[])
{
  const cs_field_t* f = cs_field_by_id(field_id);

  const int coupling_key_id = cs_field_key_id("coupling_entity");
  int coupling_id = cs_field_get_key_int(f,
                                         coupling_key_id);
  const cs_internal_coupling_t  *cpl
    = cs_internal_coupling_by_id(coupling_id);
  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t *faces_local = cpl->faces_local;

  cs_real_t *local = nullptr;
  CS_MALLOC(local, n_local, cs_real_t);

  /* Exchange distant data by face id */
  cs_internal_coupling_exchange_by_face_id(cpl,
                                           stride,
                                           tab_distant,
                                           local);

  /* Save by face id */
  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = faces_local[ii];
    for (cs_lnum_t jj = 0; jj < stride; jj++)
      tab_local[stride * face_id + jj] = local[stride * ii + jj];
  }

  CS_FREE(local);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volumes as internal coupling zones.
 *
 * These zones will be separated from the rest of the domain using automatically
 * defined thin walls.
 *
 * \deprecated
 * move contents to\ref cs_user_internal_coupling instead.
 *
 * \param[in, out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_add_volumes(cs_mesh_t  *mesh)
{
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volumes from separated meshes as internal coupling zones.
 *
 * These zones must be disjoint and the face selection criteria must be specified.
 *
 * \deprecated
 * move contents to\ref cs_user_internal_coupling instead.
 *
 * \param[in, out]  mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_from_disjoint_meshes(cs_mesh_t  *mesh)
{
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*============================================================================
 * Public C++ function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update scalar boundary condition coefficients for internal coupling.
 *
 * \param[in]     ctx              reference to dispatch context
 * \param[in]     bc_coeffs        associated BC coefficients structure
 * \param[in]     cpl              structure associated with internal coupling
 * \param[in]     halo_type        halo type
 * \param[in]     w_stride         stride for weighting coefficient
 * \param[in]     clip_coeff       clipping coefficient
 * \param[in]     var              gradient's base variable
 * \param[in]     c_weight         weighted gradient coefficient variable,
 *                                 or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_update_bc_coeffs_s
(
 cs_dispatch_context           &ctx,
 const cs_field_bc_coeffs_t    *bc_coeffs,
 const cs_internal_coupling_t  *cpl,
 cs_halo_type_t                 halo_type,
 int                            w_stride,
 double                         clip_coeff,
 const cs_real_t               *var,
 const cs_real_t               *c_weight
)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* For internal coupling, exchange local variable
     with its associated distant value */

  cs_real_t *hintp = bc_coeffs->hint;
  cs_real_t *rcodcl2p = bc_coeffs->rcodcl2;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;
  const cs_lnum_t *faces_local = cpl->faces_local;

  cs_real_t *var_ext = nullptr, *var_distant = nullptr;
  CS_MALLOC(var_ext, n_local, cs_real_t);
  CS_MALLOC(var_distant, n_distant, cs_real_t);

  const cs_lnum_t *restrict b_face_cells = mesh->b_face_cells;
  cs_real_t *bc_coeff_a = bc_coeffs->a;
  cs_real_t *bc_coeff_b = bc_coeffs->b;

  /* For cases with a stronger gradient normal to the coupling than tangential
     to the coupling, assuming a homogeneous Neuman boundary condition at the
     coupled faces for the reconstruction at I' rather than the value at I on
     non-orthogonal meshes (such as tetrahedral meshes) can actually degrade
     performance, because the only adjacent mesh locations contributing
     information are not in the plane tangential to the face and containing II'.
     So we use an iterative process here to initialize BC coefficients with
     a non-reconstructed value and refine them with a reconstructed value.
     This is actually only necessary when combining a gradient tangential to the
     coupled surface and a non-orthogonal mesh at the wall (not recommended for
     wall law modeling), so we limit this to a single iteration and do not
     provide user setting for this now. */

  int n_iter_max = 2;
  for (int iter = 0; iter < n_iter_max; iter++) {

    if (iter > 0) {
      if (w_stride <= 1)
        cs_gradient_boundary_iprime_lsq_s(ctx,
                                          mesh,
                                          cs_glob_mesh_quantities,
                                          n_distant,
                                          faces_distant,
                                          halo_type,
                                          clip_coeff,
                                          bc_coeffs,
                                          c_weight,
                                          var,
                                          var_distant);
      else {
        assert(w_stride == 6);
        cs_gradient_boundary_iprime_lsq_s_ani(ctx,
                                              mesh,
                                              cs_glob_mesh_quantities,
                                              n_distant,
                                              faces_distant,
                                              clip_coeff,
                                              bc_coeffs,
                                              (const cs_real_6_t *)c_weight,
                                              var,
                                              var_distant);
      }
    }
    else {
      for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
        cs_lnum_t face_id = faces_distant[ii];
        cs_lnum_t cell_id = b_face_cells[face_id];
        var_distant[ii] = var[cell_id];
      }
    }

    cs_internal_coupling_exchange_var(cpl,
                                      1,
                                      (cs_real_t *)var_distant,
                                      (cs_real_t *)var_ext);

    /* For internal coupling, update BC coeffs */

    for (cs_lnum_t ii = 0; ii < n_local; ii++) {
      cs_lnum_t face_id = faces_local[ii];
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t hint = hintp[face_id];
      cs_real_t hext = rcodcl2p[face_id];

      cs_real_t m_a = hext / (hint + hext);
      cs_real_t m_b = hint / (hint + hext);

      bc_coeff_a[face_id] =   m_a * var_ext[ii]
                            + m_b * var[cell_id];
      bc_coeff_b[face_id] = 0;
    }
  }

  /* Last pass to ensure face values are identical on each side:
     use the mean of local and distant reconstructed values */

  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    cs_lnum_t face_id = faces_distant[ii];
    var_distant[ii] = bc_coeff_a[face_id];
  }

  cs_internal_coupling_exchange_var(cpl, 1, var_distant, var_ext);

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = faces_local[ii];
    bc_coeff_a[face_id] = 0.5*(bc_coeff_a[face_id] + var_ext[ii]);
  }

  CS_FREE(var_ext);
  CS_FREE(var_distant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update vector boundary condition coefficients for internal coupling.
 *
 * \param[in]     ctx              reference to dispatch context
 * \param[in]     bc_coeffs_v      boundary condition structure
 * \param[in]     cpl              structure associated with internal coupling
 * \param[in]     halo_type        halo type
 * \param[in]     clip_coeff       clipping coefficient
 * \param[in]     df_limiter       diffusion limiter array
 * \param[in]     var              gradient's base variable
 * \param[in]     c_weight         weighted gradient coefficient variable,
 *                                 or nullptr
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_internal_coupling_update_bc_coeffs_strided
(
 cs_dispatch_context           &ctx,
 const cs_field_bc_coeffs_t    *bc_coeffs_v,
 const cs_internal_coupling_t  *cpl,
 cs_halo_type_t                 halo_type,
 double                         clip_coeff,
 cs_real_t                     *df_limiter,
 const cs_real_t                var[][stride],
 const cs_real_t               *c_weight
)
{
  using var_t = cs_real_t[stride];
  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* For internal coupling, exchange local variable
     with its associated distant value */

  cs_real_t *hintp = bc_coeffs_v->hint;
  cs_real_t *rcodcl2p = bc_coeffs_v->rcodcl2;

  const cs_lnum_t n_local = cpl->n_local;
  const cs_lnum_t n_distant = cpl->n_distant;
  const cs_lnum_t *faces_distant = cpl->faces_distant;
  const cs_lnum_t *faces_local = cpl->faces_local;

  var_t *var_ext = nullptr, *var_ext_lim = nullptr;
  var_t *var_distant = nullptr, *var_distant_lim = nullptr;
  CS_MALLOC(var_ext, n_local, var_t);
  CS_MALLOC(var_distant, n_distant, var_t);

  if (df_limiter != nullptr) {
    CS_MALLOC(var_ext_lim, n_local, var_t);
    CS_MALLOC(var_distant_lim, n_distant, var_t);
  }

  const cs_lnum_t *restrict b_face_cells = mesh->b_face_cells;
  cs_real_3_t  *bc_coeff_a = (cs_real_3_t  *)bc_coeffs_v->a;
  cs_real_33_t *bc_coeff_b = (cs_real_33_t *)bc_coeffs_v->b;

  /* For cases with a stronger gradient normal to the coupling than tangential
     to the coupling, assuming a homogeneous Neuman boundary condition at the
     coupled faces for the reconstruction at I' rather than the value at I on
     non-orthogonal meshes (such as tetrahedral meshes) can actually degrade
     performance, because the only adjacent mesh locations contributing
     information are not in the plane tangential to the face and containing II'.
     So we use an iterative process here to initialize BC coefficients with
     a non-reconstructed value and refine them with a reconstructed value.
     This is actually only necessary whan combining a gradient tangential to the
     coupled surface and a non-orthogonal mesh at the wall (not recommended for
     wall law modeling), so we limit this to a single iteration and do not
     provide user setting for this now. */

  int n_iter_max = 2;
  for (int iter = 0; iter < n_iter_max; iter++) {

    if (iter > 0) {

      cs_gradient_boundary_iprime_lsq_strided<stride>(ctx,
                                                      mesh,
                                                      cs_glob_mesh_quantities,
                                                      n_distant,
                                                      faces_distant,
                                                      halo_type,
                                                      clip_coeff,
                                                      df_limiter,
                                                      bc_coeffs_v,
                                                      c_weight,
                                                      var,
                                                      var_distant,
                                                      var_distant_lim);
    }
    else {
      for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
        cs_lnum_t face_id = faces_distant[ii];
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t kk = 0; kk < stride; kk++) {
          var_distant[ii][kk] = var[cell_id][kk];
          if (df_limiter != nullptr)
            var_distant_lim[ii][kk] = var[cell_id][kk];
        }
      }
    }

    cs_internal_coupling_exchange_var(cpl,
                                      stride,
                                      (cs_real_t *)var_distant,
                                      (cs_real_t *)var_ext);

    if (df_limiter != nullptr)
      cs_internal_coupling_exchange_var(cpl,
                                        stride,
                                        (cs_real_t *)var_distant_lim,
                                        (cs_real_t *)var_ext_lim);

    /* For internal coupling, update BC coeffs */

    cs_real_3_t  *bc_coeff_af = (cs_real_3_t  *)bc_coeffs_v->af;
    cs_real_33_t *bc_coeff_bf = (cs_real_33_t *)bc_coeffs_v->bf;

    for (cs_lnum_t ii = 0; ii < n_local; ii++) {
      cs_lnum_t face_id = faces_local[ii];
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t hint = hintp[face_id];
      cs_real_t hext = rcodcl2p[face_id];

      cs_real_t m_a = hext / (hint + hext);
      cs_real_t m_b = hint / (hint + hext);

      cs_real_t heq = hext * m_b;

      for (cs_lnum_t kk = 0; kk < stride; kk++) {
        bc_coeff_a[face_id][kk] =   m_a * var_ext[ii][kk]
                                  + m_b * var[cell_id][kk];

        if (df_limiter != nullptr)
          bc_coeff_af[face_id][kk] = - heq * var_ext_lim[ii][kk];
        else
          bc_coeff_af[face_id][kk] = - heq * var_ext[ii][kk];

        for (cs_lnum_t ll = 0; ll < stride; ll++) {
          bc_coeff_b[face_id][kk][ll] = 0.;
          bc_coeff_bf[face_id][kk][ll] = 0.;
        }

        bc_coeff_bf[face_id][kk][kk] = heq;
      }
    }
  }

  /* Last pass to ensure face values are identical on each side:
     use the mean of local and distant reconstructed values */

  for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
    cs_lnum_t face_id = faces_distant[ii];
    for (cs_lnum_t kk = 0; kk < stride; kk++)
      var_distant[ii][kk] = bc_coeff_a[face_id][kk];
  }

  cs_internal_coupling_exchange_var(cpl,
                                    stride,
                                    (cs_real_t *)var_distant,
                                    (cs_real_t *)var_ext);

  for (cs_lnum_t ii = 0; ii < n_local; ii++) {
    cs_lnum_t face_id = faces_local[ii];
    for (cs_lnum_t kk = 0; kk < stride; kk++) {
      bc_coeff_a[face_id][kk] = 0.5*(  bc_coeff_a[face_id][kk]
                                     + var_ext[ii][kk]);
    }
  }

  CS_FREE(var_ext);
  CS_FREE(var_distant);
  CS_FREE(var_ext_lim);
  CS_FREE(var_distant_lim);
}

// Force instanciation

template void
cs_internal_coupling_update_bc_coeffs_strided
(cs_dispatch_context           &ctx,
 const cs_field_bc_coeffs_t    *bc_coeffs_v,
 const cs_internal_coupling_t  *cpl,
 cs_halo_type_t                 halo_type,
 double                         clip_coeff,
 cs_real_t                     *df_limiter,
 const cs_real_t                var[][3],
 const cs_real_t               *c_weight
);

template void
cs_internal_coupling_update_bc_coeffs_strided
(cs_dispatch_context           &ctx,
 const cs_field_bc_coeffs_t    *bc_coeffs_v,
 const cs_internal_coupling_t  *cpl,
 cs_halo_type_t                 halo_type,
 double                         clip_coeff,
 cs_real_t                     *df_limiter,
 const cs_real_t                var[][6],
 const cs_real_t               *c_weight
);

/*----------------------------------------------------------------------------*/
