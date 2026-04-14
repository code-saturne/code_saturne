/*============================================================================
 * Repartitioning and redistribution mesh data and fields.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <ctime>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_assert.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_ext_neighborhood.h"
#include "base/cs_field.h"
#include "base/cs_mem.h"
#include "base/cs_renumber.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_coherency.h"

#include "alge/cs_cell_to_vertex.h"
#include "alge/cs_gradient.h"
#include "alge/cs_matrix_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_renumber_update.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Apply reordering to integer arrays.
 *
 * parameters:
 *   n_b_faces    <-- number of boundary faces
 *   stride       <-- number of sub values per sent element
 *   buffer       --- work array
 *   data         <-> values to write
 *   b_face_n2o   <-- boundary faces new-to-old mapping
 *----------------------------------------------------------------------------*/

template <typename T>
static void
_reorder_bc_vals(cs_lnum_t         n_b_faces,
                 int               stride,
                 T                *buffer,
                 T                *data,
                 const cs_lnum_t   b_face_n2o[])
{
  if (b_face_n2o == nullptr || data == nullptr)
    return;

  memcpy(buffer, data, n_b_faces*stride*sizeof(T));

  for (cs_lnum_t new_id = 0; new_id < n_b_faces; new_id++) {
    cs_lnum_t old_id = b_face_n2o[new_id];
    const T *orig = buffer + stride*old_id;
    T *dest = data + stride*new_id;
    for (int i = 0; i < stride; i++)
      dest[i] = orig[i];
  }
}

/*----------------------------------------------------------------------------
 * Distribute BC types.
 *
 * parameters:
 *   bfd           <-- pointer to parallel distributor
 *   n_b_faces     <-- number of boundary faces
 *   b_face_n2o    <-- boundary faces new-to-old mapping
 *----------------------------------------------------------------------------*/

static void
_reorder_bc_type(const cs_lnum_t   b_face_n2o[])
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_boundary_condition_pm_info_t *pm = cs_glob_bc_pm_info;

  int *buffer;
  CS_MALLOC_HD(buffer, mesh->n_b_faces, int, cs_alloc_mode);

  _reorder_bc_vals(mesh->n_b_faces, 1, buffer, pm->izfppp, b_face_n2o);
  _reorder_bc_vals(mesh->n_b_faces, 1, buffer, pm->iautom, b_face_n2o);

  int **bc_type = cs_boundary_conditions_get_bc_type_addr();
  _reorder_bc_vals(mesh->n_b_faces, 1, buffer, *bc_type, b_face_n2o);

  CS_FREE(buffer);
}

/*----------------------------------------------------------------------------
 * Redistribute fields and bc_coeffs based on the mesh distributors.
 *
 * parameters:
 *   cd             <-- pointer to cells distributor
 *   n_cells_ini    <-- number of cells to send
 *   cell_order     <-- cells ordering by associated global number
 *   cell_n2o       <-- cells new-to-old mapping
 *   bfd            <-- pointer to boundary faces distributor
 *   n_b_faces_ini  <-- number of boundary faces to send
 *   b_face_order   <-- boundary face ordering by associated global number
 *   b_face_n2o     <-- boundary faces new-to-old mapping
 *   ifd            <-- pointer to internal faces distributor
 *   n_i_faces_ini  <-- number of internal faces to send
 *   i_face_lst     <-- list of internal faces to send
 *   i_face_order   <-- internal faces ordering by associated global number
 *   i_face_n2o     <-- internal faces new-to-old mapping
 *----------------------------------------------------------------------------*/

static void
_renumber_update_fields(cs_mesh_t        *mesh,
                        const cs_lnum_t   cell_n2o[],
                        const cs_lnum_t   i_face_n2o[],
                        const cs_lnum_t   b_face_n2o[],
                        const cs_lnum_t   vtx_n2o[])
{
  // Cell-centered fields.

  int n_fields = cs_field_n_fields();

  cs_lnum_t buffer_size = 0;
  cs_real_t *buffer = nullptr;

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *field = cs_field_by_id(i);

    const cs_lnum_t *n2o = nullptr;
    cs_lnum_t n_elts = 0;

    cs_mesh_location_type_t mlt = cs_mesh_location_get_type(field->location_id);
    switch(mlt) {
    case CS_MESH_LOCATION_NONE:
      break;
    case CS_MESH_LOCATION_CELLS:
      n2o = cell_n2o;
      n_elts = mesh->n_cells;
      break;
    case CS_MESH_LOCATION_INTERIOR_FACES:
      n2o = i_face_n2o;
      n_elts = mesh->n_i_faces;
      break;
    case CS_MESH_LOCATION_BOUNDARY_FACES:
      n2o = i_face_n2o;
      n_elts = mesh->n_b_faces;
      break;
    case CS_MESH_LOCATION_VERTICES:
      n2o = vtx_n2o;
      n_elts = mesh->n_vertices;
      break;
    default:
      cs_assert(0); // Unhandled case
    }

    /* If the field's mesh location is a subset of a "main" location,
       values must be projected to that base location first.
       As long as this is not handled, throw an error... */
    cs_assert(mlt == (cs_mesh_location_type_t)field->location_id);

    if (n2o != nullptr) {
      cs_lnum_t n_vals = n_elts * (cs_lnum_t)field->dim;
      if (buffer_size < n_vals) {
        CS_FREE(buffer);
        buffer_size = n_vals;
        CS_MALLOC_HD(buffer, buffer_size, cs_real_t, cs_alloc_mode);
      }

      for (int j = 0; j < field->n_time_vals; j++) {
        cs_real_t *vals = field->vals[j];
        memcpy(buffer, vals, n_vals*sizeof(cs_real_t));

        cs_lnum_t stride = field->dim;
        for (cs_lnum_t k = 0; k < n_elts; k++) {
          cs_lnum_t src_id = n2o[k];
          for (cs_lnum_t ic = 0; ic < stride; ic++)
            vals[k*stride + ic] = buffer[src_id*stride + ic];
        }
      }
    }

    /* Handle bc coefficients */

    cs_field_bc_coeffs_t *bcc = field->bc_coeffs;
    if (   bcc != nullptr
        && mlt == CS_MESH_LOCATION_CELLS
        && b_face_n2o != nullptr) {

      cs_lnum_t n_b_faces = mesh->n_b_faces;

      int i_mult, a_mult, b_mult;
      cs_field_get_bc_coeff_mult(field, &i_mult, &a_mult, &b_mult);

      assert(a_mult <= b_mult);
      cs_lnum_t n_vals = n_b_faces * b_mult;
      if (buffer_size < n_vals) {
        CS_FREE(buffer);
        buffer_size = n_vals;
        CS_MALLOC_HD(buffer, buffer_size, cs_real_t, cs_alloc_mode);
      }

      _reorder_bc_vals(n_b_faces, i_mult,
                       (int *)buffer, bcc->icodcl, b_face_n2o);

      // Note: the rcodcl family of buffers is non-interleaved.

      for (int dim = 0; dim < field->dim; dim++) {
        cs_real_t *rc1 = bcc->rcodcl1 + (cs_lnum_t)dim*mesh->n_b_faces;
        _reorder_bc_vals(n_b_faces, 1, buffer, rc1, b_face_n2o);
      }

      for (int dim = 0; dim < field->dim; dim++) {
        cs_real_t *rc2 = bcc->rcodcl2 + (cs_lnum_t)dim*mesh->n_b_faces;
        _reorder_bc_vals(n_b_faces, 1, buffer, rc2, b_face_n2o);
      }

      for (int dim = 0; dim < field->dim; dim++) {
        cs_real_t *rc3 = bcc->rcodcl3 + (cs_lnum_t)dim*mesh->n_b_faces;
        _reorder_bc_vals(n_b_faces, 1, buffer, rc3, b_face_n2o);
      }

      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->a, b_face_n2o);
      _reorder_bc_vals(n_b_faces, b_mult, buffer, bcc->b, b_face_n2o);

      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->af, b_face_n2o);
      _reorder_bc_vals(n_b_faces, b_mult, buffer, bcc->bf, b_face_n2o);

      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->ad, b_face_n2o);
      _reorder_bc_vals(n_b_faces, b_mult, buffer, bcc->bd, b_face_n2o);
      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->ac, b_face_n2o);
      _reorder_bc_vals(n_b_faces, b_mult, buffer, bcc->bc, b_face_n2o);

      _reorder_bc_vals(n_b_faces, 1, buffer, bcc->h_int_tot, b_face_n2o);

      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->val_f, b_face_n2o);
      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->val_f_pre, b_face_n2o);

      _reorder_bc_vals(n_b_faces, a_mult, buffer, bcc->flux_diff, b_face_n2o);

    }
  }

  CS_FREE(buffer);
}

/*----------------------------------------------------------------------------
 * Free field gradients.
 *
 * Gradients are assumed to be computed on-the-fly when needed.
 *----------------------------------------------------------------------------*/

static void
_free_field_gradients(void)
{
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->grad)
      CS_FREE(f->grad);
  }
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber mesh and update all associated data.
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_update(void)
{
  cs_mesh_t *mesh = cs_glob_mesh;

  /* Free some elements which will need to be rebuilt when/if used.
     Halos are needed at some interior steps so destroyed later. */

  cs_mesh_quantities_free_all(cs_glob_mesh_quantities);

  _free_field_gradients();
  cs_gradient_free_quantities();
  cs_cell_to_vertex_free();

  cs_lnum_t *cell_n2o = nullptr, *vtx_n2o= nullptr;
  cs_lnum_t *i_face_n2o = nullptr, *b_face_n2o = nullptr;

  cs_renumber_mesh(mesh, &cell_n2o, &i_face_n2o, &b_face_n2o, &vtx_n2o);

  /* Re-compute mesh related quantities */

  cs_mesh_update_auxiliary(mesh);
  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);

  cs_mesh_init_group_classes(mesh);

  cs_mesh_update_selectors(mesh);
  cs_mesh_location_build(mesh, -1);
  cs_volume_zone_build_all(true);
  cs_boundary_zone_build_all(true);

#if defined(DEBUG)
  cs_mesh_coherency_check();
#endif

  cs_mesh_adjacencies_update_mesh();
  cs_matrix_update_mesh();

  /* Re-distribute fields */

  _renumber_update_fields(mesh,
                          cell_n2o,
                          i_face_n2o,
                          b_face_n2o,
                          vtx_n2o);

  _reorder_bc_type(b_face_n2o);

  CS_FREE(vtx_n2o);
  CS_FREE(b_face_n2o);
  CS_FREE(i_face_n2o);
  CS_FREE(cell_n2o);
}

/*----------------------------------------------------------------------------*/
