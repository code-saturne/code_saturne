/*============================================================================
 * Interpolate data defined on a nodal representation associated with a mesh.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"

#include "cs_assert.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_interpolate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Interpolate vertex-based values to points located relative to a mesh.
 *
 * Least-squared based interpolation is used for now.
 *
 * parameters:
 *   this_nodal         <-- pointer to nodal mesh representation structure
 *   entity_dim         <-- 3 for location on cells, 2 for faces, 1 for edges,
 *                          and 0 for vertices
 *   data_dim           <-- data dimension
 *   n_points           <-- number of points to locate
 *   location_id        <-- id of element (with concatenated sections)
 *                          in which each point is located
 *   point_coords       <-- point coordinates
 *   src_data           <-- source data (interleaved)
 *   dest_data          <-> destination data (interleaved)
 *----------------------------------------------------------------------------*/

void
fvm_interpolate_vtx_data(const fvm_nodal_t       *this_nodal,
                         int                      entity_dim,
                         int                      data_dim,
                         cs_lnum_t                n_points,
                         const cs_lnum_t          location_id[],
                         const cs_coord_t         point_coords[],
                         const cs_real_t          src_data[],
                         cs_real_t                dest_data[])
{
  if (this_nodal == nullptr || n_points == 0)
    return;

  /* Sanity checks */
  assert(   point_coords != nullptr  && location_id != nullptr
         && src_data != nullptr && dest_data != nullptr);
  assert(this_nodal->dim == 3);

  if (this_nodal->dim != 3)
    return;

  const cs_coord_t  *vtx_coords = this_nodal->vertex_coords;

  /* Loop on points */

# pragma omp parallel for if (n_points > CS_THR_MIN)
  for (cs_lnum_t p_id = 0; p_id < n_points; p_id++) {

    cs_lnum_t elt_id = location_id[p_id];

    if (elt_id < 0)  /* Unlocated point */
      continue;

    for (int sec_id = 0; sec_id < this_nodal->n_sections; sec_id++) {

      const fvm_nodal_section_t  *section = this_nodal->sections[sec_id];

      if (section->entity_dim != entity_dim)
        continue;

      if (section->n_elements < elt_id) {
        elt_id -= section->n_elements;
        continue;
      }

      /* Location is now given by elt_id relative to current section */

      if (section->type == FVM_CELL_POLY) {

        cs_assert(data_dim <= 9);

        cs_real_t  a[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        cs_real_t  r[9][4];

        for (int i = 0; i < data_dim; i++) {
          r[i][0] = 0;
          r[i][1] = 0;
          r[i][2] = 0;
          r[i][3] = 0;
        }

        const cs_coord_t  *p_coo = point_coords + 3*p_id;

        for (cs_lnum_t j = section->face_index[elt_id];
             j < section->face_index[elt_id + 1]; j++) {

          cs_lnum_t  f_id = CS_ABS(section->face_num[j]) - 1;
          for (cs_lnum_t k = section->vertex_index[f_id];
               k < section->vertex_index[f_id + 1]; k++) {

            cs_lnum_t  v_id = section->vertex_num[k] - 1;
            if (this_nodal->parent_vertex_id != nullptr)
              v_id = this_nodal->parent_vertex_id[v_id];
            const cs_coord_t  *v_coo = vtx_coords + 3*v_id;

            cs_real_t r_coo[3]
              = {v_coo[0]-p_coo[0], v_coo[1]-p_coo[1], v_coo[2]-p_coo[2]};

            a[0] += r_coo[0] * r_coo[0]; // a00
            a[1] += r_coo[1] * r_coo[0]; // a10
            a[2] += r_coo[1] * r_coo[1]; // a11
            a[3] += r_coo[2] * r_coo[0]; // a20
            a[4] += r_coo[2] * r_coo[1]; // a21
            a[5] += r_coo[2] * r_coo[2]; // a22
            a[6] += r_coo[0];            // a30
            a[7] += r_coo[1];            // a31
            a[8] += r_coo[2];            // a32
            a[9] += 1;                   // a33

            const cs_real_t *s_var = src_data + data_dim*v_id;
            for (cs_lnum_t l = 0; l < data_dim; l++) {
              r[l][0] += r_coo[0] * s_var[l];
              r[l][1] += r_coo[1] * s_var[l];
              r[l][2] += r_coo[2] * s_var[l];
              r[l][3] += s_var[l];
            }

          }

        }

        cs_math_sym_44_factor_ldlt(a);

        for (cs_lnum_t k = 0; k < data_dim; k++) {
          dest_data[p_id*data_dim + k] = cs_math_sym_44_partial_solve_ldlt(a, r[k]);
        }

      }

      else { /* Polygonal or regular elements */

        cs_lnum_t v_s_id, v_e_id;

        if (section->type == FVM_FACE_POLY) {
          v_s_id = section->vertex_index[elt_id];
          v_e_id = section->vertex_index[elt_id+1];
        }
        else {
          const cs_lnum_t  stride = section->stride;
          v_s_id = elt_id * stride;
          v_e_id = (elt_id+1) * stride;
        }

        cs_assert(data_dim <= 9);

        cs_real_t  a[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        cs_real_t  r[9][4];

        for (int i = 0; i < data_dim; i++) {
          r[i][0] = 0;
          r[i][1] = 0;
          r[i][2] = 0;
          r[i][3] = 0;
        }

        const cs_coord_t  *p_coo = point_coords + 3*p_id;

        for (cs_lnum_t j = v_s_id; j < v_e_id; j++) {

          cs_lnum_t  v_id = section->vertex_num[j] - 1;
          if (this_nodal->parent_vertex_id != nullptr)
            v_id = this_nodal->parent_vertex_id[v_id];

          const cs_coord_t  *v_coo = vtx_coords + 3*v_id;

          cs_real_t r_coo[3]
            = {v_coo[0]-p_coo[0], v_coo[1]-p_coo[1], v_coo[2]-p_coo[2]};

          a[0] += r_coo[0] * r_coo[0]; // a00
          a[1] += r_coo[1] * r_coo[0]; // a10
          a[2] += r_coo[1] * r_coo[1]; // a11
          a[3] += r_coo[2] * r_coo[0]; // a20
          a[4] += r_coo[2] * r_coo[1]; // a21
          a[5] += r_coo[2] * r_coo[2]; // a22
          a[6] += r_coo[0];            // a30
          a[7] += r_coo[1];            // a31
          a[8] += r_coo[2];            // a32
          a[9] += 1;                   // a33

          const double *s_var = src_data + data_dim*v_id;
          for (cs_lnum_t k = 0; k < data_dim; k++) {
            r[k][0] += r_coo[0] * s_var[k];
            r[k][1] += r_coo[1] * s_var[k];
            r[k][2] += r_coo[2] * s_var[k];
            r[k][3] += s_var[k];
          }

        }

        cs_math_sym_44_factor_ldlt(a);

        for (cs_lnum_t k = 0; k < data_dim; k++) {
          dest_data[p_id*data_dim + k] = cs_math_sym_44_partial_solve_ldlt(a, r[k]);
        }

      } /* Section type */

    } /* Loop on sections */

  } /* Loop on points */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
