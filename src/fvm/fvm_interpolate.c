/*============================================================================
 * Interpolate data defined on a nodal representation associated with a mesh.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in, out]  ldlt  pointer to matrix coefficients:
 *                        (m00, m10, m11, m20, m21, m22, m30, m31, m32, m33) in
 *                        (f00, l10, f11, l20, l21, f22, l30, l31, l32, f33) out
 */
/*----------------------------------------------------------------------------*/

static inline void
_sym_44_factor_ldlt(cs_real_t  ldlt[10])
{
  /* Factorization */

  // j=0
  const cs_real_t  d00 = ldlt[0]; // m00
  assert(fabs(d00) > cs_math_zero_threshold);

  const cs_real_t  f00 = 1. / d00;
  const cs_real_t  l10 = ldlt[1] * f00; // m01
  const cs_real_t  l20 = ldlt[3] * f00; // m02
  const cs_real_t  l30 = ldlt[6] * f00; // m03

  // j=1
  const cs_real_t  d11 = ldlt[2] - l10*l10*d00; // m11
  assert(fabs(d11) > cs_math_zero_threshold);
  const cs_real_t  f11 = 1. / d11;
  const cs_real_t  l21 = (ldlt[4] - l20*d00*l10) * f11; // m12
  const cs_real_t  l31 = (ldlt[7] - l30*d00*l10) * f11; // m13

  // j=2
  const cs_real_t  d22 = ldlt[5] - l20*d00*l20 - l21*d11*l21;  // m22
  assert(fabs(d22) > cs_math_zero_threshold);
  const cs_real_t  f22 = 1. / d22;
  const cs_real_t  l32 = (ldlt[8] - l30*d00*l20 - l31*d11*l21) * f22;  // m23

  // j=3
  const cs_real_t  d33 = ldlt[9] - l30*d00*l30 - l31*d11*l31 - l32*d22*l32; // m33
  assert(fabs(d33) > cs_math_zero_threshold);
  const cs_real_t  f33 = 1. / d33;

  ldlt[0] = f00;
  ldlt[1] = l10;
  ldlt[2] = f11;
  ldlt[3] = l20;
  ldlt[4] = l21;
  ldlt[5] = f22;
  ldlt[6] = l30;
  ldlt[7] = l31;
  ldlt[8] = l32;
  ldlt[9] = f33;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * Here we only need to use the last element of the solution vector,
 * so we return that value only and simplify the computation.
 *
 * \param[in, out]  ldlt  pointer to matrix coefficients:
 *                        (f00, l10, f11, l20, l21, f22, l30, l31, l32, f33)
 * \param[in]       rhs   right-hand side
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_sym_44_partial_solve_ldlt(const cs_real_t  ldlt[10],
                           const cs_real_t  rhs[4])
{
  /* f00, f11, f22, f33, l32, l31, l30, l21, l20, l10
     0    1    2    3    4    5    6    7    8    9   */

  cs_real_t x[4]; /* solution */

  x[0] = rhs[0];
  x[1] = rhs[1] - x[0]*ldlt[1];
  x[2] = rhs[2] - x[0]*ldlt[3] - x[1]*ldlt[4];
  x[3] = rhs[3] - x[0]*ldlt[6] - x[1]*ldlt[7] - x[2]*ldlt[8];

  x[3] = x[3]*ldlt[9];

  return x[3];

  /*
    x[2] = x[2]*ldlt[5] - ldlt[8]*x[3];
    x[1] = x[1]*ldlt[2] - ldlt[7]*x[3] - ldlt[4]*x[2];
    x[0] = x[0]*ldlt[0] - ldlt[6]*x[3] - ldlt[3]*x[2] - ldlt[1]*x[1];
  */
}

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
  if (this_nodal == NULL || n_points == 0)
    return;

  /* Sanity checks */
  assert(   point_coords != NULL  && location_id != NULL
         && src_data != NULL && dest_data != NULL);
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
            if (this_nodal->parent_vertex_id != NULL)
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

        _sym_44_factor_ldlt(a);

        for (cs_lnum_t k = 0; k < data_dim; k++) {
          dest_data[p_id*data_dim + k] = _sym_44_partial_solve_ldlt(a, r[k]);
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
          if (this_nodal->parent_vertex_id != NULL)
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

        _sym_44_factor_ldlt(a);

        for (cs_lnum_t k = 0; k < data_dim; k++) {
          dest_data[p_id*data_dim + k] = _sym_44_partial_solve_ldlt(a, r[k]);
        }

      } /* Section type */

    } /* Loop on sections */

  } /* Loop on points */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
