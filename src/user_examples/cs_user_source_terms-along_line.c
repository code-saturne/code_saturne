/*============================================================================
 * Additional user-defined source terms for variable equations.
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_source_terms-along_line.c
 *
 * \brief Base examples for additional right-hand side source terms for
 *   variable equations (momentum, scalars, turbulence...).
 *
 * See the reference \ref cs_user_source_terms.c for documentation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define source terms.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  static cs_lnum_t n_elts = -1;    /* >= 0 after initialization */
  static cs_lnum_t *elt_ids = NULL;
  static cs_real_t *seg_c_len = NULL;
  cs_real_3_t *seg_c_cen = NULL;

  /* x, y, z of origin, and x, y, z of destination */
  /* From North to South at the middle of the first elevation*/
  cs_real_3_t *point_coords;
  cs_lnum_t n_points = 2;
  BFT_MALLOC(point_coords, n_points, cs_real_3_t); /* Segment: 2 points */

  point_coords[0][0] = 0.; point_coords[0][1] = 15.; point_coords[0][2] = 0.5;
  point_coords[1][0] = 0.; point_coords[1][1] = 0.;  point_coords[1][2] = 0.5;

  /* Initialize on first pass */
  if (n_elts < 0) {

    /* Select cells and count length
     * Note: elt_ids and seg_c_len are allocated here, and
     * should be deallocated at the end of the calculation */
    cs_mesh_intersect_polyline_cell_select(point_coords,
                                           n_points,
                                           &n_elts,
                                           &elt_ids,
                                           &seg_c_len,
                                           &seg_c_cen);

    cs_real_t len = 0.;
    /* To visualize selected cells */
    cs_field_t *length = cs_field_by_name_try("line_length_per_cell");
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      len += seg_c_len[i];
      cs_lnum_t cell_id = elt_ids[i];
      if (length != NULL)
        length->val[cell_id] = seg_c_len[i];
    }

    cs_gnum_t n_g_elts = n_elts;

    cs_parall_sum(1, CS_REAL_TYPE, &len);
    cs_parall_sum(1, CS_GNUM_TYPE, &n_g_elts);

    bft_printf("BRIN [%f, %f, %f] -> [%f, %f, %f] n_cells=%ld, length=%f\n",
               point_coords[0][0], point_coords[0][1], point_coords[0][2],
               point_coords[1][0], point_coords[1][1], point_coords[1][2],
               n_g_elts, len);
  }

  /* Source term, proportional to the length of the segment in each cell
   * for scalar1 */
  cs_field_t *fld = cs_field_by_name_try("scalar1");
  if (fld->id == f_id) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t cell_id = elt_ids[i];
      st_exp[cell_id] = seg_c_len[i];
    }
  }

  /* Free memory */
  BFT_FREE(seg_c_cen);
  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max) {
    BFT_FREE(elt_ids);
    BFT_FREE(seg_c_len);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
