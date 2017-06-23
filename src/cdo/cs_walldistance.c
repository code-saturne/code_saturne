/*============================================================================
 * Compute the wall distance using the CDO framework
 *============================================================================*/

/* VERS */

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

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_post.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_cdo.h"
#include "cs_param.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_walldistance.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

/*============================================================================
 * Private variables
 *============================================================================*/

static cs_equation_t  *cs_walldistance_eq = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the wall distance for a face-based scheme
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      eq        pointer to the associated cs_equation_t structure
 * \param[in]      field     pointer to a cs_field_t structure
 * \param[in, out] dist      array storing the wall distance to compute
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cdofb(const cs_cdo_connect_t     *connect,
               const cs_cdo_quantities_t  *cdoq,
               const cs_equation_t        *eq,
               const cs_field_t           *field,
               cs_real_t                   dist[])
{
  cs_lnum_t  i, k;

  const cs_real_t  *c_var = field->val;
  const cs_real_t  *f_var = cs_equation_get_face_values(eq);

  /* Loop on cells */
  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    cs_real_3_t  cell_gradient = {0., 0., 0.};
    cs_real_t  inv_cell_vol = 1/cdoq->cell_vol[c_id];

    for (i = connect->c2f->idx[c_id]; i < connect->c2f->idx[c_id+1]; i++) {

      cs_lnum_t  f_id = connect->c2f->col_id[i];
      cs_quant_t  fq = cdoq->face[f_id];
      int  sgn = connect->c2f->sgn[i];
      cs_real_t  dualedge_contrib = fq.meas*sgn*(f_var[f_id] - c_var[c_id]);

      for (k = 0; k < 3; k++)
        cell_gradient[k] += dualedge_contrib*fq.unitv[k];

    } // Loop on cell faces

    for (k = 0; k < 3; k++)
      cell_gradient[k] *= inv_cell_vol;

    /* Compute the distance from the wall at this cell center */
    cs_real_t  tmp = _dp3(cell_gradient, cell_gradient) + 2*c_var[c_id];
    assert(tmp >= 0.); // Sanity check

    dist[c_id] = sqrt(tmp) - cs_math_3_norm(cell_gradient);

  } // Loop on cells

  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    field->name,
                    1,               // dim
                    false,           // interlace
                    true,            // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    dist,            // values on cells
                    NULL,            // values at internal faces
                    NULL,            // values at border faces
                    NULL);           // time step management structure

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the wall distance for a vertex-based scheme
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      field     pointer to a cs_field_t structure
 * \param[in, out] dist      array storing the wall distance to compute
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cdovb(const cs_cdo_connect_t     *connect,
               const cs_cdo_quantities_t  *cdoq,
               const cs_field_t           *field,
               cs_real_t                   dist[])
{
  /* Initialize arrays */
  cs_real_t  *dualcell_vol = NULL;
  cs_real_3_t  *vtx_gradient = NULL;

  BFT_MALLOC(vtx_gradient, cdoq->n_vertices, cs_real_3_t);
  BFT_MALLOC(dualcell_vol, cdoq->n_vertices, cs_real_t);

# pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
    vtx_gradient[i][0] = vtx_gradient[i][1] = vtx_gradient[i][2] = 0.;
    dualcell_vol[i] = 0.;
  }

  /* Reconstruct gradient at vertices from gradient at cells */
  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_real_t  *var = field->val;

  for (cs_lnum_t  c_id = 0; c_id < cdoq->n_cells; c_id++) {

    cs_real_3_t  grd_cell;

    cs_reco_grd_cell_from_pv(c_id, connect, cdoq, var, grd_cell);

    for (cs_lnum_t i = c2v->idx[c_id]; i < c2v->idx[c_id+1]; i++) {

      cs_lnum_t  v_id = c2v->ids[i];

      dualcell_vol[v_id] += cdoq->dcell_vol[i];
      for (int k = 0; k < 3; k++)
        vtx_gradient[v_id][k] += cdoq->dcell_vol[i]*grd_cell[k];

    } // Loop on cell vertices

  } // Loop on cells

# pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
    cs_real_t  inv_dualcell_vol = 1/dualcell_vol[i];
    for (int k = 0; k < 3; k++)
      vtx_gradient[i][k] *= inv_dualcell_vol;
  }

  /* Compute now wall distance at each vertex */
# pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
    cs_real_t  tmp = _dp3(vtx_gradient[i], vtx_gradient[i]) + 2*var[i];
    assert(tmp >= 0); // Sanity check
    dist[i] = sqrt(tmp) - cs_math_3_norm(vtx_gradient[i]);
  }

  /* Post-processing */
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_ALL_ASSOCIATED,
                           field->name,
                           1,               // dim
                           false,           // interlace
                           true,            // true = original mesh
                           CS_POST_TYPE_cs_real_t,
                           dist,            // values on vertices
                           NULL);           // time step management structure

  /* Free memory */
  BFT_FREE(dualcell_vol);
  BFT_FREE(vtx_gradient);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if the computation of the wall distance is activated
 *
 * \return true if the wall distance computation is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_walldistance_is_activated(void)
{
  if (cs_walldistance_eq != NULL)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the future computation of the wall distance
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_activate(void)
{
  /* Sanity check */
  assert(cs_walldistance_eq == NULL);

  cs_walldistance_eq =
    cs_equation_add("WallDistance",              // equation name
                    "WallDistance",              // variable name
                    CS_EQUATION_TYPE_PREDEFINED, // type of the equation
                    1,                           // dimension of the variable
                    CS_PARAM_BC_HMG_NEUMANN);    // default BC
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the equation related to the wall distance
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_setup(void)
{
  cs_equation_t  *eq = cs_walldistance_eq;

  if (cs_walldistance_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Stop setting the wall distance equation.\n"
              " The wall distance computation has not been activated.");

  /* Unity is a property defined by default */
  cs_equation_link(eq, "diffusion", cs_property_by_name("unity"));

  /* Add boundary conditions */
  cs_real_t  zero_value = 0.;
  const char *bc_zone_name =
    cs_param_get_boundary_domain_name(CS_PARAM_BOUNDARY_WALL);

  cs_equation_add_bc_by_value(eq,
                              CS_PARAM_BC_HMG_DIRICHLET,
                              bc_zone_name,
                              &zero_value);

  /* Add source term */
  const char *st_zone_name =
    cs_mesh_location_get_name(CS_MESH_LOCATION_CELLS);
  cs_real_t  unity = 1.0;

  cs_equation_add_source_term_by_val(eq,
                                     st_zone_name,   // zone name
                                     &unity);        // value to set

  /* Enforcement of the Dirichlet boundary conditions */
  cs_equation_set_param(eq, CS_EQKEY_BC_ENFORCEMENT, "penalization");

  /* System to solve is SPD by construction */
  cs_equation_set_param(eq, CS_EQKEY_ITSOL, "cg");

#if defined(HAVE_PETSC)  /* Modify the default settings */
  cs_equation_set_param(eq, CS_EQKEY_SOLVER_FAMILY, "petsc");
  cs_equation_set_param(eq, CS_EQKEY_PRECOND, "amg");
#else
  cs_equation_set_param(eq, CS_EQKEY_PRECOND, "jacobi");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the wall distance
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_compute(const cs_mesh_t              *mesh,
                        const cs_time_step_t         *time_step,
                        double                        dt_cur,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *cdoq)
{
  /* First step:
     Solve the equation related to the definition of the wall distance. */

  cs_equation_t  *eq = cs_walldistance_eq;

  /* Sanity check */
  assert(cs_equation_is_steady(eq));

  /* Define the algebraic system */
  cs_equation_build_system(mesh, time_step, dt_cur, eq);

  /* Solve the algebraic system */
  cs_equation_solve(eq);

  /* Second step:
     Compute the wall distance. */

  cs_field_t  *field = cs_equation_get_field(eq);

  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(field->location_id);

  /* Sanity checks */
  assert(field->is_owner);
  assert(field->dim == 1);

  /* Initialize dist array */
  cs_real_t  *dist = NULL;
  BFT_MALLOC(dist, n_elts[0], cs_real_t);
# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts[0]; i++)
    dist[i] = 0;

  cs_space_scheme_t  space_scheme = cs_equation_get_space_scheme(eq);
  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    assert(n_elts[0] == cdoq->n_vertices);
    _compute_cdovb(connect, cdoq, field, dist);
    break;

  case CS_SPACE_SCHEME_CDOFB:
    assert(n_elts[0] == cdoq->n_cells);
    _compute_cdofb(connect, cdoq, eq, field, dist);
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    bft_error(__FILE__, __LINE__, 0,
              " CDO Vertex+Cell-based is not yet implemented to compute"
              " the wall distance.");
    break;

  default:
    assert(0);
    break;
  }

  /* Replace field values by dist */
# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts[0]; i++)
    field->val[i] = dist[i];

  /* Free memory */
  BFT_FREE(dist);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
