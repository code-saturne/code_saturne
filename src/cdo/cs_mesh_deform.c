/*============================================================================
 * Mesh deformation.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_boundary_zone.h"
#include "cs_domain.h"
#include "cs_domain_setup.h"
#include "cs_equation.h"
#include "cs_mesh_builder.h"
#include "cs_mesh_group.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_deform.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool  _active = false;

/* Displacement to prescribe */

static int           _n_b_zones = 0;
static int          *_b_zone_ids = NULL;

static cs_lnum_t     _vd_size = 0;
static cs_real_3_t  *_vd = NULL;
static cs_lnum_t     _cs_comp_id[] = {0, 1, 2};

static bool          _fixed_vertices = false;
static cs_lnum_t     _n_fixed_vertices = 0;
static cs_lnum_t    *_fixed_vtx_ids = NULL;
static cs_real_3_t  *_fixed_vtx_values = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the dirichlet boundary conditions
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time       when ?
 * \param[in]      n_pts      number of points to consider
 * \param[in]      pt_ids     list of points ids (to access coords and fill)
 * \param[in]      coords     where ?
 * \param[in]      compact    true:no indirection, false:indirection for filling
 * \param[in]      input      pointer to a structure cast on-the-fly
 *                            (may be NULL)
 * \param[in, out] res        result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_define_displ_bcs(cs_real_t           time,
                  cs_lnum_t           n_pts,
                  const cs_lnum_t    *pt_ids,
                  const cs_real_t    *xyz,
                  bool                compact,
                  void               *input,
                  cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(n_pts);
  CS_UNUSED(xyz);
  CS_UNUSED(compact);

  const cs_lnum_t c_id = *((const cs_lnum_t *)input);

  assert(pt_ids != NULL);

  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t  id = pt_ids[p];
      res[id] = _vd[id][c_id];
    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t  id = pt_ids[p];
      res[p] = _vd[id][c_id];
    }

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if mesh deformation is activated
 *
 * \return true if mesh deformation computation is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_mesh_deform_is_activated(void)
{
  if (_active)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the future mesh deformation
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_activate(void)
{
  if (_active)
    return;

  _active = true;

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};

  for (int i = 0; i < 3; i++) {

    cs_equation_t  *eq =
      cs_equation_add(eq_name[i], // equation name
                      eq_name[i], // associated variable field name
                      CS_EQUATION_TYPE_PREDEFINED,
                      1,                 // dimension of the unknown
                      CS_PARAM_BC_HMG_NEUMANN); // default boundary

    cs_equation_param_t  *eqp = cs_equation_get_param(eq);

    /* System to solve is SPD by construction */
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "amg");

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the boundary zones on which mesh deformation is prescribed.
 *
 * Only those values at vertices matching boundary zones with prescribed
 * displacement will really be used.
 *
 * \param[in]  n_boundary_zones   number of boundary zones at which to
 *                                prescribe displacements
 * \param[in]  boundary_zone_ids  ids of boundary zones at which to
 *                                prescribe displacements
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_define_dirichlet_bc_zones(cs_lnum_t  n_boundary_zones,
                                         const int  boundary_zone_ids[])
{
  if (n_boundary_zones != _n_b_zones) {
    _n_b_zones = n_boundary_zones;
    BFT_REALLOC(_b_zone_ids, _n_b_zones, int);
    memcpy(_b_zone_ids, boundary_zone_ids, sizeof(int)*_n_b_zones);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup the equations related to mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_setup(cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  cs_property_t  *conductivity = cs_property_by_name("unity");

  /* Retrieve the equation parameter to set
     >> cs_equation_param_t  *eqp = cs_equation_param_by_name("eq_name");  */

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};

  for (int i = 0; i < 3; i++) {

    cs_equation_param_t  *eqp = cs_equation_param_by_name(eq_name[i]);

    for (int j = 0; j < _n_b_zones; j++) {
      const cs_zone_t *z = cs_boundary_zone_by_id(_b_zone_ids[j]);
      cs_equation_add_bc_by_analytic(eqp,
                                     CS_PARAM_BC_DIRICHLET,
                                     z->name,
                                     _define_displ_bcs,
                                     _cs_comp_id + i);
    }

    if (_fixed_vertices) {
      cs_real_t  *fixed_vtx_values;
      BFT_MALLOC(fixed_vtx_values, _n_fixed_vertices, cs_real_t);
      if (_fixed_vtx_values != NULL) {
#       pragma omp parallel for if (_n_fixed_vertices > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < _n_fixed_vertices; j++)
          fixed_vtx_values[j] = _fixed_vtx_values[j][i];
      }
      else {
#       pragma omp parallel for if (_n_fixed_vertices > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < _n_fixed_vertices; j++)
          fixed_vtx_values[j] = 0;
      }

      cs_equation_enforce_vertex_dofs(eqp,
                                      _n_fixed_vertices,
                                      _fixed_vtx_ids,
                                      NULL, /* reference value */
                                      fixed_vtx_values);

      BFT_FREE(fixed_vtx_values);
    }

    cs_equation_add_diffusion(eqp, conductivity);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prescribe the displacement vector for a set of vertices.
 *
 * Only those values at vertices matching boundary zones with prescribed
 * displacement will really be used, as defined by
 * \ref cs_mesh_deform_define_dirichlet_bc_zones.
 *
 * When calling this function multiple times for different vertex sets,
 * the most recently defined values are used for vertices belonging to
 * multiple sets.
 *
 * \param[in]  n_vertices         number of vertices at which to prescribe
 *                                displacements
 * \param[in]  vertex_ids         ids of vertices at which to prescribe
 *                                displacements, or NULL if
 *                                [0, ... n_vertices-1]
 * \param[in]  displacement       pointer to prescribed displacements
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_prescribe_displacement(cs_lnum_t          n_vertices,
                                      const cs_lnum_t    vertex_ids[],
                                      const cs_real_3_t  displacement[])
{
  const cs_mesh_t *m = cs_glob_mesh;

  if (_vd_size != m->n_vertices) {
    _vd_size = m->n_vertices;
    BFT_REALLOC(_vd, _vd_size, cs_real_3_t);
#   pragma omp parallel for if (_vd_size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < _vd_size; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _vd[i][j] = 0.;
    }
  }

  if (vertex_ids != NULL) {
#   pragma omp parallel for if (_vd_size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++) {
      cs_lnum_t v_id = vertex_ids[i];
      for (cs_lnum_t j = 0; j < 3; j++)
        _vd[v_id][j] = displacement[i][j];
    }
  }

  else { /* if (vertex_ids == NULL) */
#   pragma omp parallel for if (_vd_size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _vd[i][j] = displacement[i][j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a fixed displacement vector for given vertices.
 *
 * This displacement is enforced at all given vertices, including interior
 * vertices.
 *
 * If this function is called multiple times, the previous definitions
 * are overwritten, so all displacements of this type must be defined
 * in a single call to this function.
 *
 * \param[in]  n_vertices         number of vertices at which to prescribe
 *                                displacements
 * \param[in]  vertex_ids         ids of vertices at which to prescribe
 *                                displacements, or NULL if
 *                                [0, ... n_vertices-1]
 * \param[in]  displacement       pointer to prescribed displacements,
 *                                or NULL for no displacement
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_force_displacements(cs_lnum_t          n_vertices,
                                   const cs_lnum_t    vertex_ids[],
                                   const cs_real_3_t  displacement[])
{
  BFT_REALLOC(_fixed_vtx_ids, n_vertices, cs_lnum_t);

  if (displacement != NULL)
    BFT_REALLOC(_fixed_vtx_values, n_vertices, cs_real_3_t);
  else
    BFT_FREE(_fixed_vtx_values);

  _n_fixed_vertices = n_vertices;

  cs_gnum_t n_g_v = n_vertices;
  cs_parall_counter(&n_g_v, 1);
  _fixed_vertices = (n_g_v > 0) ? true : false;

  if (vertex_ids != NULL) {
#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++)
      _fixed_vtx_ids[i] = vertex_ids[i];
  }
  else {
#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++)
      _fixed_vtx_ids[i] = i;
  }

  if (displacement != NULL) {
#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _fixed_vtx_values[i][j] = displacement[i][j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute displacement for mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_solve_displacement(cs_domain_t  *domain)
{
  /* Build and solve equations related to the deformation */

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};

  for (int i = 0; i < 3; i++) {

    cs_equation_t *eq = cs_equation_by_name(eq_name[i]);

    if (cs_equation_uses_new_mechanism(eq))
      cs_equation_solve_steady_state(domain->mesh, eq);

    else { /* Deprecated */

      /* Sanity check */
      assert(cs_equation_is_steady(eq));

      /* Define the algebraic system */
      cs_equation_build_system(domain->mesh, eq);

      /* Solve the algebraic system */
      cs_equation_solve_deprecated(eq);

    }

  } /* Loop on Cartesian component */

  {
    cs_field_t *fx = cs_field_by_name("mesh_deform_x");
    cs_field_t *fy = cs_field_by_name("mesh_deform_y");
    cs_field_t *fz = cs_field_by_name("mesh_deform_z");

    cs_mesh_t *m = cs_glob_mesh;

    assert(_vd_size == m->n_vertices);

#   pragma omp parallel for if (_vd_size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      _vd[i][0] = fx->val[i];
      _vd[i][1] = fy->val[i];
      _vd[i][2] = fz->val[i];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to current mesh displacement vector.
 *
 * \return  pointer to current displacement vector
 */
/*----------------------------------------------------------------------------*/

const cs_real_3_t *
cs_mesh_deform_get_displacement(void)
{
  return (const cs_real_3_t *)_vd;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free structures used fo mesh deformation.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_finalize(void)
{
  BFT_FREE(_b_zone_ids);
  _n_b_zones = 0;

  BFT_FREE(_vd);
  _vd_size = 0;

  BFT_FREE(_fixed_vtx_ids);
  BFT_FREE(_fixed_vtx_values);
  _n_fixed_vertices = 0;
  _fixed_vertices = false;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
