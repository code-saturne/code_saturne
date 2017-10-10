/*============================================================================
 * Mesh deformation.
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

/* Displacement to prescribe */

static cs_lnum_t     _vd_size = 0;
static cs_real_3_t  *_vd = NULL;
static cs_lnum_t     _cs_comp_id[] = {0, 1, 2};

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
 * \param[in]      input      pointer to a structure cast on-the-fly (may be NULL)
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
      res[id] = - _vd[id][c_id];
    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t  id = pt_ids[p];
      res[p] = - _vd[id][c_id];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for the computational domain:
 *         -- which type of boundaries closed the computational domain
 *         -- the settings for the time step
 *         -- activate predefined equations or modules
 *         -- add user-defined properties and/or advection fields
 *         -- add user-defined equations
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_cdo_init_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* Choose a boundary by default.
     Valid choice is CS_PARAM_BOUNDARY_WALL or CS_PARAM_BOUNDARY_SYMMETRY */

  cs_domain_set_default_boundary(domain, CS_PARAM_BOUNDARY_SYMMETRY);

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};

  for (int i = 0; i < 3; i++) {

    cs_equation_add(eq_name[i],        // equation name
                    eq_name[i],        // associated variable field name
                    CS_EQUATION_TYPE_PREDEFINED,
                    1,                 // dimension of the unknown
                    CS_PARAM_BC_HMG_NEUMANN); // default boundary

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify the elements such as properties, advection fields,
 *         user-defined equations and modules which have been previously
 *         added.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 * \param[in]        n_zones    number of zones with prescribed deformation
 * \param[in]        zones_ids  ids of zones with prescribed deformation
 */
/*----------------------------------------------------------------------------*/

static void
_cdo_finalize_setup(cs_domain_t   *domain,
                    int            n_zones,
                    const int      zone_ids[])
{
  CS_UNUSED(domain);

  cs_property_t  *conductivity = cs_property_by_name("unity");

  /* Retrieve the equation to set
     >> cs_equation_t  *eq = cs_equation_by_name("eq_name");  */

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};

  for (int i = 0; i < 3; i++) {

    cs_equation_t  *eq = cs_equation_by_name(eq_name[i]);

    for (int j = 0; j < n_zones; j++) {
      const cs_boundary_zone_t *z = cs_boundary_zone_by_id(zone_ids[j]);
      cs_equation_add_bc_by_analytic(eq,
                                     CS_PARAM_BC_DIRICHLET,
                                     z->name,
                                     _define_displ_bcs,
                                     _cs_comp_id + i);
    }

    cs_equation_link(eq, "diffusion", conductivity);

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 * \param[in]        n_zones    number of zones with prescribed deformation
 * \param[in]        zones_ids  ids of zones with prescribed deformation
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_initialize(cs_domain_t   *domain,
                          int            n_zones,
                          const int      zone_ids[])
{
  if (cs_glob_domain == NULL)
    cs_glob_domain = cs_domain_create();

  _cdo_init_setup(cs_glob_domain);

  _cdo_finalize_setup(domain,
                      n_zones,
                      zone_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prescribe the displacement vector for a set of vertices.
 *
 * Only those values at vertices matching boundary zones with prescribed
 * displacement will really be used.
 *
 * This function may be called multiple times before displacements are
 * are computed, so displacements may be prescribed for different zones
 * separately or simulataneously.
 *
 * \param[in]  n_vertices    number of vertices at which to prescribe
 *                           displacements
 * \param[in]  vertex_ids    ids of vertices at which to prescribe
 *                           displacements, or NULL if [0, ... n_vertices-1]
 * \param[in]  displacement  pointer to prescribed displacements
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_prescribe_displacement(cs_lnum_t          n_vertices,
                                      const cs_lnum_t    vertex_ids[],
                                      const cs_real_3_t  displacement[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  if (_vd_size != m->n_vertices) {
    BFT_REALLOC(_vd, m->n_vertices, cs_real_3_t);
#   pragma omp parallel for if (_vd_size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _vd[i][j] = 0.;
    }
    _vd_size = m->n_vertices;
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
 * \brief Compute displacement for mesh deformation.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_solve_displacement(void)
{
  /*  Build high-level structures and create algebraic systems */
  cs_domain_t  *domain = cs_glob_domain;

  cs_domain_initialize_systems(domain);

  /* Define the current time step */
  // cs_domain_define_current_time_step(domain);

  /* Build and solve equations related to the deformation */

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};

  for (int i = 0; i < 3; i++) {

    cs_equation_t *eq = cs_equation_by_name(eq_name[i]);

    /* Define the algebraic system */
    cs_equation_build_system(domain->mesh,
                             domain->time_step,
                             domain->dt_cur,
                             eq);

    /* Solve the algebraic system */
    cs_equation_solve(eq);

  }

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
  BFT_FREE(_vd);
  _vd_size = 0;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
