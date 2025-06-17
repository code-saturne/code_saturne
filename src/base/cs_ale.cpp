/*============================================================================
 * Functions associated to ALE formulation
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_dispatch.h"
#include "alge/cs_cell_to_vertex.h"
#include "cdo/cs_cdo_quantities.h"
#include "cdo/cs_cdo_connect.h"
#include "cdo/cs_cdo_main.h"
#include "base/cs_dispatch.h"
#include "alge/cs_divergence.h"
#include "cdo/cs_domain.h"
#include "cdo/cs_domain_setup.h"
#include "cdo/cs_equation.h"
#include "base/cs_equation_iterative_solve.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_operator.h"
#include "gui/cs_gui_mobile_mesh.h"
#include "base/cs_interface.h"
#include "base/cs_log.h"
#include "base/cs_physical_constants.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_bad_cells.h"
#include "base/cs_parall.h"
#include "base/cs_post.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"
#include "base/cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_ale.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

static cs_ale_data_t  _cs_glob_ale_data = {.impale = nullptr,
                                           .bc_type = nullptr,
                                           .ale_iteration = 0,
                                           .i_mass_flux_ale = nullptr,
                                           .b_mass_flux_ale = nullptr};

cs_ale_type_t cs_glob_ale = CS_ALE_NONE;

cs_ale_data_t  *cs_glob_ale_data = &_cs_glob_ale_data;

/*! Number of iterations for fluid flow initialization */
int cs_glob_ale_n_ini_f = 0;

/*! Indicate whether an iteration to initialize ALE is required */
int cs_glob_ale_need_init = -999;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  /* Array of values to set the boundary conditions (array allocated to all
     mesh vertices since there is no dedicated numbering to identify easily
     vertices lying on the boundary */
  cs_real_t    *vtx_values;

  /* Pointer on a list of arrays of selected vertices for each zone associated
     to a definition by array \ref CS_XDEF_BY_ARRAY for
     CS_BOUNDARY_ALE_IMPOSED_VEL, CS_BOUNDARY_ALE_IMPOSED_DISP,
     CS_BOUNDARY_ALE_INTERNAL_COUPLING and CS_BOUNDARY_ALE_EXTERNAL_COUPLING
     (definitions by value or by a sliding condition are excluded. The position
     of the array in the list is given implicitly and is related to the order
     in which the boundary conditions are defined. (cf. \ref
     cs_ale_setup_boundaries definition for more details). The aim of this list
     of arrays is to speed-up the settings of the boundary conditions by
     avoiding doing several times the same enforcement. */

  int           n_selections;   /* Number of selections */
  cs_lnum_t    *n_vertices;     /* Number of vertices in each selections  */
  cs_lnum_t   **vtx_select;     /* List of vertices for each selection */

} cs_ale_cdo_bc_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_real_3_t  *_vtx_coord0 = nullptr;
static cs_ale_cdo_bc_t  *_cdo_bc = nullptr;

static bool cs_ale_active = false;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute center of gravity and surface of a set of nb_points
 *
 * \param[in]      n_nodes    nodes number
 * \param[in]      xn         nodes coordinates
 * \param[inout]   cog        face center of gravity
 * \param[inout]   normal     normal ported by his surface
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cog_surf_from_nodes(int           n_nodes,
                             cs_real_t     xn[][3],
                             cs_real_t     cog[3],
                             cs_real_t     normal[3])
{
  cs_real_t cog0[3];

  for (int k = 0; k < 3; k++) {
    /* Add first node as last node */
    xn[n_nodes][k] = xn[0][k];

    /* Init cog and normal */
    cog[k] = 0.;
    normal[k] = 0.;
    cog0[k] = 0.;
  }

  /* First estimation of center */
  for (int inod = 0; inod < n_nodes; inod++)
    for (int k = 0; k < 3; k++)
      cog0[k] += xn[inod][k];

  for (int k = 0; k < 3; k++)
    cog0[k] /= (double)(n_nodes);

  cs_real_t surftot = 0.;

  /* loop on edges */
  for (int edge = 0; edge < n_nodes; edge++) {

    cs_real_t x1[3], x2[3];

    for (int k = 0; k < 3; k++) {
      x1[k] = xn[edge][k] - cog0[k];
      x2[k] = xn[edge+1][k] - cog0[k];
    }

    cs_real_t pvec[3];
    pvec[0] = x1[1]*x2[2] - x1[2]*x2[1];
    pvec[1] = x1[2]*x2[0] - x1[0]*x2[2];
    pvec[2] = x1[0]*x2[1] - x1[1]*x2[0];

    for (int k = 0; k < 3; k++)
      normal[k] += 0.5 * pvec[k];

    cs_real_t surfn = (pvec[0]*pvec[0] + pvec[1]*pvec[1] + pvec[2]*pvec[2]);
    surftot += surfn;

    for (int k = 0; k < 3; k++)
      cog[k] += surfn * (x1[k] + x2[k]);
  }

  for (int k = 0; k < 3; k++)
    cog[k] = cog[k]/(3. * cs_math_fmax(surftot,1.e-20)) + cog0[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the values of mesh vertices lying on a free surface boundary
 *         condition.
 *
 * \param[in]  domain      pointer to a \ref cs_domain_t structure
 * \param[in]  z           pointer to a \ref cs_zone_t structure
 * \param[in]  select_id   id in the list of selected vertices
 */
/*----------------------------------------------------------------------------*/

static void
_free_surface(const cs_domain_t  *domain,
              const cs_zone_t    *z,
              const int           select_id)
{
  const cs_real_t *grav = cs_glob_physical_constants->gravity;
  const cs_mesh_t *m = domain->mesh;
  const cs_real_3_t *restrict vtx_coord = (const cs_real_3_t *)m->vtx_coord;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;
  const cs_real_t *restrict b_face_surf = mq->b_face_surf;
  const cs_nreal_3_t *restrict b_face_u_normal = mq->b_face_u_normal;
  const cs_real_3_t *restrict b_face_cog = mq->b_face_cog;

  /* Boundary mass flux */
  int iflmab = cs_field_get_key_int(CS_F_(vel),
                                    cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Transform face flux to vertex displacement */
  cs_real_3_t *_mesh_vel = nullptr;

  /* Dual surface associated to vertices */
  cs_real_t *_v_surf = nullptr;

  /* Squared sum of partial surface associated to a vertex */
  cs_real_t *_is_loc_min = nullptr;
  cs_real_t *_is_loc_max = nullptr;

  cs_real_3_t normal;
  /* Normal direction is given by the gravity */
  cs_math_3_normalize(grav, normal);

  const cs_real_t  invdt = 1./domain->time_step->dt_ref; /* JB: dt[0] ? */

  CS_MALLOC(_mesh_vel, m->n_vertices, cs_real_3_t);
  CS_MALLOC(_v_surf, m->n_vertices, cs_real_t);
  CS_MALLOC(_is_loc_min, m->n_vertices, cs_real_t);
  CS_MALLOC(_is_loc_max, m->n_vertices, cs_real_t);

  cs_array_real_fill_zero(m->n_vertices, _v_surf);
  cs_array_real_fill_zero(3 * m->n_vertices, (cs_real_t *)_mesh_vel);
  cs_array_real_set_scalar(m->n_vertices, 1., _is_loc_min);
  cs_array_real_set_scalar(m->n_vertices, 1., _is_loc_max);

  /* First Loop over boundary faces
   * to compute if there is a local min and max elevation */
  for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {

    const cs_lnum_t face_id = z->elt_ids[elt_id];

    const cs_lnum_t s = m->b_face_vtx_idx[face_id];
    const cs_lnum_t e = m->b_face_vtx_idx[face_id+1];

    /* Compute the portion of surface associated to v_id_1 */
    for (cs_lnum_t k = s; k < e; k++) {

      const cs_lnum_t v_id = m->b_face_vtx_lst[k];

      cs_real_3_t v_cog = {
        b_face_cog[face_id][0] - vtx_coord[v_id][0],
        b_face_cog[face_id][1] - vtx_coord[v_id][1],
        b_face_cog[face_id][2] - vtx_coord[v_id][2],
      };

      /* g . (x_f - x_N) S_fN  */
      cs_real_t dz_fn = cs_math_3_dot_product(normal, v_cog);

      /* v1 is higher than x_f */
      if (dz_fn > 0.)
        _is_loc_min[v_id] = 0.; /* not a min */

      /* x_f is higher than v1 */
      if (dz_fn < 0.)
        _is_loc_max[v_id] = 0.; /* not a max */

    }

  } /* Loop on selected border faces */

  /* Handle parallelism */
  if (m->vtx_interfaces != nullptr) {

    cs_interface_set_min(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         _is_loc_min);

    cs_interface_set_min(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         _is_loc_max);
  }

  cs_gnum_t _f_count_filter = 0;
  cs_gnum_t _f_n_elts = z->n_elts;

  /* Loop over boundary faces */
  for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {

    const cs_lnum_t face_id = z->elt_ids[elt_id];

    const cs_real_t mf_inv_rho_n_dot_s = b_mass_flux[face_id] /
      (  cs_math_3_dot_product(normal, b_face_u_normal[face_id])
       * b_face_surf[face_id]
       * CS_F_(rho_b)->val[face_id]);

    cs_real_3_t f_vel = {
      normal[0] * mf_inv_rho_n_dot_s,
      normal[1] * mf_inv_rho_n_dot_s,
      normal[2] * mf_inv_rho_n_dot_s
    };

    const cs_lnum_t s = m->b_face_vtx_idx[face_id];
    const cs_lnum_t e = m->b_face_vtx_idx[face_id+1];

    cs_real_t f_has_min = 0;
    cs_real_t f_has_max = 0;

    /* First loop to determine if filtering is needed:
     * - if face has a local min and a local max
     * - if the slope is greater than pi/4
     *   */
    int f_need_filter = 0;

    for (cs_lnum_t k = s; k < e; k++) {
      const cs_lnum_t k_1 = (k+1 < e) ? k+1 : k+1-(e-s);
      const cs_lnum_t v_id0 = m->b_face_vtx_lst[k];
      const cs_lnum_t v_id1 = m->b_face_vtx_lst[k_1];
      /* Edge to CoG vector */
      cs_real_3_t e_cog = {
        0.5 * (b_face_cog[face_id][0] - vtx_coord[v_id0][0])
          + 0.5 * (b_face_cog[face_id][0] - vtx_coord[v_id1][0]),
        0.5 * (b_face_cog[face_id][1] - vtx_coord[v_id0][1])
          + 0.5 * (b_face_cog[face_id][1] - vtx_coord[v_id1][1]),
        0.5 * (b_face_cog[face_id][2] - vtx_coord[v_id0][2])
          + 0.5 * (b_face_cog[face_id][2] - vtx_coord[v_id1][2])
      };
      cs_real_t dz = cs::abs(cs_math_3_dot_product(normal, e_cog));
      cs_real_3_t e_cog_hor;
      cs_math_3_orthogonal_projection(normal, e_cog, e_cog_hor);
      cs_real_t dx = cs_math_3_norm(e_cog_hor);
      /* Too high slope */
      if (dz > dx)
        f_need_filter = 1;

      f_has_max = cs::max(f_has_max, _is_loc_max[v_id0]);
      f_has_min = cs::max(f_has_min, _is_loc_min[v_id0]);
    }

    if ((f_has_max > 0.5 && f_has_min > 0.5) || f_need_filter == 1) {
      f_need_filter = 1;
      _f_count_filter++;
    }

    /* Compute the portion of surface associated to v_id_1 */
    for (cs_lnum_t k = s; k < e; k++) {

      const cs_lnum_t k_1 = (k+1 < e) ? k+1 : k+1-(e-s);
      const cs_lnum_t k_2 = (k+2 < e) ? k+2 : k+2-(e-s);
      const cs_lnum_t v_id0 = m->b_face_vtx_lst[k];
      const cs_lnum_t v_id1 = m->b_face_vtx_lst[k_1];
      const cs_lnum_t v_id2 = m->b_face_vtx_lst[k_2];

      cs_real_3_t v0v1 = {
        vtx_coord[v_id1][0] - vtx_coord[v_id0][0],
        vtx_coord[v_id1][1] - vtx_coord[v_id0][1],
        vtx_coord[v_id1][2] - vtx_coord[v_id0][2]
      };
      cs_real_3_t v1v2 = {
        vtx_coord[v_id2][0] - vtx_coord[v_id1][0],
        vtx_coord[v_id2][1] - vtx_coord[v_id1][1],
        vtx_coord[v_id2][2] - vtx_coord[v_id1][2]
      };
      cs_real_3_t v1_cog = {
        b_face_cog[face_id][0] - vtx_coord[v_id1][0],
        b_face_cog[face_id][1] - vtx_coord[v_id1][1],
        b_face_cog[face_id][2] - vtx_coord[v_id1][2]
      };

      /* Portion of the surface associated to the vertex
       * projected in the normal direction */
      cs_real_t portion_surf
        = -0.25 * (   cs_math_3_triple_product(v0v1, v1_cog, normal)
                   + cs_math_3_triple_product(v1v2, v1_cog, normal));

      _v_surf[v_id1] += portion_surf;

      /* g . (x_f - x_N) S_fN  */
      cs_real_t dz_fn = cs_math_3_dot_product(normal, v1_cog);

      cs_real_t coeff = invdt * f_need_filter * dz_fn;
      for (int i = 0; i < 3; i++)
        _mesh_vel[v_id1][i] += (f_vel[i] + coeff * normal[i]) * portion_surf;
    }

  } /* Loop on selected border faces */

  const cs_equation_param_t *eqp =
    cs_field_get_equation_param_const(CS_F_(mesh_u));

  if (eqp->verbosity >= 1) {
    cs_gnum_t _sum[2] = { _f_count_filter, _f_n_elts };
    cs_parall_sum(2, CS_GNUM_TYPE, _sum);
    _f_count_filter = _sum[0];
    _f_n_elts       = _sum[1];
    bft_printf("Free surface condition %d: %f percents of limited face\n",
               select_id,
               (cs_real_t)_f_count_filter / (cs_real_t)_f_n_elts);
  }


  /* Handle parallelism */
  if (m->vtx_interfaces != nullptr) {
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         _mesh_vel);

    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         _v_surf);
  }

  /* Loop on selected border vertices */
  for (cs_lnum_t i = 0; i < _cdo_bc->n_vertices[select_id]; i++) {

    const cs_lnum_t v_id = _cdo_bc->vtx_select[select_id][i];
    const double  invsurf = 1./_v_surf[v_id];

    cs_real_t  *_val = _cdo_bc->vtx_values + 3*v_id;

    for (int k = 0; k < 3; k++) {
      _mesh_vel[v_id][k] *= invsurf;
      _val[k] = _mesh_vel[v_id][k];
    }

  } /* Loop on selected vertices */

  /* Free memory */
  CS_FREE(_mesh_vel);
  CS_FREE(_v_surf);
  CS_FREE(_is_loc_min);
  CS_FREE(_is_loc_max);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new selection of vertices related to a definition by array
 *        Update the cs_ale_cdo_bc_t structure.
 *
 * \param[in]      mesh     pointer to a \ref cs_mesh_t structure
 * \param[in]      z        pointer to a cs_zone_t structure
 * \param[in, out] vtag     tag array on vertices
 */
/*----------------------------------------------------------------------------*/

static void
_update_bc_list(const cs_mesh_t   *mesh,
                const cs_zone_t   *z,
                bool               vtag[])
{
  const cs_lnum_t  *bf2v_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t  *bf2v_lst = mesh->b_face_vtx_lst;
  const cs_lnum_t  n_vertices = mesh->n_vertices;
  const int  id = _cdo_bc->n_selections;

  cs_lnum_t  counter = 0;

  _cdo_bc->n_selections++;
  CS_REALLOC(_cdo_bc->n_vertices, _cdo_bc->n_selections, cs_lnum_t);
  CS_REALLOC(_cdo_bc->vtx_select, _cdo_bc->n_selections, cs_lnum_t *);

  /* Reset vtag */
  cs_array_bool_fill_false(n_vertices, vtag);

  /* Count the number of vertices to select */
  for (cs_lnum_t i = 0; i < z->n_elts; i++) {

    const cs_lnum_t bf_id = z->elt_ids[i];
    const cs_lnum_t *idx = bf2v_idx + bf_id;
    const cs_lnum_t *lst = bf2v_lst + idx[0];

    /* Loop on face vertices */
    for (cs_lnum_t j = 0; j < idx[1]-idx[0]; j++) {
      cs_lnum_t  v_id = lst[j];
      if (!vtag[v_id]) {  /* Not already selected */
        vtag[v_id] = true;
        counter++;
      }
    }

  } /* Loop on selected border faces */

  _cdo_bc->n_vertices[id] = counter;
  CS_MALLOC(_cdo_bc->vtx_select[id], counter, cs_lnum_t);

  /* Fill the list of selected vertices */
  cs_array_bool_fill_false(n_vertices, vtag);
  counter = 0;
  for (cs_lnum_t i = 0; i < z->n_elts; i++) {

    const cs_lnum_t bf_id = z->elt_ids[i];
    const cs_lnum_t *idx = bf2v_idx + bf_id;
    const cs_lnum_t *lst = bf2v_lst + idx[0];

    /* Loop on face vertices */
    for (cs_lnum_t j = 0; j < idx[1]-idx[0]; j++) {
      cs_lnum_t v_id = lst[j];
      if (!vtag[v_id]) {  /* Not already selected */
        vtag[v_id] = true;
        _cdo_bc->vtx_select[id][counter] = v_id;
        counter++;
      }
    }

  } /* Loop on selected border faces */

  assert(counter == _cdo_bc->n_vertices[id]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function updates all ALE BCs for CDO except free surface.
 *        These BCs are required for the fluid BCs and therefore needs to be
 *        updated before.
 *
 * \param[in]   domain       domain quantities
 * \param[out]  ale_bc_type  ALE boundary condition type
 * \param[out]  b_fluid_vel  boundary fluid velocity
 */
/*----------------------------------------------------------------------------*/

static void
_update_bcs(const cs_domain_t  *domain,
            int                 ale_bc_type[],
            cs_real_3_t         b_fluid_vel[])
{
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  const cs_real_3_t *restrict vtx_coord  = (const cs_real_3_t *)m->vtx_coord;
  const cs_real_3_t *restrict b_face_cog = mq->b_face_cog;

  cs_field_t *f_displ = cs_field_by_name("mesh_displacement");

  /* Only a part of the boundaries has to be updated */
  int  select_id = 0;
  for (int b_id = 0; b_id < domain->ale_boundaries->n_boundaries; b_id++) {

    const int z_id = domain->ale_boundaries->zone_ids[b_id];
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);

    switch(domain->ale_boundaries->types[b_id]) {

    case CS_BOUNDARY_ALE_IMPOSED_VEL:
      {
        /* Retrieve the velocities to enforce  on faces */
        cs_real_t *bc_vals = cs_gui_mobile_mesh_get_fixed_velocity(z->name);

        assert(select_id < _cdo_bc->n_selections);

#if 0 //TODO when we will have meg on vertices
        /* Loop on selected border vertices */
        for (cs_lnum_t i = 0; i < _cdo_bc->n_vertices[select_id]; i++) {

          cs_real_t  *_val
            = _cdo_bc->vtx_values + 3*_cdo_bc->vtx_select[select_id][i];
          _val[0] = vel[0];
          _val[1] = vel[1];
          _val[2] = vel[2];
        }
#endif
        /* Dual surface associated to vertices */
        cs_real_t *_v_surf = nullptr;

        /* Transform face flux to vertex velocities */
        cs_real_3_t *_mesh_vel = nullptr;

        CS_MALLOC(_v_surf, m->n_vertices, cs_real_t);
        CS_MALLOC(_mesh_vel, m->n_vertices, cs_real_3_t);

        /* Initialize */
        cs_array_real_fill_zero(3 * m->n_vertices, (cs_real_t *)_mesh_vel);
        cs_array_real_fill_zero(m->n_vertices, _v_surf);

        /* First Loop over boundary faces
         * to compute face portion associated to the vertex
         * and add the portion of face velocity to the vertex velocity */
        for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {

          const cs_lnum_t face_id = z->elt_ids[elt_id];

          cs_real_t normal[3]; // copy in case cs_real_t != cs_nreal_t
          for (cs_lnum_t i = 0; i < 3; i++)
            normal[i] = b_face_u_normal[face_id][i];

          const cs_lnum_t s = m->b_face_vtx_idx[face_id];
          const cs_lnum_t e = m->b_face_vtx_idx[face_id+1];

          /* Compute the portion of surface associated to v_id_1 */
          for (cs_lnum_t k = s; k < e; k++) {

            const cs_lnum_t k_1 = (k+1 < e) ? k+1 : k+1-(e-s);
            const cs_lnum_t k_2 = (k+2 < e) ? k+2 : k+2-(e-s);
            const cs_lnum_t v_id0 = m->b_face_vtx_lst[k];
            const cs_lnum_t v_id1 = m->b_face_vtx_lst[k_1];
            const cs_lnum_t v_id2 = m->b_face_vtx_lst[k_2];

            cs_real_3_t v0v1 = {
              vtx_coord[v_id1][0] - vtx_coord[v_id0][0],
              vtx_coord[v_id1][1] - vtx_coord[v_id0][1],
              vtx_coord[v_id1][2] - vtx_coord[v_id0][2]
            };
            cs_real_3_t v1v2 = {
              vtx_coord[v_id2][0] - vtx_coord[v_id1][0],
              vtx_coord[v_id2][1] - vtx_coord[v_id1][1],
              vtx_coord[v_id2][2] - vtx_coord[v_id1][2]
            };
            cs_real_3_t v1_cog = {
              b_face_cog[face_id][0] - vtx_coord[v_id1][0],
              b_face_cog[face_id][1] - vtx_coord[v_id1][1],
              b_face_cog[face_id][2] - vtx_coord[v_id1][2]
            };

            /* Portion of the surface associated to the vertex
             * projected in the normal direction */
            cs_real_t portion_surf = 0.25 * (
                  cs_math_3_triple_product(v0v1, v1_cog, normal)
                + cs_math_3_triple_product(v1v2, v1_cog, normal));

            _v_surf[v_id1] += portion_surf;

            /* Portion of the face velocity is added to the node velocity
             * Warning: bc_vals is not interleaved */
            for (int i = 0; i < 3; i++)
              _mesh_vel[v_id1][i] += bc_vals[elt_id + i * z->n_elts] * portion_surf;

          }

        } /* Loop on selected border faces */

        /* Handle parallelism */
        if (m->vtx_interfaces != nullptr) {
          cs_interface_set_sum(m->vtx_interfaces,
                               m->n_vertices,
                               3,
                               true,
                               CS_REAL_TYPE,
                               _mesh_vel);

          cs_interface_set_sum(m->vtx_interfaces,
                               m->n_vertices,
                               1,
                               true,
                               CS_REAL_TYPE,
                               _v_surf);
        }

        /* Loop on selected border vertices */
        for (cs_lnum_t i = 0; i < _cdo_bc->n_vertices[select_id]; i++) {

          const cs_lnum_t v_id = _cdo_bc->vtx_select[select_id][i];
          const double  invsurf = 1./_v_surf[v_id];

          cs_real_t  *_val = _cdo_bc->vtx_values + 3*v_id;

          for (int k = 0; k < 3; k++) {
            _mesh_vel[v_id][k] *= invsurf;
            _val[k] = _mesh_vel[v_id][k];
          }
        }

        /* Loop over boundary faces and impose node mesh velocity */
        for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {
          const cs_lnum_t face_id = z->elt_ids[elt_id];

          /* fluid velocity BC */
          ale_bc_type[face_id] = CS_BOUNDARY_ALE_IMPOSED_VEL;
          for (int d = 0; d < 3; d++)
            b_fluid_vel[face_id][d] = bc_vals[elt_id + d * z->n_elts];
        }

        /* Free memory */
        CS_FREE(bc_vals);
        CS_FREE(_mesh_vel);
        CS_FREE(_v_surf);

        select_id++;
      }
      break;

    case CS_BOUNDARY_ALE_IMPOSED_DISP:
    case CS_BOUNDARY_ALE_INTERNAL_COUPLING:
    case CS_BOUNDARY_ALE_EXTERNAL_COUPLING: {
      const cs_real_3_t *restrict disale = (const cs_real_3_t *)(f_displ->val);
      const cs_real_t invdt = 1. / domain->time_step->dt_ref; /* JB: dt[0] ? */

      assert(select_id < _cdo_bc->n_selections);

      /* Loop on selected border vertices */
      for (cs_lnum_t i = 0; i < _cdo_bc->n_vertices[select_id]; i++) {
        const cs_lnum_t  v_id           = _cdo_bc->vtx_select[select_id][i];
        const cs_real_t *_dpl           = disale[v_id];
        const cs_real_t *restrict _xyz  = vtx_coord[v_id];
        const cs_real_t *restrict _xyz0 = _vtx_coord0[v_id];

        cs_real_t *_val = _cdo_bc->vtx_values + 3 * v_id;

        for (int k = 0; k < 3; k++)
          _val[k] = (_dpl[k] + _xyz0[k] - _xyz[k]) * invdt;

      } /* Loop on selected vertices */

        /* Loop over boundary faces for the fluid velocity */
        for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {
          const cs_lnum_t face_id = z->elt_ids[elt_id];

          ale_bc_type[face_id] = CS_BOUNDARY_ALE_IMPOSED_VEL;

          cs_real_t normal[3]; // copy in case cs_real_t != cs_nreal_t
          for (cs_lnum_t i = 0; i < 3; i++)
            normal[i] = b_face_u_normal[face_id][i];

          /* Normal direction is given by the gravity */
          const cs_real_t dsurf = 1. / mq->b_face_surf[face_id];

          b_fluid_vel[face_id][0] = 0.;
          b_fluid_vel[face_id][1] = 0.;
          b_fluid_vel[face_id][2] = 0.;

          const cs_lnum_t s = m->b_face_vtx_idx[face_id];
          const cs_lnum_t e = m->b_face_vtx_idx[face_id+1];

          /* Compute the portion of surface associated to v_id_1 */
          for (cs_lnum_t k = s; k < e; k++) {

            const cs_lnum_t k_1 = (k+1 < e) ? k+1 : k+1-(e-s);
            const cs_lnum_t k_2 = (k+2 < e) ? k+2 : k+2-(e-s);
            const cs_lnum_t v_id0 = m->b_face_vtx_lst[k];
            const cs_lnum_t v_id1 = m->b_face_vtx_lst[k_1];
            const cs_lnum_t v_id2 = m->b_face_vtx_lst[k_2];

            cs_real_3_t v0v1 = {
              vtx_coord[v_id1][0] - vtx_coord[v_id0][0],
              vtx_coord[v_id1][1] - vtx_coord[v_id0][1],
              vtx_coord[v_id1][2] - vtx_coord[v_id0][2],
            };
            cs_real_3_t v1v2 = {
              vtx_coord[v_id2][0] - vtx_coord[v_id1][0],
              vtx_coord[v_id2][1] - vtx_coord[v_id1][1],
              vtx_coord[v_id2][2] - vtx_coord[v_id1][2],
            };
            cs_real_3_t v1_cog = {
              b_face_cog[face_id][0] - vtx_coord[v_id1][0],
              b_face_cog[face_id][1] - vtx_coord[v_id1][1],
              b_face_cog[face_id][2] - vtx_coord[v_id1][2],
            };

            /* Portion of the surface associated to the vertex
             * projected in the normal direction */
            cs_real_t portion_surf = 0.25 * dsurf  * (
                cs_math_3_triple_product(v0v1, v1_cog, normal)
                + cs_math_3_triple_product(v1v2, v1_cog, normal));

            cs_real_t *_val = _cdo_bc->vtx_values + 3*v_id1;

            b_fluid_vel[face_id][0] += portion_surf * _val[0];
            b_fluid_vel[face_id][1] += portion_surf * _val[1];
            b_fluid_vel[face_id][2] += portion_surf * _val[1];
          }

        }

        select_id++;
      }
      break;

      /* Treated elsewhere, only increased selected_id */
      case CS_BOUNDARY_ALE_FREE_SURFACE:
        assert(select_id < _cdo_bc->n_selections);
        select_id++;
        break;

    default:
      break; /* Nothing to do */
    }

  } /* Loop on ALE boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function updates free surface ALE BCs for CDO.
 *        These BCs are required after the fluid is solved.
 *
 * \param[in]     domain        domain quantities
 */
/*----------------------------------------------------------------------------*/

static void
_update_bcs_free_surface(const cs_domain_t  *domain)
{
  /* Only a part of the boundaries has to be updated */
  int  select_id = 0;
  for (int b_id = 0; b_id < domain->ale_boundaries->n_boundaries; b_id++) {

    const int z_id = domain->ale_boundaries->zone_ids[b_id];
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);

    switch(domain->ale_boundaries->types[b_id]) {

      /* Treated elsewhere, only increased selected_id */
      case CS_BOUNDARY_ALE_IMPOSED_VEL:
        select_id++;
        break;

      /* Treated elsewhere, only increased selected_id */
      case CS_BOUNDARY_ALE_IMPOSED_DISP:
      case CS_BOUNDARY_ALE_INTERNAL_COUPLING:
      case CS_BOUNDARY_ALE_EXTERNAL_COUPLING:
        select_id++;
        break;

      case CS_BOUNDARY_ALE_FREE_SURFACE:
        assert(select_id < _cdo_bc->n_selections);
        _free_surface(domain, z, select_id);
        select_id++;
        break;

      default:
        break; /* Nothing to do */
    }

  } /* Loop on ALE boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This subroutine performs the solving of a Poisson equation
 *        on the mesh velocity for ALE module. It also updates the mesh
 *        displacement so that it can be used to update mass fluxes (due to
 *        mesh displacement).
 *
 * \param[in]     domain        domain quantities
 * \param[in]     impale        Indicator for fixed node displacement
 */
/*----------------------------------------------------------------------------*/

static void
_ale_solve_poisson_cdo(const cs_domain_t  *domain,
                       const int           impale[])
{
  const cs_mesh_t  *m = domain->mesh;

  /* Build and solve equation on the mesh velocity */

  cs_equation_t  *eq = cs_equation_by_name("mesh_velocity");

  /* Update the values of boundary mesh vertices */

  _update_bcs_free_surface(domain);

  /* Compute the Poisson equation on the original mesh */

  cs_real_3_t *vtx_coord = (cs_real_3_t *)m->vtx_coord;
  cs_real_3_t *vtx_coord0 = (cs_real_3_t *)(cs_field_by_name("vtx_coord0")->val);
  const cs_lnum_t n_vertices = m->n_vertices;
  cs_mesh_quantities_t *mq = domain->mesh_quantities;

  /* Back to original mesh */

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    for (cs_lnum_t idim = 0; idim < 3; idim++) {
      vtx_coord[v_id][idim] = vtx_coord0[v_id][idim];
    }
  }

  cs_ale_update_mesh_quantities(&(mq->min_vol), &(mq->max_vol), &(mq->tot_vol));

  /* Solve the algebraic system */

  cs_equation_solve_steady_state(m, eq);

  /* Retrieving fields */

  cs_field_t  *f_displ = cs_field_by_name("mesh_displacement");
  cs_real_3_t *disale = (cs_real_3_t *)(f_displ->val);
  cs_real_3_t *disala = (cs_real_3_t *)(f_displ->val_pre);
  cs_real_3_t *m_vel = (cs_real_3_t *)(cs_field_by_name("mesh_velocity")->val);

  /* Back to mesh at time n */

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    for (cs_lnum_t idim = 0; idim < 3; idim++) {
      vtx_coord[v_id][idim] = vtx_coord0[v_id][idim] + disale[v_id][idim];
    }
  }

  cs_ale_update_mesh_quantities(&(mq->min_vol), &(mq->max_vol), &(mq->tot_vol));

  for (cs_lnum_t v = 0; v < m->n_vertices; v++) {
    if (impale[v] == 0) {
      for (int k = 0; k < 3; k++)
        disale[v][k] = disala[v][k] + m_vel[v][k] * domain->time_step->dt_ref;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a Poisson equation on the mesh velocity in ALE framework.
 *
 * It also updates the mesh displacement
 * so that it can be used to update mass fluxes (due to mesh displacement).
 *
 * \param[in]       domain           domain quantities
 * \param[in]       iterns           Navier-Stokes iteration number
 * \param[in]       impale           Indicator for fixed node displacement
 * \param[in]       ale_bc_type      Type of boundary for ALE
 * \param[in]       b_rho            density at boundary
 * \param[in]       i_mass_flux      interior faces mass flux
 * \param[in]       b_mass_flux      boundary faces mass flux
 */
/*----------------------------------------------------------------------------*/

static void
_ale_solve_poisson_legacy(const cs_domain_t *domain,
                          const int          iterns,
                          const int         *impale,
                          const int         *ale_bc_type,
                          const cs_real_t    b_rho[],
                          const cs_real_t    i_mass_flux[],
                          const cs_real_t    b_mass_flux[])
{
  const cs_mesh_t *m = domain->mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = (const cs_lnum_t *)m->b_face_cells;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;
  const cs_real_t *b_dist = mq->b_dist;
  const cs_real_t *restrict b_face_surf = mq->b_face_surf;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_lnum_t need_solve = 0;
  for (cs_lnum_t v = 0; v < n_vertices; v++) {
    if (impale[v] == 0) {
      need_solve += 1;
    }
  }

  cs_parall_sum(1, CS_LNUM_TYPE, &need_solve);

  if (need_solve == 0) {
    bft_printf("Mesh velocity is not computed as the displacement "
               "is imposed for each vertice (%d)\n", need_solve);
    cs_array_real_fill_zero(3*n_cells_ext, (cs_real_t *)CS_F_(mesh_u)->val);
    return;
  }

  /* 1. Initialization */
  cs_real_3_t rinfiv = { cs_math_infinite_r,
                         cs_math_infinite_r,
                         cs_math_infinite_r };

  cs_real_3_t *smbr = nullptr;
  cs_real_33_t *fimp = nullptr;
  CS_MALLOC_HD(smbr, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(fimp, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  cs_real_3_t *mshvel = (cs_real_3_t *)CS_F_(mesh_u)->val;
  cs_real_3_t *mshvela = (cs_real_3_t *)CS_F_(mesh_u)->val_pre;

  cs_field_t  *f_displ = cs_field_by_name("mesh_displacement");
  cs_real_3_t *disale = (cs_real_3_t *)(f_displ->val);
  cs_real_3_t *disala = (cs_real_3_t *)(f_displ->val_pre);

  cs_equation_param_t *eqp =
    cs_field_get_equation_param(CS_F_(mesh_u));

  if (eqp->verbosity >= 1)
    bft_printf("\n   ** SOLVING MESH VELOCITY\n"
               "      ---------------------\n");

  /* We compute the boundary condition on the mesh velocity at the free surface
   * from the new mass flux. */

  cs_field_bc_coeffs_t *bc_coeffs = CS_F_(mesh_u)->bc_coeffs;

  int idftnp = eqp->idften;

  /* The mesh moves in the direction of the gravity in case of free-surface */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (ale_bc_type[face_id] == CS_FREE_SURFACE) {
      cs_lnum_t cell_id = b_face_cells[face_id];
      cs_real_t distbf = b_dist[face_id];

      cs_real_6_t hintt = {0., 0., 0., 0., 0., 0.};
      if (idftnp & CS_ISOTROPIC_DIFFUSION) {
        for (int isou = 0; isou < 3; isou++)
          hintt[isou] = CS_F_(vism)->val[cell_id] / distbf;
      }
      else if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
        for (int isou = 0; isou < 6; isou++)
          hintt[isou] = CS_F_(vism)->val[6*cell_id+isou] / distbf;
      }

      cs_real_t prosrf =   cs_math_3_dot_product(grav, b_face_u_normal[face_id])
                         * b_face_surf[face_id];

      cs_real_t pimpv[3];
      for (int i = 0; i < 3; i++)
        pimpv[i] = grav[i]*b_mass_flux[face_id]/(b_rho[face_id]*prosrf);

      cs_boundary_conditions_set_dirichlet_vector_aniso(face_id,
                                                        bc_coeffs,
                                                        pimpv,
                                                        hintt,
                                                        rinfiv);
    }
  }

  /* 2. Solving of the mesh velocity equation */

  if (eqp->verbosity >= 1)
    bft_printf("\n\n           SOLVING VARIABLE %s\n\n",
               CS_F_(mesh_u)->name);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (int isou = 0; isou < 3; isou++) {
      smbr[cell_id][isou] = 0.;
      for (int jsou = 0; jsou < 3; jsou++)
        fimp[cell_id][jsou][isou] = 0.;
    }
  }

  cs_real_t *i_visc = nullptr, *b_visc = nullptr;

  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);

  if (idftnp & CS_ISOTROPIC_DIFFUSION) {
    CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);

    cs_face_viscosity(m,
                      mq,
                      eqp->imvisf,
                      CS_F_(vism)->val,
                      i_visc,
                      b_visc);

  }
  else if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
    CS_MALLOC_HD(i_visc, 9 * n_i_faces, cs_real_t, cs_alloc_mode);

    cs_face_anisotropic_viscosity_vector(m,
                                         mq,
                                         eqp->imvisf,
                                         (cs_real_6_t *)CS_F_(vism)->val,
                                         (cs_real_33_t *)i_visc,
                                         b_visc);
  }

  cs_equation_iterative_solve_vector(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     CS_F_(mesh_u)->id,
                                     CS_F_(mesh_u)->name,
                                     0, /* ivisep */
                                     0, /* iescap */
                                     eqp,
                                     (const cs_real_3_t *)mshvela,
                                     (const cs_real_3_t *)mshvela,
                                     bc_coeffs,
                                     i_mass_flux,
                                     b_mass_flux,
                                     i_visc,
                                     b_visc,
                                     i_visc,
                                     b_visc,
                                     nullptr, /* i_secvis */
                                     nullptr, /* b_secvis */
                                     nullptr, /* viscel */
                                     nullptr, /* weighf */
                                     nullptr, /* weighb */
                                     0,    /* icvflv */
                                     nullptr, /* icvfli */
                                     (cs_real_33_t *)fimp,
                                     smbr,
                                     mshvel,
                                     nullptr); /* eswork */

  /* Free memory */
  CS_FREE_HD(smbr);
  CS_FREE_HD(fimp);
  CS_FREE_HD(i_visc);
  CS_FREE_HD(b_visc);

  /* 3. Update nodes displacement */

  cs_real_3_t *dproj;
  cs_real_33_t *gradm;

  /* Allocate a temporary array */
  CS_MALLOC_HD(dproj, n_vertices, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(gradm, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  bool use_previous_t = false;
  int inc = 1;

  cs_field_gradient_vector(CS_F_(mesh_u),
                           use_previous_t,
                           inc,
                           gradm);

  cs_ale_project_displacement(ale_bc_type,
                              (const cs_real_3_t *)mshvel,
                              (const cs_real_33_t *)gradm,
                              (const cs_real_3_t *)bc_coeffs->a,
                              (const cs_real_33_t *)bc_coeffs->b,
                              (const cs_real_t *)CS_F_(dt)->val,
                              dproj);

  /* FIXME : warning if nterup > 1, use itrale ? */
  /* Update mesh displacement only where it is not
     imposed by the user (ie when impale <> 1) */
  for (cs_lnum_t v = 0; v < n_vertices; v++) {
    if (impale[v] == 0) {
      for (int isou = 0; isou < 3; isou++)
        disale[v][isou] = disala[v][isou] + dproj[v][isou];
    }
  }

  /* Free memory */
  CS_FREE_HD(dproj);
  CS_FREE_HD(gradm);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free ALE boundary condition mappings.
 */
/*----------------------------------------------------------------------------*/

static void
_ale_free(void)
{
  CS_FREE(cs_glob_ale_data->impale);
  CS_FREE(cs_glob_ale_data->bc_type);
  CS_FREE(cs_glob_ale_data->b_mass_flux_ale);
  CS_FREE(cs_glob_ale_data->i_mass_flux_ale);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute of ALE volumic flux from displacement and deduced volume
 *        at time step n+1.
 *
 * \param[in, out] domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_compute_volume_from_displacement(cs_domain_t *domain) {
  if (cs_glob_ale == CS_ALE_NONE)
    return;

  const cs_lnum_t *i_face_vtx_idx = domain->mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx_lst = domain->mesh->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = domain->mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx_lst = domain->mesh->b_face_vtx_lst;

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)(domain->mesh->vtx_coord);
  const cs_real_3_t *disale
    = (const cs_real_3_t *)cs_field_by_name("mesh_displacement")->val;
  const cs_real_3_t *xyzno0
    = (const cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

  const cs_lnum_t n_i_faces = domain->mesh->n_i_faces;
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

  cs_real_t *i_mass_flux_ale = cs_glob_ale_data->i_mass_flux_ale;
  cs_real_t *b_mass_flux_ale = cs_glob_ale_data->b_mass_flux_ale;

  /* Computation of volumic flux
     --------------------------- */

  /* Computation of delta_volume due to nodes displacement
     for each internal faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_lnum_t s_id = i_face_vtx_idx[f_id];
    const cs_lnum_t e_id = i_face_vtx_idx[f_id+1];
    const cs_lnum_t n_face_vertices = e_id - s_id;

    cs_real_t *_vc[10][3];
    cs_real_3_t *vc = (cs_real_3_t *)_vc;
    if (n_face_vertices > 10)
      CS_MALLOC(vc, n_face_vertices, cs_real_3_t);

    /* Init flux */
    cs_real_t dvol = 0.;

    /* Compute center of gravity and surface of the initial surface */
    for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++) {
      cs_lnum_t j = i_v - s_id;
      cs_lnum_t inod = i_face_vtx_lst[i_v];

      for (cs_lnum_t k = 0; k < 3; k++)
        vc[j][k] = vtx_coord[inod][k];
    }

    cs_real_t cogloc[3], surfloc[3];
    _compute_cog_surf_from_nodes(n_face_vertices, vc, cogloc, surfloc);

    /* Flux contribution */
    dvol -= cs_math_3_dot_product(cogloc, surfloc);

    /* Compute center of gravity and surface of the final surface
       after displacement */
    for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++) {
      cs_lnum_t j = i_v - s_id;
      cs_lnum_t inod = i_face_vtx_lst[i_v];

      for (cs_lnum_t k = 0; k < 3; k++)
        vc[j][k] = xyzno0[inod][k] + disale[inod][k];
    }
    _compute_cog_surf_from_nodes(n_face_vertices, vc, cogloc, surfloc);

    /* Flux contribution */
    dvol += cs_math_3_dot_product(cogloc, surfloc);

    /* Compute center of gravity and surface of the annex surfaces */
    for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++) {
      cs_lnum_t inod1, inod2;

      if (i_v < e_id-1) {
        inod1 = i_face_vtx_lst[i_v];
        inod2 = i_face_vtx_lst[i_v+1];
      }
      else {
        inod1 = i_face_vtx_lst[i_v];
        inod2 = i_face_vtx_lst[s_id];
      }
      int nb_node_lat = 4;
      for (cs_lnum_t k = 0; k < 3; k++) {
        vc[0][k] = vtx_coord[inod1][k];
        vc[1][k] = xyzno0[inod1][k] + disale[inod1][k];
        vc[2][k] = xyzno0[inod2][k] + disale[inod2][k];
        vc[3][k] = vtx_coord[inod2][k];
      }

      _compute_cog_surf_from_nodes(nb_node_lat, vc, cogloc, surfloc);

      /* Flux contribution */
      cs_real_t contrib_lat = cs_math_3_dot_product(cogloc, surfloc);
      dvol -= contrib_lat;
    }

    /* Finalization of the variation of volume for the internal face */
    dvol /= 3.;
    i_mass_flux_ale[f_id] = dvol;

    if (vc != (cs_real_3_t *)_vc)
      CS_FREE(vc);
  }

  /* Computation of the delta_volume due to nodes displacement
     for each boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_lnum_t s_id = b_face_vtx_idx[f_id];
    const cs_lnum_t e_id = b_face_vtx_idx[f_id+1];
    const cs_lnum_t n_face_vertices = e_id - s_id;

    cs_real_t *_vc[10][3];
    cs_real_3_t *vc = (cs_real_3_t *)_vc;
    if (n_face_vertices > 10)
      CS_MALLOC(vc, n_face_vertices, cs_real_3_t);

    /* init flux */
    cs_real_t dvol = 0.;

    /* Compute center of gravity and surface of the initial surface */
    for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++) {
      cs_lnum_t j = i_v - s_id;
      cs_lnum_t inod = b_face_vtx_lst[i_v];
      for (cs_lnum_t k = 0; k < 3; k++)
        vc[j][k] = vtx_coord[inod][k];
    }

    cs_real_t cogloc[3], surfloc[3];
    _compute_cog_surf_from_nodes(n_face_vertices, vc, cogloc, surfloc);

    /* Flux contribution */
    dvol -= cs_math_3_dot_product(cogloc, surfloc);

    /* Compute center of gravity and surface of the final surface */
    for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++) {
      cs_lnum_t j = i_v - s_id;
      cs_lnum_t inod = b_face_vtx_lst[i_v];
      for (cs_lnum_t k = 0; k < 3; k++)
        vc[j][k] = xyzno0[inod][k] + disale[inod][k];
    }
    _compute_cog_surf_from_nodes(n_face_vertices, vc, cogloc, surfloc);

    /* Flux contribution */
    dvol += cs_math_3_dot_product(cogloc, surfloc);

    /* Compute center of gravity and surface of annex surfaces */
    cs_lnum_t inod1, inod2;
    for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++) {
      if (i_v < e_id-1) {
        inod1 = b_face_vtx_lst[i_v];
        inod2 = b_face_vtx_lst[i_v+1];
      }
      else {
        inod1 = b_face_vtx_lst[i_v];
        inod2 = b_face_vtx_lst[s_id];
      }
      int nb_node_lat = 4;
      for (cs_lnum_t k = 0; k < 3; k++) {
        vc[0][k] = vtx_coord[inod1][k];
        vc[1][k] = xyzno0[inod1][k] + disale[inod1][k];
        vc[2][k] = xyzno0[inod2][k] + disale[inod2][k];
        vc[3][k] = vtx_coord[inod2][k];
      }

      _compute_cog_surf_from_nodes(nb_node_lat, vc, cogloc, surfloc);

      /* Flux contribution */
      cs_real_t contrib_lat = cs_math_3_dot_product(cogloc, surfloc);
      dvol -= contrib_lat;
    }

    /* Finalization of the volume variation for the boundary face */
    dvol /= 3.;
    b_mass_flux_ale[f_id] = dvol;

    if (vc != (cs_real_3_t *)_vc)
      CS_FREE(vc);
  }

  /* Compute new volume(n+1) from volumic flux (user2) */

  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)domain->mesh->i_face_cells;

  const cs_real_t *voln = domain->mesh_quantities->cell_vol;
  cs_real_t *volnp1 = cs_field_by_name("cell_vol")->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    volnp1[c_id] = voln[c_id];

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];
    cs_real_t dvol = i_mass_flux_ale[f_id];
    volnp1[c_id1] += dvol;
    volnp1[c_id2] -= dvol;
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t dvol = b_mass_flux_ale[f_id];
    volnp1[c_id] += dvol;
  }

  cs_halo_sync_var(domain->mesh->halo, CS_HALO_EXTENDED, volnp1);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief In the ALE framework, update mass flux by adding mesh velocity.
 *
 * \param[in]      m       pointer to associated mesh structure
 * \param[in]      mq      pointer to associated mesh quantities structure
 * \param[in]      dt      time step at cells
 * \param[in]      crom    density at cells
 * \param[in]      brom    density at boundary faces
 * \param[in, out] imasfl  interior face mass flux
 * \param[in, out] bmasfl  boundary face mass flux
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_velocity_mass_flux(const cs_mesh_t             *m,
                           const cs_mesh_quantities_t  *mq,
                           const cs_real_t              dt[],
                           const cs_real_t              crom[],
                           const cs_real_t              brom[],
                           cs_real_t                    imasfl[],
                           cs_real_t                    bmasfl[])
{
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  const cs_lnum_t *i_face_vtx_idx = m->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx_lst = m->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = m->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx_lst = m->b_face_vtx_lst;

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)(m->vtx_coord);
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;
  const cs_real_t *b_face_surf = mq->b_face_surf;
  const cs_real_t *i_face_surf = mq->i_face_surf;

  const cs_real_3_t *mshvel = (const cs_real_3_t *)CS_F_(mesh_u)->val;

  const cs_real_3_t *xyzno0
    = (const cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

  const cs_real_3_t *disale
    = (const cs_real_3_t *)cs_field_by_name("mesh_displacement")->val;

  cs_dispatch_context ctx;

  if (cs_glob_space_disc->iflxmw > 0) {

    /* One temporary array needed for internal faces,
     * in case some internal vertices are moved directly by the user */

    cs_real_t *intflx = nullptr, *bouflx = nullptr;
    CS_MALLOC_HD(intflx, n_i_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bouflx, n_b_faces, cs_real_t, cs_alloc_mode);

    cs_field_bc_coeffs_t *bc_coeffs_ale = CS_F_(mesh_u)->bc_coeffs;

    const cs_equation_param_t *eqp_mesh
      = cs_field_get_equation_param_const(CS_F_(mesh_u));

    cs_mass_flux(m,
                 mq,
                 CS_F_(mesh_u)->id,
                 1,  /* itypfl */
                 1,  /* iflmb0 */
                 1,  /* init */
                 1,  /* inc */
                 eqp_mesh->imrgra,
                 eqp_mesh->nswrgr,
                 static_cast<cs_gradient_limit_t>(eqp_mesh->imligr),
                 eqp_mesh->verbosity,
                 eqp_mesh->epsrgr,
                 eqp_mesh->climgr,
                 crom, brom,
                 mshvel,
                 bc_coeffs_ale,
                 intflx, bouflx);

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      bmasfl[face_id] -= bouflx[face_id];
    });

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      imasfl[face_id] -= intflx[face_id];
    });

    ctx.wait();

    CS_FREE_HD(intflx);
    CS_FREE_HD(bouflx);
  }

  /* Here we need of the opposite of the mesh velocity. */

  else { /* if (cs_glob_space_disc->iflxmw == 0) */

    /* Compute the mass flux using the nodes displacement */

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      cs_real_t disp_fac[3] = {0, 0, 0};
      const cs_lnum_t s_id = b_face_vtx_idx[face_id];
      const cs_lnum_t e_id = b_face_vtx_idx[face_id+1];
      const cs_lnum_t icpt = e_id - s_id;
      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t inod = b_face_vtx_lst[ii];

        for (cs_lnum_t jj = 0; jj < 3; jj++)
          disp_fac[jj] +=   disale[inod][jj]
                          - (vtx_coord[inod][jj] - xyzno0[inod][jj]);
      }
      const cs_lnum_t c_id = b_face_cells[face_id];

      const cs_real_t *n_u = b_face_u_normal[face_id];
      const cs_real_t dot
        = b_face_surf[face_id] * cs_math_3_dot_product(disp_fac, n_u);

      bmasfl[face_id] -= brom[face_id] * dot / dt[c_id] / icpt;
    });

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      cs_real_t disp_fac[3] = {0, 0, 0};
      const cs_lnum_t s_id = i_face_vtx_idx[face_id];
      const cs_lnum_t e_id = i_face_vtx_idx[face_id+1];
      const cs_lnum_t icpt = e_id - s_id;
      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t inod = i_face_vtx_lst[ii];

        for (cs_lnum_t jj = 0; jj < 3; jj++)
          disp_fac[jj] +=   disale[inod][jj]
                          - (vtx_coord[inod][jj] - xyzno0[inod][jj]);
      }

      /* For inner vertices, the mass flux due to the mesh displacement is
       * recomputed from the nodes displacement */
      const cs_lnum_t c_id1 = i_face_cells[face_id][0];
      const cs_lnum_t c_id2 = i_face_cells[face_id][1];
      const cs_real_t dtfac = 0.5*(dt[c_id1] + dt[c_id2]);
      const cs_real_t rhofac = 0.5*(crom[c_id1] + crom[c_id2]);

      const cs_real_t *n_u = i_face_u_normal[face_id];
      const cs_real_t dot
        = i_face_surf[face_id] * cs_math_3_dot_product(disp_fac, n_u);

      imasfl[face_id] -= rhofac * dot / dtfac / icpt;
    });

    ctx.wait();
  }
}

/*----------------------------------------------------------------------------
 * Compute boundary condition code for ALE
 *----------------------------------------------------------------------------*/

void
cs_boundary_condition_ale_type(const cs_mesh_t            *m,
                               const cs_mesh_quantities_t *mq,
                               const bool                  init,
                               const cs_real_t             dt[],
                               const int                   bc_type[],
                               cs_real_t                  *rcodcl1_vel)
{
  const cs_lnum_t  n_b_faces    = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;

  /* Initialization
   * -------------- */

  int *impale      = cs_glob_ale_data->impale;
  int *ale_bc_type = cs_glob_ale_data->bc_type;

  cs_real_3_t *disale =
    (cs_real_3_t *)cs_field_by_name("mesh_displacement")->val;

  int       *icodcl_mesh_u   = nullptr;
  cs_real_t *_rcodcl1_mesh_u = nullptr;
  cs_real_t *rcodcl1_mesh_u  = nullptr;

  if (CS_F_(mesh_u)->bc_coeffs != nullptr) {
    icodcl_mesh_u  = CS_F_(mesh_u)->bc_coeffs->icodcl;
    rcodcl1_mesh_u = CS_F_(mesh_u)->bc_coeffs->rcodcl1;
  }

  if (cs_glob_ale == CS_ALE_CDO) {
    const int size_uma = CS_F_(mesh_u)->dim * n_b_faces;
    CS_MALLOC(_rcodcl1_mesh_u, size_uma, cs_real_t);
    rcodcl1_mesh_u = _rcodcl1_mesh_u;

    cs_real_3_t *b_fluid_vel = nullptr;
    CS_MALLOC(b_fluid_vel, n_b_faces, cs_real_3_t);

    cs_array_real_fill_zero(3 * n_b_faces, (cs_real_t *)b_fluid_vel);

    cs_ale_update_bcs(ale_bc_type, b_fluid_vel);

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        rcodcl1_mesh_u[n_b_faces * ii + face_id] = b_fluid_vel[face_id][ii];
    }

    CS_FREE(b_fluid_vel);
  }

  assert(rcodcl1_mesh_u != nullptr);
#pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      if (rcodcl1_mesh_u[n_b_faces * ii + face_id] > cs_math_infinite_r * 0.5)
        rcodcl1_mesh_u[n_b_faces * ii + face_id] = 0.;
    }
  }

  /* Check the consistency of BC types
   * --------------------------------- */

  /* When using CDO solver, no need for checks. */
  if (!(cs_glob_ale == CS_ALE_CDO)) {

    int ierror = 0;

#pragma omp parallel for if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (ale_bc_type[face_id] != 0 &&
          ale_bc_type[face_id] != CS_BOUNDARY_ALE_FIXED &&
          ale_bc_type[face_id] != CS_BOUNDARY_ALE_SLIDING &&
          ale_bc_type[face_id] != CS_BOUNDARY_ALE_FREE_SURFACE &&
          ale_bc_type[face_id] != CS_BOUNDARY_ALE_IMPOSED_VEL) {
        if (ale_bc_type[face_id] > 0) {
          ale_bc_type[face_id] = -ale_bc_type[face_id];
        }
        ierror++;
      }
    }

    cs_parall_max(1, CS_INT_TYPE, &ierror);

    if (ierror != 0) {
      bft_printf("ALE METHOD\n"
                 "At least one boundary face has an unknown boundary type.\n"
                 "  The calculation will not be run."
                 "Check boundary condition definitions.");
      cs_boundary_conditions_error(ale_bc_type, nullptr);
    }

    /* Conversion into BC and values
     *-------------------------------*/

    cs_real_3_t *xyzno0 = (cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;
    const cs_real_3_t *xyzno = (const cs_real_3_t *)m->vtx_coord;

    /* If all the nodes of a face have an imposed displacement, rcodcl is
       computed or overwritten, ale bc type is therefore
       CS_BOUNDARY_ALE_IMPOSED_VEL */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      int iecrw = 0, icpt = 0;
      cs_real_t ddep_xyz[3] = { 0., 0., 0. };

      const cs_lnum_t s = m->b_face_vtx_idx[face_id];
      const cs_lnum_t e = m->b_face_vtx_idx[face_id + 1];

      for (cs_lnum_t ii = s; ii < e; ii++) {
        const cs_lnum_t vtx_id = m->b_face_vtx_lst[ii];
        if (impale[vtx_id] == 0) {
          iecrw++;
        }
        icpt++;
        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          ddep_xyz[jj] +=
            disale[vtx_id][jj] + xyzno0[vtx_id][jj] - xyzno[vtx_id][jj];
        }
      }

      if (iecrw == 0 && ale_bc_type[face_id] != CS_BOUNDARY_ALE_SLIDING) {
        const cs_lnum_t c_id = b_face_cells[face_id];
        cs_real_t coeff = 1. / (dt[c_id] * icpt);

        ale_bc_type[face_id] = CS_BOUNDARY_ALE_IMPOSED_VEL;
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          rcodcl1_mesh_u[n_b_faces * ii + face_id] = ddep_xyz[ii] * coeff;
      }
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      int icpt = 0;

      /* Fixed faces: we impose that nodes are also fixed */
      if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_FIXED) {
        icpt = 0;
        if (icodcl_mesh_u[face_id] == 0) {
          icpt++;
          icodcl_mesh_u[face_id] = 1;
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            rcodcl1_mesh_u[n_b_faces * ii + face_id] = 0.;
        }

        if (icpt) {
          const cs_lnum_t s = m->b_face_vtx_idx[face_id];
          const cs_lnum_t e = m->b_face_vtx_idx[face_id + 1];
          for (cs_lnum_t ii = s; ii < e; ii++) {
            const cs_lnum_t vtx_id = m->b_face_vtx_lst[ii];
            if (impale[vtx_id] == 0) {
              impale[vtx_id] = 1;
              for (cs_lnum_t jj = 0; jj < 3; jj++)
                disale[vtx_id][jj] = 0;
            }
          }
        }
      }

      /* Sliding face */
      else if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_SLIDING) {
        if (icodcl_mesh_u[face_id] == 0)
          icodcl_mesh_u[face_id] = 4;
      }

      /* Imposed mesh velocity face */
      else if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_IMPOSED_VEL) {
        if (icodcl_mesh_u[face_id] == 0)
          icodcl_mesh_u[face_id] = 1;
      }

      /* Free surface face: the mesh velocity is imposed by the mass flux */
      else if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_FREE_SURFACE) {
        if (icodcl_mesh_u[face_id] == 0)
          icodcl_mesh_u[face_id] = 1;
      }
    }

    /* Check icodcl consistency
     * ------------------------ */

    int irkerr = -1;
    int icoder[2] = {-1, -1};
    int ierr = 0;

    /* When called before time loop, some values are not yet available. */
    if (init == true)
      return;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   icodcl_mesh_u[face_id] < 1
          || icodcl_mesh_u[face_id] > 4) {
        if (ale_bc_type[face_id] > 0) {
          ale_bc_type[face_id] = -ale_bc_type[face_id];
        }
        ierr++;
        bft_printf("Error face %d with bc_type = %d, ale_bc_type %d "
                   "and icodcl_mesh_u %d.\n",
                   face_id,
                   bc_type[face_id],
                   ale_bc_type[face_id],
                   icodcl_mesh_u[face_id]);
      }
      else {
#if 0
        bft_printf("Good result for face %d with bc_type = %d, ale_bc_type %d "
                   "and icodcl_mesh_u %d.\n",
                   face_id,
                   bc_type[face_id],
                   ale_bc_type[face_id],
                   icodcl_mesh_u[face_id]);
#endif
      }

      if (ale_bc_type[face_id] < 0) {
        irkerr    = cs_glob_rank_id;
        icoder[0] = -ale_bc_type[face_id];
        for (cs_lnum_t ii = 0; ii < 3; ii++) {
          icoder[1] = icodcl_mesh_u[face_id];
        }
      }
    }

    cs_parall_max(1, CS_INT_TYPE, &ierr);
    bool stop = false;

    if (ierr > 0) {
      stop = true;
      bft_printf(_("ALE method:\n\n"
                   "At least one boundary face has a boundary condition type\n"
                   "which is not recognized for mesh velocity:\n\n"
                   "  rank_id: %d,\n"
                   "  ale_bc_type: %d\n"
                   "  mesh_u->bc_coeffs->icodcl: %d,\n\n"
                   "The only allowed values for icodcl are\n"
                   "  1: Dirichlet\n"
                   "  2: Convective outlet\n"
                   "  3: Neumann\n"
                   "  4: Slip\n\n"
                   "Check boundary condition definitions.\n"),
                 irkerr, icoder[0], icoder[1]);
    }

    if (stop) {
      bft_printf(_("ALE method:\n\n"
                   "Inconsistency in boundary condition types"
                   " for mesh velocity:\n"
                   "  (cf. message(s) above)\n"));

      cs_boundary_conditions_error(ale_bc_type, nullptr);
    }

  } /* if (cs_glob_ale != CS_ALE_CDO) */

  /* Fluid velocity BCs for walls and symmetries
   * (due to mesh movement)
   * ------------------------------------------- */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (ale_bc_type[face_id] != CS_BOUNDARY_ALE_IMPOSED_VEL)
      continue;

    if (bc_type[face_id] == CS_SYMMETRY)
      for (int ii = 0; ii < 3; ii++)
        rcodcl1_vel[n_b_faces * ii + face_id] =
          rcodcl1_mesh_u[n_b_faces * ii + face_id];

    if (bc_type[face_id] == CS_SMOOTHWALL || bc_type[face_id] == CS_ROUGHWALL) {
      /* If nothing is set by the user, then the wall is supposed
       * to move with the sliding wall */
      if (rcodcl1_vel[n_b_faces * 0 + face_id] > cs_math_infinite_r * 0.5 &&
          rcodcl1_vel[n_b_faces * 1 + face_id] > cs_math_infinite_r * 0.5 &&
          rcodcl1_vel[n_b_faces * 2 + face_id] > cs_math_infinite_r * 0.5) {
        for (int ii = 0; ii < 3; ii++)
          rcodcl1_vel[n_b_faces * ii + face_id] =
            rcodcl1_mesh_u[n_b_faces * ii + face_id];
      }
      else {
        /* Otherwise, if the user has set the fluid velocity
         * Then the normal part is the one of the mesh velocity
         * and the tangential part is the one of the user (completed with 0
         * for non set values) */
        for (int ii = 0; ii < 3; ii++) {
          if (rcodcl1_vel[n_b_faces * ii + face_id] > cs_math_infinite_r * 0.5)
            rcodcl1_vel[n_b_faces * ii + face_id] = 0.;
        }

        const cs_nreal_t *rnxyz = b_face_u_normal[face_id];

        const cs_real_t rcodcxyz[3] = { rcodcl1_vel[n_b_faces * 0 + face_id],
                                        rcodcl1_vel[n_b_faces * 1 + face_id],
                                        rcodcl1_vel[n_b_faces * 2 + face_id] };
        cs_real_t rcodsn = 0.;
        for (int ii = 0; ii < 3; ii++)
          rcodsn += (rcodcl1_mesh_u[n_b_faces * ii + face_id] - rcodcxyz[ii]) *
                    rnxyz[ii];

        for (int ii = 0; ii < 3; ii++)
          rcodcl1_vel[n_b_faces * ii + face_id] =
            rcodcxyz[ii] + rcodsn * rnxyz[ii];
      }
    }
  }

  CS_FREE(_rcodcl1_mesh_u);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocation of ialtyb and impale for the ALE structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_allocate(void)
{
  if (cs_glob_ale <= CS_ALE_NONE)
    return;

  cs_mesh_t *m = cs_glob_mesh;
  CS_MALLOC(cs_glob_ale_data->impale, m->n_vertices, int);
  CS_MALLOC(cs_glob_ale_data->bc_type, m->n_b_faces, int);

  cs_arrays_set_value<int, 1>(m->n_b_faces,
                              0,
                              cs_glob_ale_data->bc_type);
  cs_arrays_set_value<int, 1>(m->n_vertices,
                              0,
                              cs_glob_ale_data->impale);

  CS_MALLOC_HD(cs_glob_ale_data->b_mass_flux_ale, m->n_b_faces,
               cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(cs_glob_ale_data->i_mass_flux_ale, m->n_i_faces,
               cs_real_t, cs_alloc_mode);

  cs_arrays_set_value<cs_real_t, 1>(m->n_b_faces,
                                    0.,
                                    cs_glob_ale_data->b_mass_flux_ale);

  cs_arrays_set_value<cs_real_t, 1>(m->n_i_faces,
                                    0.,
                                    cs_glob_ale_data->i_mass_flux_ale);

  cs_base_at_finalize(_ale_free);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell and face centers of gravity, cell volumes
 *         and update bad cells.
 *
 * \param[out]       min_vol        Minimum cell volume
 * \param[out]       max_vol        Maximum cell volume
 * \param[out]       tot_vol        Total cell volume
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_mesh_quantities(cs_real_t  *min_vol,
                              cs_real_t  *max_vol,
                              cs_real_t  *tot_vol)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_gradient_free_quantities();
  cs_cell_to_vertex_free();
  cs_mesh_quantities_compute(m, mq);
  cs_mesh_bad_cells_detect(m, mq);

  *min_vol = mq->min_vol;
  *max_vol = mq->max_vol;
  *tot_vol = mq->tot_vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project the displacement on mesh vertices (solved on cell center).
 *
 * \param[in]       ale_bc_type   Type of boundary for ALE
 * \param[in]       meshv         Mesh velocity
 * \param[in]       gradm         Mesh velocity gradient
 *                                (du_i/dx_j : gradv[][i][j])
 * \param[in]       claale        Boundary conditions A
 * \param[in]       clbale        Boundary conditions B
 * \param[in]       dt            Time step
 * \param[out]      disp_proj     Displacement projected on vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_project_displacement(const int           ale_bc_type[],
                            const cs_real_3_t  *meshv,
                            const cs_real_33_t  gradm[],
                            const cs_real_3_t  *claale,
                            const cs_real_33_t *clbale,
                            const cs_real_t    *dt,
                            cs_real_3_t        *disp_proj)
{
  bool *vtx_interior_indicator = nullptr;
  cs_real_t *vtx_counter = nullptr;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_vertices = m->n_vertices;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const int dim = m->dim;
  const cs_real_3_t *restrict vtx_coord  = (const cs_real_3_t *)m->vtx_coord;
  const cs_real_3_t *restrict cell_cen   = mq->cell_cen;
  const cs_real_3_t *restrict b_face_cog = mq->b_face_cog;

  CS_MALLOC(vtx_counter, n_vertices, cs_real_t);
  CS_MALLOC(vtx_interior_indicator, n_vertices, bool);

  cs_array_bool_fill_true(n_vertices, vtx_interior_indicator);
  cs_array_real_fill_zero(n_vertices, vtx_counter);
  cs_array_real_fill_zero(3 * n_vertices, (cs_real_t *)disp_proj);

  /* All nodes wich belong to a boundary face where the
     displacement is imposed (that is all faces except sliding BCs)
     are boundary nodes, the others are interior nodes. */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (ale_bc_type[face_id] != CS_BOUNDARY_ALE_SLIDING) {

      for (cs_lnum_t j = m->b_face_vtx_idx[face_id];
           j < m->b_face_vtx_idx[face_id+1];
           j++) {

        const cs_lnum_t  vtx_id = m->b_face_vtx_lst[j];
        vtx_interior_indicator[vtx_id] = false;

      } /* End of loop on vertices of the face */

    } /* Not sliding */

  } /* End of loop on border faces */

  /* Interior face and nodes treatment */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    const cs_lnum_t  cell_id1 = m->i_face_cells[face_id][0];
    const cs_lnum_t  cell_id2 = m->i_face_cells[face_id][1];
    const cs_real_t  dvol1 = 1. / mq->cell_vol[cell_id1];
    const cs_real_t  dvol2 = 1. / mq->cell_vol[cell_id2];
    const cs_real_t  dt_dvol1 = dt[cell_id1] * dvol1;
    const cs_real_t  dt_dvol2 = dt[cell_id2] * dvol2;

    if (cell_id1 < n_cells) { /* Test to take into account face only once */

      for (cs_lnum_t j = m->i_face_vtx_idx[face_id];
           j < m->i_face_vtx_idx[face_id+1];
           j++) {

        /* Get the vertex id */

        const cs_lnum_t  vtx_id = m->i_face_vtx_lst[j];

        if (vtx_interior_indicator[vtx_id]) {

          /* Get the vector from the cell center to the node */

          cs_real_3_t cen1_node, cen2_node;
          for (int i = 0; i < 3; i++) {
            cen1_node[i] = vtx_coord[vtx_id][i] - cell_cen[cell_id1][i];
            cen2_node[i] = vtx_coord[vtx_id][i] - cell_cen[cell_id2][i];
          }

          for (int i = 0; i < 3; i++) {
            disp_proj[vtx_id][i] +=
              dt_dvol1 *
                (meshv[cell_id1][i] +
                 cs_math_3_dot_product(gradm[cell_id1][i], cen1_node)) +
              dt_dvol2 * (meshv[cell_id2][i] +
                          cs_math_3_dot_product(gradm[cell_id2][i], cen2_node));
          }

          vtx_counter[vtx_id] += dvol1 + dvol2;

        } /* End of Interior nodes */

      }

    }

  } /* End of loop on internal faces */

  /* Border face treatment.
     only border face contribution */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    const cs_lnum_t  cell_id = m->b_face_cells[face_id];

    for (cs_lnum_t j = m->b_face_vtx_idx[face_id];
         j < m->b_face_vtx_idx[face_id+1]
         && ale_bc_type[face_id] != CS_BOUNDARY_ALE_SLIDING; j++) {

      const cs_lnum_t  vtx_id = m->b_face_vtx_lst[j];

      if (!vtx_interior_indicator[vtx_id]) {

        /* Get the vector from the face center to the node*/

        cs_real_3_t face_node;
        for (int i = 0; i < 3; i++)
          face_node[i] = vtx_coord[vtx_id][i] - b_face_cog[face_id][i];

        /* 1st order extrapolation of the mesh velocity at the face center
         * to the node */

        cs_real_3_t vel_node;
        for (int i = 0; i < 3; i++)
          vel_node[i] = claale[face_id][i] +
                        cs_math_3_dot_product(gradm[cell_id][i], face_node);

        const cs_real_t dsurf = 1./mq->b_face_surf[face_id];

        for (int i = 0; i < 3; i++)
          disp_proj[vtx_id][i] +=
            dsurf * dt[cell_id] *
            (vel_node[i] +
             cs_math_3_dot_product(clbale[face_id][i], meshv[cell_id]));

        vtx_counter[vtx_id] += dsurf;

      } /* End of boundary nodes */

    } /* End of loop on vertices of the face */

  } /* End of loop on border faces */

  if (m->vtx_interfaces != nullptr) {

    cs_interface_set_sum(m->vtx_interfaces,
                         n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         disp_proj);

    cs_interface_set_sum(m->vtx_interfaces,
                         n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         vtx_counter);
  }

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    /* Might be null for a vertex belonging to sliding faces only */
    if (vtx_counter[v_id] > 0.) {
      for (int i = 0; i < dim; i++)
        disp_proj[v_id][i] /= vtx_counter[v_id];
    }
  }

  /* If the boundary face IS a sliding face.
     We project the displacment paralelly to the face. */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_SLIDING) {

      for (cs_lnum_t j = m->b_face_vtx_idx[face_id];
           j < m->b_face_vtx_idx[face_id+1]; j++) {

        const cs_lnum_t  vtx_id = m->b_face_vtx_lst[j];

        disp_proj[vtx_id][0] =
          cs_math_3_dot_product(clbale[face_id][0], disp_proj[vtx_id]);
        disp_proj[vtx_id][1] =
          cs_math_3_dot_product(clbale[face_id][1], disp_proj[vtx_id]);
        disp_proj[vtx_id][2] =
          cs_math_3_dot_product(clbale[face_id][2], disp_proj[vtx_id]);

      } /* End of loop on vertices of the face */

    } /* Sliding condition */

  } /* End of loop on border faces */

  CS_FREE(vtx_counter);
  CS_FREE(vtx_interior_indicator);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update mesh in the ALE framework.
 *
 * \param[in]       itrale        number of the current ALE iteration
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_mesh(int  itrale)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const int  ndim = m->dim;
  const cs_lnum_t  n_vertices = m->n_vertices;

  cs_real_3_t *vtx_coord = (cs_real_3_t *)m->vtx_coord;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  cs_time_step_t *ts = cs_get_glob_time_step();

  /* Initialization */

  const cs_equation_param_t *eqp =
    cs_field_get_equation_param_const(CS_F_(mesh_u));

  if (eqp->verbosity >= 1)
    bft_printf("\n ---------------------------------------------------"
               "---------\n\n"
               "  Update mesh (ALE)\n"
               "  =================\n\n");

  /* Retrieving fields */

  cs_field_t  *f_displ = cs_field_by_name("mesh_displacement");
  cs_real_3_t *disale = (cs_real_3_t *)(f_displ->val);
  cs_real_3_t *disala = (cs_real_3_t *)(f_displ->val_pre);
  cs_real_3_t *xyzno0 = (cs_real_3_t *)(cs_field_by_name("vtx_coord0")->val);

  /* Update geometry */

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    for (cs_lnum_t idim = 0; idim < ndim; idim++) {
      vtx_coord[v_id][idim] = xyzno0[v_id][idim] + disale[v_id][idim];
      disala[v_id][idim]    = disale[v_id][idim];
    }
  }

  cs_ale_update_mesh_quantities(&(mq->min_vol), &(mq->max_vol), &(mq->tot_vol));

  /* Abort at the end of the current time-step if there is a negative volume */

  if (mq->min_vol <= 0.)
    ts->nt_max = ts->nt_cur;

  /* The mesh velocity is reverted to its initial value if the current time step
     is the initialization time step */

  if (itrale == 0) {

    cs_field_t *f = cs_field_by_name("mesh_velocity");
    const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(f->location_id)[2];

    if (ndim == 3)
      cs_array_copy(3*n_elts, f->val_pre, f->val);
    else {
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++)
        for (int idim = 0; idim < ndim; idim++)
          f->val[3*e_id + idim] = f->val_pre[3*e_id + idim];
    }

  } /* itrale == 0 */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update ALE BCs for required for the fluid
 *
 * \param[out]      ale_bc_type   type of ALE bcs
 * \param[out]      b_fluid_vel   Fluid velocity at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_bcs(int         *ale_bc_type,
                  cs_real_3_t *b_fluid_vel)
{
  assert(cs_glob_ale == CS_ALE_CDO);

  /* First update ALE BCs expect free surface that will be updated after */
  _update_bcs(cs_glob_domain, ale_bc_type, b_fluid_vel);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a Poisson equation on the mesh velocity in ALE framework.
 *
 * It also updates the mesh displacement
 * so that it can be used to update mass fluxes (due to mesh displacement).
 *
 * \param[in]  iterns  Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_solve_mesh_velocity(const int        iterns,
                           const cs_real_t  b_rho[],
                           const cs_real_t  i_mass_flux[],
                           const cs_real_t  b_mass_flux[])
{
  const int  *impale = cs_glob_ale_data->impale;
  const int  *ale_bc_type = cs_glob_ale_data->bc_type;

  if (cs_glob_ale == CS_ALE_LEGACY)
    _ale_solve_poisson_legacy(cs_glob_domain, iterns, impale, ale_bc_type,
                              b_rho, i_mass_flux, b_mass_flux);

  else if (cs_glob_ale == CS_ALE_CDO)
    _ale_solve_poisson_cdo(cs_glob_domain, impale);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the mesh velocity solving with CDO
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_activate(void)
{
  if (cs_ale_active)
    return;

  cs_ale_active = true;

  cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_WITH_FV);

  cs_equation_t  *eq =
    cs_equation_add("mesh_velocity", /* equation name */
                    "mesh_velocity", /* associated variable field name */
                    CS_EQUATION_TYPE_PREDEFINED,
                    3,                        /* dimension of the unknown */
                    CS_BC_SYMMETRY); /* default boundary */

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  /* Predefined settings */

  cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");

  /* BC settings */

  cs_equation_param_set(eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");

  /* System to solve is SPD by construction */

  cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "cg");
  cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "jacobi");
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RESNORM_TYPE, "filtered");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if mesh velocity solving with CDO is activated
 *
 * \return true if mesh velocity solving with CDO is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_ale_is_activated(void)
{
  if (cs_ale_active)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add "property" fields dedicated to the ALE model.
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_add_property_fields(void)
{
  assert(cs_glob_ale != CS_ALE_NONE);

  int  log_key_id = cs_field_key_id("log");
  int  post_key_id = cs_field_key_id("post_vis");

  /* Add the mesh displacement
     ------------------------- */

  cs_field_t *f_md = cs_field_create("mesh_displacement",
                                     CS_FIELD_PROPERTY,
                                     CS_MESH_LOCATION_VERTICES,
                                     3,      /* dim */
                                     true);  /* has previous */

  /* Add a label for this field */

  cs_field_set_key_str(f_md, cs_field_key_id("label"), "Mesh displacement");

  /* Output related to this field */

  const int post_flag = CS_POST_ON_LOCATION;

  cs_field_set_key_int(f_md, log_key_id, 1);
  cs_field_set_key_int(f_md, post_key_id, post_flag);

  /* Handle the restart settings */

  cs_field_set_key_int(f_md,
                       cs_field_key_id("restart_file"),
                       CS_RESTART_AUXILIARY);
  cs_field_set_key_int(f_md,
                       cs_field_key_id("restart_n_values"),
                       2);

  /* Add the initial vertex coordinates
     ---------------------------------- */

  cs_field_t *f_xyz0 = cs_field_create("vtx_coord0",
                                       CS_FIELD_PROPERTY,
                                       CS_MESH_LOCATION_VERTICES,
                                       3,       /* dim */
                                       false);  /* has previous */

  /* No output related to this field */

  cs_field_set_key_int(f_xyz0, log_key_id, 0);
  cs_field_set_key_int(f_xyz0, post_key_id, 0);

  /* Add the cell volume as a owner field to have vol^n and vol^{n+1}
     ----------------------------------------------------------------  */

  int field_type = CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY;
  cs_field_t *f_cellvol = cs_field_create("cell_vol",
                                          field_type,
                                          CS_MESH_LOCATION_CELLS,
                                          1,       /* dim */
                                          true);   /* has previous */

   /* Add output related to this field */

  cs_field_set_key_int(f_cellvol, log_key_id, 0);
  cs_field_set_key_int(f_cellvol, post_key_id, 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize volume fields dedicated to the ALE model.
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_initialize_volume_fields(void)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_real_t *mq_cell_vol = mq->cell_vol;

  cs_field_t *f_cellvol = cs_field_by_name("cell_vol");
  cs_real_t *cell_vol = f_cellvol->val;
  cs_real_t *cell_vol_pre = f_cellvol->val_pre;

  /* Initialize f_cell_vol by the initial mq->cell_vol */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cell_vol[c_id] = mq_cell_vol[c_id];
    cell_vol_pre[c_id] = mq_cell_vol[c_id];
  }

  cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, cell_vol);
  cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, cell_vol_pre);

  /* In neptune, volnp1 is computed in cs_ale_compute_volume_from_displacement
     and is stored in cs_field_by_name("cell_vol"). The most used Vol^n is
     stored in cell_vol_pre.
     In saturne, we always use cell_vol_pre wich is updated in cs_ale_mesh_update
     and cs_field_by_name("cell_vol")->val is never used. As the f->type is not
     CS_FIELD_VARIABLE, we have not the current_to previous val_pre = val.
  */

  mq->cell_vol = cell_vol_pre;

  if (cs_glob_ale == CS_ALE_CDO) {
    cs_glob_domain->cdo_quantities->cell_vol = cell_vol_pre;
  }

  /* Free _cell_vol which is not owner anymore */
  CS_FREE(mq->_cell_vol);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup the equations solving the mesh velocity when CDO is activated
 *
 * \param[in, out] domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_init_setup(cs_domain_t *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /* Mesh viscosity (iso or ortho)
   * TODO declare it before: add in activate, def here...  */

  int dim = cs_field_by_name("mesh_viscosity")->dim;

  cs_property_t  *mesh_visc = cs_property_by_name("mesh_viscosity");

  if (mesh_visc == nullptr)  {      /* Not already added */

    cs_property_type_t  type = 0;

    if (dim == 1)
      type = CS_PROPERTY_ISO;
    else if (dim == 3)
      type = CS_PROPERTY_ORTHO;
    else if (dim == 6)
      type = CS_PROPERTY_ANISO_SYM;
    else if (dim == 9)
      type = CS_PROPERTY_ANISO;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid dimension (=%d) for the mesh viscosity.",
                __func__, dim);

    /* Add and define this property */

    mesh_visc = cs_property_add("mesh_viscosity", type);
    cs_property_def_by_field(mesh_visc, cs_field_by_name("mesh_viscosity"));

  }

  /* Update the equation related to the mesh displacement */

  cs_equation_t  *eq = cs_equation_by_name("mesh_velocity");
  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  assert(mesh_visc != nullptr);
  cs_equation_add_diffusion(eqp, mesh_visc);

  /* Add the variable field */

  cs_equation_predefined_create_field(1, eq); /* Always has_previous */
}

/*----------------------------------------------------------------------------
 *!
 * \brief Print the ALE options to setup.log.
 *
 *----------------------------------------------------------------------------*/

void
cs_ale_log_setup(void)
{
  if (cs_glob_ale < CS_ALE_NONE || cs_glob_ale > CS_ALE_CDO)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: invalid value for cs_glob_ale (%d)."),
              __func__, (int)cs_glob_ale);

  const char *ale_status_str[]
    = {N_("Inactive (CS_ALE_NONE)"),
       N_("Active (CS_ALE_LEGACY)"),
       N_("Active (CS_ALE_CDO)")};

  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "ALE method (moving mesh)\n"
                 "------------------------\n\n"
                 "  %s\n"),
                _(ale_status_str[cs_glob_ale]));

  if (cs_glob_ale == CS_ALE_NONE) {
    return;
  }

  cs_log_printf(CS_LOG_SETUP,
                _("   cs_glob_ale_n_ini_f: %d "
                  "(iterations for flow initialization)\n"),
                cs_glob_ale_n_ini_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup the equations solving the mesh velocity
 *
 * \param[in]   domain       pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_setup_boundaries(const cs_domain_t   *domain)
{
  const cs_mesh_t  *mesh = domain->mesh;
  const cs_lnum_t  n_vertices = mesh->n_vertices;

  cs_equation_param_t *eqp = cs_equation_param_by_name("mesh_velocity");

  if (_cdo_bc == nullptr) {

    CS_MALLOC(_cdo_bc, 1, cs_ale_cdo_bc_t);
    CS_MALLOC(_cdo_bc->vtx_values, 3*n_vertices, cs_real_t);

    cs_array_real_fill_zero(3*n_vertices, _cdo_bc->vtx_values);

    _cdo_bc->n_selections = 0;
    _cdo_bc->n_vertices = nullptr;
    _cdo_bc->vtx_select = nullptr;

  }

  bool   *vtag = nullptr;
  CS_MALLOC(vtag, n_vertices, bool);

  for (int b_id = 0;  b_id < domain->ale_boundaries->n_boundaries; b_id++) {

    const int z_id = domain->ale_boundaries->zone_ids[b_id];
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);

    switch(domain->ale_boundaries->types[b_id]) {

    case CS_BOUNDARY_ALE_FIXED:
      {
        cs_real_t  bc_value[3] = {0., 0., 0.};

        cs_equation_add_bc_by_value(eqp,
                                    CS_BC_HMG_DIRICHLET,
                                    z->name,
                                    bc_value);
      }
      break;

    case CS_BOUNDARY_ALE_IMPOSED_VEL:
      cs_equation_add_bc_by_array(eqp,
                                  CS_BC_DIRICHLET,
                                  z->name,
                                  cs_flag_primal_vtx,
                                  _cdo_bc->vtx_values,
                                  false, /* Do not transfer ownership */
                                  true); /* full length */

      /* Add a new list of selected vertices. */
      _update_bc_list(mesh, z, vtag);
      break;

    case CS_BOUNDARY_ALE_IMPOSED_DISP:
    case CS_BOUNDARY_ALE_INTERNAL_COUPLING:
    case CS_BOUNDARY_ALE_EXTERNAL_COUPLING:
      cs_equation_add_bc_by_array(eqp,
                                  CS_BC_DIRICHLET,
                                  z->name,
                                  cs_flag_primal_vtx,
                                  _cdo_bc->vtx_values,
                                  false, /* Do not transfer ownership */
                                  true); /* full length */

      /* Add a new list of selected vertices. */
      _update_bc_list(mesh, z, vtag);
      break;

    case CS_BOUNDARY_ALE_FREE_SURFACE:
      cs_equation_add_bc_by_array(eqp,
                                  CS_BC_DIRICHLET,
                                  z->name,
                                  cs_flag_primal_vtx,
                                  _cdo_bc->vtx_values,
                                  false, /* Do not transfer ownership */
                                  true); /* full length */

      /* Add a new list of selected vertices. */
      _update_bc_list(mesh, z, vtag);
      break;

    case CS_BOUNDARY_ALE_SLIDING:
      cs_equation_add_sliding_condition(eqp, z->name);
      break;

    default:
      bft_error(__FILE__,
                __LINE__,
                0,
                _(" %s: Boundary type %d for ALE not allowed for %s."),
                __func__,
                domain->ale_boundaries->types[b_id],
                z->name);
    }

  }

  /* Free memory */
  CS_FREE(vtag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for the equation of the mesh velocity
 *
 * \param[in, out] domain       pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_finalize_setup(cs_domain_t     *domain)
{
  const cs_mesh_t  *m = domain->mesh;

  if (_vtx_coord0 == nullptr) {

    CS_MALLOC(_vtx_coord0, m->n_vertices, cs_real_3_t);

    memcpy(_vtx_coord0, m->vtx_coord, m->n_vertices * sizeof(cs_real_3_t));

  }

  /* Set the boundary conditions */
  cs_ale_setup_boundaries(domain);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the ALE mesh velocity solving
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_destroy_all(void)
{
  CS_FREE(_vtx_coord0);

  if (_cdo_bc != nullptr) {
    CS_FREE(_cdo_bc->vtx_values);

    for (int i = 0; i < _cdo_bc->n_selections; i++)
      CS_FREE(_cdo_bc->vtx_select[i]);
    CS_FREE(_cdo_bc->vtx_select);
    CS_FREE(_cdo_bc->n_vertices);

    CS_FREE(_cdo_bc);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read ALE data from restart file.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_restart_read(cs_restart_t  *r)
{
  if (cs_glob_ale == CS_ALE_NONE)
    return;

  cs_field_t *f_displ = cs_field_by_name("mesh_displacement");

  int retcode = cs_restart_read_field_vals(r, f_displ->id, 0);

  if (retcode != CS_RESTART_SUCCESS) {
    cs_real_3_t *displ = (cs_real_3_t *)f_displ->val;
    retcode = cs_restart_read_real_3_t_compat(r,
                                              "vertex_displacement",
                                              "deplact_x_no",
                                              "deplact_y_no",
                                              "deplact_z_no",
                                              CS_MESH_LOCATION_VERTICES,
                                              displ);
  }

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_field_vals(r, f_displ->id, 1);
  else
    bft_error
      (__FILE__, __LINE__, 0,
       "%s: reading mesh vertices displacement in %s.",
       __func__, cs_restart_get_name(r));

  if (retcode != CS_RESTART_SUCCESS)
    cs_field_current_to_previous(f_displ);

  /* Geometric parameters must be recalculated */

  const cs_lnum_t n_vtx = cs_glob_mesh->n_vertices;

  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  cs_real_3_t *disale = (cs_real_3_t *)f_displ->val;
  cs_real_3_t *coordp = (cs_real_3_t *)cs_glob_mesh->vtx_coord;
  cs_real_3_t *coord0 = (cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

  cs_dispatch_context ctx;
  ctx.parallel_for(n_vtx, [=] CS_F_HOST_DEVICE (cs_lnum_t v_id) {
    for (cs_lnum_t idim = 0; idim < 3; idim++) {
      coordp[v_id][idim] = coord0[v_id][idim] + disale[v_id][idim];
    }
  });
  ctx.wait();

  cs_ale_update_mesh_quantities(&(mq->min_vol), &(mq->max_vol), &(mq->tot_vol));

  /* Abort at the end of the current time-step if there is a negative volume */
  if (mq->min_vol <= 0)
    cs_time_step_define_nt_max(cs_glob_time_step->nt_cur);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write ALE data from restart file.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_restart_write(cs_restart_t  *r)
{
  if (cs_glob_ale == CS_ALE_NONE)
    return;

  cs_field_t *f_displ = cs_field_by_name("mesh_displacement");

  cs_restart_write_field_vals(r, f_displ->id, 0);
  cs_restart_write_field_vals(r, f_displ->id, 1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
