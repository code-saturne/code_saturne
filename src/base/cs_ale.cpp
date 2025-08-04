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
#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "base/cs_boundary_zone.h"
#include "alge/cs_cell_to_vertex.h"
#include "cdo/cs_cdo_quantities.h"
#include "cdo/cs_cdo_connect.h"
#include "cdo/cs_cdo_main.h"
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
                                           .bc_type = nullptr};

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
     CS_BOUNDARY_ALE_IMPOSED_VEL and CS_BOUNDARY_ALE_IMPOSED_DISP (definitions
     by value or by a sliding condition are excluded. The position of the array
     in the list is given implicitly and is related to the order in which the
     boundary conditions are defined. (cf. \ref cs_ale_setup_boundaries
     definition for more details). The aim of this list of arrays is to speed-up
     the settings of the boundary conditions by avoiding doing several times the
     same enforcement. */

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
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *)mq->b_face_normal;
  const cs_real_3_t *restrict b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

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
      (cs_math_3_dot_product(normal, b_face_normal[face_id])
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
      cs_real_t dz = CS_ABS(cs_math_3_dot_product(normal, e_cog));
      cs_real_3_t e_cog_hor;
      cs_math_3_orthogonal_projection(normal, e_cog, e_cog_hor);
      cs_real_t dx = cs_math_3_norm(e_cog_hor);
      /* Too high slope */
      if (dz > dx)
        f_need_filter = 1;

      f_has_max = CS_MAX(f_has_max, _is_loc_max[v_id0]);
      f_has_min = CS_MAX(f_has_min, _is_loc_min[v_id0]);
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
  const cs_mesh_t  *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_real_3_t *b_face_normal = (const cs_real_3_t *)mq->b_face_normal;
  const cs_real_3_t *restrict vtx_coord  = (const cs_real_3_t *)m->vtx_coord;
  const cs_real_3_t *restrict b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

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

          cs_real_3_t normal;
          cs_math_3_normalize(b_face_normal[face_id], normal);

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
      {
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

          cs_real_3_t normal;
          /* Normal direction is given by the gravity */
          cs_math_3_normalize(b_face_normal[face_id], normal);
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
 * \param[in]       domain        domain quantities
 * \param[in]       iterns        Navier-Stokes iteration number
 * \param[in]       impale        Indicator for fixed node displacement
 * \param[in]       ale_bc_type   Type of boundary for ALE
 */
/*----------------------------------------------------------------------------*/

static void
_ale_solve_poisson_legacy(const cs_domain_t *domain,
                          const int          iterns,
                          const int         *impale,
                          const int         *ale_bc_type)
{
  const cs_mesh_t *m = domain->mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = (const cs_lnum_t *)m->b_face_cells;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;
  const cs_real_t *b_dist = (const cs_real_t *)mq->b_dist;
  const cs_real_3_t *b_face_normal = (const cs_real_3_t *)mq->b_face_normal;
  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  /* The mass flux is necessary to call cs_equation_iterative_solve_vector
     but not used (iconv = 0), except for the free surface, where it is used
     as a boundary condition */

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *i_massflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf))->val;
  const cs_real_t *b_massflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

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

  /* Density at the boundary */
  cs_real_t *brom = CS_F_(rho_b)->val;

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
      } else if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
        for (int isou = 0; isou < 6; isou++)
          hintt[isou] = CS_F_(vism)->val[6*cell_id+isou] / distbf;
      }

      cs_real_t prosrf = cs_math_3_dot_product(grav, b_face_normal[face_id]);

      cs_real_3_t pimpv;
      for (int i = 0; i < 3; i++)
        pimpv[i] = grav[i]*b_massflux[face_id]/(brom[face_id]*prosrf);

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
                                     i_massflux,
                                     b_massflux,
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
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

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

  CS_MALLOC(cs_glob_ale_data->impale, cs_glob_mesh->n_vertices, int);
  CS_MALLOC(cs_glob_ale_data->bc_type, cs_glob_mesh->n_b_faces, int);

  cs_arrays_set_value<int, 1>(cs_glob_mesh->n_b_faces,
                              0,
                              cs_glob_ale_data->bc_type);
  cs_arrays_set_value<int, 1>(cs_glob_mesh->n_vertices,
                              0,
                              cs_glob_ale_data->impale);

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
  const cs_real_3_t *restrict cell_cen   = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_3_t *restrict b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

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
cs_ale_solve_mesh_velocity(int  iterns)
{
  const int  *impale = cs_glob_ale_data->impale;
  const int  *ale_bc_type = cs_glob_ale_data->bc_type;

  if (cs_glob_ale == CS_ALE_LEGACY)
    _ale_solve_poisson_legacy(cs_glob_domain, iterns, impale, ale_bc_type);

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

  cs_equation_t *eq =
    cs_equation_add("mesh_velocity", /* equation name */
                    "mesh_velocity", /* associated variable field name */
                    CS_EQUATION_TYPE_PREDEFINED,
                    3,                  /* dimension of the unknown */
                    CS_BC_HMG_NEUMANN); /* default boundary */

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
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Boundary for ALE not allowed  %s."),
                __func__, z->name);
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

  #pragma omp parallel for if (n_vtx > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_vtx; v_id++) {
    for (cs_lnum_t idim = 0; idim < 3; idim++) {
      coordp[v_id][idim] = coord0[v_id][idim] + disale[v_id][idim];
    }
  }

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
