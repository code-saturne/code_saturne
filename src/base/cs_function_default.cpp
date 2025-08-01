/*============================================================================
 * Base predefined function objects.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_assert.h"
#include "alge/cs_balance_by_zone.h"
#include "base/cs_boundary_zone.h"
#include "elec/cs_elec_model.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_function.h"
#include "base/cs_function_default.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_math.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "cdo/cs_property.h"
#include "base/cs_post.h"
#include "base/cs_rotation.h"
#include "base/cs_thermal_model.h"
#include "base/cs_turbomachinery.h"
#include "base/cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_function_default.cpp
        Base predefined function objects.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the associated rank based on a mesh location's range set.
 *
 * \param[in]       rs           pointer to range set structure, or null
 * \param[in]       location_id  base associated mesh location idà
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_range_set_mpi_rank_id(const cs_range_set_t  *rs,
                       int                    location_id,
                       cs_lnum_t              n_elts,
                       const cs_lnum_t       *elt_ids,
                       void                  *vals)
{
  const cs_lnum_t n_loc_elts = cs_mesh_location_get_n_elts(location_id)[0];

  auto vals_i = static_cast<int *>(vals);

  int *e_rank_id;

  if (n_elts != n_loc_elts || elt_ids != nullptr)
    CS_MALLOC(e_rank_id, n_loc_elts, int);
  else
    e_rank_id = vals_i;

  if (rs != nullptr) {
    for (cs_lnum_t i = 0; i < rs->n_elts[0]; i++)
      e_rank_id[i] = cs_glob_rank_id;
    for (cs_lnum_t i = rs->n_elts[0]; i < n_loc_elts; i++)
      e_rank_id[i] = 0;

    cs_range_set_scatter(rs,
                         CS_INT_TYPE,
                         1,
                         e_rank_id,
                         e_rank_id);

    if (rs->ifs != nullptr)
      cs_interface_set_max(rs->ifs,
                           n_loc_elts,
                           1,
                           true,  /* interlace */
                           CS_INT_TYPE,
                           e_rank_id);
  }
  else {
    for (cs_lnum_t i = 0; i <  n_loc_elts; i++)
      e_rank_id[i] = cs_glob_rank_id;
  }

  if (e_rank_id != vals_i) {
    int *_vals = vals_i;
    for (cs_lnum_t i = 0; i < n_elts; i++)
      _vals[i] = e_rank_id[elt_ids[i]];

    CS_FREE(e_rank_id);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the associated rank at a given mesh location.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_location_mpi_rank_id(int               location_id,
                      cs_lnum_t         n_elts,
                      const cs_lnum_t  *elt_ids,
                      void             *input,
                      void             *vals)
{
  switch(location_id) {
    case CS_MESH_LOCATION_INTERIOR_FACES: {
      cs_mesh_t *m = static_cast<cs_mesh_t *>(input);

      cs_gnum_t *g_i_face_num = m->global_i_face_num;
      if (g_i_face_num == nullptr) {
        cs_lnum_t n_i_faces = m->n_i_faces;
        CS_MALLOC(g_i_face_num, n_i_faces, cs_gnum_t);
        for (cs_lnum_t i = 0; i < n_i_faces; i++)
          g_i_face_num[i] = (cs_gnum_t)i + 1;
      }

      cs_interface_set_t *face_interfaces
        = cs_interface_set_create(m->n_i_faces,
                                  nullptr,
                                  g_i_face_num,
                                  m->periodicity,
                                  0,
                                  nullptr,
                                  nullptr,
                                  nullptr);

      if (m->global_i_face_num != g_i_face_num)
        CS_FREE(g_i_face_num);

      cs_range_set_t *rs = cs_range_set_create(face_interfaces,
                                               nullptr,
                                               m->n_i_faces,
                                               false, /* balance */
                                               2,  /* tr_ignore */
                                               0); /* g_id_base */

      _range_set_mpi_rank_id(rs, location_id, n_elts, elt_ids, vals);

      cs_range_set_destroy(&rs);
      cs_interface_set_destroy(&face_interfaces);
    } break;

    case CS_MESH_LOCATION_VERTICES: {
      cs_mesh_t      *m  = static_cast<cs_mesh_t *>(input);
      cs_range_set_t *rs = m->vtx_range_set;

      if (rs == nullptr)
        rs = cs_range_set_create(m->vtx_interfaces,
                                 nullptr,
                                 m->n_vertices,
                                 false, /* balance */
                                 2,  /* tr_ignore */
                                 0); /* g_id_base */

      _range_set_mpi_rank_id(rs, location_id, n_elts, elt_ids, vals);

      if (rs != m->vtx_range_set)
        cs_range_set_destroy(&rs);
    } break;

  default:
    {
    int *_vals = static_cast<int *>(vals);
    for (cs_lnum_t i = 0; i < n_elts; i++)
      _vals[i] = cs_glob_rank_id;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the refinement level at a given mesh location.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *)
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_location_r_gen(int               location_id,
                cs_lnum_t         n_elts,
                const cs_lnum_t  *elt_ids,
                void             *input,
                void             *vals)
{
  const cs_mesh_t *m       = static_cast<const cs_mesh_t *>(input);
  const char *e_r_gen = nullptr;

  char *c_r_gen = nullptr;
  int  *r_gen   = static_cast<int *>(vals);

  /* Locations with stored refinement generations */

  if (location_id == CS_MESH_LOCATION_INTERIOR_FACES)
    e_r_gen = m->i_face_r_gen;
  else if (location_id == CS_MESH_LOCATION_VERTICES)
    e_r_gen = m->vtx_r_gen;

  /* other base locations */

  if (   location_id == CS_MESH_LOCATION_CELLS
      || location_id == CS_MESH_LOCATION_BOUNDARY_FACES) {

    CS_MALLOC(c_r_gen, m->n_cells_with_ghosts, char);

    for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++)
      c_r_gen[i] = 0;

    /* Note: when mesh/face adjacencies are available,
       a cell-based loop would allow easier threading. */

    const cs_lnum_2_t *restrict i_face_cells =
      (const cs_lnum_2_t *)m->i_face_cells;

    for (cs_lnum_t i = 0; i < m->n_i_faces; i++) {
      char f_r_gen = m->i_face_r_gen[i];
      for (cs_lnum_t j = 0; j < 2; j++) {
        cs_lnum_t c_id = i_face_cells[i][j];
        if (c_r_gen[c_id] < f_r_gen)
          c_r_gen[c_id] = f_r_gen;
      }
    }

    if (location_id == CS_MESH_LOCATION_CELLS)
      e_r_gen = c_r_gen;

    else if (location_id == CS_MESH_LOCATION_BOUNDARY_FACES) {
      if (elt_ids != nullptr) {
        for (cs_lnum_t idx = 0; idx < n_elts; idx++) {
          cs_lnum_t i = elt_ids[idx];
          cs_lnum_t c_id = m->b_face_cells[i];
          r_gen[idx] = c_r_gen[c_id];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t c_id = m->b_face_cells[i];
          r_gen[i] = c_r_gen[c_id];
        }
      }
      CS_FREE(c_r_gen);
      return;

    }

  }

  /* Now apply stored or generated refinement generation
     (note that for boundary faces, processing has already been done
     and we have returned prior to reaching this point) */

  if (e_r_gen != nullptr) {
    if (elt_ids != nullptr) {
      for (cs_lnum_t idx = 0; idx < n_elts; idx++) {
        cs_lnum_t i = elt_ids[idx];
        r_gen[idx] = e_r_gen[i];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        r_gen[i] = e_r_gen[i];
      }
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++)
      r_gen[i] = 0;
  }

  CS_FREE(c_r_gen);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the absolute pressure associated to the given cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_absolute_pressure_f(int               location_id,
                     cs_lnum_t         n_elts,
                     const cs_lnum_t  *elt_ids,
                     void             *input,
                     void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_rotation_t *r = cs_glob_rotation;
  const int _c_r_num[1] = {0};
  const int *c_r_num = _c_r_num;
  cs_lnum_t c_r_step = 0;

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE) {
    c_r_num = cs_turbomachinery_get_cell_rotor_num();
    c_r_step = 1;
  }

  cs_real_t *p_abs = static_cast<cs_real_t *>(vals);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_t *cvar_pr = cs_field_by_name("pressure")->val;
  const cs_real_t *cpro_rho = cs_field_by_name("density")->val;

  if (elt_ids != nullptr) {
    for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
      cs_lnum_t i = elt_ids[idx];
      int r_num = c_r_num[i * c_r_step];
      cs_real_t vr[3];
      cs_rotation_velocity(r + r_num, cell_cen[i], vr);
      p_abs[idx] = cvar_pr[i] + cpro_rho[i] * 0.5 * cs_math_3_square_norm(vr);
    }
  }

  else {
    for (cs_lnum_t i = 0; i <  n_elts; i++) {
      int r_num = c_r_num[i * c_r_step];
      cs_real_t vr[3];
      cs_rotation_velocity(r + r_num, cell_cen[i], vr);
      p_abs[i] = cvar_pr[i] + cpro_rho[i] * 0.5 * cs_math_3_square_norm(vr);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the absolute velocity associated to the given cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_absolute_velocity_f(int               location_id,
                     cs_lnum_t         n_elts,
                     const cs_lnum_t  *elt_ids,
                     void             *input,
                     void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_rotation_t *r = cs_glob_rotation;
  const int _c_r_num[1] = {0};
  const int *c_r_num = _c_r_num;
  cs_lnum_t c_r_step = 0;

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE) {
    c_r_num = cs_turbomachinery_get_cell_rotor_num();
    c_r_step = 1;
  }

  cs_real_3_t *v_abs = static_cast<cs_real_3_t *>(vals);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_3_t *cvar_vel
    = (const cs_real_3_t *)cs_field_by_name("velocity")->val;

  if (elt_ids != nullptr) {
    for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
      cs_lnum_t i = elt_ids[idx];
      int r_num = c_r_num[i * c_r_step];
      cs_real_t vr[3];
      cs_rotation_velocity(r + r_num, cell_cen[i], vr);
      for (cs_lnum_t j = 0; j < 3; j++)
        v_abs[idx][j] = cvar_vel[i][j] + vr[j];
    }
  }

  else {
    for (cs_lnum_t i = 0; i <  n_elts; i++) {
      int r_num = c_r_num[i * c_r_step];
      cs_real_t vr[3];
      cs_rotation_velocity(r + r_num, cell_cen[i], vr);
      for (cs_lnum_t j = 0; j < 3; j++)
        v_abs[i][j] = cvar_vel[i][j] + vr[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access function objects specific to
 *        turbomachinery models (absolute_pressure, absolute_velocity).
 *
 * \return  pointer to the associated function object if turbomachinery model
 *          is active, or null
 */
/*----------------------------------------------------------------------------*/

static void
_define_coriolis_functions(void)
{
  assert(cs_glob_physical_constants->icorio > 0);

  /* Absolute pressure */

  {
    cs_function_t *f
      = cs_function_define_by_func("absolute_pressure",
                                   CS_MESH_LOCATION_CELLS,
                                   1,
                                   true,
                                   CS_REAL_TYPE,
                                   _absolute_pressure_f,
                                   nullptr);

    const char label[] = "Abs Pressure";
    CS_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }

  /* Absolute velocity */

  {
    cs_function_t *f
      = cs_function_define_by_func("absolute_velocity",
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true,
                                   CS_REAL_TYPE,
                                   _absolute_velocity_f,
                                   nullptr);

    const char label[] = "Abs Velocity";
    CS_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ensure the "boundary_stress" field is present.
 */
/*----------------------------------------------------------------------------*/

static void
_ensure_boundary_stress_is_present(void)
{
  const char name[] = "boundary_stress";
  cs_field_t *bf = cs_field_by_name_try(name);
  if (bf == nullptr) {
    int type = CS_FIELD_INTENSIVE | CS_FIELD_POSTPROCESS;
    int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;

    bf = cs_field_create(name, type, location_id, 3, false);
    cs_field_set_key_int(bf, cs_field_key_id("log"), 0);
    cs_field_set_key_int(bf, cs_field_key_id("post_vis"), 0);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Private function to output property values at cells.
 */
/*----------------------------------------------------------------------------*/

static void
_output_cells_property
(
  int               location_id,
  cs_lnum_t         n_elts,
  const cs_lnum_t  *elt_ids,
  void             *input,
  void             *vals
)
{
  assert(location_id == CS_MESH_LOCATION_CELLS);

  const char *_ppty_name = static_cast<const char *>(input);

  cs_real_t *_vals = static_cast<cs_real_t *>(vals);

  cs_property_t *ppty = cs_property_by_name(_ppty_name);

  cs_real_t _t_cur = cs_glob_time_step->t_cur;
  if (elt_ids != nullptr) {
    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = elt_ids[e_id];
      _vals[e_id] = cs_property_get_cell_value(c_id, _t_cur, ppty);
    }
  }
  else {
    cs_property_eval_at_cells(_t_cur,
                              ppty,
                              _vals);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Private function which computes the vorticity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_vorticity
(
  int              location_id, /*!<[in] location id (mesh location) */
  cs_lnum_t        n_elts,      /*!<[in] number of cells for computation */
  const cs_lnum_t *elt_ids,     /*!<[in] cells' ids indirection */
  void            *input,       /*!<[in] input pointer to cast into cs_field_t * */
  void            *vals         /*!<[in] array of values */
)
{
  /* Sanity check */
  assert(location_id == CS_MESH_LOCATION_CELLS);

  /* Check for single of multiphase case */
  cs_field_t *f = static_cast<cs_field_t *>(input);

  cs_real_33_t *gradv = nullptr;

  /* Use field gradient if available, otherwise recompute */
  if (f->grad != nullptr)
    gradv = reinterpret_cast<cs_real_33_t *>(f->grad);
  else {
    const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
    CS_MALLOC(gradv, n_cells_ext, cs_real_33_t);

    int inc = 1;
    cs_field_gradient_vector(f,
                             false, /* use previous */
                             inc,
                             gradv);
  }

  cs_real_t *vorticity = static_cast<cs_real_t *>(vals);

  if (elt_ids != nullptr) {
    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = elt_ids[e_id];
      cs_real_t *vort = vorticity + 3 * e_id;

      vort[0] = gradv[c_id][2][1] - gradv[c_id][1][2];
      vort[1] = gradv[c_id][0][2] - gradv[c_id][2][0];
      vort[2] = gradv[c_id][1][0] - gradv[c_id][0][1];
    }
  }
  else {
    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_real_t *vort = vorticity + 3 * e_id;
      vort[0] = gradv[e_id][2][1] - gradv[e_id][1][2];
      vort[1] = gradv[e_id][0][2] - gradv[e_id][2][0];
      vort[2] = gradv[e_id][1][0] - gradv[e_id][0][1];
    }
  }

  /* Free pointer if needed */
  if (f->grad == nullptr)
    CS_FREE(gradv);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define functions based on code_saturne case setup.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_default_define(void)
{
  if (cs_glob_n_ranks > 1) {
    cs_function_define_mpi_rank_id(CS_MESH_LOCATION_CELLS);
    cs_function_define_mpi_rank_id(CS_MESH_LOCATION_BOUNDARY_FACES);
    /* cs_function_define_mpi_rank_id(CS_MESH_LOCATION_VERTICES); */
  }

  cs_function_define_by_func("boundary_zone_class_id",
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             1,
                             false,
                             CS_INT_TYPE,
                             cs_function_class_or_zone_id,
                             nullptr);

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE)
    cs_turbomachinery_define_functions();

  if (cs_glob_physical_constants->icorio > 0)
    _define_coriolis_functions();

  if (   cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0
      || cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0)
    cs_elec_define_functions();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access a function for evaluation of element's MPI rank id.
 *
 * \param[in]   location_id  base associated mesh location id
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_mpi_rank_id(cs_mesh_location_type_t  location_id)
{
  const char base_name[] = "mpi_rank_id";
  const char *loc_name = cs_mesh_location_get_name(location_id);

  size_t l_name = strlen(loc_name) + strlen(base_name) + 1;
  char *name;
  CS_MALLOC(name, l_name + 1, char);
  snprintf(name, l_name, "%s_%s", base_name, loc_name);

  cs_function_t *f
    = cs_function_define_by_func(name,
                                 location_id,
                                 1,
                                 false,
                                 CS_INT_TYPE,
                                 _location_mpi_rank_id,
                                 cs_glob_mesh);

  CS_FREE(name);

  /* Use a different label for vertex data and element data, to avoid
     conflicts when outputting values with some writer formats,
     which do not accept 2 fields of the same name on different locations */

  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(location_id);
  if (loc_type != CS_MESH_LOCATION_VERTICES) {
    CS_MALLOC(f->label, strlen(base_name) + 1, char);
    strcpy(f->label, base_name);
  }
  else {
    const char base_name_v[] = "mpi_rank_id_v";
    CS_MALLOC(f->label, strlen(base_name_v) + 1, char);
    strcpy(f->label, base_name_v);
  }

  f->type = 0;
  if (cs_glob_mesh->time_dep < CS_MESH_TRANSIENT_CONNECT)
    f->type |= CS_FUNCTION_TIME_INDEPENDENT;

  // Before activating for cells and boundary faces, remove
  // post_mesh->post_domain feature from cs_post.c.

  if (   location_id != CS_MESH_LOCATION_CELLS
      && location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
    f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access a function for evaluation of mesh element's
 *        refinement generation (i.e. level).
 *
 * \param[in]   location_id  base associated mesh location id
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_refinement_generation(cs_mesh_location_type_t  location_id)
{
  const char base_name[] = "r_gen";
  const char *loc_name = cs_mesh_location_get_name(location_id);

  size_t l_name = strlen(loc_name) + strlen(base_name) + 1;
  char *name;
  CS_MALLOC(name, l_name + 1, char);
  snprintf(name, l_name, "%s_%s", base_name, loc_name);

  cs_function_t *f
    = cs_function_define_by_func(name,
                                 location_id,
                                 1,
                                 false,
                                 CS_INT_TYPE,
                                 _location_r_gen,
                                 cs_glob_mesh);

  CS_FREE(name);

  /* Use a different label for vertex data and element data, to avoid
     conflicts when outputting values with some writer formats,
     which do not accept 2 fields of the same name on different locations */

  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(location_id);
  if (loc_type != CS_MESH_LOCATION_VERTICES) {
    CS_MALLOC(f->label, strlen(base_name) + 1, char);
    strcpy(f->label, base_name);
  }
  else {
    const char base_name_v[] = "r_gen_v";
    CS_MALLOC(f->label, strlen(base_name_v) + 1, char);
    strcpy(f->label, base_name_v);
  }

  f->type = 0;
  if (cs_glob_mesh->time_dep < CS_MESH_TRANSIENT_CONNECT)
    f->type |= CS_FUNCTION_TIME_INDEPENDENT;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of normal boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress_normal(void)
{
  /* Create appropriate fields if needed. */
  _ensure_boundary_stress_is_present();

  cs_function_t *f = nullptr;

  f = cs_function_define_by_func("boundary_stress_normal",
                                 CS_MESH_LOCATION_BOUNDARY_FACES,
                                 1,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_boundary_stress_normal,
                                 nullptr);

  cs_function_set_label(f, "Normal Stress");

  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of tangential boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress_tangential(void)
{
  /* Create appropriate fields if needed. */
  _ensure_boundary_stress_is_present();

  cs_function_t *f = nullptr;

  f = cs_function_define_by_func("boundary_stress_tangential",
                                 CS_MESH_LOCATION_BOUNDARY_FACES,
                                 3,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_boundary_stress_tangential,
                                 cs_glob_mesh);

  cs_function_set_label(f, "Shear Stress");

  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of boundary thermal flux.
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_thermal_flux(void)
{
  cs_function_t *f = nullptr;

  /* Create appropriate fields if needed;
     check that a thermal variable is present first */

  cs_field_t *f_t = cs_thermal_model_field();

  if (f_t == nullptr)
    return f;

  const char *names[] = {"tplus", "tstar"};

  for (int i = 0; i < 2; i++) {

    cs_field_t *bf = cs_field_by_name_try(names[i]);
    if (bf == nullptr) {
      int type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
      int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;

      bf = cs_field_create(names[i], type, location_id, 1, false);
      cs_field_set_key_int(bf, cs_field_key_id("log"), 0);
      cs_field_set_key_int(bf, cs_field_key_id("post_vis"), 0);
    }
  }

  f = cs_function_define_by_func("boundary_thermal_flux",
                                 CS_MESH_LOCATION_BOUNDARY_FACES,
                                 1,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_boundary_thermal_flux,
                                 nullptr);

  cs_function_set_label(f, "Input thermal flux");

  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function for computation of boundary layer Nusselt.
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_nusselt(void)
{
  cs_function_t *f = nullptr;

  /* Create appropriate fields if needed;
     check that a thermal variable is present first */

  cs_field_t *f_t = cs_thermal_model_field();

  if (f_t == nullptr)
    return f;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_t);

  if (eqp->idiff != 0) {

    const char *names[] = {"tplus", "tstar"};

    for (int i = 0; i < 2; i++) {

      cs_field_t *bf = cs_field_by_name_try(names[i]);
      if (bf == nullptr) {
        int type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
        int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;

        bf = cs_field_create(names[i], type, location_id, 1, false);
        cs_field_set_key_int(bf, cs_field_key_id("log"), 0);
        cs_field_set_key_int(bf, cs_field_key_id("post_vis"), 0);
      }
    }

    f = cs_function_define_by_func("boundary_layer_nusselt",
                                   CS_MESH_LOCATION_BOUNDARY_FACES,
                                   1,
                                   true,
                                   CS_REAL_TYPE,
                                   cs_function_boundary_nusselt,
                                   nullptr);

    cs_function_set_label(f, "Dimensionless heat flux");

    f->type = CS_FUNCTION_INTENSIVE;

    f->post_vis = CS_POST_ON_LOCATION;

  }

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function for computation of cell Q criterion.
 *
 * \return  pointer to the associated function object in case of success,
 *          or null in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_q_criterion(void)
{
  cs_function_t *f
    = cs_function_define_by_func("q_criterion",
                                 CS_MESH_LOCATION_CELLS,
                                 1,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_q_criterion,
                                 nullptr);

  cs_function_set_label(f, "Q criterion");

  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract optional boundary face class of element zone id.
 *
 * For boundary faces, if no face classes have been defined by
 * \ref cs_boundary_zone_face_class_id the highest boundary face zone id is
 *
 * For cells, the highest cell volume zone id is used.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        pointer to field
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_class_or_zone_id(int               location_id,
                             cs_lnum_t         n_elts,
                             const cs_lnum_t  *elt_ids,
                             void             *input,
                             void             *vals)
{
  CS_UNUSED(input);

  int       *z_id     = static_cast<int *>(vals);
  const int *elt_z_id = nullptr;;

  if (location_id == CS_MESH_LOCATION_CELLS)
    elt_z_id = cs_volume_zone_cell_zone_id();
  else if (location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
    elt_z_id = cs_boundary_zone_face_class_or_zone_id();

  if (elt_z_id != nullptr) {
    if (elt_ids != nullptr) {
      for (cs_lnum_t i = 0; i < n_elts; i++)
        z_id[i] = elt_z_id[elt_ids[i]];
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++)
        z_id[i] = elt_z_id[i];
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++)
      z_id[i] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute non-reconstructed cell-based field values at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        pointer to field
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_field_boundary_nr(int               location_id,
                              cs_lnum_t         n_elts,
                              const cs_lnum_t  *elt_ids,
                              void             *input,
                              void             *vals)
{
  cs_assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_field_t *f            = static_cast<const cs_field_t *>(input);
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  /* For higher dimensions, check if components are coupled */

  int coupled = 0;
  if (   f->type & CS_FIELD_VARIABLE
      && f->dim > 1
      && f->bc_coeffs != nullptr) {
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > -1)
      coupled = cs_field_get_key_int(f, coupled_key_id);
  }

  cs_real_t       *_vals  = static_cast<cs_real_t *>(vals);
  const cs_real_t *c_vals = f->val;

  /* Scalar fields
     ------------- */

  if (f->dim == 1) {

    if (f->bc_coeffs != nullptr) { /* Variable field with BC definitions */

      const cs_real_t *coefa = f->bc_coeffs->a;
      const cs_real_t *coefb = f->bc_coeffs->b;

      if (elt_ids != nullptr) {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t face_id = elt_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];
          _vals[i] = coefa[face_id] + coefb[face_id]*c_vals[cell_id];
        }
      }
      else {
        for (cs_lnum_t face_id = 0; face_id < n_elts; face_id++) {
          cs_lnum_t cell_id = b_face_cells[face_id];
          _vals[face_id] = coefa[face_id] + coefb[face_id]*c_vals[cell_id];
        }
      }

    }

    else { /* Simply use cell values when no BC's available */

      if (elt_ids != nullptr) {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t face_id = elt_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];
          _vals[i] = c_vals[cell_id];
        }
      }
      else {
        for (cs_lnum_t face_id = 0; face_id < n_elts; face_id++) {
          cs_lnum_t cell_id = b_face_cells[face_id];
          _vals[face_id] = c_vals[cell_id];
        }
      }
    }

  }

  /* Coupled vector and tensor fields (with BC definitions)
   -------------------------------------------------------- */

  else if (coupled) {

    const cs_real_t *coefa = f->bc_coeffs->a;
    const cs_real_t *coefb = f->bc_coeffs->b;

    const cs_lnum_t dim = f->dim;
    const cs_lnum_t dim2 = f->dim * f->dim;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t face_id = (elt_ids != nullptr) ? elt_ids[i] : i;
      cs_lnum_t cell_id = b_face_cells[face_id];
      for (cs_lnum_t j = 0; j < 3; j++)
        _vals[i*dim + j] = coefa[face_id*dim + j];
      for (cs_lnum_t j = 0; j < 3; j++) {
        for (cs_lnum_t k = 0; k < 3; k++)
          _vals[i*dim + j] +=   coefb[face_id*dim2 + k*dim + j]
                              * c_vals[cell_id*dim + k];
      }
    }

  }

  /* Uncoupled vector and tensor fields (with BC definitions)
     -------------------------------------------------------- */

  else if (f->bc_coeffs != nullptr) {

    const cs_real_t *coefa = f->bc_coeffs->a;
    const cs_real_t *coefb = f->bc_coeffs->b;

    const cs_lnum_t dim = f->dim;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t face_id = (elt_ids != nullptr) ? elt_ids[i] : i;
      cs_lnum_t cell_id = b_face_cells[face_id];
      for (cs_lnum_t j = 0; j < 3; j++)
        _vals[i*dim + j] =   coefa[face_id*dim + j]
                           + (  coefb[face_id*dim + j]
                              * c_vals[cell_id*dim + j]);
    }

  }

  /* Vector and tensor fields without BC definitions
     ----------------------------------------------- */

  else {

    const cs_lnum_t dim = f->dim;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t face_id = (elt_ids != nullptr) ? elt_ids[i] : i;
      cs_lnum_t cell_id = b_face_cells[face_id];
      for (cs_lnum_t j = 0; j < 3; j++)
        _vals[i*dim + j] = c_vals[cell_id*dim + j];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute normal stress at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_stress_normal(int               location_id,
                                   cs_lnum_t         n_elts,
                                   const cs_lnum_t  *elt_ids,
                                   void             *input,
                                   void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_t *stress = static_cast<cs_real_t *>(vals);

  const cs_real_3_t *b_stress
    = (const cs_real_3_t *)(cs_field_by_name("boundary_stress")->val);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_nreal_3_t *restrict b_face_u_normal = mq->b_face_u_normal;

  if (elt_ids != nullptr) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t e_id = elt_ids[i];
      stress[i] = cs_math_3_dot_product(b_stress[e_id],
                                        b_face_u_normal[e_id]);
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      stress[i] = cs_math_3_dot_product(b_stress[i],
                                        b_face_u_normal[i]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute tangential stress at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_stress_tangential(int               location_id,
                                       cs_lnum_t         n_elts,
                                       const cs_lnum_t  *elt_ids,
                                       void             *input,
                                       void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_3_t *stress = static_cast<cs_real_3_t *>(vals);

  const cs_real_3_t *b_stress
    = (const cs_real_3_t *)(cs_field_by_name("boundary_stress")->val);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_nreal_3_t *restrict b_face_u_normal = mq->b_face_u_normal;

  if (elt_ids != nullptr) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t e_id = elt_ids[i];
      const cs_nreal_t *u_n = b_face_u_normal[e_id];
      cs_real_t s_nor = cs_math_3_dot_product(b_stress[e_id], u_n);
      for (cs_lnum_t j = 0; j < 3; j++)
        stress[i][j] = b_stress[e_id][j] - s_nor*u_n[j];
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_nreal_t *u_n = b_face_u_normal[i];
      cs_real_t s_nor = cs_math_3_dot_product(b_stress[i], u_n);
      for (cs_lnum_t j = 0; j < 3; j++)
        stress[i][j] = (b_stress[i][j] - s_nor*u_n[j]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute thermal flux at boundary (in \f$ W\,m^{-2} \f$),
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_thermal_flux(int               location_id,
                                  cs_lnum_t         n_elts,
                                  const cs_lnum_t  *elt_ids,
                                  void             *input,
                                  void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_t *b_face_flux = static_cast<cs_real_t *>(vals);

  cs_field_t *f_t = cs_thermal_model_field();

  if (f_t != nullptr) {

    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
    const cs_real_t *restrict b_face_surf = fvq->b_face_surf;

    cs_real_t normal[] = {0, 0, 0};

    cs_flux_through_surface(f_t->name,
                            normal,
                            n_elts,
                            0,
                            elt_ids,
                            nullptr,
                            nullptr,
                            b_face_flux,
                            nullptr);

    if (elt_ids != nullptr) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t e_id = elt_ids[i];
        b_face_flux[i] = (b_face_surf[e_id] > 0.) ?
                          b_face_flux[i] / b_face_surf[e_id] : 0.;
      }
    }
    else {
      for (cs_lnum_t f_id = 0; f_id < n_elts; f_id++) {
        b_face_flux[f_id] = (b_face_surf[f_id] > 0.) ?
                             b_face_flux[f_id] / b_face_surf[f_id] : 0.;
      }
    }

  }

  else { /* Default if not available */

    for (cs_lnum_t i = 0; i < n_elts; i++)
      b_face_flux[i] = 0.;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute local Nusselt number near boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_nusselt(int               location_id,
                             cs_lnum_t         n_elts,
                             const cs_lnum_t  *elt_ids,
                             void             *input,
                             void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_t *bnussl = static_cast<cs_real_t *>(vals);

  /* Remarks:
   *
   * This function uses local "boundary-only" reconstruction when possible,
   * to reduce computational cost.
   *
   * A more general solution would be to compute the boundary thermal flux
   * using cs_flux_through_surface(), dividing by the face surfaces,
   * then removing the convective part if present:
   * (b_mass_flux/face_surf)*(coefa + coefb*t_cel)
   *
   * And finally multiplying by:
   * b_face_dist / (xvsl * tplus * tstar)
   *
   * Where xvsl is the thermal diffusivity (uniform or not) of the adjancent cell.
   *
   * This would present the advantage of factoring more code, but in this case,
   * a boundary-only version of the cs_flux_through_surface function could be useful
   * so as to allow computing gradients only at the boundary (at least when
   * least-squares are used).
   */

  /* T+ and T* if saved */

  const cs_field_t *f_tplus = cs_field_by_name_try("tplus");
  const cs_field_t *f_tstar = cs_field_by_name_try("tstar");

  if (f_tplus != nullptr && f_tstar != nullptr) {

    cs_field_t *f_t = cs_thermal_model_field();
    cs_real_t *tscalp = f_t->val_pre;

    cs_real_t *tplus = f_tplus->val;
    cs_real_t *tstar = f_tstar->val;

    /* Boundary condition pointers for diffusion */

    cs_real_t *cofaf = f_t->bc_coeffs->af;
    cs_real_t *cofbf = f_t->bc_coeffs->bf;

    /* Boundary condition pointers for diffusion with coupling */

    cs_real_t *rcodcl2 = f_t->bc_coeffs->rcodcl2;
    cs_real_t *hint = f_t->bc_coeffs->hint;

    /* Compute variable values at boundary faces */

    const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

    const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_t);

    const bool *is_coupled = nullptr;

    cs_real_t *theipb = nullptr, *dist_theipb = nullptr;
    CS_MALLOC(theipb, n_elts, cs_real_t);

    /* Reconstructed fluxes */

    if (   eqp->ircflu > 0 && eqp->b_diff_flux_rc > 0
        && cs_glob_space_disc->itbrrb == 1) {

      cs_field_gradient_boundary_iprime_scalar(f_t,
                                               false, /* use_previous_t */
                                               n_elts,
                                               elt_ids,
                                               theipb);

      /* In previous (Fortran) version, we used the previous value for
         the thermal scalar, with the current gradient. This might be
         an error, but for now, add term to obtain similar behavior... */

      cs_real_t *tscal = f_t->val;

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = (elt_ids != nullptr) ? elt_ids[i] : i;
        cs_lnum_t k = b_face_cells[j];
        theipb[i] += tscalp[k] - tscal[k];
      }

    }
    else {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = (elt_ids != nullptr) ? elt_ids[i] : i;
        theipb[i] = tscalp[b_face_cells[j]];
      }

    }

    /* Special case for internal coupling */

    if (eqp->icoupl > 0) {

      cs_real_t *loc_theipb;
      CS_MALLOC(loc_theipb, n_b_faces, cs_real_t);
      CS_MALLOC(dist_theipb, n_b_faces, cs_real_t);

      for (cs_lnum_t i = 0; i < n_b_faces; i++)
        loc_theipb[i] = 0;

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = (elt_ids != nullptr) ? elt_ids[i] : i;
        loc_theipb[j] = theipb[i];
      }

      const int coupling_key_id = cs_field_key_id("coupling_entity");
      int coupling_id = cs_field_get_key_int(f_t, coupling_key_id);
      const cs_internal_coupling_t  *cpl
        = cs_internal_coupling_by_id(coupling_id);

      is_coupled = cpl->coupled_faces;

      cs_ic_field_dist_data_by_face_id(f_t->id, 1, loc_theipb, dist_theipb);

      CS_FREE(loc_theipb);

    }

    /* Physical property pointers */

    const int kivisl = cs_field_key_id("diffusivity_id");
    const int diff_id = cs_field_get_key_int(f_t, kivisl);

    cs_real_t visls_0 = -1;
    cs_lnum_t viscl_step = 0;

    const cs_real_t *cviscl = &visls_0;

    if (diff_id > -1) {
      cviscl = cs_field_by_id(diff_id)->val;
      viscl_step = 1;
    }
    else {
      const int kvisls0 = cs_field_key_id("diffusivity_ref");
      visls_0 = cs_field_get_key_double(f_t, kvisls0);
    }

    /* Compute using reconstructed temperature value in boundary cells */

    const cs_real_t *b_dist = fvq->b_dist;

    bool have_coupled = (   eqp->icoupl > 0
                         && (  cs_glob_time_step->nt_cur
                             > cs_glob_time_step->nt_prev));

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_lnum_t face_id = (elt_ids != nullptr) ? elt_ids[i] : i;
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t numer =   (cofaf[face_id] + cofbf[face_id]*theipb[i])
                        * b_dist[face_id];

      /* here numer = 0 if current face is coupled.
         FIXME exchange coefs not computed at start of calculation */

      if (have_coupled) {
        if (is_coupled[face_id]) {
          cs_real_t heq =   rcodcl2[face_id] * hint[face_id]
                          / (rcodcl2[face_id] + hint[face_id]);
          numer = heq*(theipb[i]-dist_theipb[face_id]) * b_dist[face_id];
        }
      }

      cs_real_t xvsl = cviscl[cell_id * viscl_step];

      cs_real_t denom = xvsl * tplus[face_id] * tstar[face_id];

      if (fabs(denom) > 1e-30)
        bnussl[i] = numer / denom;
      else
        bnussl[i] = 0.;

    }

    CS_FREE(dist_theipb);
    CS_FREE(theipb);

  }
  else { /* Default if not computable */

    for (cs_lnum_t i = 0; i < n_elts; i++)
      bnussl[i] = -1.;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Q-criterion from Hunt et. al over each cell of a specified
 *        volume region.
 *
 * \f[
 *    Q = \tens{\Omega}:\tens{\Omega} -
 *    \deviator{ \left(\tens{S} \right)}:\deviator{ \left(\tens{S} \right)}
 * \f]
 * where \f$\tens{\Omega}\f$ is the vorticity tensor and
 * \f$\deviator{ \left(\tens{S} \right)}\f$ the deviatoric of the rate of strain
 * tensor.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or null if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_q_criterion(int               location_id,
                        cs_lnum_t         n_elts,
                        const cs_lnum_t  *elt_ids,
                        void             *input,
                        void             *vals)
{
  CS_UNUSED(input);

  cs_assert(location_id == CS_MESH_LOCATION_CELLS);

  cs_real_t *q_crit = static_cast<cs_real_t *>(vals);

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_real_33_t *gradv;
  CS_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(cs_field_by_name("velocity"),
                           false,  /* use_previous_t */
                           1,      /* not on increment */
                           gradv);

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t c_id = (elt_ids != nullptr) ? elt_ids[i] : i;
    q_crit[i] = -1./6. * (   cs_math_sq(gradv[c_id][0][0])
                          +  cs_math_sq(gradv[c_id][1][1])
                          +  cs_math_sq(gradv[c_id][2][2]))
                - gradv[c_id][0][1]*gradv[c_id][1][0]
                - gradv[c_id][0][2]*gradv[c_id][2][0]
                - gradv[c_id][1][2]*gradv[c_id][2][1];
  }

  CS_FREE(gradv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define output function of a property (which is not a field).
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_property_cells
(
  char *property_name /*!<[in] name of the property */
)
{
  /* sanity check */
  cs_property_t *ppty = cs_property_by_name(property_name);
  if (ppty == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: property \"%s\" does not exist.\n",
              __func__, property_name);

  int _dim = cs_property_get_dim(ppty);
  cs_function_t *f
    = cs_function_define_by_func(property_name,
                                 CS_MESH_LOCATION_CELLS,
                                 _dim,
                                 true,
                                 CS_REAL_TYPE,
                                 _output_cells_property,
                                 property_name);

  cs_function_set_label(f, property_name);
  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define output function for the vorticity.
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_vorticity
(
  cs_field_t *velocity_field /*!<[in] Pointer to velocity field */
)
{
  if (velocity_field == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Null pointer provided for the velocity field.\n", __func__);

  std::string f_name = "vorticity::" + std::string(velocity_field->name);
  cs_function_t *f
    = cs_function_define_by_func(f_name.c_str(),
                                 CS_MESH_LOCATION_CELLS,
                                 3,
                                 true,
                                 CS_REAL_TYPE,
                                 _compute_vorticity,
                                 velocity_field);

  cs_function_set_label(f, f_name.c_str());
  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
