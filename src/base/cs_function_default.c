/*============================================================================
 * Base predefined function objects.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_log.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_assert.h"
#include "cs_balance_by_zone.h"
#include "cs_boundary_zone.h"
#include "cs_elec_model.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_function.h"
#include "cs_function_default.h"
#include "cs_mesh_quantities.h"
#include "cs_math.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_rotation.h"
#include "cs_thermal_model.h"
#include "cs_turbomachinery.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_function_default.c
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
 * \param[in]       rs           pointer to range set structure, or NULL
 * \param[in]       location_id  base associated mesh location idà
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  int *e_rank_id;

  if (n_elts != n_loc_elts || elt_ids != NULL)
    BFT_MALLOC(e_rank_id, n_loc_elts, int);
  else
    e_rank_id = vals;

  if (rs != NULL) {
    for (cs_lnum_t i = 0; i < rs->n_elts[0]; i++)
      e_rank_id[i] = cs_glob_rank_id;
    for (cs_lnum_t i = rs->n_elts[0]; i < n_loc_elts; i++)
      e_rank_id[i] = 0;

    cs_range_set_scatter(rs,
                         CS_INT_TYPE,
                         1,
                         e_rank_id,
                         e_rank_id);

    if (rs->ifs != NULL)
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

  if (e_rank_id != vals) {
    int *_vals = vals;
    for (cs_lnum_t i = 0; i < n_elts; i++)
      _vals[i] = e_rank_id[elt_ids[i]];

    BFT_FREE(e_rank_id);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the associated rank at a given mesh location.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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
  case CS_MESH_LOCATION_INTERIOR_FACES:
    {
      cs_mesh_t *m = input;

      cs_gnum_t *g_i_face_num = m->global_i_face_num;
      if (g_i_face_num == NULL) {
        cs_lnum_t n_i_faces = m->n_i_faces;
        BFT_MALLOC(g_i_face_num, n_i_faces, cs_gnum_t);
        for (cs_lnum_t i = 0; i < n_i_faces; i++)
          g_i_face_num[i] = (cs_gnum_t)i+1;
      }

      cs_interface_set_t *face_interfaces
        = cs_interface_set_create(m->n_i_faces,
                                  NULL,
                                  g_i_face_num,
                                  m->periodicity,
                                  0,
                                  NULL,
                                  NULL,
                                  NULL);

      if (m->global_i_face_num != g_i_face_num)
        BFT_FREE(g_i_face_num);

      cs_range_set_t *rs = cs_range_set_create(face_interfaces,
                                               NULL,
                                               m->n_i_faces,
                                               false, /* balance */
                                               2,  /* tr_ignore */
                                               0); /* g_id_base */

      _range_set_mpi_rank_id(rs, location_id, n_elts, elt_ids, vals);

      cs_range_set_destroy(&rs);
      cs_interface_set_destroy(&face_interfaces);
    }
    break;

  case CS_MESH_LOCATION_VERTICES:
    {
      cs_mesh_t *m = input;
      cs_range_set_t *rs = m->vtx_range_set;

      if (rs == NULL)
        rs = cs_range_set_create(m->vtx_interfaces,
                                 NULL,
                                 m->n_vertices,
                                 false, /* balance */
                                 2,  /* tr_ignore */
                                 0); /* g_id_base */

      _range_set_mpi_rank_id(rs, location_id, n_elts, elt_ids, vals);

      if (rs != m->vtx_range_set)
        cs_range_set_destroy(&rs);
    }
    break;

  default:
    {
      int *_vals = vals;
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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
  const cs_mesh_t *m = input;
  const char *e_r_gen = NULL;

  char *c_r_gen = NULL;
  int *r_gen = vals;

  /* Locations with stored refinement generations */

  if (location_id == CS_MESH_LOCATION_INTERIOR_FACES)
    e_r_gen = m->i_face_r_gen;
  else if (location_id == CS_MESH_LOCATION_VERTICES)
    e_r_gen = m->vtx_r_gen;

  /* other base locations */

  if (   location_id == CS_MESH_LOCATION_CELLS
      || location_id == CS_MESH_LOCATION_BOUNDARY_FACES) {

    BFT_MALLOC(c_r_gen, m->n_cells_with_ghosts, char);

    for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++)
      c_r_gen[i] = 0;

    /* Note: when mesh/face adjacencies are available,
       a cell-based loop would allow easier threading. */

    const cs_lnum_2_t *restrict i_face_cells
      = (const cs_lnum_2_t *restrict)m->i_face_cells;

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
      if (elt_ids != NULL) {
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
      BFT_FREE(c_r_gen);
      return;

    }

  }

  /* Now apply stored or generated refinement generation
     (note that for boundary faces, processing has already been done
     and we have returned prior to reaching this point) */

  if (e_r_gen != NULL) {
    if (elt_ids != NULL) {
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

  BFT_FREE(c_r_gen);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the absolute pressure associated to the given cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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
  const int *c_r_num = NULL;
  cs_lnum_t c_r_step = 0;

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE) {
    c_r_num = cs_turbomachinery_get_cell_rotor_num();
    c_r_step = 1;
  }

  cs_real_t *p_abs = vals;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)(mq->cell_cen);
  const cs_real_t *cvar_pr = cs_field_by_name("pressure")->val;
  const cs_real_t *cpro_rho = cs_field_by_name("density")->val;

  if (elt_ids != NULL) {
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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
  const int *c_r_num = NULL;
  cs_lnum_t c_r_step = 0;

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE) {
    c_r_num = cs_turbomachinery_get_cell_rotor_num();
    c_r_step = 1;
  }

  cs_real_3_t *v_abs = vals;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)(mq->cell_cen);
  const cs_real_3_t *cvar_vel
    = (const cs_real_3_t *)cs_field_by_name("velocity")->val;

  if (elt_ids != NULL) {
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
 *          is active, or NULL
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
                                   NULL);

    const char label[] = "Abs Pressure";
    BFT_MALLOC(f->label, strlen(label) + 1, char);
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
                                   NULL);

    const char label[] = "Abs Velocity";
    BFT_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ensure the "boundary_forces" fiels is present.
 */
/*----------------------------------------------------------------------------*/

static void
_ensure_boundary_forces_are_present(void)
{
  const char name[] = "boundary_forces";
  cs_field_t *bf = cs_field_by_name_try(name);
  if (bf == NULL) {
    int type = CS_FIELD_EXTENSIVE | CS_FIELD_POSTPROCESS;
    int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;

    bf = cs_field_create(name, type, location_id, 3, false);
    cs_field_set_key_int(bf, cs_field_key_id("log"), 0);
    cs_field_set_key_int(bf, cs_field_key_id("post_vis"), 0);
  }
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
                             NULL);

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
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_mpi_rank_id(cs_mesh_location_type_t  location_id)
{
  const char base_name[] = "mpi_rank_id";
  const char *loc_name = cs_mesh_location_get_name(location_id);

  size_t l_name = strlen(loc_name) + strlen(base_name) + 1;
  char *name;
  BFT_MALLOC(name, l_name + 1, char);
  snprintf(name, l_name, "%s_%s", base_name, loc_name);

  cs_function_t *f
    = cs_function_define_by_func(name,
                                 location_id,
                                 1,
                                 false,
                                 CS_INT_TYPE,
                                 _location_mpi_rank_id,
                                 cs_glob_mesh);

  BFT_FREE(name);

  /* Use a different label for vertex data and element data, to avoid
     conflicts when outputting values with some writer formats,
     which do not accept 2 fields of the same name on different locations */

  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(location_id);
  if (loc_type != CS_MESH_LOCATION_VERTICES) {
    BFT_MALLOC(f->label, strlen(base_name) + 1, char);
    strcpy(f->label, base_name);
  }
  else {
    const char base_name_v[] = "mpi_rank_id_v";
    BFT_MALLOC(f->label, strlen(base_name_v) + 1, char);
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
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_refinement_generation(cs_mesh_location_type_t  location_id)
{
  const char base_name[] = "r_gen";
  const char *loc_name = cs_mesh_location_get_name(location_id);

  size_t l_name = strlen(loc_name) + strlen(base_name) + 1;
  char *name;
  BFT_MALLOC(name, l_name + 1, char);
  snprintf(name, l_name, "%s_%s", base_name, loc_name);

  cs_function_t *f
    = cs_function_define_by_func(name,
                                 location_id,
                                 1,
                                 false,
                                 CS_INT_TYPE,
                                 _location_r_gen,
                                 cs_glob_mesh);

  BFT_FREE(name);

  /* Use a different label for vertex data and element data, to avoid
     conflicts when outputting values with some writer formats,
     which do not accept 2 fields of the same name on different locations */

  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(location_id);
  if (loc_type != CS_MESH_LOCATION_VERTICES) {
    BFT_MALLOC(f->label, strlen(base_name) + 1, char);
    strcpy(f->label, base_name);
  }
  else {
    const char base_name_v[] = "r_gen_v";
    BFT_MALLOC(f->label, strlen(base_name_v) + 1, char);
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
 * \brief Define function object for computation of boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress(void)
{
  /* Create appropriate fields if needed. */
  _ensure_boundary_forces_are_present();

  cs_function_t *f = NULL;

  f = cs_function_define_by_func("boundary_stress",
                                 CS_MESH_LOCATION_BOUNDARY_FACES,
                                 3,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_boundary_stress,
                                 cs_glob_mesh);

  cs_function_set_label(f, "Stress");

  f->type = CS_FUNCTION_INTENSIVE;

  f->post_vis = CS_POST_ON_LOCATION;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of normal boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress_normal(void)
{
  /* Create appropriate fields if needed. */
  _ensure_boundary_forces_are_present();

  cs_function_t *f = NULL;

  f = cs_function_define_by_func("boundary_stress_normal",
                                 CS_MESH_LOCATION_BOUNDARY_FACES,
                                 1,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_boundary_stress_normal,
                                 NULL);

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
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress_tangential(void)
{
  /* Create appropriate fields if needed. */
  _ensure_boundary_forces_are_present();

  cs_function_t *f = NULL;

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
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_thermal_flux(void)
{
  cs_function_t *f = NULL;

  /* Create appropriate fields if needed;
     check that a thermal variable is present first */

  cs_field_t *f_t = cs_thermal_model_field();

  if (f_t == NULL)
    return f;

  f = cs_function_define_by_func("boundary_thermal_flux",
                                 CS_MESH_LOCATION_BOUNDARY_FACES,
                                 1,
                                 true,
                                 CS_REAL_TYPE,
                                 cs_function_boundary_thermal_flux,
                                 NULL);

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
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_nusselt(void)
{
  cs_function_t *f = NULL;

  /* Create appropriate fields if needed;
     check that a thermal variable is present first */

  cs_field_t *f_t = cs_thermal_model_field();

  if (f_t == NULL)
    return f;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_t);

  if (eqp->idiff != 0) {

    const char *names[] = {"tplus", "tstar"};

    for (int i = 0; i < 2; i++) {

      cs_field_t *bf = cs_field_by_name_try(names[i]);
      if (bf == NULL) {
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
                                   NULL);

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
 *          or NULL in case of error
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
                                 NULL);

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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  int *z_id = vals;
  const int *elt_z_id = NULL;;

  if (location_id == CS_MESH_LOCATION_CELLS)
    elt_z_id = cs_volume_zone_cell_zone_id();
  else if (location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
    elt_z_id = cs_boundary_zone_face_class_or_zone_id();

  if (elt_z_id != NULL) {
    if (elt_ids != NULL) {
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  const cs_field_t *f = input;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  /* For higher dimensions, check if components are coupled */

  int coupled = 0;
  if (   f->type & CS_FIELD_VARIABLE
      && f->dim > 1
      && f->bc_coeffs != NULL) {
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > -1)
      coupled = cs_field_get_key_int(f, coupled_key_id);
  }

  cs_real_t *_vals = vals;
  const cs_real_t *c_vals = f->val;

  /* Scalar fields
     ------------- */

  if (f->dim == 1) {

    if (f->bc_coeffs != NULL) { /* Variable field with BC definitions */

      const cs_real_t *coefa = f->bc_coeffs->a;
      const cs_real_t *coefb = f->bc_coeffs->b;

      if (elt_ids != NULL) {
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

      if (elt_ids != NULL) {
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
      cs_lnum_t face_id = (elt_ids != NULL) ? elt_ids[i] : i;
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

  else if (f->bc_coeffs != NULL) {

    const cs_real_t *coefa = f->bc_coeffs->a;
    const cs_real_t *coefb = f->bc_coeffs->b;

    const cs_lnum_t dim = f->dim;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t face_id = (elt_ids != NULL) ? elt_ids[i] : i;
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
      cs_lnum_t face_id = (elt_ids != NULL) ? elt_ids[i] : i;
      cs_lnum_t cell_id = b_face_cells[face_id];
      for (cs_lnum_t j = 0; j < 3; j++)
        _vals[i*dim + j] = c_vals[cell_id*dim + j];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute stress at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_stress(int               location_id,
                            cs_lnum_t         n_elts,
                            const cs_lnum_t  *elt_ids,
                            void             *input,
                            void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_3_t *stress = vals;

  const cs_real_3_t *b_forces
    = (const cs_real_3_t *)(cs_field_by_name("boundary_forces")->val);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)mq->b_face_surf;

  if (elt_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t e_id = elt_ids[i];
      cs_real_t s_mult = 1./b_face_surf[e_id];
      for (cs_lnum_t j = 0; j < 3; j++)
        stress[i][j] = b_forces[e_id][j] * s_mult;
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_real_t s_mult = 1./b_face_surf[i];
      for (cs_lnum_t j = 0; j < 3; j++)
        stress[i][j] = b_forces[i][j] * s_mult;
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  cs_real_t *stress = vals;

  const cs_real_3_t *b_forces
    = (const cs_real_3_t *)(cs_field_by_name("boundary_forces")->val);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)mq->b_face_surf;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;

  if (elt_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t e_id = elt_ids[i];
      cs_real_t s_mult_2 = 1. / cs_math_sq(b_face_surf[e_id]);
      cs_real_t f_nor = cs_math_3_dot_product(b_forces[e_id],
                                              b_face_normal[e_id]);
      stress[i] = f_nor * s_mult_2;
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_real_t s_mult_2 = 1. / cs_math_sq(b_face_surf[i]);
      cs_real_t f_nor = cs_math_3_dot_product(b_forces[i],
                                              b_face_normal[i]);
      stress[i] = f_nor * s_mult_2;
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  cs_real_3_t *stress = vals;

  const cs_real_3_t *b_forces
    = (const cs_real_3_t *)(cs_field_by_name("boundary_forces")->val);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)mq->b_face_surf;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;

  if (elt_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t e_id = elt_ids[i];
      cs_real_t s_mult = 1./b_face_surf[e_id];
      cs_real_t s_nor[3] = {b_face_normal[e_id][0] * s_mult,
                            b_face_normal[e_id][1] * s_mult,
                            b_face_normal[e_id][2] * s_mult};
      cs_real_t f_nor = cs_math_3_dot_product(b_forces[e_id], s_nor);
      for (cs_lnum_t j = 0; j < 3; j++)
        stress[i][j] = (b_forces[e_id][j] - f_nor*s_nor[j]) * s_mult;
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_real_t s_mult = 1./b_face_surf[i];
      cs_real_t s_nor[3] = {b_face_normal[i][0] * s_mult,
                            b_face_normal[i][1] * s_mult,
                            b_face_normal[i][2] * s_mult};
      cs_real_t f_nor = cs_math_3_dot_product(b_forces[i], s_nor);
      for (cs_lnum_t j = 0; j < 3; j++)
        stress[i][j] = (b_forces[i][j] - f_nor*s_nor[j]) * s_mult;
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  cs_real_t *b_face_flux = vals;

  cs_field_t *f_t = cs_thermal_model_field();

  if (f_t != NULL) {

    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
    const cs_real_t *restrict b_face_surf = fvq->b_face_surf;

    cs_real_t normal[] = {0, 0, 0};

    cs_flux_through_surface(f_t->name,
                            normal,
                            n_elts,
                            0,
                            elt_ids,
                            NULL,
                            NULL,
                            b_face_flux,
                            NULL);

    if (elt_ids != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t e_id = elt_ids[i];
        b_face_flux[i] /= b_face_surf[e_id];
      }
    }
    else {
      for (cs_lnum_t f_id = 0; f_id < n_elts; f_id++) {
        b_face_flux[f_id] /= b_face_surf[f_id];
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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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

  cs_real_t *bnussl = vals;

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

  if (f_tplus != NULL && f_tstar != NULL) {

    cs_field_t *f_t = cs_thermal_model_field();
    cs_real_t *tscalp = f_t->val_pre;

    cs_real_t *tplus = f_tplus->val;
    cs_real_t *tstar = f_tstar->val;

    /* Boundary condition pointers for diffusion */

    cs_real_t *cofaf = f_t->bc_coeffs->af;
    cs_real_t *cofbf = f_t->bc_coeffs->bf;

    /* Boundary condition pointers for diffusion with coupling */

    cs_real_t *hext = f_t->bc_coeffs->hext;
    cs_real_t *hint = f_t->bc_coeffs->hint;

    /* Compute variable values at boundary faces */

    const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

    const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_t);

    const bool *is_coupled = NULL;

    cs_real_t *theipb = NULL, *dist_theipb = NULL;
    BFT_MALLOC(theipb, n_elts, cs_real_t);

    /* Reconstructed fluxes */

    if (eqp->ircflu > 0 && cs_glob_space_disc->itbrrb == 1) {

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
        cs_lnum_t j = (elt_ids != NULL) ? elt_ids[i] : i;
        cs_lnum_t k = b_face_cells[j];
        theipb[i] += tscalp[k] - tscal[k];
      }

    }
    else {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = (elt_ids != NULL) ? elt_ids[i] : i;
        theipb[i] = tscalp[b_face_cells[j]];
      }

    }

    /* Special case for internal coupling */

    if (eqp->icoupl > 0) {

      cs_real_t *loc_theipb;
      BFT_MALLOC(loc_theipb, n_b_faces, cs_real_t);
      BFT_MALLOC(dist_theipb, n_b_faces, cs_real_t);

      for (cs_lnum_t i = 0; i < n_b_faces; i++)
        loc_theipb[i] = 0;

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = (elt_ids != NULL) ? elt_ids[i] : i;
        loc_theipb[j] = theipb[i];
      }

      const int coupling_key_id = cs_field_key_id("coupling_entity");
      int coupling_id = cs_field_get_key_int(f_t, coupling_key_id);
      const cs_internal_coupling_t  *cpl
        = cs_internal_coupling_by_id(coupling_id);

      is_coupled = cpl->coupled_faces;

      cs_ic_field_dist_data_by_face_id(f_t->id, 1, loc_theipb, dist_theipb);

      BFT_FREE(loc_theipb);

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
    const cs_real_t *srfbn = fvq->b_f_face_surf;

    bool have_coupled = (   eqp->icoupl > 0
                         && (  cs_glob_time_step->nt_cur
                             > cs_glob_time_step->nt_prev));

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_lnum_t face_id = (elt_ids != NULL) ? elt_ids[i] : i;
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t numer =   (cofaf[face_id] + cofbf[face_id]*theipb[i])
                        * b_dist[face_id];

      /* here numer = 0 if current face is coupled.
         FIXME exchange coefs not computed at start of calculation */

      if (have_coupled) {
        if (is_coupled[face_id]) {
          cs_real_t heq =   hext[face_id] * hint[face_id]
                          / ((hext[face_id] + hint[face_id]) * srfbn[face_id]);
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

    BFT_FREE(dist_theipb);
    BFT_FREE(theipb);

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
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
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
  cs_assert(location_id == CS_MESH_LOCATION_CELLS);

  cs_real_t  *q_crit = vals;

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_real_33_t *gradv;
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(cs_field_by_name("velocity"),
                           false,  /* use_previous_t */
                           1,      /* not on increment */
                           gradv);

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t c_id = (elt_ids != NULL) ? elt_ids[i] : i;
    q_crit[i] = -1./6. * (   cs_math_sq(gradv[c_id][0][0])
                          +  cs_math_sq(gradv[c_id][1][1])
                          +  cs_math_sq(gradv[c_id][2][2]))
                - gradv[c_id][0][1]*gradv[c_id][1][0]
                - gradv[c_id][0][2]*gradv[c_id][2][0]
                - gradv[c_id][1][2]*gradv[c_id][2][1];
  }

  BFT_FREE(gradv);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
