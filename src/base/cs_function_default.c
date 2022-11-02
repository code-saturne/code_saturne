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

#include "cs_function.h"
#include "cs_function_default.h"
#include "cs_mesh_quantities.h"
#include "cs_math.h"
#include "cs_post.h"
#include "cs_physical_constants.h"
#include "cs_rotation.h"
#include "cs_turbomachinery.h"

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
                                   false,
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
                                   false,
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

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE)
    cs_turbomachinery_define_functions();

  if (cs_glob_physical_constants->icorio > 0)
    _define_coriolis_functions();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access a function whose data values will be computed
 *        using the a predefined evaluation function.
 *
 * \param[in]   location_id  base associated mesh location id

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

END_C_DECLS
