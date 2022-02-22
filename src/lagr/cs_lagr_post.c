/*============================================================================
 * Lagrangian module postprocessing
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

#include <stddef.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_mesh_location.h"

#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_parall.h"
#include "cs_post.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_post.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local types and structures
 *============================================================================*/

/* Structure associated with postprocessing options */
/*--------------------------------------------------*/

typedef struct {

  /*! \anchor particle_attr
    flag for activation of output for each possible particle attribute:
    0:   not active
    1:   active
    > 1: for atributes with multiple components, number of components
    postprocessed as separate scalars */
  int  attr_output[CS_LAGR_N_ATTRIBUTES];

} cs_lagr_post_options_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool                    _lagr_post_options_is_set = false;

/* Postprocessing options structure and associated pointer */

static cs_lagr_post_options_t  _lagr_post_options
= {.attr_output[0] = -1};

const cs_lagr_post_options_t *cs_glob_lagr_post_options = &_lagr_post_options;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define which Lagragian variables should be postprocessed
 *----------------------------------------------------------------------------*/

static void
_activate_particle_output(void)
{
  /* No output if nothing initialized by now */

  if (_lagr_post_options.attr_output[0] == -1) {
    for (cs_lagr_attribute_t i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
      _lagr_post_options.attr_output[i] = 0;
    }
  }

  else {
    for (cs_lagr_attribute_t i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
      if (_lagr_post_options.attr_output[i]) {
        int count = 0;
        cs_lagr_get_attr_info(cs_glob_lagr_particle_set,
                              0,
                              i,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              &count);

        if (count == 3) {
          switch(i) {
          case CS_LAGR_COORDS:
          case CS_LAGR_VELOCITY:
          case CS_LAGR_VELOCITY_SEEN:
          case CS_LAGR_PRED_VELOCITY:
          case CS_LAGR_PRED_VELOCITY_SEEN:
          case CS_LAGR_ORIENTATION:
          case CS_LAGR_RADII:
          case CS_LAGR_ANGULAR_VEL:
            count = 1;
            break;
          default:
            break;
          }
        }

        _lagr_post_options.attr_output[i] = count;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Default additional particle output of mesh and time-dependent variables
 * for the call to pstvar / cs_post_write_vars.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * parameters:
 *   input       <-> pointer to optional (untyped) value or structure;
 *                   here, we should point to _lagr_post_options.
 *   mesh_id     <-- id of the output mesh for the current call
 *   cat_id      <-- category id of the output mesh for the current call
 *   ts          <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_write_particle_vars(cs_lagr_post_options_t  *options,
                     int                      mesh_id,
                     const cs_time_step_t    *ts)
{
  cs_lagr_attribute_t attr_id;

  char var_name[64];
  int  component_id;
  char var_name_component[96];

  for (attr_id = 0; attr_id < CS_LAGR_N_ATTRIBUTES; attr_id++) {

    if (options->attr_output[attr_id] > 0) {

      /* build name */

      snprintf(var_name,
               63,
               "particle_%s",
               cs_lagr_attribute_name[attr_id]);
      var_name[63] = '\0';

      /* Output values */

      if (options->attr_output[attr_id] == 1)
        cs_post_write_particle_values(mesh_id,
                                      CS_POST_WRITER_ALL_ASSOCIATED,
                                      attr_id,
                                      var_name,
                                      -1,
                                      ts);
      else {
        /* Create one output per component */
        for (component_id = 0;
             component_id < options->attr_output[attr_id];
             component_id++) {
          snprintf(var_name_component,
                   95,
                   "%s_layer_%2.2i",
                   var_name,
                   component_id+1);
          var_name_component[95] = '\0';
          cs_post_write_particle_values(mesh_id,
                                        CS_POST_WRITER_ALL_ASSOCIATED,
                                        attr_id,
                                        var_name_component,
                                        component_id,
                                        ts);
        }
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Function for additional postprocessing of Lagrangian data.
 *
 * This function should match the prototype of a
 * (cs_post_time_mesh_dep_output_t) function, and be registered using
 * cs_post_add_time_mesh_dep_output().
 *
 * parameters:
 *   input       <-> pointer to optional (untyped) value or structure.
 *   mesh_id     <-- id of the output mesh for the current call
 *   cat_id      <-- category id of the output mesh for the current call
 *   ent_flag    <-- indicate global presence of cells (ent_flag[0]), interior
 *                   faces (ent_flag[1]), boundary faces (ent_flag[2]),
 *                   particles (ent_flag[3]) or probes (ent_flag[4])
 *   n_cells     <-- local number of cells of post_mesh
 *   n_i_faces   <-- local number of interior faces of post_mesh
 *   n_b_faces   <-- local number of boundary faces of post_mesh
 *   cell_ids    <-- list of cells (0 to n-1) of post-processing mesh
 *   i_face_ids  <-- list of interior faces (0 to n-1) of post-processing mesh
 *   b_face_ids  <-- list of boundary faces (0 to n-1) of post-processing mesh
 *   ts          <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_lagr_post(void                  *input,
              int                    mesh_id,
              int                    cat_id,
              int                    ent_flag[5],
              cs_lnum_t              n_cells,
              cs_lnum_t              n_i_faces,
              cs_lnum_t              n_b_faces,
              const cs_lnum_t        cell_ids[],
              const cs_lnum_t        i_face_ids[],
              const cs_lnum_t        b_face_ids[],
              const cs_time_step_t  *ts)
{
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);

  /* Specific handling for particle meshes */

  if (cat_id == -3) {
    _write_particle_vars(input, mesh_id, ts);
    return;
  }

  /* Boundary statistics */
  /*---------------------*/

  else if (cat_id == -2 && cs_glob_lagr_time_scheme->iilagr > 0) {

    cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;

    const cs_lagr_boundary_interactions_t
      *lagr_b = cs_glob_lagr_boundary_interactions;

    cs_real_t  threshold = cs_glob_lagr_stat_options->threshold;

    cs_real_t *val;
    BFT_MALLOC(val, nfabor, cs_real_t);

    for (int irf = 0; irf < cs_glob_lagr_dim->n_boundary_stats; irf++) {

      const char *var_name = lagr_b->nombrd[irf];

      const cs_real_t *_b_stats = bound_stat + nfabor*irf;

      const cs_real_t *_f_count = bound_stat + nfabor*lagr_b->inbr;

      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        cs_lnum_t f_id = b_face_ids[i];
        if (_f_count[f_id] > threshold)
          val[i] = _b_stats[f_id];
        else
          val[i] = 0.;
      }

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        var_name,
                        1,       /* var_dim */
                        true,    /* interlace */
                        false,   /* use_parent */
                        CS_POST_TYPE_cs_real_t,
                        NULL,
                        NULL,
                        val,
                        cs_glob_time_step);

    }

    BFT_FREE(val);

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize Lagrangian postprocessing.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_post_init(void)
{
  _activate_particle_output();

  cs_post_add_time_mesh_dep_output(_cs_lagr_post, &_lagr_post_options);
  _lagr_post_options_is_set = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate or deactive postprocessing for a given particle attribute.
 *
 * \param[in]  attr_id  associated attribute id
 *
 * \return     true if output of given attribute is active, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_lagr_post_get_attr(cs_lagr_attribute_t  attr_id)
{
  /* Initialize if not done yet */

  if (_lagr_post_options.attr_output[0] == -1) {
    for (cs_lagr_attribute_t i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
      _lagr_post_options.attr_output[i] = 0;
    }
  }

  bool retval = false;
  if (_lagr_post_options.attr_output[attr_id] > 0)
    retval = true;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate or deactive postprocessing for a given particle attribute.
 *
 * \param[in]  attr_id  associated attribute id
 * \param[in]  active   true if postprocessing is required, false otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_post_set_attr(cs_lagr_attribute_t  attr_id,
                      bool                 active)
{
  if (_lagr_post_options_is_set)
    bft_error(__FILE__, __LINE__, 0,
              _("%s should not be called after %s."),
              __func__, "cs_lagr_post_init");

  /* Initialize if not done yet */

  if (_lagr_post_options.attr_output[0] == -1) {
    for (cs_lagr_attribute_t i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
      _lagr_post_options.attr_output[i] = 0;
    }
  }

  cs_lagr_particle_attr_in_range(attr_id);

  if (active == false)
    _lagr_post_options.attr_output[attr_id] = 0;
  else
    _lagr_post_options.attr_output[attr_id] = 1;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
