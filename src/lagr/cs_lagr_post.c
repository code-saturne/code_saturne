/*============================================================================
 * Lagrangian module postprocessing
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "cs_prototypes.h"

#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_parall.h"
#include "cs_post.h"

#include "cs_lagr.h"
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

/* Structure used to pass postprocessing options */
/*-----------------------------------------------*/

typedef struct {

  bool      particle_attr[CS_LAGR_N_ATTRIBUTES];
  cs_int_t  particle_multicomponent_export[CS_LAGR_N_ATTRIBUTES];

} cs_post_lagr_input_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default output format and options */

static cs_post_lagr_input_t  _lagr_input;
static bool                  _lagr_input_is_set = false;


/* Postprocessing options structure and associated pointer */

static cs_lagr_post_options_t  _lagr_post_options
= {.iensi3 = 0,
   .ivisv1 = 0,
   .ivisv2 = 0,
   .ivistp = 0,
   .ivisdm = 0,
   .iviste = 0,
   .ivismp = 0,
   .ivisdk = 0,
   .ivisch = 0,
   .ivisck = 0,
   .iviswat = 0};

const cs_lagr_post_options_t *cs_glob_lagr_post_options = &_lagr_post_options;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

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
 *                   here, we should point to _default_input.
 *   mesh_id     <-- id of the output mesh for the current call
 *   cat_id      <-- category id of the output mesh for the current call
 *   ts          <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_write_particle_vars(cs_post_lagr_input_t   *input,
                     int                     mesh_id,
                     const cs_time_step_t   *ts)
{
  cs_lagr_attribute_t attr_id;

  char var_name[64];
  int  component_id;
  char var_name_component[64];

  for (attr_id = 0; attr_id < CS_LAGR_N_ATTRIBUTES; attr_id++) {

    if (input->particle_attr[attr_id]) {

      /* build name */

      int i;
      int l = snprintf(var_name,
                       63,
                       "particle_%s",
                       cs_lagr_attribute_name[attr_id] + strlen("cs_lagr_"));
      var_name[63] = '\0';
      for (i = 0; i < l; i++)
        var_name[i] = tolower(var_name[i]);

      /* Output values */

      if (input->particle_multicomponent_export[attr_id] == -1)
        cs_post_write_particle_values(mesh_id,
                                      attr_id,
                                      var_name,
                                      input->particle_multicomponent_export[attr_id],
                                      ts);
      else {
        /* Create one output per component */
        for (component_id = 0;
             component_id < input->particle_multicomponent_export[attr_id];
             component_id++) {
          snprintf(var_name_component,
                   63,
                   "%s_layer_%2.2i",
                   var_name,
                   component_id+1);
          var_name_component[63] = '\0';
          cs_post_write_particle_values(mesh_id,
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
 * Define which Lagragian variables should be postprocessed
 *----------------------------------------------------------------------------*/

static void
_activate_particle_output(void)
{
  cs_lagr_attribute_t attr_id;

  for (attr_id = 0; attr_id < CS_LAGR_N_ATTRIBUTES; attr_id++) {
    _lagr_input.particle_attr[attr_id] = false;
    _lagr_input.particle_multicomponent_export[attr_id] = -1;
  }

  if (cs_glob_lagr_post_options->ivisv1)
    _lagr_input.particle_attr[CS_LAGR_VELOCITY] = true;

  if (cs_glob_lagr_post_options->ivisv2)
    _lagr_input.particle_attr[CS_LAGR_VELOCITY_SEEN] = true;

  if (cs_glob_lagr_post_options->ivistp)
    _lagr_input.particle_attr[CS_LAGR_RESIDENCE_TIME] = true;

  if (cs_glob_lagr_post_options->ivisdm)
    _lagr_input.particle_attr[CS_LAGR_DIAMETER] = true;

  if (cs_glob_lagr_post_options->iviste) {
    _lagr_input.particle_attr[CS_LAGR_TEMPERATURE] = true;
    if (cs_glob_lagr_model->n_temperature_layers > 1)
      _lagr_input.particle_multicomponent_export[CS_LAGR_TEMPERATURE]
        = cs_glob_lagr_model->n_temperature_layers;
  }

  if (cs_glob_lagr_post_options->ivismp)
    _lagr_input.particle_attr[CS_LAGR_MASS] = true;

  if (cs_glob_lagr_post_options->ivisdk)
    _lagr_input.particle_attr[CS_LAGR_SHRINKING_DIAMETER] = true;

  if (cs_glob_lagr_post_options->iviswat)
    _lagr_input.particle_attr[CS_LAGR_WATER_MASS] = true;

  if (cs_glob_lagr_post_options->ivisch) {
    _lagr_input.particle_attr[CS_LAGR_COAL_MASS] = true;
    if (cs_glob_lagr_model->n_temperature_layers > 1)
      _lagr_input.particle_multicomponent_export[CS_LAGR_COAL_MASS]
        = cs_glob_lagr_model->n_temperature_layers;
  }

  if (cs_glob_lagr_post_options->ivisck) {
    _lagr_input.particle_attr[CS_LAGR_COKE_MASS] = true;
    if (cs_glob_lagr_model->n_temperature_layers > 1)
      _lagr_input.particle_multicomponent_export[CS_LAGR_COKE_MASS]
        = cs_glob_lagr_model->n_temperature_layers;
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
 *   cell_list   <-- list of cells (1 to n) of post-processing mesh
 *   i_face_list <-- list of interior faces (1 to n) of post-processing mesh
 *   b_face_list <-- list of boundary faces (1 to n) of post-processing mesh
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
              const cs_lnum_t        cell_list[],
              const cs_lnum_t        i_face_list[],
              const cs_lnum_t        b_face_list[],
              const cs_time_step_t  *ts)
{
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(cell_list);
  CS_UNUSED(i_face_list);

  /* Specific handling for particle meshes */

  if (cat_id == -3) {
    _write_particle_vars(input, mesh_id, ts);
    return;
  }

  /* Boundary statistics */
  /*---------------------*/

  else if (   cat_id == -2
           && cs_glob_lagr_time_scheme->iilagr > 0
           && cs_glob_lagr_post_options->iensi3 > 0) {

    cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;

    const cs_lagr_boundary_interactions_t
      *lagr_b = cs_glob_lagr_boundary_interactions;

    int  nvisbr = cs_glob_lagr_dim->nvisbr;
    cs_real_t  seuilf = cs_glob_lagr_stat_options->threshold;

    cs_real_t *val;
    BFT_MALLOC(val, nfabor, cs_real_t);

    for (int irf = 0; irf < nvisbr; irf++) {

      const char *var_name = lagr_b->nombrd[irf];

      const cs_real_t *_b_stats = bound_stat + nfabor*irf;

      if (lagr_b->imoybr[irf] > 1) {

        const cs_real_t *_f_count = NULL;

        if (lagr_b->imoybr[irf] == 3)
          _f_count = bound_stat + nfabor*lagr_b->iencnb;

        else if (lagr_b->imoybr[irf] == 2)
          _f_count = bound_stat + nfabor*lagr_b->inbr;

        for (cs_lnum_t i = 0; i < n_b_faces; i++) {
          cs_lnum_t f_id = b_face_list[i] - 1;
          if (_f_count[f_id] > seuilf)
            val[f_id] = _b_stats[f_id] / _f_count[f_id];
          else
            val[f_id] = 0.;
        }

      }
      else if (lagr_b->imoybr[irf] == 1) {

        cs_real_t  tstatp = lagr_b->tstatp;

        for (cs_lnum_t i = 0; i < n_b_faces; i++) {
          cs_lnum_t f_id = b_face_list[i] - 1;
          val[f_id] = _b_stats[f_id] / tstatp;
        }

      }

      else {

        const cs_real_t *_f_count = bound_stat + nfabor*lagr_b->inbr;

        for (cs_lnum_t i = 0; i < n_b_faces; i++) {
          cs_lnum_t f_id = b_face_list[i] - 1;
          if (_f_count[f_id] > seuilf)
            val[f_id] = _b_stats[f_id];
          else
            val[f_id] = 0.;
        }

      }

      cs_post_write_var(mesh_id,
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

  /* Add boundary zone ids */

  if (cat_id == -2) {
    cs_lagr_bdy_condition_t  *bdy_cond = cs_glob_lagr_bdy_conditions;
    cs_post_write_var(mesh_id,
                      "lagrangian_boundary_zones",
                      1,       /* var_dim */
                      true,    /* interlace */
                      true,    /* use_parent */
                      CS_POST_TYPE_int,
                      NULL,
                      NULL,
                      bdy_cond->b_face_zone_id,
                      cs_glob_time_step);

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

  cs_post_add_time_mesh_dep_output(_cs_lagr_post, &_lagr_input);
  _lagr_input_is_set = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the global cs_lagr_post_options_t structure.
 *
 * This pointer allows write access to the structure.
 *
 * \return pointer to cs_glob_lagr_post_options
 */
/*----------------------------------------------------------------------------*/

cs_lagr_post_options_t *
cs_lagr_post_get_options(void)
{
  return &_lagr_post_options;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
