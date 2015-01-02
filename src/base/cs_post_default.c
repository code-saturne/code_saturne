/*============================================================================
 * Management of the post-processing
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_lagr_tracking.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_time_step.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local types and structures
 *============================================================================*/

/* Structure used to pass Fortran array pointer arguments */
/*--------------------------------------------------------*/

typedef struct {

  const cs_int_t   *nvar;
  const cs_int_t   *nscal;
  const cs_int_t   *nvlsta;
  const cs_int_t   *nvisbr;

  /* Lagrangian variables */

  bool      particle_attr[CS_LAGR_N_ATTRIBUTES];
  cs_int_t  particle_multicomponent_export[CS_LAGR_N_ATTRIBUTES];

} cs_post_default_input_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default output format and options */

static cs_post_default_input_t  _default_input;
static bool                     _default_input_is_set = false;

/*============================================================================
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
_write_particle_vars(cs_post_default_input_t  *input,
                     int                       mesh_id,
                     const cs_time_step_t     *ts)
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
 * Default additional output of mesh and time-dependent variables for the
 * call to pstvar / cs_post_write_vars.
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
 *   ent_flag    <-- indicate global presence of cells (ent_flag[0]), interior
 *                   faces (ent_flag[1]), boundary faces (ent_flag[2]),
 *                   or particles (ent_flag[3])
 *   n_cells     <-- local number of cells of post_mesh
 *   n_i_faces   <-- local number of interior faces of post_mesh
 *   n_b_faces   <-- local number of boundary faces of post_mesh
 *   cell_list   <-- list of cells (1 to n) of post-processing mesh
 *   i_face_list <-- list of interior faces (1 to n) of post-processing mesh
 *   b_face_list <-- list of boundary faces (1 to n) of post-processing mesh
 *   ts          <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_write_additional_vars(void                  *input,
                       int                    mesh_id,
                       int                    cat_id,
                       int                    ent_flag[4],
                       cs_lnum_t              n_cells,
                       cs_lnum_t              n_i_faces,
                       cs_lnum_t              n_b_faces,
                       const cs_lnum_t        cell_list[],
                       const cs_lnum_t        i_face_list[],
                       const cs_lnum_t        b_face_list[],
                       const cs_time_step_t  *ts)
{
  /* Local variables */

  cs_post_default_input_t  *_input = input;

  int i;
  cs_int_t   itypps[4];
  cs_int_t   nummai = mesh_id;
  cs_int_t   numtyp = cat_id;

  cs_real_t  *var_trav = NULL;
  cs_real_t  *cel_vals = NULL;
  cs_real_t  *b_face_vals = NULL;

  /* Specific handling for particle meshes */

  if (cat_id == -3) {
    _write_particle_vars(_input, mesh_id, ts);
    return;
  }

  /* Basic initialization */

  for (i = 0; i < 4; i++)
    itypps[i] = ent_flag[i];

 /* Allocate work array to build variables */

  BFT_MALLOC(var_trav,
             (n_cells + n_i_faces + n_b_faces) * 3,
             cs_real_t);

  /* Pointers to variable assembly arrays, set to NULL if unused
     (so as to provoke an immediate error in case of incorrect use) */

  cel_vals = var_trav;
  b_face_vals = cel_vals + (n_cells * 3);

  if (n_cells == 0)
    cel_vals = NULL;
  if (n_b_faces == 0)
    b_face_vals = NULL;

  /* Add specific outputs for Code_Saturne */

  if (cat_id < 0)
    CS_PROCF(dvvpst, DVVPST) (&nummai, &numtyp,
                              _input->nvar, _input->nscal,
                              _input->nvlsta, _input->nvisbr,
                              &n_cells, &n_b_faces,
                              cell_list, b_face_list,
                              cel_vals, b_face_vals);

  /* Free work array */

  BFT_FREE(var_trav);

  /* Call to user subroutine for additional post-processing */

  CS_PROCF(usvpst, USVPST) (&nummai,
                            _input->nvar, _input->nscal, _input->nvlsta,
                            &n_cells, &n_i_faces, &n_b_faces,
                            itypps,
                            cell_list, i_face_list, b_face_list);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * Fortran interface:
 *
 * subroutine pstgeo
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstgeo, PSTGEO)
(
 void
)
{
  cs_post_write_meshes(cs_glob_time_step);
}

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * Fortran interface:
 *
 * subroutine pstvar
 * *****************
 *                  ( ntcabs,
 *                    nvar,   nscal,  nvlsta, nvisbr )
 *
 * integer          ntcabs      : --> : current time step number
 * integer          nvar        : <-- : number of variables
 * integer          nscal       : <-- : number of scalars
 * integer          nvlsta      : <-- : number of statistical variables (lagr)
 * integer          nvisbr      : <-- : number of boundary stat. variables (lagr)
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *ntcabs,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr
)
{
  /* Define or update map of variables */

  _default_input.nvar = nvar;
  _default_input.nscal = nscal;

  _default_input.nvlsta = nvlsta;
  _default_input.nvisbr = nvisbr;

  /* Register function for first pass */

  if (_default_input_is_set == false) {

    cs_post_add_time_mesh_dep_output(_write_additional_vars,
                                     &_default_input);

    _default_input_is_set = true;
  }

  /* Call main post-processing function */

  if (*ntcabs > -1)
    cs_post_write_vars(cs_glob_time_step);
  else
    cs_post_write_vars(NULL);
}

/*----------------------------------------------------------------------------
 * Define which Lagragian variables should be postprocessed
 *
 * Fortran interface:
 *
 * subroutine lagpvr
 * *****************
 *                  ( ivisv1, ivisv2, ivistp,  ivisdm, iviste,
 *                    ivismp, ivisdk, iviswat, ivisch, ivisck )
 *
 * integer          ivisv1      : <-- : display of variable 'fluid velocity'
 * integer          ivisv2      : <-- : display of variable 'particles velocity'
 * integer          ivistp      : <-- : display of variable 'resident time'
 * integer          ivisdm      : <-- : display of variable 'particle diameter'
 * integer          iviste      : <-- : display of variable 'particle temperature'
 * integer          ivismp      : <-- : display of variable 'particle mass'
 * integer          ivisdk      : <-- : display of variable 'core diameter of part.'
 * integer          iviswat     : <-- : display of variable 'mass of water in coal'
 * integer          ivisch      : <-- : display of variable 'mass of reactive coal'
 * integer          ivisck      : <-- : display of variable 'mass of char'
 *----------------------------------------------------------------------------*/

void CS_PROCF (lagpvr, LAGPVR)
(
 const cs_int_t  *ivisv1,
 const cs_int_t  *ivisv2,
 const cs_int_t  *ivistp,
 const cs_int_t  *ivisdm,
 const cs_int_t  *iviste,
 const cs_int_t  *ivismp,
 const cs_int_t  *ivisdk,
 const cs_int_t  *iviswat,
 const cs_int_t  *ivisch,
 const cs_int_t  *ivisck
)
{
  cs_lagr_attribute_t attr_id;

  for (attr_id = 0; attr_id < CS_LAGR_N_ATTRIBUTES; attr_id++) {
    _default_input.particle_attr[attr_id] = false;
    _default_input.particle_multicomponent_export[attr_id] = -1;
  }

  if (*ivisv1)
    _default_input.particle_attr[CS_LAGR_VELOCITY] = true;

  if (*ivisv2)
    _default_input.particle_attr[CS_LAGR_VELOCITY_SEEN] = true;

  if (*ivistp)
    _default_input.particle_attr[CS_LAGR_RESIDENCE_TIME] = true;

  if (*ivisdm)
    _default_input.particle_attr[CS_LAGR_DIAMETER] = true;

  if (*iviste) {
    _default_input.particle_attr[CS_LAGR_TEMPERATURE] = true;
    if (cs_glob_lagr_params->n_temperature_layers > 1)
      _default_input.particle_multicomponent_export[CS_LAGR_TEMPERATURE]
        = cs_glob_lagr_params->n_temperature_layers;
  }

  if (*ivismp)
    _default_input.particle_attr[CS_LAGR_MASS] = true;

  if (*ivisdk)
    _default_input.particle_attr[CS_LAGR_SHRINKING_DIAMETER] = true;

  if (*iviswat)
    _default_input.particle_attr[CS_LAGR_WATER_MASS] = true;

  if (*ivisch) {
    _default_input.particle_attr[CS_LAGR_COAL_MASS] = true;
    if (cs_glob_lagr_params->n_temperature_layers > 1)
      _default_input.particle_multicomponent_export[CS_LAGR_COAL_MASS]
        = cs_glob_lagr_params->n_temperature_layers;
  }

  if (*ivisck) {
    _default_input.particle_attr[CS_LAGR_COKE_MASS] = true;
    if (cs_glob_lagr_params->n_temperature_layers > 1)
      _default_input.particle_multicomponent_export[CS_LAGR_COKE_MASS]
        = cs_glob_lagr_params->n_temperature_layers;
  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
