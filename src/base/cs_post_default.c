/*============================================================================
 * Management of the post-processing
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
#include "cs_timer_stats.h"

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
 *                   particles (ent_flag[3]) or probes (ent_flag[4])
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
                       int                    ent_flag[5],
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
                              _input->nvar,
                              &n_cells, &n_b_faces,
                              cell_list, b_face_list,
                              cel_vals, b_face_vals);

  /* Free work array */

  BFT_FREE(var_trav);

  /* Call to user subroutine for additional post-processing */

  CS_PROCF(usvpst, USVPST) (&nummai,
                            _input->nvar, _input->nscal, 0,
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
  int t_top_id
    = cs_timer_stats_switch(cs_timer_stats_id_by_name("postprocessing_stage"));

  cs_post_write_meshes(cs_glob_time_step);

  cs_timer_stats_switch(t_top_id);
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
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *ntcabs,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal
)
{
  /* Define or update map of variables */

  _default_input.nvar = nvar;
  _default_input.nscal = nscal;

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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
