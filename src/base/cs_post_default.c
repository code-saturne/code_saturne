/*============================================================================
 * Management of the post-processing
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_prototypes.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

  const cs_real_t  *dt;
  const cs_real_t  *rtpa;
  const cs_real_t  *rtp;
  const cs_real_t  *propce;
  const cs_real_t  *propfa;
  const cs_real_t  *propfb;
  const cs_real_t  *coefa;
  const cs_real_t  *coefb;
  const cs_real_t  *statce;
  const cs_real_t  *stativ;
  const cs_real_t  *statfb;

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
 *                   faces (ent_flag[1]), or boundary faces (ent_flag[2])
 *   n_cells     <-- local number of cells of post_mesh
 *   n_i_faces   <-- local number of interior faces of post_mesh
 *   n_b_faces   <-- local number of boundary faces of post_mesh
 *   cell_list   <-- list of cells (1 to n) of post-processing mesh
 *   i_face_list <-- list of interior faces (1 to n) of post-processing mesh
 *   b_face_list <-- list of boundary faces (1 to n) of post-processing mesh
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- current physical time
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

static void
_write_additional_vars(void             *input,
                       int               mesh_id,
                       int               cat_id,
                       int               ent_flag[3],
                       cs_lnum_t         n_cells,
                       cs_lnum_t         n_i_faces,
                       cs_lnum_t         n_b_faces,
                       const cs_lnum_t   cell_list[],
                       const cs_lnum_t   i_face_list[],
                       const cs_lnum_t   b_face_list[],
                       int               nt_cur_abs,
                       cs_real_t         t_cur_abs)
{
  /* Local variables */

  cs_post_default_input_t  *_input = input;

  int i;
  cs_int_t   itypps[3];
  cs_int_t   nummai = mesh_id;
  cs_int_t   numtyp = cat_id;

  cs_real_t  *var_trav = NULL;
  cs_real_t  *cel_vals = NULL;
  cs_real_t  *i_face_vals = NULL;
  cs_real_t  *b_face_vals = NULL;

  /* Basic initialization */

  for (i = 0; i < 3; i++)
    itypps[i] = ent_flag[i];

 /* Allocate work array to build variables */

  BFT_MALLOC(var_trav,
             (n_cells + n_i_faces + n_b_faces) * 3,
             cs_real_t);

  /* Pointers to variable assembly arrays, set to NULL if unused
     (so as to provoke an immediate error in case of incorrect use) */

  cel_vals = var_trav;
  i_face_vals = cel_vals + (n_cells * 3);
  b_face_vals = i_face_vals + (n_i_faces * 3);

  if (n_cells == 0)
    cel_vals = NULL;
  if (n_i_faces == 0)
    i_face_vals = NULL;
  if (n_b_faces == 0)
    b_face_vals = NULL;

  /* Add specific outputs for Code_Saturne */

  if (cat_id < 0)
    CS_PROCF(dvvpst, DVVPST) (&nummai, &numtyp,
                              _input->nvar, _input->nscal,
                              _input->nvlsta, _input->nvisbr,
                              &n_cells, &n_i_faces, &n_b_faces,
                              itypps,
                              cell_list, i_face_list, b_face_list,
                              _input->dt,
                              _input->rtpa, _input->rtp,
                              _input->propce, _input->propfa, _input->propfb,
                              _input->coefa, _input->coefb,
                              _input->statce, _input->stativ, _input->statfb,
                              cel_vals, i_face_vals, b_face_vals);

  /* Call to user subroutine for additional post-processing */

  CS_PROCF(usvpst, USVPST) (&nummai,
                            _input->nvar, _input->nscal, _input->nvlsta,
                            &n_cells, &n_i_faces, &n_b_faces,
                            itypps,
                            cell_list, i_face_list, b_face_list,
                            _input->dt,
                            _input->rtpa, _input->rtp,
                            _input->propce, _input->propfa, _input->propfb,
                            _input->statce,
                            cel_vals, i_face_vals, b_face_vals);

  /* Free work array */

  BFT_FREE(var_trav);
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * Fortran interface:
 *
 * subroutine pstvar
 * *****************
 *                  ( ntcabs,
 *                    nvar,   nscal,  nvlsta, nvisbr,
 *                    ttcabs,
 *                    dt,     rtpa,   rtp,    propce, propfa, propfb,
 *                    coefa,  coefb,
 *                    statce, stativ, statfb )
 *
 * integer          ntcabs      : --> : current time step number
 * integer          nvar        : <-- : number of variables
 * integer          nscal       : <-- : number of scalars
 * integer          nvlsta      : <-- : number of statistical variables (lagr)
 * integer          nvisbr      : <-- : number of boundary stat. variables (lagr)
 * double precision ttcabs      : <-- : current physical time
 * double precision dt          : <-- : local time step
 * double precision rtpa        : <-- : cell variables at previous time step
 * double precision rtp         : <-- : cell variables
 * double precision propce      : <-- : cell physical properties
 * double precision propfa      : <-- : interior face physical properties
 * double precision propfb      : <-- : boundary face physical properties
 * double precision coefa       : <-- : boundary conditions array
 * double precision coefb       : <-- : boundary conditions array
 * double precision statce      : <-- : cell statistics (lagrangian)
 * double precision stativ      : <-- : cell variance statistics (lagrangian)
 * double precision statfb      : <-- : boundary face statistics (lagrangian)
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *ntcabs,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr,
 const cs_real_t  *ttcabs,
 const cs_real_t   dt[],
 const cs_real_t   rtpa[],
 const cs_real_t   rtp[],
 const cs_real_t   propce[],
 const cs_real_t   propfa[],
 const cs_real_t   propfb[],
 const cs_real_t   coefa[],
 const cs_real_t   coefb[],
 const cs_real_t   statce[],
 const cs_real_t   stativ[],
 const cs_real_t   statfb[]
)
{
  /* Define or update map of variables */

  _default_input.nvar = nvar;
  _default_input.nscal = nscal;

  _default_input.nvlsta = nvlsta;
  _default_input.nvisbr = nvisbr;

  _default_input.dt = dt;
  _default_input.rtpa = rtpa;
  _default_input.rtp = rtp;
  _default_input.propce = propce;
  _default_input.propfa = propfa;
  _default_input.propfb = propfb;
  _default_input.coefa = coefa;
  _default_input.coefb = coefb;
  _default_input.statce = statce;
  _default_input.stativ = stativ;
  _default_input.statfb = statfb;

  /* Register function for first pass */

  if (_default_input_is_set == false) {

    cs_post_add_time_mesh_dep_output(_write_additional_vars,
                                     &_default_input);

    _default_input_is_set = true;
  }

  /* Call main post-processing function */

  cs_post_write_vars(*ntcabs, *ttcabs);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
