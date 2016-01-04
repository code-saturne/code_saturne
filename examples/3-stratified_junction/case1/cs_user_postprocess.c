/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_mesh.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Example function for selection of cells with a temperature below 21
 * degrees.
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_cells  --> number of selected cells
 *   cell_ids --> array of selected cell ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_t_lt_21_select(void        *input,
                cs_lnum_t   *n_cells,
                cs_lnum_t  **cell_ids)
{
  cs_lnum_t i;

  cs_lnum_t _n_cells = 0;
  cs_lnum_t *_cell_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  cs_field_t *f = cs_field_by_name("temperature"); /* Get access to field */

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "No field with name \"temperature\" defined");

  /* Before time loop, field is defined, but has no values yet,
     so ignore that case (postprocessing mesh will be initially empty) */

  if (f->val != NULL) {

    BFT_MALLOC(_cell_ids, m->n_cells, cs_lnum_t); /* Allocate selection list */

    for (i = 0; i < m->n_cells; i++) {
      if (f->val[i] < 21) {
        _cell_ids[_n_cells] = i;
        _n_cells += 1;
      }
    }

    BFT_REALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Adjust size (good practice,
                                                    but not required) */

  }

  /* Set return values */

  *n_cells = _n_cells;
  *cell_ids = _cell_ids;
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void)
{
  /* Every writer has a a strictly positive or negative id. Negative ids
   * are for predefined writers, positive ids for user writers.
   * All predefined writers use the settings from writer -1, and
   * redefining that writer here allows changing from the default or GUI
   * settings.
   *
   * Defining or configuring a writer is done by calling the
   * cs_post_define_writer() function, whose arguments are:
   *   writer_id     <-- number of writer to create (< 0 reserved, > 0 for user)
   *   case_name     <-- associated case name
   *   dir_name      <-- associated directory name
   *   fmt_name      <-- associated format name
   *   fmt_opts      <-- associated format options string
   *   time_dep      <-- FVM_WRITER_FIXED_MESH if mesh definitions are fixed,
   *                     FVM_WRITER_TRANSIENT_COORDS if coordinates change,
   *                     FVM_WRITER_TRANSIENT_CONNECT if connectivity changes
   *   output_at_end <-- force output at calculation end if not 0
   *   frequency_n   <-- default output frequency in time-steps, or < 0
   *   frequency_t   <-- default output frequency in seconds, or < 0
   *                     (has priority over frequency_n)
   *
   * Allowed output format names: "EnSight Gold", "MED", or "CGNS".
   * (EnSight output is built-in; MED or CGNS are only available if the
   * code was built with these optional libraries)
   *
   * An output options string may contain options (separated by whitespace
   * or commas) from the following list:
   *   'text'              (text format, for EnSight)
   *   'big_endian'        (forces binary EnSight output to 'big-endian' mode)
   *   'adf'               (use ADF file type, for CGNS)
   *   'hdf5'              (force HDF5 file type, usual the default for CGNS)
   *   'discard_polygons'  (ignore polygon-type faces)
   *   'discard_polyhedra' (ignore polyhedron-type cells)
   *   'divide_polygons'   (subdivides polygon-type faces)
   *   'divide_polyhedra'  (subdivides polyhedron-type cells)
   *   'split_tensors'     (writes tensors as separate scalars) */

  /* Define additional writers */
  /* ------------------------- */

  cs_post_define_writer(1,                            /* writer_id */
                        "user",                       /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format name */
                        "",                           /* format options */
                        FVM_WRITER_TRANSIENT_CONNECT, /* time dependency */
                        true,                         /* output at end */
                        5,                            /* time step frequency */
                        -1);                          /* Time value frequency */
}

/*----------------------------------------------------------------------------
 * Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void)
{
  /* Advanced volume mesh element selection is possible using
   * cs_post_define_volume_mesh_by_func(), which allows defining
   * meshes using user-defined element lists.
   *
   * parameters for cs_post_define_volume_mesh_by_func():
   *   mesh_id           <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   mesh_name         <-- associated mesh name
   *   cell_select_func  <-- pointer to cells selection function
   *   cell_select_input <-> pointer to optional input data for the cell
   *                         selection function, or NULL
   *   time_varying      <-- if true, try to redefine mesh at each output time
   *   add_groups        <-- if true, add group information if present
   *   auto_variables    <-- if true, automatic output of main variables
   *   n_writers         <-- number of associated writers
   *   writer_ids          <-- ids of associated writers */

  /* Build a (time varying) volume mesh containing cells
     with values of field named "temperature" > < 21 */

  const int n_writers = 1;
  const int writer_ids[] = {1};  /* Associate to writer 1 */

  /* Define postprocessing mesh */

  cs_post_define_volume_mesh_by_func(1,               /* mesh id */
                                     "T_lt_21",
                                     _t_lt_21_select,
                                     NULL,            /* _t_lt_21_select input */
                                     true,            /* time varying */
                                     false,           /* add_groups */
                                     false,           /* auto_variables */
                                     n_writers,
                                     writer_ids);
}

/*----------------------------------------------------------------------------
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * parameters:
 *   nt_max_abs <-- maximum time step number
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
