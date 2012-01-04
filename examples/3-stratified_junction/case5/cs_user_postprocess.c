
/* Code_Saturne version 2.1.0-alpha1 */

/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2012 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Define (conforming or non-conforming) mesh joinings.
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvm_writer.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
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
   *   'discard_polygons'  (ignore polygon-type faces)
   *   'discard_polyhedra' (ignore polyhedron-type cells)
   *   'divide_polygons'   (subdivides polygon-type faces)
   *   'divide_polyhedra'  (subdivides polyhedron-type cells)
   *   'split_tensors'     (writes tensors as separate scalars) */

  /* Default writer time dependency */

  fvm_writer_time_dep_t   time_dep = FVM_WRITER_FIXED_MESH;

  /* Default time step or physical time based output frequencies */

  cs_bool_t  output_at_end = true;
  int        ntchr = -1;
  double     frchr = -1.0;

  /* Default output format and options */

  const char format_name[] = "EnSight Gold";
  const char format_options[] = "";

  /* Define additional writers */
  /* ------------------------- */

  ntchr = 5;

  cs_post_define_writer(1,                            /* writer_id */
                        "tinf21",                     /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format name */
                        "discard_polygons, discard_polyhedra",
                        FVM_WRITER_TRANSIENT_CONNECT,
                        output_at_end,
                        ntchr,
                        frchr);
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
  /* Post-processing meshes may be defined using one of several functions,
   * whose protypes are defined in cs_post.h; these functions are:
   *
   * Functions simplest to use are cs_post_define_volume_mesh() and
   * cs_post_define_surface_mesh(), which allow defining volume or surface
   * post-processing meshes using selection criteria.
   *
   * parameters for cs_post_define_volume_mesh():
   *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   mesh_name      <-- associated mesh name
   *   cell_criteria  <-- selection criteria for cells
   *   add_groups     <-- if true, add group information if present
   *   auto_variables <-- if true, automatic output of main variables
   *   n_writers      <-- number of associated writers
   *   writer_ids     <-- ids of associated writers
   *
   * parameters for cs_post_define_surface_mesh():
   *   mesh_id         <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   mesh_name       <-- associated mesh name
   *   i_face_criteria <-- selection criteria for interior faces
   *   b_face_criteria <-- selection criteria for boundary faces
   *   add_groups      <-- if true, add group information if present
   *   auto_variables  <-- if true, automatic output of main variables
   *   n_writers       <-- number of associated writers
   *   writer_ids      <-- ids of associated writers
   *
   * If no writer is associated to a mesh, it is not output, and its
   * construction may be avoided altogether (at least when defined
   * by one of the above functions).
   *
   * More advanced functions are described along with examples below. */

  /*--------------------------------------------------------------------------*/

  /* Select interior faces with y = 0. */

  if (true) {

    const int n_writers = 2;
    const int writer_ids[] = {-1};  /* Associate to default writer */

    /* Select cells with y = 0 */
    const char *interior_criteria = "plane[0, -1, 0, 0.0, "
                                    "epsilon = 0.0001]";
    const char *boundary_criteria = NULL;

    cs_post_define_surface_mesh(1,               /* mesh id */
                                "Cut 1",
                                interior_criteria,
                                boundary_criteria,
                                false, /* add_groups */
                                false, /* auto_variables */
                                n_writers,
                                writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  /* Select no cells, will select cells with T < 21 degrees in 'usmpst' */

  if (true) {

    const int n_writers = 1;
    const int writer_ids[] = {1};  /* Associate to writer 1 */
  
    cs_post_define_volume_mesh(2,                 /* mesh id */
                               "celTinf21",
                               NULL,
                               false,             /* add_groups */
                               true,              /* auto_variables */
                               n_writers,
                               writer_ids);
  
  }

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
  /* Use the cs_post_activate_writer() function to force the
   * "active" or "inactive" flag for a specific writer or for all
   * writers for the current time step.

   * the parameters for cs_post_activate_writer() are:
   *   writer_id <-- writer id, or 0 for all writers
   *   activate  <-- false to deactivate, true to activate */

  if (false) { /* example: deactivate all output before time step 1000 */

    if (nt_max_abs < 1000) {
      int writer_id = 0; /* 0: all writers */
      cs_post_activate_writer(writer_id, false);
    }

  }
}

/*============================================================================
 * Fortran-callable wrapper for user function definitions (do not remove).
 *============================================================================*/

/*----------------------------------------------------------------------------
 * User override of default frequency or calculation end based output.
 *
 * Fortran interface:
 *
 * subroutine pstusn (ntmabs, ntcabs, ttcabs)
 * *****************
 *
 * integer          ntmabs      : <-- : maximum time step number
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : absolute time at the current time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstusn, PSTUSN)
(
 const cs_int_t  *ntmabs,
 const cs_int_t  *ntcabs,
 const cs_real_t *ttcabs
)
{
  cs_user_postprocess_activate(*ntmabs, *ntcabs, *ttcabs);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
