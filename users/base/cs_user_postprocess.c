/*============================================================================
 * Define (conforming or non-conforming) mesh joinings.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Define default writer */
  /* --------------------- */

  if (false)
    cs_post_define_writer(-1,                /* writer_id */
                          "results",         /* writer name */
                          "postprocessing",  /* directory name */
                          format_name,
                          format_options,
                          time_dep,
                          output_at_end,
                          ntchr,
                          frchr);

  /* Define additional writers */
  /* ------------------------- */

  if (false)
    cs_post_define_writer(1,                            /* writer_id */
                          "user",                       /* writer name */
                          "postprocessing",             /* directory name */
                          "EnSight Gold",               /* format name */
                          "binary, discard_polygons, discard_polyhedra",
                          time_dep,
                          output_at_end,
                          ntchr,
                          frchr);

  if (false)
    cs_post_define_writer(2,                            /* writer_id */
                          "user_txt",                   /* writer name */
                          "postprocessing",             /* directory name */
                          "ensight",                    /* format name */
                          "text, divide_polyhedra", 
                          time_dep,                     /* modification flag */
                          output_at_end,
                          ntchr,
                          frchr);

  if (false)
    cs_post_define_writer(3,                            /* writer_id */
                          "modif",                      /* writer name */
                          "postprocessing",             /* directory name */
                          "ensight",                    /* format name */
                          "discard_polyhedra", 
                          FVM_WRITER_TRANSIENT_CONNECT,
                          false,
                          3,
                          frchr);       

  if (false)
    cs_post_define_writer(4,                            /* writer_id */
                          "exchange",                   /* writer name */
                          "postprocessing",             /* directory name */
                          "MED",                        /* format name  */
                          "",                           /* format options */
                          FVM_WRITER_TRANSIENT_COORDS,  /* modification flag */
                          false,
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
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */
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

  /* Reconfigure predefined meshes (mesh_id -1 for volume, -2 for boundary */

  if (false) {

    /* De-activate boundary mesh output by redefining it with no writer
       association (default is:
       int n_writers = 1;
       const int writer_ids[] = {-1});
    */

    int n_writers = 0;
    const int *writer_ids = NULL;

    cs_post_define_surface_mesh(-2,          /* mesh_id of main boundary mesh */
                                "Boundary",  /* mesh name */
                                NULL,        /* interior face selection criteria */
                                "all[]",     /* boundary face selection criteria */
                                true,        /* add_groups */
                                true,        /* automatic variables output */
                                n_writers,
                                writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  if (false) {

    /* Example: select interior faces with y = 0.5 */

    const int n_writers = 2;
    const int writer_ids[] = {1, 4};  /* Associate to writers 1 and 4 */

    const char *interior_criteria = "plane[0, -1, 0, 0.5, "
                                    "epsilon = 0.0001]";
    const char *boundary_criteria = NULL;

    cs_post_define_surface_mesh(1,               /* mesh id */
                                "Median plane",
                                interior_criteria,
                                boundary_criteria,
                                false, /* add_groups */
                                false, /* auto_variables */
                                n_writers,
                                writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  if (false) {

    /* Example: select all cells, and will modify the selection in 'usmpst' */

    const int n_writers = 1;
    const int writer_ids[] = {3};  /* Associate to writer 3 */

    cs_post_define_volume_mesh(2,                 /* mesh id */
                               "Volume v > 0.5",
                               "all[]",
                               false,             /* add_groups */
                               true,              /* auto_variables */
                               n_writers,
                               writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  if (false) {

    /* Example: select all boundary faces, and will modify the selection
       in 'usmpst'. */

    const int n_writers = 1;
    const int writer_ids[] = {3};  /* Associate to writer 3 */

    const char *interior_criteria = NULL;
    const char *boundary_criteria = "all[]";

    cs_post_define_surface_mesh(3,                /* mesh id */
                                "Surface iso v",
                                interior_criteria,
                                boundary_criteria,
                                false,            /* add_groups */
                                true,             /* auto_variables */
                                n_writers,
                                writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  /* The same variables will be output through all writers associated
   * to a mesh. In cases where different variables of a same mesh should
   * be output throught different writers, the solution is to define one or
   * several "aliases" of that mesh, allowing to assign a different id,
   * writers, and variables to each secondary copy of the mesh, without the
   * overhead of a full copy. The cs_post_define_alias_mesh() function
   * may be used for such a purpose.
   *
   * parameters for cs_post_define_alias_mesh():
   *   mesh_id         <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   aliased_mesh_id <-- id of aliased mesh
   *   auto_variables  <-- if true, automatic output of main variables
   *   n_writers       <-- number of associated writers
   *   writer_ids      <-- ids of associated writers  */

  if (false) {

    /* Example: define an alias of the main surface mesh */

    const int n_writers = 1;
    const int writer_ids[] = {4};  /* Associate to writer 4 */

    cs_post_define_alias_mesh(4,                 /* mesh id */
                              -2,                /* aliased mesh id */
                              false,             /* auto_variables */
                              n_writers,
                              writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  /* More advanced mesh element selection is possible using
   * cs_post_define_volume_mesh_by_list() or
   * cs_post_define_surface_mesh_by_list(), which allow defining
   * volume or surface meshes using user-defined element lists.
   *
   * parameters for cs_post_define_volume_mesh_by_list():
   *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   mesh_name      <-- associated mesh name
   *   n_cells        <-- number of selected cells
   *   cell_list      <-> list of selected cells (1 to n numbering)
   *   add_groups     <-- if true, add group information if present
   *   auto_variables <-- if true, automatic output of main variables
   *   n_writers      <-- number of associated writers
   *   writer_ids     <-- ids of associated writers
   *
   * parameters for cs_post_define_surface_mesh_by_list():
   *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   mesh_name      <-- associated mesh name
   *   n_i_faces      <-- number of associated interior faces
   *   n_b_faces      <-- number of associated boundary faces
   *   i_face_list    <-> list of associated interior faces (1 to n numbering)
   *   b_face_list    <-> list of associated boundary faces (1 to n numbering)
   *   add_groups     <-- if true, add group information if present
   *   auto_variables <-- if true, automatic output of main variables
   *   n_writers      <-- number of associated writers
   *   writer_ids     <-- ids of associated writers */

  if (false) {

    /* Advanced example:
       Build a surface mesh containing interior faces separating cells of group "2"
       from those of group "3", (assuming no cell has both colors), as well as
       boundary faces of group "4". */

    fvm_lnum_t n_i_faces = 0, n_b_faces = 0;
    fvm_lnum_t *i_face_list = NULL, *b_face_list = NULL;

    const int n_writers = 1;
    const int writer_ids[] = {1};  /* Associate to writer 1 */

    const cs_mesh_t *m = cs_glob_mesh;

    /* Allocate selection lists */

    BFT_MALLOC(i_face_list, m->n_i_faces, fvm_lnum_t);
    BFT_MALLOC(b_face_list, m->n_b_faces, fvm_lnum_t);

    /* Select interior faces */

    {
      fvm_lnum_t i, face_id;
      fvm_lnum_t n_families = 0;
      cs_int_t *family_list = NULL;
      int *family_mask = NULL;

      /* Build mask on families matching groups "2" (1), "3" (2) */

      BFT_MALLOC(family_list, m->n_families, cs_int_t);
      BFT_MALLOC(family_mask, m->n_families, int);

      for (i = 0; i < m->n_families; i++)
        family_mask[i] = 0;

      cs_selector_get_family_list("2",  &n_families, family_list);

      for (i = 0; i < n_families; i++)
        family_mask[family_list[i] - 1] += 1;

      cs_selector_get_family_list("3",  &n_families, family_list);

      for (i = 0; i < n_families; i++)
        family_mask[family_list[i] - 1] += 2;

      BFT_FREE(family_list);

      /* Now that mask is built, test for adjacency */

      for (face_id = 0; face_id < m->n_i_faces; face_id++) {

        /* Adjacent cells  and flags */

        fvm_lnum_t c1 = m->i_face_cells[face_id*2] - 1;
        fvm_lnum_t c2 = m->i_face_cells[face_id*2 + 1] - 1;

        int iflag1 = family_mask[m->cell_family[c1]];
        int iflag2 = family_mask[m->cell_family[c2]];

        /* Should the face belong to the extracted mesh ? */

        if ((iflag1 == 1 && iflag2 == 2) || (iflag1 == 2 && iflag2 == 1)) {
          i_face_list[n_i_faces] = face_id + 1;
          n_i_faces += 1;
        }

      }

      BFT_FREE(family_mask);
    }

    /* Select boundary faces */

    cs_selector_get_b_face_list("4", &n_b_faces, b_face_list);

    /* Define postprocessing mesh */

    cs_post_define_surface_mesh_by_list(5,               /* mesh id */
                                        "Mixed surface",
                                        n_i_faces,
                                        n_b_faces,
                                        i_face_list,
                                        b_face_list,
                                        false,           /* add_groups */
                                        false,           /* auto_variables */
                                        n_writers,
                                        writer_ids);

    /* Free selection lists */

    BFT_FREE(i_face_list);
    BFT_FREE(b_face_list);

  }

  /*--------------------------------------------------------------------------*/
 
  /* In cases where a mesh containing polygonal elements is output through
   * a writer configured to divide polygons into triangles (for example when
   * visualization tools do not support polygons, or when higly non convex
   * faces lead to visualization artifacts), it may be useful to extract a
   * mesh containing the edges of the original mesh so as to view the polygon
   * boundaries as an overlay.
   *
   * This is possible using the cs_post_define_edges_mesh() function,
   * whose parameters are:
   *   mesh_id       <-- id of edges mesh to create (< 0 reserved, > 0 for user)
   *   base_mesh_id  <-- id of existing mesh (< 0 reserved, > 0 for user)
   *   n_writers     <-- number of associated writers
   *   writer_ids    <-- ids of associated writers */

  if (false) {

    /* Example of edges mesh extraction */

    const int n_writers = 1;
    const int writer_ids[] = {4};  /* Associate to writer 4 */

    cs_post_define_edges_mesh(6, /* mesh_id */
                              1, /* base_mesh_id */
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
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

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
