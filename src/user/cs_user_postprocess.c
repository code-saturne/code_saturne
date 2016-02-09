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
#include "cs_probe.h"
#include "cs_post.h"
#include "cs_time_plot.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Example function for advanced selection of interior faces.
 *
 * Selects interior faces separating cells of group "2" from those
 * of group "3", (assuming no cell has both colors).
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_i_faces_select_example(void         *input,
                        cs_lnum_t    *n_faces,
                        cs_lnum_t   **face_ids)
{
  cs_lnum_t i, face_id;
  cs_lnum_t n_families = 0;
  cs_int_t *family_list = NULL;
  int *family_mask = NULL;

  cs_lnum_t n_i_faces = 0;
  cs_lnum_t *i_face_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  /* Allocate selection list */

  BFT_MALLOC(i_face_ids, m->n_i_faces, cs_lnum_t);

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

    cs_lnum_t c1 = m->i_face_cells[face_id][0];
    cs_lnum_t c2 = m->i_face_cells[face_id][1];

    int iflag1 = family_mask[m->cell_family[c1]];
    int iflag2 = family_mask[m->cell_family[c2]];

    /* Should the face belong to the extracted mesh ? */

    if ((iflag1 == 1 && iflag2 == 2) || (iflag1 == 2 && iflag2 == 1)) {
      i_face_ids[n_i_faces] = face_id;
      n_i_faces += 1;
    }

  }

  /* Free memory */

  BFT_FREE(family_mask);
  BFT_REALLOC(i_face_ids, n_i_faces, cs_lnum_t);

  /* Set return values */

  *n_faces = n_i_faces;
  *face_ids = i_face_ids;
}

/*----------------------------------------------------------------------------
 * Example function for selection of boundary faces.
 *
 * selects boundary faces of group "4".
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_b_faces_select_example(void         *input,
                        cs_lnum_t    *n_faces,
                        cs_lnum_t   **face_ids)
{
  cs_lnum_t n_b_faces = 0;
  cs_lnum_t *b_face_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  /* Allocate selection list */

  BFT_MALLOC(b_face_ids, m->n_b_faces, cs_lnum_t);

  /* Use simple selection function */

  cs_selector_get_b_face_list("4", &n_b_faces, b_face_ids);

  /* Adjust array to final size (cleaner, but not required) */

  BFT_REALLOC(b_face_ids, n_b_faces, cs_lnum_t);

  /* Set return values */

  *n_faces = n_b_faces;
  *face_ids = b_face_ids;
}

/*----------------------------------------------------------------------------
 * Example function for selection of cells with scalar field values above
 * a certain threshold.
 *
 * In this example, the selection is base on the value of a scalar field
 * named "he_fraction" being above above 0.05.
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_cells  --> number of selected cells
 *   cell_ids --> array of selected cell ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_he_fraction_05_select(void        *input,
                       cs_lnum_t   *n_cells,
                       cs_lnum_t  **cell_ids)
{
  cs_lnum_t _n_cells = 0;
  cs_lnum_t *_cell_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  cs_field_t *f = cs_field_by_name("He_fraction"); /* Get access to field */

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "No field with name \"He_fraction\" defined");

  /* Before time loop, field is defined, but has no values yet,
     so ignore that case (postprocessing mesh will be initially empty) */

  if (f->val != NULL) {

    BFT_MALLOC(_cell_ids, m->n_cells, cs_lnum_t); /* Allocate selection list */

    for (cs_lnum_t i = 0; i < m->n_cells; i++) {
      if (f->val[i] > 5.e-2) {
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
   *   'hdf5'              (force HDF5 file type, usually the default for CGNS)
   *   'discard_polygons'  (ignore polygon-type faces)
   *   'discard_polyhedra' (ignore polyhedron-type cells)
   *   'divide_polygons'   (subdivides polygon-type faces)
   *   'divide_polyhedra'  (subdivides polyhedron-type cells)
   *   'split_tensors'     (writes tensors as separate scalars)
   *   'separate_meshes' */

  /* Set time plot file writer flush behavior defaults.
   *   flush_wtime     <-- elapsed time interval between file flushes;
   *                       if < 0, no forced flush
   *   n_buffer_steps  <-- number of time steps in output buffer if
   *                       file is not to be kept open */

  if (false)
    cs_time_plot_set_flush_default(1800, /* flush_wtime */
                                   -1);  /* n_buffer_steps */

  /* Default writer time dependency */

  fvm_writer_time_dep_t   time_dep = FVM_WRITER_FIXED_MESH;

  /* Default time step or physical time based output frequencies */

  bool       output_at_end = true;
  int        frequency_n = -1;
  double     frequency_t = -1.0;

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
                          frequency_n,
                          frequency_t);

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
                          frequency_n,
                          frequency_t);

  if (false)
    cs_post_define_writer(2,                            /* writer_id */
                          "user_txt",                   /* writer name */
                          "postprocessing",             /* directory name */
                          "ensight",                    /* format name */
                          "text, divide_polyhedra",
                          time_dep,                     /* modification flag */
                          output_at_end,
                          frequency_n,
                          frequency_t);

  if (false)
    cs_post_define_writer(3,                            /* writer_id */
                          "modif",                      /* writer name */
                          "postprocessing",             /* directory name */
                          "ensight",                    /* format name */
                          "discard_polyhedra",
                          FVM_WRITER_TRANSIENT_CONNECT,
                          false,
                          3,
                          frequency_t);

  if (false)
    cs_post_define_writer(4,                            /* writer_id */
                          "exchange",                   /* writer name */
                          "postprocessing",             /* directory name */
                          "MED",                        /* format name  */
                          "",                           /* format options */
                          FVM_WRITER_TRANSIENT_COORDS,  /* modification flag */
                          false,
                          frequency_n,
                          frequency_t);

  if (false)
          cs_post_define_writer(5,               /* writer_id */
                             "profile",        /* writer name */
                         "postprocessing", /* directory name */
                         "plot",           /* format name */
                         "",               /* format options */
                         FVM_WRITER_FIXED_MESH,
                         false,            // output_at_end
                         100,              // nt_freq
                         -1.0);            // dt_freq

}

/*-----------------------------------------------------------------------------
 * Define monitoring probes and profiles. A profile is seen as a set of probes.
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_probes(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Post-processing for probes and profiles relies on the following functions
   * whose prototypes are defined in cs_probe.h
   *
   * cs_probe_set_t structure handles the management of probes for monitoring
   * or profile.
   *
   * Call cs_probe_set_create_from_array() for defining a list of probes whose
   * parameters are:
   * name      <--  name of the set of probes
   * n_probes  <--  number of probes in coords and labels
   * coords    <--  list of coordinates related to each probe
   * labels    <--  list of label related to each probe (optional)
   *
   * Call cs_probe_set_create_from_segment() for defining a profile
   * name         <--  name of the set of probes
   * n_probes     <--  number of probes defining the profile
   * start_coord  <--  coordinate of the starting point
   * start_label  <--  NULL or label of the starting point (optional)
   * end_coord    <--  coordinate of the ending point
   * end_label    <--  NULL or label of the ending point (optional)
   *
   * For a more advanced usage, the functions cs_probe_set_create() and
   * cs_probe_set_add_probe() can be called

   * parameters for cs_probe_set_create():
   * name    <--  name of the set of probes
   *
   * parameters for cs_probe_set_add_probe():
   * pset    <->  set of probes (pointer to a cs_probe_set_t structure)
   * xyz     <--  coordinates of the point to add
   * label   <--  NULL or the name of the point (optional)
   *
   * All cs_probe_set_create_* functions return a pointer to a cs_probe_set_t
   * structure which can be used to specify additional options thanks to
   * cs_probe_set_option() with parameters:
   * pset     <->  set of probes (pointer to a cs_probe_set_t structure)
   * keyname  <--  name of the keyword related to the parameter to set
   * keyval   <--  value of the keyword to set
   *
   * Available keynames are the following:
   *
   * "activated" where keyval is either "true" or "false" (default)
   * "boundary"           where keyval is either "true" or "false" (default)
   * "mode"               where keyval is "exact", "nearest_vertex" or
   *                      "nearest_cell_center (default)
   * "profile"            where keyval is either "true" or "false"
   * "selection_criteria" where keyval is a string like "x > 0.5"
   * "tolerance"          where keyval is for instance "0.05" (default "0.10")
   * "moving_probes"      where keyval is either "true" or "false"
   *
   * More advanced functions are described along with examples below. */

  if (false) { /* Define monitoring probes */

    cs_real_3_t  m1 = {0.25, 0.025, 0.025};
    cs_real_3_t  m2 = {0.50, 0.025, 0.025};
    cs_real_3_t  m3 = {0.75, 0.025, 0.025};
    cs_probe_set_t  *pset = cs_probe_set_create("Monitoring");

    cs_probe_set_add_probe(pset, m1, "M1");
    cs_probe_set_add_probe(pset, m2, "M2");
    cs_probe_set_add_probe(pset, m3, "M3");

    /* A writer (id = -3) using the format "time_plot" is associated by default
       to a set of monitoring probes. This is not the case for a profile. */

  }

  if (true) { /* Add a first profile */

    cs_coord_3_t  start = {0., 0.025, 0.025};
    cs_coord_3_t  end = {1., 0.025, 0.025};
    int  writer_ids[] = {2};

    cs_probe_set_t  *pset =
      cs_probe_set_create_from_segment("Prof1", // name
                                       11,      // n_probes
                                       start,   // start coordinate
                                       NULL,    // start label
                                       end,     // end coordinate
                                       NULL);   // end label

    //    cs_probe_set_associate_writers(pset, 1, writer_ids);
  }

  if (true)  {  /* Add a second profile attached to border vertices */

    cs_coord_3_t  start = {0., 0., 0.};
    cs_coord_3_t  end = {1., 0., 0.};

    cs_probe_set_create_from_segment("P1",    // name
                                     11,      // n_probes
                                     start,   // start coordinate
                                     NULL,    // start label
                                     end,     // end coordinate
                                     NULL);   // end label

  }

  if (false)  {  /* Add a second profile attached to border vertices */

    cs_coord_3_t  start = {0., 0., 0.};
    cs_coord_3_t  end = {1., 0., 0.};
    int  writer_ids[] = {5};

    cs_probe_set_t *pset =
      cs_probe_set_create_from_segment("P2",     // name
                                       21,       // n_probes
                                       start,    // start coordinate
                                       NULL,     // start label
                                       end,      // end coordinate
                                       NULL);    // end label

    cs_probe_set_associate_writers(pset, 1, writer_ids);

    cs_probe_set_option(pset, "boundary", "true");
    cs_probe_set_option(pset, "mode", "nearest_vertex");

  }

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

    cs_post_define_alias_mesh(2,                 /* mesh id */
                              -2,                /* aliased mesh id */
                              false,             /* auto_variables */
                              n_writers,
                              writer_ids);

  }

  /*--------------------------------------------------------------------------*/

  /* More advanced mesh element selection is possible using
   * cs_post_define_volume_mesh_by_func() or
   * cs_post_define_surface_mesh_by_func(), which allow defining
   * volume or surface meshes using user-defined element lists.
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
   *   writer_ids        <-- ids of associated writers
   *
   * parameters for cs_post_define_surface_mesh_by_func():
   *   mesh_id             <-- id of mesh to define (< 0 reserved, > 0 for user)
   *   mesh_name           <-- associated mesh name
   *   i_face_select_func  <-- pointer to interior faces selection function
   *   b_face_select_func  <-- pointer to boundary faces selection function
   *   i_face_select_input <-> pointer to optional input data for the interior
   *                           faces selection function, or NULL
   *   b_face_select_input <-> pointer to optional input data for the boundary
   *                           faces selection function, or NULL
   *   time_varying        <-- if true, try to redefine mesh at each output time
   *   add_groups          <-- if true, add group information if present
   *   auto_variables      <-- if true, automatic output of main variables
   *   n_writers           <-- number of associated writers
   *   writer_ids          <-- ids of associated writers */

  if (false) {

    /* Advanced example:
       Build a surface mesh containing interior faces separating cells of group "2"
       from those of group "3", (assuming no cell has both colors), as well as
       boundary faces of group "4". */

    const int n_writers = 1;
    const int writer_ids[] = {1};  /* Associate to writer 1 */

    /* Define postprocessing mesh */

    cs_post_define_surface_mesh_by_func(3,               /* mesh id */
                                        "Mixed surface",
                                        _i_faces_select_example,
                                        _b_faces_select_example,
                                        NULL,            /* i_faces_sel_input */
                                        NULL,            /* b_faces_sel_input */
                                        false,           /* time varying */
                                        false,           /* add_groups */
                                        false,           /* auto_variables */
                                        n_writers,
                                        writer_ids);
  }

  /* Advanced example:
     Build a (time varying) volume mesh containing cells
     with values of field named "He_fraction" > 0.05 */

  if (false) {

    const int n_writers = 1;
    const int writer_ids[] = {2};  /* Associate to writer 2 */

    /* Define postprocessing mesh */

    cs_post_define_volume_mesh_by_func(4,               /* mesh id */
                                       "He_fraction_05",
                                       _he_fraction_05_select,
                                       NULL,            /* _c_05_select_input */
                                       true,            /* time varying */
                                       false,           /* add_groups */
                                       false,           /* auto_variables */
                                       n_writers,
                                       writer_ids);
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

    cs_post_define_edges_mesh(5, /* mesh_id */
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
