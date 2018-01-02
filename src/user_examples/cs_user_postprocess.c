/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_geom.h"
#include "cs_interpolate.h"
#include "cs_mesh.h"
#include "cs_selector.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_probe.h"
#include "cs_time_plot.h"

#include "cs_field_pointer.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

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

/*! [post_select_func_1] */
static void
_i_faces_select_example(void         *input,
                        cs_lnum_t    *n_faces,
                        cs_lnum_t   **face_ids)
{
  CS_UNUSED(input);

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
/*! [post_select_func_1] */

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

/*! [post_select_func_2] */
static void
_b_faces_select_example(void         *input,
                        cs_lnum_t    *n_faces,
                        cs_lnum_t   **face_ids)
{
  CS_UNUSED(input);

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
/*! [post_select_func_2] */

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

/*! [post_select_func_3] */
static void
_he_fraction_05_select(void        *input,
                       cs_lnum_t   *n_cells,
                       cs_lnum_t  **cell_ids)
{
  CS_UNUSED(input);

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
/*! [post_select_func_3] */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void)
{
  /* Set time plot file writer flush behavior defaults. */

  /*! [post_set_tp_flush] */
  cs_time_plot_set_flush_default(1800, /* flush_wtime */
                                 -1);  /* n_buffer_steps */
  /*! [post_set_tp_flush] */

  /* Default output format and options */

  /* Redefine default writer */
  /* ----------------------- */

  /*! [post_define_writer_m1] */
  cs_post_define_writer(CS_POST_WRITER_DEFAULT,       /* writer_id */
                        "results",                    /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format_name */
                        "",                           /* format_options */
                        FVM_WRITER_FIXED_MESH,
                        false,                        /* output_at_start */
                        true,                         /* output_at_end */
                        -1,                           /* frequency_n */
                        -1.0);                        /* frequency_t */
  /*! [post_define_writer_m1] */

  /* Define additional writers */
  /* ------------------------- */

  /* Common parameters for all writers */

  /*! [post_define_writer_freq] */
  double frequency_n = -1.0;
  double frequency_t = -1.0;
  /*! [post_define_writer_freq] */

  /*! [post_define_writer_1] */
  cs_post_define_writer(1,                            /* writer_id */
                        "user_txt",                   /* writer name */
                        "postprocessing",             /* directory name */
                        "MED",                        /* format name */
                        "divide_polyhedra",
                        FVM_WRITER_FIXED_MESH,
                        false,                        /* output_at_start */
                        true,                         /* output_at_end */
                        -1,                           /* frequency_n */
                        -1.0);                        /* frequency_t */
  /*! [post_define_writer_1] */

  /*! [post_define_writer_2] */
  cs_post_define_writer(2,                            /* writer_id */
                        "modif",                      /* writer name */
                        "postprocessing",             /* directory name */
                        "ensight",                    /* format name */
                        "text",
                        FVM_WRITER_TRANSIENT_CONNECT,
                        false,                        /* output_at_start */
                        false,                        /* output_at_end */
                        3,
                        frequency_t);
  /*! [post_define_writer_2] */

  /*! [post_define_writer_3] */
  cs_post_define_writer(3,                /* writer_id */
                        "profile",        /* writer name */
                        "postprocessing", /* directory name */
                        "plot",           /* format name */
                        "",               /* format options */
                        FVM_WRITER_FIXED_MESH,
                        false,            /* output_at_start */
                        false,            /* output_at_end */
                        100,              /* nt_freq */
                        -1.0);            /* dt_freq */
  /*! [post_define_writer_3] */

  /*! [post_define_writer_4] */
  cs_post_define_writer(6,                        /* writer_id */
                        "Histogram",              /* writer name */
                        "histograms",             /* directory name */
                        "histogram",              /* format name */
                        "10 tex",                 /* format options */
                        FVM_WRITER_FIXED_MESH,
                        false,                    /* output_at_start */
                        true,                     /* output at end */
                        -1,                       /* time step frequency */
                        -1.0);                    /* time value frequency */
  /*! [post_define_writer_4] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void)
{
  /* Reconfigure predefined meshes (mesh_id -1 for volume, -2 for boundary */

  /* De-activate boundary mesh output by redefining it with no writer
     association (default is:
     int n_writers = 1;
     const int writer_ids[] = {CS_POST_WRITER_DEFAULT});
  */

  /*! [post_define_mesh_m2] */
  {
    int n_writers = 0;
    const int *writer_ids = NULL;

    cs_post_define_surface_mesh(CS_POST_MESH_BOUNDARY,  /* mesh_id */
                                "Boundary",  /* mesh name */
                                NULL,        /* interior face selection criteria */
                                "all[]",     /* boundary face selection criteria */
                                true,        /* add_groups */
                                true,        /* automatic variables output */
                                n_writers,
                                writer_ids);
  }
  /*! [post_define_mesh_m2] */

  /*--------------------------------------------------------------------------*/

  /* Example: select interior faces with y = 0.5 */

  /*! [post_define_mesh_1] */
  {
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
  /*! [post_define_mesh_1] */

  /*--------------------------------------------------------------------------*/

  /* Advanced example:
     Build a surface mesh containing interior faces separating cells of group "2"
     from those of group "3", (assuming no cell has both colors), as well as
     boundary faces of group "4". */

  /*! [post_define_mesh_3] */
  {
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
  /*! [post_define_mesh_3] */

  /* Advanced example:
     Build a (time varying) volume mesh containing cells
     with values of field named "He_fraction" > 0.05 */

  /*! [post_define_mesh_4] */
  {
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
  /*! [post_define_mesh_4] */

  /*--------------------------------------------------------------------------*/

  /* Example: extract face edges of another mesh */

  /*! [post_define_mesh_5] */
  {
    const int n_writers = 1;
    const int writer_ids[] = {4};  /* Associate to writer 4 */

    cs_post_define_edges_mesh(5, /* mesh_id */
                              1, /* base_mesh_id */
                              n_writers,
                              writer_ids);
  }
  /*! [post_define_mesh_5] */

  /*--------------------------------------------------------------------------*/

  /* Example: attach default txt histogram writer on boundary mesh */

  /*! [post_attach_mesh_1] */
  {
    cs_post_mesh_attach_writer(CS_POST_MESH_BOUNDARY, CS_POST_WRITER_HISTOGRAMS);
  }
  /*! [post_attach_mesh_1] */

  /*--------------------------------------------------------------------------*/

  /* Example: attach user tex histogram writer of id 6 on volume mesh */

  /*! [post_attach_mesh_2] */
  {
    cs_post_mesh_attach_writer(CS_POST_MESH_VOLUME, 6);
  }
  /*! [post_attach_mesh_2] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define monitoring probes and profiles.
 *
 * Profiles are defined as sets of probes.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_probes(void)
{
  /* Define monitoring probes */

  /* A writer (id = CS_POST_WRITER_PROBES) using the format "time_plot" is
     associated by default to a set of monitoring probes.
     This is not the case for a profile. */

  /*! [post_define_probes_1] */
  {
    cs_probe_set_t  *pset = cs_probe_set_create("Monitoring");

    cs_probe_set_add_probe(pset, 0.25, 0.025, 0.025, "M1");
    cs_probe_set_add_probe(pset, 0.50, 0.025, 0.025, "M2");
    cs_probe_set_add_probe(pset, 0.75, 0.025, 0.025, "M3");

  }
  /*! [post_define_probes_1] */

  /* Add a first profile */
  {
    cs_coord_3_t  start = {0., 0.025, 0.025};
    cs_coord_3_t  end = {1., 0.025, 0.025};
    int  writer_ids[] = {2};

    cs_probe_set_t  *pset =
      cs_probe_set_create_from_segment("Prof1", // name
                                       11,      // n_probes
                                       start,   // start coordinates
                                       end);    // end coordinates

    cs_probe_set_associate_writers(pset, 1, writer_ids);
  }

  /* Add a second profile attached to boundary vertices */
  {
    cs_coord_3_t  start = {0., 0., 0.};
    cs_coord_3_t  end = {1., 0., 0.};

    cs_probe_set_create_from_segment("P1",    // name
                                     11,      // n_probes
                                     start,   // start coordinates
                                     end);    // end coordinates
  }

  /* Add a second profile attached to boundary vertices */
  {
    cs_coord_3_t  start = {0., 0., 0.};
    cs_coord_3_t  end = {1., 0., 0.};
    int  writer_ids[] = {5};

    cs_probe_set_t *pset =
      cs_probe_set_create_from_segment("P2",     // name
                                       21,       // n_probes
                                       start,    // start coordinate
                                       end);     // end coordinate

    cs_probe_set_associate_writers(pset, 1, writer_ids);

    cs_probe_set_option(pset, "boundary", "true");
    cs_probe_set_snap_mode(pset, CS_PROBE_SNAP_VERTEX);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, NULL otherwise
 * \param[in]       n_cells      local number of cells of post_mesh
 * \param[in]       n_i_faces    local number of interior faces of post_mesh
 * \param[in]       n_b_faces    local number of boundary faces of post_mesh
 * \param[in]       n_vertices   local number of vertices faces of post_mesh
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing
 *                               mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertices,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts)
{
  /* Output of k=1/2(R11+R22+R33) for the Rij-epsilon model
     ------------------------------------------------------ */

  /*< [postprocess_values_ex_1] */
  if (cat_id == CS_POST_MESH_VOLUME && cs_glob_turb_model->itytur == 3) {

    cs_real_t *s_cell;
    BFT_MALLOC(s_cell, n_cells, cs_real_t);

    if (cs_glob_turb_rans_model->irijco) {
      const cs_real_6_t *cvar_r = (const cs_real_6_t *)(CS_F_(rij)->val);
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        cs_lnum_t cell_id = cell_list[i];
        s_cell[i] = 0.5* (  cvar_r[cell_id][0]
                          + cvar_r[cell_id][1]
                          + cvar_r[cell_id][2]);
      }
    }

    else {
      const cs_real_t *cvar_r11 = CS_F_(r11)->val;
      const cs_real_t *cvar_r22 = CS_F_(r22)->val;
      const cs_real_t *cvar_r33 = CS_F_(r33)->val;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        cs_lnum_t cell_id = cell_list[i];
        s_cell[i] = 0.5* (  cvar_r11[cell_id]
                          + cvar_r22[cell_id]
                          + cvar_r33[cell_id]);
      }
    }

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "Turb energy",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    BFT_FREE(s_cell);
  }
  /*< [postprocess_values_ex_1] */

}

/*----------------------------------------------------------------------------*/
/*!
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * \param  nt_max_abs  maximum time step number
 * \param  nt_cur_abs  current time step number
 * \param  t_cur_abs   absolute time at the current time step
 */
/*----------------------------------------------------------------------------*/

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

  /* Example: deactivate all output before time step 1000 */

  /*! [post_activate] */
  if (nt_max_abs < 1000) {
    int writer_id = 0; /* 0: all writers */
    cs_post_activate_writer(writer_id, false);
  }
  /*! [post_activate] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
