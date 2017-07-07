/*============================================================================
 * Management of the post-processing
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "fvm_nodal.h"
#include "fvm_nodal_append.h"
#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_file.h"
#include "cs_lagr_extract.h"
#include "cs_log.h"
#include "cs_lagr_query.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_post.c

  \brief Post-processing management.

  \var  CS_POST_ON_LOCATION
        postprocess variables on their base location (volume for variables)
  \var  CS_POST_BOUNDARY_NR
        postprocess boundary without reconstruction

  \enum cs_post_type_t

  \brief Postprocessing input variable type

  \var CS_POST_TYPE_cs_int_t
       Fortran integer
  \var CS_POST_TYPE_cs_real_t
       Fortran double precision
  \var CS_POST_TYPE_int
       integer
  \var CS_POST_TYPE_float
       single precision floating-point value
  \var CS_POST_TYPE_double
       double precision floating-point value

  \typedef  cs_post_elt_select_t

  \brief  Function pointer to elements selection definition

  Each function of this sort may be used to select a given type of element,
  usually cells, interior faces, or boundary faces.

  If non-empty and not containing all elements, a list of elements of the
  main mesh should be allocated (using BFT_MALLOC) and defined by this
  function when called. This list's lifecycle is then managed by the
  postprocessing subsystem.

   Note: if the input pointer is non-NULL, it must point to valid data
   when the selection function is called, so either:
   - that value or structure should not be temporary (i.e. local);
   - post-processing output must be ensured using cs_post_write_meshes()
   with a fixed-mesh writer before the data pointed to goes out of scope;

  \param[in, out]  input     pointer to optional (untyped) value or structure
  \param[out]      n_elts    number of selected elements
  \param[out]      elt_list  list of selected elements (0 to n-1 numbering)

  \typedef  cs_post_time_dep_output_t

  Function pointer associated with a specific post-processing output.

  Such functions are registered using the \ref cs_post_add_time_dep_output,
  and all registered functions are automatically called by
  \ref cs_post_write_vars.

  Note: if the input pointer is non-NULL, it must point to valid data
  when the output function is called, so either:
  - that value or structure should not be temporary (i.e. local);
  - post-processing output must be ensured using cs_post_write_var()
  or similar before the data pointed to goes out of scope.

  \param[in, out]  input       pointer to optional (untyped) value or structure
  \param[in]       nt_cur_abs  current time step number
  \param[in]       t_cur_abs   absolute time at the current time step

  \typedef cs_post_time_mesh_dep_output_t

  Function pointer associated with a specific post-processing output
  on multiple meshes.

  Such functions are registered using the cs_post_add_time_mesh_dep_vars(),
  and all registered functions are automatically called by
  cs_post_write_vars().

  Note: if the input pointer is non-NULL, it must point to valid data
  when the output function is called, so either:
  - that value or structure should not be temporary (i.e. local);
  - post-processing output must be ensured using cs_post_write_var()
  or similar before the data pointed to goes out of scope.

  \param[in, out]  input        pointer to optional (untyped) value or structure
  \param[in]       mesh_id      id of the output mesh for the current call
  \param[in]       cat_id       category id of the output mesh for the
                                current call
  \param[in]       ent_flag     indicate global presence of cells
                                (ent_flag[0]), interior faces (ent_flag[1]),
                                boundary faces (ent_flag[2]), particles
                                (ent_flag[3]) or probes (ent_flag[4])
  \param[in]       n_cells      local number of cells of post_mesh
  \param[in]       n_i_faces    local number of interior faces of post_mesh
  \param[in]       n_b_faces    local number of boundary faces of post_mesh
  \param[in]       cell_list    list of cells (1 to n) of post-processing mesh
  \param[in]       i_face_list  list of interior faces (1 to n) of
                                post-processing mesh
  \param[in]       b_face_list  list of boundary faces (1 to n) of
                                post-processing mesh
  \param[in]       nt_cur_abs   current time step number
  \param[in]       t_cur_abs    current physical time
  \param[in]       nt_cur_abs   current time step number
  \param[in]       t_cur_abs    absolute time at the current time step
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define _MIN_RESERVED_MESH_ID    CS_POST_MESH_PROBES
#define _MIN_RESERVED_WRITER_ID  CS_POST_WRITER_PROFILES

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Specific (forced) writer output times */
/*---------------------------------------*/

typedef struct {

  int      n_t_steps_max ;   /* Max. number of forced time steps */
  int      n_t_vals_max;     /* Max. number of forced time values */

  int      n_t_steps;        /* Number of forced time steps */
  int      n_t_vals;         /* Number of forced time values */

  int     *t_steps;          /* Forced output time steps (unordered) */
  double  *t_vals;           /* Forced output time values (unordered) */

} cs_post_writer_times_t;

/* Writer structure definition parameters */
/*----------------------------------------*/

typedef struct {

  fvm_writer_time_dep_t   time_dep;     /* Time dependency */
  int                     fmt_id;       /* format id */
  char                   *case_name;    /* Case (writer) name */
  char                   *dir_name;     /* Associated directory name */
  char                   *fmt_opts;     /* Format options */

} cs_post_writer_def_t;

/* Mesh location type */
/*--------------------*/

typedef enum {

  CS_POST_LOCATION_CELL,         /* Values located at cells */
  CS_POST_LOCATION_I_FACE,       /* Values located at interior faces */
  CS_POST_LOCATION_B_FACE,       /* Values located at boundary faces */
  CS_POST_LOCATION_VERTEX,       /* Values located at vertices */
  CS_POST_LOCATION_PARTICLE,     /* Values located at particles */

} cs_post_location_t;

/* Writer structure */
/*------------------*/

/* This object is based on a choice of a case, directory, and format,
   as well as a flag for associated mesh's time dependency, and the default
   output frequency for associated variables. */

typedef struct {

  int            id;            /* Identifier (< 0 for "reservable" writer,
                                 * > 0 for user writer */
  int            output_start;  /* Output at start of calculation if nonzero */
  int            output_end;    /* Output at end of calculation if nonzero */
  int            frequency_n;   /* Default output frequency in time-steps */
  double         frequency_t;   /* Default output frequency in seconds */

  int            active;        /* -1 if blocked at this stage,
                                   0 if no output at current time step,
                                   1 in case of output */
  int            n_last;        /* Time step number for the last
                                   activation (-1 before first output) */
  double         t_last;        /* Time value number for the last
                                   activation (0.0 before first output) */

  cs_post_writer_times_t  *ot;  /* Specific output times */
  cs_post_writer_def_t    *wd;  /* Associated writer definition */

  fvm_writer_t  *writer;        /* Associated FVM writer */

} cs_post_writer_t;

/* Post-processing mesh structure */
/*--------------------------------*/

/* This object manages the link between an exportable mesh and
   associated writers. */

typedef struct {

  int                     id;            /* Identifier (< 0 for "reservable"
                                            mesh, > 0 for user mesh */

  char                   *name;          /* Mesh name */
  char                   *criteria[5];   /* Base selection criteria for
                                            cells, interior faces,
                                            boundary faces, and particles
                                            respectively */
  cs_post_elt_select_t   *sel_func[5];   /* Advanced selection functions for
                                            cells, interior faces,
                                            boundary faces, particles and
                                            probes respectively */
  void                   *sel_input[5];  /* Advanced selection input for
                                            matching selection functions */
  int                     ent_flag[5];   /* Presence of cells (ent_flag[0],
                                            interior faces (ent_flag[1]),
                                            boundary faces (ent_flag[2]),
                                            or particles (ent_flag[3] = 1
                                            for particles, 2 for trajectories),
                                            probes (ent_flag[4] = 1 for
                                            monitoring probes, 2 for profile)
                                            on one processor at least */

  int                     cat_id;        /* Optional category id as regards
                                            variables output (CS_POST_MESH_...,
                                            0 by default) */

  int                     edges_ref;     /* Base mesh for edges mesh */
  int                     locate_ref;    /* Base mesh for location mesh */

  bool                    add_groups;    /* Add group information if present */
  bool                    post_domain;   /* Output domain number in parallel
                                            if true */
  bool                    time_varying;  /* Time varying if associated writers
                                            allow it */

  int                     n_writers;     /* Number of associated writers */
  int                    *writer_id;     /* Array of associated writer ids */
  int                     nt_last;       /* Time step number for the last
                                            output (-2 before first output,
                                            -1 for time-indepedent output) */

  cs_lnum_t               n_i_faces;     /* N. associated interior faces */
  cs_lnum_t               n_b_faces;     /* N. associated boundary faces */

  double                  density;       /* Particles density in case
                                            of particle mesh */

  const fvm_nodal_t      *exp_mesh;      /* Associated exportable mesh */
  fvm_nodal_t            *_exp_mesh;     /* Associated exportable mesh,
                                            if owner */

  fvm_writer_time_dep_t   mod_flag_min;  /* Minimum mesh time dependency */
  fvm_writer_time_dep_t   mod_flag_max;  /* Maximum mesh time dependency */

} cs_post_mesh_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default output format and options */

static int _cs_post_default_format_id = 0;
static char *_cs_post_default_format_options = NULL;

/* Backup of initial vertex coordinates */

static bool        _cs_post_deformable = false;
static cs_real_t  *_cs_post_ini_vtx_coo = NULL;

/* Minimum global mesh time dependency */

fvm_writer_time_dep_t  _cs_post_mod_flag_min = FVM_WRITER_FIXED_MESH;

/* Flag for stable numbering of particles */

static bool        _number_particles_by_coord = false;

/* Array of exportable meshes associated with post-processing;
   free ids start under the last CS_POST_MESH_* definition,
   currently at -5) */

static int              _cs_post_min_mesh_id = _MIN_RESERVED_MESH_ID;
static int              _cs_post_n_meshes = 0;
static int              _cs_post_n_meshes_max = 0;
static cs_post_mesh_t  *_cs_post_meshes = NULL;

/* Array of writers for post-processing; */
/* writers CS_POST_WRITER_... are reserved */

static int                _cs_post_min_writer_id = _MIN_RESERVED_WRITER_ID;
static int                _cs_post_n_writers = 0;
static int                _cs_post_n_writers_max = 0;
static cs_post_writer_t  *_cs_post_writers = NULL;

/* Array of registered variable output functions and instances */

static int                _cs_post_n_output_tp = 0;
static int                _cs_post_n_output_tp_max = 0;

static int                _cs_post_n_output_mtp = 0;
static int                _cs_post_n_output_mtp_max = 0;

static cs_post_time_dep_output_t  **_cs_post_f_output_tp = NULL;
static void                       **_cs_post_i_output_tp = NULL;

static cs_post_time_mesh_dep_output_t  **_cs_post_f_output_mtp = NULL;
static void                            **_cs_post_i_output_mtp = NULL;

/* Default directory name */

static const char  _cs_post_dirname[] = "postprocessing";

/* Timer statistics */

static int  _post_out_stat_id = -1;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_post_activate_by_time_step(void);

void
cs_f_post_write_var(int               mesh_id,
                    const char       *var_name,
                    int               var_dim,
                    bool              interlace,
                    bool              use_parent,
                    int               nt_cur_abs,
                    double            t_cur_abs,
                    const cs_real_t  *cel_vals,
                    const cs_real_t  *i_face_vals,
                    const cs_real_t  *b_face_vals);

/*============================================================================
 * Prototypes for user functions called only by functions from this module.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * User function for output of values on a post-processing mesh.
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
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertice,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Clear temporary writer definition information.
 *
 * parameters:
 *   writer <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_destroy_writer_def(cs_post_writer_t  *writer)
{
  assert(writer != NULL);
  if (writer->wd != NULL) {
    cs_post_writer_def_t  *wd = writer->wd;
    BFT_FREE(wd->case_name);
    BFT_FREE(wd->dir_name);
    BFT_FREE(wd->fmt_opts);
    BFT_FREE(writer->wd);
  }
}

/*----------------------------------------------------------------------------
 * Print writer information to log file
 *----------------------------------------------------------------------------*/

static void
_writer_info(void)
{
  if (cs_glob_rank_id < 1) {

    int i;

    bft_printf(_("\n"
                 "Postprocessing output writers:\n"
                 "------------------------------\n\n"));

    for (i = 0; i < _cs_post_n_writers; i++) {

      int fmt_id = 0, n_fmt_str = 0;
      fvm_writer_time_dep_t   time_dep = FVM_WRITER_FIXED_MESH;
      const char  *fmt_name, *fmt_opts = NULL;
      const char  *case_name = NULL, *dir_name = NULL;
      const char empty[] = "";
      char frequency_s[80] = "";

      const cs_post_writer_t  *writer = _cs_post_writers + i;

      if (writer->wd != NULL) {
        const cs_post_writer_def_t *wd = writer->wd;
        fmt_id = wd->fmt_id;
        time_dep = wd->time_dep;
        fmt_opts = wd->fmt_opts;
        case_name = wd->case_name;
        dir_name = wd->dir_name;
      }
      else if (writer->writer != NULL) {
        const fvm_writer_t *w = writer->writer;
        fmt_id = fvm_writer_get_format_id(fvm_writer_get_format(w));
        time_dep = fvm_writer_get_time_dep(w);
        case_name = fvm_writer_get_name(w);
        fmt_opts = fvm_writer_get_options(w);
        dir_name = fvm_writer_get_path(w);
      }
      if (fmt_opts == NULL)
        fmt_opts = empty;

      n_fmt_str = fvm_writer_n_version_strings(fmt_id);
      if (n_fmt_str == 0)
        fmt_name = fvm_writer_format_name(fmt_id);
      else
        fmt_name = fvm_writer_version_string(fmt_id, 0, 0);

      if (writer->output_end != 0) {
        if (writer->frequency_t > 0)
          snprintf(frequency_s, 79,
                   _("every %12.5e s and at calculation end"),
                   writer->frequency_t);
        else if (writer->frequency_n >= 0)
          snprintf(frequency_s, 79,
                   _("every %d time steps and at calculation end"),
                   writer->frequency_n);
        else
          snprintf(frequency_s, 79, _("at calculation end"));
      }
      else {
        if (writer->frequency_t > 0)
          snprintf(frequency_s, 79, _("every %12.5e s"),
                   writer->frequency_t);
        else if (writer->frequency_n >= 0)
          snprintf(frequency_s, 79, _("every %d time steps"),
                   writer->frequency_n);
      }
      frequency_s[79] = '\0';

      bft_printf(_("  %2d: name: %s\n"
                   "      directory: %s\n"
                   "      format: %s\n"
                   "      options: %s\n"
                   "      time dependency: %s\n"
                   "      output: %s\n\n"),
                 writer->id, case_name, dir_name, fmt_name, fmt_opts,
                 _(fvm_writer_time_dep_name[time_dep]), frequency_s);
    }
  }
}

/*----------------------------------------------------------------------------
 * Initialize a writer; this creates the FVM writer structure, and
 * clears the temporary writer definition information.
 *
 * parameters:
 *   writer <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_init_writer(cs_post_writer_t  *writer)
{
  assert(writer != NULL);

  if (writer->writer == NULL) {

    cs_post_writer_def_t  *wd = writer->wd;

    /* Sanity checks */
    assert(writer->wd != NULL);
    if (wd->fmt_id >= fvm_writer_n_formats())
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid format name for writer (case: %s, dirname: %s)."),
                wd->case_name, wd->dir_name);

    writer->writer = fvm_writer_init(wd->case_name,
                                     wd->dir_name,
                                     fvm_writer_format_name(wd->fmt_id),
                                     wd->fmt_opts,
                                     wd->time_dep);
    _destroy_writer_def(writer);

  }

}

/*----------------------------------------------------------------------------
 * Free a writer's forced output time values.
 *
 * parameters:
 *   w <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_free_writer_times(cs_post_writer_t  *w)
{
  assert(w != NULL);

  if (w->ot == NULL) {
    BFT_FREE(w->ot->t_vals);
    BFT_FREE(w->ot->t_steps);
    BFT_FREE(w->ot);
  }
}

/*----------------------------------------------------------------------------
 * Create a specific writer output times structure.
 *
 * returns:
 *   structure for handling of specific output times
 *----------------------------------------------------------------------------*/

static cs_post_writer_times_t *
_writer_times_create(void)
{
  cs_post_writer_times_t  *ot;
  BFT_MALLOC(ot, 1, cs_post_writer_times_t);

  ot->n_t_steps_max = 0;
  ot->n_t_vals_max = 0;

  ot->n_t_steps = 0;
  ot->n_t_vals = 0;

  ot->t_steps = NULL;
  ot->t_vals = NULL;

  return ot;
}

/*----------------------------------------------------------------------------
 * Add an activation time step for a specific writer.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   nt        <-- time step value to add (or remove)
 *----------------------------------------------------------------------------*/

static void
_add_writer_ts(cs_post_writer_t  *w,
               int                nt)
{
  int prev_id;
  int nt_abs = CS_ABS(nt);

  if (w->ot == NULL)
    w->ot = _writer_times_create();

  /* Search for previous value */

  for (prev_id = 0; prev_id < w->ot->n_t_steps; prev_id++) {
    if (w->ot->t_steps[prev_id] == nt_abs)
      break;
  }

  /* If value already present */

  if (prev_id < w->ot->n_t_steps) {

    /* Remove previous value from unsorted list (swap with last, remove last) */

    if (nt < 0) {
      w->ot->t_steps[prev_id] = w->ot->t_steps[w->ot->n_t_steps - 1];
      w->ot->n_t_steps -= 1;
    }

  }

  /* If values not already present */

  else if (nt > -1) {

    if (w->ot->n_t_steps_max < w->ot->n_t_steps + 1) {
      if (w->ot->n_t_steps_max == 0)
        w->ot->n_t_steps_max = 1;
      else
        w->ot->n_t_steps_max *= 2;
      BFT_REALLOC(w->ot->t_steps, w->ot->n_t_steps_max, int);
    }

    w->ot->t_steps[w->ot->n_t_steps] = nt;
    w->ot->n_t_steps += 1;

  }
}

/*----------------------------------------------------------------------------
 * Add an activation time value for a specific writer.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   t         <-- time value to add (or remove)
 *----------------------------------------------------------------------------*/

static void
_add_writer_tv(cs_post_writer_t  *w,
               double             t)
{
  int prev_id;
  double t_abs = CS_ABS(t);

  if (w->ot == NULL)
    w->ot = _writer_times_create();

  /* Search for previous value */

  for (prev_id = 0; prev_id < w->ot->n_t_steps; prev_id++) {
    double td = w->ot->t_vals[prev_id] - t_abs;
    if (td > -1.e-35 && td < 1.e-35)
      break;
  }

  /* If value already present */

  if (prev_id < w->ot->n_t_vals) {

    /* Remove previous value from unsorted list (swap with last, remove last) */

    if (t < 0.) {
      w->ot->t_vals[prev_id] = w->ot->t_vals[w->ot->n_t_vals - 1];
      w->ot->n_t_vals -= 1;
    }

  }

  /* If values not already present */

  else if (t >= 0.) {

    if (w->ot->n_t_vals_max < w->ot->n_t_vals + 1) {
      if (w->ot->n_t_vals_max == 0)
        w->ot->n_t_vals_max = 1;
      else
        w->ot->n_t_vals_max *= 2;
      BFT_REALLOC(w->ot->t_vals, w->ot->n_t_vals_max, double);
    }

    w->ot->t_vals[w->ot->n_t_vals] = t;
    w->ot->n_t_vals += 1;

  }
}

/*----------------------------------------------------------------------------
 * Update "active" or "inactive" flag of a writer based on specified
 * output lists.
 *
 * parameters:
 *   w  <-> pointer to writer structure
 *   ts <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_activate_if_listed(cs_post_writer_t      *w,
                    const cs_time_step_t  *ts)
{
  int  i;
  bool force_status = false;
  int prev_status = w->active;

  cs_post_writer_times_t *ot = w->ot;

  /* If no output times list is provided, nothing to do */

  if (ot == NULL)
    return;

  /* In case of previous calls for a given time step,
     do not change status (which must have been forced otherwise),
     but update lists so as not to provoke an output at the next
     time step (so as to be consistent with the forcing that must have
     been done prior to entering here for this situation to exist). */

  if (w->n_last == ts->nt_cur)
    force_status = true;

  /* Test for listed time steps */

  i = 0;
  while (i < ot->n_t_steps) {
    /* Activate, then remove current or previous time steps from list */
    if (ot->t_steps[i] <= ts->nt_cur) {
      if (w->active > -1)
        w->active = 1;
      ot->t_steps[i] = ot->t_steps[ot->n_t_steps - 1];
      ot->n_t_steps -= 1;
    }
    else
      i++;
  }

  /* Test for listed time values */

  i = 0;
  while (i < ot->n_t_vals) {
    /* Activate, then remove current or previous time values from list */
    if (ot->t_vals[i] <= ts->t_cur) {
      if (w->active > -1)
        w->active = 1;
      ot->t_vals[i] = ot->t_vals[ot->n_t_steps - 1];
      ot->n_t_vals -= 1;
    }
    else
      i++;
  }

  if (force_status)
    w->active = prev_status;
}

/*----------------------------------------------------------------------------
 * Convert cs_post_type_t datatype to cs_datatype_t.
 *
 * parameters:
 *   type_cs <-- Code_Saturne data type
 *
 * returns
 *   corresponding FVM datatype
 *----------------------------------------------------------------------------*/

static cs_datatype_t
_cs_post_cnv_datatype(cs_post_type_t  type_cs)
{
  cs_datatype_t type_fvm = CS_DATATYPE_NULL;

  switch(type_cs) {

  case CS_POST_TYPE_cs_int_t:
    if (sizeof(cs_int_t) == 4)
      type_fvm = CS_INT32;
    else if (sizeof(cs_int_t) == 8)
      type_fvm = CS_INT64;
    break;

  case CS_POST_TYPE_cs_real_t:
    if (sizeof(cs_real_t) == sizeof(double))
      type_fvm = CS_DOUBLE;
    else if (sizeof(cs_real_t) == sizeof(float))
      type_fvm = CS_FLOAT;
    break;

  case CS_POST_TYPE_int:
    if (sizeof(int) == 4)
      type_fvm = CS_INT32;
    else if (sizeof(int) == 8)
      type_fvm = CS_INT64;
    break;

  case CS_POST_TYPE_float:
    type_fvm = CS_FLOAT;
    break;

  case CS_POST_TYPE_double:
    type_fvm = CS_DOUBLE;
    break;

  default:
    assert(0);
  }

  return type_fvm;
}

/*----------------------------------------------------------------------------
 * Search for position in the array of writers of a writer with a given id.
 *
 * parameters:
 *   writer_id <-- id of writer
 *
 * returns:
 *   position in the array of writers
 *----------------------------------------------------------------------------*/

static int
_cs_post_writer_id(const int  writer_id)
{
  int  id;

  cs_post_writer_t  *writer = NULL;

  /* Search for requested writer */

  for (id = 0; id < _cs_post_n_writers; id++) {
    writer = _cs_post_writers + id;
    if (writer->id == writer_id)
      break;
  }
  if (id >= _cs_post_n_writers)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested post-processing writer number\n"
                "%d is not defined.\n"), (int)(writer_id));

  return id;
}

/*----------------------------------------------------------------------------
 * Search for position in the array of writers of a writer with a given id,
 * allowing the writer not to be present.
 *
 * parameters:
 *   writer_id <-- id of writer
 *
 * returns:
 *   position in the array of writers, or -1
 *----------------------------------------------------------------------------*/

static int
_cs_post_writer_id_try(const int  writer_id)
{
  int  id;

  cs_post_writer_t  *writer = NULL;

  /* Search for requested writer */

  for (id = 0; id < _cs_post_n_writers; id++) {
    writer = _cs_post_writers + id;
    if (writer->id == writer_id)
      break;
  }
  if (id >= _cs_post_n_writers)
    id = -1;

  return id;
}

/*----------------------------------------------------------------------------
 * Search for position in the array of meshes of a mesh with a given id.
 *
 * parameters:
 *   mesh_id <-- id of mesh
 *
 * returns:
 *   position in the array of meshes
 *----------------------------------------------------------------------------*/

static int
_cs_post_mesh_id(int  mesh_id)
{
  int id;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_meshes; id++) {
    post_mesh = _cs_post_meshes + id;
    if (post_mesh->id == mesh_id)
      break;
  }
  if (id >= _cs_post_n_meshes)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested post-processing mesh number\n"
                "%d is not defined.\n"), (int)mesh_id);

  return id;
}

/*----------------------------------------------------------------------------
 * Search for position in the array of meshes of a mesh with a given id,
 * allowing the id not to be present
 *
 * parameters:
 *   mesh_id <-- id of mesh
 *
 * returns:
 *   position in the array of meshes, or -1
 *----------------------------------------------------------------------------*/

static int
_cs_post_mesh_id_try(int  mesh_id)
{
  int id;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_meshes; id++) {
    post_mesh = _cs_post_meshes + id;
    if (post_mesh->id == mesh_id)
      break;
  }
  if (id >= _cs_post_n_meshes)
    id = -1;

  return id;
}

/*----------------------------------------------------------------------------
 * Return indicator base on Lagrangian calculation status:
 *
 * parameters:
 *   ts            <-- time step structure, or NULL
 *
 * returns:
 *   0 if Lagrangian model is not active
 *   1 if Lagrangian model is active but no particle data is ready
 *   2 if current but not previous particle data is present
 *   3 if current and previous particle data is present
 *----------------------------------------------------------------------------*/

static int
_lagrangian_needed(const cs_time_step_t  *ts)
{
  int retval = 0;

  int _model = cs_lagr_model_type();

  if (_model != 0) {

    retval = 1;

    if (ts != NULL) {
      int _restart = cs_lagr_particle_restart();
      int _nt_start = (_restart) ? ts->nt_prev : ts->nt_prev + 1;
      if (ts->nt_cur == _nt_start)
        retval = 2;
      else if (ts->nt_cur > _nt_start)
        retval = 3;
    }

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Update mesh attributes related to writer association.
 *
 * parameters:
 *   post_mesh  <-> pointer to postprocessing mesh
 *----------------------------------------------------------------------------*/

static void
_update_mesh_writer_associations(cs_post_mesh_t  *post_mesh)
{
  int i, j;

  /* Minimum and maximum time dependency flags initially inverted,
     will be recalculated after mesh - writer associations */

  if (post_mesh->time_varying)
    post_mesh->mod_flag_min = FVM_WRITER_TRANSIENT_CONNECT;
  else
    post_mesh->mod_flag_min = _cs_post_mod_flag_min;
  post_mesh->mod_flag_max = FVM_WRITER_FIXED_MESH;

  int   n_writers = post_mesh->n_writers;

  if (post_mesh->ent_flag[3] == 0) { /* Non-Lagrangian mesh */

    for (i = 0; i < n_writers; i++) {

      fvm_writer_time_dep_t mod_flag;
      const int _writer_id = post_mesh->writer_id[i];
      cs_post_writer_t  *writer = _cs_post_writers + _writer_id;

      if (writer->wd != NULL)
        mod_flag = writer->wd->time_dep;
      else
        mod_flag = fvm_writer_get_time_dep(writer->writer);

      if (mod_flag < post_mesh->mod_flag_min)
        post_mesh->mod_flag_min = mod_flag;
      if (mod_flag > post_mesh->mod_flag_max)
        post_mesh->mod_flag_max = mod_flag;

    }

  }
  else { /* Lagrangian mesh: post_mesh->ent_flag[3] != 0 */

    int mode = post_mesh->ent_flag[3];
    fvm_writer_time_dep_t mod_type = (mode == 2) ?
      FVM_WRITER_FIXED_MESH : FVM_WRITER_TRANSIENT_CONNECT;

    post_mesh->mod_flag_min = FVM_WRITER_TRANSIENT_CONNECT;
    post_mesh->mod_flag_max = FVM_WRITER_TRANSIENT_CONNECT;

    for (i = 0, j = 0; i < n_writers; i++) {

      fvm_writer_time_dep_t mod_flag;
      const int _writer_id = post_mesh->writer_id[i];
      cs_post_writer_t  *writer = _cs_post_writers + _writer_id;

      if (writer->wd != NULL)
        mod_flag = writer->wd->time_dep;
      else
        mod_flag = fvm_writer_get_time_dep(writer->writer);

      if (mod_flag == mod_type)
        post_mesh->writer_id[j++] = _writer_id;

    }

    if (j < n_writers) {
      post_mesh->n_writers = j;
      BFT_REALLOC(post_mesh->writer_id, j, int);
    }

  }
}

/*----------------------------------------------------------------------------
 * Add or select a post-processing mesh, do basic initialization, and return
 * a pointer to the associated structure.
 *
 * parameters:
 *   mesh_id      <-- requested mesh id
 *   time_varying <-- if true, mesh may be redefined over time if associated
 *                    writers allow it
 *   mode         <-- 0 for standard mesh, 1 for particles, 2 for trajectories,
 *                    3 for probes, 4 for profiles
 *   n_writers    <-- number of associated writers
 *   writer_ids   <-- ids of associated writers
 *
 * returns:
 *   pointer to associated structure
 *----------------------------------------------------------------------------*/

static cs_post_mesh_t *
_predefine_mesh(int        mesh_id,
                bool       time_varying,
                int        mode,
                int        n_writers,
                const int  writer_ids[])
{
  /* local variables */

  int  i, j;

  cs_post_mesh_t  *post_mesh = NULL;

  /* Check that the requested mesh is available */

  if (mesh_id == 0)
      bft_error(__FILE__, __LINE__, 0,
                _("The requested post-processing mesh number\n"
                  "must be < 0 (reserved) or > 0 (user).\n"));

  for (i = 0; i < _cs_post_n_meshes; i++) {
    if ((_cs_post_meshes + i)->id == mesh_id) {

      post_mesh = _cs_post_meshes + i;

      BFT_FREE(post_mesh->name);
      for (j = 0; j < 5; j++)
        BFT_FREE(post_mesh->criteria[j]);
      BFT_FREE(post_mesh->writer_id);

      post_mesh->exp_mesh = NULL;
      if (post_mesh->_exp_mesh != NULL)
        post_mesh->_exp_mesh = fvm_nodal_destroy(post_mesh->_exp_mesh);

      break;

    }
  }

  if (i == _cs_post_n_meshes) {

    /* Resize global array of exportable meshes */

    if (_cs_post_n_meshes == _cs_post_n_meshes_max) {
      if (_cs_post_n_meshes_max == 0)
        _cs_post_n_meshes_max = 8;
      else
        _cs_post_n_meshes_max *= 2;
      BFT_REALLOC(_cs_post_meshes,
                  _cs_post_n_meshes_max,
                  cs_post_mesh_t);
    }

    post_mesh = _cs_post_meshes + i;

    _cs_post_n_meshes += 1;
  }

  if (mesh_id < _cs_post_min_mesh_id)
    _cs_post_min_mesh_id = mesh_id;

  /* Assign newly created mesh to the structure */

  post_mesh->id = mesh_id;
  post_mesh->name = NULL;
  post_mesh->cat_id = mesh_id;
  post_mesh->edges_ref = -1;
  post_mesh->locate_ref = -1;

  post_mesh->n_writers = 0;
  post_mesh->writer_id = NULL;

  post_mesh->nt_last = -2;

  post_mesh->add_groups = false;
  if (post_mesh->id == -1 || post_mesh->id == -2)
    post_mesh->post_domain = true;
  else
    post_mesh->post_domain = false;

  post_mesh->time_varying = time_varying;

  for (j = 0; j < 5; j++) {
    post_mesh->criteria[j] = NULL;
    post_mesh->sel_func[j] = NULL;
    post_mesh->sel_input[j] = NULL;
    post_mesh->ent_flag[j] = 0;
  }

  post_mesh->n_i_faces = 0;
  post_mesh->n_b_faces = 0;

  post_mesh->density = 1.;

  post_mesh->exp_mesh = NULL;
  post_mesh->_exp_mesh = NULL;

  /* Minimum and maximum time dependency flags initially inverted,
     will be recalculated after mesh - writer associations */

  if (post_mesh->time_varying)
    post_mesh->mod_flag_min = FVM_WRITER_TRANSIENT_CONNECT;
  else
    post_mesh->mod_flag_min = _cs_post_mod_flag_min;
  post_mesh->mod_flag_max = FVM_WRITER_FIXED_MESH;

  /* Associate mesh with writers */

  post_mesh->n_writers = n_writers;
  BFT_MALLOC(post_mesh->writer_id, n_writers, int);

  for (i = 0; i < n_writers; i++)
    post_mesh->writer_id[i] = _cs_post_writer_id(writer_ids[i]);

  if (mode == 1 || mode == 2)          /* Lagrangian mesh */
    post_mesh->ent_flag[3] = mode;

  else if (mode == 3 || mode == 4)     /* Probe or profile mesh */
    post_mesh->ent_flag[4] = mode - 2; /* 1 = probe monitoring,
                                          2 = profile */

  _update_mesh_writer_associations(post_mesh);

  return post_mesh;
}

/*----------------------------------------------------------------------------
 * Free a postprocessing mesh's data.
 *
 * parameters:
 *   _mesh_id <-- local id of mesh to remove
 *----------------------------------------------------------------------------*/

static void
_free_mesh(int _mesh_id)
{
  int i;
  cs_post_mesh_t  *post_mesh = _cs_post_meshes + _mesh_id;

  if (post_mesh->_exp_mesh != NULL)
    post_mesh->_exp_mesh = fvm_nodal_destroy(post_mesh->_exp_mesh);

  BFT_FREE(post_mesh->writer_id);
  post_mesh->n_writers = 0;

  for (i = 0; i < 4; i++)
    BFT_FREE(post_mesh->criteria[i]);

  BFT_FREE(post_mesh->name);

  /* Shift remaining meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->locate_ref > _mesh_id)
      post_mesh->locate_ref -= 1;
    else if (post_mesh->locate_ref == _mesh_id)
      post_mesh->locate_ref = -1;
    if (post_mesh->edges_ref >= _mesh_id) {
      assert(post_mesh->edges_ref != _mesh_id);
      post_mesh->edges_ref -= 1;
    }
  }

  for (i = _mesh_id + 1; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    _cs_post_meshes[i-1] = _cs_post_meshes[i];
  }
  _cs_post_n_meshes -= 1;
}

/*----------------------------------------------------------------------------
 * Check and possibly fix postprocessing mesh category id once we have
 * knowledge of entity counts.
 *
 * parameters:
 *   post_mesh <-> pointer to partially initialized post-processing mesh
 *----------------------------------------------------------------------------*/

static void
_check_mesh_cat_id(cs_post_mesh_t  *post_mesh)
{
  if (   post_mesh->cat_id == CS_POST_MESH_VOLUME
      || post_mesh->cat_id == CS_POST_MESH_BOUNDARY) {
    const int *ef = post_mesh->ent_flag;
    if (ef[0] == 1 && ef[1] == 0 && ef[2] == 0)
      post_mesh->cat_id = CS_POST_MESH_VOLUME;
    else if (ef[0] == 0 && ef[1] == 0 && ef[2] == 1)
      post_mesh->cat_id = CS_POST_MESH_BOUNDARY;
  }
}

/*----------------------------------------------------------------------------
 * Create a post-processing mesh; lists of cells or faces to extract are
 * sorted upon exit, whether they were sorted upon calling or not.
 *
 * The list of associated cells is only necessary if the number of cells
 * to extract is strictly greater than 0 and less than the number of cells
 * of the computational mesh.
 *
 * Lists of faces are ignored if the number of extracted cells is nonzero;
 * otherwise, if the number of boundary faces to extract is equal to the
 * number of boundary faces in the computational mesh, and the number of
 * interior faces to extract is zero, then we extract by default the boundary
 * mesh, and the list of associated boundary faces is thus not necessary.
 *
 * parameters:
 *   post_mesh   <-> pointer to partially initialized post-processing mesh
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-> list of associated cells
 *   i_face_list <-> list of associated interior faces
 *   b_face_list <-> list of associated boundary faces
 *----------------------------------------------------------------------------*/

static void
_define_export_mesh(cs_post_mesh_t  *post_mesh,
                    cs_lnum_t        n_cells,
                    cs_lnum_t        n_i_faces,
                    cs_lnum_t        n_b_faces,
                    cs_lnum_t        cell_list[],
                    cs_lnum_t        i_face_list[],
                    cs_lnum_t        b_face_list[])
{
  /* local variables */

  fvm_nodal_t  *exp_mesh = NULL;

  /* Create associated structure */

  if (post_mesh->ent_flag[0] == 1) {

    if (n_cells >= cs_glob_mesh->n_cells)
      exp_mesh = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                post_mesh->name,
                                                post_mesh->add_groups,
                                                cs_glob_mesh->n_cells,
                                                NULL);
    else
      exp_mesh = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                post_mesh->name,
                                                post_mesh->add_groups,
                                                n_cells,
                                                cell_list);

  }
  else {

    if (   n_b_faces >= cs_glob_mesh->n_b_faces
        && n_i_faces == 0)
      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                post_mesh->name,
                                                post_mesh->add_groups,
                                                0,
                                                cs_glob_mesh->n_b_faces,
                                                NULL,
                                                NULL);
    else
      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                post_mesh->name,
                                                post_mesh->add_groups,
                                                n_i_faces,
                                                n_b_faces,
                                                i_face_list,
                                                b_face_list);

  }

  /* Fix category id now that we have knowledge of entity counts */

  _check_mesh_cat_id(post_mesh);

  /* Local dimensions */

  post_mesh->n_i_faces = n_i_faces;
  post_mesh->n_b_faces = n_b_faces;

  /* As faces might be split, ensure the number of faces is correct */

  /* Link to newly created mesh */

  post_mesh->exp_mesh = exp_mesh;
  post_mesh->_exp_mesh = exp_mesh;
}

/*----------------------------------------------------------------------------
 * Create a particles post-processing mesh;
 *
 * parameters:
 *   post_mesh     <-> pointer to partially initialized post-processing mesh
 *   n_particles   <-- number of associated particles
 *   particle_list <-> list of associated particles
 *   ts            <-- time step structure
 *----------------------------------------------------------------------------*/

static void
_define_particle_export_mesh(cs_post_mesh_t        *post_mesh,
                             cs_lnum_t              n_particles,
                             cs_lnum_t              particle_list[],
                             const cs_time_step_t  *ts)
{
  /* local variables */

  fvm_nodal_t  *exp_mesh = NULL;

  assert(ts != NULL);

  /* Create associated structure */

  {
    cs_gnum_t *global_num = NULL;
    cs_coord_3_t *coords = NULL;
    fvm_io_num_t  *io_num = NULL;

    cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();

    /* Particle positions */

    if (post_mesh->ent_flag[3] == 1) {

      assert(ts->nt_cur > -1);

      exp_mesh = fvm_nodal_create(post_mesh->name, 3);

      BFT_MALLOC(coords, n_particles, cs_coord_3_t);

      cs_lagr_get_particle_values(p_set,
                                  CS_LAGR_COORDS,
                                  CS_REAL_TYPE,
                                  3,
                                  -1,
                                  n_particles,
                                  particle_list,
                                  coords);

      fvm_nodal_define_vertex_list(exp_mesh, n_particles, NULL);
      fvm_nodal_transfer_vertices(exp_mesh, (cs_coord_t *)coords);

    }

    /* Particle trajectories */

    else if (post_mesh->ent_flag[3] == 2) {

      cs_lnum_t i;
      cs_lnum_t  *vertex_num;
      char *mesh_name;

      assert(ts->nt_cur > 0);

      BFT_MALLOC(mesh_name, strlen(post_mesh->name) + 32, char);
      sprintf(mesh_name, "%s_%05d", post_mesh->name, ts->nt_cur);

      exp_mesh = fvm_nodal_create(mesh_name, 3);

      BFT_FREE(mesh_name);

      BFT_MALLOC(vertex_num, n_particles*2, cs_lnum_t);

      for (i = 0; i < n_particles*2; i++)
        vertex_num[i] = i+1;

      BFT_MALLOC(coords, n_particles*2, cs_coord_3_t);

      cs_lagr_get_trajectory_values(p_set,
                                    CS_LAGR_COORDS,
                                    CS_REAL_TYPE,
                                    3,
                                    -1,
                                    n_particles,
                                    particle_list,
                                    coords);

      fvm_nodal_append_by_transfer(exp_mesh,
                                   n_particles,
                                   FVM_EDGE,
                                   NULL,
                                   NULL,
                                   NULL,
                                   vertex_num,
                                   NULL);

      fvm_nodal_transfer_vertices(exp_mesh, (cs_coord_t *)coords);

      if (post_mesh->nt_last < ts->nt_cur)
        post_mesh->nt_last = -2;
    }

    /* Build global numbering if required */

    if (_number_particles_by_coord)
      io_num = fvm_io_num_create_from_sfc((const cs_coord_t *)coords,
                                          3,
                                          n_particles,
                                          FVM_IO_NUM_SFC_MORTON_BOX);
    else if (cs_glob_n_ranks > 1)
      io_num = fvm_io_num_create_from_scan(n_particles);

    if (io_num != NULL) {

      global_num = fvm_io_num_transfer_global_num(io_num);
      fvm_io_num_destroy(io_num);

      if (post_mesh->ent_flag[3] == 1) {

        fvm_nodal_init_io_num(exp_mesh, global_num, 0);
        BFT_FREE(global_num);

      }
      else if (post_mesh->ent_flag[3] == 2) {

        cs_lnum_t i;
        cs_gnum_t *g_coord_num;

        fvm_nodal_init_io_num(exp_mesh, global_num, 1);

        BFT_MALLOC(g_coord_num, n_particles*2, cs_gnum_t);
        for (i = 0; i < n_particles; i++) {
          g_coord_num[i*2] = global_num[i]*2 - 1;
          g_coord_num[i*2+1] = global_num[i]*2;
        }
        BFT_FREE(global_num);

        fvm_nodal_init_io_num(exp_mesh, g_coord_num, 0);

        BFT_FREE(g_coord_num);

      }

    }

    /* Drop trajectory sub-mesh if it is empty
       (otherwise, EnSight is OK, but ParaView can't display variables) */

    if (   post_mesh->ent_flag[3] == 2
        && fvm_nodal_get_n_g_elements(exp_mesh, FVM_EDGE) == 0)
      exp_mesh = fvm_nodal_destroy(exp_mesh);

  }

  /* Fix category id */

  if (post_mesh->cat_id < 0)
    post_mesh->cat_id = CS_POST_MESH_PARTICLES;

  /* Link to newly created mesh */

  post_mesh->exp_mesh = exp_mesh;
  post_mesh->_exp_mesh = exp_mesh;
}

/*----------------------------------------------------------------------------
 * Initialize a volume or surface post-processing mesh based on its
 * selection criteria or selection functions.
 *
 * parameters:
 *   post_mesh <-> pointer to partially initialized post-processing mesh
 *   ts        <-- time step structure
 *----------------------------------------------------------------------------*/

static void
_define_regular_mesh(cs_post_mesh_t  *post_mesh)
{
  const cs_mesh_t *mesh = cs_glob_mesh;

  assert(post_mesh != NULL);

  assert(post_mesh->exp_mesh == NULL);

  cs_lnum_t i;
  cs_lnum_t n_cells = 0, n_i_faces = 0, n_b_faces = 0;
  cs_lnum_t *cell_list = NULL, *i_face_list = NULL, *b_face_list = NULL;

  /* Define element lists based on selection criteria */

  if (post_mesh->criteria[0] != NULL) {
    const char *criteria = post_mesh->criteria[0];
    if (!strcmp(criteria, "all[]"))
      n_cells = mesh->n_cells;
    else {
      BFT_MALLOC(cell_list, mesh->n_cells, cs_lnum_t);
      cs_selector_get_cell_num_list(criteria, &n_cells, cell_list);
    }
  }
  else if (post_mesh->sel_func[0] != NULL) {
    cs_post_elt_select_t *sel_func = post_mesh->sel_func[0];
    sel_func(post_mesh->sel_input[0], &n_cells, &cell_list);
    for (i = 0; i < n_cells; i++)
      cell_list[i] += 1;
  }

  if (post_mesh->criteria[1] != NULL) {
    const char *criteria = post_mesh->criteria[1];
    if (!strcmp(criteria, "all[]"))
      n_i_faces = mesh->n_i_faces;
    else {
      BFT_MALLOC(i_face_list, mesh->n_i_faces, cs_lnum_t);
      cs_selector_get_i_face_num_list(criteria, &n_i_faces, i_face_list);
    }
  }
  else if (post_mesh->sel_func[1] != NULL) {
    cs_post_elt_select_t *sel_func = post_mesh->sel_func[1];
    sel_func(post_mesh->sel_input[1], &n_i_faces, &i_face_list);
    for (i = 0; i < n_i_faces; i++)
      i_face_list[i] += 1;
  }

  if (post_mesh->criteria[2] != NULL) {
    const char *criteria = post_mesh->criteria[2];
    if (!strcmp(criteria, "all[]"))
      n_b_faces = mesh->n_b_faces;
    else {
      BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);
      cs_selector_get_b_face_num_list(criteria, &n_b_faces, b_face_list);
    }
  }
  else if (post_mesh->sel_func[2] != NULL) {
    cs_post_elt_select_t *sel_func = post_mesh->sel_func[2];
    sel_func(post_mesh->sel_input[2], &n_b_faces, &b_face_list);
    for (i = 0; i < n_b_faces; i++)
      b_face_list[i] += 1;
  }

  /* Define mesh based on current arguments */

  _define_export_mesh(post_mesh,
                      n_cells,
                      n_i_faces,
                      n_b_faces,
                      cell_list,
                      i_face_list,
                      b_face_list);

  BFT_FREE(cell_list);
  BFT_FREE(i_face_list);
  BFT_FREE(b_face_list);
}

/*----------------------------------------------------------------------------
 * Create a post-processing mesh for probes
 *
 * parameters:
 *   post_mesh     <-> pointer to partially initialized post-processing mesh
 *----------------------------------------------------------------------------*/

static void
_define_probe_export_mesh(cs_post_mesh_t  *post_mesh)
{
  /* Sanity checks */
  assert(post_mesh != NULL);

  cs_probe_set_t     *pset = (cs_probe_set_t *)post_mesh->sel_input[4];
  fvm_nodal_t        *exp_mesh = NULL;
  const fvm_nodal_t  *location_mesh = NULL;

  /* First step: locate probes and update their coordinates */

  if (post_mesh->locate_ref > -1) {
    cs_post_mesh_t *post_mesh_loc = _cs_post_meshes + post_mesh->locate_ref;
    if (post_mesh_loc->exp_mesh == NULL)
      _define_regular_mesh(post_mesh_loc);
    location_mesh = post_mesh_loc->exp_mesh;
  }

  cs_probe_set_locate(pset, location_mesh);

  /* Create associated structure */

  exp_mesh = cs_probe_set_export_mesh(pset, cs_probe_set_get_name(pset));

  /* Link to newly created mesh */

  post_mesh->exp_mesh = exp_mesh;
  post_mesh->_exp_mesh = exp_mesh;

  /* Unassign matching location mesh ids for non-transient probe sets
     to allow freeing them */

  bool time_varying, is_profile, on_boundary, auto_variables;

  int  n_writers = 0;
  int  *writer_id = NULL;

  cs_probe_set_get_post_info(pset,
                             &time_varying,
                             &on_boundary,
                             &is_profile,
                             &auto_variables,
                             &n_writers,
                             &writer_id);

  if (time_varying == false)
    post_mesh->locate_ref = -1;
}

/*----------------------------------------------------------------------------
 * Initialize a post-processing mesh based on its selection criteria
 * or selection functions.
 *
 * parameters:
 *   post_mesh <-> pointer to partially initialized post-processing mesh
 *   ts        <-- time step structure
 *----------------------------------------------------------------------------*/

static void
_define_mesh(cs_post_mesh_t        *post_mesh,
             const cs_time_step_t  *ts)

{
  const cs_mesh_t *mesh = cs_glob_mesh;

  assert(post_mesh != NULL);

  assert(post_mesh->exp_mesh == NULL);

  /* Edges mesh */

  if (post_mesh->edges_ref > -1) {

    fvm_nodal_t *exp_edges = NULL;
    cs_post_mesh_t *post_base
      = _cs_post_meshes + _cs_post_mesh_id(post_mesh->edges_ref);

    /* if base mesh structure is not built yet, force its build now */

    if (post_base->exp_mesh == NULL)
      _define_mesh(post_base, ts);

    /* Copy mesh edges to new mesh structure */

    exp_edges = fvm_nodal_copy_edges(post_mesh->name, post_mesh->exp_mesh);

    /* Create mesh and assign to structure */

    post_mesh->exp_mesh = exp_edges;
    post_mesh->_exp_mesh = exp_edges;
  }

  /* Particle (Lagrangian) mesh */

  else if (post_mesh->ent_flag[3] != 0 && ts != NULL) {

    cs_lnum_t n_post_particles = 0, n_particles = cs_lagr_get_n_particles();
    cs_lnum_t *particle_list = NULL;

    if (post_mesh->criteria[3] != NULL) {

      cs_lnum_t n_cells = 0;
      cs_lnum_t *cell_list = NULL;
      const char *criteria = post_mesh->criteria[3];

      if (!strcmp(criteria, "all[]"))
        n_cells = mesh->n_cells;
      else {
        BFT_MALLOC(cell_list, mesh->n_cells, cs_lnum_t);
        cs_selector_get_cell_num_list(criteria, &n_cells, cell_list);
      }
      if (n_cells < mesh->n_cells || post_mesh->density < 1.) {
        BFT_MALLOC(particle_list, n_particles, cs_lnum_t);
        cs_lagr_get_particle_list(n_cells,
                                  cell_list,
                                  post_mesh->density,
                                  &n_post_particles,
                                  particle_list);
        BFT_REALLOC(particle_list, n_post_particles, cs_lnum_t);
      }
      else
        n_post_particles = n_particles;
      BFT_FREE(cell_list);
    }

    else if (post_mesh->sel_func[3] != NULL) {
      cs_post_elt_select_t *sel_func = post_mesh->sel_func[3];
      sel_func(post_mesh->sel_input[0], &n_post_particles, &particle_list);
    }

    _define_particle_export_mesh(post_mesh,
                                 n_post_particles,
                                 particle_list,
                                 ts);

    BFT_FREE(particle_list);
  }

  /* Probe mesh */

  else if (post_mesh->ent_flag[4] != 0) {

    _define_probe_export_mesh(post_mesh);

  }

  /* Standard (non-particle) meshes */

  else
    _define_regular_mesh(post_mesh);
}

/*----------------------------------------------------------------------------
 * Modify an existing post-processing mesh.
 *
 * It is not necessary to use this function if a mesh is simply deformed.
 *
 * parameters:
 *   post_mesh <-- pointer to postprocessing mesh structure
 *   ts        <-- time step structure
 *----------------------------------------------------------------------------*/

static void
_redefine_mesh(cs_post_mesh_t        *post_mesh,
               const cs_time_step_t  *ts)
{
  /* Remove previous base structure (return if we do not own the mesh) */

  if (post_mesh->exp_mesh != NULL) {
    if (post_mesh->_exp_mesh == NULL)
      return;
    else
      post_mesh->_exp_mesh = fvm_nodal_destroy(post_mesh->_exp_mesh);
  }
  post_mesh->exp_mesh = NULL;

  /* Define new mesh */

  _define_mesh(post_mesh, ts);
}

/*----------------------------------------------------------------------------
 * Remove meshes which are associated with no writer
 *----------------------------------------------------------------------------*/

static void
_clear_unused_meshes(void)
{
  int  i;
  int *discard = NULL;

  cs_post_mesh_t  *post_mesh;

  /* Mark used meshes */

  BFT_MALLOC(discard, _cs_post_n_meshes, int);

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->n_writers == 0)
      discard[i] = 1;
    else
      discard[i] = 0;
  }

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->locate_ref > -1) {
      if (post_mesh->n_writers > 0)
        discard[post_mesh->locate_ref] = 0;
    }
  }

  /* Discard meshes not required, compacting array */

  for (i = _cs_post_n_meshes - 1; i >= 0; i--) {
    if (discard[i] == 1)
      _free_mesh(i);  /* shifts other meshes and reduces _cs_post_n_meshes */
  }

  BFT_FREE(discard);
}

/*----------------------------------------------------------------------------
 * Divide polygons or polyhedra in simpler elements if necessary.
 *
 * parameters:
 *   post_mesh <-> pointer to post-processing mesh
 *   writer    <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

static void
_divide_poly(cs_post_mesh_t    *post_mesh,
             cs_post_writer_t  *writer)
{
  if (fvm_writer_needs_tesselation(writer->writer,
                                   post_mesh->exp_mesh,
                                   FVM_CELL_POLY) > 0)
    fvm_nodal_tesselate(post_mesh->_exp_mesh, FVM_CELL_POLY, NULL);

  if (fvm_writer_needs_tesselation(writer->writer,
                                   post_mesh->exp_mesh,
                                   FVM_FACE_POLY) > 0)
    fvm_nodal_tesselate(post_mesh->_exp_mesh, FVM_FACE_POLY, NULL);
}

/*----------------------------------------------------------------------------
 * Write parallel domain (rank) number to post-processing mesh
 *
 * parameters:
 *   writer     <-- FVM writer
 *   exp_mesh   <-- exportable mesh
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *---------------------------------------------------------------------------*/

static void
_cs_post_write_domain(fvm_writer_t       *writer,
                      const fvm_nodal_t  *exp_mesh,
                      int                 nt_cur_abs,
                      double              t_cur_abs)
{
  int  dim_ent;
  cs_lnum_t  i, n_elts;

  cs_lnum_t  parent_num_shift[1]  = {0};
  int32_t  *domain = NULL;

  int _nt_cur_abs = -1;
  double _t_cur_abs = 0.;

  const int32_t   *var_ptr[1] = {NULL};

  if (cs_glob_n_ranks < 2)
    return;

  dim_ent = fvm_nodal_get_max_entity_dim(exp_mesh);
  n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ent);

  /* Prepare domain number */

  BFT_MALLOC(domain, n_elts, int32_t);

  for (i = 0; i < n_elts; i++)
    domain[i] = cs_glob_rank_id;

  /* Prepare post-processing output */

  var_ptr[0] = domain;

  if (fvm_writer_get_time_dep(writer) != FVM_WRITER_FIXED_MESH) {
    _nt_cur_abs = nt_cur_abs;
    _t_cur_abs = t_cur_abs;
  }

  fvm_writer_export_field(writer,
                          exp_mesh,
                          "mpi_rank_id",
                          FVM_WRITER_PER_ELEMENT,
                          1,
                          CS_INTERLACE,
                          0,
                          parent_num_shift,
                          CS_INT32,
                          _nt_cur_abs,
                          _t_cur_abs,
                          (const void * *)var_ptr);

  /* Free memory */

  BFT_FREE(domain);
}

/*----------------------------------------------------------------------------
 * Output fixed zone information if needed.
 *
 * This function is called when we know some writers are active
 *
 * parameters:
 *   writer     <-- FVM writer
 *   post_mesh  <-- postprocessing mesh
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

static void
_cs_post_write_fixed_zone_info(fvm_writer_t          *writer,
                               const cs_post_mesh_t  *post_mesh,
                               int                    nt_cur_abs,
                               double                 t_cur_abs)
{
  assert(post_mesh->exp_mesh != NULL);

  bool output = false;
  const int   *var_ptr[1] = {NULL};

  if (post_mesh->id == CS_POST_MESH_VOLUME) {

    /* Ignore cases where all zones include all cells */

    int n_zones = cs_volume_zone_n_zones();
    int z_id = 0;

    for (z_id = 0; z_id < n_zones; z_id++) {
      const cs_volume_zone_t  *z = cs_volume_zone_by_id(z_id);
      if (z->location_id != CS_MESH_LOCATION_CELLS)
        break;
    }
    if (z_id >= n_zones)
      return;

    const int *zone_id = cs_volume_zone_cell_zone_id();

    if (cs_volume_zone_n_zones_time_varying() == 0) {
      output = true;
      var_ptr[0] = zone_id;
    }

  }

  else if (post_mesh->id == CS_POST_MESH_BOUNDARY) {

    /* Ignore cases where all zones include all boundary faces */

    int n_zones = cs_boundary_zone_n_zones();
    int z_id = 0;

    for (z_id = 0; z_id < n_zones; z_id++) {
      const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);
      if (z->location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
        break;
    }
    if (z_id >= n_zones)
      return;

    const int *zone_id = cs_boundary_zone_face_zone_id();

    if (cs_boundary_zone_n_zones_time_varying() == 0) {
      output = true;
      var_ptr[0] = zone_id;
    }

  }

  if (output) {

    cs_lnum_t  parent_num_shift[1]  = {0};
    int _nt_cur_abs = -1;
    double _t_cur_abs = 0.;

    if (fvm_writer_get_time_dep(writer) != FVM_WRITER_FIXED_MESH) {
      _nt_cur_abs = nt_cur_abs;
      _t_cur_abs = t_cur_abs;
    }

    fvm_writer_export_field(writer,
                            post_mesh->exp_mesh,
                            "boundary zone id",
                            FVM_WRITER_PER_ELEMENT,
                            1,
                            CS_INTERLACE,
                            1,
                            parent_num_shift,
                            CS_INT_TYPE,
                            _nt_cur_abs,
                            _t_cur_abs,
                            (const void * *)var_ptr);

  }
}

/*----------------------------------------------------------------------------
 * Output varying zone information if needed.
 *
 * This function is called when we know some writers are active
 *
 * parameters:
 *   ts <-- time step structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_post_write_transient_zone_info(const cs_post_mesh_t  *post_mesh,
                                   const cs_time_step_t  *ts)
{
  if (post_mesh->id == CS_POST_MESH_VOLUME) {
    if (cs_volume_zone_n_zones_time_varying() > 0) {
      cs_post_write_var(post_mesh->id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "volume zone id",
                        1,       /* var_dim */
                        true,    /* interlace */
                        true,    /* use_parent */
                        CS_POST_TYPE_int,
                        cs_volume_zone_cell_zone_id(),
                        NULL,
                        NULL,
                        ts);
    }
  }

  else if (post_mesh->id == CS_POST_MESH_BOUNDARY) {
    if (cs_boundary_zone_n_zones_time_varying() > 0)
      cs_post_write_var(post_mesh->id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "boundary zone id",
                        1,       /* var_dim */
                        true,    /* interlace */
                        true,    /* use_parent */
                        CS_POST_TYPE_int,
                        NULL,
                        NULL,
                        cs_boundary_zone_face_zone_id(),
                        ts);
  }
}

/*----------------------------------------------------------------------------
 * Output a post-processing mesh using associated writers.
 *
 * If the time step structure argument passed is NULL, a time-independent
 * output will be assumed.
 *
 * parameters:
 *   post_mesh  <-> pointer to post-processing mesh
 *   ts         <-- time step structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_post_write_mesh(cs_post_mesh_t        *post_mesh,
                    const cs_time_step_t  *ts)
{
  int  j;
  fvm_writer_time_dep_t  time_dep;
  bool  write_mesh = false;

  cs_post_writer_t *writer = NULL;

  const int nt_cur = (ts != NULL) ? ts->nt_cur : -1;
  const double t_cur = (ts != NULL) ? ts->t_cur : 0.;

  /* Special case: particle trajectories can not be
     output before time stepping starts */

  if (post_mesh->ent_flag[3] == 2 && nt_cur < 1)
    return;

  /* Loop on writers */

  for (j = 0; j < post_mesh->n_writers; j++) {

    writer = _cs_post_writers + post_mesh->writer_id[j];

    if (writer->wd != NULL)
      time_dep = writer->wd->time_dep;
    else
      time_dep = fvm_writer_get_time_dep(writer->writer);

    write_mesh = false;

    if (   time_dep == FVM_WRITER_FIXED_MESH
        && writer->active > -1
        && post_mesh->ent_flag[3] != 2) {
      if (post_mesh->nt_last < -1)
        write_mesh = true;
    }
    else {
      if (post_mesh->nt_last < nt_cur && writer->active == 1)
        write_mesh = true;
    }

    if (write_mesh == true) {

      if (writer->writer == NULL)
        _init_writer(writer);

      if (post_mesh->exp_mesh == NULL)
        _define_mesh(post_mesh, ts);

      if (post_mesh->exp_mesh == NULL)
        continue;

      _divide_poly(post_mesh, writer);

      fvm_writer_set_mesh_time(writer->writer, nt_cur, t_cur);
      fvm_writer_export_nodal(writer->writer, post_mesh->exp_mesh);

      if (nt_cur >= 0) {
        writer->n_last = nt_cur;
        writer->t_last = t_cur;
      }

    }

    if (write_mesh == true) {

      if (post_mesh->post_domain)
        _cs_post_write_domain(writer->writer,
                              post_mesh->exp_mesh,
                              nt_cur,
                              t_cur);

      _cs_post_write_fixed_zone_info(writer->writer,
                                     post_mesh,
                                     nt_cur,
                                     t_cur);

      if (nt_cur >= 0) {
        writer->n_last = nt_cur;
        writer->t_last = t_cur;
      }

    }

  }

  if (write_mesh == true)
    post_mesh->nt_last = nt_cur;
}

/*----------------------------------------------------------------------------
 * Assemble variable values defined on a mix of interior and boundary
 * faces (with no indirection) into an array defined on a single faces set.
 *
 * The resulting variable is not interlaced.
 *
 * parameters:
 *   exp_mesh    <-- exportable mesh
 *   n_i_faces   <-- number of interior faces
 *   n_b_faces   <-- number of boundary faces
 *   var_dim     <-- variable dimension
 *   interlace   <-- for vector, interlace if 1, no interlace if 0
 *   i_face_vals <-- values at interior faces
 *   b_face_vals <-- values at boundary faces
 *   var_tmp[]   --> assembled values
 *----------------------------------------------------------------------------*/

static void
_cs_post_assmb_var_faces(const fvm_nodal_t  *exp_mesh,
                         cs_lnum_t           n_i_faces,
                         cs_lnum_t           n_b_faces,
                         int                 var_dim,
                         cs_interlace_t      interlace,
                         const cs_real_t     i_face_vals[],
                         const cs_real_t     b_face_vals[],
                         cs_real_t           var_tmp[])
{
  cs_lnum_t  i, j, stride_1, stride_2;

  cs_lnum_t  n_elts = n_i_faces + n_b_faces;

  assert(exp_mesh != NULL);

  /* The variable is defined on interior and boundary faces of the
     post-processing mesh, and has been built using values
     at the corresponding interior and boundary faces */

  /* Boundary faces contribution */

  if (interlace == CS_INTERLACE) {
    stride_1 = var_dim;
    stride_2 = 1;
  }
  else {
    stride_1 = 1;
    stride_2 = n_b_faces;
  }

  for (i = 0; i < n_b_faces; i++) {
    for (j = 0; j < var_dim; j++)
      var_tmp[i + j*n_elts] = b_face_vals[i*stride_1 + j*stride_2];
  }

  /* Interior faces contribution */

  if (interlace == CS_INTERLACE) {
    stride_1 = var_dim;
    stride_2 = 1;
  }
  else {
    stride_1 = 1;
    stride_2 = n_i_faces;
  }

  for (i = 0; i < n_i_faces; i++) {
    for (j = 0; j < var_dim; j++)
      var_tmp[i + n_b_faces + j*n_elts] = i_face_vals[i*stride_1 + j*stride_2];
  }

}

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

static void
_cs_post_write_displacements(int     nt_cur_abs,
                             double  t_cur_abs)
{
  int i, j;
  cs_lnum_t  k, nbr_val;
  cs_datatype_t datatype;

  cs_lnum_t  parent_num_shift[1]  = {0};
  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_writer_t  *writer = NULL;
  cs_real_t  *deplacements = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_real_t   *var_ptr[1] = {NULL};

  /* Loop on writers to check if something must be done */
  /*----------------------------------------------------*/

  if (_cs_post_deformable == false)
    return;

  for (j = 0; j < _cs_post_n_writers; j++) {
    writer = _cs_post_writers + j;
    if (writer->active == 1)
      break;
  }
  if (j == _cs_post_n_writers)
    return;

  /* Compute main deformation field */
  /*--------------------------------*/

  nbr_val = mesh->n_vertices * 3;

  BFT_MALLOC(deplacements, nbr_val, cs_real_t);

  assert(mesh->n_vertices == 0 || _cs_post_ini_vtx_coo != NULL);

  for (k = 0; k < nbr_val; k++)
    deplacements[k] = mesh->vtx_coord[k] - _cs_post_ini_vtx_coo[k];

  /* Prepare post-processing */
  /*-------------------------*/

  if (sizeof(cs_real_t) == sizeof(double))
    datatype = CS_DOUBLE;
  else if (sizeof(cs_real_t) == sizeof(float))
    datatype = CS_FLOAT;

  var_ptr[0] = deplacements;

  /* Loop on meshes to output displacements */
  /*----------------------------------------*/

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    for (j = 0; j < post_mesh->n_writers; j++) {

      writer = _cs_post_writers + post_mesh->writer_id[j];

      if (writer->active == 1) {

        fvm_writer_export_field(writer->writer,
                                post_mesh->exp_mesh,
                                _("displacement"),
                                FVM_WRITER_PER_NODE,
                                3,
                                CS_INTERLACE,
                                1,
                                parent_num_shift,
                                datatype,
                                (int)nt_cur_abs,
                                (double)t_cur_abs,
                                (const void * *)var_ptr);

        if (nt_cur_abs >= 0) {
          writer->n_last = nt_cur_abs;
          writer->t_last = t_cur_abs;
        }

      }

    }

  }

  /* Free memory */

  BFT_FREE(deplacements);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if post-processing is activated and then update post-processing
 *        of meshes if there is a time-dependent mesh
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

static void
_update_meshes(const cs_time_step_t  *ts)
{
  bool  active;

  /* Loop on writers to check if something must be done */
  /*----------------------------------------------------*/

  int  w;
  for (w = 0; w < _cs_post_n_writers; w++) {
    cs_post_writer_t  *writer = _cs_post_writers + w;
    if (writer->active == 1)
      break;
  }
  if (w == _cs_post_n_writers)
    return;

  int t_top_id = cs_timer_stats_switch(_post_out_stat_id);

  const int nt_cur = (ts != NULL) ? ts->nt_cur : -1;
  const double t_cur = (ts != NULL) ? ts->t_cur : 0.;

  /* Possible modification of post-processing meshes */
  /*-------------------------------------------------*/

  for (int i = 0; i < _cs_post_n_meshes; i++) {

    cs_post_mesh_t  *post_mesh = _cs_post_meshes + i;

    active = false;

    for (int j = 0; j < post_mesh->n_writers; j++) {
      cs_post_writer_t  *writer = _cs_post_writers + post_mesh->writer_id[j];
      if (writer->active == 1)
        active = true;
    }

    /* Modifiable user mesh, active at this time step */

    if (   active == true
        && post_mesh->mod_flag_min == FVM_WRITER_TRANSIENT_CONNECT)
      _redefine_mesh(post_mesh, ts);

  }

  /* Output of meshes or vertex displacement field if necessary */
  /*------------------------------------------------------------*/

  cs_post_write_meshes(ts);

  if (_cs_post_deformable == true)
    _cs_post_write_displacements(nt_cur, t_cur);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Generate global group flags array from local family flags
 *
 * The caller should free the returned array once it is no longer needed.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   fam_flag <-- flag values (size: mesh->n_families + 1)
 *
 * returns:
 *   group flag (size: mesh->n_groups)
 *----------------------------------------------------------------------------*/

static char *
_build_group_flag(const cs_mesh_t  *mesh,
                  int              *fam_flag)
{
  int i, j;

  char *group_flag = NULL;

  BFT_MALLOC(group_flag, mesh->n_groups, char);
  memset(group_flag, 0, mesh->n_groups);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    int *_fam_flag = NULL;
    BFT_MALLOC(_fam_flag, mesh->n_families + 1, int);
    MPI_Allreduce(fam_flag, _fam_flag, mesh->n_families + 1,
                  MPI_INT, MPI_MAX, cs_glob_mpi_comm);
    memcpy(fam_flag, _fam_flag, (mesh->n_families + 1)*sizeof(int));
    BFT_FREE(_fam_flag);
  }
#endif /* defined(HAVE_MPI) */

  for (i = 0; i < mesh->n_families; i++) {
    if (fam_flag[(i+1)] != 0) {
      char mask = fam_flag[i+1];
      for (j = 0; j < mesh->n_max_family_items; j++) {
        int g_id = - mesh->family_item[mesh->n_families*j + i] - 1;
        if (g_id >= 0)
          group_flag[g_id] = group_flag[g_id] | mask;
      }
    }
  }

  return group_flag;
}

/*----------------------------------------------------------------------------
 * Set a family flags array to 1 for families containg a given group,
 * and to 0 for others.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   g_id     <-- group id
 *   fam_flag --> flag values (size: mesh->n_families)
 *----------------------------------------------------------------------------*/

static void
_set_fam_flags(const cs_mesh_t  *mesh,
               int               g_id,
               int              *fam_flag)
{
  int j, k;
  memset(fam_flag, 0, mesh->n_families*sizeof(int));
  for (j = 0; j < mesh->n_families; j++) {
    for (k = 0; k < mesh->n_max_family_items; k++) {
      int _g_id = - mesh->family_item[mesh->n_families*k + j] - 1;
      if (_g_id == g_id)
        fam_flag[j] = 1;
    }
  }
}

/*----------------------------------------------------------------------------
 * Output volume sub-meshes by group
 *
 * parameters:
 *   mesh      <-- base mesh
 *   fmt_name  <-- format name
 *   fmt_opts  <-- format options
 *---------------------------------------------------------------------------*/

static void
_vol_submeshes_by_group(const cs_mesh_t  *mesh,
                        const char       *fmt_name,
                        const char       *fmt_opts)
{
  cs_lnum_t i, j;
  cs_lnum_t n_cells, n_i_faces, n_b_faces;
  char part_name[81];
  int max_null_family = 0;
  int *fam_flag = NULL;
  char *group_flag = NULL;
  cs_lnum_t *cell_list = NULL, *i_face_list = NULL, *b_face_list = NULL;
  fvm_writer_t *writer = NULL;
  fvm_nodal_t *exp_mesh = NULL;

  if (mesh->n_families == 0)
    return;

  /* Families should be sorted, so if a nonzero family is empty,
     it is family 1 */

  if (mesh->family_item[0] == 0)
    max_null_family = 1;

  if (mesh->n_families <= max_null_family)
    return;

  /* Get writer info */

  /* Default values */

  /* Create default writer */

  writer = fvm_writer_init("mesh_groups",
                           _cs_post_dirname,
                           fmt_name,
                           fmt_opts,
                           FVM_WRITER_FIXED_MESH);

  /* Now detect which groups may be referenced */

  BFT_MALLOC(fam_flag, (mesh->n_families + 1), int);
  memset(fam_flag, 0, (mesh->n_families + 1) * sizeof(int));

  if (mesh->cell_family != NULL) {
    for (i = 0; i < mesh->n_cells; i++)
      fam_flag[mesh->cell_family[i]]
        = fam_flag[mesh->cell_family[i]] | 1;
  }
  if (mesh->i_face_family != NULL) {
    for (i = 0; i < mesh->n_i_faces; i++)
      fam_flag[mesh->i_face_family[i]]
        = fam_flag[mesh->i_face_family[i]] | 2;
  }
  if (mesh->b_face_family != NULL) {
    for (i = 0; i < mesh->n_b_faces; i++)
      fam_flag[mesh->b_face_family[i]]
        = fam_flag[mesh->b_face_family[i]] | 4;
  }

  group_flag = _build_group_flag(mesh, fam_flag);

  /* Now extract volume elements by groups.
     Note that selector structures may not have been initialized yet,
     so to avoid issue, we use a direct selection here. */

  BFT_REALLOC(fam_flag, mesh->n_families, int);

  BFT_MALLOC(cell_list, mesh->n_cells, cs_lnum_t);

  for (i = 0; i < mesh->n_groups; i++) {

    if (group_flag[i] & '\1') {

      const char *g_name = mesh->group + mesh->group_idx[i];

      _set_fam_flags(mesh, i, fam_flag);

      for (j = 0, n_cells = 0; j < mesh->n_cells; j++) {
        int f_id = mesh->cell_family[j];
        if (f_id > 0 && fam_flag[f_id - 1])
          cell_list[n_cells++] = j + 1;
      }
      strcpy(part_name, "vol: ");
      strncat(part_name, g_name, 80 - strlen(part_name));
      exp_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                                part_name,
                                                false,
                                                n_cells,
                                                cell_list);

      if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_CELL_POLY) > 0)
        fvm_nodal_tesselate(exp_mesh, FVM_CELL_POLY, NULL);

      fvm_writer_set_mesh_time(writer, -1, 0);
      fvm_writer_export_nodal(writer, exp_mesh);

      exp_mesh = fvm_nodal_destroy(exp_mesh);
    }

  }

  /* Now export cells with no groups */

  if (mesh->cell_family != NULL) {
    for (j = 0, n_cells = 0; j < mesh->n_cells; j++) {
      if (mesh->cell_family[j] <= max_null_family)
        cell_list[n_cells++] = j + 1;
    }
  }
  else {
    for (j = 0, n_cells = 0; j < mesh->n_cells; j++)
      cell_list[n_cells++] = j + 1;
  }

  i = n_cells;
  cs_parall_counter_max(&i, 1);

  if (i > 0) {
    exp_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                              "vol: no_group",
                                              false,
                                              n_cells,
                                              cell_list);

    if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_CELL_POLY) > 0)
      fvm_nodal_tesselate(exp_mesh, FVM_CELL_POLY, NULL);

    fvm_writer_set_mesh_time(writer, -1, 0);
    fvm_writer_export_nodal(writer, exp_mesh);

    exp_mesh = fvm_nodal_destroy(exp_mesh);
  }

  BFT_FREE(cell_list);

  /* Now extract faces by groups */

  BFT_MALLOC(i_face_list, mesh->n_i_faces, cs_lnum_t);
  BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);

  for (i = 0; i < mesh->n_groups; i++) {

    if ((group_flag[i] & '\2') || (group_flag[i] & '\4')) {

      const char *g_name = mesh->group + mesh->group_idx[i];

      _set_fam_flags(mesh, i, fam_flag);

      n_i_faces = 0;
      if (mesh->i_face_family != NULL) {
        for (j = 0; j < mesh->n_i_faces; j++) {
          int f_id = mesh->i_face_family[j];
          if (f_id > 0 && fam_flag[f_id - 1])
            i_face_list[n_i_faces++] = j + 1;
        }
      }
      n_b_faces = 0;
      if (mesh->b_face_family != NULL) {
        for (j = 0; j < mesh->n_b_faces; j++) {
          int f_id = mesh->b_face_family[j];
          if (f_id > 0 && fam_flag[f_id - 1])
            b_face_list[n_b_faces++] = j + 1;
        }
      }

      strcpy(part_name, "surf: ");
      strncat(part_name, g_name, 80 - strlen(part_name));
      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                part_name,
                                                false,
                                                n_i_faces,
                                                n_b_faces,
                                                i_face_list,
                                                b_face_list);

      if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_FACE_POLY) > 0)
        fvm_nodal_tesselate(exp_mesh, FVM_FACE_POLY, NULL);

      fvm_writer_set_mesh_time(writer, -1, 0);
      fvm_writer_export_nodal(writer, exp_mesh);

      exp_mesh = fvm_nodal_destroy(exp_mesh);
    }

  }

  writer = fvm_writer_finalize(writer);

  BFT_FREE(b_face_list);
  BFT_FREE(i_face_list);

  BFT_FREE(fam_flag);
  BFT_FREE(group_flag);
}

/*----------------------------------------------------------------------------
 * Output boundary sub-meshes by group, if it contains multiple groups.
 *
 * parameters:
 *   mesh        <-- base mesh
 *   fmt_name    <-- format name
 *   fmt_opts    <-- format options
 *---------------------------------------------------------------------------*/

static void
_boundary_submeshes_by_group(const cs_mesh_t   *mesh,
                             const char        *fmt_name,
                             const char        *fmt_opts)
{
  cs_lnum_t i, j;
  cs_lnum_t n_b_faces;
  cs_gnum_t n_no_group = 0;
  int max_null_family = 0;
  int *fam_flag = NULL;
  char *group_flag = NULL;
  cs_lnum_t *b_face_list = NULL;
  fvm_writer_t *writer = NULL;
  fvm_nodal_t *exp_mesh = NULL;

  if (mesh->n_families == 0)
    return;

  /* Families should be sorted, so if a nonzero family is empty,
     it is family 1 */

  if (mesh->family_item[0] == 0)
    max_null_family = 1;

  if (mesh->n_families <= max_null_family)
    return;

  /* Check how many boundary faces belong to no group */

  if (mesh->b_face_family != NULL) {
    for (j = 0, n_b_faces = 0; j < mesh->n_b_faces; j++) {
      if (mesh->b_face_family[j] <= max_null_family)
        n_no_group += 1;
    }
  }
  else
    n_no_group = mesh->n_b_faces;

  cs_parall_counter(&n_no_group, 1);

  if (n_no_group == mesh->n_g_b_faces)
    return;

  /* Get writer info */

  /* Default values */

  /* Create default writer */

  writer = fvm_writer_init("boundary_groups",
                           _cs_post_dirname,
                           fmt_name,
                           fmt_opts,
                           FVM_WRITER_FIXED_MESH);

  /* Now detect which groups may be referenced */

  BFT_MALLOC(fam_flag, mesh->n_families + 1, int);
  memset(fam_flag, 0, (mesh->n_families + 1)*sizeof(int));

  if (mesh->b_face_family != NULL) {
    for (i = 0; i < mesh->n_b_faces; i++)
      fam_flag[mesh->b_face_family[i]] = 1;
  }

  group_flag = _build_group_flag(mesh, fam_flag);

  /* Now extract boundary faces by groups.
     Note that selector structures may not have been initialized yet,
     so to avoid issue, we use a direct selection here. */

  BFT_REALLOC(fam_flag, mesh->n_families, int);

  BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);

  for (i = 0; i < mesh->n_groups; i++) {

    if (group_flag[i] != 0) {

      const char *g_name = mesh->group + mesh->group_idx[i];

      _set_fam_flags(mesh, i, fam_flag);

      n_b_faces = 0;
      if (mesh->b_face_family != NULL) {
        for (j = 0; j < mesh->n_b_faces; j++) {
          int f_id = mesh->b_face_family[j];
          if (f_id > 0 && fam_flag[f_id - 1])
            b_face_list[n_b_faces++] = j + 1;
        }
      }

      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                g_name,
                                                false,
                                                0,
                                                n_b_faces,
                                                NULL,
                                                b_face_list);

      if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_FACE_POLY) > 0)
        fvm_nodal_tesselate(exp_mesh, FVM_FACE_POLY, NULL);

      fvm_writer_set_mesh_time(writer, -1, 0);
      fvm_writer_export_nodal(writer, exp_mesh);

      exp_mesh = fvm_nodal_destroy(exp_mesh);
    }

  }

  /* Output boundary faces belonging to no group */

  if (n_no_group > 0) {

    if (mesh->b_face_family != NULL) {
      for (j = 0, n_b_faces = 0; j < mesh->n_b_faces; j++) {
        if (mesh->b_face_family[j] <= max_null_family)
          b_face_list[n_b_faces++] = j + 1;
      }
    }
    else {
      for (j = 0, n_b_faces = 0; j < mesh->n_b_faces; j++)
        b_face_list[n_b_faces++] = j + 1;
    }

    exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                              "no_group",
                                              false,
                                              0,
                                              n_b_faces,
                                              NULL,
                                              b_face_list);

    if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_FACE_POLY) > 0)
      fvm_nodal_tesselate(exp_mesh, FVM_FACE_POLY, NULL);

    fvm_writer_set_mesh_time(writer, -1, 0);
    fvm_writer_export_nodal(writer, exp_mesh);

    exp_mesh = fvm_nodal_destroy(exp_mesh);
  }

  BFT_FREE(b_face_list);

  writer = fvm_writer_finalize(writer);

  BFT_FREE(fam_flag);
  BFT_FREE(group_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interpolate values from fields defined on cells or faces.
 *
 * \param[in, out]  input           pointer to optional (untyped) value
 *                                  or structure.
 * \param[in]       datatype        associated datatype
 * \param[in]       val_dim         dimension of data values
 * \param[in]       n_points        number of interpolation points
 * \param[in]       point_location  location of points in mesh elements
 * \param[in]       point_coords    point coordinates
 * \param[in]       location_vals   values at mesh location
 * \param[out]      point_vals      interpolated values at points
 */
/*----------------------------------------------------------------------------*/

static void
_field_interpolate(void                *input,
                   cs_datatype_t        datatype,
                   int                  val_dim,
                   cs_lnum_t            n_points,
                   const cs_lnum_t      point_location[],
                   const cs_real_3_t    point_coords[],
                   const void          *location_vals,
                   void                *point_vals)
{
  CS_UNUSED(input);

  /* int f_id = *((int *)input); */

  /* TODO add options here */

  cs_interpolate_from_location_p0(NULL,
                                  datatype,
                                  val_dim,
                                  n_points,
                                  point_location,
                                  point_coords,
                                  location_vals,
                                  point_vals);
}

/*----------------------------------------------------------------------------
 * Main post-processing output of variables.
 *
 * parameters:
 *   post_mesh   <-- pointer to post-processing mesh structure
 *   ts          <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_post_output_fields(cs_post_mesh_t        *post_mesh,
                       const cs_time_step_t  *ts)
{
  /* Output for cell and boundary meshes */
  /*-------------------------------------*/

  if (   post_mesh->cat_id == CS_POST_MESH_VOLUME
      || post_mesh->cat_id == CS_POST_MESH_BOUNDARY) {

    const char *name;

    cs_mesh_location_type_t  mesh_cat_type = CS_MESH_LOCATION_NONE;

    if (post_mesh->cat_id == CS_POST_MESH_VOLUME)
      mesh_cat_type = CS_MESH_LOCATION_CELLS;
    else if (post_mesh->cat_id == CS_POST_MESH_BOUNDARY)
      mesh_cat_type = CS_MESH_LOCATION_BOUNDARY_FACES;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of mesh for the generic postprocessing of"
                  " fields.\n"
                  " Requested mesh is neither volumic nor surfacic."));

    const int n_fields = cs_field_n_fields();
    const int vis_key_id = cs_field_key_id("post_vis");
    const int label_key_id = cs_field_key_id("label");

    /* Loop on fields */

    for (int f_id = 0; f_id < n_fields; f_id++) {

      bool  use_parent = true;

      const cs_field_t  *f = cs_field_by_id(f_id);

      const cs_mesh_location_type_t field_loc_type
         = cs_mesh_location_get_type(f->location_id);

      if (mesh_cat_type == CS_MESH_LOCATION_CELLS) {
        if (field_loc_type != CS_MESH_LOCATION_CELLS &&
            field_loc_type != CS_MESH_LOCATION_VERTICES)
          continue;
      }

      else if (mesh_cat_type == CS_MESH_LOCATION_BOUNDARY_FACES) {
        if (field_loc_type != CS_MESH_LOCATION_BOUNDARY_FACES &&
            field_loc_type != CS_MESH_LOCATION_VERTICES)
          continue;
      }

      if (! (cs_field_get_key_int(f, vis_key_id) & CS_POST_ON_LOCATION))
        continue;

      name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;

      if (   field_loc_type == CS_MESH_LOCATION_CELLS
          || field_loc_type == CS_MESH_LOCATION_BOUNDARY_FACES) {

        const cs_real_t *cell_val = NULL, *b_face_val = NULL;

        if (field_loc_type == CS_MESH_LOCATION_CELLS)
          cell_val = f->val;
        else /* if (field_loc_type == CS_MESH_LOCATION_BOUNDARY_FACES) */
          b_face_val = f->val;

        cs_post_write_var(post_mesh->id,
                          CS_POST_WRITER_ALL_ASSOCIATED,
                          name,
                          f->dim,
                          true,
                          use_parent,
                          CS_POST_TYPE_cs_real_t,
                          cell_val,
                          NULL,
                          b_face_val,
                          ts);

      }

      else if (field_loc_type == CS_MESH_LOCATION_VERTICES)
        cs_post_write_vertex_var(post_mesh->id,
                                 CS_POST_WRITER_ALL_ASSOCIATED,
                                 name,
                                 f->dim,
                                 true,
                                 use_parent,
                                 CS_POST_TYPE_cs_real_t,
                                 f->val,
                                 ts);

    } /* End of loop on fields */

  } /* End of main output for cell or boundary mesh or submesh */

  /* Output for probes and profile meshes */
  /*--------------------------------------*/

  else if (post_mesh->cat_id == CS_POST_MESH_PROBES) {

    const int n_fields = cs_field_n_fields();
    const int vis_key_id = cs_field_key_id("post_vis");
    const int label_key_id = cs_field_key_id("label");

    /* Loop on fields */

    for (int f_id = 0; f_id < n_fields; f_id++) {

      const cs_field_t  *f = cs_field_by_id(f_id);

      const cs_mesh_location_type_t field_loc_type
         = cs_mesh_location_get_type(f->location_id);

      if (   field_loc_type != CS_MESH_LOCATION_CELLS
          && field_loc_type == CS_MESH_LOCATION_BOUNDARY_FACES
          && field_loc_type != CS_MESH_LOCATION_VERTICES)
        continue;

      if (! (cs_field_get_key_int(f, vis_key_id) & CS_POST_MONITOR))
        continue;

      const char *name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;

      cs_post_write_probe_values(post_mesh->id,
                                 CS_POST_WRITER_ALL_ASSOCIATED,
                                 name,
                                 f->dim,
                                 CS_POST_TYPE_cs_real_t,
                                 f->location_id,
                                 _field_interpolate,
                                 &f_id,
                                 f->val,
                                 ts);

    } /* End of loop on fields */

  } /* End of main output for cell or boundary mesh or submesh */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a post-processing mesh for probes (set of probes should have
 *        been already defined)
 *
 * \param[in] mesh_id        id of mesh to define (<0 reserved, >0 for user)
 * \param[in] pset           pointer to a cs_probe_set_t structure
 * \param[in] time_varying   true if probe coords may change during computation
 * \param[in] is_profile     true if probe set is related to a profile
 * \param[in] on_boundary    true if probes are located on boundary
 * \param[in] auto_variable  true if the set of variables to post is predefined
 * \param[in] n_writers      number of associated writers
 * \param[in] writer_ids     ids of associated writers
 */
/*----------------------------------------------------------------------------*/

static void
_cs_post_define_probe_mesh(int                    mesh_id,
                           cs_probe_set_t        *pset,
                           bool                   time_varying,
                           bool                   is_profile,
                           bool                   on_boundary,
                           bool                   auto_variable,
                           int                    n_writers,
                           const int              writer_ids[])
{
  assert(pset != NULL); /* Sanity check */

  /* Common initializations */

  int  mode = (is_profile == true) ? 4 : 3;
  cs_post_mesh_t *post_mesh = _predefine_mesh(mesh_id, time_varying, mode,
                                              n_writers, writer_ids);

  /* Define mesh based on current arguments */

  const char  *mesh_name = cs_probe_set_get_name(pset);
  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  post_mesh->sel_func[4] = NULL;
  post_mesh->sel_input[4] = pset;

  post_mesh->add_groups = false;

  if (auto_variable) {
    if (is_profile) {
      if (on_boundary)
        post_mesh->cat_id = CS_POST_MESH_BOUNDARY;
      else
        post_mesh->cat_id = CS_POST_MESH_VOLUME;
    }
    else
      post_mesh->cat_id = CS_POST_MESH_PROBES;
  }

  /* Try to assign probe location mesh */

  const char _select_all[] = "all[]";
  const char *sel_criteria = cs_probe_set_get_location_criteria(pset);
  if (sel_criteria == NULL)
    sel_criteria = _select_all;

  /* Check for existing meshes with the same selection criteria */

  int  ent_flag_id = (on_boundary) ? 2 : 0;
  int  match_partial[2] = {-1, -1};

  for (int i = 0; i < _cs_post_n_meshes; i++) {

    cs_post_mesh_t *post_mesh_cmp = _cs_post_meshes + i;
    if (post_mesh_cmp->criteria[ent_flag_id] != NULL) {
      if (strcmp(sel_criteria, post_mesh_cmp->criteria[ent_flag_id]) == 0) {
        if (time_varying == false)
          post_mesh->locate_ref = i;
        break;
      }
      else {
        if (post_mesh_cmp->n_writers == 0)
          match_partial[1] = i;
        else {
          for (int j = 0; j < n_writers && match_partial[0] == -1; j++) {
            for (int k = 0; k < post_mesh_cmp->n_writers; k++)
              if (writer_ids[j] == post_mesh_cmp->writer_id[k])
                match_partial[0] = i;
          }
        }
      }
    }
  }

  if (post_mesh->locate_ref < 0) {
    if (match_partial[0] >= 0)
      post_mesh->locate_ref = match_partial[0];
    else if (match_partial[1] >= 0)
      post_mesh->locate_ref = match_partial[1];
  }

  /* Add (define) location mesh if none found */

  if (post_mesh->locate_ref == -1) {
    int new_id = cs_post_get_free_mesh_id();
    if (on_boundary)
      cs_post_define_surface_mesh(new_id,
                                  "probe_set_location_mesh",
                                  NULL,
                                  sel_criteria,
                                  false,
                                  false,
                                  0,
                                  NULL);
    else
      cs_post_define_volume_mesh(new_id,
                                 "probe_set_location_mesh",
                                 sel_criteria,
                                 false,
                                 false,
                                 0,
                                 NULL);

    post_mesh->locate_ref = _cs_post_mesh_id(new_id);
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update "active" or "inactive" flag of writers based on the time step.
 *
 * Writers are activated if their output frequency is a divisor of the
 * current time step, or if their optional time step and value output lists
 * contain matches for the current time step.
 *----------------------------------------------------------------------------*/

void
cs_f_post_activate_by_time_step(void)
{
  cs_post_activate_by_time_step(cs_glob_time_step);
}

/*----------------------------------------------------------------------------
 * Output a floating point variable defined at cells or faces of a
 * post-processing mesh using associated writers.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   var_name    <-- name of variable to output
 *   var_dim     <-- 1 for scalar, 3 for vector
 *   interlace   <-- if a vector, true for interlaced values, false otherwise
 *   use_parent  <-- true if values are defined on "parent" mesh,
 *                   false if values are defined on post-processing mesh
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- current physical time
 *   cel_vals    <-- cell values
 *   i_face_vals <-- interior face values
 *   b_face_vals <-- boundary face values
 *----------------------------------------------------------------------------*/

void
cs_f_post_write_var(int               mesh_id,
                    const char       *var_name,
                    int               var_dim,
                    bool              interlace,
                    bool              use_parent,
                    int               nt_cur_abs,
                    double            t_cur_abs,
                    const cs_real_t  *cel_vals,
                    const cs_real_t  *i_face_vals,
                    const cs_real_t  *b_face_vals)
{
  CS_UNUSED(t_cur_abs);

  cs_post_type_t var_type
    = (sizeof(cs_real_t) == 8) ? CS_POST_TYPE_double : CS_POST_TYPE_float;

  const cs_time_step_t  *ts = cs_glob_time_step;

  if (nt_cur_abs < 0) /* Allow forcing of time-independent output */
    ts = NULL;

  cs_post_write_var(mesh_id,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    var_name,
                    var_dim,
                    interlace,
                    use_parent,
                    var_type,
                    cel_vals,
                    i_face_vals,
                    b_face_vals,
                    ts);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transform time-independent values into time-dependent
 *        values for transient meshes.
 *
 * \param[in]       writer  pointer to associated writer
 * \param[in, out]  nt_cur  associated time step (-1 initially for
 *                          time-independent values)
 * \param[in, out]  t_cur   associated time value
 */
/*----------------------------------------------------------------------------*/

static inline void
_check_non_transient(const cs_post_writer_t  *writer,
                     int                     *nt_cur,
                     double                  *t_cur)
{
  assert(writer->active > 0);

  fvm_writer_time_dep_t time_dep = fvm_writer_get_time_dep(writer->writer);

  if (time_dep == FVM_WRITER_TRANSIENT_CONNECT) {
    *nt_cur = writer->n_last;
    *t_cur = writer->t_last;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a writer; this objects manages a case's name, directory,
 *        and format, as well as associated mesh's time dependency, and the
 *        default output frequency for associated variables.
 *
 * This function must be called before the time loop. If a writer with a
 * given id is defined multiple times, the last definition supercedes the
 * previous ones.
 *
 * Current reserved ids are the following: CS_POST_WRITER_DEFAULT
 * for main/default output, CS_POST_WRITER_ERRORS for error visualization,
 * CS_POST_WRITER_PROBES for main probes, CS_POST_WRITER_PARTICLES for
 * particles, CS_POST_WRITER_TRAJECTORIES for trajectories. Other negative
 * ids may be dynamically reserved by the code depending on options.
 * Positive ids identify user-defined writers.
 *
 * \warning depending on the chosen format, the \em case_name may be
 * shortened (maximum number of characters: 32 for \em MED, 19 for \em EnSight,
 * or modified automatically (white-space or forbidden characters will be
 * replaced by "_").
 *
 * The \c \b format_name argument is used to choose the output format, and the
 * following values are allowed (assuming the matching
 * support was built):
 *
 * - \c \b EnSight \c \b Gold (\c \b EnSight also accepted)
 * - \c \b MED
 * - \c \b CGNS
 * - \c \b CCM (only for the full volume and boundary meshes)
 * - \c \b Catalyst (in-situ visualization)
 * - \c \b MEDCoupling (in-memory structure, to be used from other code)
 * - \c \b plot (comma or whitespace separated 2d plot files)
 * - \c \b time_plot (comma or whitespace separated time plot files)
 *
 * The format name is case-sensitive, so \c \b ensight or \c \b cgns are also valid.
 *
 * The optional \c \b fmt_opts character string contains a list of options related
 * to the format, separated by spaces or commas; these options include:
 *
 * - \c \b binary for a binary format version (default)
 * - \c \b big_endian to force outputs to be in \c \b big-endian mode
 *         (for \c \b EnSight).
 * - \c \b text for a text format version (for \c \b EnSight).
 * - \c \b adf for ADF file type (for \c \b CGNS).
 * - \c \b hdf5 for HDF5 file type (for \c \b CGNS, normally the default if
 *         HDF5 support is available).
 * - \c \b discard_polygons to prevent from exporting faces with more than
 *         four edges (which may not be recognized by some post-processing
 *         tools); such faces will therefore not appear in the post-processing
 *         mesh.
 * - \c \b discard_polyhedra to prevent from exporting elements which are
 *         neither tetrahedra, prisms, pyramids nor hexahedra (which may not
 *         be recognized by some post-processing tools); such elements will
 *         therefore not appear in the post-processing mesh.
 * - \c \b divide_polygons to divide faces with more than four edges into
 *         triangles, so that any post-processing tool can recognize them.
 * - \c \b divide_polyhedra} to divide elements which are neither tetrahedra,
 *         prisms, pyramids nor hexahedra into simpler elements (tetrahedra and
 *         pyramids), so that any post-processing tool can recognize them.
 * - \c \b separate_meshes to multiple meshes and associated fields to
 *         separate outputs.
 *
 * Note that the white-spaces in the beginning or in the end of the
 * character strings given as arguments here are suppressed automatically.
 *
 * \param[in]  writer_id        id of writer to create. (< 0 reserved,
 *                              > 0 for user); eveb for reserved ids,
 *                              the matching writer's options
 *                              may be redifined by calls to this function
 * \param[in]  case_name        associated case name
 * \param[in]  dir_name         associated directory name
 * \param[in]  fmt_name         associated format name
 * \param[in]  fmt_opts         associated format options string
 * \param[in]  time_dep         \ref FVM_WRITER_FIXED_MESH if mesh definitions
 *                              are fixed, \ref FVM_WRITER_TRANSIENT_COORDS if
 *                              coordinates change,
 *                              \ref FVM_WRITER_TRANSIENT_CONNECT if
 *                              connectivity changes
 * \param[in]  output_at_start  force output at calculation start if true
 * \param[in]  output_at_end    force output at calculation end if true
 * \param[in]  frequency_n      default output frequency in time-steps, or < 0
 * \param[in]  frequency_t      default output frequency in seconds, or < 0
 *                              (has priority over frequency_n)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_writer(int                     writer_id,
                      const char             *case_name,
                      const char             *dir_name,
                      const char             *fmt_name,
                      const char             *fmt_opts,
                      fvm_writer_time_dep_t   time_dep,
                      bool                    output_at_start,
                      bool                    output_at_end,
                      int                     frequency_n,
                      double                  frequency_t)
{
  /* local variables */

  int    i;

  cs_post_writer_t  *w = NULL;
  cs_post_writer_def_t  *wd = NULL;

  /* Initialize timer statistics if necessary */

  if (_post_out_stat_id < 0)
    _post_out_stat_id =  cs_timer_stats_id_by_name("postprocessing_output");

  /* Check if the required writer already exists */

  if (writer_id == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested post-processing writer number\n"
                "must be < 0 (reserved) or > 0 (user).\n"));

  for (i = 0; i < _cs_post_n_writers; i++) {
    if ((_cs_post_writers + i)->id == writer_id) {
      w = _cs_post_writers + i;
      BFT_FREE(w->ot);
      wd = w->wd;
      assert(wd != NULL);
      BFT_FREE(wd->case_name);
      BFT_FREE(wd->dir_name);
      BFT_FREE(wd->fmt_opts);
      break;
    }
  }

  if (i == _cs_post_n_writers) { /* New definition */

    /* Resize global writers array */

    if (_cs_post_n_writers == _cs_post_n_writers_max) {
      if (_cs_post_n_writers_max == 0)
        _cs_post_n_writers_max = 4;
      else
        _cs_post_n_writers_max *= 2;
      BFT_REALLOC(_cs_post_writers,
                  _cs_post_n_writers_max,
                  cs_post_writer_t);
    }

    if (writer_id < _cs_post_min_writer_id)
      _cs_post_min_writer_id = writer_id;
    _cs_post_n_writers += 1;

    w = _cs_post_writers + i;
    BFT_MALLOC(w->wd, 1, cs_post_writer_def_t);
    wd = w->wd;

  }

  /* Assign writer definition to the structure */

  w->id = writer_id;
  w->output_start = false;
  w->output_start = output_at_start;
  w->output_end = output_at_end;
  w->frequency_n = frequency_n;
  w->frequency_t = frequency_t;
  w->active = 0;
  w->n_last = -2;
  w->t_last = cs_glob_time_step->t_prev;
  w->ot = NULL;

  wd->time_dep = time_dep;

  BFT_MALLOC(wd->case_name, strlen(case_name) + 1, char);
  strcpy(wd->case_name, case_name);

  BFT_MALLOC(wd->dir_name, strlen(dir_name) + 1, char);
  strcpy(wd->dir_name, dir_name);

  wd->fmt_id = fvm_writer_get_format_id(fmt_name);

  if (fmt_opts != NULL) {
    BFT_MALLOC(wd->fmt_opts, strlen(fmt_opts) + 1, char);
    strcpy(wd->fmt_opts, fmt_opts);
  }
  else {
    BFT_MALLOC(wd->fmt_opts, 1, char);
    wd->fmt_opts[0] = '\0';
  }

  w->writer = NULL;

  /* If writer is the default writer (id -1), update defaults */

  if (writer_id == -1) {
    _cs_post_default_format_id = wd->fmt_id;
    if (wd->fmt_opts != NULL) {
      BFT_REALLOC(_cs_post_default_format_options,
                  strlen(wd->fmt_opts)+ 1,
                  char);
      strcpy(_cs_post_default_format_options, wd->fmt_opts);
    }
    else
      BFT_FREE(_cs_post_default_format_options);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a volume post-processing mesh.
 *
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name       associated mesh name
 * \param[in]  cell_criteria   selection criteria for cells
 * \param[in]  add_groups      if true, add group information if present
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_volume_mesh(int          mesh_id,
                           const char  *mesh_name,
                           const char  *cell_criteria,
                           bool         add_groups,
                           bool         auto_variables,
                           int          n_writers,
                           const int    writer_ids[])
{
  /* Call common initialization */

  cs_post_mesh_t *post_mesh = NULL;

  post_mesh = _predefine_mesh(mesh_id, true, 0, n_writers, writer_ids);

  /* Define mesh based on current arguments */

  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  if (cell_criteria != NULL) {
    BFT_MALLOC(post_mesh->criteria[0], strlen(cell_criteria) + 1, char);
    strcpy(post_mesh->criteria[0], cell_criteria);
  }
  post_mesh->ent_flag[0] = 1;

  post_mesh->add_groups = (add_groups) ? true : false;
  if (auto_variables)
    post_mesh->cat_id = CS_POST_MESH_VOLUME;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a volume post-processing mesh using a selection function.
 *
 * The selection may be updated over time steps if both the time_varying
 * flag is set to true and the mesh is only associated with writers defined
 * with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * Note: if the cell_select_input pointer is non-NULL, it must point
 * to valid data when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * \param[in]  mesh_id            id of mesh to define
 *                                (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name          associated mesh name
 * \param[in]  cell_select_func   pointer to cells selection function
 * \param[in]  cell_select_input  pointer to optional input data for the cell
 *                                selection function, or NULL
 * \param[in]  time_varying       if true, try to redefine mesh at each
 *                                output time
 * \param[in]  add_groups         if true, add group information if present
 * \param[in]  auto_variables     if true, automatic output of main variables
 * \param[in]  n_writers          number of associated writers
 * \param[in]  writer_ids         ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_volume_mesh_by_func(int                    mesh_id,
                                   const char            *mesh_name,
                                   cs_post_elt_select_t  *cell_select_func,
                                   void                  *cell_select_input,
                                   bool                   time_varying,
                                   bool                   add_groups,
                                   bool                   auto_variables,
                                   int                    n_writers,
                                   const int              writer_ids[])
{
  /* Call common initialization */

  cs_post_mesh_t *post_mesh = NULL;

  post_mesh = _predefine_mesh(mesh_id, time_varying, 0, n_writers, writer_ids);

  /* Define mesh based on current arguments */

  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  post_mesh->sel_func[0] = cell_select_func;
  post_mesh->sel_input[0] = cell_select_input;
  post_mesh->ent_flag[0] = 1;

  post_mesh->add_groups = (add_groups) ? true : false;
  if (auto_variables)
    post_mesh->cat_id = CS_POST_MESH_VOLUME;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a surface post-processing mesh.
 *
 * \param[in]  mesh_id          id of mesh to define
 *                              (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name        associated mesh name
 * \param[in]  i_face_criteria  selection criteria for interior faces
 * \param[in]  b_face_criteria  selection criteria for boundary faces
 * \param[in]  add_groups       if true, add group information if present
 * \param[in]  auto_variables   if true, automatic output of main variables
 * \param[in]  n_writers        number of associated writers
 * \param[in]  writer_ids       ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_surface_mesh(int          mesh_id,
                            const char  *mesh_name,
                            const char  *i_face_criteria,
                            const char  *b_face_criteria,
                            bool         add_groups,
                            bool         auto_variables,
                            int          n_writers,
                            const int    writer_ids[])
{
  /* Call common initialization */

  cs_post_mesh_t *post_mesh = NULL;

  post_mesh = _predefine_mesh(mesh_id, true, 0, n_writers, writer_ids);

  /* Define mesh based on current arguments */

  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  if (i_face_criteria != NULL) {
    BFT_MALLOC(post_mesh->criteria[1], strlen(i_face_criteria) + 1, char);
    strcpy(post_mesh->criteria[1], i_face_criteria);
    post_mesh->ent_flag[1] = 1;
  }

  if (b_face_criteria != NULL) {
    BFT_MALLOC(post_mesh->criteria[2], strlen(b_face_criteria) + 1, char);
    strcpy(post_mesh->criteria[2], b_face_criteria);
    post_mesh->ent_flag[2] = 1;
  }

  post_mesh->add_groups = (add_groups != 0) ? true : false;
  if (auto_variables && post_mesh->ent_flag[1] == 0)
    post_mesh->cat_id = CS_POST_MESH_BOUNDARY;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a surface post-processing mesh using selection functions.
 *
 * The selection may be updated over time steps if both the time_varying
 * flag is set to true and the mesh is only associated with writers defined
 * with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * Note: if i_face_select_input or b_face_select_input pointer is non-NULL,
 * it must point to valid data when the selection function is called,
 * so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * \param[in]  mesh_id              id of mesh to define
 *                                  (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name            associated mesh name
 * \param[in]  i_face_select_func   pointer to interior faces selection function
 * \param[in]  b_face_select_func   pointer to boundary faces selection function
 * \param[in]  i_face_select_input  pointer to optional input data for the
 *                                  interior faces selection function, or NULL
 * \param[in]  b_face_select_input  pointer to optional input data for the
 *                                  boundary faces selection function, or NULL
 * \param[in]  time_varying         if true, try to redefine mesh at each
 *                                  output time
 * \param[in]  add_groups           if true, add group information if present
 * \param[in]  auto_variables       if true, automatic output of main variables
 * \param[in]  n_writers            number of associated writers
 * \param[in]  writer_ids          ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_surface_mesh_by_func(int                    mesh_id,
                                    const char            *mesh_name,
                                    cs_post_elt_select_t  *i_face_select_func,
                                    cs_post_elt_select_t  *b_face_select_func,
                                    void                  *i_face_select_input,
                                    void                  *b_face_select_input,
                                    bool                   time_varying,
                                    bool                   add_groups,
                                    bool                   auto_variables,
                                    int                    n_writers,
                                    const int              writer_ids[])
{
  /* Call common initialization */

  cs_post_mesh_t *post_mesh = NULL;

  post_mesh = _predefine_mesh(mesh_id, time_varying, 0, n_writers, writer_ids);

  /* Define mesh based on current arguments */

  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  post_mesh->sel_func[1] = i_face_select_func;
  post_mesh->sel_func[2] = b_face_select_func;

  post_mesh->sel_input[1] = i_face_select_input;
  post_mesh->sel_input[2] = b_face_select_input;

  post_mesh->add_groups = (add_groups != 0) ? true : false;

  if (post_mesh->sel_func[1] != NULL)
    post_mesh->ent_flag[1] = 1;
  if (post_mesh->sel_func[2] != NULL)
    post_mesh->ent_flag[2] = 1;

  if (auto_variables)
    post_mesh->cat_id = CS_POST_MESH_BOUNDARY;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particles post-processing mesh.
 *
 * Such a mesh is always time-varying, and will only be output by writers
 * defined with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * If the trajectory_mode argument is set to true, this logic is reversed,
 * and output will only occur for writers defined with the
 * FVM_WRITER_FIXED_MESH option. In this case, a submesh consisting of
 * trajectory segments for the current time step will be added to
 * the output at each output time step.
 *
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name       associated mesh name
 * \param[in]  cell_criteria   selection criteria for cells containing
 *                             particles, or NULL.
 * \param[in]  density         fraction of the particles in the selected area
 *                             which should be output (0 < density <= 1)
 * \param[in]  trajectory      if true, activate trajectory mode
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_particles_mesh(int          mesh_id,
                              const char  *mesh_name,
                              const char  *cell_criteria,
                              double       density,
                              bool         trajectory,
                              bool         auto_variables,
                              int          n_writers,
                              const int    writer_ids[])
{
  /* Call common initialization */

  int flag = (trajectory) ? 2 : 1;
  cs_post_mesh_t *post_mesh = NULL;

  post_mesh = _predefine_mesh(mesh_id, true, flag, n_writers, writer_ids);

  /* Define mesh based on current arguments */

  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  if (cell_criteria != NULL) {
    BFT_MALLOC(post_mesh->criteria[3], strlen(cell_criteria) + 1, char);
    strcpy(post_mesh->criteria[3], cell_criteria);
  }

  post_mesh->add_groups = false;

  post_mesh->density = CS_MIN(density, 1.);
  post_mesh->density = CS_MAX(post_mesh->density, 0.);

  if (auto_variables)
    post_mesh->cat_id = CS_POST_MESH_VOLUME;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particles post-processing mesh using a selection function.
 *
 * The selection may be updated over time steps.
 *
 * Such a mesh is always time-varying, and will only be output by writers
 * defined with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * If the trajectory_mode argument is set to true, this logic is reversed,
 * and output will only occur for writers defined with the
 * FVM_WRITER_FIXED_MESH option. In this case, a submesh consisting of
 * trajectory segments for the current time step will be added to
 * the output at each output time step.
 *
 * Note: if the p_select_input pointer is non-NULL, it must point
 * to valid data when the selection function is called, so
 * that value or structure should not be temporary (i.e. local);
 *
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name       associated mesh name
 * \param[in]  p_select_func   pointer to particles selection function
 * \param[in]  p_select_input  pointer to optional input data for the particles
 *                             selection function, or NULL
 * \param[in]  trajectory      if true, activate trajectory mode
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_particles_mesh_by_func(int                    mesh_id,
                                      const char            *mesh_name,
                                      cs_post_elt_select_t  *p_select_func,
                                      void                  *p_select_input,
                                      bool                   trajectory,
                                      bool                   auto_variables,
                                      int                    n_writers,
                                      const int              writer_ids[])
{
  /* Call common initialization */

  int flag = (trajectory) ? 2 : 1;
  cs_post_mesh_t *post_mesh = NULL;

  post_mesh = _predefine_mesh(mesh_id, true, flag, n_writers, writer_ids);

  /* Define mesh based on current arguments */

  BFT_MALLOC(post_mesh->name, strlen(mesh_name) + 1, char);
  strcpy(post_mesh->name, mesh_name);

  post_mesh->sel_func[3] = p_select_func;
  post_mesh->sel_input[3] = p_select_input;
  post_mesh->ent_flag[3] = 1;

  post_mesh->add_groups = false;

  post_mesh->density = 1.;

  if (auto_variables)
    post_mesh->cat_id = CS_POST_MESH_PARTICLES;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a post-processing mesh associated with an existing exportable
 * mesh representation.
 *
 * If the exportable mesh is not intended to be used elsewhere, one can choose
 * to transfer its property to the post-processing mesh, which will then
 * manage its lifecycle based on its own requirements.
 *
 * If the exportable mesh must still be shared, one must be careful to
 * maintain consistency between this mesh and the post-processing output.
 *
 * The mesh in exportable dimension may be of a lower dimension than
 * its parent mesh, if it has been projected. In this case, a
 * dim_shift value of 1 indicates that parent cells are mapped to
 * exportable faces, and faces to edges, while a dim_shift value of 2
 * would indicate that parent cells are mapped to edges.
 * This is important when variables values are exported.
 *
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  exp_mesh        mesh in exportable representation
 *                             (i.e. fvm_nodal_t)
 * \param[in]  dim_shift       nonzero if exp_mesh has been projected
 * \param[in]  transfer        if true, ownership of exp_mesh is transferred
 *                             to the post-processing mesh
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_existing_mesh(int           mesh_id,
                             fvm_nodal_t  *exp_mesh,
                             int           dim_shift,
                             bool          transfer,
                             bool          auto_variables,
                             int           n_writers,
                             const int     writer_ids[])
{
  /* local variables */

  int        i;
  int        glob_flag[3];
  cs_lnum_t  b_f_num_shift, ind_fac;

  int    loc_flag[3] = {1, 1, 1};  /* Flags 0 to 2 "inverted" compared
                                      to others so as to use a single
                                      call to
                                      MPI_Allreduce(..., MPI_MIN, ...) */

  int         dim_ent = 0;
  int         dim_ext_ent = 0;
  bool        maj_ent_flag = false;
  cs_lnum_t   n_elts = 0;

  cs_lnum_t       *num_ent_parent = NULL;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Initialization of base structure */

  post_mesh = _predefine_mesh(mesh_id, true, 0, n_writers, writer_ids);

  /* Assign mesh to structure */

  post_mesh->exp_mesh = exp_mesh;

  if (transfer == true)
    post_mesh->_exp_mesh = exp_mesh;

  /* Compute number of cells and/or faces */

  dim_ext_ent = fvm_nodal_get_max_entity_dim(exp_mesh);
  dim_ent = dim_ext_ent + dim_shift;
  n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ext_ent);

  if (dim_ent == 3 && n_elts > 0)
    loc_flag[0] = 0;

  else if (dim_ent == 2 && n_elts > 0) {

    BFT_MALLOC(num_ent_parent, n_elts, cs_lnum_t);

    fvm_nodal_get_parent_num(exp_mesh, dim_ext_ent, num_ent_parent);

    b_f_num_shift = cs_glob_mesh->n_b_faces;
    for (ind_fac = 0; ind_fac < n_elts; ind_fac++) {
      if (num_ent_parent[ind_fac] > b_f_num_shift)
        post_mesh->n_i_faces += 1;
      else
        post_mesh->n_b_faces += 1;
    }

    BFT_FREE(num_ent_parent);

    if (post_mesh->n_i_faces > 0)
      loc_flag[1] = 0;
    else if (post_mesh->n_b_faces > 0)
      loc_flag[2] = 0;

  }

  for (i = 0; i < 3; i++)
    glob_flag[i] = loc_flag[i];

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Allreduce (loc_flag, glob_flag, 3, MPI_INT, MPI_MIN,
                   cs_glob_mpi_comm);
#endif

  /* Global indicators of mesh entity type presence;
     updated only if the mesh is not totally empty (for time-depending
     meshes, empty at certain times, we want to know the last type
     of entity used) */

  for (i = 0; i < 3; i++) {
    if (glob_flag[i] == 0)
      maj_ent_flag = true;
  }

  if (maj_ent_flag == true) {
    for (i = 0; i < 3; i++) {
      if (glob_flag[i] == 0)         /* Inverted glob_flag 0 to 2 logic */
        post_mesh->ent_flag[i] = 1;  /* (c.f. remark above) */
      else
        post_mesh->ent_flag[i] = 0;
    }
  }

  if (auto_variables) {
    post_mesh->cat_id = CS_POST_MESH_VOLUME;
    _check_mesh_cat_id(post_mesh);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * \param[in]  mesh_id       id of edges mesh to create
 *                           (< 0 reserved, > 0 for user)
 * \param[in]  base_mesh_id  id of existing mesh (< 0 reserved, > 0 for user)
 * \param[in]  n_writers     number of associated writers
 * \param[in]  writer_ids    ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_edges_mesh(int        mesh_id,
                          int        base_mesh_id,
                          int        n_writers,
                          const int  writer_ids[])
{
  /* local variables */

  cs_post_mesh_t *post_mesh = NULL;

  cs_post_mesh_t *post_base
    = _cs_post_meshes + _cs_post_mesh_id(base_mesh_id);

  /* Add and initialize base structure */

  post_mesh = _predefine_mesh(mesh_id, true, 0, n_writers, writer_ids);

  BFT_MALLOC(post_mesh->name,
             strlen(post_base->name) + strlen(_(" edges")) + 1,
             char);
  strcpy(post_base->name, post_base->name);
  strcat(post_mesh->name, _(" edges"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a writer to a postprocessing mesh.
 *
 * This function must be called during the postprocessing output definition
 * stage, before any output actually occurs.
 *
 * If called with a non-existing mesh or writer id, or if the writer is
 * already associated, no setting is changed, and this function
 * returns silently.
 *
 * \param[in]  mesh_id      id of mesh to define
 *                          (< 0 reserved, > 0 for user)
 * \param[in]  writer_id    id of writer to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_attach_writer(int  mesh_id,
                           int  writer_id)
{
  int _mesh_id = _cs_post_mesh_id_try(mesh_id);
  int _writer_id = _cs_post_writer_id_try(writer_id);

  if (_mesh_id < 0 || _writer_id < 0)
    return;

  cs_post_mesh_t *post_mesh = _cs_post_meshes + _mesh_id;

  /* Check we have not output this mesh yet */

  if (post_mesh->nt_last > -2)
    bft_error(__FILE__, __LINE__, 0,
              _("Error associating writer %d with mesh %d:"
                "output has already been done for this mesh, "
                "so mesh-writer association is locked."),
              writer_id, mesh_id);

  /* Ignore if writer id already associated */

  for (int i = 0; i < post_mesh->n_writers; i++) {
    if (post_mesh->writer_id[i] == _writer_id)
      return;
  }

  BFT_REALLOC(post_mesh->writer_id, post_mesh->n_writers + 1, int);
  post_mesh->writer_id[post_mesh->n_writers] = _writer_id;
  post_mesh->n_writers += 1;

  _update_mesh_writer_associations(post_mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief De-associate a writer from a postprocessing mesh.
 *
 * This function must be called during the postprocessing output definition
 * stage, before any output actually occurs.
 *
 * If called with a non-existing mesh or writer id, or if the writer was not
 * previously associated, no setting is changed, and this function
 * returns silently.
 *
 * \param[in]  mesh_id      id of mesh to define
 *                          (< 0 reserved, > 0 for user)
 * \param[in]  writer_id    id of writer to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_detach_writer(int  mesh_id,
                           int  writer_id)
{
  int _mesh_id = _cs_post_mesh_id_try(mesh_id);
  int _writer_id = _cs_post_writer_id_try(writer_id);

  if (_mesh_id < 0 || _writer_id < 0)
    return;

  cs_post_mesh_t *post_mesh = _cs_post_meshes + _mesh_id;

  /* Check we have not output this mesh yet */

  if (post_mesh->nt_last > -2)
    bft_error(__FILE__, __LINE__, 0,
              _("Error unassociating writer %d from mesh %d:"
                "output has already been done for this mesh, "
                "so mesh-writer association is locked."),
              writer_id, mesh_id);

  /* Ignore if writer id already associated */

  int i, j;
  for (i = 0, j = 0; i < post_mesh->n_writers; i++) {
    if (post_mesh->writer_id[i] != _writer_id) {
      post_mesh->writer_id[j] = post_mesh->writer_id[i];
      j++;
    }
  }

  if (j < post_mesh->n_writers) {

    post_mesh->n_writers = j;
    BFT_REALLOC(post_mesh->writer_id, post_mesh->n_writers, int);

    _update_mesh_writer_associations(post_mesh);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing meshes entity presence flag.
 *
 * This flag is an array of 3 integers, indicating the presence of elements
 * of given types on at least one subdomain (i.e. rank):
 *   0: presence of cells
 *   1: presence of interior faces
 *   2: presence of boundary faces
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  pointer to entity presence flag
 */
/*----------------------------------------------------------------------------*/

const int *
cs_post_mesh_get_ent_flag(int  mesh_id)
{
  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  return mesh->ent_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of cells
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of cells of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_cells(int  mesh_id)
{
  cs_lnum_t retval = 0;

  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL)
    retval = fvm_nodal_get_n_entities(mesh->exp_mesh, 3);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of cells
 *
 * The array of cell ids must be of at least size
 * cs_post_mesh_get_n_cells(mesh_id).
 *
 * \param[in]   mesh_id   postprocessing mesh id
 * \param[out]  cell_ids  array of associated cell ids (0 to n-1 numbering,
 *                        relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_cell_ids(int         mesh_id,
                          cs_lnum_t  *cell_ids)
{
  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL) {
    cs_lnum_t i;
    cs_lnum_t n_cells = fvm_nodal_get_n_entities(mesh->exp_mesh, 3);
    fvm_nodal_get_parent_num(mesh->exp_mesh, 3, cell_ids);
    for (i = 0; i < n_cells; i++)
      cell_ids[i] -= 1;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of interior faces
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of cells of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_i_faces(int  mesh_id)
{
  cs_lnum_t retval = 0;

  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL)
    retval = mesh->n_i_faces;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of boundary faces.
 *
 * The array of boundary face ids must be of at least size
 * cs_post_mesh_get_n_b_faces(mesh_id).
 *
 * \param[in]   mesh_id     postprocessing mesh id
 * \param[out]  i_face_ids  array of associated interior faces ids
 *                          (0 to n-1 numbering, relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_i_face_ids(int        mesh_id,
                            cs_lnum_t  i_face_ids[])
{
  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL) {
    cs_lnum_t i;
    cs_lnum_t n_faces = fvm_nodal_get_n_entities(mesh->exp_mesh, 2);
    const cs_lnum_t num_shift = cs_glob_mesh->n_b_faces + 1;
    if (mesh->n_b_faces == 0) {
      fvm_nodal_get_parent_num(mesh->exp_mesh, 3, i_face_ids);
      for (i = 0; i < n_faces; i++)
        i_face_ids[i] -= num_shift;
    }
    else {
      cs_lnum_t n_i_faces = 0;
      cs_lnum_t *tmp_ids = NULL;
      BFT_MALLOC(tmp_ids, n_faces, cs_lnum_t);
      fvm_nodal_get_parent_num(mesh->exp_mesh, 3, tmp_ids);
      for (i = 0; i < n_faces; i++) {
        if (tmp_ids[i] > cs_glob_mesh->n_b_faces)
          i_face_ids[n_i_faces++] = tmp_ids[i] - num_shift;
      }
      BFT_FREE(tmp_ids);
    }
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of boundary faces
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of cells of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_b_faces(int  mesh_id)
{
  cs_lnum_t retval = 0;

  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL)
    retval = mesh->n_b_faces;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of boundary faces.
 *
 * The array of boundary face ids must be of at least size
 * cs_post_mesh_get_n_b_faces(mesh_id).
 *
 * \param[in]   mesh_id     postprocessing mesh id
 * \param[out]  b_face_ids  array of associated boundary faces ids
 *                          (0 to n-1 numbering, relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_b_face_ids(int        mesh_id,
                            cs_lnum_t  b_face_ids[])
{
  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL) {
    cs_lnum_t i;
    cs_lnum_t n_faces = fvm_nodal_get_n_entities(mesh->exp_mesh, 2);
    if (mesh->n_i_faces == 0) {
      fvm_nodal_get_parent_num(mesh->exp_mesh, 3, b_face_ids);
      for (i = 0; i < n_faces; i++)
        b_face_ids[i] -= 1;
    }
    else {
      cs_lnum_t n_b_faces = 0;
      cs_lnum_t *tmp_ids = NULL;
      BFT_MALLOC(tmp_ids, n_faces, cs_lnum_t);
      fvm_nodal_get_parent_num(mesh->exp_mesh, 3, tmp_ids);
      for (i = 0; i < n_faces; i++) {
        if (tmp_ids[i] > cs_glob_mesh->n_b_faces)
          b_face_ids[n_b_faces++] = tmp_ids[i] - 1;
      }
      BFT_FREE(tmp_ids);
    }
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of vertices
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of vertices of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_vertices(int  mesh_id)
{
  cs_lnum_t retval = 0;

  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL)
    retval = fvm_nodal_get_n_entities(mesh->exp_mesh, 0);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of vertices
 *
 * The array of vertex ids must be of at least size
 * cs_post_mesh_get_n_vertices(mesh_id).
 *
 * \param[in]   mesh_id     postprocessing mesh id
 * \param[out]  vertex_ids  array of associated vertex ids (0 to n-1 numbering,
 *                          relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_vertex_ids(int         mesh_id,
                            cs_lnum_t  *vertex_ids)
{
  const cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  if (mesh->exp_mesh != NULL) {
    cs_lnum_t i;
    cs_lnum_t n_vertices = fvm_nodal_get_n_entities(mesh->exp_mesh, 0);
    fvm_nodal_get_parent_num(mesh->exp_mesh, 0, vertex_ids);
    for (i = 0; i < n_vertices; i++)
      vertex_ids[i] -= 1;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s called before post-processing meshes are built."),
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set whether postprocessing mesh's parallel domain should be output.
 *
 * \param[in]  mesh_id      postprocessing mesh id
 * \param[in]  post_domain  true if parallel domain should be output,
 *                          false otherwise.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_set_post_domain(int   mesh_id,
                             bool  post_domain)
{
  cs_post_mesh_t  *mesh = _cs_post_meshes + _cs_post_mesh_id(mesh_id);

  mesh->post_domain = post_domain;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove a post-processing mesh.
 *
 * No further post-processing output will be allowed on this mesh,
 * so the associated structures may be freed.
 *
 * A post-processing mesh that has been associated with a time-varying
 * writer may not be removed.
 *
 * \param[in]  mesh_id  postprocessing mesh id
 */
/*----------------------------------------------------------------------------*/

void
cs_post_free_mesh(int  mesh_id)
{
  int i;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Search for requested mesh */

  int _mesh_id = _cs_post_mesh_id(mesh_id);

  /* Check if mesh was referenced for probe location */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;
    if (post_mesh->locate_ref == _mesh_id)
      bft_error(__FILE__, __LINE__, 0,
                _("Post-processing mesh number %d has been referenced\n"
                  "by probe set mesh %d, so it may not be freed.\n"),
                mesh_id, post_mesh->id);
  }

  /* Now set pointer to mesh and check for time dependency */

  post_mesh = _cs_post_meshes + _mesh_id;

  for (i = 0; i < post_mesh->n_writers; i++) {

    cs_post_writer_t *writer = _cs_post_writers + post_mesh->writer_id[i];

    fvm_writer_time_dep_t time_dep = fvm_writer_get_time_dep(writer->writer);

    if (post_mesh->nt_last > -2 && time_dep != FVM_WRITER_FIXED_MESH)
      bft_error(__FILE__, __LINE__, 0,
                _("Post-processing mesh number %d has been associated\n"
                  "to writer %d which allows time-varying meshes, so\n"
                  "it may not be freed.\n"),
                mesh_id, writer->id);
  }

  /* Remove mesh if allowed */

  _free_mesh(_mesh_id);

  /* Finally, update free mesh ids */

  int min_id = _MIN_RESERVED_MESH_ID;
  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->id < min_id)
      min_id = post_mesh->id;
  }
  _cs_post_min_mesh_id = min_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for the existence of a writer of the given id.
 *
 * \param[in]  writer_id  writer id to check
 *
 * \return  true if writer with this id exists, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_post_writer_exists(int  writer_id)
{
  /* local variables */

  int id;
  cs_post_writer_t  *writer = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_writers; id++) {
    writer = _cs_post_writers + id;
    if (writer->id == writer_id)
      return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for the existence of a post-processing mesh of the given id.
 *
 * \param[in]  mesh_id  mesh id to check
 *
 * \return  true if mesh with this id exists, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_post_mesh_exists(int  mesh_id)
{
  int id;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_meshes; id++) {
    post_mesh = _cs_post_meshes + id;
    if (post_mesh->id == mesh_id)
      return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the default writer format name
 *
 * \return  name of the default writer format
 */
/*----------------------------------------------------------------------------*/

const char *
cs_post_get_default_format(void)
{
  return (fvm_writer_format_name(_cs_post_default_format_id));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the default writer format options
 *
 * \return  default writer format options string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_post_get_default_format_options(void)
{
  return (_cs_post_default_format_options);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the next "reservable" (i.e. non-user) writer id available.
 *
 * \return  the smallest negative integer present, -1
 */
/*----------------------------------------------------------------------------*/

int
cs_post_get_free_writer_id(void)
{
  return (_cs_post_min_writer_id - 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the next "reservable" (i.e. non-user) mesh id available.
 *
 * \return  the smallest negative integer present, -1
 */
/*----------------------------------------------------------------------------*/

int
cs_post_get_free_mesh_id(void)
{
  return (_cs_post_min_mesh_id - 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update "active" or "inactive" flag of writers based on the time step.
 *
 * Writers are activated if their output frequency is a divisor of the
 * current time step, or if their optional time step and value output lists
 * contain matches for the current time step.
 *
 * \param[in]  ts  time step status structure
 */
/*----------------------------------------------------------------------------*/

void
cs_post_activate_by_time_step(const cs_time_step_t  *ts)
{
  int  i;
  cs_post_writer_t  *writer;

  for (i = 0; i < _cs_post_n_writers; i++) {

    writer = _cs_post_writers + i;

    if (writer->active < 0)
      continue;

    /* In case of previous calls for a given time step,
       a writer's status may not be changed */

    if (writer->n_last == ts->nt_cur) {
      writer->active = 1;
      continue;
    }

    /* Activation based on frequency */

    writer->active = 0;

    if (writer->frequency_t > 0) {
      double  delta_t = ts->t_cur - writer->t_last;
      if (delta_t >= writer->frequency_t*(1-1e-6))
        writer->active = 1;
      /* delta_t = ts->t_cur - ts->t_prev; */
      /* if (delta_t < writer->frequency_t*(1-1e-6)) */
      /*   writer->active = 0; */
    }
    else if (writer->frequency_n > 0) {
      if (   ts->nt_cur % (writer->frequency_n) == 0
          && ts->nt_cur > 0
          && ts->nt_cur != ts->nt_prev)
        writer->active = 1;
    }

    if (ts->nt_cur == ts->nt_prev && writer->output_start)
      writer->active = 1;

    if (ts->nt_cur == ts->nt_max && writer->output_end)
      writer->active = 1;

    /* Activation based on time step lists */

    _activate_if_listed(writer, ts);

    /* Do not activate transient writers for time-independent stages */

    if (ts->nt_cur < 0) {
      fvm_writer_time_dep_t  time_dep;
      if (writer->writer)
        time_dep = fvm_writer_get_time_dep(writer->writer);
      else
        time_dep = writer->wd->time_dep;
      if (time_dep != FVM_WRITER_FIXED_MESH)
        writer->active = 0;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 * \param[in]  activate   false to deactivate, true to activate
 */
/*----------------------------------------------------------------------------*/

void
cs_post_activate_writer(int   writer_id,
                        bool  activate)
{
  int i;
  cs_post_writer_t  *writer;

  if (writer_id != 0) {
    i = _cs_post_writer_id(writer_id);
    writer = _cs_post_writers + i;
    writer->active = (activate) ? 1 : 0;
  }
  else {
    for (i = 0; i < _cs_post_n_writers; i++) {
      writer = _cs_post_writers + i;
      writer->active = (activate) ? 1 : 0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Disable specific writer or all writers not currently active until
 *        \ref cs_post_enable_writer or \ref cs_post_activate_writer
 *        is called for those writers.
 *
 * For each call to this function for a given writer, the same number
 * of calls to \ref cs_post_enable_writer or a single call to
 * \ref cs_post_activate_writer is required to re-enable the writer.
 *
 * This is useful to disable output even of fixed meshes in preprocessing
 * stages.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_disable_writer(int   writer_id)
{
  int i;
  cs_post_writer_t  *writer;

  if (writer_id != 0) {
    i = _cs_post_writer_id(writer_id);
    writer = _cs_post_writers + i;
    if (writer->active < 1)
      writer->active -= 1;
  }
  else {
    for (i = 0; i < _cs_post_n_writers; i++) {
      writer = _cs_post_writers + i;
      if (writer->active < 1)
        writer->active -= 1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enable a specific writer or all writers currently disabled by
 *        previous calls to \ref cs_post_disable_writer.
 *
 * For each previous call to \ref cs_post_disable_writer for a given writer,
 * a call to this function (or a single call to \ref cs_post_activate_writer)
 * is required to re-enable the writer.
 *
 * This is useful to disable output even of fixed meshes in preprocessing
 * stages.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_enable_writer(int   writer_id)
{
  int i;
  cs_post_writer_t  *writer;

  if (writer_id != 0) {
    i = _cs_post_writer_id(writer_id);
    writer = _cs_post_writers + i;
    if (writer->active < 0)
      writer->active += 1;
  }
  else {
    for (i = 0; i < _cs_post_n_writers; i++) {
      writer = _cs_post_writers + i;
      if (writer->active < 0)
        writer->active += 1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the FVM writer associated to a writer_id.
 *
 * \param[in]  writer_id  associated writer id
 *
 * \return  a pointer to a fvm_writer_t structure
 */
/*----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(int  writer_id)
{
  int  id;
  cs_post_writer_t  *writer = NULL;

  id = _cs_post_writer_id(writer_id);
  writer = _cs_post_writers + id;

  if (writer->writer == NULL)
    _init_writer(writer);

  return writer->writer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return time dependency associated to a writer_id.
 *
 * \param[in]  writer_id  associated writer id
 *
 * \return  associated writer's time dependency
 */
/*----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
cs_post_get_writer_time_dep(int  writer_id)
{
  int  id;
  cs_post_writer_t  *writer = NULL;

  fvm_writer_time_dep_t   time_dep = FVM_WRITER_FIXED_MESH;

  id = _cs_post_writer_id(writer_id);
  writer = _cs_post_writers + id;

  if (writer->wd != NULL)
    time_dep = writer->wd->time_dep;
  else if (writer->writer != NULL)
    time_dep = fvm_writer_get_time_dep(writer->writer);

  return time_dep;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an activation time step for a specific writer or for all writers.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 * \param[in]  nt         time step value to add (or remove)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_writer_t_step(int  writer_id,
                          int  nt)
{
  int i;

  if (writer_id != 0) {
    i = _cs_post_writer_id(writer_id);
    _add_writer_ts(_cs_post_writers + i, nt);
  }
  else {
    for (i = 0; i < _cs_post_n_writers; i++)
      _add_writer_ts(_cs_post_writers + i, nt);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an activation time value for a specific writer or for all writers.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 * \param[in]  t          time value to add (or remove)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_writer_t_value(int     writer_id,
                           double  t)
{
  int i;

  if (writer_id != 0) {
    i = _cs_post_writer_id(writer_id);
    _add_writer_tv(_cs_post_writers + i, t);
  }
  else {
    for (i = 0; i < _cs_post_n_writers; i++)
      _add_writer_tv(_cs_post_writers + i, t);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output post-processing meshes using associated writers.
 *
 * If the time step structure argument passed is NULL, a time-independent
 * output will be assumed.
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_meshes(const cs_time_step_t  *ts)
{
  int  i;
  cs_post_mesh_t  *post_mesh;

  int t_top_id = cs_timer_stats_switch(_post_out_stat_id);

  /* First loop on meshes, for probes and profiles (which must not be
     "reduced" afer first output, as coordinates may be required for
     interpolation, and also and share volume or surface location meshes) */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->ent_flag[4] != 0)
      _cs_post_write_mesh(post_mesh, ts);
  }

  /* Main Loops on meshes and writers for output */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->ent_flag[4] != 0)
      continue;
    _cs_post_write_mesh(post_mesh, ts);
    /* reduce mesh definitions if not required anymore */
    if (   post_mesh->mod_flag_max == FVM_WRITER_FIXED_MESH
        && post_mesh->_exp_mesh != NULL)
      fvm_nodal_reduce(post_mesh->_exp_mesh, 0);
  }

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a variable defined at cells or faces of a post-processing mesh
 *        using associated writers.
 *
 * \param[in]  mesh_id      id of associated mesh
 * \param[in]  writer_id    id of specified associated writer,
 *                          or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  var_name     name of variable to output
 * \param[in]  var_dim      1 for scalar, 3 for vector, 6 for symmetric tensor,
 *                          9 for non-symmetric tensor
 * \param[in]  interlace    if a vector, true for interlaced values,
 *                          false otherwise
 * \param[in]  use_parent   true if values are defined on "parent" mesh,
 *                          false if values are defined on post-processing mesh
 * \param[in]  var_type     variable's data type
 * \param[in]  cel_vals     cell values
 * \param[in]  i_face_vals  interior face values
 * \param[in]  b_face_vals  boundary face values
 * \param[in]  ts           time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_var(int                    mesh_id,
                  int                    writer_id,
                  const char            *var_name,
                  int                    var_dim,
                  bool                   interlace,
                  bool                   use_parent,
                  cs_post_type_t         var_type,
                  const void            *cel_vals,
                  const void            *i_face_vals,
                  const void            *b_face_vals,
                  const cs_time_step_t  *ts)
{
  cs_lnum_t  i;
  int        _mesh_id;

  cs_interlace_t  _interlace;
  cs_datatype_t    datatype;

  size_t       dec_ptr = 0;
  int          n_parent_lists = 0;
  cs_lnum_t    parent_num_shift[2]  = {0, 0};
  cs_real_t   *var_tmp = NULL;
  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_writer_t    *writer = NULL;

  const void  *var_ptr[2*9] = {NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL};

  int nt_cur = (ts != NULL) ? ts->nt_cur : -1;
  double t_cur = (ts != NULL) ? ts->t_cur : 0.;

  /* Initializations */

  _mesh_id = _cs_post_mesh_id_try(mesh_id);

  if (_mesh_id < 0)
    return;

  post_mesh = _cs_post_meshes + _mesh_id;

  if (interlace)
    _interlace = CS_INTERLACE;
  else
    _interlace = CS_NO_INTERLACE;

  datatype =  _cs_post_cnv_datatype(var_type);

  /* Assign appropriate array to FVM for output */

  /* Case of cells */
  /*---------------*/

  if (post_mesh->ent_flag[CS_POST_LOCATION_CELL] == 1) {

    if (use_parent) {
      n_parent_lists = 1;
      parent_num_shift[0] = 0;
    }
    else
      n_parent_lists = 0;

    var_ptr[0] = cel_vals;
    if (interlace == false) {
      if (use_parent)
        dec_ptr = cs_glob_mesh->n_cells_with_ghosts;
      else
        dec_ptr = fvm_nodal_get_n_entities(post_mesh->exp_mesh, 3);
      dec_ptr *= cs_datatype_size[datatype];
      for (i = 1; i < var_dim; i++)
        var_ptr[i] = ((const char *)cel_vals) + i*dec_ptr;
    }
  }

  /* Case of faces */
  /*---------------*/

  else if (   post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1
           || post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] == 1) {

    /* In case of indirection, all that is necessary is to set pointers */

    if (use_parent) {

      n_parent_lists = 2;
      parent_num_shift[0] = 0;
      parent_num_shift[1] = cs_glob_mesh->n_b_faces;

      if (post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] == 1) {
        if (interlace == false) {
          dec_ptr = cs_glob_mesh->n_b_faces * cs_datatype_size[datatype];
          for (i = 0; i < var_dim; i++)
            var_ptr[i] = ((const char *)b_face_vals) + i*dec_ptr;
        }
        else
          var_ptr[0] = b_face_vals;
      }

      if (post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1) {
        if (interlace == false) {
          dec_ptr = cs_glob_mesh->n_i_faces * cs_datatype_size[datatype];
          for (i = 0; i < var_dim; i++)
            var_ptr[var_dim + i] = ((const char *)i_face_vals) + i*dec_ptr;
        }
        else
          var_ptr[1] = i_face_vals;
      }

    }

    /* With no indirection, we must switch to a variable defined on two
       lists of faces to a variable defined on one list */

    else {

      n_parent_lists = 0;

      if (post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] == 1) {

        /* Case where a variable is defined both on boundary and
           interior faces: we must switch to a single list, as
           indirection is not used */

        if (post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1) {

          BFT_MALLOC(var_tmp,
                     (   post_mesh->n_i_faces
                      +  post_mesh->n_b_faces) * var_dim,
                     cs_real_t);

          _cs_post_assmb_var_faces(post_mesh->exp_mesh,
                                   post_mesh->n_i_faces,
                                   post_mesh->n_b_faces,
                                   var_dim,
                                   _interlace,
                                   i_face_vals,
                                   b_face_vals,
                                   var_tmp);

          _interlace = CS_NO_INTERLACE;

          dec_ptr = cs_datatype_size[datatype] * (  post_mesh->n_i_faces
                                                  + post_mesh->n_b_faces);

          for (i = 0; i < var_dim; i++)
            var_ptr[i] = ((char *)var_tmp) + i*dec_ptr;

        }

        /* Case where we only have boundary faces */

        else {

          if (interlace == false) {
            dec_ptr = cs_datatype_size[datatype] * post_mesh->n_b_faces;
            for (i = 0; i < var_dim; i++)
              var_ptr[i] = ((const char *)b_face_vals) + i*dec_ptr;
          }
          else
            var_ptr[0] = b_face_vals;
        }

      }

      /* Case where we only have interior faces */

      else if (post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1) {

        if (interlace == false) {
          dec_ptr = cs_datatype_size[datatype] * post_mesh->n_i_faces;
          for (i = 0; i < var_dim; i++)
            var_ptr[i] = ((const char *)i_face_vals) + i*dec_ptr;
        }
        else
          var_ptr[0] = i_face_vals;
      }

    }

  }

  /* Effective output: loop on writers */
  /*-----------------------------------*/

  for (i = 0; i < post_mesh->n_writers; i++) {

    writer = _cs_post_writers + post_mesh->writer_id[i];

    if (writer->id != writer_id && writer_id != CS_POST_WRITER_ALL_ASSOCIATED)
      continue;

    if (writer->active == 1) {

      _check_non_transient(writer, &nt_cur, &t_cur);

      fvm_writer_export_field(writer->writer,
                              post_mesh->exp_mesh,
                              var_name,
                              FVM_WRITER_PER_ELEMENT,
                              var_dim,
                              _interlace,
                              n_parent_lists,
                              parent_num_shift,
                              datatype,
                              nt_cur,
                              t_cur,
                              (const void * *)var_ptr);

      if (nt_cur >= 0) {
        writer->n_last = nt_cur;
        writer->t_last = t_cur;
      }

    }

  }

  /* Free memory (if both interior and boundary faces present) */

  if (var_tmp != NULL)
    BFT_FREE(var_tmp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a variable defined at vertices of a post-processing mesh using
 *        associated writers.
 *
 * \param[in]  mesh_id     id of associated mesh
 * \param[in]  writer_id   id of specified associated writer,
 *                         or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  var_name    name of variable to output
 * \param[in]  var_dim     1 for scalar, 3 for vector, 6 for symmetric tensor,
 *                         9 for non-symmetric tensor
 * \param[in]  interlace   if a vector, true for interlaced values,
 *                         false otherwise
 * \param[in]  use_parent  true if values are defined on "parent" mesh,
 *                         false if values are defined on post-processing mesh
 * \param[in]  var_type    variable's data type
 * \param[in]  vtx_vals    vertex values
 * \param[in]  ts          time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_vertex_var(int                    mesh_id,
                         int                    writer_id,
                         const char            *var_name,
                         int                    var_dim,
                         bool                   interlace,
                         bool                   use_parent,
                         cs_post_type_t         var_type,
                         const void            *vtx_vals,
                         const cs_time_step_t  *ts)
{
  cs_lnum_t  i;
  int        _mesh_id;

  cs_post_mesh_t  *post_mesh;
  cs_post_writer_t  *writer;
  cs_interlace_t  _interlace;
  cs_datatype_t  datatype;

  size_t       dec_ptr = 0;
  int          n_parent_lists = 0;
  cs_lnum_t    parent_num_shift[1]  = {0};

  const void  *var_ptr[9] = {NULL, NULL, NULL,
                             NULL, NULL, NULL,
                             NULL, NULL, NULL};

  int nt_cur = (ts != NULL) ? ts->nt_cur : -1;
  double t_cur = (ts != NULL) ? ts->t_cur : 0.;

  /* Initializations */

  _mesh_id = _cs_post_mesh_id_try(mesh_id);

  if (_mesh_id < 0)
    return;

  post_mesh = _cs_post_meshes + _mesh_id;

  if (interlace)
    _interlace = CS_INTERLACE;
  else
    _interlace = CS_NO_INTERLACE;

  assert(   sizeof(cs_real_t) == sizeof(double)
         || sizeof(cs_real_t) == sizeof(float));

  datatype =  _cs_post_cnv_datatype(var_type);

  /* Assign appropriate array to FVM for output */

  if (use_parent)
    n_parent_lists = 1;
  else
    n_parent_lists = 0;

  var_ptr[0] = vtx_vals;
  if (interlace == false) {
    if (use_parent)
      dec_ptr = cs_glob_mesh->n_vertices;
    else
      dec_ptr =   fvm_nodal_get_n_entities(post_mesh->exp_mesh, 0)
                * cs_datatype_size[datatype];
    for (i = 1; i < var_dim; i++)
      var_ptr[i] = ((const char *)vtx_vals) + i*dec_ptr;
  }

  /* Effective output: loop on writers */
  /*-----------------------------------*/

  for (i = 0; i < post_mesh->n_writers; i++) {

    writer = _cs_post_writers + post_mesh->writer_id[i];

    if (writer->id != writer_id && writer_id != CS_POST_WRITER_ALL_ASSOCIATED)
      continue;

    if (writer->active == 1) {

      _check_non_transient(writer, &nt_cur, &t_cur);

      fvm_writer_export_field(writer->writer,
                              post_mesh->exp_mesh,
                              var_name,
                              FVM_WRITER_PER_NODE,
                              var_dim,
                              _interlace,
                              n_parent_lists,
                              parent_num_shift,
                              datatype,
                              nt_cur,
                              t_cur,
                              (const void * *)var_ptr);

      if (nt_cur >= 0) {
        writer->n_last = nt_cur;
        writer->t_last = t_cur;
      }

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output an existing lagrangian particle attribute at particle
 *        positions or trajectory endpoints of a particle mesh using
 *        associated writers.
 *
 * \param[in]  mesh_id       id of associated mesh
 * \param[in]  writer_id     id of specified associated writer,
 *                           or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  attr_id       associated particle attribute id
 * \param[in]  var_name      name of variable to output
 * \param[in]  component_id  if -1 : extract the whole attribute
 *                           if >0 : id of the component to extract
 * \param[in]  ts            time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_particle_values(int                    mesh_id,
                              int                    writer_id,
                              int                    attr_id,
                              const char            *var_name,
                              int                    component_id,
                              const cs_time_step_t  *ts)
{
  int  _mesh_id, i, _length;
  int _stride_export_field = 1;
  cs_post_mesh_t  *post_mesh;
  cs_post_writer_t  *writer;

  cs_lagr_attribute_t  attr = attr_id;

  cs_lnum_t    n_particles = 0, n_pts = 0;
  cs_lnum_t    parent_num_shift[1]  = {0};
  cs_lnum_t   *particle_list = NULL;

  size_t  extents, size;
  ptrdiff_t  displ;
  cs_datatype_t datatype;
  int  stride;
  unsigned char *vals = NULL;

  int nt_cur = (ts != NULL) ? ts->nt_cur : -1;
  double t_cur = (ts != NULL) ? ts->t_cur : 0.;

  const void  *var_ptr[1] = {NULL};

  /* Initializations */

  _mesh_id = _cs_post_mesh_id_try(mesh_id);

  if (_mesh_id < 0)
    return;

  post_mesh = _cs_post_meshes + _mesh_id;

  if (post_mesh->ent_flag[3] == 0 || post_mesh->exp_mesh == NULL)
    return;

  n_particles = cs_lagr_get_n_particles();

  const cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();

  assert(p_set != NULL);

  /* Get attribute values info, returning if not present */

  cs_lagr_get_attr_info(p_set, 0, attr,
                        &extents, &size, &displ, &datatype, &stride);

  if (stride == 0)
    return;
  else {
    if (component_id == -1) {
      _length = size;
      _stride_export_field = stride;
    }
    else {
      _length = size/stride;
      _stride_export_field = 1;
     }
  }

  assert(ts->nt_cur > -1);

  /* Allocate work arrays */

  n_pts = fvm_nodal_get_n_entities(post_mesh->exp_mesh, 0);

  BFT_MALLOC(vals, n_pts*_length, unsigned char);

  var_ptr[0] = vals;

  if (n_pts != n_particles) {
    int parent_dim = (post_mesh->ent_flag[3] == 2) ? 1 : 0;
    BFT_MALLOC(particle_list, n_particles, cs_lnum_t);
    fvm_nodal_get_parent_num(post_mesh->exp_mesh, parent_dim, particle_list);
  }

  /* Particle values */

  if (post_mesh->ent_flag[3] == 1)
    cs_lagr_get_particle_values(p_set,
                                attr,
                                datatype,
                                stride,
                                component_id,
                                n_pts,
                                particle_list,
                                vals);

  else if (post_mesh->ent_flag[3] == 2) {
    nt_cur = -1; t_cur = 0.;
    cs_lagr_get_trajectory_values(p_set,
                                  attr,
                                  datatype,
                                  stride,
                                  component_id,
                                  n_pts/2,
                                  particle_list,
                                  vals);
  }

  BFT_FREE(particle_list);

  /* Effective output: loop on writers */
  /*-----------------------------------*/

  for (i = 0; i < post_mesh->n_writers; i++) {

    writer = _cs_post_writers + post_mesh->writer_id[i];

    if (writer->id != writer_id && writer_id != CS_POST_WRITER_ALL_ASSOCIATED)
      continue;

    if (writer->active == 1) {

      fvm_writer_export_field(writer->writer,
                              post_mesh->exp_mesh,
                              var_name,
                              FVM_WRITER_PER_NODE,
                              _stride_export_field,
                              CS_INTERLACE,
                              0, /* n_parent_lists, */
                              parent_num_shift,
                              datatype,
                              nt_cur,
                              t_cur,
                              (const void * *)var_ptr);

      if (nt_cur >= 0) {
        writer->n_last = nt_cur;
        writer->t_last = t_cur;
      }

    }

  }

  BFT_FREE(vals);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a variable defined at cells or faces of a post-processing mesh
 *        using associated writers.
 *
 * \param[in]  mesh_id              id of associated mesh
 * \param[in]  writer_id            id of specified associated writer,
 *                                  or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  var_name             name of variable to output
 * \param[in]  var_dim              1 for scalar, 3 for vector, 6 for symmetric
 *                                  tensor, 9 for non-symmetric tensor
 * \param[in]  var_type             variable's data type
 * \param[in]  parent_location_id   asociated values mesh location, or 0
 *                                  if values are passed directly
 * \param[in]  interpolate_func     pointer to interpolation function,
 *                                  or NULL for default
 * \param[in]  interpolate_input    pointer to optional interpolation input
 *                                  data, or NULL for default
 * \param[in]  vals                 variable's values
 * \param[in]  ts                   time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_probe_values(int                              mesh_id,
                           int                              writer_id,
                           const char                      *var_name,
                           int                              var_dim,
                           cs_post_type_t                   var_type,
                           int                              parent_location_id,
                           cs_interpolate_from_location_t  *interpolate_func,
                           void                            *interpolate_input,
                           const void                      *vals,
                           const cs_time_step_t            *ts)
{
  int nt_cur = (ts != NULL) ? ts->nt_cur : -1;
  double t_cur = (ts != NULL) ? ts->t_cur : 0.;

  bool on_boundary, on_curve;

  /* Initializations */

  int _mesh_id = _cs_post_mesh_id_try(mesh_id);

  if (_mesh_id < 0)
    return;

  cs_post_mesh_t  *post_mesh = _cs_post_meshes + _mesh_id;
  cs_probe_set_t  *pset = (cs_probe_set_t *)post_mesh->sel_input[4];

  cs_probe_set_get_post_info(pset,
                             NULL,
                             &on_boundary,
                             &on_curve,
                             NULL,
                             NULL,
                             NULL);

  cs_datatype_t datatype = _cs_post_cnv_datatype(var_type);

  const void  *var_ptr[1] = {vals};
  unsigned char *_vals = NULL;

  /* Extract or interpolate values */

  if (parent_location_id > 0) {

    const cs_lnum_t  n_points = fvm_nodal_get_n_entities(post_mesh->exp_mesh, 0);
    const cs_lnum_t *elt_ids = cs_probe_set_get_elt_ids(pset,
                                                        parent_location_id);

    cs_interpolate_from_location_t  *_interpolate_func = interpolate_func;

    cs_coord_t  *point_coords = NULL;

    if (_interpolate_func == NULL)
      _interpolate_func = cs_interpolate_from_location_p0;

    BFT_MALLOC(_vals,
               n_points*cs_datatype_size[datatype]*var_dim,
               unsigned char);

    if (_interpolate_func != cs_interpolate_from_location_p0) {
      BFT_MALLOC(point_coords, n_points*3, cs_coord_t);
      fvm_nodal_get_vertex_coords(post_mesh->exp_mesh,
                                  CS_INTERLACE,
                                  point_coords);
    }

    _interpolate_func(interpolate_input,
                      datatype,
                      var_dim,
                      n_points,
                      elt_ids,
                      (const cs_real_3_t *)point_coords,
                      vals,
                      _vals);
    var_ptr[0] = _vals;

    BFT_FREE(point_coords);
  }

  /* Effective output: loop on writers */
  /*-----------------------------------*/

  for (int i = 0; i < post_mesh->n_writers; i++) {

    cs_post_writer_t *writer = _cs_post_writers + post_mesh->writer_id[i];

    if (writer->id != writer_id && writer_id != CS_POST_WRITER_ALL_ASSOCIATED)
      continue;

    if (writer->active == 1) {

      cs_lnum_t  parent_num_shift[1] = {0};

      fvm_writer_export_field(writer->writer,
                              post_mesh->exp_mesh,
                              var_name,
                              FVM_WRITER_PER_NODE,
                              var_dim,
                              CS_INTERLACE,
                              0, /* n_parent_lists */
                              parent_num_shift,
                              datatype,
                              nt_cur,
                              t_cur,
                              (const void **)var_ptr);

      if (nt_cur >= 0) {
        writer->n_last = nt_cur;
        writer->t_last = t_cur;
      }

    }

  }

  BFT_FREE(_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update references to parent mesh of post-processing meshes in case of
 * computational mesh cell renumbering.
 *
 * This function may be called only once, after possible renumbering of cells,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * \param[in]  init_cell_num  initial cell numbering (new -> old)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_renum_cells(const cs_lnum_t  init_cell_num[])
{
  int        i;
  cs_lnum_t  icel;
  cs_lnum_t  n_elts;

  cs_lnum_t  *renum_ent_parent = NULL;

  bool  need_doing = false;

  cs_post_mesh_t   *post_mesh;
  const cs_mesh_t  *mesh = cs_glob_mesh;

  if (init_cell_num == NULL)
    return;

  /* Loop on meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    if (post_mesh->ent_flag[CS_POST_LOCATION_CELL] > 0)
      need_doing = true;
  }

  if (need_doing == true) {

    /* Prepare renumbering */

    n_elts = mesh->n_cells;

    BFT_MALLOC(renum_ent_parent, n_elts, cs_lnum_t);

    for (icel = 0; icel < mesh->n_cells; icel++)
      renum_ent_parent[init_cell_num[icel]] = icel + 1;

    /* Effective modification */

    for (i = 0; i < _cs_post_n_meshes; i++) {

      post_mesh = _cs_post_meshes + i;

      if (   post_mesh->_exp_mesh != NULL
          && post_mesh->ent_flag[CS_POST_LOCATION_CELL] > 0) {

        fvm_nodal_change_parent_num(post_mesh->_exp_mesh,
                                    renum_ent_parent,
                                    3);

      }

    }

    BFT_FREE(renum_ent_parent);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update references to parent mesh of post-processing meshes in case of
 * computational mesh interior and/or boundary faces renumbering.
 *
 * This function may be called only once, after possible renumbering of faces,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * \param[in]  init_i_face_num  initial interior numbering (new -> old)
 * \param[in]  init_b_face_num  initial boundary numbering (new -> old)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_renum_faces(const cs_lnum_t  init_i_face_num[],
                    const cs_lnum_t  init_b_face_num[])
{
  int       i;
  cs_lnum_t  ifac;
  cs_lnum_t  n_elts;

  cs_lnum_t  *renum_ent_parent = NULL;

  bool  need_doing = false;

  cs_post_mesh_t   *post_mesh;
  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Loop on meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    if (   post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] > 0
        || post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] > 0) {
      need_doing = true;
    }

  }

  if (need_doing == true) {

    /* Prepare renumbering */

    n_elts = mesh->n_i_faces + mesh->n_b_faces;

    BFT_MALLOC(renum_ent_parent, n_elts, cs_lnum_t);

    if (init_b_face_num == NULL) {
      for (ifac = 0; ifac < mesh->n_b_faces; ifac++)
        renum_ent_parent[ifac] = ifac + 1;
    }
    else {
      for (ifac = 0; ifac < mesh->n_b_faces; ifac++)
        renum_ent_parent[init_b_face_num[ifac]] = ifac + 1;
    }

    if (init_i_face_num == NULL) {
      for (ifac = 0, i = mesh->n_b_faces;
           ifac < mesh->n_i_faces;
           ifac++, i++)
        renum_ent_parent[mesh->n_b_faces + ifac]
          = mesh->n_b_faces + ifac + 1;
    }
    else {
      for (ifac = 0, i = mesh->n_b_faces;
           ifac < mesh->n_i_faces;
           ifac++, i++)
        renum_ent_parent[mesh->n_b_faces + init_i_face_num[ifac]]
          = mesh->n_b_faces + ifac + 1;
    }

    /* Effective modification */

    for (i = 0; i < _cs_post_n_meshes; i++) {

      post_mesh = _cs_post_meshes + i;

      if (   post_mesh->_exp_mesh != NULL
          && (   post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] > 0
              || post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] > 0)) {

        fvm_nodal_change_parent_num(post_mesh->_exp_mesh,
                                    renum_ent_parent,
                                    2);

      }

    }

    BFT_FREE(renum_ent_parent);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Configure the post-processing output so that a mesh displacement
 *        field may be output automatically for meshes based on the global
 *        volume mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_set_deformable(void)
{
  _cs_post_deformable = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Configure the post-processing output so that mesh connectivity
 * may be automatically updated.
 *
 * This is done for meshes defined using selection criteria or functions.
 * The behavior of Lagrangian meshes is unchanged.
 *
 * To be effective, this function should be called before defining
 * postprocessing meshes.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_set_changing_connectivity(void)
{
  _cs_post_mod_flag_min = FVM_WRITER_TRANSIENT_CONNECT;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize post-processing writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_init_writers(void)
{
  /* Ensure default is defined */

  if (!cs_post_writer_exists(CS_POST_WRITER_DEFAULT))
    cs_post_define_writer(CS_POST_WRITER_DEFAULT,   /* writer_id */
                          "results",                /* writer name */
                          _cs_post_dirname,
                          "EnSight Gold",           /* format name */
                          "",                       /* format options */
                          FVM_WRITER_FIXED_MESH,
                          false,                    /* output_at_start */
                          true,                     /* output at end */
                          -1,                       /* time step frequency */
                          -1.0);                    /* time value frequency */

  /* Additional writers for Lagrangian output */

  if (_lagrangian_needed(NULL)) {

    /* Particles */

    if (!cs_post_writer_exists(CS_POST_WRITER_PARTICLES))
      cs_post_define_writer(CS_POST_WRITER_PARTICLES,  /* writer_id */
                            "particles",            /* writer name */
                            _cs_post_dirname,
                            "EnSight Gold",         /* format name */
                            "",                     /* format options */
                            FVM_WRITER_TRANSIENT_CONNECT,
                            false,                  /* output_at_start */
                            true,                   /* output at end */
                            -1,                     /* time step frequency */
                            -1.0);                  /* time value frequency */

    if (!cs_post_writer_exists(CS_POST_WRITER_TRAJECTORIES))
      cs_post_define_writer(CS_POST_WRITER_TRAJECTORIES, /* writer_id */
                            "trajectories",         /* writer name */
                            _cs_post_dirname,
                            "EnSight Gold",         /* format name */
                            "",                     /* format options */
                            FVM_WRITER_FIXED_MESH,
                            false,                  /* output_at_start */
                            true,                   /* output at end */
                            1,                      /* time step frequency */
                            -1.0);                  /* time value frequency */

  }

  /* Additional writers for probe monitoring and profiles */

  if (!cs_post_writer_exists(CS_POST_WRITER_PROBES))
    cs_post_define_writer(CS_POST_WRITER_PROBES,    /* writer_id */
                          "",                       /* writer name */
                          "monitoring",
                          "time_plot",              /* format name */
                          "",                       /* format options */
                          FVM_WRITER_FIXED_MESH,
                          false,                    /* output_at_start */
                          false,                    /* output at end */
                          1,                        /* time step frequency */
                          -1.0);                    /* time value frequency */

  if (!cs_post_writer_exists(CS_POST_WRITER_PROFILES))
    cs_post_define_writer(CS_POST_WRITER_PROFILES,  /* writer_id */
                          "",                       /* writer name */
                          "profiles",
                          "plot",                   /* format name */
                          "",                       /* format options */
                          FVM_WRITER_FIXED_MESH,
                          false,                    /* output_at_start */
                          true,                     /* output at end */
                          -1,                       /* time step frequency */
                          -1.0);                    /* time value frequency */

  /* Print info on writers */

  _writer_info();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize main post-processing meshes
 *
 * The check_flag variable is a mask, used for additionnal post-processing:
 *
 *  - If (check_flag & 1), volume submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 *  - If (check_flag & 2), boundary submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 * It is recommended that post-processing meshes be defined before calling
 * this function, though specific "automatic" meshes (for example those
 * related to couplings) may be defined between this call and a time loop.
 *
 * \param[in]  check_mask  mask used for additional output
 */
/*----------------------------------------------------------------------------*/

void
cs_post_init_meshes(int check_mask)
{
  { /* Definition of default post-processing meshes if this has not been
       done yet */

    const int  writer_ids[] = {CS_POST_WRITER_DEFAULT};

    if (!cs_post_mesh_exists(CS_POST_MESH_VOLUME))
      cs_post_define_volume_mesh(CS_POST_MESH_VOLUME,
                                 _("Fluid domain"),
                                 "all[]",
                                 true,
                                 true,
                                 1,
                                 writer_ids);

    if (!cs_post_mesh_exists(CS_POST_MESH_BOUNDARY))
      cs_post_define_surface_mesh(CS_POST_MESH_BOUNDARY,
                                  _("Boundary"),
                                  NULL,
                                  "all[]",
                                  true,
                                  true,
                                  1,
                                  writer_ids);

  }

  /* Additional writers for Lagrangian output */

  if (_lagrangian_needed(NULL)) {
    if (!cs_post_mesh_exists(CS_POST_MESH_PARTICLES)) {
      const int writer_ids[] = {CS_POST_WRITER_PARTICLES};
      cs_post_define_particles_mesh(CS_POST_MESH_PARTICLES,
                                    _("Particles"),
                                    "all[]",
                                    1.0,    /* density */
                                    false,  /* trajectory */
                                    true,   /* auto_variables */
                                    1,
                                    writer_ids);
    }
  }

  /* Define probe meshes if needed */

  int n_probe_sets = cs_probe_get_n_sets();

  for (int pset_id = 0; pset_id < n_probe_sets; pset_id++) {

    bool  time_varying, is_profile, on_boundary, auto_variables;

    int  n_writers = 0;
    int  *writer_ids = NULL;
    cs_probe_set_t  *pset = cs_probe_set_get_by_id(pset_id);
    int  post_mesh_id = cs_post_get_free_mesh_id();

    cs_probe_set_get_post_info(pset,
                               &time_varying,
                               &on_boundary,
                               &is_profile,
                               &auto_variables,
                               &n_writers,
                               &writer_ids);

    if (is_profile) { /* User has to define an associated writer */

      _cs_post_define_probe_mesh(post_mesh_id,
                                 pset,
                                 time_varying,
                                 is_profile,
                                 on_boundary,
                                 auto_variables,
                                 n_writers,
                                 writer_ids);

    }
    else { /* Monitoring points and profiles*/

      /* Handle default writer assignment */

      if (n_writers < 0) {

        if (! is_profile) {
          const int  default_writer_ids[] = {CS_POST_WRITER_PROBES};
          cs_probe_set_associate_writers(pset, 1, default_writer_ids);
        }

        else
          cs_probe_set_associate_writers(pset, 0, NULL);

        cs_probe_set_get_post_info(pset, NULL, NULL, NULL, NULL,
                                   &n_writers, &writer_ids);

      }

      /* Now define associated mesh if active */

      if (n_writers > 0)
        _cs_post_define_probe_mesh(post_mesh_id,
                                   pset,
                                   time_varying,
                                   is_profile,
                                   on_boundary,
                                   auto_variables,
                                   n_writers,
                                   writer_ids);

    }

  } /* Loop on sets of probes */

  /* Remove meshes which are associated with no writer */

  _clear_unused_meshes();

  /* Add group parts if necessary (EnSight format) */

  if (check_mask & 1) {
    const char *fmt_name = fvm_writer_format_name(_cs_post_default_format_id);
    if (!strcmp(fmt_name, "EnSight Gold")) {
      _vol_submeshes_by_group(cs_glob_mesh,
                              fmt_name,
                              _cs_post_default_format_options);
      _boundary_submeshes_by_group(cs_glob_mesh,
                                   fmt_name,
                                   _cs_post_default_format_options);
    }
  }

#if 0
  /* Compute connectivity if not already done for delayed definitions */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    cs_post_mesh_t  *post_mesh = _cs_post_meshes + i;
    if (post_mesh->exp_mesh == NULL)
      _define_mesh(post_mesh, NULL);
  }
#endif

  /* If we must compute the vertices displacement field, we need
     to save the initial vertex coordinates */

  if (_cs_post_deformable && _cs_post_ini_vtx_coo == NULL) {
    cs_mesh_t *mesh = cs_glob_mesh;
    if (mesh->n_vertices > 0) {
      BFT_MALLOC(_cs_post_ini_vtx_coo,
                 mesh->n_vertices * 3,
                 cs_real_t);
      memcpy(_cs_post_ini_vtx_coo,
             mesh->vtx_coord,
             mesh->n_vertices * 3 * sizeof(cs_real_t));
    }
  }

  /* Initial output */

  cs_post_write_meshes(NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if post-processing is activated and then update post-processing
 *        of meshes if there is a need to update time-dependent meshes
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_time_step_begin(const cs_time_step_t  *ts)
{
  assert(ts != NULL); /* Sanity check */

  /* Activation or not of each writer according to the time step */

  cs_post_activate_by_time_step(ts);

  /* User-defined activation of writers for a fine-grained control */

  cs_user_postprocess_activate(ts->nt_max,
                               ts->nt_cur,
                               ts->t_cur);

  /* Possible modification of post-processing meshes */
  /*-------------------------------------------------*/

  _update_meshes(ts);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop on post-processing meshes to output variables.
 *
 * This handles all default fields output, as well as all
 * registered output functions and outputs defined in
 * \ref cs_user_postprocess_values
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_time_step_output(const cs_time_step_t  *ts)
{
  int  j;
  bool  active;

  /* Loop on writers to check if something must be done */
  /*----------------------------------------------------*/

  for (j = 0; j < _cs_post_n_writers; j++) {
    cs_post_writer_t  *writer = _cs_post_writers + j;
    if (writer->active == 1)
      break;
  }
  if (j == _cs_post_n_writers)
    return;

  int t_top_id = cs_timer_stats_switch(_post_out_stat_id);

    /* Output of variables by registered function instances */
  /*------------------------------------------------------*/

  for (int i = 0; i < _cs_post_n_output_tp; i++)
    _cs_post_f_output_tp[i](_cs_post_i_output_tp[i], ts);

  /* Output of variables associated with post-processing meshes */
  /*------------------------------------------------------------*/

  /* n_elts_max already initialized before and during the
     eventual modification of post-processing mesh definitions,
     and parent_ids allocated if n_elts_max > 0 */

  cs_lnum_t  *parent_ids = NULL;
  cs_lnum_t  n_elts_max = 0;

  /* Main loop on post-processing meshes */

  for (int i = 0; i < _cs_post_n_meshes; i++) {

    cs_post_mesh_t  *post_mesh = _cs_post_meshes + i;

    active = false;

    for (j = 0; j < post_mesh->n_writers; j++) {
      cs_post_writer_t  *writer = _cs_post_writers + post_mesh->writer_id[j];
      if (writer->active == 1)
        active = true;
    }

    /* If the mesh is active at this time step */
    /*-----------------------------------------*/

    if (active == true) {

      const fvm_nodal_t  *exp_mesh = post_mesh->exp_mesh;

      if (exp_mesh == NULL)
        continue;

      int  dim_ent = fvm_nodal_get_max_entity_dim(exp_mesh);
      cs_lnum_t  n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ent);

      if (n_elts > n_elts_max) {
        n_elts_max = n_elts;
        BFT_REALLOC(parent_ids, n_elts_max, cs_int_t);
      }

      /* Get corresponding element ids */

      fvm_nodal_get_parent_num(exp_mesh, dim_ent, parent_ids);

      for (cs_lnum_t k = 0; k < n_elts; k++)
        parent_ids[k] -= 1;

      /* We can output variables for this time step */
      /*--------------------------------------------*/

      cs_lnum_t  n_cells = 0, n_i_faces = 0, n_b_faces = 0;
      cs_lnum_t  *cell_ids = NULL, *i_face_ids = NULL, *b_face_ids = NULL;

      /* Here list sizes are adjusted, and we point to the array filled
         by fvm_nodal_get_parent_num() if possible. */

      if (dim_ent == 3) {
        n_cells = n_elts;
        cell_ids = parent_ids;
      }

      /* The numbers of "parent" interior faces known by FVM
         are shifted by the total number of boundary faces */

      else if (dim_ent == 2 && n_elts > 0) {

        cs_lnum_t  b_f_num_shift = cs_glob_mesh->n_b_faces;

        for (cs_lnum_t ind_fac = 0; ind_fac < n_elts; ind_fac++) {
          if (parent_ids[ind_fac] >= b_f_num_shift)
            n_i_faces++;
          else
            n_b_faces++;
        }

        /* boundary faces only: parent FVM face numbers unchanged */
        if (n_i_faces == 0) {
          b_face_ids = parent_ids;
        }

        /* interior faces only: parents FVM face numbers shifted */
        else if (n_b_faces == 0) {
          for (cs_lnum_t ind_fac = 0; ind_fac < n_elts; ind_fac++)
            parent_ids[ind_fac] -= b_f_num_shift;
          i_face_ids = parent_ids;
        }

        /* interior and boundary faces: numbers must be separated */

        else {

          BFT_MALLOC(i_face_ids, n_i_faces, cs_int_t);
          BFT_MALLOC(b_face_ids, n_b_faces, cs_int_t);

          n_i_faces = 0, n_b_faces = 0;

          for (cs_lnum_t ind_fac = 0; ind_fac < n_elts; ind_fac++) {
            if (parent_ids[ind_fac] >= b_f_num_shift)
              i_face_ids[n_i_faces++] = parent_ids[ind_fac] - b_f_num_shift;
            else
              b_face_ids[n_b_faces++] = parent_ids[ind_fac];
          }

        }

        /* In all cases, update the number of interior and boundary faces
           (useful in case of splitting of FVM mesh elements) for functions
           called by this one */

        post_mesh->n_i_faces = n_i_faces;
        post_mesh->n_b_faces = n_b_faces;

      }

      /* Output of zone information if necessary */

      _cs_post_write_transient_zone_info(post_mesh, ts);

      /* Standard post-processing */

      if (post_mesh->cat_id < 0)
        _cs_post_output_fields(post_mesh, ts);

      /* Output of variables by registered function instances */

      for (j = 0; j < _cs_post_n_output_mtp; j++)
        _cs_post_f_output_mtp[j](_cs_post_i_output_mtp[j],
                                 post_mesh->id,
                                 post_mesh->cat_id,
                                 post_mesh->ent_flag,
                                 n_cells,
                                 n_i_faces,
                                 n_b_faces,
                                 cell_ids,
                                 i_face_ids,
                                 b_face_ids,
                                 ts);

      /* User-defined output */

      cs_lnum_t  n_vertices = cs_post_mesh_get_n_vertices(post_mesh->id);

      if (post_mesh->sel_input[4] == NULL) {

        cs_lnum_t *vertex_ids;
        BFT_MALLOC(vertex_ids, n_vertices, cs_lnum_t);
        cs_post_mesh_get_vertex_ids(post_mesh->id, vertex_ids);

        cs_user_postprocess_values(post_mesh->name,
                                   post_mesh->id,
                                   post_mesh->cat_id,
                                   NULL,
                                   n_cells,
                                   n_i_faces,
                                   n_b_faces,
                                   n_vertices,
                                   cell_ids,
                                   i_face_ids,
                                   b_face_ids,
                                   vertex_ids,
                                   ts);

        BFT_FREE(vertex_ids);

        /* In case of mixed interior and boundary faces, free
           additional arrays */

        if (i_face_ids != NULL && b_face_ids != NULL) {
          BFT_FREE(i_face_ids);
          BFT_FREE(b_face_ids);
        }

      }

      else { /* if (post_mesh->sel_input[4] != NULL) */

        bool on_boundary = false;
        const cs_lnum_t *_cell_ids = NULL, *_b_face_ids = NULL;
        const cs_lnum_t *vertex_ids = NULL;
        cs_probe_set_t  *pset = (cs_probe_set_t *)post_mesh->sel_input[4];
        const char *mesh_name = cs_probe_set_get_name(pset);

        cs_probe_set_get_post_info(pset,
                                   NULL,
                                   &on_boundary,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL);

        if (on_boundary) {
          n_b_faces = n_vertices; /* n_probes */
          _b_face_ids = cs_probe_set_get_elt_ids(pset,
                                                 CS_MESH_LOCATION_BOUNDARY_FACES);
        }
        else {
          n_cells = n_vertices; /* n_probes */
          _cell_ids = cs_probe_set_get_elt_ids(pset,
                                               CS_MESH_LOCATION_CELLS);
        }
        vertex_ids = cs_probe_set_get_elt_ids(pset,
                                              CS_MESH_LOCATION_VERTICES);

        cs_user_postprocess_values(mesh_name,
                                   post_mesh->id,
                                   post_mesh->cat_id,
                                   pset,
                                   n_cells,
                                   0,
                                   n_b_faces,
                                   n_vertices,
                                   _cell_ids,
                                   NULL,
                                   _b_face_ids,
                                   vertex_ids,
                                   ts);

      }

    }

  }

  /* Free memory */

  BFT_FREE(parent_ids);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flush writers and free time-varying and Lagragian mesh if needed
 *        of meshes if there is a time-dependent mesh
 */
/*----------------------------------------------------------------------------*/

void
cs_post_time_step_end(void)
{
  int t_top_id = cs_timer_stats_switch(_post_out_stat_id);

  /* Flush writers if necessary */

  for (int i = 0; i < _cs_post_n_writers; i++) {
    cs_post_writer_t  *writer = _cs_post_writers + i;
    if (writer->active == 1) {
      if (writer->writer != NULL)
        fvm_writer_flush(writer->writer);
    }
  }

  /* Free time-varying and Lagrangian meshes unless they
     are mapped to an existing mesh */

  for (int i = 0; i < _cs_post_n_meshes; i++) {
    cs_post_mesh_t  *post_mesh = _cs_post_meshes + i;
    if (post_mesh->_exp_mesh != NULL) {
      if (   post_mesh->ent_flag[3]
          || post_mesh->mod_flag_min == FVM_WRITER_TRANSIENT_CONNECT) {
        post_mesh->exp_mesh = NULL;
        post_mesh->_exp_mesh = fvm_nodal_destroy(post_mesh->_exp_mesh);
      }
    }
  }

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop on post-processing meshes to output variables.
 *
 * This handles all default fields output, as well as all
 * registered output functions and outputs defined in
 * \ref cs_user_postprocess_values
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_vars(const cs_time_step_t  *ts)
{
  /* Output meshes if needed */

  _update_meshes(ts);

  /* Loop on post-processing meshes to output variables */

  cs_post_time_step_output(ts);

  /* Flush writers and free time-varying and Lagragian mesh if needed */

  cs_post_time_step_end();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all structures associated with post-processing
 */
/*----------------------------------------------------------------------------*/

void
cs_post_finalize(void)
{
  int i, j;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Timings */

  for (i = 0; i < _cs_post_n_writers; i++) {
    cs_timer_counter_t m_time, f_time, a_time;
    fvm_writer_t *writer = (_cs_post_writers + i)->writer;
    CS_TIMER_COUNTER_INIT(m_time);
    CS_TIMER_COUNTER_INIT(f_time);
    CS_TIMER_COUNTER_INIT(a_time);
    if (writer != NULL) {
      fvm_writer_get_times(writer,
                           &m_time, &f_time, &a_time);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("\n"
                      "Writing of \"%s\" (%s) summary:\n"
                      "\n"
                      "  CPU time for meshes:              %12.3f\n"
                      "  CPU time for variables:           %12.3f\n"
                      "  CPU time forcing output:          %12.3f\n"
                      "\n"
                      "  Elapsed time for meshes:          %12.3f\n"
                      "  Elapsed time for variables:       %12.3f\n"
                      "  Elapsed time forcing output:      %12.3f\n"),
                    fvm_writer_get_name(writer),
                    fvm_writer_get_format(writer),
                    m_time.cpu_nsec*1e-9,
                    f_time.cpu_nsec*1e-9,
                    a_time.cpu_nsec*1e-9,
                    m_time.wall_nsec*1e-9,
                    f_time.wall_nsec*1e-9,
                    a_time.wall_nsec*1e-9);
    }
  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  /* Initial coordinates (if mesh is deformable) */

  if (_cs_post_ini_vtx_coo != NULL)
    BFT_FREE(_cs_post_ini_vtx_coo);

  /* Exportable meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->_exp_mesh != NULL)
      fvm_nodal_destroy(post_mesh->_exp_mesh);
    BFT_FREE(post_mesh->name);
    for (j = 0; j < 4; j++)
      BFT_FREE(post_mesh->criteria[j]);
    BFT_FREE(post_mesh->writer_id);
  }

  BFT_FREE(_cs_post_meshes);

  _cs_post_min_mesh_id = _MIN_RESERVED_MESH_ID;
  _cs_post_n_meshes = 0;
  _cs_post_n_meshes_max = 0;

  /* Writers */

  for (i = 0; i < _cs_post_n_writers; i++) {
    cs_post_writer_t  *writer = _cs_post_writers + i;
    if (writer->ot != NULL)
      _free_writer_times(writer);
    if (writer->wd != NULL)
      _destroy_writer_def(writer);
    if (writer->writer != NULL)
      fvm_writer_finalize((_cs_post_writers + i)->writer);
  }

  BFT_FREE(_cs_post_writers);

  _cs_post_n_writers = 0;
  _cs_post_n_writers_max = 0;

  /* Registered processings if necessary */

  if (_cs_post_n_output_tp_max > 0) {
    BFT_FREE(_cs_post_f_output_tp);
    BFT_FREE(_cs_post_i_output_tp);
  }

  if (_cs_post_n_output_mtp_max > 0) {
    BFT_FREE(_cs_post_f_output_mtp);
    BFT_FREE(_cs_post_i_output_mtp);
  }

  /* Options */

  BFT_FREE(_cs_post_default_format_options);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Postprocess free (isolated) faces of the current global mesh
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_free_faces(void)
{
  cs_lnum_t i, j;
  cs_lnum_t n_f_faces = 0;
  cs_gnum_t n_no_group = 0;
  int max_null_family = 0;
  cs_lnum_t *f_face_list = NULL;

  fvm_writer_t *writer = NULL;
  fvm_nodal_t *exp_mesh = NULL;

  bool  generate_submeshes = false;
  cs_mesh_t *mesh = cs_glob_mesh;
  const char *fmt_name = fvm_writer_format_name(_cs_post_default_format_id);

  if (mesh->n_g_free_faces == 0)
    return;

  /* Create default writer */

  writer = fvm_writer_init("isolated_faces",
                           _cs_post_dirname,
                           fmt_name,
                           _cs_post_default_format_options,
                           FVM_WRITER_FIXED_MESH);

  /* Build list of faces to extract */

  BFT_MALLOC(f_face_list, mesh->n_b_faces, cs_lnum_t);

  for (i = 0; i < mesh->n_b_faces; i++) {
    if (mesh->b_face_cells[i] < 0)
      f_face_list[n_f_faces++] = i + 1;
  }

  /* Extract and output mesh of isolated faces */

  exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            "isolated faces",
                                            true,
                                            0,
                                            n_f_faces,
                                            NULL,
                                            f_face_list);

  if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_FACE_POLY) > 0)
    fvm_nodal_tesselate(exp_mesh, FVM_FACE_POLY, NULL);

  fvm_writer_set_mesh_time(writer, -1, 0);
  fvm_writer_export_nodal(writer, exp_mesh);

  exp_mesh = fvm_nodal_destroy(exp_mesh);

  /* Now check if we should generate additional meshes (EnSight Gold format) */

  if (!strcmp(fmt_name, "EnSight Gold") && mesh->n_families > 0) {

    generate_submeshes = true;

    /* Families should be sorted, so if a nonzero family is empty,
       it is family 1 */
    if (mesh->family_item[0] == 0)
      max_null_family = 1;
    if (mesh->n_families <= max_null_family)
      generate_submeshes = false;

    /* Check how many boundary faces belong to no group */

    if (mesh->b_face_family != NULL) {
      for (j = 0; j < n_f_faces; j++) {
        if (mesh->b_face_family[f_face_list[j] - 1] <= max_null_family)
          n_no_group += 1;
      }
    }
    else
      n_no_group = n_f_faces;

    cs_parall_counter(&n_no_group, 1);

    if (n_no_group == mesh->n_g_free_faces)
      generate_submeshes = false;
  }

  /* Generate submeshes if necessary */

  if (generate_submeshes) {

    cs_lnum_t n_b_faces;
    int *fam_flag = NULL;
    char *group_flag = NULL;
    cs_lnum_t *b_face_list = NULL;
    char part_name[81];

    /* Now detect which groups may be referenced */

    BFT_MALLOC(fam_flag, mesh->n_families + 1, int);
    memset(fam_flag, 0, (mesh->n_families + 1)*sizeof(int));

    if (mesh->b_face_family != NULL) {
      for (i = 0; i < n_f_faces; i++)
        fam_flag[mesh->b_face_family[f_face_list[i] - 1]] = 1;
    }

    group_flag = _build_group_flag(mesh, fam_flag);

    /* Now extract isolated faces by groups.
       Selector structures may not have been initialized yet,
       so we use a direct selection here. */

    BFT_REALLOC(fam_flag, mesh->n_families, int);

    BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);

    for (i = 0; i < mesh->n_groups; i++) {

      if (group_flag[i] != 0) {

        const char *g_name = mesh->group + mesh->group_idx[i];

        _set_fam_flags(mesh, i, fam_flag);

        n_b_faces = 0;
        if (mesh->b_face_family != NULL) {
          for (j = 0; j < n_f_faces; j++) {
            cs_lnum_t face_id = f_face_list[j] - 1;
            int fam_id = mesh->b_face_family[face_id];
            if (fam_id > 0 && fam_flag[fam_id - 1])
              b_face_list[n_b_faces++] = face_id + 1;
          }
        }

        strcpy(part_name, "isolated: ");
        strncat(part_name, g_name, 80 - strlen(part_name));

        exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                  part_name,
                                                  false,
                                                  0,
                                                  n_b_faces,
                                                  NULL,
                                                  b_face_list);

        if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_FACE_POLY) > 0)
          fvm_nodal_tesselate(exp_mesh, FVM_FACE_POLY, NULL);

        fvm_writer_set_mesh_time(writer, -1, 0);
        fvm_writer_export_nodal(writer, exp_mesh);

        exp_mesh = fvm_nodal_destroy(exp_mesh);
      }

    }

    /* Output boundary faces belonging to no group */

    if (n_no_group > 0) {

      if (mesh->b_face_family != NULL) {
        for (j = 0, n_b_faces = 0; j < n_f_faces; j++) {
          cs_lnum_t face_id = f_face_list[j] - 1;
          if (mesh->b_face_family[face_id] <= max_null_family)
            b_face_list[n_b_faces++] = face_id + 1;
        }
      }
      else {
        for (j = 0, n_b_faces = 0; j < n_f_faces; j++) {
          cs_lnum_t face_id = f_face_list[j] - 1;
          b_face_list[n_b_faces++] = face_id + 1;
        }
      }

      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                "isolated: no_group",
                                                false,
                                                0,
                                                n_b_faces,
                                                NULL,
                                                b_face_list);

      if (fvm_writer_needs_tesselation(writer, exp_mesh, FVM_FACE_POLY) > 0)
        fvm_nodal_tesselate(exp_mesh, FVM_FACE_POLY, NULL);

      fvm_writer_set_mesh_time(writer, -1, 0);
      fvm_writer_export_nodal(writer, exp_mesh);

      exp_mesh = fvm_nodal_destroy(exp_mesh);
    }

    BFT_FREE(b_face_list);

    BFT_FREE(fam_flag);
    BFT_FREE(group_flag);

  } /* End of submeshes generation */

  /* Free memory */

  writer = fvm_writer_finalize(writer);

  BFT_FREE(f_face_list);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * troubleshoot errors.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_init_error_writer(void)
{
  /* Default values */

  int writer_id = CS_POST_WRITER_ERRORS;
  if (cs_post_writer_exists(writer_id))
    return;

  /* Create default writer */

  int default_format_id = _cs_post_default_format_id;
  const char *default_format_options = _cs_post_default_format_options;
  const char null_str[] = "";

  /* Special case for Catalyst: if matching co-processing script is
     not available, revert to EnSight Gold format */

  if (default_format_id == fvm_writer_get_format_id("Catalyst")) {
    if (! cs_file_isreg("error.py")) {
      default_format_id = fvm_writer_get_format_id("EnSight Gold");
      default_format_options = null_str;
    }
  }

  cs_post_define_writer(writer_id,
                        "error",
                        _cs_post_dirname,
                        fvm_writer_format_name(default_format_id),
                        default_format_options,
                        FVM_WRITER_FIXED_MESH, /* No time dependency here */
                        false,
                        true,
                        -1,
                        -1.0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, and associate
 * and output global volume mesh.
 *
 * This is intended to help troubleshoot errors using fields based
 * on cells.
 *
 * \return  id of error output mesh (< 0), or 0 if all writers are deactivated
 */
/*----------------------------------------------------------------------------*/

int
cs_post_init_error_writer_cells(void)
{
  int mesh_id = 0;

  const cs_mesh_t *mesh = cs_glob_mesh;

  /* If post-processing is active, output info */
  /*-------------------------------------------*/

  if (mesh->i_face_vtx_idx != NULL || mesh->b_face_vtx_idx != NULL) {

    const int writer_id = CS_POST_WRITER_ERRORS;
    const char *mesh_name = N_("Calculation domain");
    cs_post_mesh_t *post_mesh = NULL;

    cs_post_init_error_writer();
    cs_post_activate_writer(writer_id, 1);

    mesh_id = cs_post_get_free_mesh_id();

    post_mesh = _predefine_mesh(mesh_id, false, 0, 1, &writer_id);

    /* Define mesh based on current arguments */

    BFT_MALLOC(post_mesh->name, strlen(_(mesh_name)) + 1, char);
    strcpy(post_mesh->name, _(mesh_name));

    post_mesh->add_groups = false;

    _define_export_mesh(post_mesh,
                        mesh->n_cells,
                        0,
                        0,
                        NULL,
                        NULL,
                        NULL);

    _cs_post_write_mesh(post_mesh, NULL);

  }

  return mesh_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Register a processing of time-dependent variables to the call to
 * cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * \param[in]       function  function to register
 * \param[in, out]  input     pointer to optional (untyped) value or structure
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_time_dep_output(cs_post_time_dep_output_t  *function,
                            void                       *input)
{
  /* Resize array of registered post-processings if necessary */

  if (_cs_post_n_output_tp >= _cs_post_n_output_tp_max) {
    if (_cs_post_n_output_tp_max == 0)
      _cs_post_n_output_tp_max = 8;
    else
      _cs_post_n_output_tp_max *= 2;
    BFT_REALLOC(_cs_post_f_output_tp,
                _cs_post_n_output_tp_max,
                cs_post_time_dep_output_t *);
    BFT_REALLOC(_cs_post_i_output_tp, _cs_post_n_output_tp_max, void *);
  }

  /* Add a post-processing */

  _cs_post_f_output_tp[_cs_post_n_output_tp] = function;
  _cs_post_i_output_tp[_cs_post_n_output_tp] = input;

  _cs_post_n_output_tp += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Register a processing of time-dependent variables than can be output
 * on different meshes to the call to cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * \param[in]       function  function to register
 * \param[in, out]  input     pointer to optional (untyped) value or structure
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_time_mesh_dep_output(cs_post_time_mesh_dep_output_t  *function,
                                 void                            *input)
{
  /* Resize array of registered post-processings if necessary */

  if (_cs_post_n_output_mtp >= _cs_post_n_output_mtp_max) {
    if (_cs_post_n_output_mtp_max == 0)
      _cs_post_n_output_mtp_max = 8;
    else
      _cs_post_n_output_mtp_max *= 2;
    BFT_REALLOC(_cs_post_f_output_mtp,
                _cs_post_n_output_mtp_max,
                cs_post_time_mesh_dep_output_t *);
    BFT_REALLOC(_cs_post_i_output_mtp, _cs_post_n_output_mtp_max, void *);
  }

  /* Add a post-processing */

  _cs_post_f_output_mtp[_cs_post_n_output_mtp] = function;
  _cs_post_i_output_mtp[_cs_post_n_output_mtp] = input;

  _cs_post_n_output_mtp += 1;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
