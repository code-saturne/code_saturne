/*============================================================================
 * Handle export of mesh and fields.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2008  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_file.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_config_defs.h"
#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_writer.h"
#include "fvm_writer_priv.h"

/* Headers for available writers (could be replaced by plugin system) */

#include "fvm_to_cgns.h"
#include "fvm_to_med.h"
#include "fvm_to_ensight.h"
#include "fvm_to_ensight_v1.h"
#include "fvm_to_text.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 * Static and constant variables
 *============================================================================*/

/* Number and status of defined formats */

static const int _fvm_writer_n_formats = 5;

static fvm_writer_format_t _fvm_writer_format_list[5] = {

  /* Built-in text writer */
  {
    "text",
    VERSION,
    (  FVM_WRITER_FORMAT_HAS_POLYGON
     | FVM_WRITER_FORMAT_HAS_POLYHEDRON),
    FVM_WRITER_TRANSIENT_CONNECT,
    NULL,                              /* n_version_strings_func */
    NULL,                              /* version_string_func */
    fvm_to_text_init_writer,           /* init_func */
    fvm_to_text_finalize_writer,       /* finalize_func */
    NULL,                              /* set_mesh_time_func */
    NULL,                              /* needs_tesselation_func */
    fvm_to_text_export_nodal,          /* export_nodal_func */
    NULL,                              /* export_field_func */
    NULL                               /* flush_func */
  },

  /* Built-in EnSight Gold writer */
  {
    "EnSight Gold",
    "7.4 +",
    (  FVM_WRITER_FORMAT_HAS_POLYGON
     | FVM_WRITER_FORMAT_HAS_POLYHEDRON),
    FVM_WRITER_TRANSIENT_CONNECT,
    NULL,                              /* n_version_strings_func */
    NULL,                              /* version_string_func */
    fvm_to_ensight_init_writer,        /* init_func */
    fvm_to_ensight_finalize_writer,    /* finalize_func */
    fvm_to_ensight_set_mesh_time,      /* set_mesh_time_func */
    fvm_to_ensight_needs_tesselation,  /* needs_tesselation_func */
    fvm_to_ensight_export_nodal,       /* export_nodal_func */
    fvm_to_ensight_export_field,       /* export_field_func */
    NULL                               /* flush_func */
  },

  /* Built-in EnSight Gold writer, older version */
  {
    "EnSight Gold (v1 driver)",
    "7.4 +",
    (  FVM_WRITER_FORMAT_HAS_POLYGON
     | FVM_WRITER_FORMAT_HAS_POLYHEDRON),
    FVM_WRITER_TRANSIENT_CONNECT,
    NULL,                                 /* n_version_strings_func */
    NULL,                                 /* version_string_func */
    fvm_to_ensight_v1_init_writer,        /* init_func */
    fvm_to_ensight_v1_finalize_writer,    /* finalize_func */
    fvm_to_ensight_v1_set_mesh_time,      /* set_mesh_time_func */
    fvm_to_ensight_v1_needs_tesselation,  /* needs_tesselation_func */
    fvm_to_ensight_v1_export_nodal,       /* export_nodal_func */
    fvm_to_ensight_v1_export_field,       /* export_field_func */
    NULL                                  /* flush_func */
  },

  /* MED_fichier 2.3 writer */
  {
    "MED_fichier",
    "2.3 +",
    (  FVM_WRITER_FORMAT_USE_EXTERNAL
     | FVM_WRITER_FORMAT_HAS_POLYGON
     | FVM_WRITER_FORMAT_HAS_POLYHEDRON),
    FVM_WRITER_FIXED_MESH,
#if defined(HAVE_MED)
    NULL,                              /* n_version_strings_func */
    NULL,                              /* version_string_func */
    fvm_to_med_init_writer,            /* init_func */
    fvm_to_med_finalize_writer,        /* finalize_func */
    fvm_to_med_set_mesh_time,          /* set_mesh_time_func */
    fvm_to_med_needs_tesselation,      /* needs_tesselation_func */
    fvm_to_med_export_nodal,           /* export_nodal_func */
    fvm_to_med_export_field,           /* export_field_func */
    NULL                               /* flush_func */
#else
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
#endif
  },

  /* CGNS writer */
  {
    "CGNS",
    "2.5 +",
    (  FVM_WRITER_FORMAT_USE_EXTERNAL
     | FVM_WRITER_FORMAT_HAS_POLYGON),
    FVM_WRITER_FIXED_MESH,
#if defined(HAVE_CGNS)
    NULL,                              /* n_version_strings_func */
    NULL,                              /* version_string_func */
    fvm_to_cgns_init_writer,           /* init_func */
    fvm_to_cgns_finalize_writer,       /* finalize_func */
    fvm_to_cgns_set_mesh_time,         /* set_mesh_time_func */
    fvm_to_cgns_needs_tesselation,     /* needs_tesselation_func */
    fvm_to_cgns_export_nodal,          /* export_nodal_func */
    fvm_to_cgns_export_field,          /* export_field_func */
    NULL                               /* flush_func */
#else
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
#endif
  }

};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Find the format with the closest name to the argument,
 *
 * parameters:
 *   format_name <-- name of desired format
 *
 * returns:
 *   index of the format with the closest name to the argument,
 *   or _fvm_writer_n_formats if no name is similar enough.
 *----------------------------------------------------------------------------*/

static int
_fvm_writer_format_closest_name(const char  *const format_name)
{
  char  tmp_name[32], closest_name[32];
  int i, l;

  if (format_name == NULL)
    return _fvm_writer_n_formats;

  l = strlen(format_name);

  /* Transform format name to lowercase, whitespace as underscore */

  strncpy(tmp_name, format_name, 32);
  tmp_name[31] = '\0';
  for (i = 0 ; i < l ; i++) {
    tmp_name[i] = tolower(tmp_name[i]);
    if (tmp_name[i] == ' ' || tmp_name[i] == '\t')
      tmp_name[i] = '_';
  }

  /* Try "known" names */

  if (strncmp(tmp_name, "ensight", 7) == 0)
    strcpy(closest_name, "EnSight Gold");
  else if (strncmp(tmp_name, "med_mem", 7) == 0)
    strcpy(closest_name, "MED_memory");
  else if (strncmp(tmp_name, "med", 3) == 0)
    strcpy(closest_name, "MED_fichier");
  else if (strncmp(tmp_name, "cgns", 4) == 0)
    strcpy(closest_name, "CGNS");

  /* Find name in list */

  for (i = 0 ; i < _fvm_writer_n_formats ; i++)
    if (strcmp(closest_name, _fvm_writer_format_list[i].name) == 0)
      break;

  return i;
}

/*----------------------------------------------------------------------------
 * Transform a string containing a list of options to lowercase with
 * whitespace separators.
 * The new list is dynamically allocated, and should be freed when no
 * longer needed.
 *
 * parameters:
 *   option_list <-- options string (case-independent, whitespace,
 *                   semicolon, or comma separated list)
 *
 * returns:
 *   single-whitespace separated option string in lowercase.
 *----------------------------------------------------------------------------*/

static char *
_fvm_writer_option_list(const char  *const option_list)
{
  char *ret_list;
  int i, j, l;

  if (option_list == NULL)
    return NULL;

  l = strlen(option_list);

  BFT_MALLOC(ret_list, l + 1, char);

  /* Transform format name to lowercase, single whitespace separated */

  for (i = 0, j = 0 ; i < l ; i++) {
    ret_list[j] = tolower(option_list[i]);
    if (ret_list[j] == ',' || ret_list[j] == ';' || ret_list[j] == '\t')
      ret_list[j] = ' ';
    if (ret_list[j] != ' ' || (j > 0 && ret_list[j-1] != ' '))
      j++;
  }
  if (j > 0 && ret_list[j-1] == ' ')
    j--;

  ret_list[j] = '\0';

  return ret_list;
}

/*============================================================================
 * Semi-private function definitions (prototypes in fvm_writer_priv.h)
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute recommended buffer sizes to input or output a nodal mesh
 * definition by slices. This is especially useful when gathering the mesh for
 * output by slices using standard I/O in parallel mode.
 *
 * The global number of vertices and elements of each slice may also
 * be returned, if the pointers n_g_vertices and n_g_elements_section
 * are non-NULL respectively.
 *
 * The number of slices indicated is a minimum, and only a target;
 * computation is based primarily on cell and face connectivity, and the
 * target should be met for strided connectivities on those types of elements
 * only. Using an "optimistic" (i.e. small) mean number of vertices per
 * polyhedra or polygon will typically lead to requiring more slices, as
 * the connectivity slice size returned will be smaller than that truly
 * required for the corresponding slice size.
 * Slice sizes required for edges connectivity will meet the target only
 * when the global numbers of cells and faces given are zero, so as to
 * avoid generating too large connectivity slice sizes for cells should a mesh
 * contain both (as for example a hexahedral connectivity slice is 8 times
 * larger than the corresponding slice size, while an edges connectivity is
 * only 2 times as large).
 *
 * parameters:
 *   this_nodal                 <-- pointer to nodal mesh structure
 *   n_slices                   <-- target number of slices required
 *   n_polyhedron_vertices_mean <-- estimate of the mean number of vertices
 *                                  per polyhedron
 *   n_polygon_vertices_mean    <-- estimate of the mean number of vertices
 *                                  per polygon
 *   n_g_vertices               --> global number of vertices (or NULL)
 *   n_g_elements_section       --> array for global number of elements per
 *                                  section (or NULL)
 *   global_s_size              --> maximum number of entities defined per slice
 *   global_connect_s_size      --> maximum number of connectivity values
 *                                  per slice
 *----------------------------------------------------------------------------*/

void
fvm_writer_def_nodal_buf_size(const fvm_nodal_t  *this_nodal,
                              int                 n_slices,
                              int                 n_polyhedron_vertices_mean,
                              int                 n_polygon_vertices_mean,
                              fvm_gnum_t         *n_g_vertices,
                              fvm_gnum_t          n_g_elements_section[],
                              fvm_gnum_t         *global_s_size,
                              fvm_gnum_t         *global_connect_s_size)
{
  int  i;
  fvm_gnum_t  n_g_elements = 0;
  fvm_gnum_t  n_g_cells = 0, n_g_faces = 0, n_g_edges = 0;
  fvm_gnum_t  _n_g_vertices = 0;

  fvm_gnum_t  connect_size = 0;
  fvm_gnum_t  *_n_g_elements_section = NULL;

  const fvm_nodal_section_t  *section = NULL;

  assert(this_nodal != NULL);

  if (n_g_elements_section == NULL)
    BFT_MALLOC(_n_g_elements_section, this_nodal->n_sections, fvm_gnum_t);
  else
    _n_g_elements_section = n_g_elements_section;

  /* Vertices */

  if (this_nodal->global_vertex_num != NULL)
    _n_g_vertices = fvm_io_num_get_global_count(this_nodal->global_vertex_num);
  else
    _n_g_vertices = this_nodal->n_vertices;

  if (n_g_vertices != NULL)
    *n_g_vertices = _n_g_vertices;

  /* Edges, faces, and cells */

  for (i = 0; i < this_nodal->n_sections; i++) {

    section = this_nodal->sections[i];

    n_g_elements = fvm_nodal_section_n_g_elements(section);

    switch(section->entity_dim) {
    case 1:
      n_g_edges += n_g_elements;
      break;
    case 2:
      n_g_faces += n_g_elements;
      break;
    default:
      n_g_cells += n_g_elements;
    }

    _n_g_elements_section[i] = n_g_elements;

  }

  /* Define global slice size as (global_size / n_ranks) +  1 */

  *global_s_size = FVM_MAX(n_g_cells, n_g_faces);

  if (*global_s_size == 0)
    *global_s_size = n_g_edges;

  if (*global_s_size == 0) /* Should the mesh contain only vertices */
    *global_s_size = _n_g_vertices;

  *global_s_size = (*global_s_size / n_slices) + 1;

  /* Now compute best size for connectivity buffer,
     using estimated "mean" numbers of nodes per polyhedron and per polygon */

  *global_connect_s_size = 0;

  for (i = 0; i < this_nodal->n_sections; i++) {

    section = this_nodal->sections[i];

    switch(section->type) {
    case FVM_FACE_POLY:
    case FVM_CELL_POLY:
      if (section->type == FVM_FACE_POLY)
        connect_size =   FVM_MIN(_n_g_elements_section[i], *global_s_size)
                       * n_polygon_vertices_mean;
      else if (section->type == FVM_CELL_POLY)
        connect_size =   FVM_MIN(_n_g_elements_section[i], *global_s_size)
                       * n_polyhedron_vertices_mean;

      if (section->tesselation != NULL) {

        int  i_type;

        const fvm_tesselation_t  *tesselation = section->tesselation;
        const int  n_sub_types = fvm_tesselation_n_sub_types(tesselation);

        for (i_type = 0; i_type < n_sub_types; i_type++) {

          int  stride;
          fvm_lnum_t  n_sub_elements_max;

          fvm_element_t  sub_type = fvm_tesselation_sub_type(tesselation,
                                                             i_type);

          fvm_tesselation_get_global_size(tesselation,
                                          sub_type,
                                          NULL,
                                          &n_sub_elements_max);

          stride = fvm_nodal_n_vertices_element[sub_type];


          connect_size = FVM_MAX(connect_size,
                                 (fvm_gnum_t)n_sub_elements_max * stride);

        }

      }

      break;
      break;
    default:
      connect_size =   FVM_MIN(_n_g_elements_section[i], *global_s_size)
                     * section->stride;
    }

    /* Buffer size is that required for the entity type
       requiring the largest buffer size */

    *global_connect_s_size = FVM_MAX(*global_connect_s_size,
                                     connect_size);

  }

  if (n_g_elements_section != _n_g_elements_section)
    BFT_FREE(_n_g_elements_section);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Returns number of known formats.
 *----------------------------------------------------------------------------*/

int
fvm_writer_n_formats(void)
{
  return _fvm_writer_n_formats;
}

/*----------------------------------------------------------------------------
 * Returns name of a known format.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *
 * returns:
 *   pointer to constant string containing the format's name
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_format_name(int format_index)
{
  if (format_index >= 0 && format_index < _fvm_writer_n_formats)
    return _fvm_writer_format_list[format_index].name;

  else
    return NULL;
}

/*----------------------------------------------------------------------------
 * Returns availability of a known format.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *
 * returns:
 *   1 if the format is available, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
fvm_writer_format_available(int format_index)
{
  int retval = 0;

  if (format_index >= 0 && format_index < _fvm_writer_n_formats) {
    if (_fvm_writer_format_list[format_index].init_func != NULL)
      retval = 1;
  }
  return retval;
}

/*----------------------------------------------------------------------------
 * Returns number of library version strings associated with a given format.
 *
 * For writers requiring an external library, the first associated
 * version string should correspond to that library, with possible
 * additional version strings for its dependencies.
 *
 * For writers only requiring standard libraries (libc, MPI, MPI-IO),
 * this function should return 0.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *
 * returns:
 *   number of library version strings associated with a given format.
 *----------------------------------------------------------------------------*/

int
fvm_writer_n_version_strings(int format_index)
{
  int retval = 0;
  fvm_writer_n_version_strings_t  *n_version_strings_func = NULL;

  if (format_index >= 0 && format_index < _fvm_writer_n_formats) {
    n_version_strings_func
      = _fvm_writer_format_list[format_index].n_version_strings_func;
    if (n_version_strings_func != NULL)
      retval = n_version_strings_func(format_index);
  }
  return retval;
}

/*----------------------------------------------------------------------------
 * Returns a library version string associated with a given format.
 *
 * We must have string_index < fvm_writer_n_version_strings(format_index).
 *
 * In certain cases, when using dynamic libraries, fvm may be compiled
 * with one library version, and linked with another. If both run-time
 * and compile-time version information is available, this function
 * will return the run-time version string by default.
 *
 * Setting the compile_time flag to 1, the compile-time version string
 * will be returned if this is different from the run-time version.
 * If the version is the same, or only one of the 2 version strings are
 * available, a NULL character string will be returned with this flag set.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *   string_index <-- index in format's version string list (0 to n-1)
 *   compile_time <-- 0 by default, 1 if we want the compile-time version
 *                    string, if different from the run-time version.
 *
 * returns:
 *   pointer to constant string containing the library's version.
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_version_string(int format_index,
                          int string_index,
                          int compile_time_version)
{
  const char * retval = NULL;
  fvm_writer_version_string_t  *version_string_func = NULL;

  if (format_index >= 0 && format_index < _fvm_writer_n_formats) {
    version_string_func
      = _fvm_writer_format_list[format_index].version_string_func;
    if (version_string_func != NULL)
      retval = version_string_func(format_index,
                                   string_index,
                                   compile_time_version);
  }
  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize FVM mesh and field output writer.
 *
 * Allowed options depend on what is applicable to a given format. Those
 * not relevant to a given writer are ignored. Possible options include:
 *   text                output text files
 *   binary              output binary files (default)
 *   big_endian          force binary files to big-endian
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   divide_polyhedra    tesselate polyhedra with tetrahedra and pyramids
 *                       (adding a vertex near each polyhedron's center)
 *   split_tensors       write tensor values as separate scalars
 *
 * parameters:
 *   name            <-- base name of output
 *   path            <-- optional directory name for output
 *                       (directory automatically created if necessary)
 *   format_name     <-- name of selected format (case-independent)
 *   format_options  <-- options for the selected format (case-independent,
 *                       whitespace or comma separated list)
 *   time_dependency <-- indicates if and how meshes will change with time
 *
 * returns:
 *   pointer to mesh and field output writer
 *----------------------------------------------------------------------------*/

fvm_writer_t *
fvm_writer_init(const char             *name,
                const char             *path,
                const char             *format_name,
                const char             *format_options,
                fvm_writer_time_dep_t   time_dependency)
{
  int  i;
  char   local_dir[] = ".";
  char  *tmp_path = NULL;
  char  *tmp_options = NULL;
  fvm_writer_t  *this_writer = NULL;
  fvm_writer_init_t  *init_func = NULL;

  /* Find corresponding format and check coherency */

  for (i = 0 ; i < _fvm_writer_n_formats ; i++)
    if (strcmp(format_name, _fvm_writer_format_list[i].name) == 0)
      break;

  if (i >= _fvm_writer_n_formats)
    i = _fvm_writer_format_closest_name(format_name);

  if (i >= _fvm_writer_n_formats)
    bft_error(__FILE__, __LINE__, 0,
              _("Format type \"%s\" required for case \"%s\" is unknown"),
              format_name, name);

  if (!fvm_writer_format_available(i))
    bft_error(__FILE__, __LINE__, 0,
              _("Format type \"%s\" required for case \"%s\" is not available"),
              format_name, name);


  /* Create directory */

  if (path != NULL) {

    int l = strlen(path);

    if (l > 0) {
      BFT_MALLOC(tmp_path, l + 2, char);
      strcpy(tmp_path, path);
      if (tmp_path[l - 1] == FVM_DIR_SEPARATOR)
        tmp_path[l - 1] = '\0';
      if (bft_file_mkdir_default(path) == 1)
        tmp_path[0] = '\0';
      else {
        l = strlen(tmp_path);
        tmp_path[l]   = FVM_DIR_SEPARATOR;
        tmp_path[l+1] = '\0';
      }
    }

  }
  else
    tmp_path = local_dir;

  tmp_options = _fvm_writer_option_list(format_options);

  /* Also, temporary option to use legacy EnSight driver
     (we assume this driver comes just after the new EnSight
     driver in the formats list) */

  if (   tmp_options != NULL
      && !strcmp(_fvm_writer_format_list[i].name, "EnSight Gold")) {

    int i1, i2, l_opt;
    int l_tot = strlen(tmp_options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {
      for (i2 = i1; i2 < l_tot && tmp_options[i2] != ' '; i2++);
      l_opt = i2 - i1;
      if ((l_opt == 2) && (strncmp(tmp_options + i1, "v1", l_opt) == 0)) {
        i = i+1;
        break;
      }
      for (i1 = i2 + 1; i1 < l_tot && tmp_options[i1] == ' '; i1++);
    }

  }

  /* Initialize writer */

  BFT_MALLOC(this_writer, 1, fvm_writer_t);

  BFT_MALLOC(this_writer->name, strlen(name) + 1, char);
  strcpy(this_writer->name, name);

  this_writer->format = &(_fvm_writer_format_list[i]);

  this_writer->time_dep = FVM_MIN(time_dependency,
                                  this_writer->format->max_time_dep);

  this_writer->mesh_wtime = 0.;
  this_writer->mesh_cpu_time = 0.;
  this_writer->field_wtime = 0.;
  this_writer->field_cpu_time = 0.;

  /* Initialize format-specific writer */

  init_func = this_writer->format->init_func;

  if (init_func != NULL) {
#if defined(HAVE_MPI)
    this_writer->format_writer = init_func(name,
                                           tmp_path,
                                           tmp_options,
                                           this_writer->time_dep,
                                           fvm_parall_get_mpi_comm());
#else
    this_writer->format_writer = init_func(name,
                                           tmp_path,
                                           tmp_options,
                                           this_writer->time_dep);
#endif
  }
  else
    this_writer->format_writer = NULL;

  if (tmp_options != NULL)
    BFT_FREE(tmp_options);

  if (tmp_path != local_dir)
    BFT_FREE(tmp_path);

  /* Return pointer to initialized writer */

  return this_writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM mesh and field output writer.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_writer_t *
fvm_writer_finalize(fvm_writer_t  *this_writer)
{
  fvm_writer_finalize_t  *finalize_func = NULL;

  assert(this_writer != NULL);
  assert(this_writer->format != NULL);

  BFT_FREE(this_writer->name);

  finalize_func = this_writer->format->finalize_func;

  if (finalize_func != NULL)
    this_writer->format_writer = finalize_func(this_writer->format_writer);
  else
    this_writer->format_writer = NULL;

  BFT_FREE(this_writer);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Return a writer's name.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   pointer to base name of output associated with the writer
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_get_name(const fvm_writer_t  *this_writer)
{
  return this_writer->name;
}

/*----------------------------------------------------------------------------
 * Return a writer's associated format name.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   pointer to output format name associated with the writer
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_get_format(const fvm_writer_t  *this_writer)
{
  return this_writer->format->name;
}

/*----------------------------------------------------------------------------
 * Return geometry time dependency status of a writer.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
fvm_writer_get_time_dep(const fvm_writer_t  *this_writer)
{
  return this_writer->time_dep;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a mesh.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_writer_set_mesh_time(fvm_writer_t  *this_writer,
                         int            time_step,
                         double         time_value)
{
  fvm_writer_set_mesh_time_t  *set_mesh_time_func = NULL;

  assert(this_writer != NULL);
  assert(this_writer->format != NULL);

  set_mesh_time_func = this_writer->format->set_mesh_time_func;

  if (set_mesh_time_func != NULL)
    set_mesh_time_func(this_writer->format_writer,
                       time_step,
                       time_value);
}

/*----------------------------------------------------------------------------
 * Query if elements of a given type will need to be tesselated
 * for use of a nodal mesh with an output writer.
 *
 * This function should be called before any fvm_writer_export_...()
 *
 * parameters:
 *   this_writer  <-- pointer to mesh and field output writer
 *   mesh         <-- pointer to nodal mesh
 *   element_type <-- type of element
 *
 * returns:
 *   0 if no tesselation is necessary, 1 if tesselation is necessary.
 *----------------------------------------------------------------------------*/

int
fvm_writer_needs_tesselation(fvm_writer_t       *this_writer,
                             const fvm_nodal_t  *mesh,
                             fvm_element_t       element_type)
{
  int retval = 0;
  fvm_writer_needs_tesselation_t  *needs_tesselation_func = NULL;

  needs_tesselation_func = this_writer->format->needs_tesselation_func;
  if (needs_tesselation_func != NULL)
    retval = needs_tesselation_func(this_writer->format_writer,
                                    mesh,
                                    element_type);
  return retval;
}

/*----------------------------------------------------------------------------
 * Export FVM nodal mesh.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *   mesh        <-- pointer to nodal mesh
 *----------------------------------------------------------------------------*/

void
fvm_writer_export_nodal(fvm_writer_t        *this_writer,
                        const fvm_nodal_t   *mesh)
{
  double w_start, w_end, cpu_start, cpu_end;

  fvm_writer_export_nodal_t  *export_nodal_func = NULL;

  assert(this_writer != NULL);
  assert(this_writer->format != NULL);

  w_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  export_nodal_func = this_writer->format->export_nodal_func;

  if (export_nodal_func != NULL)
    export_nodal_func(this_writer->format_writer,
                      mesh);

  w_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  this_writer->mesh_wtime += (w_end - w_start);
  this_writer->mesh_cpu_time += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------
 * Export field associated with a nodal mesh.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer      <-- pointer to mesh and field output writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_writer_export_field(fvm_writer_t                 *this_writer,
                        const fvm_nodal_t            *mesh,
                        const char                   *name,
                        fvm_writer_var_loc_t          location,
                        int                           dimension,
                        fvm_interlace_t               interlace,
                        int                           n_parent_lists,
                        const fvm_lnum_t              parent_num_shift[],
                        fvm_datatype_t                datatype,
                        int                           time_step,
                        double                        time_value,
                        const void             *const field_values[])
{
  double w_start, w_end, cpu_start, cpu_end;

  fvm_writer_export_field_t  *export_field_func = NULL;

  assert(this_writer != NULL);
  assert(this_writer->format != NULL);

  w_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  export_field_func = this_writer->format->export_field_func;

  if (export_field_func != NULL)
    export_field_func(this_writer->format_writer,
                      mesh,
                      name,
                      location,
                      dimension,
                      interlace,
                      n_parent_lists,
                      parent_num_shift,
                      datatype,
                      time_step,
                      time_value,
                      field_values);

  w_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  this_writer->field_wtime += (w_end - w_start);
  this_writer->field_cpu_time += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * parameters:
 *   this_writer      <-- pointer to mesh and field output writer
 *----------------------------------------------------------------------------*/

void
fvm_writer_flush(fvm_writer_t  *this_writer)
{
  fvm_writer_flush_t  *flush_func = NULL;

  assert(this_writer != NULL);
  assert(this_writer->format != NULL);

  flush_func = this_writer->format->flush_func;

  if (flush_func != NULL)
    flush_func(this_writer->format_writer);
}

/*----------------------------------------------------------------------------
 * Return accumulated wall-clock and CPU times associated with mesh and
 * field exports for a given writer.
 *
 * parameters:
 *   this_writer      <-- pointer to mesh and field output writer
 *   mesh_wtime       --> Meshes output Wall-clock time (or NULL)
 *   mesh_cpu_time    --> Meshes output CPU time (or NULL)
 *   field_wtime      --> Fields output Wall-clock time (or NULL)
 *   field_cpu_time   --> Fields output CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_writer_get_times(fvm_writer_t  *this_writer,
                     double        *mesh_wtime,
                     double        *mesh_cpu_time,
                     double        *field_wtime,
                     double        *field_cpu_time)
{
  assert(this_writer != NULL);

  if (mesh_wtime != NULL)
    *mesh_wtime = this_writer->mesh_wtime;
  if (mesh_cpu_time != NULL)
    *mesh_cpu_time = this_writer->mesh_cpu_time;
  if (field_wtime != NULL)
    *field_wtime = this_writer->field_wtime;
  if (field_cpu_time != NULL)
    *field_cpu_time = this_writer->field_cpu_time;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
