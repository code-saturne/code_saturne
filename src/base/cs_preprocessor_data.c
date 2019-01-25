/*============================================================================
 * Manage the exchange of data between Code_Saturne and the pre-processor
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_mesh_from_builder.h"
#include "cs_mesh_group.h"
#include "cs_parall.h"
#include "cs_partition.h"
#include "cs_io.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_preprocessor_data.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  const char         *filename;   /* File name */
  cs_file_off_t       offset;     /* File offsets for re-opening */
  const double       *matrix;     /* Coordinate transformation matrix */

  size_t              n_group_renames;
  const char  *const *old_group_names;
  const char  *const *new_group_names;

  /* Single allocation for all data */

  size_t              data_size;
  unsigned char      *data;

} _mesh_file_info_t;

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  /* File info */

  int                 n_files;
  _mesh_file_info_t  *file_info;

  /* Temporary dimensions necessary for multiple inputs */

  int         *gc_id_shift;

  int          n_perio_read;
  cs_lnum_t    n_cells_read;
  cs_lnum_t    n_faces_read;
  cs_lnum_t    n_faces_connect_read;
  cs_lnum_t    n_vertices_read;

  cs_gnum_t    n_g_cells_read;
  cs_gnum_t    n_g_faces_read;
  cs_gnum_t    n_g_faces_connect_read;
  cs_gnum_t    n_g_vertices_read;

} _mesh_reader_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static _mesh_reader_t *_cs_glob_mesh_reader = NULL;

#if defined(WIN32) || defined(_WIN32)
static const char _dir_separator = '\\';
#else
static const char _dir_separator = '/';
#endif

/* Definitions of file to read */

int _n_mesh_files = 0;
int _n_max_mesh_files = 0;
_mesh_file_info_t  *_mesh_file_info = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return minimum size required to align data with the base pointer size.
 *
 * parameters:
 *   min_size <-- minimum data size
 *
 * returns:
 *   the data size required including alignment
 *----------------------------------------------------------------------------*/

static inline size_t
_align_size(size_t  min_size)
{
  const size_t align = (sizeof(void *))/sizeof(unsigned char);
  return (min_size + (align-1) - ((min_size - 1) % align));
}

/*----------------------------------------------------------------------------
 * Define defaul input data in nothing has been specified by the user.
 *----------------------------------------------------------------------------*/

static void
_set_default_input_if_needed(void)
{
  const char input_default[] = "mesh_input";
  const char cp_input_default[] = "restart/mesh_input";

  if (_n_mesh_files == 0) {

    if (cs_file_isreg(input_default))
      cs_preprocessor_data_add_file(input_default, 0, NULL, NULL);

    else if (cs_file_isdir(input_default)) {
      int i;
      char **dir_files = cs_file_listdir(input_default);
      for (i = 0; dir_files[i] != NULL; i++) {
        char *tmp_name = NULL;
        BFT_MALLOC(tmp_name,
                   strlen(input_default) + 1 + strlen(dir_files[i]) + 1,
                   char);
        sprintf(tmp_name, "%s%c%s",
                input_default, _dir_separator, dir_files[i]);
        if (cs_file_isreg(tmp_name))
          cs_preprocessor_data_add_file(tmp_name, 0, NULL, NULL);
        BFT_FREE(tmp_name);
        BFT_FREE(dir_files[i]);
      }
      BFT_FREE(dir_files);
    }

    else if (cs_file_isreg(cp_input_default))
      cs_preprocessor_data_add_file(cp_input_default, 0, NULL, NULL);

    else
      bft_error(__FILE__, __LINE__, 0,
                _("No \"%s\" file or directory found."), input_default);
  }
}

/*----------------------------------------------------------------------------
 * Create an empty mesh reader helper structure.
 *
 * Property of the matching mesh file information structure is transferred
 * to the mesh reader.
 *
 * parameters:
 *   n_mesh_files   <-> number of associated mesh files
 *   mesh_file_info <-> array of mesh file information structures
 *
 * returns:
 *   A pointer to a mesh reader helper structure
 *----------------------------------------------------------------------------*/

static _mesh_reader_t *
_mesh_reader_create(int                 *n_mesh_files,
                    _mesh_file_info_t  **mesh_file_info)
{
  int i;
  _mesh_reader_t  *mr = NULL;

  BFT_MALLOC(mr, 1, _mesh_reader_t);

  memset(mr, 0, sizeof(_mesh_reader_t));

  /* Transfer ownership of mesh file info */

  mr->n_files = *n_mesh_files;
  mr->file_info = *mesh_file_info;

  BFT_REALLOC(mr->file_info, mr->n_files, _mesh_file_info_t);

  /* Setup remaining structure fields */

  *n_mesh_files = 0;
  *mesh_file_info = NULL;

  BFT_MALLOC(mr->gc_id_shift, mr->n_files, int);
  for (i = 0; i < mr->n_files; i++)
    mr->gc_id_shift[i] = 0;

  mr->n_perio_read = 0;
  mr->n_cells_read = 0;
  mr->n_faces_read = 0;
  mr->n_faces_connect_read = 0;
  mr->n_vertices_read = 0;

  mr->n_g_cells_read = 0;
  mr->n_g_faces_read = 0;
  mr->n_g_faces_connect_read = 0;

  return mr;
}

/*----------------------------------------------------------------------------
 * Destroy a mesh reader helper structure
 *
 * parameters:
 *   mr <-> pointer to a mesh reader helper
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static void
_mesh_reader_destroy(_mesh_reader_t  **mr)
{
  int i;
  _mesh_reader_t *_mr = *mr;

  for (i = 0; i < _mr->n_files; i++) {
    _mesh_file_info_t  *f = _mr->file_info + i;
    BFT_FREE(f->data);
  }
  BFT_FREE(_mr->file_info);

  BFT_FREE(_mr->gc_id_shift);

  BFT_FREE(*mr);
}

/*----------------------------------------------------------------------------
 * Add a periodicity to mesh->periodicities (fvm_periodicity_t *) structure.
 *
 * parameters:
 *   mesh       <-> mesh
 *   perio_type <-- periodicity type
 *   perio_num  <-- periodicity number (identifier)
 *   matrix     <-- transformation matrix using homogeneous coordinates
 *----------------------------------------------------------------------------*/

static void
_add_periodicity(cs_mesh_t *mesh,
                 cs_int_t   perio_type,
                 cs_int_t   perio_num,
                 cs_real_t  matrix[3][4])
{
  cs_int_t  i, j;
  double  _matrix[3][4];

  fvm_periodicity_type_t _perio_type = (fvm_periodicity_type_t)perio_type;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++)
      _matrix[i][j] = matrix[i][j];
  }

  fvm_periodicity_add_by_matrix(mesh->periodicity,
                                perio_num,
                                _perio_type,
                                _matrix);
}

/*----------------------------------------------------------------------------
 * Set block ranges for parallel reads
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   mb   <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

static void
_set_block_ranges(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb)
{
  int i;

  int min_block_size = 0;
  int rank_id = cs_glob_rank_id;
  int n_ranks = cs_glob_n_ranks;

#if defined(HAVE_MPI)
  int block_rank_step = 1;
  cs_file_get_default_comm(&block_rank_step, &min_block_size, NULL, NULL);
  mb->min_rank_step = block_rank_step;
#endif

  /* Always build per_face_range in case of periodicity */

  if (mb->n_perio > 0) {
    BFT_REALLOC(mb->per_face_bi, mb->n_perio, cs_block_dist_info_t);
    memset(mb->per_face_bi, 0, sizeof(cs_block_dist_info_t)*mb->n_perio);
  }

  /* Set block sizes and ranges (useful for parallel mode) */

  mb->cell_bi = cs_block_dist_compute_sizes(rank_id,
                                            n_ranks,
                                            mb->min_rank_step,
                                            min_block_size/sizeof(cs_gnum_t),
                                            mesh->n_g_cells);

  mb->face_bi = cs_block_dist_compute_sizes(rank_id,
                                            n_ranks,
                                            mb->min_rank_step,
                                            min_block_size/(sizeof(cs_gnum_t)*2),
                                            mb->n_g_faces);

  mb->vertex_bi = cs_block_dist_compute_sizes(rank_id,
                                              n_ranks,
                                              mb->min_rank_step,
                                              min_block_size/(sizeof(cs_real_t)*3),
                                              mesh->n_g_vertices);

  for (i = 0; i < mb->n_perio; i++)
    mb->per_face_bi[i]
      = cs_block_dist_compute_sizes(rank_id,
                                    n_ranks,
                                    mb->min_rank_step,
                                    min_block_size/sizeof(cs_gnum_t),
                                    mb->n_g_per_face_couples[i]);
}

/*----------------------------------------------------------------------------
 * Rename groups in a mesh.
 *
 * parameters:
 *   mesh            <-> mesh being modified
 *   start_id        <-- id of first group to rename
 *   n_group_renames <-- number of group rename couples
 *   old_group_names <-- old group names (size: n_group_renames)
 *   new_group_names <-- new group names (size: n_group_renames)
 *----------------------------------------------------------------------------*/

static void
_mesh_groups_rename(cs_mesh_t          *mesh,
                    size_t              start_id,
                    size_t              n_group_renames,
                    const char  *const *old_group_names,
                    const char  *const *new_group_names)
{
  size_t  i, j, k;
  int  have_rename = 0;
  size_t  end_id = mesh->n_groups;
  size_t  n_ids = 0;
  int  *new_group_id = NULL;

  if (end_id > start_id)
    n_ids = end_id - start_id;

  /* Check for rename matches */

  BFT_MALLOC(new_group_id, n_ids, int);

  for (i = 0, j = start_id; i < n_ids; i++, j++) {
    const char *g_name = mesh->group + mesh->group_idx[j];
    new_group_id[i] = -1;
    for (k = 0; k < n_group_renames; k++) {
      if (strcmp(g_name, old_group_names[k]) == 0) {
        new_group_id[i] = k;
        have_rename = 1;
        break;
      }
    }
  }

  /* Now rename matched groups */

  if (have_rename) {

    size_t new_size = 0;
    size_t old_size = mesh->group_idx[end_id] - mesh->group_idx[start_id];
    int *saved_idx = NULL;
    char *saved_names = NULL;

    BFT_MALLOC(saved_idx, n_ids + 1, int);
    BFT_MALLOC(saved_names, old_size, char);

    for (i = 0; i < n_ids+1; i++)
      saved_idx[i] = mesh->group_idx[start_id + i] - mesh->group_idx[start_id];
    memcpy(saved_names,
           mesh->group + mesh->group_idx[start_id],
           old_size);

    /* Update index */

    for (i = 0; i < n_ids; i++) {
      const char *new_src = NULL;
      new_size += 1;
      if (new_group_id[i] > -1)
        new_src = new_group_names[new_group_id[i]];
      else
        new_src = saved_names + saved_idx[i];
      if (new_src != NULL)
        new_size += strlen(new_src);
      mesh->group_idx[start_id + i + 1]
        = mesh->group_idx[start_id] + new_size;
    }

    BFT_REALLOC(mesh->group, mesh->group_idx[mesh->n_groups], char);

    for (i = 0, j = start_id; i < n_ids; i++, j++) {
      char *new_dest = mesh->group + mesh->group_idx[j];
      const char *new_src = NULL;
      if (new_group_id[i] > -1)
        new_src = new_group_names[new_group_id[i]];
      else
        new_src = saved_names + saved_idx[i];
      if (new_src != NULL)
        strcpy(new_dest, new_src);
      else
        strcpy(new_dest, "");
    }

    BFT_FREE(saved_names);
    BFT_FREE(saved_idx);

    /* Set mesh modification flag */

    mesh->modified = 1;

  }

  BFT_FREE(new_group_id);
}

/*----------------------------------------------------------------------------
 * Read sections from the pre-processor about the dimensions of mesh
 *
 * This function updates the information in the mesh and mesh reader
 * structures relative to the data in the given file.
 *
 * parameters
 *   mesh <-> pointer to mesh structure
 *   n_gc <-- number of group classes last read
 *----------------------------------------------------------------------------*/

static void
_colors_to_groups(cs_mesh_t  *mesh,
                  int         n_gc)
{
  cs_int_t  i, j;
  int  n_colors = 0;
  int  color_names_size = 0;

  /* Counting pass */

  for (j = 0; j < mesh->n_max_family_items; j++) {
    for (i = mesh->n_families - n_gc; i < mesh->n_families; i++) {
      if (mesh->family_item[mesh->n_families*j + i] > 0) {
        int color_id = mesh->family_item[mesh->n_families*j + i];
        int name_size = 1;
        while (color_id > 0) {
          color_id /= 10;
          name_size += 1;
        }
        n_colors += 1;
        color_names_size += name_size;
      }
    }
  }

  /* Reallocation */

  if (n_colors > 0) {
    if (mesh->n_groups > 0) {
      BFT_REALLOC(mesh->group_idx, mesh->n_groups + n_colors + 1, cs_int_t);
      BFT_REALLOC(mesh->group,
                  mesh->group_idx[mesh->n_groups] + color_names_size,
                  char);
    }
    else {
      BFT_MALLOC(mesh->group_idx, n_colors + 1, cs_int_t);
      BFT_MALLOC(mesh->group, color_names_size, char);
      mesh->group_idx[0] = 0;
    }
  }

  /* Assignment */

  for (j = 0; j < mesh->n_max_family_items; j++) {
    for (i = mesh->n_families - n_gc; i < mesh->n_families; i++) {
      if (mesh->family_item[mesh->n_families*j + i] > 0) {
        int color_id = mesh->family_item[mesh->n_families*j + i];
        int name_size = 1;
        cs_int_t group_end = mesh->group_idx[mesh->n_groups];
        sprintf(mesh->group + group_end, "%d", color_id);
        while (color_id > 0) {
          color_id /= 10;
          name_size += 1;
        }
        mesh->group_idx[mesh->n_groups + 1]
          = mesh->group_idx[mesh->n_groups] + name_size;
        mesh->n_groups += 1;
        mesh->family_item[mesh->n_families*j + i] = - mesh->n_groups;
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Read sections from the pre-processor about the dimensions of mesh
 *
 * This function updates the information in the mesh and mesh reader
 * structures relative to the data in the given file.
 *
 * parameters
 *   filename <-- file name
 *
 * returns:
 *   on rank 0, 1 if periodicity is present, 2 if rotational periodicity
 *   is present; 0 in all other cases
 *----------------------------------------------------------------------------*/

static int
_read_perio_info(const char  *filename)
{
  cs_io_sec_header_t  header;

  cs_io_t   *pp_in = NULL;

  int retval = 0;

  /* Initialize reading of Preprocessor output */

  bft_printf(_(" Checking metadata from file: \"%s\"\n"), filename);

#if defined(HAVE_MPI)
  pp_in = cs_io_initialize(filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           CS_FILE_STDIO_SERIAL,
                           CS_IO_ECHO_NONE,
                           MPI_INFO_NULL,
                           MPI_COMM_NULL,
                           MPI_COMM_NULL);
#else
  pp_in = cs_io_initialize(filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           CS_FILE_STDIO_SERIAL,
                           CS_IO_ECHO_NONE);
#endif

  /* Loop on read sections */

  while (true) {

    /* Receive headers and clean header names */

    cs_io_read_header(pp_in, &header);

    /* Treatment according to the header name */

    if (strncmp(header.sec_name, "EOF", CS_IO_NAME_LEN) == 0)
      break;

    /* Additional sections for periodicity. */

    else {

      if (strncmp(header.sec_name, "n_periodic_directions",
                  CS_IO_NAME_LEN) == 0)
        retval = 1;
      else if (strncmp(header.sec_name, "n_periodic_rotations",
                       CS_IO_NAME_LEN) == 0) {
        retval = 2;
        break;
      }
      else if (strncmp(header.sec_name, "end_block:dimensions",
                       CS_IO_NAME_LEN) == 0)
        break;

      cs_io_skip(&header, pp_in);

    }

  }

  cs_io_finalize(&pp_in);
  pp_in = NULL;

  /* Close file */

  cs_io_finalize(&pp_in);

  /* Return test values */

  return retval;
}

/*----------------------------------------------------------------------------
 * Read sections from the pre-processor about the dimensions of mesh
 *
 * This function updates the information in the mesh and mesh reader
 * structures relative to the data in the given file.
 *
 * parameters
 *   mesh     <-> pointer to mesh structure
 *   mb       <-> pointer to mesh builder helper structure
 *   mr       <-> pointer to mesh reader structure
 *   file_id  <-- file id in mesh reader
 *----------------------------------------------------------------------------*/

static void
_read_dimensions(cs_mesh_t          *mesh,
                 cs_mesh_builder_t  *mb,
                 _mesh_reader_t     *mr,
                 int                 file_id)
{
  cs_int_t  i, j;
  cs_io_sec_header_t  header;

  cs_gnum_t  n_elts = 0;
  int        n_gc = 0;
  int        n_gc_props_max = 0;
  int        n_groups = 0;
  int        n_init_perio = 0;
  bool       dim_read = false;
  bool       end_read = false;
  cs_io_t   *pp_in = NULL;

  _mesh_file_info_t  *f = NULL;

  const char  *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                   "unexpected or of incorrect size");

  if (file_id < 0 || file_id >= mr->n_files)
    return;

  f = mr->file_info + file_id;
  f->offset = 0;

  mr->gc_id_shift[file_id] = mesh->n_families;

  /* Initialize reading of Preprocessor output */

  bft_printf(_(" Reading metadata from file: \"%s\"\n"), f->filename);

#if defined(HAVE_MPI)
  pp_in = cs_io_initialize(f->filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           CS_FILE_STDIO_SERIAL,
                           CS_IO_ECHO_NONE,
                           MPI_INFO_NULL,
                           MPI_COMM_NULL,
                           cs_glob_mpi_comm);
#else
  pp_in = cs_io_initialize(f->filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           CS_FILE_STDIO_SERIAL,
                           CS_IO_ECHO_NONE);
#endif

  /* Loop on read sections */

  while (end_read == false) {

    /* Receive headers and clean header names */

    cs_io_read_header(pp_in, &header);

    /* Treatment according to the header name */

    if (strncmp(header.sec_name, "EOF", CS_IO_NAME_LEN)
        == 0) {
      cs_io_finalize(&pp_in);
      pp_in = NULL;
    }

    if (strncmp(header.sec_name, "start_block:dimensions",
                CS_IO_NAME_LEN) == 0) {

      if (dim_read == false)
        dim_read = true;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }
    else if (strncmp(header.sec_name, "end_block:dimensions",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read == true) {
        dim_read = false;
        end_read = true;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }

    /* Receive dimensions from the pre-processor */

    else if (strncmp(header.sec_name, "n_cells",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_cells;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_cells, pp_in);
        mesh->n_g_cells += _n_g_cells;
      }

    }
    else if (strncmp(header.sec_name, "n_faces",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_faces;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_faces, pp_in);
        mb->n_g_faces += _n_g_faces;
      }

    }
    else if (strncmp(header.sec_name, "n_vertices",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_vertices;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_vertices, pp_in);
        mesh->n_g_vertices += _n_g_vertices;
      }

    }
    else if (strncmp(header.sec_name, "face_vertices_size",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_face_connect_size;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_face_connect_size, pp_in);
        mb->n_g_face_connect_size += _n_g_face_connect_size;
      }

    }
    else if (strncmp(header.sec_name, "n_group_classes",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_gc, pp_in);
        mesh->n_families += n_gc;
      }

    }
    else if (strncmp(header.sec_name, "n_group_class_props_max",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_gc_props_max, pp_in);
        if (n_gc_props_max > mesh->n_max_family_items) {
          /* Update (pad) previous definitions */
          BFT_REALLOC(mesh->family_item,
                      mesh->n_families*n_gc_props_max,
                      cs_int_t);
          for (i = mesh->n_max_family_items;
               i < n_gc_props_max;
               i++) {
            for (j = 0; j < mesh->n_families - n_gc; j++)
              mesh->family_item[(mesh->n_families - n_gc)*i + j] = 0;
          }
          mesh->n_max_family_items = n_gc_props_max;
        }
      }

    }
    else if (strncmp(header.sec_name, "n_groups",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_groups, pp_in);
        mesh->n_groups += n_groups;
      }

    }
    else if (strncmp(header.sec_name, "group_name_index",
                     CS_IO_NAME_LEN) == 0) {

      if ((cs_int_t)header.n_vals != n_groups + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_lnum_t *_group_idx;
        BFT_MALLOC(_group_idx, n_groups + 1, cs_lnum_t);
        cs_io_read_global(&header, _group_idx, pp_in);
        if (mesh->group_idx == NULL) {
          BFT_MALLOC(mesh->group_idx, mesh->n_groups + 1, int);
          for (i = 0; i < n_groups+1; i++)
            mesh->group_idx[i] = _group_idx[i] - 1;
        }
        else {
          BFT_REALLOC(mesh->group_idx, mesh->n_groups + 1, int);
          for (i = 0, j = mesh->n_groups - n_groups; i < n_groups; i++, j++)
            mesh->group_idx[j + 1] = (   mesh->group_idx[j]
                                      + _group_idx[i+1] - _group_idx[i]);
        }
        BFT_FREE(_group_idx);
      }

    }
    else if (strncmp(header.sec_name, "group_name",
                     CS_IO_NAME_LEN) == 0) {

      if (   mesh->group_idx == NULL
          || (  (cs_int_t)header.n_vals
              != (  mesh->group_idx[mesh->n_groups]
                  - mesh->group_idx[mesh->n_groups - n_groups])))
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        i = mesh->group_idx[mesh->n_groups - n_groups] - mesh->group_idx[0];
        BFT_REALLOC(mesh->group, i + header.n_vals + 1, char);
        cs_io_read_global(&header, mesh->group + i, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "group_class_properties",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_families * mesh->n_max_family_items;
      if (dim_read != true || header.n_vals != n_gc*n_gc_props_max)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        if (mesh->family_item == NULL)
          BFT_MALLOC(mesh->family_item, n_elts, cs_int_t);
        if (mesh->n_families == n_gc)
          cs_io_read_global(&header, mesh->family_item, pp_in);
        else {
          cs_int_t *_family_item = NULL;
          BFT_REALLOC(mesh->family_item, n_elts, cs_int_t);
          BFT_MALLOC(_family_item, header.n_vals, cs_int_t);
          cs_io_read_global(&header, _family_item, pp_in);
          /* Shift previous data */
          for (j = mesh->n_max_family_items - 1; j > 0; j--) {
            for (i = mesh->n_families - n_gc - 1; i > -1; i--)
              mesh->family_item[mesh->n_families*j + i]
                = mesh->family_item[(mesh->n_families - n_gc)*j + i];
          }
          for (i = 0; i < n_gc; i++) {
            /* Copy group class data, shifting group names if necessary */
            for (j = 0; j < n_gc_props_max; j++) {
              int _family_item_j = _family_item[n_gc*j + i];
              if (_family_item_j < 0)
                _family_item_j -= (mesh->n_groups - n_groups);
              mesh->family_item[mesh->n_families*j + (mesh->n_families-n_gc+i)]
                = _family_item_j;
            }
            /* Pad if necessary */
            for (j = n_gc_props_max; j < mesh->n_max_family_items; j++)
              mesh->family_item[mesh->n_families*j + (mesh->n_families-n_gc+i)]
                = 0;
          }
          BFT_FREE(_family_item);
        }
        /* Transform colors to group names if present */
        _colors_to_groups(mesh, n_gc);
      }

      if (f->n_group_renames > 0)
        _mesh_groups_rename(mesh,
                            mesh->n_groups - n_groups,
                            f->n_group_renames,
                            f->old_group_names,
                            f->new_group_names);

    }

    /* Additional sections for periodicity. */

    else if (strncmp(header.sec_name, "n_periodic_directions",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {

        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_init_perio, pp_in);
        mesh->n_init_perio += n_init_perio;

        assert(n_init_perio > 0);

        if (mesh->periodicity == NULL)
          mesh->periodicity = fvm_periodicity_create(0.001);

        BFT_REALLOC(mb->n_per_face_couples, mesh->n_init_perio, cs_lnum_t);
        BFT_REALLOC(mb->n_g_per_face_couples, mesh->n_init_perio, cs_gnum_t);
        BFT_REALLOC(mb->per_face_couples, mesh->n_init_perio, cs_gnum_t *);

        mb->n_perio = mesh->n_init_perio;

        for (i = mesh->n_init_perio - n_init_perio;
             i < mesh->n_init_perio;
             i++) {
          mb->n_per_face_couples[i] = 0;
          mb->n_g_per_face_couples[i] = 0;
          mb->per_face_couples[i] = NULL;
        }
      }

    }
    else if (strncmp(header.sec_name, "n_periodic_rotations",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_int_t n_rot_perio;
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_rot_perio, pp_in);
        if (n_rot_perio > 0)
          mesh->have_rotation_perio = 1;
      }

    }
    else if (strncmp(header.sec_name, "n_periodic_faces",
                     CS_IO_NAME_LEN) == 0) {

      if ((cs_int_t)header.n_vals != n_init_perio)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        size_t dest_offset = mesh->n_init_perio - n_init_perio;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header,
                          mb->n_g_per_face_couples + dest_offset,
                          pp_in);
        for (i = dest_offset; i < mesh->n_init_perio; i++)
          mb->n_g_per_face_couples[i] /= 2;
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

  } /* End of test on headers */

  /* Close file */

  f->offset = cs_io_get_offset(pp_in);
  cs_io_finalize(&pp_in);
}

/*----------------------------------------------------------------------------
 * Compute value range information for a given data section.
 *
 * parameters:
 *   header          <-- pointer to current file section header data
 *   pp_in           <-- pointer to current file
 *   n_g_elts        <-> global number of elements
 *   n_g_elts_read   <-> global number of elements already read
 *   n_location_vals <-- number of values for each location
 *   is_index        <-- 1 if data is an index, 0 otherwise
 *   gnum_range      <-- global number range for all elements on this rank
 *   gnum_range_cur  --> global number range for elements in current file
 *   n_g_elts_cur    --> global number of elements in current file
 *   n_vals          --> expected number of local values from all files
 *   n_vals_cur      --> number of local values from current file
 *----------------------------------------------------------------------------*/

static void
_data_range(cs_io_sec_header_t  *header,
            const cs_io_t       *pp_in,
            cs_gnum_t            n_g_elts,
            cs_gnum_t            n_g_elts_read,
            size_t               n_location_vals,
            size_t               is_index,
            const cs_gnum_t      gnum_range[2],
            cs_gnum_t            gnum_range_cur[2],
            cs_gnum_t           *n_g_elts_cur,
            cs_lnum_t           *n_vals,
            cs_lnum_t           *n_vals_cur)
{
  size_t i;

  /* Initialization */

  gnum_range_cur[0] = gnum_range[0];
  gnum_range_cur[1] = gnum_range[1];

  *n_g_elts_cur = (header->n_vals - is_index) / n_location_vals;
  *n_vals = (gnum_range[1] - gnum_range[0]) * n_location_vals;
  *n_vals_cur = 0;

  if (*n_g_elts_cur + n_g_elts_read > n_g_elts)
    bft_error(__FILE__, __LINE__, 0,
              _("Section of type <%s> on <%s>\n"
                "has incorrect size (current: %llu, read: %llu, total: %llu."),
              header->sec_name, cs_io_get_name(pp_in),
              (unsigned long long)(*n_g_elts_cur),
              (unsigned long long)n_g_elts_read,
              (unsigned long long)n_g_elts);

  else if (header->n_location_vals != n_location_vals)
    bft_error(__FILE__, __LINE__, 0,
              _("Section of type <%s> on <%s>\n"
                "has incorrect number of values per location."),
              header->sec_name, cs_io_get_name(pp_in));

  else {

    /* Compute range for this file (based on range for all files,
       and parts of this range already read) */

    for (i = 0; i < 2; i++) {
      if (gnum_range_cur[i] <= n_g_elts_read)
        gnum_range_cur[i] = 1;
      else
        gnum_range_cur[i] -= n_g_elts_read;
      if (gnum_range_cur[i] > *n_g_elts_cur)
        gnum_range_cur[i] = *n_g_elts_cur + 1;
    }

    if (gnum_range[1] > gnum_range[0])
      *n_vals_cur = (gnum_range_cur[1] - gnum_range_cur[0]) * n_location_vals;
  }

  /* Index adds past-the-last value */
  if (is_index == 1) {
    *n_vals += 1;
    *n_vals_cur += 1;
  }
}

/*----------------------------------------------------------------------------
 * Transform coordinates.
 *
 * parameters:
 *   n_coords <-- number of coordinates.
 *   coords   <-> array of coordinates
 *   matrix   <-- matrix of the transformation in homogeneous coord.
 *                last line = [0; 0; 0; 1] (Not used here)
 *----------------------------------------------------------------------------*/

static void
_transform_coords(size_t         n_coords,
                  cs_real_t     *coords,
                  const double   matrix[])

{
  size_t  i;
  double  _matrix[3][4];
  memcpy(_matrix, matrix, 12*sizeof(double));

  for (i = 0; i < n_coords; i++) {

    size_t j, coo_shift;
    cs_real_t  tmp_coord[3];

    coo_shift = i*3;

    for (j = 0; j < 3; j++)
      tmp_coord[j] = coords[coo_shift + j];

    for (j = 0; j < 3; j++) {
      coords[coo_shift + j] = (  _matrix[j][0] * tmp_coord[0]
                               + _matrix[j][1] * tmp_coord[1]
                               + _matrix[j][2] * tmp_coord[2]
                               + _matrix[j][3]);
    }

  }
}

/*----------------------------------------------------------------------------
 * Invert a homogeneous transformation matrix.
 *
 * parameters:
 *   a  <-- matrix
 *   b  --> inverse matrix
 *----------------------------------------------------------------------------*/

static void
_inverse_transf_matrix(double  a[3][4],
                       double  b[3][4])
{
  int i, j, k, k_pivot;

  double abs_pivot, abs_a_ki, factor;
  double _a[3][4];
  double _a_swap[3], _b_swap[3];

  /* Copy array (avoid modifying a), and initialize inverse */

  memcpy(_a, a, sizeof(double)*12);

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++)
      b[i][j] = 0.0;
  }
  for (i = 0; i < 3; i++)
    b[i][i] = 1.0;

  /* Forward elimination */

  for (i = 0; i < 2; i++) {

    /* Seek best pivot */

    k_pivot = i;
    abs_pivot = fabs(_a[i][i]);

    for (k = i+1; k < 3; k++) {
      abs_a_ki = fabs(_a[k][i]);
      if (abs_a_ki > abs_pivot) {
        abs_pivot = abs_a_ki;
        k_pivot = k;
      }
    }

    /* Swap if necessary */

    if (k_pivot != i) {
      for (j = 0; j < 4; j++) {
        _a_swap[j] = _a[i][j];
        _a[i][j] = _a[k_pivot][j];
        _a[k_pivot][j] = _a_swap[j];
        _b_swap[j] = b[i][j];
        b[i][j] = b[k_pivot][j];
        b[k_pivot][j] = _b_swap[j];
      }
    }

    if (abs_pivot < 1.0e-18)
      bft_error(__FILE__, __LINE__, 0,
                _("User-defined Transformation matrix is not invertible:\n"
                  "  [[%12.5f %12.5f %12.5f %12.5f]\n"
                  "   [%12.5f %12.5f %12.5f %12.5f]\n"
                  "   [%12.5f %12.5f %12.5f %12.5f]\n"
                  "   [%12.5f %12.5f %12.5f %12.5f]]\n"),
                a[0][0], a[0][1], a[0][2], a[0][3],
                a[1][0], a[1][1], a[1][2], a[1][3],
                a[2][0], a[2][1], a[2][2], a[2][3],
                0., 0., 0., 1.);

    /* Eliminate values */

    for (k = i+1; k < 3; k++) {
      factor = _a[k][i] / _a[i][i];
      _a[k][i] = 0.0;
      for (j = i+1; j < 4; j++) {
        _a[k][j] -= _a[i][j]*factor;
        b[k][j] -= b[i][j]*factor;
      }
    }

  }

  /* Normalize lines */

  for (i = 0; i < 3; i++) {
    factor = _a[i][i];
    for (j = 0; j < 4; j++) {
      _a[i][j] /= factor;
      b[i][j] /= factor;
    }
  }

  /* Eliminate last column */

  for (j = 3; j > 0; j--) {
    for (i = 0; i < j; i++) {
      b[i][j] -= _a[i][j];
      _a[i][j] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Combine transformation matrixes (c = a.b)
 *
 * parameters:
 *   a <-- first transformation matrix
 *   b <-- second transformation matrix
 *   c --> combined transformation matrix
 *---------------------------------------------------------------------------*/

static void
_combine_tr_matrixes(double  a[3][4],
                     double  b[3][4],
                     double  c[3][4])
{
  c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
  c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
  c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
  c[0][3] = a[0][0]*b[0][3] + a[0][1]*b[1][3] + a[0][2]*b[2][3] + a[0][3];

  c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
  c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
  c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
  c[1][3] = a[1][0]*b[0][3] + a[1][1]*b[1][3] + a[1][2]*b[2][3] + a[1][3];

  c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
  c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
  c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
  c[2][3] = a[2][0]*b[0][3] + a[2][1]*b[1][3] + a[2][2]*b[2][3] + a[2][3];
}

/*----------------------------------------------------------------------------
 * Transform a periodicity matrix in case of coordinates transformation
 *
 * The new transformation is the combination of the following:
 *  - switch to initial coordinates
 *  - apply initial periodicity
 *  - switch to new coordinates
 *
 * parameters:
 *   matrix       <-- matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   perio_matrix <-- matrix of the periodic matrix transformation
 *                    in homogeneous coord.
 *----------------------------------------------------------------------------*/

static void
_transform_perio_matrix(const double  matrix[],
                        double        perio_matrix[3][4])

{
  double  _m[3][4];
  double  _inv_m[3][4], _tmp_m[3][4];

  memcpy(_m, matrix, 12*sizeof(double));

  _inverse_transf_matrix(_m, _inv_m);

  _combine_tr_matrixes(perio_matrix, _inv_m, _tmp_m);
  _combine_tr_matrixes(_m, _tmp_m, perio_matrix);
}

/*----------------------------------------------------------------------------
 * Read pre-processor mesh data for a given mesh and finalize input.
 *
 * parameters:
 *   file_id <-- id of file handled by mesh builder
 *   mesh    <-> pointer to mesh structure
 *   mr      <-> pointer to mesh reader structure
 *   echo    <-- echo (verbosity) level
 *----------------------------------------------------------------------------*/

static void
_read_data(int                 file_id,
           cs_mesh_t          *mesh,
           cs_mesh_builder_t  *mb,
           _mesh_reader_t     *mr,
           long                echo)
{
  cs_int_t  perio_id, perio_type;
  cs_file_access_t  method;
  cs_io_sec_header_t  header;

  cs_real_t  perio_matrix[3][4];

  cs_int_t  perio_num = -1;
  bool       end_read = false;
  bool       data_read = false;
  cs_io_t  *pp_in = NULL;

  int gc_id_shift = mr->gc_id_shift[file_id];

  int n_perio_read = 0;
  cs_lnum_t n_cells = 0;
  cs_lnum_t n_faces = 0;
  cs_lnum_t n_vertices = 0;
  cs_lnum_t n_face_connect_size = 0;
  cs_gnum_t n_g_cells = 0;
  cs_gnum_t n_g_faces = 0;
  cs_gnum_t n_g_vertices = 0;
  cs_gnum_t n_g_face_connect_size = 0;

  cs_gnum_t face_vtx_range[2] = {0, 0};
  _mesh_file_info_t  *f = NULL;

  const char  *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                   "unexpected or of incorrect size.");

  f = mr->file_info + file_id;

#if defined(HAVE_MPI)
  {
    MPI_Info           hints;
    MPI_Comm           block_comm, comm;
    cs_file_get_default_access(CS_FILE_MODE_READ, &method, &hints);
    cs_file_get_default_comm(NULL, NULL, &block_comm, &comm);
    assert(comm == cs_glob_mpi_comm || comm == MPI_COMM_NULL);
    pp_in = cs_io_initialize(f->filename,
                             "Face-based mesh definition, R0",
                             CS_IO_MODE_READ,
                             method,
                             echo,
                             hints,
                             block_comm,
                             comm);
  }
#else
  {
    cs_file_get_default_access(CS_FILE_MODE_READ, &method);
    pp_in = cs_io_initialize(f->filename,
                             "Face-based mesh definition, R0",
                             CS_IO_MODE_READ,
                             method,
                             echo);
  }
#endif

  cs_io_set_offset(pp_in, f->offset);

  echo = cs_io_get_echo(pp_in);

  /* Loop on sections read */

  while (end_read == false) {

    /* Receive header and clean header name */

    cs_io_read_header(pp_in, &header);

    /* Process according to the header name */

    if (strncmp(header.sec_name, "EOF", CS_IO_NAME_LEN)
        == 0) {
      cs_io_finalize(&pp_in);
      pp_in = NULL;
    }

    if (strncmp(header.sec_name, "start_block:data",
                CS_IO_NAME_LEN) == 0) {

      if (data_read == false)
        data_read = true;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }
    else if (strncmp(header.sec_name, "end_block:data",
                     CS_IO_NAME_LEN) == 0) {

      if (data_read == true) {
        data_read = false;
        end_read = true;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }

    /* Read data from the pre-processor output file */

    else {

      cs_gnum_t gnum_range_cur[2];

      cs_lnum_t n_vals = 0;
      cs_lnum_t n_vals_cur = 0;
      cs_lnum_t val_offset_cur = 0;

      if (data_read != true)
        bft_error(__FILE__, __LINE__, 0,
                  _("Section of type <%s> on <%s>\n"
                    "unexpected."), header.sec_name, cs_io_get_name(pp_in));

      /* Face-cells connectivity */

      if (strncmp(header.sec_name, "face_cells", CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mb->n_g_faces,
                    mr->n_g_faces_read,
                    2,
                    0,
                    mb->face_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_faces,
                    &n_vals,
                    &n_vals_cur);

        n_faces = n_vals_cur / 2;
        val_offset_cur = mr->n_faces_read * 2;

        /* Allocate for first file read */
        if (mb->face_cells == NULL)
          BFT_MALLOC(mb->face_cells, n_vals, cs_gnum_t);

        /* Read data */
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mb->face_cells + val_offset_cur, pp_in);

        /* Shift referenced cell numbers in case of appended data */
        if (mr->n_g_cells_read > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mb->face_cells[val_offset_cur + ii] != 0)
              mb->face_cells[val_offset_cur + ii] += mr->n_g_cells_read;
          }
        }
      }

      /* Cell group class values */

      else if (strncmp(header.sec_name, "cell_group_class_id",
                     CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mesh->n_g_cells,
                    mr->n_g_cells_read,
                    1,
                    0,
                    mb->cell_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_cells,
                    &n_vals,
                    &n_vals_cur);

        n_cells = n_vals_cur;
        val_offset_cur = mr->n_cells_read;

        /* Allocate for first file read */
        if (mb->cell_gc_id == NULL)
          BFT_MALLOC(mb->cell_gc_id, n_vals, cs_int_t);

        /* Read data */
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mb->cell_gc_id + val_offset_cur, pp_in);

        /* Shift referenced numbers in case of appended data */
        if (gc_id_shift > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mb->cell_gc_id[val_offset_cur + ii] != 0)
              mb->cell_gc_id[val_offset_cur + ii] += gc_id_shift;
          }
        }
      }

      /* Face group class values */

      else if (strncmp(header.sec_name, "face_group_class_id",
                       CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mb->n_g_faces,
                    mr->n_g_faces_read,
                    1,
                    0,
                    mb->face_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_faces,
                    &n_vals,
                    &n_vals_cur);

        n_faces = n_vals_cur;
        val_offset_cur = mr->n_faces_read;

        /* Allocate for first file read */
        if (mb->face_gc_id == NULL)
          BFT_MALLOC(mb->face_gc_id, n_vals, cs_int_t);

        /* Read data */
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mb->face_gc_id + val_offset_cur, pp_in);

        /* Shift referenced numbers in case of appended data */
        if (gc_id_shift > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mb->face_gc_id[val_offset_cur + ii] != 0)
              mb->face_gc_id[val_offset_cur + ii] += gc_id_shift;
          }
        }
      }

      /* Face -> vertices connectivity */

      else if (strncmp(header.sec_name, "face_vertices_index",
                       CS_IO_NAME_LEN) == 0) {

        cs_lnum_t ii;
        cs_lnum_t idx_offset_shift = 0;
        cs_gnum_t idx_gnum_shift = 0;
        cs_gnum_t *_g_face_vertices_idx = NULL;

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mb->n_g_faces,
                    mr->n_g_faces_read,
                    1,
                    1,
                    mb->face_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_faces,
                    &n_vals,
                    &n_vals_cur);

        n_faces = n_vals_cur - 1;
        val_offset_cur = mr->n_faces_read;

        /* Allocate for first file read */
        if (mb->face_vertices_idx == NULL)
          BFT_MALLOC(mb->face_vertices_idx, n_vals, cs_lnum_t);

        if (val_offset_cur > 0)
          idx_offset_shift = mb->face_vertices_idx[val_offset_cur];

        /* Read data */
        cs_io_set_cs_gnum(&header, pp_in);
        BFT_MALLOC(_g_face_vertices_idx, n_vals_cur+1, cs_gnum_t);
        cs_io_read_index_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                               _g_face_vertices_idx, pp_in);

        /* save start and end values for next read */
        face_vtx_range[1] = _g_face_vertices_idx[n_vals_cur - 1];
        face_vtx_range[0] = _g_face_vertices_idx[0];

        /* Shift cell numbers in case of appended data */
        idx_gnum_shift = _g_face_vertices_idx[0];
        for (ii = 0; ii < n_vals_cur; ii++) {
          cs_gnum_t _face_vtx_idx = _g_face_vertices_idx[ii] - idx_gnum_shift;
          mb->face_vertices_idx[ii + val_offset_cur]
            = _face_vtx_idx + idx_offset_shift;
        }

        BFT_FREE(_g_face_vertices_idx);
      }

      else if (strncmp(header.sec_name, "face_vertices",
                       CS_IO_NAME_LEN) == 0) {

        n_vals = 0;
        n_g_face_connect_size = header.n_vals;

        if (  (cs_gnum_t)(header.n_vals) + mr->n_g_faces_connect_read
            > mb->n_g_face_connect_size)
          bft_error
            (__FILE__, __LINE__, 0,
             _("Section of type <%s> on <%s>\n"
               "has incorrect size (current: %llu, read: %llu, total: %llu."),
             header.sec_name, cs_io_get_name(pp_in),
             (unsigned long long)(header.n_vals),
             (unsigned long long)mr->n_g_faces_connect_read,
             (unsigned long long)mb->n_g_face_connect_size);

        else if (header.n_location_vals != 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("Section of type <%s> on <%s>\n"
                      "has incorrect number of values per location."),
                    header.sec_name, cs_io_get_name(pp_in));

        n_vals_cur = face_vtx_range[1] - face_vtx_range[0];
        val_offset_cur = mr->n_faces_connect_read;

        /* Reallocate for each read, as size of indexed array
           cannot be determined before reading the previous section
           (and is thus not yet known for future files). */
        BFT_REALLOC(mb->face_vertices,
                    mr->n_faces_connect_read + n_vals_cur,
                    cs_gnum_t);

        /* Read data */
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_block(&header, face_vtx_range[0], face_vtx_range[1],
                         mb->face_vertices + val_offset_cur, pp_in);

        /* Shift referenced vertex numbers in case of appended data */
        if (mr->n_g_vertices_read > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mb->face_vertices[val_offset_cur + ii] != 0)
              mb->face_vertices[val_offset_cur + ii] += mr->n_g_vertices_read;
          }
        }

        mr->n_faces_connect_read += n_vals_cur;
      }

      else if (strncmp(header.sec_name, "vertex_coords",
                       CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mesh->n_g_vertices,
                    mr->n_g_vertices_read,
                    3,
                    0,
                    mb->vertex_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_vertices,
                    &n_vals,
                    &n_vals_cur);

        n_vertices = n_vals_cur / 3;
        val_offset_cur = mr->n_vertices_read * 3;

        /* Allocate for first file read */
        if (mb->vertex_coords == NULL)
          BFT_MALLOC(mb->vertex_coords, n_vals, cs_real_t);

        /* Read data */
        cs_io_assert_cs_real(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mb->vertex_coords + val_offset_cur, pp_in);

        /* Transform coordinates if necessary */

        if (f->matrix != NULL) {
          cs_gnum_t range_size = gnum_range_cur[1] - gnum_range_cur[0];
          _transform_coords(range_size,
                            mb->vertex_coords + val_offset_cur,
                            f->matrix);
          mesh->modified = 1;
        }
      }

      /* Additional buffers for periodicity */

      else if (strncmp(header.sec_name, "periodicity_type_",
                       strlen("periodicity_type_")) == 0) {

        if (data_read != true || header.n_vals != 1)
          bft_error(__FILE__, __LINE__, 0,
                    _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
        else {
          perio_num = atoi(header.sec_name + strlen("periodicity_type_"));
          n_perio_read = CS_MAX(n_perio_read, perio_num);
          cs_io_read_global(&header, &perio_type, pp_in);
        }

      }
      else if (strncmp(header.sec_name, "periodicity_matrix_",
                       strlen("periodicity_matrix_")) == 0) {

        n_vals = 12; /* 3x4 */
        if (data_read != true || header.n_vals != n_vals)
          bft_error(__FILE__, __LINE__, 0,
                    _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
        else {
          assert(   perio_num
                 == atoi(header.sec_name + strlen("periodicity_matrix_")));
          cs_io_assert_cs_real(&header, pp_in);
          cs_io_read_global(&header, perio_matrix, pp_in);

          /* Add a periodicity to mesh->periodicities */

          if (f->matrix != NULL)
            _transform_perio_matrix(f->matrix, perio_matrix);

          _add_periodicity(mesh,
                           perio_type,
                           perio_num + mr->n_perio_read,
                           perio_matrix);

        }

      }
      else if (strncmp(header.sec_name, "periodicity_faces_",
                       strlen("periodicity_faces_")) == 0) {

        perio_id = atoi(header.sec_name
                        + strlen("periodicity_faces_")) - 1
                        + mr->n_perio_read;

        n_vals = mb->n_g_per_face_couples[perio_id] * 2;

        if (data_read != true || header.n_vals != n_vals)
          bft_error(__FILE__, __LINE__, 0,
                    _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
        else {

          if ((mb->per_face_bi[perio_id]).gnum_range[0] > 0)
            mb->n_per_face_couples[perio_id]
              = (  (mb->per_face_bi[perio_id]).gnum_range[1]
                 - (mb->per_face_bi[perio_id]).gnum_range[0]);
          else
            mb->n_per_face_couples[perio_id] = 0;

          cs_io_set_cs_gnum(&header, pp_in);
          n_vals = mb->n_per_face_couples[perio_id]*2;
          BFT_MALLOC(mb->per_face_couples[perio_id], n_vals, cs_gnum_t);
          assert(header.n_location_vals == 2);
          cs_io_read_block(&header,
                           (mb->per_face_bi[perio_id]).gnum_range[0],
                           (mb->per_face_bi[perio_id]).gnum_range[1],
                           mb->per_face_couples[perio_id],
                           pp_in);

          /* Shift referenced face numbers in case of appended data */
          if (mr->n_g_faces_read > 0) {
            cs_lnum_t ii;
            for (ii = 0; ii < n_vals; ii++) {
              if (mb->per_face_couples[perio_id][ii] != 0)
                mb->per_face_couples[perio_id][ii] += mr->n_g_faces_read;
            }
          }
        }
      }
    }

  } /* End of loop on sections */

  mr->n_perio_read += n_perio_read;
  mr->n_cells_read += n_cells;
  mr->n_faces_read += n_faces;
  mr->n_faces_connect_read += n_face_connect_size;
  mr->n_vertices_read += n_vertices;
  mr->n_g_cells_read += n_g_cells;
  mr->n_g_faces_read += n_g_faces;
  mr->n_g_faces_connect_read += n_g_face_connect_size;
  mr->n_g_vertices_read += n_g_vertices;

  /* Finalize pre-processor input */
  /*------------------------------*/

  f->offset = 0;
  cs_io_finalize(&pp_in);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Pass information relative to mesh metadata to the Fortran API
 *
 * Fortran Interface:
 *
 * subroutine ledevi(iperio, iperot)
 * *****************
 *
 * integer          iperio      : <-- : Periodicity indicator
 * integer          iperot      : <-- : Number of rotation periodicities
 *----------------------------------------------------------------------------*/

void
CS_PROCF(ledevi, LEDEVI)(cs_int_t   *iperio,
                         cs_int_t   *iperot)
{
  cs_mesh_t  *m = cs_glob_mesh;

  /* Initialize parameter values */

  if (m != NULL) {
    if (m->n_init_perio > 0)
      *iperio = 1;
    if (m->have_rotation_perio > 0)
      *iperot = 1;
  }
  else {
    int retval = cs_preprocessor_check_perio();
    if (retval > 0) {
      *iperio = 1;
      if (retval > 1)
        *iperot = 1;
    }
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define input mesh file to read.
 *
 * If this function is never called, the default file is read.
 * The first time this function is called,  this default is overriden by the
 * defined file, and all subsequent calls define additional meshes to read.
 *
 * parameters:
 *   file_name       <-- name of file to read
 *   n_group_renames <-- number of groups to rename
 *   group_rename    <-- old (group_rename[i*2]) to new (group_rename[i*2 + 1])
 *                       group names array (size: n_group_renames*2)
 *   transf_matrix   <-- coordinate transformation matrix (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_add_file(const char     *file_name,
                              size_t          n_group_renames,
                              const char    **group_rename,
                              const double    transf_matrix[3][4])
{
  size_t  i, l;
  size_t  data_size = 0;
  char  **_old_group_names = NULL, **_new_group_names = NULL;
  _mesh_file_info_t  *f = NULL;

  /* Compute data size */

  data_size = _align_size(strlen(file_name) + 1);

  if (transf_matrix != NULL)
    data_size += _align_size(12*sizeof(double));

  data_size += (_align_size(n_group_renames * sizeof(char *)) * 2);

  assert(n_group_renames == 0 || group_rename != NULL);

  for (i = 0; i < n_group_renames; i++) {
    data_size += _align_size(strlen(group_rename[i*2]) + 1);
    if (group_rename[i*2+1] != NULL)
      data_size += _align_size(strlen(group_rename[i*2+1]) + 1);
  }

  /* Allocate data (reallocate mesh file info array if necesary) */

  if (_n_max_mesh_files == 0) {
    _n_max_mesh_files = 1;
    BFT_MALLOC(_mesh_file_info, 1, _mesh_file_info_t);
  }

  if (_n_mesh_files + 1 > _n_max_mesh_files) {
    _n_max_mesh_files *= 2;
    BFT_REALLOC(_mesh_file_info, _n_max_mesh_files, _mesh_file_info_t);
  }

  f = _mesh_file_info + _n_mesh_files;
  _n_mesh_files += 1;

  /* Setup base structure fields */

  f->offset = 0;
  f->data_size = data_size;
  BFT_MALLOC(f->data, f->data_size, unsigned char);
  memset(f->data, 0, f->data_size);

  /* Setup data */

  data_size = 0;

  l = strlen(file_name) + 1;
  memcpy(f->data, file_name, l);
  f->filename = (const char *)(f->data);

  data_size = _align_size(l);

  if (transf_matrix != NULL) {
    l = 12*sizeof(double);
    memcpy(f->data + data_size, transf_matrix, l);
    f->matrix = (const double *)(f->data + data_size);
    data_size += _align_size(l);
  }
  else
    f->matrix = NULL;

  f->n_group_renames = n_group_renames;
  f->old_group_names = NULL;
  f->new_group_names = NULL;

  if (n_group_renames > 0) {

    _old_group_names = (char **)(f->data + data_size);
    f->old_group_names = (const char *const *)_old_group_names;
    data_size += _align_size(n_group_renames * sizeof(char *));

    _new_group_names = (char **)(f->data + data_size);
    f->new_group_names = (const char *const *)_new_group_names;
    data_size += _align_size(n_group_renames * sizeof(char *));

  }

  for (i = 0; i < n_group_renames; i++) {
    l = strlen(group_rename[i*2]) + 1;
    _old_group_names[i] = (char *)(f->data + data_size);
    memcpy(_old_group_names[i], group_rename[i*2], l);
    data_size += _align_size(l);
    if (group_rename[i*2+1] != NULL) {
      l = strlen(group_rename[i*2+1]) + 1;
      _new_group_names[i] = (char *)(f->data + data_size);
      memcpy(_new_group_names[i], group_rename[i*2+1], l);
      data_size += _align_size(l);
    }
    else
      _new_group_names[i] = NULL;
  }
}

/*----------------------------------------------------------------------------
 * Check for periodicity information in mesh meta-data.
 *
 * returns:
 *   0 if no periodicity is present in mesh input,
 *   1 for translation periodicity only,
 *   2 for rotation or mixed periodicity
 *----------------------------------------------------------------------------*/

int
cs_preprocessor_check_perio(void)
{
  int retval;
  _mesh_reader_t *mr = NULL;

  int perio_flag = 0;

  /* Initialize reading of Preprocessor output */

  _set_default_input_if_needed();

  mr = _mesh_reader_create(&_n_mesh_files,
                           &_mesh_file_info);

  _n_max_mesh_files = 0;

  for (int i = 0; i < _n_mesh_files; i++) {
    retval = _read_perio_info((_mesh_file_info + i)->filename);
    perio_flag = CS_MAX(retval, perio_flag);
  }

  _mesh_reader_destroy(&mr);

  /* Return values */

  cs_parall_max(1, CS_INT_TYPE, &perio_flag);

  return perio_flag;
}

/*----------------------------------------------------------------------------
 * Read mesh meta-data.
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_headers(cs_mesh_t          *mesh,
                                  cs_mesh_builder_t  *mesh_builder)
{
  int file_id;

  _mesh_reader_t *mr = NULL;

  /* Initialize reading of Preprocessor output */

  _set_default_input_if_needed();

  _cs_glob_mesh_reader = _mesh_reader_create(&_n_mesh_files,
                                             &_mesh_file_info);

  _n_max_mesh_files = 0;

  mr = _cs_glob_mesh_reader;

  for (file_id = 0; file_id < mr->n_files; file_id++)
    _read_dimensions(mesh, mesh_builder, mr, file_id);

  /* Return values */

  assert(mesh->dim == 3);

  mesh->n_domains = cs_glob_n_ranks;
  mesh->domain_num = cs_glob_rank_id + 1;

  /* Update data in cs_mesh_t structure in serial mode */

  if (cs_glob_n_ranks == 1) {
    mesh->n_cells = mesh->n_g_cells;
    mesh->n_cells_with_ghosts = mesh->n_cells;
    mesh->domain_num = 1;
  }
  else
    mesh->domain_num = cs_glob_rank_id + 1;

  /* Clean group names */

  cs_mesh_group_clean(mesh);
}

/*----------------------------------------------------------------------------
 * Read pre-processor mesh data and finalize input.
 *
 * At this stage, ghost cells are not generated yet, so the interior
 * face connectivity is not complete near parallel domain or periodic
 * boundaries. Also, isolated faces, if present, are considered to be
 * boundary faces, as they may participate in future mesh joining
 * operations. Their matching cell number will be set to -1.
 * Remaining isolated faces should be removed before completing
 * the mesh structure (for example, just before building ghost cells).
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_mesh(cs_mesh_t          *mesh,
                               cs_mesh_builder_t  *mesh_builder)
{
  int file_id;

  cs_partition_stage_t partition_stage
    =     (cs_partition_get_preprocess())
       ?  CS_PARTITION_FOR_PREPROCESS : CS_PARTITION_MAIN;

  long  echo = CS_IO_ECHO_OPEN_CLOSE;
  _mesh_reader_t  *mr = _cs_glob_mesh_reader;

  bool pre_partitioned = false;

  /* Check for existing partitioning and cell block info (set by
     cs_mesh_to_builder_partition and valid if the global number of
     cells has not changed), in which case the existing
     partitioning may be used */

  if (mesh_builder->have_cell_rank) {

    cs_block_dist_info_t cell_bi_ref;
    memcpy(&cell_bi_ref,
           &(mesh_builder->cell_bi),
           sizeof(cs_block_dist_info_t));
    _set_block_ranges(mesh, mesh_builder);
    cs_gnum_t n_g_cells_ref = 0;
    if (cell_bi_ref.gnum_range[1] > cell_bi_ref.gnum_range[0])
      n_g_cells_ref = cell_bi_ref.gnum_range[1] - cell_bi_ref.gnum_range[0];
    cs_parall_counter(&n_g_cells_ref, 1);

    _set_block_ranges(mesh, mesh_builder);

    if (n_g_cells_ref == mesh->n_g_cells) {
      memcpy(&(mesh_builder->cell_bi),
             &cell_bi_ref,
             sizeof(cs_block_dist_info_t));
      pre_partitioned = true;
    }
    else {
      mesh_builder->have_cell_rank = false;
      BFT_FREE(mesh_builder->cell_rank);
    }

  }
  else
    _set_block_ranges(mesh, mesh_builder);

  for (file_id = 0; file_id < mr->n_files; file_id++)
    _read_data(file_id, mesh, mesh_builder, mr, echo);

  if (mr->n_files > 1)
    mesh->modified = 1;

  /* Partition data */

  if (! pre_partitioned)
    cs_partition(mesh, mesh_builder, partition_stage);

  bft_printf("\n");

  /* Now send data to the correct rank */
  /*-----------------------------------*/

  cs_mesh_from_builder(mesh, mesh_builder);

  /* Free temporary memory */

  _mesh_reader_destroy(&mr);
  _cs_glob_mesh_reader = mr;

  /* Remove family duplicates */

  cs_mesh_clean_families(mesh);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
