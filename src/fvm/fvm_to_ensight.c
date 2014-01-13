/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to EnSight Gold files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_to_ensight_case.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_ensight.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * EnSight Gold writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;               /* Writer name */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

  _Bool        text_mode;          /* true if using text output */
  _Bool        swap_endian;        /* true if binary file endianness must
                                      be changed */

  _Bool        discard_polygons;   /* Option to discard polygonal elements */
  _Bool        discard_polyhedra;  /* Option to discard polyhedral elements */

  _Bool        divide_polygons;    /* Option to tesselate polygonal elements */
  _Bool        divide_polyhedra;   /* Option to tesselate polyhedral elements */

  fvm_to_ensight_case_t  *case_info;  /* Associated case structure */

#if defined(HAVE_MPI)
  int          min_rank_step;      /* Minimum rank step */
  int          min_block_size;     /* Minimum block buffer size */
  MPI_Comm     block_comm;         /* Associated MPI block communicator */
  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

} fvm_to_ensight_writer_t;

/*----------------------------------------------------------------------------
 * Indirect file structure to handle both text and binary files
 *----------------------------------------------------------------------------*/

typedef struct {

  FILE        *tf;                 /* Text file handing structure */
  cs_file_t   *bf;                 /* Binary file handling structure */

} _ensight_file_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char  *_ensight_type_name[FVM_N_ELEMENT_TYPES] = {"bar2",
                                                               "tria3",
                                                               "quad4",
                                                               "nsided",
                                                               "tetra4",
                                                               "pyramid5",
                                                               "penta6",
                                                               "hexa8",
                                                               "nfaced"};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Open an EnSight Gold geometry or variable file
 *
 * parameters:
 *   this_writer <-- pointer to Ensight Gold writer structure.
 *   filename    <-- name of file to open.
 *   apend       <-- if true, append to file instead of overwriting
 *----------------------------------------------------------------------------*/

static _ensight_file_t
_open_ensight_file(const fvm_to_ensight_writer_t  *this_writer,
                   const char                     *filename,
                   _Bool                           append)
{
  _ensight_file_t f = {NULL, NULL};

  if (this_writer->text_mode == true) {
    if (this_writer->rank == 0) {
      if (append)
        f.tf = fopen(filename, "a");
      else
        f.tf = fopen(filename, "w");
      if (f.tf == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error opening file \"%s\":\n\n"
                    "  %s"), filename, strerror(errno));
    }
  }
  else {

    cs_file_mode_t mode = append ? CS_FILE_MODE_APPEND : CS_FILE_MODE_WRITE;
    cs_file_access_t method;

#if defined(HAVE_MPI)

    MPI_Info hints;
    cs_file_get_default_access(CS_FILE_MODE_WRITE, &method, &hints);
    f.bf = cs_file_open(filename,
                        mode,
                        method,
                        hints,
                        this_writer->block_comm,
                        this_writer->comm);

#else

    cs_file_get_default_access(CS_FILE_MODE_WRITE, &method);
    f.bf = cs_file_open(filename, mode, method);

#endif

    if (this_writer->swap_endian == true)
      cs_file_set_swap_endian(f.bf, 1);
  }

  return f;
}

/*----------------------------------------------------------------------------
 * close an EnSight Gold geometry or variable file
 *
 * parameters:
 *   f <-- pointer to file handler structure.
 *----------------------------------------------------------------------------*/

static void
_free_ensight_file(_ensight_file_t  *f)
{
  if (f->tf != NULL) {
    if (fclose(f->tf) != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error closing EnSight output file (text mode):\n\n"
                  "  %s"), strerror(errno));
    f->tf = NULL;
  }

  else if (f->bf != NULL)
    f->bf = cs_file_free(f->bf);
}

/*----------------------------------------------------------------------------
 * Count number of extra vertices when tesselations are present
 *
 * parameters:
 *   mesh               <-- pointer to nodal mesh structure
 *   divide_polyhedra   <-- true if polyhedra are tesselated
 *   n_extra_vertices_g --> global number of extra vertices (optional)
 *   n_extra_vertices   --> local number of extra vertices (optional)
 *----------------------------------------------------------------------------*/

static void
_count_extra_vertices(const fvm_nodal_t  *mesh,
                      _Bool               divide_polyhedra,
                      cs_gnum_t          *n_extra_vertices_g,
                      cs_lnum_t          *n_extra_vertices)
{
  int  i;

  const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (n_extra_vertices_g != NULL)
    *n_extra_vertices_g = 0;
  if (n_extra_vertices != NULL)
    *n_extra_vertices   = 0;

  if (divide_polyhedra) {

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == export_dim
          && section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        if (n_extra_vertices_g != NULL)
          *n_extra_vertices_g
            += fvm_tesselation_n_g_vertices_add(section->tesselation);

        if (n_extra_vertices != NULL)
          *n_extra_vertices
            += fvm_tesselation_n_vertices_add(section->tesselation);

      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Write string to a text or C binary EnSight Gold file
 *
 * parameters:
 *   f <-- file to write to
 *   s <-- string to write
 *----------------------------------------------------------------------------*/

static void
_write_string(_ensight_file_t   f,
              const char       *s)
{
  size_t  i;
  char  buf[82];

  if (f.tf != NULL) {
    strncpy(buf, s, 80);
    buf[80] = '\0';
    fprintf(f.tf, "%s\n", buf);
  }

  else if (f.bf != NULL) {
    strncpy(buf, s, 80);
    buf[80] = '\0';
    for (i = strlen(buf); i < 80; i++)
      buf[i] = ' ';
    cs_file_write_global(f.bf, buf, 1, 80);
  }
}

/*----------------------------------------------------------------------------
 * Write integer to a text or C binary EnSight Gold file
 *
 * parameters:
 *   f <-- file to write to
 *   n <-- integer value to write
 *----------------------------------------------------------------------------*/

inline static void
_write_int(_ensight_file_t  f,
           int32_t          n)
{
  if (f.tf != NULL)
    fprintf(f.tf, "%10d\n", (int)n);

  else if (f.bf != NULL) {
    int _n = n;
    cs_file_write_global(f.bf, &_n, sizeof(int32_t), 1);
  }
}

/*----------------------------------------------------------------------------
 * Write headers to an EnSight Gold geometry file.
 *
 * parameters:
 *   this_writer <-- pointer to Ensight Gold writer structure.
 *   f           <-- pointer to file to initialize.
 *----------------------------------------------------------------------------*/

static void
_write_geom_headers(fvm_to_ensight_writer_t  *this_writer,
                    _ensight_file_t           f)
{
  if (f.bf != NULL)
    _write_string(f, "C Binary");

  /* 1st description line */
  {
    char buf[81] = "";
    if (this_writer->name != NULL)
      strncpy(buf, this_writer->name, 80);
    buf[80] = '\0';
    _write_string(f, buf);
  }
  /* 2nd description line */
  _write_string(f, "Output by Code_Saturne version "VERSION);
  _write_string(f, "node id assign");
  _write_string(f, "element id assign");
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write block of a vector of floats to an EnSight Gold file.
 *
 * This variant is called in parallel mode, where the values are already
 * in a temporary block buffer and may be discarded.
 *
 * parameters:
 *   num_start <-- global number of first element for this block
 *   num_end   <-- global number of past the last element for this block
 *   values    <-- pointer to values block array
 *   comm      <-- associated MPI communicator
 *   f         <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_block_floats_g(cs_gnum_t         num_start,
                      cs_gnum_t         num_end,
                      float             values[],
                      MPI_Comm          comm,
                      _ensight_file_t   f)
{
  size_t  i;
  cs_gnum_t   j;

  /* In Binary mode, all ranks have a file structure,
     we may use use a collective call */

  if (f.bf != NULL)
    cs_file_write_block_buffer(f.bf,
                               values,
                               sizeof(float),
                               1,
                               num_start,
                               num_end);

  /* If all ranks do not have a binary file structure pointer, then
     we are using a text file, open only on rank 0 */

  else {

    float *_values = NULL;
    cs_file_serializer_t *s
      = cs_file_serializer_create(sizeof(float),
                                  1,
                                  num_start,
                                  num_end,
                                  0,
                                  values,
                                  comm);

    do {

      cs_gnum_t range[2] = {num_start, num_end};

      _values = cs_file_serializer_advance(s, range);

      if (_values != NULL) { /* only possible on rank 0 */
        assert(f.tf != NULL);
        for (i = 0, j = range[0]; j < range[1]; i++, j++)
          fprintf(f.tf, "%12.5e\n", _values[i]);
      }

    } while (_values != NULL);

    cs_file_serializer_destroy(&s);
  }

}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write block of a vector of floats to an EnSight Gold file
 *
 * parameters:
 *   n_values <-- number of values to write
 *   num_end  <-- global number of past the last element for this block
 *   values   <-- pointer to values block array
 *   f        <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_block_floats_l(size_t           n_values,
                      const float      values[],
                      _ensight_file_t  f)
{
  size_t  i;

  if (f.tf != NULL) {
    for (i = 0; i < n_values; i++)
      fprintf(f.tf, "%12.5e\n", values[i]);
  }
  else if (f.bf != NULL)
    cs_file_write_global(f.bf, values, sizeof(float), n_values);
}

/*----------------------------------------------------------------------------
 * Return extra vertex coordinates when tesselations are present
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   n_extra_vertices <-- number of extra vertices
 *
 * returns:
 *   array containing all extra vertex coordinates
 *----------------------------------------------------------------------------*/

static cs_coord_t *
_extra_vertex_coords(const fvm_nodal_t  *mesh,
                     cs_lnum_t           n_extra_vertices)
{
  int  i;
  cs_lnum_t   n_extra_vertices_section;

  size_t  coord_shift = 0;
  cs_coord_t  *coords = NULL;

  if (n_extra_vertices > 0) { /* This implies divide_polyhedra = true */

    BFT_MALLOC(coords, n_extra_vertices * 3, cs_coord_t);

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      if (   section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        n_extra_vertices_section
          = fvm_tesselation_n_vertices_add(section->tesselation);

        if (n_extra_vertices_section > 0) {

          fvm_tesselation_vertex_coords(section->tesselation,
                                        coords + coord_shift);

          coord_shift += n_extra_vertices_section * 3;

        }

      }
    }
  }

  return coords;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get extra vertex global numbers when tesselations are present
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   n_extra_vertices <-- number of extra vertices
 *   vtx_gnum         --> extra vertex global numbers (size: n_extra_vertices)
 *----------------------------------------------------------------------------*/

static void
_extra_vertex_get_gnum(const fvm_nodal_t  *mesh,
                       cs_lnum_t           n_extra_vertices,
                       cs_gnum_t           vtx_gnum[])
{
  int  i;
  cs_lnum_t   j = 0;
  cs_lnum_t   start_id = 0;
  cs_gnum_t   gnum_shift
    = fvm_io_num_get_global_count(mesh->global_vertex_num);

  if (n_extra_vertices > 0) { /* Implies divide_polyhedra */

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      if (   section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        cs_lnum_t n_extra_vertices_section
          = fvm_tesselation_n_vertices_add(section->tesselation);

        if (n_extra_vertices_section > 0) {

          const fvm_io_num_t *extra_vertex_num
            = fvm_tesselation_global_vertex_num(section->tesselation);
          const cs_gnum_t *extra_gnum
            = fvm_io_num_get_global_num(extra_vertex_num);

          for (j = 0; j < n_extra_vertices_section; j++)
            vtx_gnum[start_id + j] = extra_gnum[j] + gnum_shift;

          start_id += n_extra_vertices_section;

        }

        gnum_shift
          = fvm_tesselation_n_g_vertices_add(section->tesselation);

      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Build block info and part to block distribution helper for vertices.
 *
 * parameters:
 *   w    <-- pointer to writer structure
 *   mesh <-- pointer to nodal mesh structure
 *   bi   --> block information structure
 *   d    --> part to bloc distributor
 *----------------------------------------------------------------------------*/

static void
_vertex_part_to_block_create(const fvm_to_ensight_writer_t   *w,
                             const fvm_nodal_t               *mesh,
                             cs_block_dist_info_t            *bi,
                             cs_part_to_block_t             **d)
{
  cs_gnum_t    n_g_extra_vertices = 0, n_g_vertices_tot = 0;
  cs_lnum_t    n_extra_vertices = 0, n_vertices_tot = 0;

  cs_gnum_t   *_g_num = NULL;

  cs_block_dist_info_t  _bi;
  cs_part_to_block_t   *_d;

  size_t min_block_size = w->min_block_size / sizeof(float);

  const cs_lnum_t   n_vertices
    = fvm_io_num_get_local_count(mesh->global_vertex_num);
  cs_gnum_t   n_g_vertices
    = fvm_io_num_get_global_count(mesh->global_vertex_num);
  const cs_gnum_t   *g_num
    = fvm_io_num_get_global_num(mesh->global_vertex_num);

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(mesh,
                        w->divide_polyhedra,
                        &n_g_extra_vertices,
                        &n_extra_vertices);

  n_vertices_tot = n_vertices + n_extra_vertices;
  n_g_vertices_tot = n_g_vertices + n_g_extra_vertices;

  _bi = cs_block_dist_compute_sizes(w->rank,
                                    w->n_ranks,
                                    w->min_rank_step,
                                    min_block_size,
                                    n_g_vertices_tot);

  /* Global vertex numbers */


  if (n_extra_vertices > 0) {

    BFT_MALLOC(_g_num, n_vertices_tot, cs_gnum_t);

    memcpy(_g_num, g_num, n_vertices*sizeof(cs_gnum_t));
    _extra_vertex_get_gnum(mesh, n_extra_vertices, _g_num + n_vertices);

    g_num = _g_num;

  }

  /* Build distribution structures */

  _d = cs_part_to_block_create_by_gnum(w->comm, _bi, n_vertices_tot, g_num);

  if (n_extra_vertices > 0)
    cs_part_to_block_transfer_gnum(_d, _g_num);

  /* Return initialized structures */

  if (bi != NULL)
    *bi = _bi;

  if (d != NULL)
    *d = _d;
}

/*----------------------------------------------------------------------------
 * Write vertex coordinates to an EnSight Gold file in parallel mode
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   mesh        <-- pointer to nodal mesh structure
 *   f           <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_g(const fvm_to_ensight_writer_t  *this_writer,
                        const fvm_nodal_t              *mesh,
                        _ensight_file_t                 f)
{
  cs_lnum_t   i, j;
  size_t stride;
  cs_block_dist_info_t  bi;

  cs_gnum_t    n_g_extra_vertices = 0, n_g_vertices_tot = 0;
  cs_lnum_t    n_extra_vertices = 0, n_vertices_tot = 0;
  cs_coord_t  *extra_vertex_coords = NULL;
  float        *part_coords = NULL, *block_coords = NULL;

  cs_part_to_block_t   *d = NULL;
  size_t                block_buf_size = 0;

  const double      *vertex_coords = mesh->vertex_coords;
  const cs_lnum_t   *parent_vertex_num = mesh->parent_vertex_num;
  const cs_lnum_t   n_vertices
      = fvm_io_num_get_local_count(mesh->global_vertex_num);
  cs_gnum_t   n_g_vertices
      = fvm_io_num_get_global_count(mesh->global_vertex_num);

  /* Initialize distribution info */

  _vertex_part_to_block_create(this_writer,
                               mesh,
                               &bi,
                               &d);

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(mesh,
                        this_writer->divide_polyhedra,
                        &n_g_extra_vertices,
                        &n_extra_vertices);

  n_vertices_tot = n_vertices + n_extra_vertices;
  n_g_vertices_tot = n_g_vertices + n_g_extra_vertices;

  extra_vertex_coords = _extra_vertex_coords(mesh,
                                             n_extra_vertices);

  /* Build arrays */

  block_buf_size = (bi.gnum_range[1] - bi.gnum_range[0]);
  BFT_MALLOC(block_coords, block_buf_size, float);
  BFT_MALLOC(part_coords, n_vertices_tot, float);

  /* Vertex coordinates */
  /*--------------------*/

  stride = (size_t)(mesh->dim);

  _write_string(f, "coordinates");
  _write_int(f, (int)(n_g_vertices_tot));

  /* Loop on dimension (de-interlace coordinates, always 3D for EnSight) */

  for (j = 0; j < 3; j++) {

    if (j < mesh->dim) {

      if (parent_vertex_num != NULL) {
        for (i = 0; i < n_vertices; i++)
          part_coords[i] = vertex_coords[(parent_vertex_num[i]-1)*stride + j];
      }
      else {
        for (i = 0; i < n_vertices; i++)
          part_coords[i] = vertex_coords[i*stride + j];
      }

      for (i = 0; i < n_extra_vertices; i++)
        part_coords[n_vertices + i] = extra_vertex_coords[(i*stride) + j];
    }
    else {
      for (i = 0; i < n_vertices_tot; i++)
        part_coords[i] = 0.0;
    }

    cs_part_to_block_copy_array(d,
                                CS_FLOAT,
                                1,
                                part_coords,
                                block_coords);

    _write_block_floats_g(bi.gnum_range[0],
                          bi.gnum_range[1],
                          block_coords,
                          this_writer->comm,
                          f);

  } /* end of loop on spatial dimension */

  cs_part_to_block_destroy(&d);

  BFT_FREE(block_coords);
  BFT_FREE(part_coords);
  if (extra_vertex_coords != NULL)
    BFT_FREE(extra_vertex_coords);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write vertex coordinates to an EnSight Gold file in serial mode
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   mesh        <-- pointer to nodal mesh structure
 *   f           <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_l(const fvm_to_ensight_writer_t  *this_writer,
                        const fvm_nodal_t              *mesh,
                        _ensight_file_t                 f)
{
  cs_lnum_t    i, j;
  cs_lnum_t    n_extra_vertices = 0;
  cs_coord_t  *extra_vertex_coords = NULL;
  float        *coords_tmp = NULL;

  const cs_lnum_t    n_vertices = mesh->n_vertices;
  const double      *vertex_coords = mesh->vertex_coords;
  const cs_lnum_t   *parent_vertex_num = mesh->parent_vertex_num;

  const size_t  stride = (size_t)(mesh->dim);

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(mesh,
                        this_writer->divide_polyhedra,
                        NULL,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(mesh,
                                             n_extra_vertices);

  /* Vertex coordinates */
  /*--------------------*/

  BFT_MALLOC(coords_tmp, CS_MAX(n_vertices, n_extra_vertices), float);

  _write_string(f, "coordinates");
  _write_int(f, n_vertices + n_extra_vertices);

  /* Loop on dimension (de-interlace coordinates, always 3D for EnSight) */

  for (j = 0; j < 3; j++) {

    /* First, handle regular vertices */

    if (j < mesh->dim) {
      if (parent_vertex_num != NULL) {
        for (i = 0; i < n_vertices; i++) {
          assert(parent_vertex_num[i] != 0);
          coords_tmp[i]
            = (float)(vertex_coords[(parent_vertex_num[i]-1)*stride + j]);
        }
      }
      else {
        for (i = 0; i < n_vertices; i++)
          coords_tmp[i] = (float)(vertex_coords[i*stride + j]);
      }
    }
    else {
      for (i = 0; i < (n_vertices); i++)
        coords_tmp[i] = 0.0;
    }

    _write_block_floats_l(n_vertices,
                          coords_tmp,
                          f);

    /* Handle extra vertices (only occur with polyhedra tesselations in 3d) */

    for (i = 0; i < n_extra_vertices; i++)
      coords_tmp[i] = (float)(extra_vertex_coords[i*stride + j]);

    if (n_extra_vertices > 0)
      _write_block_floats_l(n_extra_vertices,
                            coords_tmp,
                            f);

  } /* end of loop on mesh dimension */

  BFT_FREE(coords_tmp);

  if (extra_vertex_coords != NULL)
    BFT_FREE(extra_vertex_coords);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided connectivity block to an EnSight Gold file in text mode
 *
 * parameters:
 *   stride  <-- number of vertices per element type
 *   n_elems <-- number of elements
 *   connect <-- connectivity array
 *   tf      <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_connect_block_gt(int             stride,
                        cs_lnum_t       n_elems,
                        const int32_t   connect[],
                        FILE           *tf)
{
  cs_lnum_t   i;

  assert(tf != NULL);

  switch(stride) {

  case 1: /* length */
    for (i = 0; i < n_elems; i++)
      fprintf(tf, "%10d\n",
              (int)connect[i]);
    break;

  case 2: /* edge */
    for (i = 0; i < n_elems; i++)
      fprintf(tf, "%10d%10d\n",
              (int)connect[i*2],
              (int)connect[i*2+1]);
    break;

  case 3: /* FVM_FACE_TRIA */
    for (i = 0; i < n_elems; i++)
      fprintf(tf, "%10d%10d%10d\n",
              (int)connect[i*3],
              (int)connect[i*3+1],
              (int)connect[i*3+2]);
    break;

  case 4: /* FVM_FACE_QUAD or FVM_CELL_TETRA */
    for (i = 0; i < n_elems; i++)
      fprintf(tf, "%10d%10d%10d%10d\n",
              (int)connect[i*4],
              (int)connect[i*4+1],
              (int)connect[i*4+2],
              (int)connect[i*4+3]);
    break;

  case 5: /* FVM_CELL_PYRAM */
    for (i = 0; i < n_elems; i++)
      fprintf(tf, "%10d%10d%10d%10d%10d\n",
              (int)connect[i*5],
              (int)connect[i*5+1],
              (int)connect[i*5+2],
              (int)connect[i*5+3],
              (int)connect[i*5+4]);
    break;

  case 6: /* FVM_CELL_PRISM */
    for (i = 0; i < n_elems; i++)
      fprintf(tf, "%10d%10d%10d%10d%10d%10d\n",
              (int)connect[i*6],
              (int)connect[i*6+1],
              (int)connect[i*6+2],
              (int)connect[i*6+3],
              (int)connect[i*6+4],
              (int)connect[i*6+5]);
    break;

  case 8: /* FVM_CELL_HEXA */
    for (i = 0; i < n_elems; i++)
      fprintf(tf,
              "%10d%10d%10d%10d%10d%10d%10d%10d\n",
              (int)connect[i*8],
              (int)connect[i*8+1],
              (int)connect[i*8+2],
              (int)connect[i*8+3],
              (int)connect[i*8+4],
              (int)connect[i*8+5],
              (int)connect[i*8+6],
              (int)connect[i*8+7]);
    break;

  default:
    assert(0);
  }
}

/*----------------------------------------------------------------------------
 * Write strided global connectivity block to an EnSight Gold file
 *
 * parameters:
 *   stride        <-- number of vertices per element type
 *   num_start     <-- global number of first element for this block
 *   num_end       <-- global number of past last element for this block
 *   block_connect <-> global connectivity block array
 *   comm          <-- associated MPI communicator
 *   f             <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_write_block_connect_g(int              stride,
                       cs_gnum_t        num_start,
                       cs_gnum_t        num_end,
                       int32_t          block_connect[],
                       MPI_Comm         comm,
                       _ensight_file_t  f)
{
  /* In Binary mode, all ranks have a file structure,
     we may use use a collective call */

  if (f.bf != NULL)
    cs_file_write_block(f.bf,
                        block_connect,
                        sizeof(int32_t),
                        stride,
                        num_start,
                        num_end);

  /* If all ranks do not have a binary file structure pointer, then
     we are using a text file, open only on rank 0 */

  else {

    int32_t *_block_connect = NULL;

    cs_file_serializer_t *s = cs_file_serializer_create(sizeof(int32_t),
                                                        stride,
                                                        num_start,
                                                        num_end,
                                                        0,
                                                        block_connect,
                                                        comm);

    do {
      cs_gnum_t range[2] = {num_start, num_end};

      _block_connect = cs_file_serializer_advance(s, range);

      if (_block_connect != NULL) /* only possible on rank 0 */
        _write_connect_block_gt(stride,
                                (range[1] - range[0]),
                                _block_connect,
                                f.tf);

    } while (_block_connect != NULL);

    cs_file_serializer_destroy(&s);
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write strided local connectivity to an EnSight Gold file
 *
 * parameters:
 *   stride      <-- number of vertices per element type
 *   n_elems     <-- number of elements
 *   connect     <-- connectivity array
 *   f           <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_connect_l(int                stride,
                 cs_lnum_t          n_elems,
                 const cs_lnum_t    connect[],
                 _ensight_file_t    f)
{
  cs_lnum_t   i;

  if (f.tf != NULL) { /* Text mode */

    switch(stride) {

    case 2: /* edge */
      for (i = 0; i < n_elems; i++)
        fprintf(f.tf, "%10d%10d\n",
                (int)connect[i*2],
                (int)connect[i*2+1]);
      break;

    case 3: /* FVM_FACE_TRIA */
      for (i = 0; i < n_elems; i++)
        fprintf(f.tf, "%10d%10d%10d\n",
                (int)connect[i*3],
                (int)connect[i*3+1],
                (int)connect[i*3+2]);
      break;

    case 4: /* FVM_FACE_QUAD or FVM_CELL_TETRA */
      for (i = 0; i < n_elems; i++)
        fprintf(f.tf, "%10d%10d%10d%10d\n",
                (int)connect[i*4],
                (int)connect[i*4+1],
                (int)connect[i*4+2],
                (int)connect[i*4+3]);
      break;

    case 5: /* FVM_CELL_PYRAM */
      for (i = 0; i < n_elems; i++)
        fprintf(f.tf, "%10d%10d%10d%10d%10d\n",
                (int)connect[i*5],
                (int)connect[i*5+1],
                (int)connect[i*5+2],
                (int)connect[i*5+3],
                (int)connect[i*5+4]);
      break;

    case 6: /* FVM_CELL_PRISM */
      for (i = 0; i < n_elems; i++)
        fprintf(f.tf, "%10d%10d%10d%10d%10d%10d\n",
                (int)connect[i*6],
                (int)connect[i*6+1],
                (int)connect[i*6+2],
                (int)connect[i*6+3],
                (int)connect[i*6+4],
                (int)connect[i*6+5]);
      break;

    case 8: /* FVM_CELL_HEXA */
      for (i = 0; i < n_elems; i++)
        fprintf(f.tf,
                "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                (int)connect[i*8],
                (int)connect[i*8+1],
                (int)connect[i*8+2],
                (int)connect[i*8+3],
                (int)connect[i*8+4],
                (int)connect[i*8+5],
                (int)connect[i*8+6],
                (int)connect[i*8+7]);
      break;

    default:
      assert(0);
    }

  }
  else if (f.bf != NULL) { /* Binary mode */

    size_t  j;
    size_t  k = 0;
    int32_t  *buffer = NULL;
    const size_t  n_values = n_elems*stride;
    const size_t  buffer_size = n_values >  64 ? (n_values / 8) : n_values;

    BFT_MALLOC(buffer, buffer_size, int32_t);

    while (k < n_values) {
      for (j = 0; j < buffer_size && k < n_values; j++)
        buffer[j] = (int)(connect[k++]);
      cs_file_write_global(f.bf, buffer, sizeof(int32_t), j);
    }

    BFT_FREE(buffer);
  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write "trivial" point elements to an EnSight Gold file in parallel mode
 *
 * parameters:
 *   w    <-- pointer to writer structure
 *   mesh <-- pointer to nodal mesh structure
 *   f    <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_export_point_elements_g(const fvm_to_ensight_writer_t  *w,
                         const fvm_nodal_t              *mesh,
                         _ensight_file_t                 f)
{
  const cs_gnum_t   n_g_vertices
    = fvm_io_num_get_global_count(mesh->global_vertex_num);

  _write_string(f, "point");
  _write_int(f, (int)n_g_vertices);

  if (f.tf != NULL) { /* Text mode, rank 0 only */

    cs_gnum_t   i;
    int32_t  j = 1;

    for (i = 0; i < n_g_vertices; i++)
      fprintf(f.tf, "%10d\n", j++);

  }
  else if (f.bf != NULL) { /* Binary mode */

    cs_lnum_t i;
    cs_gnum_t j;
    cs_block_dist_info_t  bi;

    size_t min_block_size = w->min_block_size / sizeof(float);
    int32_t  *connect = NULL;

    bi = cs_block_dist_compute_sizes(w->rank,
                                     w->n_ranks,
                                     w->min_rank_step,
                                     min_block_size,
                                     n_g_vertices);

    BFT_MALLOC(connect, bi.gnum_range[1] - bi.gnum_range[0], int32_t);

    for (i = 0, j = bi.gnum_range[0]; j < bi.gnum_range[1]; i++, j++)
      connect[i] = j;

    _write_block_connect_g(1,
                           bi.gnum_range[0],
                           bi.gnum_range[1],
                           connect,
                           w->comm,
                           f);

    BFT_FREE(connect);
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write "trivial" point elements to an EnSight Gold file in serial mode
 *
 * parameters:
 *   mesh <-- pointer to nodal mesh structure
 *   f    <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_export_point_elements_l(const fvm_nodal_t  *mesh,
                         _ensight_file_t     f)
{
  const cs_lnum_t   n_vertices = mesh->n_vertices;

  _write_string(f, "point");
  _write_int(f, (int)n_vertices);

  if (n_vertices == 0)
    return;

  if (f.tf != NULL) { /* Text mode */
    int i;
    for (i = 0; i < n_vertices; i++)
      fprintf(f.tf, "%10d\n", i+1);
  }

  else if (f.bf != NULL) { /* Binary mode */

    int32_t  k, j_end;
    int32_t  j = 1;
    int32_t  *buf = NULL;
    const int32_t  bufsize = n_vertices >  64 ? (n_vertices / 8) : n_vertices;

    BFT_MALLOC(buf, bufsize, int32_t);

    j_end = n_vertices + 1;
    while (j < j_end) {
      for (k = 0;  j < j_end  && k < bufsize; k++)
        buf[k] = j++;
      cs_file_write_global(f.bf, buf, sizeof(int32_t), k);
    }

    BFT_FREE(buf);
  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write indexed element lengths from a nodal mesh to an EnSight Gold
 * file in parallel mode
 *
 * parameters:
 *   w                  <-- pointer to writer structure
 *   global_element_num <-- global element numbering
 *   vertex_index       <-- pointer to element -> vertex index
 *   n_ranks            <-- number of processes in communicator
 *   f                  <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_write_lengths_g(const fvm_to_ensight_writer_t  *w,
                 const fvm_io_num_t             *global_element_num,
                 const cs_lnum_t                 vertex_index[],
                 _ensight_file_t                 f)
{
  cs_lnum_t   i;
  cs_block_dist_info_t   bi;

  int32_t  *part_lengths = NULL;
  int32_t  *block_lengths = NULL;

  cs_part_to_block_t  *d = NULL;

  const size_t min_block_size = w->min_block_size / sizeof(int32_t);
  const cs_lnum_t   n_elements
    = fvm_io_num_get_local_count(global_element_num);
  const cs_lnum_t   n_g_elements
    = fvm_io_num_get_global_count(global_element_num);
  const cs_gnum_t   *g_num
    = fvm_io_num_get_global_num(global_element_num);

  /* Allocate block buffer */

  bi = cs_block_dist_compute_sizes(w->rank,
                                   w->n_ranks,
                                   w->min_rank_step,
                                   min_block_size,
                                   n_g_elements);

  /* Build distribution structures */

  BFT_MALLOC(block_lengths, bi.gnum_range[1] - bi.gnum_range[0], int);
  BFT_MALLOC(part_lengths, n_elements, int32_t);

  for (i = 0; i < n_elements; i++)
    part_lengths[i] = vertex_index[i+1] - vertex_index[i];

  d = cs_part_to_block_create_by_gnum(w->comm, bi, n_elements, g_num);

  cs_part_to_block_copy_array(d,
                              CS_INT32,
                              1,
                              part_lengths,
                              block_lengths);

  cs_part_to_block_destroy(&d);
  BFT_FREE(part_lengths);

  /* Write to file */

  _write_block_connect_g(1,
                         bi.gnum_range[0],
                         bi.gnum_range[1],
                         block_lengths,
                         w->comm,
                         f);

  BFT_FREE(block_lengths);
}

/*----------------------------------------------------------------------------
 * Write block-distributed indexed element (polygons or polyhedra)
 * cell -> vertex connectivity to an EnSight Gold file in parallel mode.
 *
 * In text mode, zeroes may be used in place of extra vertex numbers
 * to indicate extra newlines between face -> vertex definitions.
 *
 * parameters:
 *   num_start     <-- global number of first element for this block
 *   num_end       <-- global number of past last element for this block
 *   block_index   <-- global connectivity block array
 *   block_connect <-> global connectivity block array
 *   comm          <-- associated MPI communicator
 *   f             <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_write_block_indexed(cs_gnum_t         num_start,
                     cs_gnum_t         num_end,
                     const cs_lnum_t   block_index[],
                     int32_t           block_connect[],
                     MPI_Comm          comm,
                     _ensight_file_t   f)
{
  cs_gnum_t block_size = 0, block_start = 0, block_end = 0;

  /* Prepare write to file */

  block_size = block_index[num_end - num_start];

  MPI_Scan(&block_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM, comm);
  block_end += 1;
  block_start = block_end - block_size;

  /* In Binary mode, all ranks have a file structure,
     we may use use a collective call */

  if (f.bf != NULL)
    cs_file_write_block(f.bf,
                        block_connect,
                        sizeof(int32_t),
                        1,
                        block_start,
                        block_end);

  /* If all ranks do not have a binary file structure pointer, then
     we are using a text file, open only on rank 0 */

  else {
    cs_lnum_t   i;
    int32_t *_block_vtx_num = NULL;
    cs_file_serializer_t *s = cs_file_serializer_create(sizeof(int32_t),
                                                        1,
                                                        block_start,
                                                        block_end,
                                                        0,
                                                        block_connect,
                                                        comm);

    do {
      cs_gnum_t j;
      cs_gnum_t range[2] = {block_start, block_end};
      _block_vtx_num = cs_file_serializer_advance(s, range);
      if (_block_vtx_num != NULL) { /* only possible on rank 0 */
        assert(f.tf != NULL);
        for (i = 0, j = range[0]; j < range[1]; i++, j++) {
          if (_block_vtx_num[i] != 0)
            fprintf(f.tf, "%10d", _block_vtx_num[i]);
          else
            fprintf(f.tf, "\n");
        }
      }
    } while (_block_vtx_num != NULL);

    cs_file_serializer_destroy(&s);
  }
}

/*----------------------------------------------------------------------------
 * Write indexed element (polygons or polyhedra) cell -> vertex connectivity
 * to an EnSight Gold file in parallel mode.
 *
 * In text mode, zeroes may be used in place of extra vertex numbers
 * to indicate extra newlines between face -> vertex definitions.
 *
 * parameters:
 *   w                  <-- pointer to writer structure
 *   global_vertex_num  <-- vertex global numbering
 *   global_element_num <-- global element numbering
 *   vertex_index       <-- element -> vertex index
 *   vertex_num         <-- element -> vertex number
 *   f                  <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_write_indexed_connect_g(const fvm_to_ensight_writer_t  *w,
                         const fvm_io_num_t             *global_element_num,
                         const cs_lnum_t                 vertex_index[],
                         const int32_t                   vertex_num[],
                         _ensight_file_t                 f)
{
  cs_block_dist_info_t bi;

  cs_gnum_t loc_size = 0, tot_size = 0, block_size = 0;
  cs_part_to_block_t  *d = NULL;
  cs_lnum_t   *block_index = NULL;
  int32_t  *block_vtx_num = NULL;
  size_t  min_block_size = w->min_block_size / sizeof(int32_t);

  const cs_gnum_t   n_g_elements
    = fvm_io_num_get_global_count(global_element_num);
  const cs_lnum_t   n_elements
    = fvm_io_num_get_local_count(global_element_num);
  const cs_gnum_t   *g_elt_num
    = fvm_io_num_get_global_num(global_element_num);

  /* Adjust min block size based on minimum element size */

  loc_size = vertex_index[n_elements];
  MPI_Allreduce(&loc_size, &tot_size, 1, CS_MPI_GNUM, MPI_SUM, w->comm);

  min_block_size /= (tot_size / n_g_elements);

  /* Allocate memory for additionnal indexes */

  bi = cs_block_dist_compute_sizes(w->rank,
                                   w->n_ranks,
                                   w->min_rank_step,
                                   min_block_size,
                                   n_g_elements);

  BFT_MALLOC(block_index, bi.gnum_range[1] - bi.gnum_range[0] + 1, cs_lnum_t);

  d = cs_part_to_block_create_by_gnum(w->comm, bi, n_elements, g_elt_num);

  cs_part_to_block_copy_index(d,
                              vertex_index,
                              block_index);

  block_size = block_index[bi.gnum_range[1] - bi.gnum_range[0]];

  BFT_MALLOC(block_vtx_num, block_size, int32_t);

  cs_part_to_block_copy_indexed(d,
                                CS_INT32,
                                vertex_index,
                                vertex_num,
                                block_index,
                                block_vtx_num);

  /* Write to file */

  _write_block_indexed(bi.gnum_range[0],
                       bi.gnum_range[1],
                       block_index,
                       block_vtx_num,
                       w->comm,
                       f);

  /* Free memory */

  BFT_FREE(block_vtx_num);
  cs_part_to_block_destroy(&d);
  BFT_FREE(block_index);
}

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to an EnSight Gold file in parallel mode
 *
 * parameters:
 *   w                 <-- pointer to writer structure
 *   export_section    <-- pointer to EnSight section helper structure
 *   global_vertex_num <-- pointer to vertex global numbering
 *   f                 <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polyhedra_g(const fvm_to_ensight_writer_t  *w,
                          const fvm_writer_section_t     *export_section,
                          const fvm_io_num_t             *global_vertex_num,
                          _ensight_file_t                 f)
{
  cs_lnum_t   i, j, k, l, face_id;

  cs_lnum_t   face_length, cell_length;
  cs_block_dist_info_t  bi;

  cs_part_to_block_t  *d = NULL;
  const fvm_writer_section_t  *current_section;

  /* Export number of faces per polyhedron */
  /*---------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_lengths_g(w,
                     section->global_element_num,
                     section->face_index,
                     f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Export number of vertices per face per polyhedron */
  /*---------------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    cs_gnum_t   block_size = 0, block_start = 0, block_end = 0;
    cs_lnum_t *block_index = NULL;

    size_t  min_block_size = w->min_block_size / sizeof(int32_t);
    int32_t  *part_face_len = NULL, *block_face_len = NULL;

    const fvm_nodal_section_t  *section = current_section->section;
    const cs_lnum_t   n_elements
      = fvm_io_num_get_local_count(section->global_element_num);
    const cs_gnum_t n_g_elements
      = fvm_io_num_get_global_count(section->global_element_num);
    const cs_gnum_t   *g_elt_num
      = fvm_io_num_get_global_num(section->global_element_num);

    /* Build local polyhedron face lengths information */

    BFT_MALLOC(part_face_len,
               section->face_index[section->n_elements],
               int32_t);

    k = 0;
    for (i = 0; i < section->n_elements; i++) {
      for (j = section->face_index[i]; j < section->face_index[i+1]; j++) {
        face_id = CS_ABS(section->face_num[j]) - 1;
        face_length = (  section->vertex_index[face_id+1]
                       - section->vertex_index[face_id]);
        part_face_len[k++] = face_length;
      }
    }
    assert(k == section->face_index[section->n_elements]);

    /* Prepare distribution structures */

    bi = cs_block_dist_compute_sizes(w->rank,
                                     w->n_ranks,
                                     w->min_rank_step,
                                     min_block_size,
                                     n_g_elements);

    d = cs_part_to_block_create_by_gnum(w->comm,
                                        bi,
                                        n_elements,
                                        g_elt_num);

    BFT_MALLOC(block_index, bi.gnum_range[1] - bi.gnum_range[0] + 1, cs_lnum_t);

    cs_part_to_block_copy_index(d,
                                section->face_index,
                                block_index);

    block_size = block_index[bi.gnum_range[1] - bi.gnum_range[0]];

    BFT_MALLOC(block_face_len, block_size, int32_t);

    cs_part_to_block_copy_indexed(d,
                                  CS_INT32,
                                  section->face_index,
                                  part_face_len,
                                  block_index,
                                  block_face_len);

    MPI_Scan(&block_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM, w->comm);
    block_end += 1;
    block_start = block_end - block_size;

    _write_block_connect_g(1,
                           block_start,
                           block_end,
                           block_face_len,
                           w->comm,
                           f);

    BFT_FREE(block_face_len);

    cs_part_to_block_destroy(&d);

    BFT_FREE(block_index);
    BFT_FREE(part_face_len);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Export cell->vertex connectivity by blocks */
  /*--------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    cs_lnum_t   *part_vtx_idx = NULL;
    int32_t  *part_vtx_num = NULL;

    const fvm_nodal_section_t  *section = current_section->section;
    const cs_gnum_t   *g_vtx_num
      = fvm_io_num_get_global_num(global_vertex_num);

    BFT_MALLOC(part_vtx_idx, section->n_elements + 1, cs_lnum_t);

    l = 0;

    if (f.bf != NULL) { /* In binary mode, build cell -> vertex connectivity */

      part_vtx_idx[0] = 0;
      for (i = 0; i < section->n_elements; i++) {
        cell_length = 0;
        for (j = section->face_index[i]; j < section->face_index[i+1]; j++) {
          face_id = CS_ABS(section->face_num[j]) - 1;
          face_length = (  section->vertex_index[face_id+1]
                         - section->vertex_index[face_id]);
          cell_length += face_length;
        }
        part_vtx_idx[i+1] = part_vtx_idx[i] + cell_length;
      }

    }

    /* In text mode, add zeroes to cell vertex connectivity to mark face
       bounds (so as to add newlines) */

    else { /* we are in text mode if f.bf == NULL */

      part_vtx_idx[0] = 0;
      for (i = 0; i < section->n_elements; i++) {
        cell_length = 0;
        for (j = section->face_index[i]; j < section->face_index[i+1]; j++) {
          face_id = CS_ABS(section->face_num[j]) - 1;
          face_length = (  section->vertex_index[face_id+1]
                         - section->vertex_index[face_id]);
          cell_length += face_length + 1;
        }
        part_vtx_idx[i+1] = part_vtx_idx[i] + cell_length;
      }

    }

    BFT_MALLOC(part_vtx_num, part_vtx_idx[section->n_elements], cs_lnum_t);

    l = 0;

    for (i = 0; i < section->n_elements; i++) {
      for (j = section->face_index[i]; j < section->face_index[i+1]; j++) {
        if (section->face_num[j] > 0) {
          face_id = section->face_num[j] - 1;
          for (k = section->vertex_index[face_id];
               k < section->vertex_index[face_id+1];
               k++)
            part_vtx_num[l++] = g_vtx_num[section->vertex_num[k] - 1];
        }
        else {
          face_id = -section->face_num[j] - 1;
          k = section->vertex_index[face_id];
          part_vtx_num[l++] = g_vtx_num[section->vertex_num[k] - 1];
          for (k = section->vertex_index[face_id+1] - 1;
               k > section->vertex_index[face_id];
               k--)
            part_vtx_num[l++] = g_vtx_num[section->vertex_num[k] - 1];
        }
        if (f.bf == NULL)
          part_vtx_num[l++] = 0; /* mark face limits in text mode */
      }

    }

    /* Now distribute and write cells -> vertices connectivity */

    _write_indexed_connect_g(w,
                             section->global_element_num,
                             part_vtx_idx,
                             part_vtx_num,
                             f);

    BFT_FREE(part_vtx_num);
    BFT_FREE(part_vtx_idx);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to an EnSight Gold file in serial mode
 *
 * parameters:
 *   export_section <-- pointer to EnSight section helper structure
 *   f              <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polyhedra_l(const fvm_writer_section_t  *export_section,
                          _ensight_file_t              f)
{
  int  face_sgn;
  cs_lnum_t   i, j, k, l;
  size_t  i_buf;

  cs_lnum_t   face_length, face_id;

  size_t    buffer_size = 0;
  int32_t  *buffer = NULL;

  const fvm_writer_section_t  *current_section;

  /* Write number of faces per cell */
  /*--------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    if (f.tf != NULL) { /* Text mode */
      for (i = 0; i < section->n_elements; i++)
        fprintf(f.tf, "%10d\n",
                (int)(  section->face_index[i+1]
                      - section->face_index[i]));
    }
    else { /* binary mode */

      /* First, allocate a buffer large enough so that the number of
         writes is limited, small enough so that the memory overhead is
         minimal; polyhedral connectivity is at least 4 faces x 3 vertices
         per cell, usually quite a bit more, so this is 1/3 of the minimum */

      if (buffer_size < (size_t)section->n_elements * 4) {
        buffer_size = section->n_elements * 4;
        BFT_REALLOC(buffer, buffer_size, int32_t);
      }

      /* Now fill buffer and write */

      for (i = 0, i_buf = 0; i < section->n_elements; i++) {
        if (i_buf == buffer_size) {
          cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);
          i_buf = 0;
        }
        buffer[i_buf++] = (int)(  section->face_index[i+1]
                                - section->face_index[i]);
      }
      if (i_buf > 0)
        cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);

    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Write number of vertices/face */
  /*-------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    for (i = 0, i_buf = 0; i < section->n_elements; i++) {

      /* Loop on cell faces */

      for (j = section->face_index[i];
           j < section->face_index[i+1];
           j++) {

        if (section->face_num[j] > 0)
          face_id = section->face_num[j] - 1;
        else
          face_id = -section->face_num[j] - 1;

        face_length = (  section->vertex_index[face_id+1]
                       - section->vertex_index[face_id]);

        if (f.tf != NULL)
          fprintf(f.tf, "%10d\n",
                  (int)face_length);
        else {
          if (i_buf == buffer_size) {
            cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);
            i_buf = 0;
          }
          buffer[i_buf++] = (int)face_length;
        }

      }

    }

    if (f.bf != NULL && i_buf > 0)
      cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Write cell/vertex connectivity */
  /*--------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    for (i = 0, i_buf = 0; i < section->n_elements; i++) {

      /* Loop on cell faces */

      for (j = section->face_index[i];
           j < section->face_index[i+1];
           j++) {

        /* Print face vertex numbers */

        if (section->face_num[j] > 0) {
          face_id = section->face_num[j] - 1;
          face_sgn = 1;
        }
        else {
          face_id = -section->face_num[j] - 1;
          face_sgn = -1;
        }

        face_length = (  section->vertex_index[face_id+1]
                       - section->vertex_index[face_id]);

        if (f.tf != NULL) { /* text mode */
          for (k = 0; k < face_length; k++) {
            l =   section->vertex_index[face_id]
                + (face_length + (k*face_sgn))%face_length;
            fprintf(f.tf, "%10d", (int)section->vertex_num[l]);
          }
          fprintf(f.tf, "\n");
        }
        else { /* binary mode */
          for (k = 0; k < face_length; k++) {
            l =   section->vertex_index[face_id]
                + (face_length + (k*face_sgn))%face_length;
            if (i_buf == buffer_size) {
              cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);
              i_buf = 0;
            }
            buffer[i_buf++] = (int)section->vertex_num[l];
          }
        }

      } /* End of loop on cell faces */

    } /* End of loop on polyhedral cells */

    if (f.bf != NULL && i_buf > 0)
      cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  if (buffer != NULL)
    BFT_FREE(buffer);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to an EnSight Gold file in parallel mode
 *
 * parameters:
 *   w                 <-- pointer to writer structure
 *   export_section    <-- pointer to EnSight section helper structure
 *   global_vertex_num <-- pointer to vertex global numbering
 *   f                 <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_g(const fvm_to_ensight_writer_t  *w,
                         const fvm_writer_section_t     *export_section,
                         const fvm_io_num_t             *global_vertex_num,
                         _ensight_file_t                 f)
{
  const fvm_writer_section_t  *current_section;

  /* Export number of vertices per polygon by blocks */
  /*-------------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_lengths_g(w,
                     section->global_element_num,
                     section->vertex_index,
                     f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Export face->vertex connectivity by blocks */
  /*--------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    cs_lnum_t   i, j, k;
    cs_lnum_t   *_part_vtx_idx = NULL;
    const cs_lnum_t   *part_vtx_idx = NULL;
    int32_t  *part_vtx_num = NULL;

    const fvm_nodal_section_t  *section = current_section->section;
    const cs_gnum_t   *g_vtx_num
      = fvm_io_num_get_global_num(global_vertex_num);

    if (f.bf != NULL) /* In binary mode, use existing index */
      part_vtx_idx = section->vertex_index;

    /* In text mode, add zeroes to cell vertex connectivity to mark face
       bounds (so as to add newlines) */

    else { /* we are in text mode if f.bf == NULL */

      BFT_MALLOC(_part_vtx_idx, section->n_elements + 1, cs_lnum_t);

      _part_vtx_idx[0] = 0;
      for (i = 0; i < section->n_elements; i++)
        _part_vtx_idx[i+1] = _part_vtx_idx[i] + (  section->vertex_index[i+1]
                                                 - section->vertex_index[i]) + 1;

      part_vtx_idx = _part_vtx_idx;
    }

    /* Build connectivity array */

    BFT_MALLOC(part_vtx_num, part_vtx_idx[section->n_elements], int32_t);

    if (f.bf != NULL) { /* binary mode */
      for (i = 0, k = 0; i < section->n_elements; i++) {
        for (j = section->vertex_index[i];
             j < section->vertex_index[i+1];
             j++)
          part_vtx_num[k++] = g_vtx_num[section->vertex_num[j] - 1];
      }
    }

    else { /* text mode */
      for (i = 0, k = 0; i < section->n_elements; i++) {
        for (j = section->vertex_index[i];
             j < section->vertex_index[i+1];
             j++)
          part_vtx_num[k++] = g_vtx_num[section->vertex_num[j] - 1];
        part_vtx_num[k++] = 0; /* mark face bounds in text mode */
      }
    }

    /* Now distribute and write cell -> vertices connectivity */

    _write_indexed_connect_g(w,
                             section->global_element_num,
                             part_vtx_idx,
                             part_vtx_num,
                             f);

    BFT_FREE(part_vtx_num);
    if (_part_vtx_idx != NULL)
      BFT_FREE(_part_vtx_idx);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to a text file in serial mode
 *
 * parameters:
 *   export_section <-- pointer to EnSight section helper structure
 *   f              <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
*----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_l(const fvm_writer_section_t  *export_section,
                         _ensight_file_t              f)


{
  cs_lnum_t   i, j;
  size_t  i_buf;

  size_t    buffer_size = 0;
  int32_t  *buffer = NULL;

  const fvm_writer_section_t  *current_section = NULL;

  /* Print face connectivity directly, without using extra buffers */

  /* First loop on all polygonal faces, to write number of vertices */
  /*----------------------------------------------------------------*/

  current_section = export_section;

  do { /* Loop on sections that should be grouped */

    const fvm_nodal_section_t  *section = current_section->section;

    if (f.tf != NULL) { /* Text mode */
      for (i = 0; i < section->n_elements; i++)
        fprintf(f.tf, "%10d\n", (int)(  section->vertex_index[i+1]
                                      - section->vertex_index[i]));
    }
    else { /* binary mode */

      /* First, allocate a buffer large enough so that the number of
         writes is limited, small enough so that the memory overhead is
         minimal; polygonal connectivity is at least 3 vertices per face,
         usually 5 or more, so this is 1/3 of the minimum */

      if (buffer_size < (size_t)section->n_elements) {
        buffer_size = section->n_elements;
        BFT_REALLOC(buffer, buffer_size, int32_t);
      }

      /* Now fill buffer and write */

      for (i = 0, i_buf = 0; i < section->n_elements; i++) {
        if (i_buf == buffer_size) {
          cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);
          i_buf = 0;
        }
        buffer[i_buf++] = (int)(  section->vertex_index[i+1]
                                - section->vertex_index[i]);
      }
      if (i_buf > 0)
        cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Loop on all polygonal faces */
  /*-----------------------------*/

  current_section = export_section;

  do { /* Loop on sections that should be grouped */

    const fvm_nodal_section_t  *section = current_section->section;

    for (i = 0, i_buf = 0; i < section->n_elements; i++) {

      /* Print face vertex numbers */

      if (f.tf != NULL) { /* text mode */
        for (j = section->vertex_index[i];
             j < section->vertex_index[i+1];
             j++)
          fprintf(f.tf, "%10d", (int)section->vertex_num[j]);
        fprintf(f.tf, "\n");
      }
      else { /* binary mode */
        for (j = section->vertex_index[i];
             j < section->vertex_index[i+1];
             j++) {
          if (i_buf == buffer_size) {
            cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);
            i_buf = 0;
          }
          buffer[i_buf++] = (int)section->vertex_num[j];
        }
      }

    } /* End of loop on polygonal faces */

    if (f.bf != NULL && i_buf > 0)
      cs_file_write_global(f.bf, buffer, sizeof(int32_t), i_buf);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  if (buffer != NULL)
    BFT_FREE(buffer);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write tesselated element cell -> vertex connectivity to an EnSight Gold
 * file in parallel mode.
 *
 * parameters:
 *   w                  <-- pointer to writer structure
 *   global_vertex_num  <-- vertex global numbering
 *   global_element_num <-- global element numbering
 *   tesselation        <-- element tesselation description
 *   type               <-- tesselated sub-element type
 *   extra_vertex_base  <-- starting number for added vertices
 *   f                  <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_write_tesselated_connect_g(const fvm_to_ensight_writer_t  *w,
                            const fvm_io_num_t             *global_vertex_num,
                            const fvm_io_num_t             *global_element_num,
                            const fvm_tesselation_t        *tesselation,
                            fvm_element_t                   type,
                            const cs_gnum_t                 extra_vertex_base,
                            _ensight_file_t                 f)
{
  cs_lnum_t   i;
  cs_block_dist_info_t bi;

  cs_lnum_t   part_size = 0;

  cs_lnum_t   end_id = 0;
  cs_gnum_t   n_g_sub_elements = 0, global_num_end = 0;
  cs_gnum_t   block_size = 0, block_start = 0, block_end = 0;

  cs_part_to_block_t  *d = NULL;
  cs_lnum_t   *part_index, *block_index = NULL;
  int32_t     *part_vtx_num = NULL, *block_vtx_num = NULL;
  cs_gnum_t   *part_vtx_gnum = NULL;

  size_t  min_block_size = w->min_block_size / sizeof(int32_t);

  const int  stride = fvm_nodal_n_vertices_element[type];
  const cs_lnum_t   n_elements = fvm_tesselation_n_elements(tesselation);
  const cs_gnum_t   n_g_elements
    = fvm_io_num_get_global_count(global_element_num);
  const cs_lnum_t   n_sub_elements
    = fvm_tesselation_n_sub_elements(tesselation, type);
  const cs_lnum_t   *sub_element_idx
      = fvm_tesselation_sub_elt_index(tesselation, type);
  const cs_gnum_t   *g_elt_num
    = fvm_io_num_get_global_num(global_element_num);

  /* Adjust min block size based on mean number of sub-elements */

  fvm_tesselation_get_global_size(tesselation,
                                  type,
                                  &n_g_sub_elements,
                                  NULL);

  min_block_size /= (n_g_sub_elements/n_g_elements) * stride;

  /* Decode connectivity */

  part_size = n_sub_elements * stride;
  assert(sub_element_idx[n_elements]*stride == part_size);

  global_num_end = n_g_elements + 1;

  if (n_elements > 0) {
    BFT_MALLOC(part_vtx_num, part_size, int32_t);
    BFT_MALLOC(part_vtx_gnum, part_size, cs_gnum_t);
  }

  end_id = fvm_tesselation_decode_g(tesselation,
                                    type,
                                    0,
                                    part_size,
                                    &global_num_end,
                                    global_vertex_num,
                                    extra_vertex_base,
                                    part_vtx_gnum,
                                    w->comm);

  assert(end_id == n_elements);
  assert(global_num_end == n_g_elements + 1);

  /* Convert to write type */

  if (n_elements > 0) {
    for (i = 0; i < part_size; i++)
      part_vtx_num[i] = part_vtx_gnum[i];
    BFT_FREE(part_vtx_gnum);
  }

  /* Allocate memory for additionnal indexes and decoded connectivity */

  bi = cs_block_dist_compute_sizes(w->rank,
                                   w->n_ranks,
                                   w->min_rank_step,
                                   min_block_size,
                                   n_g_elements);

  BFT_MALLOC(block_index, bi.gnum_range[1] - bi.gnum_range[0] + 1, cs_lnum_t);
  BFT_MALLOC(part_index, n_elements + 1, cs_lnum_t);

  d = cs_part_to_block_create_by_gnum(w->comm, bi, n_elements, g_elt_num);

  part_index[0] = 0;
  for (i = 0; i < n_elements; i++) {
    part_index[i+1] = part_index[i] + (  sub_element_idx[i+1]
                                       - sub_element_idx[i]) * stride;
  }

  /* Copy index */

  cs_part_to_block_copy_index(d,
                              part_index,
                              block_index);

  block_size = (block_index[bi.gnum_range[1] - bi.gnum_range[0]]);

  /* Copy connectivity */

  BFT_MALLOC(block_vtx_num, block_size, int32_t);

  cs_part_to_block_copy_indexed(d,
                                CS_INT32,
                                part_index,
                                part_vtx_num,
                                block_index,
                                block_vtx_num);

  cs_part_to_block_destroy(&d);

  BFT_FREE(part_vtx_num);
  BFT_FREE(part_index);
  BFT_FREE(block_index);

  /* Write to file */

  block_size /= stride;

  MPI_Scan(&block_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM, w->comm);
  block_end += 1;
  block_start = block_end - block_size;

  _write_block_connect_g(stride,
                         block_start,
                         block_end,
                         block_vtx_num,
                         w->comm,
                         f);

  /* Free remaining memory */

  BFT_FREE(block_vtx_num);
}

/*----------------------------------------------------------------------------
 * Write tesselated element connectivity from a nodal mesh to an EnSight Gold
 * file in parallel mode
 *
 * parameters:
 *   w                 <-- pointer to writer structure
 *   export_section    <-- pointer to EnSight section helper structure
 *   global_vertex_num <-- pointer to vertex global numbering
 *   f                 <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_tesselated_g(const fvm_to_ensight_writer_t  *w,
                           const fvm_writer_section_t     *export_section,
                           const fvm_io_num_t             *global_vertex_num,
                           _ensight_file_t                 f)
{
  const fvm_writer_section_t  *current_section;

  /* Export face->vertex connectivity by blocks */
  /*--------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_tesselated_connect_g(w,
                                global_vertex_num,
                                section->global_element_num,
                                section->tesselation,
                                current_section->type,
                                current_section->extra_vertex_base,
                                f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write tesselated element connectivity from a nodal mesh to an EnSight Gold
 * file in parallel mode
 *
 * parameters:
 *   export_section <-- pointer to EnSight section helper structure
 *   f              <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_tesselated_l(const fvm_writer_section_t  *export_section,
                           _ensight_file_t              f)
{
  const fvm_writer_section_t  *current_section;

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    cs_lnum_t   start_id, end_id;
    cs_lnum_t   n_sub_elements_max;
    cs_lnum_t   n_buffer_elements_max = section->n_elements;
    cs_lnum_t *vertex_num = NULL;

    const cs_lnum_t *sub_element_idx
      = fvm_tesselation_sub_elt_index(section->tesselation,
                                      export_section->type);

    fvm_tesselation_get_global_size(section->tesselation,
                                    export_section->type,
                                    NULL,
                                    &n_sub_elements_max);
    if (n_sub_elements_max > n_buffer_elements_max)
      n_buffer_elements_max = n_sub_elements_max;

    BFT_MALLOC(vertex_num,
               (  n_buffer_elements_max
                * fvm_nodal_n_vertices_element[export_section->type]),
               cs_lnum_t);

    for (start_id = 0;
         start_id < section->n_elements;
         start_id = end_id) {

      end_id
        = fvm_tesselation_decode(section->tesselation,
                                 current_section->type,
                                 start_id,
                                 n_buffer_elements_max,
                                 export_section->extra_vertex_base,
                                 vertex_num);

      _write_connect_l(fvm_nodal_n_vertices_element[export_section->type],
                       (  sub_element_idx[end_id]
                        - sub_element_idx[start_id]),
                       vertex_num,
                       f);

    }

    BFT_FREE(vertex_num);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided elements from a nodal mesh to an EnSight Gold file in
 * parallel mode
 *
 * parameters:
 *   w                 <-- pointer to writer structure
 *   export_section    <-- pointer to EnSight section helper structure
 *   global_vertex_num <-- pointer to vertex global numbering
 *   f                 <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_strided_g(const fvm_to_ensight_writer_t  *w,
                        const fvm_writer_section_t     *export_section,
                        const fvm_io_num_t             *global_vertex_num,
                        _ensight_file_t                 f)
{
  cs_lnum_t   i, j;

  const fvm_writer_section_t  *current_section;

  /* Export vertex connectivity */

  current_section = export_section;

  do { /* loop on sections which should be appended */

    cs_block_dist_info_t bi;

    cs_lnum_t   block_size = 0;
    cs_part_to_block_t  *d = NULL;
    int32_t  *part_vtx_num = NULL, *block_vtx_num = NULL;

    const fvm_nodal_section_t  *section = current_section->section;
    const int  stride = fvm_nodal_n_vertices_element[section->type];

    const size_t  min_block_size
      = w->min_block_size / (sizeof(int32_t) * stride);

    const cs_lnum_t   n_elements
      = fvm_io_num_get_local_count(section->global_element_num);
    const cs_gnum_t   n_g_elements
      = fvm_io_num_get_global_count(section->global_element_num);
    const cs_gnum_t   *g_elt_num
      = fvm_io_num_get_global_num(section->global_element_num);
    const cs_gnum_t   *g_vtx_num
      = fvm_io_num_get_global_num(global_vertex_num);

    /* Prepare distribution structures */

    bi = cs_block_dist_compute_sizes(w->rank,
                                     w->n_ranks,
                                     w->min_rank_step,
                                     min_block_size,
                                     n_g_elements);

    d = cs_part_to_block_create_by_gnum(w->comm,
                                        bi,
                                        n_elements,
                                        g_elt_num);

    /* Build connectivity */

    block_size = bi.gnum_range[1] - bi.gnum_range[0];

    BFT_MALLOC(block_vtx_num, block_size*stride, int32_t);
    BFT_MALLOC(part_vtx_num, n_elements*stride, int32_t);

    for (i = 0; i < n_elements; i++) {
      for (j = 0; j < stride; j++) {
        part_vtx_num[i*stride + j]
          = g_vtx_num[section->vertex_num[i*stride + j] - 1];
      }
    }

    cs_part_to_block_copy_array(d,
                                CS_INT32,
                                stride,
                                part_vtx_num,
                                block_vtx_num);

    BFT_FREE(part_vtx_num);

    _write_block_connect_g(stride,
                           bi.gnum_range[0],
                           bi.gnum_range[1],
                           block_vtx_num,
                           w->comm,
                           f);

    BFT_FREE(block_vtx_num);

    cs_part_to_block_destroy(&d);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

/*----------------------------------------------------------------------------
 * Write field values associated with nodal values of a nodal mesh to
 * an EnSight Gold file in serial mode.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   w                <-- pointer to writer structure
 *   mesh             <-- pointer to nodal mesh structure
 *   divide_polyhedra <-- true if polyhedra are tesselated
 *   input_dim        <-- input field dimension
 *   output_dim       <-- output field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- input data type (output is real)
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_export_field_values_ng(const fvm_to_ensight_writer_t  *w,
                        const fvm_nodal_t              *mesh,
                        int                             input_dim,
                        int                             output_dim,
                        cs_interlace_t                  interlace,
                        int                             n_parent_lists,
                        const cs_lnum_t                 parent_num_shift[],
                        cs_datatype_t                   datatype,
                        const void               *const field_values[],
                        _ensight_file_t                 f)
{
  int  i;
  cs_block_dist_info_t  bi;

  cs_lnum_t   part_size = 0, block_size = 0;
  float  *part_values = NULL, *block_values = NULL;
  cs_part_to_block_t  *d = NULL;

  /* Initialize distribution info */

  _vertex_part_to_block_create(w,
                               mesh,
                               &bi,
                               &d);

  part_size = cs_part_to_block_get_n_part_ents(d);
  block_size = bi.gnum_range[1] - bi.gnum_range[0];

  BFT_MALLOC(part_values, part_size, float);
  BFT_MALLOC(block_values, block_size, float);

  for (i = 0; i < output_dim; i++) {

    cs_lnum_t   start_id = 0;
    cs_lnum_t   end_id = mesh->n_vertices;

    /* Distribute partition to block values */

    if (i < input_dim) {

      int j;

      /* Main vertices */

      fvm_convert_array(input_dim,
                        i,
                        1, /* stride */
                        start_id,
                        end_id,
                        interlace,
                        datatype,
                        CS_FLOAT,
                        n_parent_lists,
                        parent_num_shift,
                        mesh->parent_vertex_num,
                        field_values,
                        part_values);

      /* Additional vertices in case of tesselation
         (end_id == part_size with no tesselation or if all tesselated
         sections have been accounted for).*/

      for (j = 0; end_id < part_size && j < mesh->n_sections; j++) {

        const fvm_nodal_section_t  *section = mesh->sections[j];

        assert(w->divide_polyhedra == true);

        if (section->type == FVM_CELL_POLY && section->tesselation != NULL) {

          cs_lnum_t   n_extra_vertices
            = fvm_tesselation_n_vertices_add(section->tesselation);

          start_id = end_id;
          end_id = start_id + n_extra_vertices;

          fvm_tesselation_vertex_values(section->tesselation,
                                        input_dim,
                                        i,
                                        1, /* stride, */
                                        0,
                                        n_extra_vertices,
                                        interlace,
                                        datatype,
                                        CS_FLOAT,
                                        n_parent_lists,
                                        parent_num_shift,
                                        mesh->parent_vertex_num,
                                        field_values,
                                        part_values + start_id);

        }

      } /* End of loops on tesselated sections */

      assert(end_id == part_size);

      cs_part_to_block_copy_array(d,
                                  CS_FLOAT,
                                  1,
                                  part_values,
                                  block_values);

    }

    /* Zero extra dimensions */

    else {
      cs_lnum_t j;
      for (j = 0; j < block_size; j++)
        block_values[j] = 0.0;
    }

    _write_block_floats_g(bi.gnum_range[0],
                          bi.gnum_range[1],
                          block_values,
                          w->comm,
                          f);
  }

  BFT_FREE(block_values);
  BFT_FREE(part_values);

  cs_part_to_block_destroy(&d);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write field values associated with nodal values of a nodal mesh to
 * an EnSight Gold file in serial mode.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   n_entities         <-- number of entities
 *   input_dim          <-- input field dimension
 *   output_dim         <-- output field dimension
 *   interlace          <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- input data type (output is real)
 *   field_values       <-- array of associated field value arrays
 *   f                  <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_export_field_values_nl(const fvm_nodal_t           *mesh,
                        fvm_writer_field_helper_t   *helper,
                        int                          input_dim,
                        cs_interlace_t               interlace,
                        int                          n_parent_lists,
                        const cs_lnum_t              parent_num_shift[],
                        cs_datatype_t                datatype,
                        const void            *const field_values[],
                        _ensight_file_t              f)
{
  int  i;
  size_t  output_size;
  float  *output_buffer;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  const size_t  output_buffer_size
    = mesh->n_vertices >  16 ? (mesh->n_vertices / 4) : mesh->n_vertices;

  BFT_MALLOC(output_buffer, output_buffer_size, float);

  for (i = 0; i < output_dim; i++) {

    while (fvm_writer_field_helper_step_n(helper,
                                          mesh,
                                          input_dim,
                                          i,
                                          interlace,
                                          n_parent_lists,
                                          parent_num_shift,
                                          datatype,
                                          field_values,
                                          output_buffer,
                                          output_buffer_size,
                                          &output_size) == 0) {

      _write_block_floats_l(output_size,
                            output_buffer,
                            f);

    }
  }

  BFT_FREE(output_buffer);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to
 * an EnSight Gold file.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   w                <-- pointer to writer structure
 *   export_section   <-- pointer to EnSight section helper structure
 *   helper           <-- pointer to general writer helper structure
 *   input_dim        <-- input field dimension
 *   output_dim       <-- output field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_field_values_eg(const fvm_to_ensight_writer_t   *w,
                        const fvm_writer_section_t      *export_section,
                        int                              input_dim,
                        int                              output_dim,
                        cs_interlace_t                   interlace,
                        int                              n_parent_lists,
                        const cs_lnum_t                  parent_num_shift[],
                        cs_datatype_t                    datatype,
                        const void                *const field_values[],
                        _ensight_file_t                  f)
{
  int  i;
  cs_lnum_t   j, k;

  cs_block_dist_info_t  bi;
  cs_part_to_block_t  *d = NULL;

  int         n_sections = 0;
  _Bool       have_tesselation = false;
  cs_lnum_t   part_size = 0, block_size = 0;
  cs_gnum_t   block_sub_size = 0, block_start = 0, block_end = 0;
  cs_gnum_t   n_g_elements = 0;

  int  *part_n_sub = NULL, *block_n_sub = NULL;
  float  *part_values = NULL, *block_values = NULL, *_block_values = NULL;

  cs_gnum_t         *_g_elt_num = NULL;
  const cs_gnum_t   *g_elt_num
    = fvm_io_num_get_global_num(export_section->section->global_element_num);

  const fvm_writer_section_t  *current_section = NULL;

  size_t  min_block_size = w->min_block_size / sizeof(int32_t);

  /* Loop on sections to count output size */

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;

    n_sections += 1;
    n_g_elements += fvm_io_num_get_global_count(section->global_element_num);
    part_size += fvm_io_num_get_local_count(section->global_element_num);
    if (current_section->type != section->type)
      have_tesselation = true;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Build global numbering if necessary */

  if (n_sections > 1) {

    cs_lnum_t start_id = 0;

    BFT_MALLOC(_g_elt_num, part_size, cs_gnum_t);
    g_elt_num = _g_elt_num;

    /* loop on sections which should be appended */

    current_section = export_section;
    do {

      const fvm_nodal_section_t  *section = current_section->section;
      const cs_lnum_t section_size
        = fvm_io_num_get_local_count(section->global_element_num);

      memcpy(_g_elt_num + start_id,
             fvm_io_num_get_global_num(section->global_element_num),
             sizeof(cs_gnum_t)*section_size);
      start_id += section_size;

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous == true);
  }

  /* Build sub-element count if necessary */

  if (have_tesselation) {

    cs_lnum_t start_id = 0;

    BFT_MALLOC(part_n_sub, part_size, int);

    current_section = export_section;
    do {

      const fvm_nodal_section_t  *section = current_section->section;
      const cs_lnum_t section_size
        = fvm_io_num_get_local_count(section->global_element_num);

      if (current_section->type != section->type) {
        const cs_lnum_t   *sub_element_idx
          = fvm_tesselation_sub_elt_index(section->tesselation,
                                          current_section->type);
        for (j = 0; j < section_size; j++)
          part_n_sub[start_id + j] = sub_element_idx[j+1] - sub_element_idx[j];
      }
      else {
        for (j = 0; j < section_size; j++)
          part_n_sub[start_id + j] = 1;
      }
      start_id += section_size;

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous == true);
  }

  /* Build distribution structures */

  bi = cs_block_dist_compute_sizes(w->rank,
                                   w->n_ranks,
                                   w->min_rank_step,
                                   min_block_size,
                                   n_g_elements);

  block_size = bi.gnum_range[1] - bi.gnum_range[0];

  d = cs_part_to_block_create_by_gnum(w->comm, bi, part_size, g_elt_num);

  if (_g_elt_num != NULL)
    cs_part_to_block_transfer_gnum(d, _g_elt_num);

  g_elt_num = NULL;
  _g_elt_num = NULL;

  /* Distribute sub-element info in case of tesselation */

  if (have_tesselation) {

    BFT_MALLOC(block_n_sub, block_size, int);
    cs_part_to_block_copy_array(d,
                                CS_INT32,
                                1,
                                part_n_sub,
                                block_n_sub);
    BFT_FREE(part_n_sub);

    for (j = 0; j < block_size; j++)
      block_sub_size += block_n_sub[j];

  }
  else
    block_sub_size = block_size;

  /* To save space, in case of tesselation, part_values and _block_n_sub
     point to the same memory space, as they are not needed simultaneously.
     Without tesselation, _block_n_sub simply points to block_n_sub */

  BFT_MALLOC(part_values,
             CS_MAX(part_size, (cs_lnum_t)block_sub_size),
             float);
  BFT_MALLOC(block_values, block_size, float);

  if (have_tesselation) {
    MPI_Scan(&block_sub_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM, w->comm);
    block_end += 1;
    block_start = block_end - block_sub_size;
    _block_values = part_values;
  }
  else {
    block_start = bi.gnum_range[0];
    block_end = bi.gnum_range[1];
    _block_values = block_values;
  }

  /* Loop on dimension (de-interlace vectors, always 3D for EnSight) */

  for (i = 0; i < output_dim; i++) {

    /* Distribute partition to block values */

    if (i < input_dim) {

      cs_lnum_t start_id = 0;
      cs_lnum_t src_shift = 0;

      /* loop on sections which should be appended */

      current_section = export_section;
      do {

        const fvm_nodal_section_t  *section = current_section->section;

        if (n_parent_lists == 0)
          src_shift = export_section->num_shift;

        fvm_convert_array(input_dim,
                          i,
                          1,
                          src_shift,
                          section->n_elements + src_shift,
                          interlace,
                          datatype,
                          CS_FLOAT,
                          n_parent_lists,
                          parent_num_shift,
                          section->parent_element_num,
                          field_values,
                          part_values + start_id);

        start_id += fvm_io_num_get_local_count(section->global_element_num);

        current_section = current_section->next;

      } while (   current_section != NULL
               && current_section->continues_previous == true);

      /* Distribute part values */

      cs_part_to_block_copy_array(d,
                                  CS_FLOAT,
                                  1,
                                  part_values,
                                  block_values);

      /* Scatter values to sub-elements in case of tesselation */

      if (have_tesselation) {
        cs_lnum_t   l = 0;
        for (j = 0; j < block_size; j++) {
          for (k = 0; k < block_n_sub[j]; k++)
            _block_values[l++] = block_values[j];
        }
      }

    }

    /* Zero extra dimensions */

    else {
      for (j = 0; j < (cs_lnum_t)block_sub_size; j++)
        block_values[j] = 0.0;
    }

    /* Write block values */

    _write_block_floats_g(block_start,
                          block_end,
                          _block_values,
                          w->comm,
                          f);

  } /* end of loop on spatial dimension */

  BFT_FREE(block_values);
  BFT_FREE(part_values);

  cs_part_to_block_destroy(&d);

  if (block_n_sub != NULL)
    BFT_FREE(block_n_sub);

  return current_section;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to
 * an EnSight Gold file.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   export_section   <-- pointer to EnSight section helper structure
 *   helper           <-- pointer to general writer helper structure
 *   input_dim        <-- input field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated file handle
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_field_values_el(const fvm_writer_section_t      *export_section,
                        fvm_writer_field_helper_t       *helper,
                        int                              input_dim,
                        cs_interlace_t                   interlace,
                        int                              n_parent_lists,
                        const cs_lnum_t                  parent_num_shift[],
                        cs_datatype_t                    datatype,
                        const void                *const field_values[],
                        _ensight_file_t                  f)
{
  int  i;
  size_t  input_size = 0, output_size = 0;
  size_t  min_output_buffer_size = 0, output_buffer_size = 0;
  float  *output_buffer = NULL;

  const fvm_writer_section_t  *current_section = NULL;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  /* Blocking for arbitrary buffer size, but should be small enough
     to add little additional memory requirement (in proportion), large
     enough to limit number of write and gather calls. */

  fvm_writer_field_helper_get_size(helper,
                                   &input_size,
                                   &output_size,
                                   NULL,
                                   &min_output_buffer_size);

  output_buffer_size = input_size / 4;
  output_buffer_size = CS_MAX(output_buffer_size, min_output_buffer_size);
  output_buffer_size = CS_MAX(output_buffer_size, 128);
  output_buffer_size = CS_MIN(output_buffer_size, output_size);

  BFT_MALLOC(output_buffer, output_buffer_size, float);

  /* Loop on dimension (de-interlace vectors, always 3D for EnSight) */

  for (i = 0; i < output_dim; i++) {

    _Bool loop_on_sections = true;

    current_section = export_section;

    while (loop_on_sections == true) {

      while (fvm_writer_field_helper_step_e(helper,
                                            current_section,
                                            input_dim,
                                            i,
                                            interlace,
                                            n_parent_lists,
                                            parent_num_shift,
                                            datatype,
                                            field_values,
                                            output_buffer,
                                            output_buffer_size,
                                            &output_size) == 0) {

        _write_block_floats_l(output_size,
                              output_buffer,
                              f);

      }

      current_section = current_section->next;

      if (   current_section == NULL
          || current_section->continues_previous == false)
        loop_on_sections = false;

    } /* while (loop on sections) */

  } /* end of loop on spatial dimension */

  BFT_FREE(output_buffer);

  return current_section;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to EnSight Gold file writer.
 *
 * Options are:
 *   text                output text files
 *   binary              output binary files (default)
 *   big_endian          force binary files to big-endian
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   divide_polyhedra    tesselate polyhedra with tetrahedra and pyramids
 *                       (adding a vertex near each polyhedron's center)
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque EnSight Gold writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_ensight_init_writer(const char             *name,
                           const char             *path,
                           const char             *options,
                           fvm_writer_time_dep_t   time_dependency,
                           MPI_Comm                comm)
#else
void *
fvm_to_ensight_init_writer(const char             *name,
                           const char             *path,
                           const char             *options,
                           fvm_writer_time_dep_t   time_dependency)
#endif
{
  fvm_to_ensight_writer_t  *this_writer = NULL;

  /* Initialize writer */

  BFT_MALLOC(this_writer, 1, fvm_to_ensight_writer_t);

  BFT_MALLOC(this_writer->name, strlen(name) + 1, char);
  strcpy(this_writer->name, name);

  this_writer->text_mode = false;
  this_writer->swap_endian = false;
  this_writer->discard_polygons = false;
  this_writer->discard_polyhedra = false;
  this_writer->divide_polygons = false;
  this_writer->divide_polyhedra = false;

  this_writer->rank = 0;
  this_writer->n_ranks = 1;

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks, min_rank_step, min_block_size;
    MPI_Comm w_block_comm, w_comm;
    this_writer->min_rank_step = 1;
    this_writer->min_block_size = 1024*1024*8;
    this_writer->block_comm = MPI_COMM_NULL;
    this_writer->comm = MPI_COMM_NULL;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag && comm != MPI_COMM_NULL) {
      this_writer->comm = comm;
      MPI_Comm_rank(this_writer->comm, &rank);
      MPI_Comm_size(this_writer->comm, &n_ranks);
      this_writer->rank = rank;
      this_writer->n_ranks = n_ranks;
      cs_file_get_default_comm(&min_rank_step, &min_block_size,
                               &w_block_comm, &w_comm);
      if (comm == w_comm) {
        this_writer->min_rank_step = min_rank_step;
        this_writer->min_block_size = min_block_size;
        this_writer->block_comm = w_block_comm;
      }
      this_writer->comm = comm;
    }
  }
#endif /* defined(HAVE_MPI) */

  /* Parse options */

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      l_opt = i2 - i1;

      if ((l_opt == 4) && (strncmp(options + i1, "text", l_opt) == 0))
        this_writer->text_mode = true;
      else if ((l_opt == 6) && (strncmp(options + i1, "binary", l_opt) == 0))
        this_writer->text_mode = false;

      else if (   (l_opt == 10)
               && (strncmp(options + i1, "big_endian", l_opt) == 0)) {
        int int_endian = 0;
        this_writer->text_mode = false;
        /* Check if system is "big-endian" or "little-endian" */
        *((char *)(&int_endian)) = '\1';
        if (int_endian == 1)
          this_writer->swap_endian = 1;
      }

      else if (   (l_opt == 16)
               && (strncmp(options + i1, "discard_polygons", l_opt) == 0))
        this_writer->discard_polygons = true;
      else if (   (l_opt == 17)
               && (strncmp(options + i1, "discard_polyhedra", l_opt) == 0))
        this_writer->discard_polyhedra = true;

      else if (   (l_opt == 15)
               && (strncmp(options + i1, "divide_polygons", l_opt) == 0))
        this_writer->divide_polygons = true;
      else if (   (l_opt == 16)
               && (strncmp(options + i1, "divide_polyhedra", l_opt) == 0))
        this_writer->divide_polyhedra = true;

      for (i1 = i2 + 1; i1 < l_tot && options[i1] == ' '; i1++);

    }

  }

  this_writer->case_info = fvm_to_ensight_case_create(name,
                                                      path,
                                                      time_dependency);

  /* Return writer */

  return this_writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to EnSight Gold file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque Ensight Gold writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_ensight_finalize_writer(void  *this_writer_p)
{
  fvm_to_ensight_writer_t  *this_writer
                             = (fvm_to_ensight_writer_t *)this_writer_p;

  BFT_FREE(this_writer->name);

  fvm_to_ensight_case_destroy(this_writer->case_info);

  BFT_FREE(this_writer);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with an EnSight geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_ensight_set_mesh_time(void          *this_writer_p,
                             const int      time_step,
                             const double   time_value)
{
  fvm_to_ensight_writer_t  *this_writer
                             = (fvm_to_ensight_writer_t *)this_writer_p;

  fvm_to_ensight_case_set_geom_time(this_writer->case_info,
                                    time_step,
                                    time_value);
}

/*----------------------------------------------------------------------------
 * Indicate if a elements of a given type in a mesh associated to a given
 * EnSight Gold file writer need to be tesselated.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *   element_type  <-- element type we are interested in
 *
 * returns:
 *   1 if tesselation of the given element type is needed, 0 otherwise
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_needs_tesselation(fvm_writer_t       *this_writer_p,
                                 const fvm_nodal_t  *mesh,
                                 fvm_element_t       element_type)
{
  int  i;
  int  retval = 0;
  fvm_to_ensight_writer_t  *this_writer
                             = (fvm_to_ensight_writer_t *)this_writer_p;

  const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (   (   element_type == FVM_FACE_POLY
          && this_writer->divide_polygons == true)
      || (   element_type == FVM_CELL_POLY
          && this_writer->divide_polyhedra == true)) {

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (section->entity_dim == export_dim) {
        if (section->type == element_type)
          retval = 1;
      }

    }

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a an EnSight Gold file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_ensight_export_nodal(void               *this_writer_p,
                            const fvm_nodal_t  *mesh)
{
  int  part_num;

  const fvm_writer_section_t  *export_section = NULL;
  fvm_writer_section_t        *export_list = NULL;
  fvm_to_ensight_writer_t     *this_writer
                                  = (fvm_to_ensight_writer_t *)this_writer_p;
  _ensight_file_t  f = {NULL, NULL};

  const int  rank = this_writer->rank;
  const int  n_ranks = this_writer->n_ranks;

  /* Initialization */
  /*----------------*/

  fvm_to_ensight_case_file_info_t  file_info;

  /* Get part number */

  part_num = fvm_to_ensight_case_get_part_num(this_writer->case_info,
                                              mesh->name);
  if (part_num == 0)
    part_num = fvm_to_ensight_case_add_part(this_writer->case_info,
                                            mesh->name);

  /* Open geometry file in append mode */

  file_info = fvm_to_ensight_case_get_geom_file(this_writer->case_info);

  f = _open_ensight_file(this_writer,
                         file_info.name,
                         file_info.queried);

  if (file_info.queried == false)
    _write_geom_headers(this_writer, f);

  /* Part header */

  _write_string(f, "part");
  _write_int(f, part_num);
  if (mesh->name != NULL)
    _write_string(f, mesh->name);
  else
    _write_string(f, "unnamed");

  /* Build list of sections that are used here, in order of output */

  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       this_writer->discard_polygons,
                                       this_writer->discard_polyhedra,
                                       this_writer->divide_polygons,
                                       this_writer->divide_polyhedra);

  /* Vertex coordinates */
  /*--------------------*/

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    _export_vertex_coords_g(this_writer, mesh, f);
#endif

  if (n_ranks == 1)
    _export_vertex_coords_l(this_writer, mesh, f);

  /* If no sections are present (i.e. we may only have vertices),
     add  "point" elements */

  if (export_list == NULL) {

#if defined(HAVE_MPI)
    if (n_ranks > 1)
      _export_point_elements_g(this_writer, mesh, f);
#endif
    if (n_ranks == 1)
      _export_point_elements_l(mesh, f);

  }

  /* Element connectivity */
  /*----------------------*/

  export_section = export_list;

  while (export_section != NULL) {

    const fvm_nodal_section_t  *section = export_section->section;

    /* Print header if start of corresponding EnSight section */

    if (export_section->continues_previous == false) {

      cs_gnum_t n_g_elements = 0;
      const fvm_writer_section_t  *next_section = export_section;

      do {

        if (next_section->section->type == export_section->type)
          n_g_elements += fvm_nodal_section_n_g_elements(next_section->section);

        else {
          cs_gnum_t n_g_sub_elements = 0;
          fvm_tesselation_get_global_size(next_section->section->tesselation,
                                          next_section->type,
                                          &n_g_sub_elements,
                                          NULL);
          n_g_elements += n_g_sub_elements;
        }

        next_section = next_section->next;

      } while (next_section != NULL && next_section->continues_previous == true);

      _write_string(f, _ensight_type_name[export_section->type]);
      _write_int(f, n_g_elements);
    }

    /* Output for strided (regular) element types */
    /*--------------------------------------------*/

    if (section->stride > 0) {

#if defined(HAVE_MPI)

      if (n_ranks > 1)
        export_section = _export_nodal_strided_g(this_writer,
                                                 export_section,
                                                 mesh->global_vertex_num,
                                                 f);

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1) { /* start of output in serial mode */

        _write_connect_l(section->stride,
                         section->n_elements,
                         section->vertex_num,
                         f);

        export_section = export_section->next;

      }

    } /* end of output for strided element types */

    /* Output for tesselated polygons or polyhedra */
    /*---------------------------------------------*/

    else if (export_section->type != section->type) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1)
        export_section = _export_nodal_tesselated_g(this_writer,
                                                    export_section,
                                                    mesh->global_vertex_num,
                                                    f);
#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)
        export_section = _export_nodal_tesselated_l(export_section,
                                                    f);

    }

    /* Output for polygons */
    /*---------------------*/

    else if (export_section->type == FVM_FACE_POLY) {
#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1)
        export_section = _export_nodal_polygons_g(this_writer,
                                                  export_section,
                                                  mesh->global_vertex_num,
                                                  f);
#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)
        export_section = _export_nodal_polygons_l(export_section,
                                                  f);

    }

    /* Output for polyhedra */
    /*----------------------*/

    else if (export_section->type == FVM_CELL_POLY) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1)
        export_section =_export_nodal_polyhedra_g(this_writer,
                                                  export_section,
                                                  mesh->global_vertex_num,
                                                  f);

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)
        export_section = _export_nodal_polyhedra_l(export_section,
                                                   f);

    }

  } /* End of loop on sections */

  BFT_FREE(export_list);

  /* Close geometry file and update case file */
  /*------------------------------------------*/

  _free_ensight_file(&f);

  fvm_to_ensight_case_write_case(this_writer->case_info, rank);
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to an EnSight Gold file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
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
fvm_to_ensight_export_field(void                  *this_writer_p,
                            const fvm_nodal_t     *mesh,
                            const char            *name,
                            fvm_writer_var_loc_t   location,
                            int                    dimension,
                            cs_interlace_t         interlace,
                            int                    n_parent_lists,
                            const cs_lnum_t        parent_num_shift[],
                            cs_datatype_t          datatype,
                            int                    time_step,
                            double                 time_value,
                            const void      *const field_values[])
{
  int   output_dim, part_num;
  fvm_to_ensight_case_file_info_t  file_info;

  const fvm_writer_section_t  *export_section = NULL;
  fvm_writer_field_helper_t  *helper = NULL;
  fvm_writer_section_t  *export_list = NULL;
  fvm_to_ensight_writer_t  *this_writer
                             = (fvm_to_ensight_writer_t *)this_writer_p;
  _ensight_file_t  f = {NULL, NULL};

  const int  rank = this_writer->rank;
  const int  n_ranks = this_writer->n_ranks;

  /* Initialization */
  /*----------------*/

  /* Dimension */

  output_dim = dimension;
  if (dimension == 2)
    output_dim = 3;
  else if (dimension > 3 && dimension != 6 && dimension != 9)
    bft_error(__FILE__, __LINE__, 0,
              _("Data of dimension %d not handled"), dimension);

  /* Get part number */

  part_num = fvm_to_ensight_case_get_part_num(this_writer->case_info,
                                              mesh->name);
  if (part_num == 0)
    part_num = fvm_to_ensight_case_add_part(this_writer->case_info,
                                            mesh->name);

  /* Open variable file */

  file_info = fvm_to_ensight_case_get_var_file(this_writer->case_info,
                                               name,
                                               output_dim,
                                               location,
                                               time_step,
                                               time_value);

  f = _open_ensight_file(this_writer, file_info.name, file_info.queried);

  if (file_info.queried == false) {

    char buf[81] = "";

    /* New files start with description line */
#if HAVE_SNPRINTF
    if (time_step > -1)
      snprintf(buf, 80, "%s (time values: %d, %g)",
               name, time_step, time_value);
    else
      strncpy(buf, name, 80);
#else
    strncpy(buf, name, 80);
#endif
    buf[80] = '\0';
    _write_string(f, buf);
  }

  /* Initialize writer helper */
  /*--------------------------*/

  /* Build list of sections that are used here, in order of output */

  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       this_writer->discard_polygons,
                                       this_writer->discard_polyhedra,
                                       this_writer->divide_polygons,
                                       this_writer->divide_polyhedra);

  if (n_ranks == 1)
    helper = fvm_writer_field_helper_create(mesh,
                                            export_list,
                                            output_dim,
                                            CS_NO_INTERLACE,
                                            CS_FLOAT,
                                            location);

  /* Part header */

  _write_string(f, "part");
  _write_int(f, part_num);

  /* Per node variable */
  /*-------------------*/

  if (location == FVM_WRITER_PER_NODE) {

    _write_string(f, "coordinates");

#if defined(HAVE_MPI)

    if (n_ranks > 1)
      _export_field_values_ng(this_writer,
                              mesh,
                              dimension,
                              output_dim,
                              interlace,
                              n_parent_lists,
                              parent_num_shift,
                              datatype,
                              field_values,
                              f);

#endif /* defined(HAVE_MPI) */

    if (n_ranks == 1)
      _export_field_values_nl(mesh,
                              helper,
                              dimension,
                              interlace,
                              n_parent_lists,
                              parent_num_shift,
                              datatype,
                              field_values,
                              f);
  }

  /* Per element variable */
  /*----------------------*/

  else if (location == FVM_WRITER_PER_ELEMENT) {

    export_section = export_list;

    while (export_section != NULL) {

      /* Print header if start of corresponding EnSight section */

      if (export_section->continues_previous == false)
        _write_string(f, _ensight_type_name[export_section->type]);

      /* Output per grouped sections */

#if defined(HAVE_MPI)

      if (n_ranks > 1)
        export_section = _export_field_values_eg(this_writer,
                                                 export_section,
                                                 dimension,
                                                 output_dim,
                                                 interlace,
                                                 n_parent_lists,
                                                 parent_num_shift,
                                                 datatype,
                                                 field_values,
                                                 f);

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)
        export_section = _export_field_values_el(export_section,
                                                 helper,
                                                 dimension,
                                                 interlace,
                                                 n_parent_lists,
                                                 parent_num_shift,
                                                 datatype,
                                                 field_values,
                                                 f);

    } /* End of loop on sections */

  } /* End for per element variable */

  /* Free helper structures */
  /*------------------------*/

  if (helper != NULL)
    helper = fvm_writer_field_helper_destroy(helper);

  BFT_FREE(export_list);

  /* Close variable file and update case file */
  /*------------------------------------------*/

  _free_ensight_file(&f);

  fvm_to_ensight_case_write_case(this_writer->case_info, rank);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
