/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to EnSight Gold files
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2007  EDF

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_file.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_config_defs.h"
#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_gather.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_parall.h"
#include "fvm_to_ensight_case.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_ensight_v1.h"

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

/*----------------------------------------------------------------------------
 * EnSight Gold writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;             /* Writer name */

  int          rank;             /* Rank of current process in communicator */
  int          n_ranks;          /* Number of processes in communicator */

  _Bool        text_mode;        /* true if using text output */
  _Bool        swap_endian;      /* true if binary file endianness must
                                    be changed */

  _Bool        discard_polygons;   /* Option to discard polygonal elements */
  _Bool        discard_polyhedra;  /* Option to discard polyhedral elements */

  _Bool        divide_polygons;    /* Option to tesselate polygonal elements */
  _Bool        divide_polyhedra;   /* Option to tesselate polyhedral elements */

  fvm_to_ensight_case_t  *case_info;  /* Associated case structure */

#if defined(HAVE_MPI)
  MPI_Comm     comm;             /* Associated MPI communicator */
#endif

} fvm_to_ensight_writer_t;

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
 * Count number of extra vertices when tesselations are present
 *
 * parameters:
 *   this_writer         <-- pointer to associated writer
 *   mesh                <-- pointer to nodal mesh structure
 *   n_extra_vertices_g  --> global number of extra vertices (optional)
 *   n_extra_vertices    --> local number of extra vertices (optional)
 *----------------------------------------------------------------------------*/

static void
_count_extra_vertices(const fvm_to_ensight_writer_t  *this_writer,
                      const fvm_nodal_t              *mesh,
                      fvm_gnum_t                     *n_extra_vertices_g,
                      fvm_lnum_t                     *n_extra_vertices)
{
  int  i;

  const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (n_extra_vertices_g != NULL)
    *n_extra_vertices_g = 0;
  if (n_extra_vertices != NULL)
    *n_extra_vertices   = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvm_nodal_section_t  *const  section = mesh->sections[i];

    /* Output if entity dimension equal to highest in mesh
       (i.e. no output of faces if cells present, or edges
       if cells or faces) */

    if (   section->entity_dim == export_dim
        && section->type == FVM_CELL_POLY
        && section->tesselation != NULL
        && this_writer->divide_polyhedra == true) {

      if (n_extra_vertices_g != NULL)
        *n_extra_vertices_g
          += fvm_tesselation_n_g_vertices_add(section->tesselation);

      if (n_extra_vertices != NULL)
        *n_extra_vertices
          += fvm_tesselation_n_vertices_add(section->tesselation);

    }

  }

}

/*----------------------------------------------------------------------------
 * Write string to a text or C binary EnSight Gold file
 *
 * parameters:
 *   f                <-- file to write to
 *   s                <-- string to write
 *----------------------------------------------------------------------------*/

static void
_write_string(bft_file_t  *f,
              const char  *s)
{
  size_t  i;
  char  buf[82];
  bft_file_type_t  type;

  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  type = bft_file_get_type(f);

  if (type == BFT_FILE_TYPE_TEXT) {
    strncpy(buf, s, 80);
    buf[80] = '\0';
    bft_file_printf(f, "%s\n", buf);
  }

  else {
    strncpy(buf, s, 80);
    buf[80] = '\0';
    for (i = strlen(buf) ; i < 80 ; i++)
      buf[i] = ' ';
    bft_file_write(buf, 1, 80, f);
  }

}

/*----------------------------------------------------------------------------
 * Write integer to a text or C binary EnSight Gold file
 *
 * parameters:
 *   f                <-- file to write to
 *   n                <-- integer value to write
 *----------------------------------------------------------------------------*/

inline static void
_write_int(bft_file_t  *f,
           const int    n)
{
  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  if (bft_file_get_type(f) == BFT_FILE_TYPE_TEXT)
    bft_file_printf(f, "%10d\n", n);

  else
    bft_file_write(&n, sizeof(int), 1, f);
}

/*----------------------------------------------------------------------------
 * Initialize FVM to EnSight Gold geometry file.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque Ensight Gold writer structure.
 *----------------------------------------------------------------------------*/

static void
_init_geom_file(fvm_to_ensight_writer_t  *this_writer)
{
  if (this_writer->rank == 0) {

    bft_file_t *f;
    bft_file_type_t file_type;
    fvm_to_ensight_case_file_info_t  file_info;

    if (this_writer->text_mode == true)
      file_type = BFT_FILE_TYPE_TEXT;
    else
      file_type = BFT_FILE_TYPE_BINARY;

    /* Write file header (future calls will be in append mode) */

    file_info = fvm_to_ensight_case_get_geom_file(this_writer->case_info);

    f = bft_file_open(file_info.name,
                      BFT_FILE_MODE_WRITE,
                      file_type);

    if (this_writer->swap_endian == true)
      bft_file_set_swap_endian(f, 1);

    if (this_writer->text_mode == false)
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

    f = bft_file_free(f);
  }

}

/*----------------------------------------------------------------------------
 * Write slice of a vector of floats to an EnSight Gold file
 *
 * parameters:
 *   num_start      <-- global number of first element for this slice
 *   num_end        <-- global number of past the last element for this slice
 *   values         <-- pointer to values slice array
 *   f              <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_slice_values(fvm_gnum_t    num_start,
                    fvm_gnum_t    num_end,
                    const float   values[],
                    bft_file_t   *f)
{
  size_t  i;
  fvm_gnum_t  j;

  /* If called by non I/O rank or no values, return */
  if (f == NULL || (num_start > num_end))
    return;

  if (bft_file_get_type(f) == BFT_FILE_TYPE_TEXT) {
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12.5e\n", values[i]);
  }
  else
    bft_file_write(values, sizeof(float), num_end - num_start, f);
}

/*----------------------------------------------------------------------------
 * Return extra vertex coordinates when tesselations are present
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *
 * returns:
 *   array containing all extra vertex coordinates
 *----------------------------------------------------------------------------*/

static fvm_coord_t *
_extra_vertex_coords(const fvm_to_ensight_writer_t  *this_writer,
                     const fvm_nodal_t              *mesh)
{
  int  i;
  fvm_lnum_t  n_extra_vertices_section;

  fvm_lnum_t  n_extra_vertices = 0;
  size_t  coord_shift = 0;
  fvm_coord_t  *coords = NULL;

  _count_extra_vertices(this_writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  if (n_extra_vertices > 0) {

    const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

    BFT_MALLOC(coords, n_extra_vertices * 3, fvm_coord_t);

    for (i = 0 ; i < mesh->n_sections ; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == export_dim
          && section->type == FVM_CELL_POLY
          && section->tesselation != NULL
          && this_writer->divide_polyhedra == true) {

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
 * Write vertex coordinates to an EnSight Gold file in serial mode
 *
 * parameters:
 *   this_writer   <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure
 *   comm          <-- associated MPI communicator
 *   global_s_size <-- global slice size
 *   f             <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_g(const fvm_to_ensight_writer_t  *this_writer,
                        const fvm_nodal_t              *mesh,
                        MPI_Comm                        comm,
                        fvm_gnum_t                      global_s_size,
                        bft_file_t                     *f)
{
  size_t      stride;
  fvm_lnum_t  i, j;

  fvm_gnum_t   global_num_start;
  fvm_gnum_t   global_num_end;

  int          rank;
  int          section_id;

  fvm_gnum_t   n_g_extra_vertices = 0;
  fvm_lnum_t   n_extra_vertices = 0;
  fvm_coord_t  *extra_vertex_coords = NULL;
  float        *coords_tmp = NULL;
  float        *global_coords_s = NULL;

  fvm_gather_slice_t  *vertices_slice = NULL;

  const double      *vertex_coords = mesh->vertex_coords;
  const fvm_lnum_t  *parent_vertex_num = mesh->parent_vertex_num;
  const fvm_lnum_t  n_vertices
                      = fvm_io_num_get_local_count(mesh->global_vertex_num);
  const fvm_gnum_t  n_g_vertices
                      = fvm_io_num_get_global_count(mesh->global_vertex_num);

  /* Get info on the current MPI communicator */

  MPI_Comm_rank(comm, &rank);

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(this_writer,
                        mesh,
                        &n_g_extra_vertices,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(this_writer,
                                             mesh);

  /* Vertex coordinates */
  /*--------------------*/

  stride = (size_t)(mesh->dim);

  BFT_MALLOC(coords_tmp, FVM_MAX(n_vertices, n_extra_vertices), float);
  BFT_MALLOC(global_coords_s, global_s_size, float);

  if (rank == 0) {
    _write_string(f, "coordinates");
    _write_int(f, (int)(n_g_vertices + n_g_extra_vertices));
  }

  vertices_slice = fvm_gather_slice_create(mesh->global_vertex_num,
                                           global_s_size,
                                           comm);

  /* Loop on dimension (de-interlace coordinates, always 3D for EnSight) */

  for (j = 0 ; j < 3 ; j++) {

    if (j < mesh->dim) {
      if (parent_vertex_num != NULL) {
        for (i = 0 ; i < n_vertices ; i++)
          coords_tmp[i]
            = (float)(vertex_coords[(parent_vertex_num[i]-1)*stride + j]);
      }
      else {
        for (i = 0 ; i < n_vertices ; i++)
          coords_tmp[i] = (float)(vertex_coords[i*stride + j]);
      }
    }
    else {
      for (i = 0 ; i < n_vertices ; i++)
        coords_tmp[i] = 0.0;
    }

    /* loop on slices in parallel mode */

    if (j > 0)
      fvm_gather_slice_reinitialize(vertices_slice);

    while (fvm_gather_slice_advance(vertices_slice,
                                    &global_num_start,
                                    &global_num_end) == 0) {

      fvm_gather_array(coords_tmp,
                       global_coords_s,
                       MPI_FLOAT,
                       1,
                       mesh->global_vertex_num,
                       comm,
                       vertices_slice);

      if (rank == 0)
        _write_slice_values(global_num_start,
                            global_num_end,
                            global_coords_s,
                            f);

    }

    /* Now handle extra vertices */
    /*---------------------------*/

    if (n_g_extra_vertices > 0) {

      fvm_lnum_t extra_vertices_count = 0;

      for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

        fvm_gather_slice_t  *extra_vertices_slice = NULL;
        const fvm_nodal_section_t  *const  section = mesh->sections[section_id];

        /* Output if entity dimension equal to highest in mesh
           (i.e. no output of faces if cells present, or edges
           if cells or faces) */

        if (   section->entity_dim == mesh->dim
            && section->type == FVM_CELL_POLY
            && section->tesselation != NULL
            && this_writer->divide_polyhedra == true) {

          const fvm_io_num_t *extra_vertex_num
            = fvm_tesselation_global_vertex_num(section->tesselation);
          const fvm_lnum_t n_extra_vertices_section
            = fvm_tesselation_n_vertices_add(section->tesselation);

          for (i = 0 ; i < n_extra_vertices ; i++)
            coords_tmp[i]
              = (float)(extra_vertex_coords[  (i+extra_vertices_count)*stride
                                            + j]);

          extra_vertices_slice = fvm_gather_slice_create(extra_vertex_num,
                                                         global_s_size,
                                                         comm);

          /* loop on slices in parallel mode */

          while (fvm_gather_slice_advance(extra_vertices_slice,
                                          &global_num_start,
                                          &global_num_end) == 0) {

            fvm_gather_array(coords_tmp,
                             global_coords_s,
                             MPI_FLOAT,
                             1,
                             extra_vertex_num,
                             comm,
                             extra_vertices_slice);

            if (rank == 0)
              _write_slice_values(global_num_start,
                                  global_num_end,
                                  global_coords_s,
                                  f);

          }

          fvm_gather_slice_destroy(extra_vertices_slice);

          extra_vertices_count += n_extra_vertices_section;

        }

      } /* end of loop on sections for extra vertices */

    } /* end handling for extra vertices */

  } /* end of loop on spatial dimension */

  fvm_gather_slice_destroy(vertices_slice);

  BFT_FREE(global_coords_s);
  BFT_FREE(coords_tmp);

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
 *   f           <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_l(const fvm_to_ensight_writer_t  *this_writer,
                        const fvm_nodal_t              *mesh,
                        bft_file_t                     *f)
{
  fvm_lnum_t   i, j;
  fvm_lnum_t   n_extra_vertices = 0;
  fvm_coord_t  *extra_vertex_coords = NULL;
  float        *coords_tmp = NULL;

  const fvm_lnum_t   n_vertices = mesh->n_vertices;
  const double      *vertex_coords = mesh->vertex_coords;
  const fvm_lnum_t  *parent_vertex_num = mesh->parent_vertex_num;

  const size_t  stride = (size_t)(mesh->dim);

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(this_writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(this_writer,
                                             mesh);

  /* Vertex coordinates */
  /*--------------------*/

  BFT_MALLOC(coords_tmp, FVM_MAX(n_vertices, n_extra_vertices), float);

  _write_string(f, "coordinates");
  _write_int(f, n_vertices + n_extra_vertices);

  /* Loop on dimension (de-interlace coordinates, always 3D for EnSight) */

  for (j = 0 ; j < 3 ; j++) {

    /* First, handle regular vertices */

    if (j < mesh->dim) {
      if (parent_vertex_num != NULL) {
        for (i = 0 ; i < n_vertices ; i++)
          coords_tmp[i]
            = (float)(vertex_coords[(parent_vertex_num[i]-1)*stride + j]);
      }
      else {
        for (i = 0 ; i < n_vertices ; i++)
          coords_tmp[i] = (float)(vertex_coords[i*stride + j]);
      }
    }
    else {
      for (i = 0 ; i < (n_vertices) ; i++)
        coords_tmp[i] = 0.0;
    }

    _write_slice_values(1,
                        (fvm_gnum_t)(n_vertices + 1),
                        coords_tmp,
                        f);

    /* Handle extra vertices (only occur with polyhedra tesselations in 3d) */

    for (i = 0 ; i < n_extra_vertices ; i++)
      coords_tmp[i] = (float)(extra_vertex_coords[i*stride + j]);

    if (n_extra_vertices > 0)
      _write_slice_values(1,
                          (fvm_gnum_t)(n_extra_vertices + 1),
                          coords_tmp,
                          f);

  } /* end of loop on mesh dimension */

  BFT_FREE(coords_tmp);

  if (extra_vertex_coords != NULL)
    BFT_FREE(extra_vertex_coords);
}

/*----------------------------------------------------------------------------
 * Write "trivial" point elements to an EnSight Gold file
 *
 * parameters:
 *   n_g_vertices     <-- number of vertices
 *   buffer_size      <-- size of write buffer
 *   buffer           --- write buffer (for binary mode)
 *   f                <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_export_point_elements_l(fvm_gnum_t   n_g_vertices,
                         size_t       buffer_size,
                         int          buffer[],
                         bft_file_t  *f)
{
  fvm_gnum_t  i;
  int  j = 1;

  _write_string(f, "point");
  _write_int(f, (int)n_g_vertices);

  if (bft_file_get_type(f) == BFT_FILE_TYPE_TEXT) { /* Text mode */

    for (i = 0 ; i < n_g_vertices ; i++)
      bft_file_printf(f, "%10d\n", j++);

  }
  else { /* Binary mode */

    size_t  k;
    int  j_end = n_g_vertices + 1;

    while (j < j_end) {
      for (k = 0 ;  j < j_end  && k < buffer_size ; k++)
        buffer[k] = j++;
      bft_file_write(buffer, sizeof(int), k, f);
    }

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided global connectivity slice to an EnSight Gold file
 *
 * parameters:
 *   stride           <-- number of vertices per element type
 *   num_start        <-- global number of first element for this slice
 *   num_end          <-- global number of past last element for this slice
 *   global_connect_s <-- global connectivity slice array
 *   buffer_size      <-- size of write buffer
 *   buffer           --- write buffer (for binary mode)
 *   f                <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_slice_connect_g(int                stride,
                       fvm_gnum_t         num_start,
                       fvm_gnum_t         num_end,
                       const fvm_gnum_t   global_connect_s[],
                       size_t             buffer_size,
                       int                buffer[],
                       bft_file_t        *f)
{
  size_t  i;
  fvm_gnum_t  j;

  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  if (bft_file_get_type(f) == BFT_FILE_TYPE_TEXT) { /* Text mode */

    switch(stride) {

    case 2: /* edge */
      for (i = 0, j = num_start ; j < num_end ; i++, j++)
        bft_file_printf(f, "%10d%10d\n",
                        (int)global_connect_s[i*2],
                        (int)global_connect_s[i*2+1]);
      break;

    case 3: /* FVM_FACE_TRIA */
      for (i = 0, j = num_start ; j < num_end ; i++, j++)
        bft_file_printf(f, "%10d%10d%10d\n",
                        (int)global_connect_s[i*3],
                        (int)global_connect_s[i*3+1],
                        (int)global_connect_s[i*3+2]);
      break;

    case 4: /* FVM_FACE_QUAD or FVM_CELL_TETRA */
      for (i = 0, j = num_start ; j < num_end ; i++, j++)
        bft_file_printf(f, "%10d%10d%10d%10d\n",
                        (int)global_connect_s[i*4],
                        (int)global_connect_s[i*4+1],
                        (int)global_connect_s[i*4+2],
                        (int)global_connect_s[i*4+3]);
      break;

    case 5: /* FVM_CELL_PYRAM */
      for (i = 0, j = num_start ; j < num_end ; i++, j++)
        bft_file_printf(f, "%10d%10d%10d%10d%10d\n",
                        (int)global_connect_s[i*5],
                        (int)global_connect_s[i*5+1],
                        (int)global_connect_s[i*5+2],
                        (int)global_connect_s[i*5+3],
                        (int)global_connect_s[i*5+4]);
      break;

    case 6: /* FVM_CELL_PRISM */
      for (i = 0, j = num_start ; j < num_end ; i++, j++)
        bft_file_printf(f, "%10d%10d%10d%10d%10d%10d\n",
                        (int)global_connect_s[i*6],
                        (int)global_connect_s[i*6+1],
                        (int)global_connect_s[i*6+2],
                        (int)global_connect_s[i*6+3],
                        (int)global_connect_s[i*6+4],
                        (int)global_connect_s[i*6+5]);
      break;

    case 8: /* FVM_CELL_HEXA */
      for (i = 0, j = num_start ; j < num_end ; i++, j++)
        bft_file_printf(f,
                        "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                        (int)global_connect_s[i*8],
                        (int)global_connect_s[i*8+1],
                        (int)global_connect_s[i*8+2],
                        (int)global_connect_s[i*8+3],
                        (int)global_connect_s[i*8+4],
                        (int)global_connect_s[i*8+5],
                        (int)global_connect_s[i*8+6],
                        (int)global_connect_s[i*8+7]);
      break;

    default:
      assert(0);
    }

  }
  else { /* Binary mode */

    size_t  k = 0;
    size_t  n_values = (num_end - num_start)*stride;

    while (k < n_values) {
      for (i = 0 ; i < buffer_size && k < n_values ; i++)
        buffer[i] = (int)(global_connect_s[k++]);
      bft_file_write(buffer, sizeof(int), i, f);
    }

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
 *   buffer_size <-- size of write buffer
 *   buffer      --- write buffer (for binary mode)
 *   f           <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_connect_l(int                stride,
                 fvm_lnum_t         n_elems,
                 const fvm_lnum_t   connect[],
                 size_t             buffer_size,
                 int                buffer[],
                 bft_file_t        *f)
{
  fvm_lnum_t  i;

  if (bft_file_get_type(f) == BFT_FILE_TYPE_TEXT) { /* Text mode */

    switch(stride) {

    case 2: /* edge */
      for (i = 0 ; i < n_elems ; i++)
        bft_file_printf(f, "%10d%10d\n",
                        (int)connect[i*2],
                        (int)connect[i*2+1]);
      break;

    case 3: /* FVM_FACE_TRIA */
      for (i = 0 ; i < n_elems ; i++)
        bft_file_printf(f, "%10d%10d%10d\n",
                        (int)connect[i*3],
                        (int)connect[i*3+1],
                        (int)connect[i*3+2]);
      break;

    case 4: /* FVM_FACE_QUAD or FVM_CELL_TETRA */
      for (i = 0 ; i < n_elems ; i++)
        bft_file_printf(f, "%10d%10d%10d%10d\n",
                        (int)connect[i*4],
                        (int)connect[i*4+1],
                        (int)connect[i*4+2],
                        (int)connect[i*4+3]);
      break;

    case 5: /* FVM_CELL_PYRAM */
      for (i = 0 ; i < n_elems ; i++)
        bft_file_printf(f, "%10d%10d%10d%10d%10d\n",
                        (int)connect[i*5],
                        (int)connect[i*5+1],
                        (int)connect[i*5+2],
                        (int)connect[i*5+3],
                        (int)connect[i*5+4]);
      break;

    case 6: /* FVM_CELL_PRISM */
      for (i = 0 ; i < n_elems ; i++)
        bft_file_printf(f, "%10d%10d%10d%10d%10d%10d\n",
                        (int)connect[i*6],
                        (int)connect[i*6+1],
                        (int)connect[i*6+2],
                        (int)connect[i*6+3],
                        (int)connect[i*6+4],
                        (int)connect[i*6+5]);
      break;

    case 8: /* FVM_CELL_HEXA */
      for (i = 0 ; i < n_elems ; i++)
        bft_file_printf(f,
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
  else { /* Binary mode */

    size_t  j;
    size_t  k = 0;
    size_t  n_values = n_elems*stride;

    while (k < n_values) {
      for (j = 0 ; j < buffer_size && k < n_values  ; j++)
        buffer[j] = (int)(connect[k++]);
      bft_file_write(buffer, sizeof(int), j, f);
    }

  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write indexed element lengths from a nodal mesh to an EnSight Gold
 * file in parallel mode
 *
 * parameters:
 *   global_element_num           <-- global element numbering
 *   vertex_index                 <-- pointer to element -> vertex index
 *   comm                         <-- associated MPI communicator
 *   n_ranks                      <-- number of processes in communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller      --- global connectivity slice provided
 *                                    by caller
 *   text_mode                    <-- true if text output, false if binary
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_write_lengths_g(const fvm_io_num_t  *global_element_num,
                 const fvm_lnum_t     vertex_index[],
                 MPI_Comm             comm,
                 fvm_gnum_t           global_s_size,
                 fvm_gnum_t           global_connect_s_size_caller,
                 fvm_gnum_t           global_connect_s_caller[],
                 _Bool                text_mode,
                 size_t               buffer_size,
                 int                  buffer[],
                 bft_file_t          *f)
{
  fvm_lnum_t  i;
  fvm_gnum_t  i_s;
  size_t      i_buf;

  int         rank;

  fvm_gnum_t  global_num_start;
  fvm_gnum_t  global_num_end;

  int        *lengths = NULL;
  int        *lengths_s = NULL;

  fvm_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvm_gnum_t  *global_connect_s = global_connect_s_caller;

  fvm_gather_slice_t   *elements_slice = NULL;

  const fvm_lnum_t  n_elements = fvm_io_num_get_local_count(global_element_num);

  /* Get info on the current MPI communicator */

  MPI_Comm_rank(comm, &rank);

  /* Build local lengths */

  BFT_MALLOC(lengths, n_elements, fvm_lnum_t);

  for (i = 0 ; i < n_elements ; i++)
    lengths[i] = vertex_index[i+1] - vertex_index[i];

  /* Use global_connect_s memory area for lengths_s; if
     necessary, increase its size */

  if (  sizeof(fvm_lnum_t)*global_s_size
      > sizeof(fvm_gnum_t)*global_connect_s_size) {
    global_connect_s_size =   (global_s_size*sizeof(fvm_lnum_t))
                            / sizeof(fvm_gnum_t);
    if (global_connect_s == global_connect_s_caller)
      global_connect_s = NULL;
    BFT_REALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);
  }
  lengths_s = (int *)global_connect_s;

  /* Loop on slices */
  /*----------------*/

  elements_slice = fvm_gather_slice_create(global_element_num,
                                           global_s_size,
                                           comm);

  while (fvm_gather_slice_advance(elements_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather number of vertices per element */

    fvm_gather_array(lengths,
                     lengths_s,
                     MPI_INT,
                     1,
                     global_element_num,
                     comm,
                     elements_slice);

    /* Do all printing for rank 0 */

    if (rank == 0) {

      /* Print number of vertices per element */

      if (text_mode == true) {
        for (i = 0, i_s = global_num_start ;
             i_s < global_num_end ;
             i++, i_s++)
          bft_file_printf(f, "%10d\n", lengths_s[i]);
      }
      else {
        for (i = 0, i_s = global_num_start, i_buf = 0 ;
             i_s < global_num_end ;
             i++, i_s++) {
          if (i_buf == buffer_size) {
            bft_file_write(buffer, sizeof(int), i_buf, f);
            i_buf = 0;
          }
          buffer[i_buf++] = (int)(lengths_s[i]);
        }
        if (i_buf > 0)
          bft_file_write(buffer, sizeof(int), i_buf, f);
      }

    }

  } /* end of loop on slices */

  fvm_gather_slice_destroy(elements_slice);

  BFT_FREE(lengths);

  if (global_connect_s != global_connect_s_caller)
    BFT_FREE(global_connect_s);

}

/*----------------------------------------------------------------------------
 * Write indexed element (polygons or polyhedra) cell -> vertex connectivity
 * to an EnSight Gold file in parallel mode.
 *
 * In text mode, zeroes may be used in place of extra vertex numbers
 * to indicate extra newlines (so as to place newline characters between
 * face -> vertex definitions for a polyhedral cell->vertex connectivity).
 *
 * parameters:
 *   global_vertex_num            <-- vertex global numbering
 *   global_element_num           <-- global element numbering
 *   vertex_index                 <-- element -> vertex index
 *   vertex_num                   <-- element -> vertex number
 *   comm                         <-- associated MPI communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller      --- global connectivity slice provided
 *                                    by caller
 *   text_mode                    <-- true if text output, false if binary
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static void
_write_indexed_connect_g(const fvm_io_num_t  *global_vertex_num,
                         const fvm_io_num_t  *global_element_num,
                         const fvm_lnum_t     vertex_index[],
                         const fvm_lnum_t     vertex_num[],
                         MPI_Comm             comm,
                         fvm_gnum_t           global_s_size,
                         fvm_gnum_t           global_connect_s_size_caller,
                         fvm_gnum_t           global_connect_s_caller[],
                         _Bool                text_mode,
                         size_t               buffer_size,
                         int                  buffer[],
                         bft_file_t          *f)
{
  fvm_lnum_t  i, j;
  fvm_gnum_t  i_s, j_s, k_s;
  size_t      i_buf;

  int         rank;

  fvm_gnum_t  global_num_start;
  fvm_gnum_t  global_num_end;

  fvm_gnum_t  *global_idx_s = NULL;

  fvm_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvm_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvm_gnum_t  *global_connect_s = global_connect_s_caller;

  fvm_gather_slice_t   *elements_slice = NULL;

  /* Get info on the current MPI communicator */

  MPI_Comm_rank(comm, &rank);

  /* Allocate memory for additionnal indexes */

  BFT_MALLOC(global_idx_s, global_s_size + 1, fvm_gnum_t);

  /* Loop on slices */
  /*----------------*/

  elements_slice = fvm_gather_slice_create(global_element_num,
                                           global_s_size,
                                           comm);

  while (fvm_gather_slice_advance(elements_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather element->vertices index */

    fvm_gather_slice_index(vertex_index,
                           global_idx_s,
                           global_element_num,
                           comm,
                           elements_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvm_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    comm,
                                    global_idx_s,
                                    elements_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_caller)
        global_connect_s = NULL;
      BFT_REALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now gather element->vertices connectivity */

    fvm_gather_indexed_numbers(vertex_index,
                               vertex_num,
                               global_connect_s,
                               global_vertex_num,
                               global_element_num,
                               comm,
                               global_idx_s,
                               elements_slice);

    /* Do all printing for cells on rank 0 */

    if (rank == 0) {

      if (text_mode == true) {

        /* Print face connectivity */

        k_s = 0;

        /* Loop on elements in slice */

        for (i = 0, i_s = global_num_start ;
             i_s < global_num_end ;
             i++, i_s++) {

          /* Print element vertex numbers */

          for (j = 0, j_s = global_idx_s[i] ;
               j_s < global_idx_s[i+1] ;
               j++, j_s++) {
            if (global_connect_s[k_s] != 0)
              bft_file_printf(f, "%10d",
                              (int)global_connect_s[k_s++]);
            else {
              if (j_s < global_idx_s[i+1] - 1)
                bft_file_printf(f, "\n");
              k_s++;
            }
          }

          bft_file_printf(f, "\n");

          assert(k_s == global_idx_s[i+1]);

        } /* End of loop elements in slice */

      } /* End of text file specific treatement */

      else { /* text_mode = false */

        size_t  i_start = global_idx_s[0];
        size_t  i_end   = global_idx_s[global_num_end - global_num_start];

        for (i_s = i_start, i_buf = 0 ; i_s < i_end ; i_s++) {
          if (i_buf == buffer_size) {
            bft_file_write(buffer, sizeof(int), i_buf, f);
            i_buf = 0;
          }
          buffer[i_buf++] = (int)(global_connect_s[i_s]);
        }
        if (i_buf > 0)
          bft_file_write(buffer, sizeof(int), i_buf, f);

      } /* End of binary-file specific treatement */

    } /* End of printing for rank 0 for this slice */

  }

  /* Free memory */

  elements_slice = fvm_gather_slice_destroy(elements_slice);

  BFT_FREE(global_idx_s);

  if (global_connect_s != global_connect_s_caller)
    BFT_FREE(global_connect_s);
}

/*----------------------------------------------------------------------------
 * Write tesselated element cell -> vertex connectivity to an EnSight Gold
 * file in parallel mode.
 *
 * parameters:
 *   global_vertex_num            <-- vertex global numbering
 *   global_element_num           <-- global element numbering
 *   tesselation                  <-- element tesselation description
 *   type                         <-- tesselated sub-element type
 *   extra_vertex_base            <-- starting number for added vertices
 *   comm                         <-- associated MPI communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller      --- global connectivity slice provided
 *                                    by caller
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static void
_write_tesselated_connect_g(const fvm_io_num_t  *global_vertex_num,
                            const fvm_io_num_t  *global_element_num,
                            const fvm_tesselation_t  *tesselation,
                            fvm_element_t        type,
                            const fvm_gnum_t     extra_vertex_base,
                            MPI_Comm             comm,
                            fvm_gnum_t           global_s_size,
                            fvm_gnum_t           global_connect_s_size_caller,
                            fvm_gnum_t           global_connect_s_caller[],
                            size_t               buffer_size,
                            int                  buffer[],
                            bft_file_t          *f)
{
  int         rank;

  fvm_lnum_t  n_sub_elements_max, local_connect_size;
  fvm_gnum_t  global_num_start, global_num_end;

  fvm_lnum_t  start_id = 0;
  fvm_lnum_t  end_id = 0;

  fvm_lnum_t  *local_idx = NULL;
  fvm_gnum_t  *global_idx_s = NULL;
  fvm_gnum_t  *sub_elt_vertex_num = NULL;

  fvm_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvm_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvm_gnum_t  *global_connect_s = global_connect_s_caller;

  fvm_gather_slice_t   *elements_slice = NULL;

  const fvm_lnum_t n_elements = fvm_tesselation_n_elements(tesselation);
  const int stride = fvm_nodal_n_vertices_element[type];

  /* Get info on the current MPI communicator */

  MPI_Comm_rank(comm, &rank);

  fvm_tesselation_get_global_size(tesselation,
                                  type,
                                  NULL,
                                  &n_sub_elements_max);

  /* Allocate memory for additionnal indexes and decoded connectivity */

  BFT_MALLOC(local_idx, n_elements + 1, fvm_lnum_t);
  BFT_MALLOC(global_idx_s, global_s_size + 1, fvm_gnum_t);

  local_connect_size = FVM_MAX(global_s_size,
                               (fvm_gnum_t)n_sub_elements_max*10);
  BFT_MALLOC(sub_elt_vertex_num, local_connect_size * stride, fvm_gnum_t);

  /* Loop on slices */
  /*----------------*/

  elements_slice = fvm_gather_slice_create(global_element_num,
                                           global_s_size,
                                           comm);

  while (fvm_gather_slice_advance(elements_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Build element->vertices index */

    end_id
      = fvm_tesselation_range_index_g(tesselation,
                                      type,
                                      fvm_nodal_n_vertices_element[type],
                                      start_id,
                                      local_connect_size,
                                      &global_num_end,
                                      local_idx,
                                      comm);

    /* Check if the maximum id returned on some ranks leads to a
       lower global_num_end than initially required (due to the
       local buffer being too small) and adjust slice if necessary */

    fvm_gather_slice_limit(elements_slice, &global_num_end);

    /* Gather element->vertices index */

    fvm_gather_slice_index(local_idx,
                           global_idx_s,
                           global_element_num,
                           comm,
                           elements_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvm_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    comm,
                                    global_idx_s,
                                    elements_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_caller)
        global_connect_s = NULL;
      BFT_REALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now decode tesselation */

    end_id = fvm_tesselation_decode_g(tesselation,
                                      type,
                                      start_id,
                                      local_connect_size,
                                      &global_num_end,
                                      global_vertex_num,
                                      extra_vertex_base,
                                      sub_elt_vertex_num,
                                      comm);

    /* No need to check if the maximum id returned on some ranks
       leads to a lower global_num_end than initially required
       (due to local buffer being full), as this was already done
       above for the local index */

    /* Now gather decoded element->vertices connectivity */

    fvm_gather_indexed(sub_elt_vertex_num,
                       global_connect_s,
                       FVM_MPI_GNUM,
                       local_idx,
                       global_element_num,
                       comm,
                       global_idx_s,
                       elements_slice);

    /* Do all printing for cells on rank 0 */

    if (rank == 0) {
      fvm_gnum_t num_start = global_idx_s[0] / stride;
      fvm_gnum_t num_end = global_idx_s[(  global_num_end
                                         - global_num_start)] / stride;
      _write_slice_connect_g(stride,
                             num_start,
                             num_end,
                             global_connect_s,
                             buffer_size,
                             buffer,
                             f);
    }

    start_id = end_id;
  }

  /* Free memory */

  elements_slice = fvm_gather_slice_destroy(elements_slice);

  BFT_FREE(sub_elt_vertex_num);
  BFT_FREE(global_idx_s);
  BFT_FREE(local_idx);

  if (global_connect_s != global_connect_s_caller)
    BFT_FREE(global_connect_s);
}

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to an EnSight Gold file in parallel mode
 *
 * parameters:
 *   export_section               <-- pointer to EnSight section helper structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller      --- global connectivity slice provided
 *                                    by caller
 *   text_mode                    <-- true if text output, false if binary
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polyhedra_g(const fvm_writer_section_t  *export_section,
                          const fvm_io_num_t          *global_vertex_num,
                          MPI_Comm            comm,
                          fvm_gnum_t          global_s_size,
                          fvm_gnum_t          global_connect_s_size_caller,
                          fvm_gnum_t          global_connect_s_caller[],
                          _Bool               text_mode,
                          size_t              buffer_size,
                          int                 buffer[],
                          bft_file_t         *f)
{
  fvm_lnum_t  i, j, k, l;
  fvm_gnum_t  i_s;
  size_t      i_buf;

  fvm_lnum_t  cell_length, face_length;
  fvm_lnum_t  face_id;
  fvm_gnum_t  global_num_start;
  fvm_gnum_t  global_num_end;

  int         rank;

  fvm_lnum_t  n_face_lengths_max = 0;

  fvm_lnum_t  *face_lengths = NULL;
  fvm_lnum_t  *cell_vtx_idx = NULL;
  fvm_lnum_t  *cell_vtx_num = NULL;
  fvm_gnum_t  *global_idx_s = NULL;

  fvm_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvm_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvm_gnum_t  *global_connect_s = global_connect_s_caller;

  fvm_gather_slice_t   *polyhedra_slice = NULL;

  const fvm_writer_section_t  *current_section;

  /* Get info on the current MPI communicator */

  MPI_Comm_rank(comm, &rank);

  /* Allocate memory for additionnal indexes */

  BFT_MALLOC(global_idx_s, global_s_size + 1, fvm_gnum_t);

  /* Export number of faces per polyhedron */
  /*---------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_lengths_g(section->global_element_num,
                     section->face_index,
                     comm,
                     global_s_size,
                     global_connect_s_size,
                     global_connect_s,
                     text_mode,
                     buffer_size,
                     buffer,
                     f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Export number of vertices per face per polyhedron */
  /*---------------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    /* Build local polyhedron face lengths information */

    if (section->face_index[section->n_elements] > n_face_lengths_max) {
      BFT_REALLOC(face_lengths,
                  section->face_index[section->n_elements],
                  fvm_lnum_t);
      n_face_lengths_max = section->face_index[section->n_elements];
    }

    l = 0;
    for (i = 0 ; i < section->n_elements ; i++) {
      for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
        face_id = FVM_ABS(section->face_num[j]) - 1;
        face_length = (  section->vertex_index[face_id+1]
                       - section->vertex_index[face_id]);
        face_lengths[l++] = face_length;
      }
    }
    assert(l == section->face_index[section->n_elements]);

    /* Loop on slices */

    polyhedra_slice = fvm_gather_slice_create(section->global_element_num,
                                              global_s_size,
                                              comm);

    while (fvm_gather_slice_advance(polyhedra_slice,
                                    &global_num_start,
                                    &global_num_end) == 0) {

      /* Now build the slice index for number of vertices per face */

      fvm_gather_slice_index(section->face_index,
                             global_idx_s,
                             section->global_element_num,
                             comm,
                             polyhedra_slice);

      /* Recompute slice size */

      fvm_gather_resize_indexed_slice(10,
                                      &global_num_end,
                                      &global_connect_s_size,
                                      comm,
                                      global_idx_s,
                                      polyhedra_slice);

      /* If the buffer passed to this function is too small, allocate a
         larger one; in this case, we may as well keep it for all slices */

      if (global_connect_s_size_prev < global_connect_s_size) {
        if (global_connect_s == global_connect_s_caller)
          global_connect_s = NULL;
        BFT_REALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);
        global_connect_s_size_prev = global_connect_s_size;
      }

      /* Now gather the number of vertices per face */

      fvm_gather_indexed_numbers(section->face_index,
                                 face_lengths,
                                 global_connect_s,
                                 NULL,
                                 section->global_element_num,
                                 comm,
                                 global_idx_s,
                                 polyhedra_slice);

      /* Print number of vertices per face on rank 0 */

      if (rank == 0) {

        size_t  i_start = global_idx_s[0];
        size_t  i_end   = global_idx_s[global_num_end - global_num_start];

        if (text_mode == true) {

          for (i_s = i_start ; i_s < i_end ; i_s++)
            bft_file_printf(f, "%10d\n", global_connect_s[i_s]);

        }
        else { /* text_mode = false */

          for (i_s = i_start, i_buf = 0 ; i_s < i_end ; i_s++) {
            if (i_buf == buffer_size) {
              bft_file_write(buffer, sizeof(int), i_buf, f);
              i_buf = 0;
            }
            buffer[i_buf++] = (int)(global_connect_s[i_s]);
          }
          if (i_buf > 0)
            bft_file_write(buffer, sizeof(int), i_buf, f);

        }

      }

    } /* end of loop on slices */

    fvm_gather_slice_destroy(polyhedra_slice);

    BFT_FREE(face_lengths);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  if (global_connect_s != global_connect_s_caller)
    BFT_FREE(global_connect_s);

  global_connect_s_size = global_connect_s_size_caller;
  global_connect_s = global_connect_s_caller;

  /* Export cell->vertex connectivity by slices */
  /*--------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;
    const fvm_io_num_t  *_global_vertex_num = global_vertex_num;

    BFT_MALLOC(cell_vtx_idx, section->n_elements + 1, fvm_lnum_t);

    l = 0;

    /* In text mode, add zeroes to cell vertex connectivity to mark face
       limits (so as to add newlines) */

    if (text_mode == true) {

      cell_vtx_idx[0] = 0;
      for (i = 0 ; i < section->n_elements ; i++) {
        cell_length = 0;
        for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
          face_id = FVM_ABS(section->face_num[j]) - 1;
          face_length = (  section->vertex_index[face_id+1]
                         - section->vertex_index[face_id]);
          cell_length += face_length + 1;
        }
        cell_vtx_idx[i+1] = cell_vtx_idx[i] + cell_length;
      }

    }
    else { /* In binary mode, simply build true cell -> vertex connectivity */

      cell_vtx_idx[0] = 0;
      for (i = 0 ; i < section->n_elements ; i++) {
        cell_length = 0;
        for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
          face_id = FVM_ABS(section->face_num[j]) - 1;
          face_length = (  section->vertex_index[face_id+1]
                         - section->vertex_index[face_id]);
          cell_length += face_length;
        }
        cell_vtx_idx[i+1] = cell_vtx_idx[i] + cell_length;
      }

    }

    BFT_MALLOC(cell_vtx_num, cell_vtx_idx[section->n_elements], fvm_lnum_t);

    l = 0;

    for (i = 0 ; i < section->n_elements ; i++) {
      for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
        if (section->face_num[j] > 0) {
          face_id = section->face_num[j] - 1;
          for (k = section->vertex_index[face_id] ;
               k < section->vertex_index[face_id+1] ;
               k++)
            cell_vtx_num[l++] = section->vertex_num[k];
        }
        else {
          face_id = -section->face_num[j] - 1;
          k = section->vertex_index[face_id] ;
          cell_vtx_num[l++] = section->vertex_num[k];
          for (k = section->vertex_index[face_id+1] - 1 ;
               k > section->vertex_index[face_id] ;
               k--)
            cell_vtx_num[l++] = section->vertex_num[k];
        }
        if (text_mode == true)
          cell_vtx_num[l++] = 0; /* mark face limits in text mode */
      }

    }

    /* In text mode, we must apply global vertex numberings here so as to add
       zeroes to cell vertex connectivity to mark face limits; in binary mode,
       we have a regular indexed connectivity, so local vertex numbers may be
       converted to global numbers by fvm_gather_...() functions */

    if (text_mode == true) {
      if (global_vertex_num != NULL) {
        const fvm_gnum_t * g_v_num
          = fvm_io_num_get_global_num(global_vertex_num);
        for (l = 0; l < cell_vtx_idx[section->n_elements]; l++) {
          if (cell_vtx_num[l] != 0)
            cell_vtx_num[l] = g_v_num[cell_vtx_num[l] - 1];
        }
      }
      _global_vertex_num = NULL;
    }

    _write_indexed_connect_g(_global_vertex_num,
                             section->global_element_num,
                             cell_vtx_idx,
                             cell_vtx_num,
                             comm,
                             global_s_size,
                             global_connect_s_size,
                             global_connect_s,
                             text_mode,
                             buffer_size,
                             buffer,
                             f);

    BFT_FREE(cell_vtx_idx);
    BFT_FREE(cell_vtx_num);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Free memory */

  BFT_FREE(global_idx_s);

  return current_section;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to an EnSight Gold file in serial mode
 *
 * parameters:
 *   export_section <-- pointer to EnSight section helper structure
 *   text_mode      <-- true if text output, false if binary
 *   buffer_size    <-- size of write buffer
 *   buffer         --- write buffer (for binary mode)
 *   f              <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polyhedra_l(const fvm_writer_section_t  *export_section,
                          _Bool                        text_mode,
                          size_t                       buffer_size,
                          int                          buffer[],
                          bft_file_t                  *f)


{
  fvm_lnum_t  i, j, k, l;
  size_t  i_buf;

  fvm_lnum_t  face_length;
  fvm_lnum_t  face_id;

  int  face_sgn;

  const fvm_writer_section_t  *current_section;

  /* Print cell connectivity directly, without using extra buffers */

  /* Write number of faces per cell */
  /*--------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    if (text_mode == true) {
      for (i = 0 ; i < section->n_elements ; i++)
        bft_file_printf(f, "%10d\n",
                        (int)(  section->face_index[i+1]
                              - section->face_index[i]));
    }
    else {
      for (i = 0, i_buf = 0 ; i < section->n_elements ; i++) {
        if (i_buf == buffer_size) {
          bft_file_write(buffer, sizeof(int), i_buf, f);
          i_buf = 0;
        }
        buffer[i_buf++] = (int)(  section->face_index[i+1]
                                - section->face_index[i]);
      }
      if (i_buf > 0)
        bft_file_write(buffer, sizeof(int), i_buf, f);
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Write number of vertices/face */
  /*-------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    for (i = 0, i_buf = 0 ; i < section->n_elements ; i++) {

      /* Loop on cell faces */

      for (j = section->face_index[i] ;
           j < section->face_index[i+1] ;
           j++) {

        if (section->face_num[j] > 0)
          face_id = section->face_num[j] - 1;
        else
          face_id = -section->face_num[j] - 1;

        face_length = (  section->vertex_index[face_id+1]
                       - section->vertex_index[face_id]);

        if (text_mode == true)
          bft_file_printf(f, "%10d\n",
                          (int)face_length);
        else {
          if (i_buf == buffer_size) {
            bft_file_write(buffer, sizeof(int), i_buf, f);
            i_buf = 0;
          }
          buffer[i_buf++] = (int)face_length;
        }

      }

    }

    if (text_mode == false && i_buf > 0)
      bft_file_write(buffer, sizeof(int), i_buf, f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Write cell/vertex connectivity */
  /*--------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    for (i = 0, i_buf = 0 ; i < section->n_elements ; i++) {

      /* Loop on cell faces */

      for (j = section->face_index[i] ;
           j < section->face_index[i+1] ;
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

        if (text_mode == true) {
          for (k = 0 ; k < face_length ; k++) {
            l =   section->vertex_index[face_id]
                + (face_length + (k*face_sgn))%face_length;
            bft_file_printf(f, "%10d", (int)section->vertex_num[l]);
          }
          bft_file_printf(f, "\n");
        }
        else { /* text_mode == false */
          for (k = 0 ; k < face_length ; k++) {
            l =   section->vertex_index[face_id]
                + (face_length + (k*face_sgn))%face_length;
            if (i_buf == buffer_size) {
              bft_file_write(buffer, sizeof(int), i_buf, f);
              i_buf = 0;
            }
            buffer[i_buf++] = (int)section->vertex_num[l];
          }
        }

      } /* End of loop on cell faces */

    } /* End of loop on polyhedral cells */

    if (text_mode == false && i_buf > 0)
      bft_file_write(buffer, sizeof(int), i_buf, f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to an EnSight Gold file in parallel mode
 *
 * parameters:
 *   export_section               <-- pointer to EnSight section
 *                                    helper structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   n_ranks                      <-- number of processes in communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size        <-- global connectivity slice size
 *   global_connect_s             --- global connectivity slice
 *   text_mode                    <-- true if text output, false if binary
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_g(const fvm_writer_section_t  *export_section,
                         const fvm_io_num_t          *global_vertex_num,
                         MPI_Comm                  comm,
                         fvm_gnum_t                global_s_size,
                         fvm_gnum_t                global_connect_s_size,
                         fvm_gnum_t                global_connect_s[],
                         _Bool                     text_mode,
                         size_t                    buffer_size,
                         int                       buffer[],
                         bft_file_t               *f)
{
  const fvm_writer_section_t  *current_section;

  /* Export number of vertices per polygon by slices */
  /*-------------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_lengths_g(section->global_element_num,
                     section->vertex_index,
                     comm,
                     global_s_size,
                     global_connect_s_size,
                     global_connect_s,
                     text_mode,
                     buffer_size,
                     buffer,
                     f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Export face->vertex connectivity by slices */
  /*--------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_indexed_connect_g(global_vertex_num,
                             section->global_element_num,
                             section->vertex_index,
                             section->vertex_num,
                             comm,
                             global_s_size,
                             global_connect_s_size,
                             global_connect_s,
                             text_mode,
                             buffer_size,
                             buffer,
                             f);

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
 *   text_mode      <-- true if text output, false if binary
 *   buffer_size    <-- size of write buffer
 *   buffer         --- write buffer (for binary mode)
 *   f              <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
*----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_l(const fvm_writer_section_t  *export_section,
                         _Bool                        text_mode,
                         size_t                       buffer_size,
                         int                          buffer[],
                         bft_file_t                  *f)


{
  fvm_lnum_t  i, j;
  size_t  i_buf;

  const fvm_writer_section_t  *current_section = NULL;

  /* Print face connectivity directly, without using extra buffers */

  /* First loop on all polygonal faces, to write number of vertices */
  /*----------------------------------------------------------------*/

  current_section = export_section;

  do { /* Loop on sections that should be grouped */

    const fvm_nodal_section_t  *section = current_section->section;

    if (text_mode == true) {
      for (i = 0 ; i < section->n_elements ; i++)
        bft_file_printf(f, "%10d\n", (int)(  section->vertex_index[i+1]
                                           - section->vertex_index[i]));
    }
    else {
      for (i = 0, i_buf = 0 ; i < section->n_elements ; i++) {
        if (i_buf == buffer_size) {
          bft_file_write(buffer, sizeof(int), i_buf, f);
          i_buf = 0;
        }
        buffer[i_buf++] = (int)(  section->vertex_index[i+1]
                                - section->vertex_index[i]);
      }
      if (i_buf > 0)
        bft_file_write(buffer, sizeof(int), i_buf, f);
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Loop on all polygonal faces */
  /*-----------------------------*/

  current_section = export_section;

  do { /* Loop on sections that should be grouped */

    const fvm_nodal_section_t  *section = current_section->section;

    for (i = 0, i_buf = 0 ; i < section->n_elements ; i++) {

      /* Print face vertex numbers */

      if (text_mode == true) {
        for (j = section->vertex_index[i] ;
             j < section->vertex_index[i+1] ;
             j++)
          bft_file_printf(f, "%10d", (int)section->vertex_num[j]);
        bft_file_printf(f, "\n");
      }
      else { /* text_mode = false */
        for (j = section->vertex_index[i] ;
             j < section->vertex_index[i+1] ;
             j++) {
          if (i_buf == buffer_size) {
            bft_file_write(buffer, sizeof(int), i_buf, f);
            i_buf = 0;
          }
          buffer[i_buf++] = (int)section->vertex_num[j];
        }
      }

    } /* End of loop on polygonal faces */

    if (text_mode == false && i_buf > 0)
      bft_file_write(buffer, sizeof(int), i_buf, f);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write tesselated element connectivity from a nodal mesh to an EnSight Gold
 * file in parallel mode
 *
 * parameters:
 *   export_section               <-- pointer to EnSight section
 *                                    helper structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller      --- global connectivity slice provided
 *                                    by caller
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_tesselated_g(const fvm_writer_section_t  *export_section,
                           const fvm_io_num_t          *global_vertex_num,
                           MPI_Comm                  comm,
                           fvm_gnum_t                global_s_size,
                           fvm_gnum_t                global_connect_s_size,
                           fvm_gnum_t                global_connect_s[],
                           size_t                    buffer_size,
                           int                       buffer[],
                           bft_file_t               *f)
{
  const fvm_writer_section_t  *current_section;

  /* Export face->vertex connectivity by slices */
  /*--------------------------------------------*/

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    _write_tesselated_connect_g(global_vertex_num,
                                section->global_element_num,
                                section->tesselation,
                                current_section->type,
                                current_section->extra_vertex_base,
                                comm,
                                global_s_size,
                                global_connect_s_size,
                                global_connect_s,
                                buffer_size,
                                buffer,
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
 *   export_section               <-- pointer to EnSight section
 *                                    helper structure
 *   tesselation                  <-- element tesselation description
 *   buffer_size                  <-- size of write buffer
 *   buffer                       --- write buffer (for binary mode)
 *   f                            <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_tesselated_l(const fvm_writer_section_t  *export_section,
                           size_t                       buffer_size,
                           int                          buffer[],
                           bft_file_t                  *f)
{
  const fvm_writer_section_t  *current_section;

  current_section = export_section;

  do { /* loop on sections which should be appended */

    const fvm_nodal_section_t  *section = current_section->section;

    fvm_lnum_t  start_id, end_id;
    fvm_lnum_t  n_sub_elements_max;
    fvm_lnum_t  n_buffer_elements_max = section->n_elements;
    fvm_lnum_t *vertex_num = NULL;

    const fvm_lnum_t *sub_element_idx
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
               fvm_lnum_t);

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
                       buffer_size,
                       buffer,
                       f);

    }

    BFT_FREE(vertex_num);

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
 *   field_values       <-- pointer to arre (output is real)
 *   rank               <-- rank of current process in communicator
 *   output_buffer_size <-- output buffer size
 *   output_buffer      --- output buffer (size output_buffer_size)
 *   f                  <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_field_values_n(const fvm_nodal_t           *mesh,
                       fvm_writer_field_helper_t   *helper,
                       int                          input_dim,
                       fvm_interlace_t              interlace,
                       int                          n_parent_lists,
                       const fvm_lnum_t             parent_num_shift[],
                       fvm_datatype_t               datatype,
                       const void            *const field_values[],
                       int                          rank,
                       size_t                       output_buffer_size,
                       float                        output_buffer[],
                       bft_file_t                  *f)
{
  int  i;
  size_t  output_size;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  for (i = 0 ; i < output_dim ; i++) {

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

      if (rank == 0)
        _write_slice_values(0,
                            output_size,
                            output_buffer,
                            f);

    }

  }

}

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to
 * an EnSight Gold file.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   export_section     <-- pointer to EnSight section helper structure
 *   helper             <-- pointer to general writer helper structure
 *   input_dim          <-- input field dimension
 *   interlace          <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   rank               <-- rank of current process in communicator
 *   output_buffer_size <-- output buffer size
 *   output_buffer      --- output buffer (size output_buffer_size)
 *   f                  <-- pointer to associated file
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_field_values_e(const fvm_writer_section_t      *export_section,
                       fvm_writer_field_helper_t       *helper,
                       int                              input_dim,
                       fvm_interlace_t                  interlace,
                       int                              n_parent_lists,
                       const fvm_lnum_t                 parent_num_shift[],
                       fvm_datatype_t                   datatype,
                       const void                *const field_values[],
                       int                              rank,
                       size_t                           output_buffer_size,
                       float                            output_buffer[],
                       bft_file_t                      *f)
{
  int  i;
  size_t  output_size;

  const fvm_writer_section_t  *current_section = NULL;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  /* Loop on dimension (de-interlace vectors, always 3D for EnSight) */

  for (i = 0 ; i < output_dim ; i++) {

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

        if (rank == 0)
          _write_slice_values(0,
                              output_size,
                              output_buffer,
                              f);

      }

      current_section = current_section->next;

      if (   current_section == NULL
          || current_section->continues_previous == false)
        loop_on_sections = false;

    } /* while (loop on sections) */

  } /* end of loop on spatial dimension */

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
fvm_to_ensight_v1_init_writer(const char             *name,
                              const char             *path,
                              const char             *options,
                              fvm_writer_time_dep_t   time_dependency,
                              MPI_Comm                comm)
#else
void *
fvm_to_ensight_v1_init_writer(const char             *name,
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
    int mpi_flag, rank, n_ranks;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag && comm != MPI_COMM_NULL) {
      this_writer->comm = comm;
      MPI_Comm_rank(this_writer->comm, &rank);
      MPI_Comm_size(this_writer->comm, &n_ranks);
      this_writer->rank = rank;
      this_writer->n_ranks = n_ranks;
    }
    else
      this_writer->comm = MPI_COMM_NULL;
  }
#endif /* defined(HAVE_MPI) */

  /* Parse options */

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1 ; i2 < l_tot && options[i2] != ' ' ; i2++);
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

      for (i1 = i2 + 1 ; i1 < l_tot && options[i1] == ' ' ; i1++);

    }

  }

  this_writer->case_info = fvm_to_ensight_case_create(name,
                                                      path,
                                                      time_dependency);

  /* Initialize geometry file name */

  if (time_dependency == FVM_WRITER_FIXED_MESH)
    _init_geom_file(this_writer);

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
fvm_to_ensight_v1_finalize_writer(void  *this_writer_p)
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
fvm_to_ensight_v1_set_mesh_time(void          *this_writer_p,
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
fvm_to_ensight_v1_needs_tesselation(fvm_writer_t       *this_writer_p,
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

    for (i = 0 ; i < mesh->n_sections ; i++) {

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
fvm_to_ensight_v1_export_nodal(void               *this_writer_p,
                               const fvm_nodal_t  *mesh)
{
  int     part_num;
  size_t  write_connect_buffer_size = 0;
  int    *write_connect_buffer = NULL;

  const fvm_writer_section_t  *export_section = NULL;
  fvm_writer_section_t  *export_list = NULL;
  fvm_to_ensight_writer_t  *this_writer
                             = (fvm_to_ensight_writer_t *)this_writer_p;
  bft_file_t  *f = NULL;

  fvm_gnum_t   global_connect_s_size, global_s_size;
  fvm_gnum_t   n_g_vertices;
  fvm_gnum_t  *n_g_elements_section = NULL;

#if defined(HAVE_MPI)

  fvm_gnum_t  *global_connect_s = NULL;
  MPI_Comm    comm = this_writer->comm;

#endif

  const int  rank = this_writer->rank;
  const int  n_ranks = this_writer->n_ranks;

  /* Initialization */
  /*----------------*/

  /* Get part number */

  part_num = fvm_to_ensight_case_get_part_num(this_writer->case_info,
                                              mesh->name);
  if (part_num == 0)
    part_num = fvm_to_ensight_case_add_part(this_writer->case_info,
                                            mesh->name);

  /* Open geometry file in append mode */

  if (rank == 0) {

    bft_file_type_t                  file_type;
    fvm_to_ensight_case_file_info_t  file_info;

    if (this_writer->text_mode == true)
      file_type = BFT_FILE_TYPE_TEXT;
    else

      file_type = BFT_FILE_TYPE_BINARY;

    file_info = fvm_to_ensight_case_get_geom_file(this_writer->case_info);

    if (file_info.queried == false)
      _init_geom_file(this_writer); /* Create, write header, and close file */

    f = bft_file_open(file_info.name,
                      BFT_FILE_MODE_APPEND,
                      file_type);

    if (this_writer->swap_endian == true)
      bft_file_set_swap_endian(f, 1);
  }

  /* Part header */

  if (rank == 0) {
    _write_string(f, "part");
    _write_int(f, part_num);
    if (mesh->name != NULL)
      _write_string(f, mesh->name);
    else
      _write_string(f, "unnamed");
  }

  /* Build list of sections that are used here, in order of output */

  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       this_writer->discard_polygons,
                                       this_writer->discard_polyhedra,
                                       this_writer->divide_polygons,
                                       this_writer->divide_polyhedra);

  /* Buffer and global sizes required in parallel mode and for binary output */
  /*-------------------------------------------------------------------------*/

  BFT_MALLOC(n_g_elements_section, mesh->n_sections, fvm_gnum_t);

  fvm_writer_def_nodal_buf_size(mesh,
                                n_ranks,
                                12,
                                5,
                                &n_g_vertices,
                                n_g_elements_section,
                                &global_s_size,
                                &global_connect_s_size);

  /* Avoid too many small communications with large processor counts */

  if (n_ranks > 1 && global_connect_s_size > 0) {
    size_t min_buffer_size =   fvm_parall_get_min_coll_buf_size()
                             / sizeof(fvm_gnum_t);
    if (min_buffer_size > global_connect_s_size) {
      fvm_gnum_t global_s_size_min = global_s_size * (  min_buffer_size
                                                      / global_connect_s_size);
      if (global_s_size_min > global_s_size) {
        global_s_size = global_s_size_min;
        global_connect_s_size = min_buffer_size;
      }
    }
  }

  /* Vertex coordinates */
  /*--------------------*/

#if defined(HAVE_MPI)

  if (n_ranks > 1 )
    _export_vertex_coords_g(this_writer,
                            mesh,
                            comm,
                            global_s_size,
                            f);

#endif

  if (n_ranks == 1)
    _export_vertex_coords_l(this_writer,
                            mesh,
                            f);

  /* Allocate connectivity buffer for use wih all types of elements */

#if defined(HAVE_MPI)

  if (n_ranks > 1)
    BFT_MALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);

#endif

  /* Buffer required in binary mode */
  /*--------------------------------*/

  if (rank == 0) {
    if (this_writer->text_mode == true) {
      write_connect_buffer_size = 0;
      write_connect_buffer = NULL;
    }
    else { /* Arbitrary write buffer size, small enough to add little
              additional memory requirement (in proportion), large enough
              to limit number of write calls */
      write_connect_buffer_size = (global_s_size + 4) / 4;
      BFT_MALLOC(write_connect_buffer, write_connect_buffer_size, int);
    }
  }

  /* If no sections are present (i.e. we may only have vertices),
     add  "point" elements */

  if (export_list == NULL && rank == 0)
    _export_point_elements_l(n_g_vertices,
                             write_connect_buffer_size,
                             write_connect_buffer,
                             f);

  /* Element connectivity */
  /*----------------------*/

  export_section = export_list;

  while (export_section != NULL) {

    const fvm_nodal_section_t  *section = export_section->section;

    /* Print header if start of corresponding EnSight section */

    if (export_section->continues_previous == false && rank == 0) {

      fvm_gnum_t n_g_elements = 0;
      const fvm_writer_section_t  *next_section = export_section;

      do {

        if (next_section->section->type == export_section->type)
          n_g_elements += fvm_nodal_section_n_g_elements(next_section->section);

        else {
          fvm_gnum_t n_g_sub_elements = 0;
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

      if (n_ranks > 1) { /* start of output in parallel mode */

        fvm_gnum_t   global_num_start;
        fvm_gnum_t   global_num_end;

        fvm_gather_slice_t  *elements_slice = NULL;

        elements_slice
          = fvm_gather_slice_create(section->global_element_num,
                                    global_s_size,
                                    comm);

        while (fvm_gather_slice_advance(elements_slice,
                                        &global_num_start,
                                        &global_num_end) == 0) {

          fvm_gather_strided_connect(section->vertex_num,
                                     global_connect_s,
                                     section->stride,
                                     mesh->global_vertex_num,
                                     section->global_element_num,
                                     comm,
                                     elements_slice);

          if (rank == 0)
            _write_slice_connect_g(section->stride,
                                   global_num_start,
                                   global_num_end,
                                   global_connect_s,
                                   write_connect_buffer_size,
                                   write_connect_buffer,
                                   f);

        }

        fvm_gather_slice_destroy(elements_slice);

      } /* end of output in parallel mode */

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1) { /* start of output in serial mode */

        _write_connect_l(section->stride,
                         section->n_elements,
                         section->vertex_num,
                         write_connect_buffer_size,
                         write_connect_buffer,
                         f);

      }

      export_section = export_section->next;

    } /* end of output for strided element types */

    /* Output for tesselated polygons or polyhedra */
    /*---------------------------------------------*/

    else if (export_section->type != section->type) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1) {

        export_section = _export_nodal_tesselated_g(export_section,
                                                    mesh->global_vertex_num,
                                                    comm,
                                                    global_s_size,
                                                    global_connect_s_size,
                                                    global_connect_s,
                                                    write_connect_buffer_size,
                                                    write_connect_buffer,
                                                    f);

      }

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1) {

        export_section = _export_nodal_tesselated_l(export_section,
                                                    write_connect_buffer_size,
                                                    write_connect_buffer,
                                                    f);

      }

    }

    /* Output for polygons */
    /*---------------------*/

    else if (export_section->type == FVM_FACE_POLY) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1) {

        export_section = _export_nodal_polygons_g(export_section,
                                                  mesh->global_vertex_num,
                                                  comm,
                                                  global_s_size,
                                                  global_connect_s_size,
                                                  global_connect_s,
                                                  this_writer->text_mode,
                                                  write_connect_buffer_size,
                                                  write_connect_buffer,
                                                  f);

      }

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1) {

        /* Export as polygons */

        export_section = _export_nodal_polygons_l(export_section,
                                                  this_writer->text_mode,
                                                  write_connect_buffer_size,
                                                  write_connect_buffer,
                                                  f);

      }

    }

    /* Output for polyhedra */
    /*----------------------*/

    else if (export_section->type == FVM_CELL_POLY) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1)

        export_section =_export_nodal_polyhedra_g(export_section,
                                                  mesh->global_vertex_num,
                                                  comm,
                                                  global_s_size,
                                                  global_connect_s_size,
                                                  global_connect_s,
                                                  this_writer->text_mode,
                                                  write_connect_buffer_size,
                                                  write_connect_buffer,
                                                  f);

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)

        export_section = _export_nodal_polyhedra_l(export_section,
                                                   this_writer->text_mode,
                                                   write_connect_buffer_size,
                                                   write_connect_buffer,
                                                   f);

    }

  } /* End of loop on sections */

  BFT_FREE(export_list);

  /* Free buffers */
  /*--------------*/

  BFT_FREE(n_g_elements_section);

  if (write_connect_buffer != NULL) {
    write_connect_buffer_size = 0;
    BFT_FREE(write_connect_buffer);
  }

  /* Free buffers required in parallel mode */

#if defined(HAVE_MPI)

  if (n_ranks > 1)
    BFT_FREE(global_connect_s);

#endif /* defined(HAVE_MPI) */

  /* Close geometry file and update case file */
  /*------------------------------------------*/

  if (rank == 0)
    bft_file_free(f);

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
fvm_to_ensight_v1_export_field(void                  *this_writer_p,
                               const fvm_nodal_t     *mesh,
                               const char            *name,
                               fvm_writer_var_loc_t   location,
                               int                    dimension,
                               fvm_interlace_t        interlace,
                               int                    n_parent_lists,
                               const fvm_lnum_t       parent_num_shift[],
                               fvm_datatype_t         datatype,
                               int                    time_step,
                               double                 time_value,
                               const void      *const field_values[])
{
  int     output_dim;
  int     part_num;
  size_t  input_size = 0;
  size_t  output_size = 0;
  size_t  min_var_buffer_size = 0;
  size_t  var_buffer_size = 0;
  float  *var_buffer = NULL;

  const fvm_writer_section_t  *export_section = NULL;
  fvm_writer_field_helper_t  *helper = NULL;
  fvm_writer_section_t  *export_list = NULL;
  fvm_to_ensight_writer_t  *this_writer
                             = (fvm_to_ensight_writer_t *)this_writer_p;
  bft_file_t  *f = NULL;

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

  if (rank == 0) {

    bft_file_type_t                  file_type;
    bft_file_mode_t                  file_mode;
    fvm_to_ensight_case_file_info_t  file_info;

    if (this_writer->text_mode == true)
      file_type = BFT_FILE_TYPE_TEXT;
    else
      file_type = BFT_FILE_TYPE_BINARY;

    file_info = fvm_to_ensight_case_get_var_file(this_writer->case_info,
                                                 name,
                                                 output_dim,
                                                 location,
                                                 time_step,
                                                 time_value);

    if (file_info.queried == true)
      file_mode = BFT_FILE_MODE_APPEND;
    else
      file_mode = BFT_FILE_MODE_WRITE;

    f = bft_file_open(file_info.name,
                      file_mode,
                      file_type);

    if (this_writer->swap_endian == true)
      bft_file_set_swap_endian(f, 1);

    /* New files start with description line */

    if (file_mode == BFT_FILE_MODE_WRITE){
      char buf[81] = "";
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

  helper = fvm_writer_field_helper_create(mesh,
                                          export_list,
                                          output_dim,
                                          FVM_NO_INTERLACE,
                                          FVM_FLOAT,
                                          location);

#if defined(HAVE_MPI)

  fvm_writer_field_helper_init_g(helper,
                                 export_list,
                                 mesh,
                                 this_writer->comm);

#endif

  /* Buffer size computation and allocation */
  /*----------------------------------------*/

  fvm_writer_field_helper_get_size(helper,
                                   &input_size,
                                   &output_size,
                                   NULL,
                                   &min_var_buffer_size);

  /* Slicing allows for arbitrary buffer size, but should be small enough
     to add little additional memory requirement (in proportion), large
     enough to limit number of write and gather calls. */

  if (n_ranks > 1) {
    size_t min_buffer_size = fvm_parall_get_min_coll_buf_size() / sizeof(float);
    var_buffer_size = input_size / n_ranks;
    if (var_buffer_size < min_buffer_size)
      var_buffer_size = min_buffer_size;
  }
  else
    var_buffer_size = input_size / 4;

  var_buffer_size = FVM_MAX(var_buffer_size, min_var_buffer_size);
  var_buffer_size = FVM_MAX(var_buffer_size, 128);
  var_buffer_size = FVM_MIN(var_buffer_size, output_size);

  BFT_MALLOC(var_buffer, var_buffer_size, float);

  /* Part header */

  if (rank == 0) {
    _write_string(f, "part");
    _write_int(f, part_num);
  }

  /* Per node variable */
  /*-------------------*/

  if (location == FVM_WRITER_PER_NODE) {

    if (rank == 0)
      _write_string(f, "coordinates");

    _export_field_values_n(mesh,
                           helper,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           rank,
                           var_buffer_size,
                           var_buffer,
                           f);

  }

  /* Per element variable */
  /*----------------------*/

  else if (location == FVM_WRITER_PER_ELEMENT) {

    export_section = export_list;

    while (export_section != NULL) {

      /* Print header if start of corresponding EnSight section */

      if (export_section->continues_previous == false && rank == 0)
        _write_string(f, _ensight_type_name[export_section->type]);

      /* Output per grouped sections */

      export_section = _export_field_values_e(export_section,
                                              helper,
                                              dimension,
                                              interlace,
                                              n_parent_lists,
                                              parent_num_shift,
                                              datatype,
                                              field_values,
                                              rank,
                                              var_buffer_size,
                                              var_buffer,
                                              f);

    } /* End of loop on sections */

  } /* End for per element variable */

  /* Free buffers and helper structures */
  /*------------------------------------*/

  BFT_FREE(var_buffer);

  helper = fvm_writer_field_helper_destroy(helper);

  BFT_FREE(export_list);

  /* Close variable file and update case file */
  /*------------------------------------------*/

  if (rank == 0)
    bft_file_free(f);

  fvm_to_ensight_case_write_case(this_writer->case_info, rank);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
