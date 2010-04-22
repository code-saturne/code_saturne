/*============================================================================
 * Write a nodal representation associated with a mesh to file
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2006  EDF

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

#include <bft_mem.h>
#include <bft_file.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_config_defs.h"
#include "fvm_defs.h"
#include "fvm_gather.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_parall.h"
#include "fvm_writer_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Text writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  bft_file_t  *file;      /* Output file */

  fvm_writer_time_dep_t   time_dependency; /* Mesh time dependency */

  int          rank;      /* Rank of current process in communicator */
  int          n_ranks;   /* Number of processes in communicator */

#if defined(HAVE_MPI)
  MPI_Comm     comm;      /* Associated MPI communicator */
#endif

} fvm_to_text_writer_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write slice of a vector of doubles to a text file
 *
 * parameters:
 *   stride         <-- number of values per element
 *   num_start      <-- global number of first element for this slice
 *   num_end        <-- global number of past the last element for this slice
 *   values         <-- pointer to values slice array
 *   f              <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_slice_vector(const size_t             stride,
                    const fvm_gnum_t         num_start,
                    const fvm_gnum_t         num_end,
                    const double             values[],
                    bft_file_t        *const f)
{
  size_t  i, k;
  fvm_gnum_t  j;

  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  switch(stride) {

  case 1:
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12.5f\n",
                      (unsigned long long)j,
                      values[i]);
    break;

  case 2:
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12.5f %12.5f\n",
                      (unsigned long long)j,
                      values[2*i], values[2*i+1]);
    break;

  case 3:
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12.5f %12.5f %12.5f\n",
                      (unsigned long long)j,
                      values[3*i], values[3*i+1], values[3*i+2]);
    break;

  default: /* Fallback, requiring more calls */
    for (i = 0, j = num_start ; j < num_end ; i++, j++) {
      bft_file_printf(f, "%12llu :", (unsigned long long)j);
      for (k = 0 ; k < stride ; k++)
        bft_file_printf(f, " %12.5f",
                        values[stride*i+k]);
      bft_file_printf(f, "\n");
    }
    break;
  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided global connectivity slice to a text file
 *
 * parameters:
 *   stride           <-- number of vertices per element type
 *   num_start        <-- global number of first element for this slice
 *   num_end          <-- global number of last element for this slice
 *   global_connect_s <-- global connectivity slice array
 *   f                <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_slice_connect_g(const int                stride,
                       const fvm_gnum_t         num_start,
                       const fvm_gnum_t         num_end,
                       const fvm_gnum_t         global_connect_s[],
                       bft_file_t        *const f)
{
  fvm_gnum_t  i, j;

  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  switch(stride) {

  case 2: /* edge */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12llu %12llu\n",
                      (unsigned long long)j,
                      (unsigned long long)global_connect_s[i*2],
                      (unsigned long long)global_connect_s[i*2+1]);
    break;

  case 3: /* FVM_FACE_TRIA */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu\n",
                      (unsigned long long)j,
                      (unsigned long long)global_connect_s[i*3],
                      (unsigned long long)global_connect_s[i*3+1],
                      (unsigned long long)global_connect_s[i*3+2]);
    break;

  case 4: /* FVM_FACE_QUAD or FVM_CELL_TETRA */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu %12llu\n",
                      (unsigned long long)j,
                      (unsigned long long)global_connect_s[i*4],
                      (unsigned long long)global_connect_s[i*4+1],
                      (unsigned long long)global_connect_s[i*4+2],
                      (unsigned long long)global_connect_s[i*4+3]);
    break;

  case 5: /* FVM_CELL_PYRAM */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu %12llu %12llu\n",
                      (unsigned long long)j,
                      (unsigned long long)global_connect_s[i*5],
                      (unsigned long long)global_connect_s[i*5+1],
                      (unsigned long long)global_connect_s[i*5+2],
                      (unsigned long long)global_connect_s[i*5+3],
                      (unsigned long long)global_connect_s[i*5+4]);
    break;

  case 6: /* FVM_CELL_PRISM */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu %12llu %12llu %12llu\n",
                      (unsigned long long)j,
                      (unsigned long long)global_connect_s[i*6],
                      (unsigned long long)global_connect_s[i*6+1],
                      (unsigned long long)global_connect_s[i*6+2],
                      (unsigned long long)global_connect_s[i*6+3],
                      (unsigned long long)global_connect_s[i*6+4],
                      (unsigned long long)global_connect_s[i*6+5]);
    break;

  case 8: /* FVM_CELL_HEXA */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bft_file_printf(f,
                      "%12llu : "
                      "%12llu %12llu %12llu %12llu\n"
                      "               "
                      "%12llu %12llu %12llu %12llu\n",
                      (unsigned long long)j,
                      (unsigned long long)global_connect_s[i*8],
                      (unsigned long long)global_connect_s[i*8+1],
                      (unsigned long long)global_connect_s[i*8+2],
                      (unsigned long long)global_connect_s[i*8+3],
                      (unsigned long long)global_connect_s[i*8+4],
                      (unsigned long long)global_connect_s[i*8+5],
                      (unsigned long long)global_connect_s[i*8+6],
                      (unsigned long long)global_connect_s[i*8+7]);
    break;

  default:
    assert(0);
  }

}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write strided local connectivity to a text file
 *
 * parameters:
 *   stride  <-- number of vertices per element type
 *   n_elems <-- number of elements
 *   connect <-- connectivity array
 *   f       <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_connect_l(const int                stride,
                 const fvm_lnum_t         n_elems,
                 const fvm_lnum_t         connect[],
                 bft_file_t        *const f)
{
  fvm_lnum_t  i;

  switch(stride) {

  case 2: /* edge */
    for (i = 0 ; i < n_elems ; i++)
      bft_file_printf(f, "%12llu : %12llu %12llu\n",
                      (unsigned long long)(i+1),
                      (unsigned long long)connect[i*2],
                      (unsigned long long)connect[i*2+1]);
    break;

  case 3: /* FVM_FACE_TRIA */
    for (i = 0 ; i < n_elems ; i++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu\n",
                      (unsigned long long)(i+1),
                      (unsigned long long)connect[i*3],
                      (unsigned long long)connect[i*3+1],
                      (unsigned long long)connect[i*3+2]);
    break;

  case 4: /* FVM_FACE_QUAD or FVM_CELL_TETRA */
    for (i = 0 ; i < n_elems ; i++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu %12llu\n",
                      (unsigned long long)(i+1),
                      (unsigned long long)connect[i*4],
                      (unsigned long long)connect[i*4+1],
                      (unsigned long long)connect[i*4+2],
                      (unsigned long long)connect[i*4+3]);
    break;

  case 5: /* FVM_CELL_PYRAM */
    for (i = 0 ; i < n_elems ; i++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu %12llu %12llu\n",
                      (unsigned long long)(i+1),
                      (unsigned long long)connect[i*5],
                      (unsigned long long)connect[i*5+1],
                      (unsigned long long)connect[i*5+2],
                      (unsigned long long)connect[i*5+3],
                      (unsigned long long)connect[i*5+4]);
    break;

  case 6: /* FVM_CELL_PRISM */
    for (i = 0 ; i < n_elems ; i++)
      bft_file_printf(f, "%12llu : %12llu %12llu %12llu %12llu %12llu %12llu\n",
                      (unsigned long long)(i+1),
                      (unsigned long long)connect[i*6],
                      (unsigned long long)connect[i*6+1],
                      (unsigned long long)connect[i*6+2],
                      (unsigned long long)connect[i*6+3],
                      (unsigned long long)connect[i*6+4],
                      (unsigned long long)connect[i*6+5]);
    break;

  case 8: /* FVM_CELL_HEXA */
    for (i = 0 ; i < n_elems ; i++)
      bft_file_printf(f,
                      "%12llu : "
                      "%12llu %12llu %12llu %12llu\n"
                      "               "
                      "%12llu %12llu %12llu %12llu\n",
                      (unsigned long long)(i+1),
                      (unsigned long long)connect[i*8],
                      (unsigned long long)connect[i*8+1],
                      (unsigned long long)connect[i*8+2],
                      (unsigned long long)connect[i*8+3],
                      (unsigned long long)connect[i*8+4],
                      (unsigned long long)connect[i*8+5],
                      (unsigned long long)connect[i*8+6],
                      (unsigned long long)connect[i*8+7]);
    break;

  default:
    assert(0);
  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to a text file in parallel mode
 *
 * parameters:
 *   section                      <-- pointer to nodal mesh section structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   rank                         <-- rank in communicator
 *   n_ranks                      <-- number of processes in communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller       -- global connectivity slice provided
 *                                    by caller
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polyhedra_g(const fvm_nodal_section_t  *const section,
                          const fvm_io_num_t         *const global_vertex_num,
                          MPI_Comm                  comm,
                          const int                 rank,
                          const fvm_gnum_t          global_s_size,
                          const fvm_gnum_t          global_connect_s_size_caller,
                          fvm_gnum_t                global_connect_s_caller[],
                          bft_file_t         *const f)

{
  fvm_lnum_t  i, j, k, l;
  fvm_gnum_t  i_s, j_s, k_s;

  fvm_lnum_t  cell_length, face_length;
  fvm_lnum_t  face_id;
  fvm_gnum_t  global_num_start = 1;
  fvm_gnum_t  global_num_end = 0;
  fvm_gnum_t  global_cell_face_idx_shift = 0 ;
  fvm_gnum_t  global_cell_vtx_idx_shift = 0 ;

  fvm_lnum_t  *face_lengths = NULL;
  fvm_lnum_t  *cell_vtx_idx = NULL;
  fvm_lnum_t  *cell_connect = NULL;
  fvm_gnum_t  *global_cell_face_idx_s = NULL;
  fvm_gnum_t  *global_face_lengths_s = NULL;
  fvm_gnum_t  *global_cell_vtx_idx_s = NULL;

  fvm_gnum_t   global_face_lengths_s_size = 0;
  fvm_gnum_t   global_face_lengths_s_size_prev = 0;

  fvm_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvm_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvm_gnum_t  *global_connect_s = global_connect_s_caller;

  fvm_gather_slice_t   *polyhedra_slice = NULL;

  const fvm_gnum_t   n_g_polyhedra = fvm_io_num_get_global_count
                                       (section->global_element_num);

  /* Allocate memory for additionnal indexes */

  BFT_MALLOC(global_cell_face_idx_s, global_s_size + 1, fvm_gnum_t);
  BFT_MALLOC(global_cell_vtx_idx_s, global_s_size + 1, fvm_gnum_t);

  /* Every face should have at least 3 vertices, so cell->vertex connectivity
     should be at least 3 times the size of the cell->face connectivity;
     So we choose 1/3 (rounded up) of the size of the cell->vertex slice
     buffer as the initial size of the face_lengths slice buffer */

  global_face_lengths_s_size = (global_connect_s_size / 3) + 1;
  global_face_lengths_s_size_prev = global_face_lengths_s_size;
  BFT_MALLOC(global_face_lengths_s, global_face_lengths_s_size, fvm_gnum_t);

  /* Build local polyhedron indexes and connectivity information */

  BFT_MALLOC(cell_vtx_idx,
             section->n_elements + 1,
             fvm_lnum_t);

  BFT_MALLOC(face_lengths,
             section->face_index[section->n_elements],
             fvm_lnum_t);

  j_s = 0;
  l = 0;

  cell_vtx_idx[0] = 0;
  for (i = 0 ; i < section->n_elements ; i++) {
    cell_length = 0;
    for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
      face_id = FVM_ABS(section->face_num[j]) - 1;
      face_length = (  section->vertex_index[face_id+1]
                     - section->vertex_index[face_id]);
      face_lengths[l++] = face_length;
      cell_length += face_length;
    }
    cell_vtx_idx[i+1] = cell_vtx_idx[i] + cell_length;
  }

  BFT_MALLOC(cell_connect,
             cell_vtx_idx[section->n_elements],
             fvm_lnum_t);

  l = 0;

  for (i = 0 ; i < section->n_elements ; i++) {
    for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
      if (section->face_num[j] > 0) {
        face_id = section->face_num[j] - 1;
        for (k = section->vertex_index[face_id] ;
             k < section->vertex_index[face_id+1] ;
             k++)
          cell_connect[l++] = section->vertex_num[k];
      }
      else {
        face_id = -section->face_num[j] - 1;
        k = section->vertex_index[face_id] ;
        cell_connect[l++] = section->vertex_num[k];
          for (k = section->vertex_index[face_id+1] - 1 ;
               k > section->vertex_index[face_id] ;
               k--)
            cell_connect[l++] = section->vertex_num[k];
      }
    }
  }

  /* Export by slices */

  polyhedra_slice = fvm_gather_slice_create(section->global_element_num,
                                            global_s_size,
                                            comm);

  while (fvm_gather_slice_advance(polyhedra_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather cell->vertices index */

    fvm_gather_slice_index(cell_vtx_idx,
                           global_cell_vtx_idx_s,
                           section->global_element_num,
                           comm,
                           polyhedra_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvm_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    comm,
                                    global_cell_vtx_idx_s,
                                    polyhedra_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_caller)
        global_connect_s = NULL;
      BFT_REALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now gather cell->vertices connectivity */

    fvm_gather_indexed_numbers(cell_vtx_idx,
                               cell_connect,
                               global_connect_s,
                               global_vertex_num,
                               section->global_element_num,
                               comm,
                               global_cell_vtx_idx_s,
                               polyhedra_slice);

    /* Now build the slice index for number of vertices per face */

    fvm_gather_slice_index(section->face_index,
                           global_cell_face_idx_s,
                           section->global_element_num,
                           comm,
                           polyhedra_slice);

    /* If the face_lengths slice buffer is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (rank == 0) {
      global_face_lengths_s_size =
        FVM_MAX(global_face_lengths_s_size,
                global_cell_face_idx_s[global_num_end - global_num_start]);
    }
    MPI_Bcast(&global_face_lengths_s_size, 1, FVM_MPI_GNUM, 0, comm);
    if (global_face_lengths_s_size_prev < global_face_lengths_s_size) {
      BFT_REALLOC(global_face_lengths_s,
                  global_face_lengths_s_size, fvm_gnum_t);
      global_face_lengths_s_size_prev = global_face_lengths_s_size;
    }

    /* Now gather the number of vertices per face */

    fvm_gather_indexed_numbers(section->face_index,
                               face_lengths,
                               global_face_lengths_s,
                               NULL,
                               section->global_element_num,
                               comm,
                               global_cell_face_idx_s,
                               polyhedra_slice);

    /* Do all printing for cells on rank 0 */

    if (rank == 0) {

      int  line_values;
      char str_num_cell[16];
      char str_num_face[16];
      char str_idx_cell[32];
      char str_idx_face[32];

      /* Print cell connectivity */

      k_s = 0;

      /* Loop on polyhedral cells in slice */

      for (i = 0, i_s = global_num_start ; i_s < global_num_end ; i++, i_s++) {

        /* Loop on cell faces */

        for (j = 0, j_s = global_cell_face_idx_s[i] ;
             j_s < global_cell_face_idx_s[i+1] ;
             j++, j_s++) {

          /* Print cell and face numbers and indexes */

          if (j_s == global_cell_face_idx_s[i]) {
            sprintf(str_num_cell, "%12llu", (unsigned long long)i_s);
            sprintf(str_idx_cell, "[%llu] :",
                    (unsigned long long)(j_s + global_cell_face_idx_shift + 1));
          }
          else {
            str_num_cell[0] = '\0';
            str_idx_cell[0] = '\0';
          }
          sprintf(str_num_face, "%5u", (unsigned)(j+1));
          sprintf(str_idx_face, "[%llu] :",
                    (unsigned long long)(k_s + global_cell_vtx_idx_shift + 1));

          bft_file_printf(f, "%12s %14s %5s %14s",
                          str_num_cell, str_idx_cell,
                          str_num_face, str_idx_face);

          /* Print face vertex numbers */

          line_values = 0;
          for (k = 0 ; k < (fvm_lnum_t)global_face_lengths_s[j_s] ; k++) {
            if (line_values > 2) {
              line_values = 0;
              bft_file_printf(f,"\n%48s", "");
            }
            bft_file_printf(f, " %12llu",
                            (unsigned long long)global_connect_s[k_s++]);
            line_values++;
          }
          bft_file_printf(f, "\n");

        } /* End of loop on cell faces */

        assert(k_s == global_cell_vtx_idx_s[i+1]);

      } /* End of loop on polyhedral cells in slice */

      if (global_num_end > n_g_polyhedra) {
        str_num_cell[0] = '\0';
        sprintf(str_idx_cell, "[%llu] :",
                (unsigned long long)(j_s + global_cell_face_idx_shift + 1));
        str_num_face[0] = '\0';
        sprintf(str_idx_face, "[%llu] :",
                (unsigned long long)(k_s + global_cell_vtx_idx_shift + 1));
        bft_file_printf(f, "%12s %14s %5s %14s",
                        str_num_cell, str_idx_cell,
                        str_num_face, str_idx_face);
      }

      /* Update shift for conversion from slice index to true index */

      global_cell_vtx_idx_shift
        += global_cell_vtx_idx_s[global_num_end - global_num_start];
      global_cell_face_idx_shift
        += global_cell_face_idx_s[global_num_end - global_num_start];

    } /* End of printing for rank 0 for this slice */

  }

  fvm_gather_slice_destroy(polyhedra_slice);

  /* Free memory */

  BFT_FREE(global_cell_face_idx_s);
  BFT_FREE(global_cell_vtx_idx_s);
  BFT_FREE(global_face_lengths_s);
  BFT_FREE(cell_vtx_idx);
  BFT_FREE(face_lengths);
  BFT_FREE(cell_connect);

  if (global_connect_s != global_connect_s_caller)
    BFT_FREE(global_connect_s);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to a text file in serial mode
 *
 * parameters:
 *   section      <-- pointer to nodal mesh section structure
 *   f            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polyhedra_l(const fvm_nodal_section_t  *const section,
                          bft_file_t                 *const f)

{
  fvm_lnum_t  i, j, k, l;

  fvm_lnum_t  connect_length, face_length;
  fvm_lnum_t  face_id;

  int  face_sgn, line_values;
  char str_num_cell[16];
  char str_num_face[16];
  char str_idx_cell[32];
  char str_idx_face[32];

  /* Print cell connectivity directly, without using extra buffers */

  connect_length = 0;
  j = 0;

  /* Loop on all polyhedral cells */

  for (i = 0 ; i < section->n_elements ; i++) {

    /* Loop on cell faces */

    for (j = section->face_index[i] ;
         j < section->face_index[i+1] ;
         j++) {

      /* Print cell and face numbers and indexes */

      if (j == section->face_index[i]) {
        sprintf(str_num_cell, "%12llu", (unsigned long long)(i+1));
        sprintf(str_idx_cell, "[%llu] :",
                (unsigned long long)(section->face_index[i] + 1));
      }
      else {
        str_num_cell[0] = '\0';
        str_idx_cell[0] = '\0';
      }
      sprintf(str_num_face, "%5u", (unsigned)(j-section->face_index[i]+1));
      sprintf(str_idx_face, "[%llu] :", (unsigned long long)(connect_length+1));

      bft_file_printf(f, "%12s %14s %5s %14s",
                      str_num_cell, str_idx_cell,
                      str_num_face, str_idx_face);

      /* Print face vertex numbers */

      if (section->face_num[j] > 0) {
        face_id = section->face_num[j] - 1;
        face_sgn = 1;
      }
      else {
        face_id = -section->face_num[j] - 1;
        face_sgn = -1;
      }

      line_values = 0;
      face_length = (  section->vertex_index[face_id+1]
                     - section->vertex_index[face_id]);
      connect_length += face_length;

      for (k = 0 ; k < face_length ; k++) {
        l =   section->vertex_index[face_id]
            + (face_length + (k*face_sgn))%face_length;
        if (line_values > 2) {
          line_values = 0;
          bft_file_printf(f,"\n%48s", "");
        }
        bft_file_printf(f, " %12llu",
                        (unsigned long long)section->vertex_num[l]);
        line_values++;
      }
      bft_file_printf(f, "\n");

    } /* End of loop on cell faces */

  } /* End of loop on polyhedral cells */

  str_num_cell[0] = '\0';
  sprintf(str_idx_cell, "[%llu] :", (unsigned long long)(j + 1));
  str_num_face[0] = '\0';
  sprintf(str_idx_face, "[%llu] :", (unsigned long long)(connect_length+1));
  bft_file_printf(f, "%12s %14s %5s %14s",
                  str_num_cell, str_idx_cell,
                  str_num_face, str_idx_face);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to a text file in parallel mode
 *
 * parameters:
 *   section                      <-- pointer to nodal mesh section structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   rank                         <-- rank in communicator
 *   n_ranks                      <-- number of processes in communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller       -- global connectivity slice provided
 *                                    by caller
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polygons_g(const fvm_nodal_section_t  *const section,
                         const fvm_io_num_t         *const global_vertex_num,
                         MPI_Comm                  comm,
                         const int                 rank,
                         const fvm_gnum_t          global_s_size,
                         const fvm_gnum_t          global_connect_s_size_caller,
                         fvm_gnum_t                global_connect_s_caller[],
                         bft_file_t         *const f)

{
  fvm_lnum_t  i, j;
  fvm_gnum_t  i_s, j_s;

  fvm_gnum_t  global_num_start;
  fvm_gnum_t  global_num_end;
  fvm_gnum_t  global_idx_shift = 0 ;

  fvm_gnum_t  *global_idx_s = NULL;

  fvm_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvm_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvm_gnum_t  *global_connect_s = global_connect_s_caller;

  fvm_gather_slice_t   *polygons_slice = NULL;

  const fvm_gnum_t   n_g_polygons = fvm_io_num_get_global_count
                                       (section->global_element_num);

  /* Allocate memory for additionnal indexes */

  BFT_MALLOC(global_idx_s, global_s_size + 1, fvm_gnum_t);

  /* Export by slices */

  polygons_slice = fvm_gather_slice_create(section->global_element_num,
                                           global_s_size,
                                           comm);

  while (fvm_gather_slice_advance(polygons_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather face->vertices index */

    fvm_gather_slice_index(section->vertex_index,
                           global_idx_s,
                           section->global_element_num,
                           comm,
                           polygons_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvm_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    comm,
                                    global_idx_s,
                                    polygons_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_caller)
        global_connect_s = NULL;
      BFT_REALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now gather face->vertices connectivity */

    fvm_gather_indexed_numbers(section->vertex_index,
                               section->vertex_num,
                               global_connect_s,
                               global_vertex_num,
                               section->global_element_num,
                               comm,
                               global_idx_s,
                               polygons_slice);

    /* Do all printing for faces on rank 0 */

    if (rank == 0) {

      int  line_values;
      char str_num_face[16];
      char str_idx_face[32];

      /* Print face connectivity */

      /* Loop on polygonal faces in slice */

      for (i = 0, i_s = global_num_start ; i_s < global_num_end ; i++, i_s++) {

        /* Print cell and face numbers and indexes */

        sprintf(str_num_face, "%12llu", (unsigned long long)i_s);
        sprintf(str_idx_face, "[%llu] :",
                (unsigned long long)(global_idx_s[i] + global_idx_shift + 1));

        bft_file_printf(f, "%12s %14s",
                        str_num_face, str_idx_face);

        /* Print face vertex numbers */

        line_values = 0;
        for (j = 0, j_s = global_idx_s[i] ;
             j_s < global_idx_s[i+1] ;
             j++, j_s++) {
          if (line_values > 2) {
            line_values = 0;
            bft_file_printf(f,"\n%27s", "");
          }
          bft_file_printf(f, " %12llu",
                          (unsigned long long)global_connect_s[j_s]);
          line_values++;
        }
        bft_file_printf(f, "\n");

      } /* End of loop on polyhedral cells in slice */

      if (global_num_end > n_g_polygons) {
        str_num_face[0] = '\0';
        sprintf(str_idx_face, "[%llu] :",
                (unsigned long long)(global_idx_s[i] + global_idx_shift + 1));
        bft_file_printf(f, "%12s %14s",
                        str_num_face, str_idx_face);
      }

      /* Update shift for conversion from slice index to true index */

      global_idx_shift
        += global_idx_s[global_num_end - global_num_start + 1];

    } /* End of printing for rank 0 for this slice */

  }

  fvm_gather_slice_destroy(polygons_slice);

  /* Free memory */

  BFT_FREE(global_idx_s);

  if (global_connect_s != global_connect_s_caller)
    BFT_FREE(global_connect_s);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to a text file in serial mode
 *
 * parameters:
 *   section      <-- pointer to nodal mesh section structure
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polygons_l(const fvm_nodal_section_t  *const section,
                         bft_file_t                 *const f)

{
  fvm_lnum_t  i, j;

  int  line_values;
  char str_num_face[16];
  char str_idx_face[32];

  /* Print face connectivity directly, without using extra buffers */

  j = 0; /* Initialize here in case section->n_elements = 0 */

  /* Loop on all polygonal faces */

  for (i = 0 ; i < section->n_elements ; i++) {

    /* Print face numbers and indexes */

    sprintf(str_num_face, "%12llu", (unsigned long long)(i+1));
    sprintf(str_idx_face, "[%llu] :",
            (unsigned long long)(section->vertex_index[i] + 1));

    bft_file_printf(f, "%12s %14s",
                    str_num_face, str_idx_face);

    /* Print face vertex numbers */

    line_values = 0;

    for (j = section->vertex_index[i] ;
         j < section->vertex_index[i+1] ;
         j++) {
      if (line_values > 2) {
        line_values = 0;
        bft_file_printf(f,"\n%27s", "");
      }
      bft_file_printf(f, " %12llu",
                      (unsigned long long)section->vertex_num[j]);
      line_values++;
    }
    bft_file_printf(f, "\n");

  } /* End of loop on polygonal faces */

  str_num_face[0] = '\0';
  sprintf(str_idx_face, "[%llu] :", (unsigned long long)(j + 1));
  bft_file_printf(f, "%12s %14s",
                  str_num_face, str_idx_face);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to text file writer.
 *
 * parameters:
 *   name    <-- base output case name.
 *   options <-- whitespace separated, lowercase options list
 *   comm    <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque text file writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_text_init_writer(const char                   *const name,
                        const char                   *const path,
                        const char                   *const options,
                        const fvm_writer_time_dep_t         time_dependency,
                        const MPI_Comm                      comm)
#else
void *
fvm_to_text_init_writer(const char                   *const name,
                        const char                   *const path,
                        const char                   *const options,
                        const fvm_writer_time_dep_t         time_dependency)
#endif

{
  fvm_to_text_writer_t  *this_writer = NULL;
  int  rank = 0;

  /* Initialize writer */

  BFT_MALLOC(this_writer, 1, fvm_to_text_writer_t);

  this_writer->time_dependency = time_dependency;

  this_writer->rank = 0;
  this_writer->n_ranks = 1;

#if defined(HAVE_MPI)
  {
    int mpi_flag, n_ranks;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag) {
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

  if (rank == 0) {

    char * file_name;
    int path_len = 0;
    const char extension[] = ".txt";

    if (path != NULL)
      path_len = strlen(path);

    BFT_MALLOC(file_name,
               path_len + strlen(name) + strlen(extension) + 1,
               char);

    if (path != NULL)
      strcpy(file_name, path);
    else
      file_name[0] = '\0';

    strcat(file_name, name);
    strcat(file_name, extension);
    this_writer->file = bft_file_open(file_name,
                                      BFT_FILE_MODE_WRITE,
                                      BFT_FILE_TYPE_TEXT);
    BFT_FREE(file_name);

  }
  else

    this_writer->file = NULL;

  /* Parse options */

  if (rank == 0 && options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1 ; i2 < l_tot && options[i2] != ' ' ; i2++);
      l_opt = i2 - i1 + 1;

      bft_file_printf(this_writer->file,
                      _("Option: %*s\n"), l_opt, options + i1);

      i1 = i2;

    }

  }

  /* Return writer */

  return this_writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to text file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque text file writer structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_text_finalize_writer(void  *this_writer_p)
{
  fvm_to_text_writer_t *this_writer = (fvm_to_text_writer_t *)this_writer_p;

  if (this_writer->file != NULL)
    this_writer->file = bft_file_free(this_writer->file);

  BFT_FREE(this_writer);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a text file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_text_export_nodal(void               *const this_writer_p,
                         const fvm_nodal_t  *const mesh)
{
  int         section_id;
  fvm_lnum_t  i, j;

  fvm_to_text_writer_t  *this_writer = (fvm_to_text_writer_t *)this_writer_p;
  bft_file_t  *f = NULL;

  fvm_gnum_t   global_connect_s_size, global_s_size;
  fvm_gnum_t   n_g_vertices = 0;
  fvm_gnum_t   n_g_edges = 0;
  fvm_gnum_t   n_g_faces = 0;
  fvm_gnum_t   n_g_cells = 0;
  fvm_gnum_t  * n_g_elements_section = NULL;

#if defined(HAVE_MPI)

  fvm_gnum_t  *global_connect_s = NULL;
  MPI_Comm    comm = this_writer->comm;

#endif

  const int  rank = this_writer->rank;
  const int  n_ranks = this_writer->n_ranks;

  /* Initialization */
  /*----------------*/

  f = this_writer->file;

  /* Buffer sizes required in parallel mode, global sizes always required */
  /*----------------------------------------------------------------------*/

  BFT_MALLOC(n_g_elements_section, mesh->n_sections, fvm_gnum_t);

  fvm_writer_def_nodal_buf_size(mesh,
                                n_ranks,
                                12,
                                5,
                                &n_g_vertices,
                                n_g_elements_section,
                                &global_s_size,
                                &global_connect_s_size);

  for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

    const fvm_nodal_section_t  *const  section = mesh->sections[section_id];

    switch(section->entity_dim) {
    case 1:
      n_g_edges += n_g_elements_section[section_id];
      break;
    case 2:
      n_g_faces += n_g_elements_section[section_id];
      break;
    case 3:
      n_g_cells += n_g_elements_section[section_id];
      break;
    default:
      assert(0);
    }

  }

  /* Global indicators */
  /*--------------------*/

  if (rank == 0) {

    if (mesh->name != NULL)
      bft_file_printf(f, _("\n"
                           "Mesh name: %s\n"),
                      mesh->name);
    else
      bft_file_printf(f, _("\n"
                           "Unnamed mesh\n"));

    bft_file_printf(f, _("\n"
                         "Mesh dimension:     %d\n"
                         "Number of domains:  %d\n"
                         "Number of sections:  %d\n"),
                    mesh->dim, mesh->n_doms, mesh->n_sections);

    bft_file_printf(f, _("\n"
                         "Number of cells:               %d\n"
                         "Number of faces:               %d\n"
                         "Number of edges:               %d\n"
                         "Number of vertices:            %d\n"),
                    n_g_cells,
                    n_g_faces,
                    n_g_edges,
                    n_g_vertices);

  }

  /* Vertex coordinates */
  /*--------------------*/

  {
    const int      stride = mesh->dim;
    const double  *local_coords;
    double        *coords_tmp = NULL;

    if (mesh->parent_vertex_num != NULL) {
      BFT_MALLOC(coords_tmp, stride * mesh->n_vertices, double);
      for (i = 0 ; i < mesh->n_vertices ; i++) {
        for (j = 0 ; j < stride ; j++)
          coords_tmp[i*stride + j]
            = mesh->vertex_coords[(mesh->parent_vertex_num[i]-1)*stride + j];
      }
      local_coords = coords_tmp;
    }
    else
      local_coords = mesh->vertex_coords;

    if (rank == 0)
      bft_file_printf(f, _("\nVertex coordinates:\n\n"));

    /* loop on slices in parallel mode, use whole array in serial mode */

#if defined(HAVE_MPI)

    if (n_ranks > 1) { /* start of output in parallel mode */

      fvm_gnum_t   global_num_start = 1;
      fvm_gnum_t   global_num_end = 0;

      fvm_gather_slice_t   *vertices_slice = NULL;
      double               *global_coords_s = NULL;

      BFT_MALLOC(global_coords_s, global_s_size * stride, double);

      vertices_slice = fvm_gather_slice_create(mesh->global_vertex_num,
                                               global_s_size,
                                               comm);

      while (fvm_gather_slice_advance(vertices_slice,
                                      &global_num_start,
                                      &global_num_end) == 0) {

        fvm_gather_array(local_coords,
                         global_coords_s,
                         MPI_DOUBLE,
                         (size_t)stride,
                         mesh->global_vertex_num,
                         comm,
                         vertices_slice);

        if (rank == 0)
          _write_slice_vector(stride,
                              global_num_start,
                              global_num_end,
                              global_coords_s,
                              f);

      }

      fvm_gather_slice_destroy(vertices_slice);

      BFT_FREE(global_coords_s);

    } /* end of output in parallel mode */

#endif /* defined(HAVE_MPI) */

    if (n_ranks == 1) { /* start of output in serial mode */

      _write_slice_vector(stride,
                          1,
                          (fvm_gnum_t)(mesh->n_vertices + 1),
                          local_coords,
                          f);

    } /* end of output in serial mode */

    if (coords_tmp != NULL)
      BFT_FREE(coords_tmp);
  }

  /* Allocate connectivity buffer for use wih all types of elements */

#if defined(HAVE_MPI)

  if (n_ranks > 1)
    BFT_MALLOC(global_connect_s, global_connect_s_size, fvm_gnum_t);

#endif

  /* Section connectivity */
  /*----------------------*/

  for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

    const fvm_nodal_section_t  *const  section = mesh->sections[section_id];

    if (rank == 0)
      bft_file_printf(f, _("\nSection: %s\n"
                           "  Number of elements: %lu\n\n"),
                      _(fvm_elements_type_name[section->type]),
                      (unsigned long long)(n_g_elements_section[section_id]));

    /* Output for strided (regular) element types */
    /*--------------------------------------------*/

    if (section->stride > 0) {

#if defined(HAVE_MPI)

      if (n_ranks > 1) { /* start of output in parallel mode */

        fvm_gnum_t   global_num_start;
        fvm_gnum_t   global_num_end;

        fvm_gather_slice_t   *elements_slice = NULL;

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
                                   f);

        }

        fvm_gather_slice_destroy(elements_slice);

      } /* end of output in parallel mode */

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1) { /* start of output in serial mode */

        _write_connect_l(section->stride,
                         section->n_elements,
                         section->vertex_num,
                         f);

      } /* end of output in serial mode */

    } /* end of output for strided element types */

    /* Output for polygons */
    /*---------------------*/

    else if (section->type == FVM_FACE_POLY) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1) {

        _export_nodal_polygons_g(section,
                                 mesh->global_vertex_num,
                                 comm,
                                 rank,
                                 global_s_size,
                                 global_connect_s_size,
                                 global_connect_s,
                                 f);

      }

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)

        _export_nodal_polygons_l(section, f);

    }

    /* Output for polyhedra */
    /*----------------------*/

    else if (section->type == FVM_CELL_POLY) {

#if defined(HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1) {

        _export_nodal_polyhedra_g(section,
                                  mesh->global_vertex_num,
                                  comm,
                                  rank,
                                  global_s_size,
                                  global_connect_s_size,
                                  global_connect_s,
                                  f);

      }

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)

        _export_nodal_polyhedra_l(section, f);

    }

  } /* End of loop on sections */

  /* Free buffers */
  /*--------------*/

  BFT_FREE(n_g_elements_section);

  /* Free buffers required in parallel mode */

#if defined(HAVE_MPI)

  if (n_ranks > 1)
    BFT_FREE(global_connect_s);

#endif /* defined(HAVE_MPI) */

  /* Close dump file */
  /*-----------------*/

  if (rank == 0)
    bft_file_flush(f);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
