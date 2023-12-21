/*============================================================================
 * Fortran interfaces of functions needing a synchronization of the extended
 * neighborhood.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_sort.h"
#include "fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ext_neighborhood.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_ext_neighborhood.c
        Extended cell neighborhood.

  \enum cs_ext_neighborhood_type_t

  \brief Type of extended neighborhood associated with the mesh

  Extended cell neighborhoods are based on cells adjacent through a shared
  vertex but with no shared face.

  \var CS_EXT_NEIGHBORHOOD_NONE
       No extended neighborhood

       \image html ext_neighborhood_none.svg "No extended neighborhood"

  \var CS_EXT_NEIGHBORHOOD_COMPLETE
       Full extended neighborhood

       \image html ext_neighborhood_complete.svg "Full extended neighborhood"

       This option should lead to the smoothest gradient, as it uses information
       from all neighbors, but is quite costly. On average, a hexahedral mesh
       has about 21 vertex neighbors per cell, a tetrahedral mesh
       around 150 vertex neighbors.

  \var CS_EXT_NEIGHBORHOOD_OPTIMIZED
       Optimized extended neighborhood based on heuristics.

       This option is similar to \ref CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE,
       keeping additional neighbors based on geometric tests and cell type.

       See the theory documentation for additional details.

  \var CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE
       Cells with centers best aligned opposite to face-adjacent cell centers

       \image html ext_neighborhood_cell_center_opposite.svg "Opposite cell centers"

       Add cells whose centers are closest to the half-line prolonging
       the [face-adjacent cell-center, cell center segment]. The number of
       additional cells in the extended neighborhood is at most equal to the
       number of cell faces.

  \var CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX
       Cells adjacent to faces whose non-orthogonality exceeds
       a given threshold (45 degrees by default)

       \image html ext_neighborhood_non_ortho_max.svg "Maximum non-orthogonality"

       Add cells adjacent to vertices of faces whose non-orthogonality exceeds
       a given threshold.

       This is the legacy reduced extended neighborhood option.

       Depending on the configuration, information may be heavily weighted to
       some sides of a cell and quite limited on other sides, which may explain
       why user feedback seems to dismiss this option.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* Vector size, in cs_real_t units (large enough for AVX512) */

#define CS_VS  8

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/* Short names for extended neighborhood types */

const char *cs_ext_neighborhood_type_name[]
= {N_("none"),
   N_("complete"),
   N_("optimized"),
   N_("opposite cell centers"),
   N_("non-orthogonality threshold")};

static cs_ext_neighborhood_type_t _ext_nbh_type
  = CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE;
static cs_real_t                  _non_ortho_max = 45;
static bool                       _full_nb_boundary = false;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]   n     size of array
 * \param[out]  s_id  start index for the current thread
 * \param[out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static void
_thread_range(cs_lnum_t   n,
              cs_lnum_t  *s_id,
              cs_lnum_t  *e_id)
{
#if defined(HAVE_OPENMP)
  int t_id = omp_get_thread_num();
  int n_t = omp_get_num_threads();
  cs_lnum_t t_n = (n + n_t - 1) / n_t;
  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, CS_CL);
  *e_id = cs_align(*e_id, CS_CL);
  if (*e_id > n) *e_id = n;
#else
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the square norm of the distance from a point to a line
 *        defined by a vector.
 *
 * \param[in]     u     vector of 3 real values aligned with line
 * \param[in]     v     vector of 3 real values from one pont one line
 *                      to given point
 *
 * \return square of distance, given by: |u x v| / |v|
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_cross_product_sq_norm(const cs_real_t u[3],
                       const cs_real_t v[3])
{
  cs_real_t uv[3];

  uv[0] = u[1]*v[2] - u[2]*v[1];
  uv[1] = u[2]*v[0] - u[0]*v[2];
  uv[2] = u[0]*v[1] - u[1]*v[0];

  cs_real_t d2n = uv[0]*uv[0] + uv[1]*uv[1] + uv[2]*uv[2];
  cs_real_t d2d = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];

  return d2n / d2d;
}

/*----------------------------------------------------------------------------
 * Extract a mesh's "cell -> internal faces" connectivity.
 *
 * parameters:
 *   mesh               <-- pointer to a cs_mesh_t structure
 *   p_cell_i_faces_idx --> pointer to the "cell -> faces" connectivity index
 *                          (0 to n-1 numbering)
 *   p_cell_i_faces_lst --> pointer to the "cell -> faces" connectivity list
 *                          (0 to n-1, no orientation)
 *----------------------------------------------------------------------------*/

static void
_get_cell_i_faces_connectivity(const cs_mesh_t          *mesh,
                               cs_lnum_t         **const p_cell_i_faces_idx,
                               cs_lnum_t         **const p_cell_i_faces_lst)
{

  cs_lnum_t  i, j, j1, j2;

  cs_lnum_t  *cell_faces_count = NULL;
  cs_lnum_t  *cell_faces_idx = NULL;
  cs_lnum_t  *cell_faces_lst = NULL;

  /* Allocate and initialize index */

  BFT_MALLOC(cell_faces_idx, mesh->n_cells + 1, cs_lnum_t);

  for (i = 0; i < mesh->n_cells + 1; i++)
    cell_faces_idx[i] = 0;

  /* Count number of faces per cell (we assign the temporary counter
   * to i and cell_faces_idx[i + 1] instead of cell_faces_idx[i]
   * to simplify the next stage) */

  /* Note: test if j < mesh->n_cells on internal faces to ignore
     parallel and/or periodic ghost cells */

  for (i = 0; i < mesh->n_i_faces; i++) {
    j1 = mesh->i_face_cells[i][0];
    j2 = mesh->i_face_cells[i][1];
    if (j1 < mesh->n_cells)
      cell_faces_idx[j1 + 1] += 1;
    if (j2 < mesh->n_cells)
      cell_faces_idx[j2 + 1] += 1;
  }

  /* Build position index */

  cell_faces_idx[0] = 0;
  for (j = 0; j < mesh->n_cells; j++)
    cell_faces_idx[j + 1] += cell_faces_idx[j];

  /* Build array of values */

  BFT_MALLOC(cell_faces_lst, cell_faces_idx[mesh->n_cells], cs_lnum_t);
  BFT_MALLOC(cell_faces_count, mesh->n_cells, cs_lnum_t);

  for (i = 0; i < mesh->n_cells; i++)
    cell_faces_count[i] = 0;

  for (i = 0; i < mesh->n_i_faces; i++) {
    j1 = mesh->i_face_cells[i][0];
    j2 = mesh->i_face_cells[i][1];
    if (j1 < mesh->n_cells) {
      cell_faces_lst[cell_faces_idx[j1] + cell_faces_count[j1]] = i;
      cell_faces_count[j1] += 1;
    }
    if (j2 < mesh->n_cells) {
      cell_faces_lst[cell_faces_idx[j2] + cell_faces_count[j2]] = i;
      cell_faces_count[j2] += 1;
    }
  }

  BFT_FREE(cell_faces_count);

  /* Set return values */

  *p_cell_i_faces_idx = cell_faces_idx;
  *p_cell_i_faces_lst = cell_faces_lst;
}

/*----------------------------------------------------------------------------
 * Extract a mesh's "cell -> boundary faces" connectivity.
 *
 * parameters:
 *   mesh               <-- pointer to a cs_mesh_t structure
 *   p_cell_b_faces_idx --> pointer to the "cell -> faces" connectivity index
 *                          (0 to n-1 numbering)
 *   p_cell_b_faces_lst --> pointer to the "cell -> faces" connectivity list
 *                          (0 to n-1, no orientation)
 *----------------------------------------------------------------------------*/

static void
_get_cell_b_faces_connectivity(const cs_mesh_t          *mesh,
                               cs_lnum_t         **const p_cell_b_faces_idx,
                               cs_lnum_t         **const p_cell_b_faces_lst)
{

  cs_lnum_t  i, j, j1;

  cs_lnum_t  *cell_faces_count = NULL;
  cs_lnum_t  *cell_faces_idx = NULL;
  cs_lnum_t  *cell_faces_lst = NULL;

  /* Allocate and initialize index */

  BFT_MALLOC(cell_faces_idx, mesh->n_cells + 1, cs_lnum_t);

  for (i = 0; i < mesh->n_cells + 1; i++)
    cell_faces_idx[i] = 0;

  /* Count number of faces per cell (we assign the temporary counter
   * to i and cell_faces_idx[i + 1] instead of cell_faces_idx[i]
   * to simplify the next stage) */

  /* Note: test if j < mesh->n_cells on internal faces to ignore
     parallel and/or periodic ghost cells */

  for (i = 0; i < mesh->n_b_faces; i++) {
    j1 = mesh->b_face_cells[i];
    if (j1 < mesh->n_cells)
      cell_faces_idx[j1 + 1] += 1;
  }

  /* Build position index */

  cell_faces_idx[0] = 0;
  for (j = 0; j < mesh->n_cells; j++)
    cell_faces_idx[j + 1] += cell_faces_idx[j];

  /* Build array of values */

  BFT_MALLOC(cell_faces_lst, cell_faces_idx[mesh->n_cells], cs_lnum_t);
  BFT_MALLOC(cell_faces_count, mesh->n_cells, cs_lnum_t);

  for (i = 0; i < mesh->n_cells; i++)
    cell_faces_count[i] = 0;

  for (i = 0; i < mesh->n_b_faces; i++) {
    j1 = mesh->b_face_cells[i];
    if (j1 < mesh->n_cells) {
      cell_faces_lst[cell_faces_idx[j1] + cell_faces_count[j1]] = i;
      cell_faces_count[j1] += 1;
    }
  }

  BFT_FREE(cell_faces_count);

  /* Set return values */

  *p_cell_b_faces_idx = cell_faces_idx;
  *p_cell_b_faces_lst = cell_faces_lst;
}

/*----------------------------------------------------------------------------
 * Create a "vertex -> cells" connectivity.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_vtx_cells_idx --> pointer to the "vtx -> cells" connectivity index
 *   p_vtx_cells_lst --> pointer to the "vtx -> cells" connectivity list
 *----------------------------------------------------------------------------*/

static void
_create_vtx_cells_connect(cs_mesh_t  *mesh,
                          cs_lnum_t  *p_vtx_cells_idx[],
                          cs_lnum_t  *p_vtx_cells_lst[])
{
  cs_lnum_t  i, j, idx;
  cs_lnum_t  vtx_id, face_id, cell_id;

  bool      already_seen;

  cs_lnum_t  vtx_cells_connect_size = 0;

  cs_lnum_t  *vtx_faces_idx = NULL, *vtx_faces_lst = NULL;
  cs_lnum_t  *vtx_cells_idx = NULL, *vtx_cells_lst = NULL;

  const cs_lnum_t  n_vertices = mesh->n_vertices;
  const cs_lnum_t  n_faces = mesh->n_i_faces;
  const cs_lnum_t  *face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t  *face_vtx_lst = mesh->i_face_vtx_lst;
  const cs_lnum_2_t  *face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);

  cs_lnum_t  vtx_cells_estimated_connect_size = 3 * n_vertices;

  BFT_MALLOC(vtx_cells_idx, n_vertices + 1, cs_lnum_t);
  BFT_MALLOC(vtx_faces_idx, n_vertices + 1, cs_lnum_t);

  for (vtx_id = 0; vtx_id < n_vertices + 1; vtx_id++) {
    vtx_cells_idx[vtx_id] = 0;
    vtx_faces_idx[vtx_id] = 0;
  }

  /* Define vtx -> faces connectivity index */

  for (face_id = 0; face_id < n_faces; face_id++) {

    for (i = face_vtx_idx[face_id]; i < face_vtx_idx[face_id+1]; i++) {
      vtx_id = face_vtx_lst[i];
      vtx_faces_idx[vtx_id + 1] += 1;
    }

  } /* End of loop on faces */

  vtx_faces_idx[0] = 0;
  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++)
    vtx_faces_idx[vtx_id + 1] += vtx_faces_idx[vtx_id];

  /* Allocation and definiton of "vtx -> faces" connectivity list */

  BFT_MALLOC(vtx_faces_lst, vtx_faces_idx[n_vertices], cs_lnum_t);

  for (face_id = 0; face_id < n_faces; face_id++) {

    for (i = face_vtx_idx[face_id]; i < face_vtx_idx[face_id+1]; i++) {

      vtx_id = face_vtx_lst[i];
      vtx_faces_lst[vtx_faces_idx[vtx_id]] = face_id;
      vtx_faces_idx[vtx_id] += 1;

    }

  } /* End of loop on faces */

  for (vtx_id = n_vertices; vtx_id > 0; vtx_id--)
    vtx_faces_idx[vtx_id] = vtx_faces_idx[vtx_id-1];
  vtx_faces_idx[0] = 0;

  /* Define "vertex -> cells" connectivity.
     Use "vertex -> faces" connectivity and "face -> cells" connectivity */

  BFT_MALLOC(vtx_cells_lst, vtx_cells_estimated_connect_size, cs_lnum_t);

  vtx_cells_idx[0] = 0;

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++) {

    for (i = vtx_faces_idx[vtx_id]; i < vtx_faces_idx[vtx_id+1]; i++) {

      face_id = vtx_faces_lst[i];

      for (j = 0; j < 2; j++) { /* For the cells sharing this face */

        cell_id = face_cells[face_id][j];

        already_seen = false;
        idx = vtx_cells_idx[vtx_id];

        while ((already_seen == false) && (idx < vtx_cells_connect_size)) {
          if (cell_id == vtx_cells_lst[idx])
            already_seen = true;
          idx++;
        }

        if (already_seen == false) {

          if (vtx_cells_connect_size + 1 > vtx_cells_estimated_connect_size) {
            vtx_cells_estimated_connect_size *= 2;
            BFT_REALLOC(vtx_cells_lst,
                        vtx_cells_estimated_connect_size, cs_lnum_t);
          }

          vtx_cells_lst[vtx_cells_connect_size] = cell_id;
          vtx_cells_connect_size += 1;

        }

      } /* End of loop on cells sharing this face */

    } /* End of loop on faces sharing this vertex */

    vtx_cells_idx[vtx_id+1] = vtx_cells_connect_size;

  } /* End of loop on vertices */

  BFT_REALLOC(vtx_cells_lst, vtx_cells_connect_size, cs_lnum_t);

  /* Free memory */

  BFT_FREE(vtx_faces_idx);
  BFT_FREE(vtx_faces_lst);

  *p_vtx_cells_idx = vtx_cells_idx;
  *p_vtx_cells_lst = vtx_cells_lst;

}

/*----------------------------------------------------------------------------
 * Create a "vertex -> cells" connectivity.
 *
 * parameters:
 *   face_id        <-- identification number for the face
 *   cell_id        <-- identification number for the cell sharing this face
 *   mesh           <-- pointer to a cs_mesh_t structure
 *   cell_cells_tag <-- pointer to a tag array
 *   vtx_cells_idx  --> pointer to the "vtx -> cells" connectivity index
 *   vtx_cells_lst  --> pointer to the "vtx -> cells" connectivity list
 *   vtx_gcells_idx --> pointer to the "vtx -> gcells" connectivity index
 *   vtx_gcells_lst --> pointer to the "vtx -> gcells" connectivity list
 *----------------------------------------------------------------------------*/

static void
_tag_cells(cs_lnum_t         face_id,
           cs_lnum_t         cell_id,
           const cs_mesh_t  *mesh,
           char              cell_cells_tag[],
           cs_lnum_t         vtx_cells_idx[],
           cs_lnum_t         vtx_cells_lst[],
           cs_lnum_t         vtx_gcells_idx[],
           cs_lnum_t         vtx_gcells_lst[])
{
  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;

  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  *face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t  *face_vtx_lst = mesh->i_face_vtx_lst;

  if (cell_id < n_cells) {

    for (cs_lnum_t i = cell_cells_idx[cell_id];
         i < cell_cells_idx[cell_id+1];
         i++) {

      /* For cell not tagged yet */

      if (cell_cells_tag[i] == 0) {

        cs_lnum_t ext_cell_id = cell_cells_lst[i];

        /* Cells sharing a vertex with the face */

        for (cs_lnum_t j = face_vtx_idx[face_id];
             j < face_vtx_idx[face_id+1];
             j++) {

          cs_lnum_t vtx_id = face_vtx_lst[j];

          if (ext_cell_id < n_cells) { /* true cell */

            for (cs_lnum_t k = vtx_cells_idx[vtx_id];
                 k < vtx_cells_idx[vtx_id+1];
                 k++) {

              cs_lnum_t loc_cell_id = vtx_cells_lst[k];

              /* Comparison and selection */

              if (loc_cell_id == ext_cell_id && cell_cells_tag[i] == 0)
                cell_cells_tag[i] = 1;

            } /* End of loop on cells sharing this vertex */

          }
          else { /* ext_cell_num >= n_cells; i.e. ghost cell */

            for (cs_lnum_t k = vtx_gcells_idx[vtx_id];
                 k < vtx_gcells_idx[vtx_id+1];
                 k++) {

              cs_lnum_t loc_cell_id = vtx_gcells_lst[k] + n_cells;

              /* Comparison and selection */

              if (loc_cell_id == ext_cell_id && cell_cells_tag[i] == 0)
                cell_cells_tag[i] = 1;

            }

          }

        } /* End of loop on vertices of the face */

      }

    } /* End of loop on cells in the extended neighborhood */

  } /* End if cell_id < n_cells */

}

/*---------------------------------------------------------------------------
 * Reverse "ghost cell -> vertex" connectivity into "vertex -> ghost cells"
 * connectivity for halo elements.
 * Build the connectivity index.
 *
 * parameters:
 *   halo            <-- pointer to a cs_halo_t structure
 *   rank_id         <-- rank number to work with
 *   checker         <-> temporary array to check vertices
 *   gcell_vtx_idx   <-- "ghost cell -> vertices" connectivity index
 *   gcell_vtx_lst   <-- "ghost cell -> vertices" connectivity list
 *   vtx_gcells_idx  <-> "vertex -> ghost cells" connectivity index
 *---------------------------------------------------------------------------*/

static void
_reverse_connectivity_idx(cs_halo_t   *halo,
                          cs_lnum_t    n_vertices,
                          cs_lnum_t    rank_id,
                          cs_lnum_t   *checker,
                          cs_lnum_t   *gcell_vtx_idx,
                          cs_lnum_t   *gcell_vtx_lst,
                          cs_lnum_t   *vtx_gcells_idx)
{
  cs_lnum_t  i, j, id, vtx_id, start_idx, end_idx;

  /* Initialize index */

  vtx_gcells_idx[0] = 0;
  for (i = 0; i < n_vertices; i++) {
    vtx_gcells_idx[i+1] = 0;
    checker[i] = -1;
  }

  if (rank_id == -1) {
    start_idx = 0;
    end_idx = halo->n_elts[CS_HALO_EXTENDED];
  }
  else { /* Call with rank_id > 1 for standard halo */
    start_idx = halo->index[2*rank_id];
    end_idx = halo->index[2*rank_id+1];
  }

  /* Define index */

  for (id = start_idx; id < end_idx; id++) {

    for (j = gcell_vtx_idx[id]; j < gcell_vtx_idx[id+1]; j++) {

      vtx_id = gcell_vtx_lst[j];

      if (checker[vtx_id] != id) {
        checker[vtx_id] = id;
        vtx_gcells_idx[vtx_id+1] += 1;
      }

    }

  } /* End of loop of ghost cells */

  for (i = 0; i < n_vertices; i++)
    vtx_gcells_idx[i+1] += vtx_gcells_idx[i];
}

/*---------------------------------------------------------------------------
 * Reverse "ghost cells -> vertex" connectivity into "vertex -> ghost cells"
 * connectivity for halo elements.
 * Build the connectivity list.
 *
 * parameters:
 *   halo            <-- pointer to a cs_halo_t structure
 *   n_vertices      <-- number of vertices
 *   rank_id         <-- rank number to work with
 *   counter         <-> temporary array to count vertices
 *   checker         <-> temporary array to check vertices
 *   gcell_vtx_idx   <-- "ghost cell -> vertices" connectivity index
 *   gcell_vtx_lst   <-- "ghost cell -> vertices" connectivity list
 *   vtx_gcells_idx  <-- "vertex -> ghost cells" connectivity index
 *   vtx_gcells_lst  <-> "vertex -> ghost cells" connectivity list
 *---------------------------------------------------------------------------*/

static void
_reverse_connectivity_lst(cs_halo_t   *halo,
                          cs_lnum_t    n_vertices,
                          cs_lnum_t    rank_id,
                          cs_lnum_t   *counter,
                          cs_lnum_t   *checker,
                          cs_lnum_t   *gcell_vtx_idx,
                          cs_lnum_t   *gcell_vtx_lst,
                          cs_lnum_t   *vtx_gcells_idx,
                          cs_lnum_t   *vtx_gcells_lst)
{
  cs_lnum_t  i, j, id, shift, vtx_id, start_idx, end_idx;

  /* Initialize buffers */

  for (i = 0; i < n_vertices; i++) {
    counter[i] = 0;
    checker[i] = -1;
  }

  if (rank_id == -1) {
    start_idx = 0;
    end_idx = halo->n_elts[CS_HALO_EXTENDED];
  }
  else {
    start_idx = halo->index[2*rank_id];
    end_idx = halo->index[2*rank_id+1];
  }

  /* Fill the connectivity list */

  for (id = start_idx; id < end_idx; id++) {

    for (j = gcell_vtx_idx[id]; j < gcell_vtx_idx[id+1]; j++) {

      vtx_id = gcell_vtx_lst[j];

      if (checker[vtx_id] != id) {

        checker[vtx_id] = id;
        shift = vtx_gcells_idx[vtx_id] + counter[vtx_id];
        vtx_gcells_lst[shift] = id;
        counter[vtx_id] += 1;

      }

    }

  } /* End of loop of ghost cells */

}

/*---------------------------------------------------------------------------
 * Create a "vertex -> ghost cells" connectivity.
 * Add mesh->n_cells to all the elements of vtx_gcells_lst to obtain
 * the local cell numbering.
 *
 * parameters:
 *   halo              <-- pointer to a cs_halo_t structure
 *   n_vertices        <-- number of vertices
 *   gcell_vtx_idx     <-- "ghost cell -> vertices" connectivity index
 *   gcell_vtx_lst     <-- "ghost cell -> vertices" connectivity list
 *   p_vtx_gcells_idx  --> pointer to "vertex -> ghost cells" index
 *   p_vtx_gcells_lst  --> pointer to "vertex -> ghost cells" list
 *---------------------------------------------------------------------------*/

static void
_create_vtx_gcells_connect(cs_halo_t    *halo,
                           cs_lnum_t     n_vertices,
                           cs_lnum_t    *gcells_vtx_idx,
                           cs_lnum_t    *gcells_vtx_lst,
                           cs_lnum_t    *p_vtx_gcells_idx[],
                           cs_lnum_t    *p_vtx_gcells_lst[])
{
  cs_lnum_t  *vtx_buffer = NULL, *vtx_counter = NULL, *vtx_checker = NULL;
  cs_lnum_t  *vtx_gcells_idx = NULL, *vtx_gcells_lst = NULL;

  BFT_MALLOC(vtx_buffer, 2*n_vertices, cs_lnum_t);
  vtx_counter = &(vtx_buffer[0]);
  vtx_checker = &(vtx_buffer[n_vertices]);

  BFT_MALLOC(vtx_gcells_idx, n_vertices + 1, cs_lnum_t);

  /* Create a vertex -> ghost cells connectivity */

  _reverse_connectivity_idx(halo,
                            n_vertices,
                            -1,
                            vtx_checker,
                            gcells_vtx_idx,
                            gcells_vtx_lst,
                            vtx_gcells_idx);

  BFT_MALLOC(vtx_gcells_lst, vtx_gcells_idx[n_vertices], cs_lnum_t);

  _reverse_connectivity_lst(halo,
                            n_vertices,
                            -1,
                            vtx_counter,
                            vtx_checker,
                            gcells_vtx_idx,
                            gcells_vtx_lst,
                            vtx_gcells_idx,
                            vtx_gcells_lst);

  *p_vtx_gcells_idx = vtx_gcells_idx;
  *p_vtx_gcells_lst = vtx_gcells_lst;

  /* Free memory */

  BFT_FREE(vtx_buffer);

}

/*---------------------------------------------------------------------------
 * Create a "vertex -> cells" connectivity.
 *
 * parameters:
 *   mesh              <-- pointer to cs_mesh_t structure
 *   cell_i_faces_idx  <-- "cell -> internal faces" connectivity index
 *   cell_i_faces_lst  <-- "cell -> internal faces" connectivity list
 *   p_vtx_cells_idx   --> pointer to "vertex -> cells" connectivity index
 *   p_vtx_cells_lst   --> pointer to "vertex -> cells" connectivity list
 *---------------------------------------------------------------------------*/

static void
_create_vtx_cells_connect2(cs_mesh_t   *mesh,
                           cs_lnum_t   *cell_i_faces_idx,
                           cs_lnum_t   *cell_i_faces_lst,
                           cs_lnum_t   *p_vtx_cells_idx[],
                           cs_lnum_t   *p_vtx_cells_lst[])
{
  cs_lnum_t  i, cell_id, fac_id, i_vtx;
  cs_lnum_t  shift, vtx_id;

  cs_lnum_t  *vtx_buffer = NULL, *vtx_count = NULL, *vtx_tag = NULL;
  cs_lnum_t  *vtx_cells_idx = NULL, *vtx_cells_lst = NULL;

  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_vertices = mesh->n_vertices;
  const cs_lnum_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  /* Initialize buffers */

  BFT_MALLOC(vtx_buffer, 2*n_vertices, cs_lnum_t);
  BFT_MALLOC(vtx_cells_idx, n_vertices + 1, cs_lnum_t);

  vtx_count = &(vtx_buffer[0]);
  vtx_tag = &(vtx_buffer[n_vertices]);

  vtx_cells_idx[0] = 0;
  for (i = 0; i < n_vertices; i++) {

    vtx_cells_idx[i + 1] = 0;
    vtx_tag[i] = -1;
    vtx_count[i] = 0;

  }

  /* Define index */

  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    for (i = cell_i_faces_idx[cell_id]; i < cell_i_faces_idx[cell_id+1]; i++) {

      fac_id = cell_i_faces_lst[i];

      for (i_vtx = fac_vtx_idx[fac_id];
           i_vtx < fac_vtx_idx[fac_id+1];
           i_vtx++) {

        vtx_id = fac_vtx_lst[i_vtx];

        if (vtx_tag[vtx_id] != cell_id) {

          vtx_cells_idx[vtx_id+1] += 1;
          vtx_tag[vtx_id] = cell_id;

        } /* Add this cell to the connectivity of the vertex */

      } /* End of loop on vertices */

    } /* End of loop on cell's faces */

  } /* End of loop on cells */

  for (i = 0; i < n_vertices; i++) {

    vtx_cells_idx[i+1] += vtx_cells_idx[i];
    vtx_tag[i] = -1;

  }

  BFT_MALLOC(vtx_cells_lst, vtx_cells_idx[n_vertices], cs_lnum_t);

  /* Fill list */

  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    for (i = cell_i_faces_idx[cell_id]; i < cell_i_faces_idx[cell_id+1]; i++) {

      fac_id = cell_i_faces_lst[i];

      for (i_vtx = fac_vtx_idx[fac_id];
           i_vtx < fac_vtx_idx[fac_id+1];
           i_vtx++) {

        vtx_id = fac_vtx_lst[i_vtx];

        if (vtx_tag[vtx_id] != cell_id) {

          shift = vtx_cells_idx[vtx_id] + vtx_count[vtx_id];
          vtx_tag[vtx_id] = cell_id;
          vtx_cells_lst[shift] = cell_id;
          vtx_count[vtx_id] += 1;

        } /* Add this cell to the connectivity of the vertex */

      } /* End of loop on vertices */

    } /* End of loop on cell's faces */

  } /* End of loop on cells */

  BFT_FREE(vtx_buffer);

  *p_vtx_cells_idx = vtx_cells_idx;
  *p_vtx_cells_lst = vtx_cells_lst;

}

/*---------------------------------------------------------------------------
 * Create a "cell -> cells" connectivity.
 *
 * parameters:
 *   mesh              <-- pointer to cs_mesh_t structure
 *   cell_i_faces_idx  <-- "cell -> faces" connectivity index
 *   cell_i_faces_lst  <-- "cell -> faces" connectivity list
 *   vtx_gcells_idx    --> "vertex -> ghost cells" connectivity index
 *   vtx_gcells_lst    --> "vertex -> ghost cells" connectivity list
 *   vtx_cells_idx     --> "vertex -> cells" connectivity index
 *   vtx_cells_lst     --> "vertex -> cells" connectivity list
 *   p_cell_cells_idx  --> pointer to "cell -> cells" connectivity index
 *   p_cell_cells_lst  --> pointer to "cell -> cells" connectivity list
 *---------------------------------------------------------------------------*/

static void
_create_cell_cells_connect(cs_mesh_t  *mesh,
                           cs_lnum_t  *cell_i_faces_idx,
                           cs_lnum_t  *cell_i_faces_lst,
                           cs_lnum_t  *vtx_gcells_idx,
                           cs_lnum_t  *vtx_gcells_lst,
                           cs_lnum_t  *vtx_cells_idx,
                           cs_lnum_t  *vtx_cells_lst,
                           cs_lnum_t  *p_cell_cells_idx[],
                           cs_lnum_t  *p_cell_cells_lst[])
{
  cs_lnum_t  i, j, i_cel, fac_id, i_vtx;
  cs_lnum_t  cell_id, vtx_id, shift;

  cs_lnum_t  *cell_buffer = NULL, *cell_tag = NULL, *cell_count = NULL;
  cs_lnum_t  *cell_cells_idx = NULL, *cell_cells_lst = NULL;

  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;
  const cs_lnum_2_t  *face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const cs_lnum_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  /* Allocate and initialize buffers */

  BFT_MALLOC(cell_cells_idx, n_cells + 1, cs_lnum_t);
  BFT_MALLOC(cell_buffer, n_cells_wghosts + n_cells, cs_lnum_t);

  cell_tag = &(cell_buffer[0]);
  cell_count = &(cell_buffer[n_cells_wghosts]);

  cell_cells_idx[0] = 0;
  for (i = 0; i < n_cells; i++) {
    cell_cells_idx[i+1] = 0;
    cell_count[i] = 0;
  }

  for (i = 0; i < n_cells_wghosts; i++)
    cell_tag[i] = -1;

  /* Define index */

  for (i_cel = 0; i_cel < n_cells; i_cel++) {

    /* First loop on faces to tag cells sharing a face */

    for (i = cell_i_faces_idx[i_cel]; i < cell_i_faces_idx[i_cel+1]; i++) {

      fac_id = cell_i_faces_lst[i];

      cell_tag[face_cells[fac_id][0]] = i_cel;
      cell_tag[face_cells[fac_id][1]] = i_cel;

    } /* End of loop on cell's faces */

    /* Second loop on faces to update index */

    for (i = cell_i_faces_idx[i_cel]; i < cell_i_faces_idx[i_cel+1]; i++) {

      fac_id = cell_i_faces_lst[i];

      for (i_vtx = fac_vtx_idx[fac_id];
           i_vtx < fac_vtx_idx[fac_id+1];
           i_vtx++) {

        vtx_id = fac_vtx_lst[i_vtx];

        /* For cells belonging to this rank, get vertex -> cells connect. */

        for (j = vtx_cells_idx[vtx_id]; j < vtx_cells_idx[vtx_id+1]; j++) {

          cell_id = vtx_cells_lst[j];

          if (cell_tag[cell_id] != i_cel) {
            cell_cells_idx[i_cel+1] += 1;
            cell_tag[cell_id] = i_cel;
          }

        }

        if (n_cells_wghosts - n_cells > 0) { /* If there are ghost cells */

          /* For ghost cells, get vertex -> ghost cells connect. */

          for (j = vtx_gcells_idx[vtx_id]; j < vtx_gcells_idx[vtx_id+1]; j++) {

            cell_id = vtx_gcells_lst[j] + n_cells;

            if (cell_tag[cell_id] != i_cel) {
              cell_cells_idx[i_cel+1] += 1;
              cell_tag[cell_id] = i_cel;
            }

          }

        } /* If there are ghost cells */

      } /* End of loop on vertices */

    } /* End of loop on cell's faces */

  } /* End of loop on cells */

  /* Create index */

  for (i = 0; i < n_cells; i++)
    cell_cells_idx[i+1] += cell_cells_idx[i];

  for (i = 0; i < n_cells_wghosts; i++)
    cell_tag[i] = -1;

  BFT_MALLOC(cell_cells_lst, cell_cells_idx[n_cells], cs_lnum_t);

  /* Fill list */

  for (i_cel = 0; i_cel < n_cells; i_cel++) {

    /* First loop on faces to tag cells sharing a face */

    for (i = cell_i_faces_idx[i_cel]; i < cell_i_faces_idx[i_cel+1]; i++) {

      fac_id = cell_i_faces_lst[i];

      cell_tag[face_cells[fac_id][0]] = i_cel;
      cell_tag[face_cells[fac_id][1]] = i_cel;

    } /* End of loop on cell's faces */

    /* Second loop on faces to update list */

    for (i = cell_i_faces_idx[i_cel]; i < cell_i_faces_idx[i_cel+1]; i++) {

      fac_id = cell_i_faces_lst[i];

      for (i_vtx = fac_vtx_idx[fac_id];
           i_vtx < fac_vtx_idx[fac_id+1];
           i_vtx++) {

        vtx_id = fac_vtx_lst[i_vtx];

        /* For cells belonging to this rank, get vertex -> cells connect. */

        for (j = vtx_cells_idx[vtx_id]; j < vtx_cells_idx[vtx_id+1]; j++) {

          cell_id = vtx_cells_lst[j];

          if (cell_tag[cell_id] != i_cel) {

            shift = cell_cells_idx[i_cel] + cell_count[i_cel];
            cell_cells_lst[shift] = cell_id;
            cell_tag[cell_id] = i_cel;
            cell_count[i_cel] += 1;

          } /* Add this cell_id */

        }

        if (n_cells_wghosts - n_cells > 0) { /* If there are ghost cells */

          /* For ghost cells, get vertex -> ghost cells connect. */

          for (j = vtx_gcells_idx[vtx_id]; j < vtx_gcells_idx[vtx_id+1]; j++) {

            cell_id = vtx_gcells_lst[j] + n_cells;

            if (cell_tag[cell_id] != i_cel) {

              shift = cell_cells_idx[i_cel] + cell_count[i_cel];
              cell_cells_lst[shift] = cell_id;
              cell_tag[cell_id] = i_cel;
              cell_count[i_cel] += 1;

            }

          }

        } /* If there are ghost cells */

      } /* End of loop on vertices */

    } /* End of loop on cell's faces */

  } /* End of loop on cells */

  /* Sort line elements by column id (for better access patterns) */

  bool unique = cs_sort_indexed(n_cells, cell_cells_idx, cell_cells_lst);

  assert(unique == true);

  *p_cell_cells_idx = cell_cells_idx;
  *p_cell_cells_lst = cell_cells_lst;

  /* Free memory */

  BFT_FREE(cell_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reduce the "cell -> cells" connectivity for the
 *        extended neighborhood using a non-orthogonality criterion.
 *
 * Note: Only cells sharing only a vertex or vertices (not a face)
 *       belong to the "cell -> cells" connectivity.
 *
 * \param[in]       mesh             pointer to mesh structure
 * \param[in]       mesh_quantities  associated mesh quantities
 * \param[in, out]  cell_cells_tag   tag adjacencies to retain
 */
/*----------------------------------------------------------------------------*/

static void
_neighborhood_reduce_anomax(cs_mesh_t             *mesh,
                            cs_mesh_quantities_t  *mesh_quantities,
                            char                   cell_cells_tag[])
{
  cs_lnum_t  i, face_id, cell_i, cell_j;

  cs_real_t  v_ij[3];
  cs_real_t  face_normal[3];
  cs_real_t  norm_ij, face_norm, cos_ij_fn;
  cs_real_t  dprod;

  cs_lnum_t  *vtx_cells_idx = NULL, *vtx_cells_lst = NULL;
  cs_lnum_t  *vtx_gcells_idx = NULL, *vtx_gcells_lst = NULL;

  /* Convert degrees to radians */
  const cs_real_t non_ortho_max = _non_ortho_max * cs_math_pi / 180;

  const cs_lnum_t  n_faces = mesh->n_i_faces;
  const cs_lnum_2_t  *face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);

  const cs_real_t  cos_ij_fn_min = cos(non_ortho_max);
  const cs_real_t  *cell_cen = mesh_quantities->cell_cen;

  enum {X, Y, Z};

  assert(mesh->dim == 3);

  if (non_ortho_max <= 0)
    return;

  /*
    For each internal face, we select in the extended neighborhood
    of the two cells sharing this face all the cells sharing a
    vertex of this face if the non-orthogonality angle is above a
    criterion.
  */

  /*
    First: re-build a "vertex -> cells" connectivity
    ------------------------------------------------
    We have to invert the "face -> vertices" connectivity and then
    we will use the "face -> cells" conectivity.
  */

  _create_vtx_cells_connect(mesh,
                            &vtx_cells_idx,
                            &vtx_cells_lst);

  if (cs_mesh_n_g_ghost_cells(mesh) > 0)
    _create_vtx_gcells_connect(mesh->halo,
                               mesh->n_vertices,
                               mesh->gcell_vtx_idx,
                               mesh->gcell_vtx_lst,
                               &vtx_gcells_idx,
                               &vtx_gcells_lst);

  for (face_id = 0; face_id < n_faces; face_id++) {

    /* We compute the cosine of the non-orthogonality angle
       of internal faces (angle between the normal of the face
       and the line between I (center of the cell I) and J (center
       of the cell J) */

    /* Vector IJ and normal of the face */

    cell_i = face_cells[face_id][0];
    cell_j = face_cells[face_id][1];
    dprod = 0;

    for (i = 0; i < 3; i++) {
      v_ij[i] = cell_cen[3*cell_j + i] - cell_cen[3*cell_i + i];
      face_normal[i] = mesh_quantities->i_face_normal[3*face_id + i];
      dprod += v_ij[i]*face_normal[i];
    }

    norm_ij = cs_math_3_norm(v_ij);
    face_norm = cs_math_3_norm(face_normal);

    assert(norm_ij > 0.);
    assert(face_norm > 0.);

    /* Dot product : norm_ij . face_norm */

    cos_ij_fn = dprod / (norm_ij * face_norm);

    /* Comparison to a predefined limit.
       This is non-orthogonal if we are below the limit and so we keep
       the cell in the extended neighborhood of the two cells sharing
       the face. (The cell is tagged (<0) then we will change the sign
       and eliminate all cells < 0 */

    if (cos_ij_fn <= cos_ij_fn_min) {

      /* For each cell sharing the face : intersection between
         cells in the extended neighborhood and cells sharing a
         vertex of the face. */

      _tag_cells(face_id, cell_i, mesh,
                 cell_cells_tag,
                 vtx_cells_idx, vtx_cells_lst,
                 vtx_gcells_idx, vtx_gcells_lst);
      _tag_cells(face_id, cell_j, mesh,
                 cell_cells_tag,
                 vtx_cells_idx, vtx_cells_lst,
                 vtx_gcells_idx, vtx_gcells_lst);

    }

  } /* End of loop on faces */

  /* Free "vertex -> cells" connectivity */

  BFT_FREE(vtx_cells_idx);
  BFT_FREE(vtx_cells_lst);

  if (cs_mesh_n_g_ghost_cells(mesh) > 0) {
    BFT_FREE(vtx_gcells_idx);
    BFT_FREE(vtx_gcells_lst);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Guess cell type based on connectivity.
 *
 * This is not always exact (i.e. a cell with 3 quadrangle and 2 triangle
 * faces might not always be a prism), but is simple and precice enough
 * for heuristics.
 *
 * \param[in]  c_id       cell id
 * \param[in]  i_f2v_idx  interior face->vertices index
 * \param[in]  i_f2v_idx  boundary face->vertices index
 * \param[in]  cell_i_faces_idx  cell->interior faces index
 * \param[in]  cell_i_faces_lst  cell->interior faces adjacency list
 * \param[in]  cell_b_faces_idx  cell->interior faces index
 * \param[in]  cell_b_faces_lst  cell->interior faces adjacency list
 *
 * \return  cell type
 */
/*----------------------------------------------------------------------------*/

static inline fvm_element_t
_guess_cell_type(cs_lnum_t             c_id,
                 const cs_lnum_t       i_f2v_idx[],
                 const cs_lnum_t       b_f2v_idx[],
                 const cs_lnum_t       cell_i_faces_idx[],
                 const cs_lnum_t       cell_i_faces_lst[],
                 const cs_lnum_t       cell_b_faces_idx[],
                 const cs_lnum_t       cell_b_faces_lst[])
{
  fvm_element_t cell_type = FVM_CELL_POLY;

  cs_lnum_t i_f_s_id = cell_i_faces_idx[c_id];
  cs_lnum_t i_f_e_id = cell_i_faces_idx[c_id+1];
  cs_lnum_t b_f_s_id = cell_b_faces_idx[c_id];
  cs_lnum_t b_f_e_id = cell_b_faces_idx[c_id+1];

  cs_lnum_t n_tria = 0;
  cs_lnum_t n_quad = 0;

  for (cs_lnum_t i = i_f_s_id; i < i_f_e_id; i++) {
    cs_lnum_t f_id = cell_i_faces_lst[i];
    cs_lnum_t n_vtx_f = i_f2v_idx[f_id+1] - i_f2v_idx[f_id];
    if (n_vtx_f == 3)
      n_tria++;
    else if (n_vtx_f == 4)
      n_quad++;
  }
  for (cs_lnum_t i = b_f_s_id; i < b_f_e_id; i++) {
    cs_lnum_t f_id = cell_b_faces_lst[i];
    cs_lnum_t n_vtx_f = b_f2v_idx[f_id+1] - b_f2v_idx[f_id];
    if (n_vtx_f == 3)
      n_tria++;
    else if (n_vtx_f == 4)
      n_quad++;
  }

  if (n_tria == 0 && n_quad == 6)
    cell_type = FVM_CELL_HEXA;
  else if (n_tria == 4 && n_quad == 0)
    cell_type = FVM_CELL_TETRA;
  else if (n_tria == 2 && n_quad == 3)
    cell_type = FVM_CELL_PRISM;
  else if (n_tria == 4 && n_quad == 1)
    cell_type = FVM_CELL_PYRAM;

  return cell_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reduce extended neighborhood based on heuristics.
 *
 * \param[in, out]  mesh            pointer to a mesh structure
 * \param[in]       mq              associated mesh quantities
 * \param[in, out]  cell_cells_tag  tag adjacencies to retain
 */
/*----------------------------------------------------------------------------*/

static void
_neighborhood_reduce_optimized(cs_mesh_t             *mesh,
                               cs_mesh_quantities_t  *mq,
                               char                   cell_cells_tag[])
{
  cs_lnum_t  *cell_i_faces_idx = NULL, *cell_i_faces_lst = NULL;
  cs_lnum_t  *cell_b_faces_idx = NULL, *cell_b_faces_lst = NULL;

  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;

  const cs_lnum_t n_cells = mesh->n_cells;

  const cs_lnum_2_t  *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_real_3_t  *cell_cen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_3_t  *b_face_cog = (const cs_real_3_t *)mq->b_face_cog;
  const cs_real_3_t  *i_face_cog = (const cs_real_3_t *)mq->i_face_cog;
  const cs_real_t    *cell_vol = (const cs_real_t *)mq->cell_vol;
  const cs_lnum_t  *i_f2v_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t  *b_f2v_idx = mesh->b_face_vtx_idx;

  const cs_real_t one_over_3 = 1./3.;

  const cs_real_t *b_face_surf = (const cs_real_t *)mq->b_face_surf;
  const cs_real_t *i_face_surf = (const cs_real_t *)mq->i_face_surf;

  /* TODO could use mq->b_sym_flag if it was available */

  /* Min number of cells retained in the ext neighborhood. Setting to 1
     would have same spirit as prev opposite adjacent cell center.
     This parameter could be differentiated for different cell
     types (ex : 10 for tetra, 16 for pyram, 15 for prism and 12 for hexa)
     but "optimal" values are derived from the case of regular
     grids... where they are not used. We only differentiate for hexa
     later to maintain symmetry of the neighborhood.*/

  const cs_lnum_t n_c_min = 10;

  /* Get "cell -> faces" connectivity for the local mesh */

  _get_cell_i_faces_connectivity(mesh,
                                 &cell_i_faces_idx,
                                 &cell_i_faces_lst);

  _get_cell_b_faces_connectivity(mesh,
                                 &cell_b_faces_idx,
                                 &cell_b_faces_lst);

# pragma omp parallel if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t t_s_id, t_e_id;
    _thread_range(n_cells, &t_s_id, &t_e_id);

    cs_real_3_t *n_c_s = NULL;
    cs_lnum_t   *i_c_selected = NULL;
    bool        *is_selected  = NULL;

    cs_lnum_t n_max_c = 0;

    /* Debugging trace */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    FILE  *dbg_trace = NULL;
    {
      char  filename[128];

#if defined(HAVE_OPENMP)
      int t_id = omp_get_thread_num();
      snprintf(filename, 127, "neighborhood_reduce_r%04d_t%02d.txt",
               CS_MAX(cs_glob_rank_id, 0), t_id);
#else
      snprintf(filename, 127, "neighborhood_reduce_r%04d.txt",
               CS_MAX(cs_glob_rank_id, 0));
#endif
      filename[127] = '\0';
      dbg_trace = fopen(filename, "w");
    }
#endif

    /* Loop on cells */

    for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++) {

      /* Build cache for extended neighbors */

      cs_lnum_t c_s_id = cell_cells_idx[c_id];
      cs_lnum_t c_e_id = cell_cells_idx[c_id+1];

      cs_lnum_t n_c = c_e_id - c_s_id; /* n. cells in ext. neighborhood only */

      cs_lnum_t n_c_selected_in_ext = 0;

      cs_lnum_t n_c_min_true = CS_MIN(n_c, n_c_min);

      if (n_c > n_max_c) {
        n_max_c = n_c*2;
        BFT_REALLOC(n_c_s, n_max_c, cs_real_3_t);
        BFT_REALLOC(is_selected, n_max_c, bool);
        BFT_REALLOC(i_c_selected, n_max_c, cs_lnum_t);
      }

      for (cs_lnum_t i = 0; i < n_c; i++) {
        cs_lnum_t c_id_n = cell_cells_lst[c_s_id + i];
        is_selected[i] = false;
        for (cs_lnum_t j = 0; j < 3; j++)
          /* list of vectors IJ after restricted neighborhood */
          n_c_s[i][j] = cell_cen[c_id_n][j] - cell_cen[c_id][j];
      }

      /* Initialize IC on restricted neighborhood (C = neighborhood CoG) */

      cs_real_t ic[3] = {0., 0., 0.};
      cs_lnum_t n_c_restricted = 0.; /* n. cells in restricted neighborhood */

      /* Identify/guess cell type */

      cs_lnum_t i_f_s_id = cell_i_faces_idx[c_id];
      cs_lnum_t i_f_e_id = cell_i_faces_idx[c_id+1];
      cs_lnum_t b_f_s_id = cell_b_faces_idx[c_id];
      cs_lnum_t b_f_e_id = cell_b_faces_idx[c_id+1];

      fvm_element_t cell_type
        = _guess_cell_type(c_id,
                           i_f2v_idx, b_f2v_idx,
                           cell_i_faces_idx, cell_i_faces_lst,
                           cell_b_faces_idx, cell_b_faces_lst);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (dbg_trace != NULL) {
        fprintf(dbg_trace, "\nc_id %d (%lf,%lf,%lf) - ",
                c_id, cell_cen[c_id][0], cell_cen[c_id][1], cell_cen[c_id][2]);

        fprintf(dbg_trace, "%s - ", fvm_element_type_name[cell_type]);
        if (b_f_e_id - b_f_s_id == 1)
          fprintf(dbg_trace, "simple b cell\n");
        else if (b_f_e_id - b_f_s_id > 1)
          fprintf(dbg_trace, "multiple b cell\n");
        else
          fprintf(dbg_trace, "interior cell\n");
      }
#endif

      /* Get regularity of the grid; only hexahedra and prisms are assumed to
         be possibly regular (if regular, the neighborhood will be restricted) */

      bool is_regular = true;
      if (cell_type != FVM_CELL_HEXA && cell_type != FVM_CELL_PRISM)
        is_regular = false;

      cs_real_t avg_dj = 0.;
      cs_real_t sum_s = 0.;

      /* Loop on cell's boundary faces */

      for (cs_lnum_t i = b_f_s_id; i < b_f_e_id; i++) {

        cs_lnum_t f_id_0 = cell_b_faces_lst[i];
        n_c_restricted += 1;
        sum_s+= b_face_surf[f_id_0];

        cs_real_t dj[3];
        for (cs_lnum_t j = 0; j < 3; j++) {
          /* 2 factor is not mandatory below */
          dj[j] = 2*(b_face_cog[f_id_0][j] - cell_cen[c_id][j]);
          ic[j] =   ic[j] * (n_c_restricted-1.)/n_c_restricted
                  + dj[j] * 1./n_c_restricted;
        }

        avg_dj =   avg_dj * (n_c_restricted-1.)/n_c_restricted
                 + cs_math_3_norm(dj) * 1./n_c_restricted;

      } /* End of loop on boundary faces */

      /* Loop on cell's interior faces */

      for (cs_lnum_t i = i_f_s_id; i < i_f_e_id; i++) {

        cs_lnum_t f_id_0 = cell_i_faces_lst[i];
        cs_lnum_t c_id_1 = i_face_cells[f_id_0][1];
        if (c_id_1 == c_id)
          c_id_1 = i_face_cells[f_id_0][0];

        sum_s+= i_face_surf[f_id_0];

        n_c_restricted += 1;

        cs_real_t dj[3];
        for (cs_lnum_t j = 0; j < 3; j++) {
          dj[j] = cell_cen[c_id_1][j] - cell_cen[c_id][j];
          ic[j] =   ic[j] * (n_c_restricted-1.)/n_c_restricted
                  + dj[j] * 1./n_c_restricted;
        }

        avg_dj =   avg_dj * (n_c_restricted-1.)/n_c_restricted
                 + cs_math_3_norm(dj) * 1./n_c_restricted;

        if (is_regular) {
          /* For prisms, only the triangle faces are considered */
          cs_lnum_t n_f_vtx = i_f2v_idx[f_id_0+1] - i_f2v_idx[f_id_0];
          if (cell_type == FVM_CELL_PRISM && n_f_vtx != 3)
            continue;

          cs_real_t v_if[3], v_jf[3];
          for (cs_lnum_t j = 0; j < 3; j++) {
            v_if[j] = i_face_cog[f_id_0][j] - cell_cen[c_id][j];
            v_jf[j] = i_face_cog[f_id_0][j] - cell_cen[c_id_1][j];
          }

          /* IF.FJ (== - IF.JF) */
          cs_real_t dot_p = - cs_math_3_dot_product(v_if, v_jf);

          cs_real_t norm_if_2 = cs_math_3_square_norm(v_if);
          cs_real_t norm_fj_2 = cs_math_3_square_norm(v_jf);

          /* if dot_p < 0.9 * sqrt(norm_if_2) * sqrt(norm_fj_2) */
          if (cs_math_pow2(dot_p) < cs_math_pow2(0.9) * norm_if_2 * norm_fj_2) {
            is_regular = false;
            /* 2 hexa at more than 60 deg */
            if (cell_type == FVM_CELL_HEXA && cs_math_pow2(dot_p) < cs_math_pow2(0.5) * norm_if_2 * norm_fj_2)
              n_c_min_true = 6;
          }
        }

      } /* End of loop on interior faces */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (dbg_trace != NULL) {
        if (is_regular)
          fprintf(dbg_trace, "regular\n");
        else
          fprintf(dbg_trace, "not regular\n");
      }
#endif

      /* At simple boundary, adding neighboors increases the risk of
         spurious modes.
         At double boundary no risk of spurious mode
         (the cell is basically alone). */

      cs_real_t l1 = pow(cell_vol[c_id], one_over_3);

      cs_real_t aspect_r = 6.*l1*l1/sum_s;

      if (b_f_e_id-b_f_s_id == 1)
        is_regular = true;
      else if (b_f_e_id-b_f_s_id > 1) {
        is_regular = false;
        n_c_min_true = n_c_restricted + 1;
      }
      else if (aspect_r < 0.1365) { /* large size / small size > 100
                                       -> use full neighborhood*/
        is_regular = false;
        for (cs_lnum_t i = c_s_id; i < c_e_id; i++)
          cell_cells_tag[i] = 1; /* keep this cell */
      }

      /* Limit to a restricted neighborhood if the grid is regular */

      if (! is_regular) {
        cs_lnum_t n_c_selected = n_c_restricted;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (dbg_trace != NULL)
          fprintf(dbg_trace, "n_c_min %d\n", n_c_min_true);
#endif

        cs_real_t min_crit = HUGE_VAL;

        for (cs_lnum_t k = 0; k < n_c; k++)
          i_c_selected[k] = -1;

        for (cs_lnum_t k = 0; k < n_c; k++) {

          cs_real_t min_loc_crit = HUGE_VAL;
          for (cs_lnum_t i = 0; i < n_c; i++) {
            if (is_selected[i] == false) {
              cs_real_t loc_crit = 0.;
              for (cs_lnum_t j = 0; j < 3; j++)
                loc_crit+= cs_math_pow2(n_c_selected * ic[j] + 2*n_c_s[i][j]);
              if (loc_crit < min_loc_crit) {
                i_c_selected[k] = i;
                min_loc_crit = loc_crit;
              }
            }
          }
          is_selected[i_c_selected[k]] = true;
          n_c_selected += 1;

          for (cs_lnum_t j = 0; j < 3; j++)
            ic[j] =   ic[j] * (n_c_selected-1.)/n_c_selected
                    + n_c_s[i_c_selected[k]][j] * 1./n_c_selected;

          cs_real_t norm_ic = cs_math_3_norm(ic);

          avg_dj =   avg_dj * (n_c_selected-1.)/n_c_selected
                   + cs_math_3_norm(n_c_s[i_c_selected[k]]) * 1./n_c_selected;

          cs_real_t crit = (norm_ic + 0.1*l1) * avg_dj;
          /* Optimize offset and, to a small extent, the average distance
             of the selected elements. */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (dbg_trace != NULL)
            fprintf(dbg_trace, "crit %d %lf (norm_ic %lf - avg_dj %lf)\n",
                    k, crit, norm_ic, avg_dj);
#endif

          if (crit < min_crit && k >= n_c_min_true-1) {
            min_crit = crit;
            n_c_selected_in_ext = k+1;
            if (norm_ic < 0.01*l1)
              break;
          }

        } /* End of crit loop (definition of the neighborhood) */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (dbg_trace != NULL)
          fprintf(dbg_trace,
                  "n_c_selected_in_ext %ld\n",
                  (long)n_c_selected_in_ext);
#endif

        for (cs_lnum_t k = 0; k < n_c_selected_in_ext; k++)
          cell_cells_tag[c_s_id + i_c_selected[k]] = 1; /* keep this cell */

      } /* End if not regular */

    } /* End of loop on cells */

    BFT_FREE(n_c_s);
    BFT_FREE(i_c_selected);
    BFT_FREE(is_selected);

  } /* End of OpenMP block */

  BFT_FREE(cell_i_faces_idx);
  BFT_FREE(cell_i_faces_lst);
  BFT_FREE(cell_b_faces_idx);
  BFT_FREE(cell_b_faces_lst);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reduce extended neighborhood based on cell center opposite.
 *
 * \param[in, out]  mesh            pointer to a mesh structure
 * \param[in]       mq              associated mesh quantities
 * \param[in, out]  cell_cells_tag  tag adjacencies to retain
 */
/*----------------------------------------------------------------------------*/

static void
_neighborhood_reduce_cell_center_opposite(cs_mesh_t             *mesh,
                                          cs_mesh_quantities_t  *mq,
                                          char                   cell_cells_tag[])
{
  cs_lnum_t  *cell_i_faces_idx = NULL, *cell_i_faces_lst = NULL;
  cs_lnum_t  *cell_b_faces_idx = NULL, *cell_b_faces_lst = NULL;

  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;

  const cs_lnum_t n_cells = mesh->n_cells;

  const cs_lnum_2_t  *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_real_3_t  *cell_cen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_3_t  *b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

  /* Get "cell -> faces" connectivity for the local mesh */

  _get_cell_i_faces_connectivity(mesh,
                                 &cell_i_faces_idx,
                                 &cell_i_faces_lst);

  _get_cell_b_faces_connectivity(mesh,
                                 &cell_b_faces_idx,
                                 &cell_b_faces_lst);

# pragma omp parallel if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n_cells, sizeof(cs_real_t), &t_s_id, &t_e_id);

    cs_real_3_t  *n_c_s = NULL;

    cs_lnum_t n_max_c = 0;

    /* Loop on cells */

    for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++) {

      /* Build cache for extended neighbors */

      cs_lnum_t c_s_id = cell_cells_idx[c_id];
      cs_lnum_t c_e_id = cell_cells_idx[c_id+1];
      cs_lnum_t n_c = c_e_id - c_s_id;

      if (n_c > n_max_c) {
        n_max_c = n_c*2;
        BFT_REALLOC(n_c_s, n_max_c, cs_real_3_t);
      }

      for (cs_lnum_t i = 0; i < n_c; i++) {
        cs_lnum_t c_id_n = cell_cells_lst[c_s_id + i];
        for (cs_lnum_t j = 0; j < 3; j++)
          n_c_s[i][j] = cell_cen[c_id_n][j] - cell_cen[c_id][j];
      }

      /* Loop on interior cell faces */

      cs_lnum_t i_f_s_id = cell_i_faces_idx[c_id];
      cs_lnum_t i_f_e_id = cell_i_faces_idx[c_id+1];

      cs_lnum_t b_f_s_id = cell_b_faces_idx[c_id];
      cs_lnum_t b_f_e_id = cell_b_faces_idx[c_id+1];

      for (cs_lnum_t i = i_f_s_id; i < i_f_e_id; i++) {

        cs_real_t i_dist_align_min = HUGE_VAL;

        cs_lnum_t f_id_0 = cell_i_faces_lst[i];
        cs_lnum_t c_id_0 = i_face_cells[f_id_0][0];
        if (c_id_0 == c_id)
          c_id_0 = i_face_cells[f_id_0][1];

        cs_real_t seg_0[3];
        for (cs_lnum_t k = 0; k < 3; k++)
          seg_0[k] = cell_cen[c_id][k] - cell_cen[c_id_0][k];

        /* Initialize using base (face->cells) connectivity) */

        for (cs_lnum_t j = i_f_s_id; j < i_f_e_id; j++) {

          if (i == j) continue;

          cs_real_t seg_1[3];

          cs_lnum_t f_id_1 = cell_i_faces_lst[j];
          cs_lnum_t c_id_1 = i_face_cells[f_id_1][0];
          if (c_id_1 == c_id)
            c_id_1 = i_face_cells[f_id_1][1];
          for (cs_lnum_t k = 0; k < 3; k++)
            seg_1[k] = cell_cen[c_id_1][k] - cell_cen[c_id][k];

          if (cs_math_3_dot_product(seg_0, seg_1) <= 0)
            continue;

          cs_real_t d2 = _cross_product_sq_norm(seg_0, seg_1);

          if (d2 < i_dist_align_min)
            i_dist_align_min = d2;

        }

        for (cs_lnum_t j = b_f_s_id; j < b_f_e_id; j++) {

          cs_real_t seg_1[3];

          cs_lnum_t f_id_1 = cell_b_faces_lst[j];
          for (cs_lnum_t k = 0; k < 3; k++)
            seg_1[k] = b_face_cog[f_id_1][k] - cell_cen[c_id][k];

          if (cs_math_3_dot_product(seg_0, seg_1) <= 0)
            continue;

          cs_real_t d2 = _cross_product_sq_norm(seg_0, seg_1);

          if (d2 < i_dist_align_min)
            i_dist_align_min = d2;

        }

        /* Compare with cells in extended neighborhood */

        cs_lnum_t i_cell_align_min = -1;
        for (cs_lnum_t j = 0; j < n_c; j++) {

          if (cs_math_3_dot_product(seg_0, n_c_s[j]) <= 0)
            continue;

          cs_real_t d2 = _cross_product_sq_norm(seg_0, n_c_s[j]);

          if (d2 < i_dist_align_min) {
            i_cell_align_min = j;
            i_dist_align_min = d2;
          }

        }

        if (i_cell_align_min > -1)
          cell_cells_tag[c_s_id + i_cell_align_min] = 1;

      } /* End of loop on interior cell faces */

      /* Loop on boundary cell faces */

      for (cs_lnum_t i = b_f_s_id; i < b_f_e_id; i++) {

        cs_real_t i_dist_align_min = HUGE_VAL;

        cs_lnum_t f_id_0 = cell_b_faces_lst[i];

        cs_real_t seg_0[3];
        for (cs_lnum_t k = 0; k < 3; k++)
          seg_0[k] = cell_cen[c_id][k] - b_face_cog[f_id_0][k];

        /* Initialize using base (face->cells) connectivity) */

        for (cs_lnum_t j = i_f_s_id; j < i_f_e_id; j++) {

          cs_real_t seg_1[3];

          cs_lnum_t f_id_1 = cell_i_faces_lst[j];
          cs_lnum_t c_id_1 = i_face_cells[f_id_1][0];
          if (c_id_1 == c_id)
            c_id_1 = i_face_cells[f_id_1][1];
          for (cs_lnum_t k = 0; k < 3; k++)
            seg_1[k] = cell_cen[c_id_1][k] - cell_cen[c_id][k];

          if (cs_math_3_dot_product(seg_0, seg_1) <= 0)
            continue;

          cs_real_t d2 = _cross_product_sq_norm(seg_0, seg_1);

          if (d2 < i_dist_align_min)
            i_dist_align_min = d2;

        }

        for (cs_lnum_t j = b_f_s_id; j < b_f_e_id; j++) {

          if (i == j) continue;

          cs_real_t seg_1[3];

          cs_lnum_t f_id_1 = cell_b_faces_lst[j];
          for (cs_lnum_t k = 0; k < 3; k++)
            seg_1[k] = b_face_cog[f_id_1][k] - cell_cen[c_id][k];

          if (cs_math_3_dot_product(seg_0, seg_1) <= 0)
            continue;

          cs_real_t d2 = _cross_product_sq_norm(seg_0, seg_1);

          if (d2 < i_dist_align_min)
            i_dist_align_min = d2;

        }

        /* Compare with cells in extended neighborhood */

        cs_lnum_t i_cell_align_min = -1;
        for (cs_lnum_t j = 0; j < n_c; j++) {

          if (cs_math_3_dot_product(seg_0, n_c_s[j]) <= 0)
            continue;

          cs_real_t d2 = _cross_product_sq_norm(seg_0, n_c_s[j]);

          if (d2 < i_dist_align_min) {
            i_cell_align_min = j;
            i_dist_align_min = d2;
          }

        }

        if (i_cell_align_min > -1)
          cell_cells_tag[c_s_id + i_cell_align_min] = 1;

      } /* End of loop on boundary cell faces */

    } /* End of loop on cells */

    BFT_FREE(n_c_s);

  } /* End of Open MP block */

  BFT_FREE(cell_i_faces_idx);
  BFT_FREE(cell_i_faces_lst);
  BFT_FREE(cell_b_faces_idx);
  BFT_FREE(cell_b_faces_lst);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mark cells adjacent to boundary for extended neighborhood selection.
 *
 * \param[in, out]  mesh            pointer to a mesh structure
 * \param[in, out]  cell_cells_tag  tag adjacencies to retain
 */
/*----------------------------------------------------------------------------*/

static void
_neighborhood_reduce_full_boundary(cs_mesh_t             *mesh,
                                   char                   cell_cells_tag[])
{
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;

  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_cells = mesh->n_b_cells;

  int *cell_b_flag;
  BFT_MALLOC(cell_b_flag, n_cells_ext, int);

  cs_mesh_tag_boundary_cells(mesh, cell_b_flag);

# pragma omp parallel if (n_b_cells > CS_THR_MIN)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n_b_cells, sizeof(cs_real_t), &t_s_id, &t_e_id);

    /* Loop on boundary cells */

    for (cs_lnum_t b_c_id = t_s_id; b_c_id < t_e_id; b_c_id++) {

      cs_lnum_t c_id = mesh->b_cells[b_c_id];

      /* Build cache for extended neighbors */

      cs_lnum_t c_s_id = cell_cells_idx[c_id];
      cs_lnum_t c_e_id = cell_cells_idx[c_id+1];

      for (cs_lnum_t i = c_s_id; i < c_e_id; i++) {
        if (cell_b_flag[cell_cells_lst[i]])
          cell_cells_tag[i] = 1;
      }

    } /* End of loop on boundary cells */

  } /* End of Open MP block */

  BFT_FREE(cell_b_flag);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the extended neighborhood type.
 *
 * \return  extended neighborhood type
 */
/*----------------------------------------------------------------------------*/

cs_ext_neighborhood_type_t
cs_ext_neighborhood_get_type(void)
{
  return _ext_nbh_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the extended neighborhood type.
 *
 * \param[in]  enh_type  extended neighborhood type
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_set_type(cs_ext_neighborhood_type_t  enh_type)
{
  if (   enh_type < CS_EXT_NEIGHBORHOOD_NONE
      || enh_type > CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX)
    enh_type = CS_EXT_NEIGHBORHOOD_NONE;

  _ext_nbh_type = enh_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query whether full extended neighborhood should be kept at boundaries.
 *
 * \return  extended neighborhood type
 */
/*----------------------------------------------------------------------------*/

bool
cs_ext_neighborhood_get_boundary_complete(void)
{
  return _full_nb_boundary;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the extended neighborhood type.
 *
 * \param[in]  keep  true if full extended neighborhood should be maintained
 *             at boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_set_boundary_complete(bool  keep)
{
  _full_nb_boundary = keep;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the non_orthogonality threshold  (in degrees) associated with the
 *        CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX neighborhood type.
 *
 * \return  non-orthogonality threshold
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ext_neighborhood_get_non_ortho_max(void)
{
  return _non_ortho_max;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the non_orthogonality threshold  (in degrees) associated with the
 *        CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX neighborhood type.
 *
 * \param[in]  non_ortho_max  non-orthogonality threshold
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_set_non_ortho_max(cs_real_t  non_ortho_max)
{
  _non_ortho_max = non_ortho_max;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reduce the "cell -> cells" connectivity for the
 *        extended neighborhood using a non-orthogonality criterion.
 *
 * Note: Only cells sharing only a vertex or vertices (not a face)
 *       belong to the "cell -> cells" connectivity.
 *
 * \param[in]  mesh             pointer to mesh structure
 * \param[in]  mesh_quantities  associated mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_reduce(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities)
{
  /* Return if there is no extended neighborhood */

  if (   mesh->cell_cells_idx == NULL
      || (mesh->halo_type == CS_HALO_STANDARD && _full_nb_boundary == false)
      || _ext_nbh_type == CS_EXT_NEIGHBORHOOD_COMPLETE)
    return;

  cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;
  cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;

  const cs_lnum_t  n_cells = mesh->n_cells;

  if (mesh->verbosity > 0) {
    char full_b_str[64] = "";
    if (_full_nb_boundary) {
      strncpy(full_b_str, _(", complete on boundary"), 63);
      full_b_str[63] = '\0';
    }
    bft_printf
      (_("\n"
         " Extended neighborhood reduction: %s%s\n"
         " --------------------------------\n\n"),
       _(cs_ext_neighborhood_type_name[_ext_nbh_type]),
       full_b_str);
  }

  /* Tag cells to eliminate (0 to eliminate, 1 to keep) */

  char *cell_cells_tag;
  BFT_MALLOC(cell_cells_tag, mesh->cell_cells_idx[n_cells], char);
  for (cs_lnum_t i = 0; i < mesh->cell_cells_idx[n_cells]; i++)
    cell_cells_tag[i] = 0;

  switch(_ext_nbh_type) {
  case CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX:
    _neighborhood_reduce_anomax(mesh,
                                mesh_quantities,
                                cell_cells_tag);
    break;
  case CS_EXT_NEIGHBORHOOD_OPTIMIZED:
    _neighborhood_reduce_optimized(mesh,
                                   mesh_quantities,
                                   cell_cells_tag);
    break;
  case CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE:
    _neighborhood_reduce_cell_center_opposite(mesh,
                                              mesh_quantities,
                                              cell_cells_tag);
    break;
  default:
    break;
  }

  if (_full_nb_boundary)
    _neighborhood_reduce_full_boundary(mesh, cell_cells_tag);

  /* Delete cells with tag == 0 */

  cs_gnum_t n_n_sum[2] = {mesh->cell_cells_idx[n_cells], 0};
  cs_lnum_t n_n_min[2] = {mesh->cell_cells_idx[n_cells],
                          mesh->cell_cells_idx[n_cells]};
  cs_lnum_t n_n_max[2] = {0, 0};

  cs_lnum_t s_id = 0, j = 0;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    cs_lnum_t e_id = cell_cells_idx[cell_id+1];

    cs_lnum_t n_c_ini = e_id - s_id;
    cs_lnum_t n_c_red = 0;
    n_n_min[0] = CS_MIN(n_n_min[0], n_c_ini);
    n_n_max[0] = CS_MAX(n_n_max[0], n_c_ini);

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      if (cell_cells_tag[i] != 0) {
        cell_cells_lst[j] = cell_cells_lst[i];
        j++;
        n_c_red++;
      }

    } /* End of loop on cells in the extended neighborhood of cell_id+1 */

    s_id = e_id;
    cell_cells_idx[cell_id+1] = j;

    n_n_min[1] = CS_MIN(n_n_min[1], n_c_red);
    n_n_max[1] = CS_MAX(n_n_max[1], n_c_red);

  } /* End of loop on cells */

  n_n_sum[1] = mesh->cell_cells_idx[n_cells];

  BFT_FREE(cell_cells_tag);

  /* Reallocation of cell_cells_lst */

  BFT_REALLOC(mesh->cell_cells_lst, cell_cells_idx[n_cells], cs_lnum_t);

  /* Output for log */

  if (mesh->verbosity > 0) {

    cs_parall_sum(2, CS_GNUM_TYPE, n_n_sum);
    cs_parall_min(2, CS_LNUM_TYPE, n_n_min);
    cs_parall_max(2, CS_LNUM_TYPE, n_n_max);

    double n_n_mean[2];
    for (int i = 0; i < 2; i++)
      n_n_mean[i] = (double)(n_n_sum[i]) / (double)(mesh->n_g_cells);

    bft_printf
      (_(" Complete extended cell-cell connectivity size: %llu\n"
         "   Mean neighbors: %3.2f\n"
         "   Min. neighbors: %d\n"
         "   Max. neighbors: %d\n"
         "\n"
         " Reduced extended cell-cell connectivity size: %llu\n"
         "   Mean neighbors: %3.2f\n"
         "   Min. neighbors: %d\n"
         "   Max. neighbors: %d\n"),
       (unsigned long long)n_n_sum[0], n_n_mean[0],
       (int)n_n_min[0], (int)n_n_max[0],
       (unsigned long long)n_n_sum[1], n_n_mean[1],
       (int)n_n_min[1], (int)n_n_max[1]);

  }

#if 0 /* For debugging purposes */
  for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
    bft_printf(" cell %d :: ", i+1);
    for (cs_lnum_t j = mesh->cell_cells_idx[i];
         j < mesh->cell_cells_idx[i+1];
         j++)
      bft_printf(" %d ", mesh->cell_cells_lst[j]);
    bft_printf("\n");
  }
#endif

  cs_sort_indexed(n_cells, mesh->cell_cells_idx, mesh->cell_cells_lst);

  cs_mesh_adjacencies_update_cell_cells_e();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the  "cell -> cells" connectivity.
 *
 * \param[in, out]  mesh  pointer to a mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_define(cs_mesh_t  *mesh)
{
  cs_lnum_t  *vtx_gcells_idx = NULL, *vtx_gcells_lst = NULL;
  cs_lnum_t  *vtx_cells_idx = NULL, *vtx_cells_lst = NULL;
  cs_lnum_t  *cell_i_faces_idx = NULL, *cell_i_faces_lst = NULL;
  cs_lnum_t  *cell_cells_idx = NULL, *cell_cells_lst = NULL;

  cs_halo_t  *halo = mesh->halo;

  /* Get "cell -> faces" connectivity for the local mesh */

  _get_cell_i_faces_connectivity(mesh,
                                 &cell_i_faces_idx,
                                 &cell_i_faces_lst);

  /* Create a "vertex -> cell" connectivity */

  _create_vtx_cells_connect2(mesh,
                             cell_i_faces_idx,
                             cell_i_faces_lst,
                             &vtx_cells_idx,
                             &vtx_cells_lst);

  if (cs_mesh_n_g_ghost_cells(mesh) > 0) {

    /* Create a "vertex -> ghost cells" connectivity */

    _create_vtx_gcells_connect(halo,
                               mesh->n_vertices,
                               mesh->gcell_vtx_idx,
                               mesh->gcell_vtx_lst,
                               &vtx_gcells_idx,
                               &vtx_gcells_lst);

  }

  /* Create the "cell -> cells" connectivity for the extended halo */

  _create_cell_cells_connect(mesh,
                             cell_i_faces_idx,
                             cell_i_faces_lst,
                             vtx_gcells_idx,
                             vtx_gcells_lst,
                             vtx_cells_idx,
                             vtx_cells_lst,
                             &cell_cells_idx,
                             &cell_cells_lst);

  mesh->cell_cells_idx = cell_cells_idx;
  mesh->cell_cells_lst = cell_cells_lst;

  /* Free memory */

  BFT_FREE(vtx_gcells_idx);
  BFT_FREE(vtx_gcells_lst);

  BFT_FREE(cell_i_faces_idx);
  BFT_FREE(cell_i_faces_lst);
  BFT_FREE(vtx_cells_idx);
  BFT_FREE(vtx_cells_lst);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
