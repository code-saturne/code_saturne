/*============================================================================
 * Initialization of a nodal connectivity definition based upon
 * a (possibly partial) descending connectivity
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Local headers associated with the current file
 *----------------------------------------------------------------------------*/

#include "fvm_nodal_from_desc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Enumeration definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Nodal connectivity reconstruction return code.
 *----------------------------------------------------------------------------*/

typedef enum {

  FVM_NODAL_FROM_DESC_SUCCESS,  /* Successful reconstruction */
  FVM_NODAL_FROM_DESC_FAILURE,  /* Reconstruction failure (connectivity pb.) */
  FVM_NODAL_FROM_DESC_WORIENT   /* Reconstruction with orientation warning
                                   (for possible orientation problem) */
} fvm_nodal_from_desc_t;

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Global static variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine if a given cell is a prism or a polyhedron
 *
 * This function should only be called on cells having 2 triangular and
 * 3 quadrangular faces; if the triangles are adjacent, the cell
 * is not a prism, and is considered to be a polyhedron. Otherwise, it
 * seems to be truly a prism (supposing it is closed and well oriented).
 *
 * parameters:
 *   cell_id,        <-- cell id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *
 * returns:
 *   type of cell defined by cell_id
 *----------------------------------------------------------------------------*/

inline static fvm_element_t
_is_prism_or_poly(const cs_lnum_t    cell_id,
                  const int          n_face_lists,
                  const cs_lnum_t    face_list_shift[],
                  const cs_lnum_t   *face_vertex_idx[],
                  const cs_lnum_t   *face_vertex[],
                  const cs_lnum_t    cell_face_idx[],
                  const cs_lnum_t    cell_face_num[])
{
  int         vtx_id_1, vtx_id_2;
  cs_lnum_t   face_id, fl;
  cs_lnum_t   idx, idx_start, idx_end;
  cs_lnum_t   n_face_vertices;
  cs_lnum_t   vtx, vertex_id_start, vertex_id_end;

  cs_lnum_t   vtx_tria[6];
  int         n_trias = 0;

  /* Extract 2 triangles */

  idx_start = cell_face_idx[cell_id]     - 1;
  idx_end   = cell_face_idx[cell_id + 1] - 1;

  for (idx = idx_start ; idx < idx_end ; idx++) {

    face_id = CS_ABS(cell_face_num[idx]) - 1;

    for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
    assert(fl > -1);
    face_id -= face_list_shift[fl];

    vertex_id_start = face_vertex_idx[fl][face_id];
    vertex_id_end   = face_vertex_idx[fl][face_id + 1];
    n_face_vertices = vertex_id_end - vertex_id_start;

    if (n_face_vertices == 3) {

      if (cell_face_num[idx] > 0) {
        for (vtx = 0 ; vtx < 3 ; vtx++)
          vtx_tria[n_trias*3 + vtx]
            = face_vertex[fl][vertex_id_start + vtx] + 1;
      }
      else {
        for (vtx = 0 ; vtx < 3 ; vtx++)
          vtx_tria[n_trias*3 + vtx]
            = face_vertex[fl][vertex_id_end - 1 - vtx] + 1;
      }

      n_trias += 1;

      if (n_trias == 2) /* Once we have found the two triangles, */
        break;          /*  we do not need to look at other faces */
    }

  }

  assert(n_trias == 2);

  /* Are the triangles adjacent ? */

  for (vtx_id_1 = 0; vtx_id_1 < 3; vtx_id_1++) {
    for (vtx_id_2 = 0; vtx_id_2 < 3; vtx_id_2++) {
      if (vtx_tria[3 + vtx_id_2] == vtx_tria[vtx_id_1]) {
        return FVM_CELL_POLY;
      }
    }
  }

  /* If no adjacent triangles were found, we have a prism */

  return FVM_CELL_PRISM;
}

/*----------------------------------------------------------------------------
 * Determination of a given cell's type.
 *
 * If the optional cell_vtx_tria[3*4] and cell_vtx_quad[4*6] arrays are given,
 * they are filled with the vertex indexes of the cell's triangle and
 * quadrangle type faces (as long as no simple polyhedral face is encountered
 * and we have no more than 4 triangles or 6 quadrangles), for easier
 * nodal cell reconstuction. Vertex indexes are ordered as usual for
 * an outward pointing normal.
 *
 * parameters:
 *   cell_id,        <-- cell id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   cell_vtx_tria   --> local triangle definitions (optional, 4 max.)
 *   cell_vtx_quad   --> local quadrangle definitions (optional, 6 max.)
 *
 * returns:
 *   type of cell defined by cell_id
 *----------------------------------------------------------------------------*/

inline static fvm_element_t
_nodal_cell_from_desc(const cs_lnum_t    cell_id,
                      const int          n_face_lists,
                      const cs_lnum_t    face_list_shift[],
                      const cs_lnum_t   *face_vertex_idx[],
                      const cs_lnum_t   *face_vertex[],
                      const cs_lnum_t    cell_face_idx[],
                      const cs_lnum_t    cell_face_num[],
                      cs_lnum_t         *cell_vtx_tria,
                      cs_lnum_t         *cell_vtx_quad)
{
  cs_lnum_t   n_cell_faces;
  cs_lnum_t   face_id, fl;
  cs_lnum_t   idx, idx_start, idx_end;
  cs_lnum_t   n_face_vertices;
  cs_lnum_t   vtx, vertex_id_start, vertex_id_end;

  fvm_element_t cell_type;

  cs_lnum_t   n_trias = 0;
  cs_lnum_t   n_quads = 0;
  cs_lnum_t   n_ngons = 0;


  /* Guess connectivity types */
  /*--------------------------*/

  n_cell_faces = cell_face_idx[cell_id + 1] - cell_face_idx[cell_id];

  /* If we have more than 6 faces, we have a general simple polyhedron */

  if (n_cell_faces > 6)
    return FVM_CELL_POLY;

  /*
    Otherwise, we probably have a "classical" element; for example,
    with 6 faces, we probably have a hexahedron, though we could
    also have a tetrahedron with one of its edges and thus two of
    its faces split. We count vertices per face to check.
  */

  idx_start = cell_face_idx[cell_id]     - 1;
  idx_end   = cell_face_idx[cell_id + 1] - 1;

  for (idx = idx_start ; idx < idx_end ; idx++) {

    face_id = CS_ABS(cell_face_num[idx]) - 1;

    for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
    assert(fl > -1);
    face_id -= face_list_shift[fl];

    vertex_id_start = face_vertex_idx[fl][face_id];
    vertex_id_end   = face_vertex_idx[fl][face_id + 1];
    n_face_vertices = vertex_id_end - vertex_id_start;

    if (n_face_vertices == 3) {

      if (cell_vtx_tria != NULL && n_trias < 4) {
        if (cell_face_num[idx] > 0) {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_tria[n_trias*3 + vtx]
              = face_vertex[fl][vertex_id_start + vtx] + 1;
        }
        else {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_tria[n_trias*3 + vtx]
              = face_vertex[fl][vertex_id_end - 1 - vtx] + 1;
        }
      }

      n_trias += 1;

    }
    else if (n_face_vertices == 4) {

      if (cell_vtx_quad != NULL && n_quads < 6) {
        if (cell_face_num[idx] > 0) {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_quad[n_quads*4 + vtx]
              = face_vertex[fl][vertex_id_start + vtx] + 1;
        }
        else {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_quad[n_quads*4 + vtx]
              = face_vertex[fl][vertex_id_end - 1 - vtx] + 1;
        }
      }

      n_quads += 1;
    }
    else

      n_ngons += 1;

  }

  /* Return element type */

  if (n_ngons > 0)
    cell_type = FVM_CELL_POLY;
  else {
    if (n_trias == 0 && n_quads == 6)
      cell_type = FVM_CELL_HEXA;
    else if (n_trias == 2 && n_quads == 3)
      cell_type = _is_prism_or_poly(cell_id,
                                    n_face_lists,
                                    face_list_shift,
                                    face_vertex_idx,
                                    face_vertex,
                                    cell_face_idx,
                                    cell_face_num);
    else if (n_trias == 4) {
      if (n_quads == 0)
        cell_type = FVM_CELL_TETRA;
      else if (n_quads == 1)
        cell_type = FVM_CELL_PYRAM;
      else
        cell_type = FVM_CELL_POLY;
    }
    else
      cell_type = FVM_CELL_POLY;
  }

  return cell_type;

}

/*----------------------------------------------------------------------------
 * Tetrahedron construction
 *
 * parameters:
 *   cell_vtx_tria  <-- triangular faces connectivity
 *   cell_vtx_tetra --> tetrahedron connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvm_nodal_from_desc_t
_nodal_from_desc_cnv_cel_tetra(const cs_lnum_t   cell_vtx_tria[],
                               cs_lnum_t   cell_vtx_tetra[])
{
  cs_lnum_t   vertex_id, face_id;
  cs_lnum_t   direction;
  cs_lnum_t   vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("face 1 : %d %d %d\n",
             cell_vtx_tria[0], cell_vtx_tria[1], cell_vtx_tria[2]);
  bft_printf("face 2 : %d %d %d\n",
             cell_vtx_tria[3], cell_vtx_tria[4], cell_vtx_tria[5]);
  bft_printf("face 3 : %d %d %d\n",
             cell_vtx_tria[6], cell_vtx_tria[7], cell_vtx_tria[8]);
  bft_printf("face 4 : %d %d %d\n",
             cell_vtx_tria[9], cell_vtx_tria[10], cell_vtx_tria[11]);
#endif

  /*
    Base : vertices of the first face; we take the vertices in opposite
    order ("bottom" face numbering of the tetrahedron with outward
    pointing normal).

    *        x 4
    *       /|\
    *      / | \
    *     /  |  \
    *  1 x- -|- -x 3
    *     \  |  /
    *      \ | /
    *       \|/
    *        x 2
  */


  cell_vtx_tetra[0] = cell_vtx_tria[2];
  cell_vtx_tetra[1] = cell_vtx_tria[1];
  cell_vtx_tetra[2] = cell_vtx_tria[0];

  /*
    We have found 3 of 4 vertices; all other triangles should share
    vertex 4, and one of those should share vertices 1 and 2 of the
    base triangle.
  */

  vtx_num_1 = cell_vtx_tetra[0];
  vtx_num_2 = cell_vtx_tetra[1];

  direction = 0;

  for (face_id = 1 ; face_id < 4 ; face_id++) {

    for (vertex_id = 0 ; vertex_id < 3 ; vertex_id++) {

      vtx_num = cell_vtx_tria[face_id*3 + vertex_id];

      if (vtx_num == vtx_num_1) {
        if (cell_vtx_tria[face_id*3 + ((vertex_id+1) % 3)] == vtx_num_2) {
          direction = 1;
          break;
        }
        else if (cell_vtx_tria[face_id*3 + ((vertex_id-1+3) % 3)] == vtx_num_2) {
          direction = -1;
          break;
        }
      }

    }

    if (direction != 0)
      break;

  }

  if (direction == -1)
    warn_orient = true;
  else if (direction == 0)
    return FVM_NODAL_FROM_DESC_FAILURE;

  cell_vtx_tetra[3]
    = cell_vtx_tria[face_id*3 + ((vertex_id + (3 + (2 * direction))) % 3)];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("tetra : %d %d %d %d\n",
             cell_vtx_tetra[0], cell_vtx_tetra[1],
             cell_vtx_tetra[2], cell_vtx_tetra[3]);
#endif

  if (warn_orient == true)
    return FVM_NODAL_FROM_DESC_WORIENT;
  else
    return FVM_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Pyramid construction
 *
 * parameters:
 *   cell_vtx_tria  <-- triangular faces connectivity
 *   cell_vtx_quad  <-- quadrangle faces connectivity
 *   cell_vtx_pyram --> pyramid connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvm_nodal_from_desc_t
_nodal_from_desc_cnv_cel_pyram(const cs_lnum_t   cell_vtx_tria[],
                               const cs_lnum_t   cell_vtx_quad[],
                               cs_lnum_t   cell_vtx_pyram[])
{
  cs_lnum_t   vertex_id, face_id;
  cs_lnum_t   direction;
  cs_lnum_t   vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("face 1 : %d %d %d %d\n",
             cell_vtx_quad[0], cell_vtx_quad[1],
             cell_vtx_quad[2], cell_vtx_quad[3]);
  bft_printf("face 2 : %d %d %d\n",
             cell_vtx_tria[0], cell_vtx_tria[1], cell_vtx_tria[2]);
  bft_printf("face 3 : %d %d %d\n",
             cell_vtx_tria[3], cell_vtx_tria[4], cell_vtx_tria[5]);
  bft_printf("face 4 : %d %d %d\n",
             cell_vtx_tria[6], cell_vtx_tria[7], cell_vtx_tria[8]);
  bft_printf("face 5 : %d %d %d\n",
             cell_vtx_tria[9], cell_vtx_tria[10], cell_vtx_tria[11]);
#endif

  /*
    Base : vertices of the quadrangle; we take the vertices in opposite
    order ("bottom" face numbering of the pyramid with outward
    pointing normal).

    *         5 x
    *          /|\
    *         //| \
    *        // |  \
    *     4 x/--|---x 3
    *      //   |  /
    *     //    | /
    *  1 x-------x 2
  */


  cell_vtx_pyram[0] = cell_vtx_quad[3];
  cell_vtx_pyram[1] = cell_vtx_quad[2];
  cell_vtx_pyram[2] = cell_vtx_quad[1];
  cell_vtx_pyram[3] = cell_vtx_quad[0];

  /*
    We have found 4 out of 5 vertices; all 4 triangles should share
    vertex 5, and one of those should share vertices 1 and 2 of the
    base quadrangle.
  */

  vtx_num_1 = cell_vtx_pyram[0];
  vtx_num_2 = cell_vtx_pyram[1];

  direction = 0;

  for (face_id = 0 ; face_id < 4 ; face_id++) {

    for (vertex_id = 0 ; vertex_id < 3 ; vertex_id++) {

      vtx_num = cell_vtx_tria[face_id*3 + vertex_id];

      if (vtx_num == vtx_num_1) {
        if (cell_vtx_tria[face_id*3 + ((vertex_id+1) % 3)] == vtx_num_2) {
          direction = 1;
          break;
        }
        else if (cell_vtx_tria[face_id*3 + ((vertex_id-1+3) % 3)] == vtx_num_2) {
          direction = -1;
          break;
        }
      }

    }

    if (direction != 0)
      break;

  }

  if (direction == -1)
    warn_orient = true;
  else if (direction == 0)
    return FVM_NODAL_FROM_DESC_FAILURE;

  cell_vtx_pyram[4]
    = cell_vtx_tria[face_id*3 + ((vertex_id + (3 + (2 * direction))) % 3)];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("pyram : %d %d %d %d %d\n",
             cell_vtx_pyram[0], cell_vtx_pyram[1], cell_vtx_pyram[2],
             cell_vtx_pyram[3], cell_vtx_pyram[4]);
#endif

  if (warn_orient == true)
    return FVM_NODAL_FROM_DESC_WORIENT;
  else
    return FVM_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Prism (pentahedron) construction
 *
 * parameters:
 *   cell_vtx_tria  <-- triangular faces connectivity
 *   cell_vtx_quad  <-- quadrangle faces connectivity
 *   cell_vtx_prism --> prism connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvm_nodal_from_desc_t
_nodal_from_desc_cnv_cel_prism(const cs_lnum_t   cell_vtx_tria[],
                               const cs_lnum_t   cell_vtx_quad[],
                               cs_lnum_t   cell_vtx_prism[])
{
  cs_lnum_t   vertex_id, face_id;
  cs_lnum_t   ipass, direction;
  cs_lnum_t   vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("face 1 : %d %d %d\n",
             cell_vtx_tria[0], cell_vtx_tria[1], cell_vtx_tria[2]);
  bft_printf("face 2 : %d %d %d\n",
             cell_vtx_tria[3], cell_vtx_tria[4], cell_vtx_tria[5]);
  bft_printf("face 3 : %d %d %d %d\n",
             cell_vtx_quad[0], cell_vtx_quad[1],
             cell_vtx_quad[2], cell_vtx_quad[3]);
  bft_printf("face 4 : %d %d %d %d\n",
             cell_vtx_quad[4], cell_vtx_quad[5],
             cell_vtx_quad[6], cell_vtx_quad[7]);
  bft_printf("face 5 : %d %d %d %d\n",
             cell_vtx_quad[8], cell_vtx_quad[9],
             cell_vtx_quad[10], cell_vtx_quad[11]);
#endif

  /*
    Base : vertices of the first triangle; we take the vertices in opposite
    order ("bottom" face numbering of the prism with outward
    pointing normal).

    *  4 x-------x 6
    *    |\     /|
    *    | \   / |
    *  1 x- \-/ -x 3
    *     \ 5x  /
    *      \ | /
    *       \|/
    *        x 2
  */


  cell_vtx_prism[0] = cell_vtx_tria[2];
  cell_vtx_prism[1] = cell_vtx_tria[1];
  cell_vtx_prism[2] = cell_vtx_tria[0];

  /*
    We have found 3 out of 6 vertices; we first seek the quadrangle sharing
    vertices 1 and 2, so as to determine vertices 4 and 5,
    then we seek the quadrangle sharing vertices 2 and 3, so as to
    determine vertices 5 and 6.
  */

  for (ipass = 0 ; ipass < 2 ; ipass++) {

    vtx_num_1 = cell_vtx_prism[    ipass];
    vtx_num_2 = cell_vtx_prism[1 + ipass];

    direction = 0;

    for (face_id = 0 ; face_id < 4 ; face_id++) {

      for (vertex_id = 0 ; vertex_id < 4 ; vertex_id++) {

        vtx_num = cell_vtx_quad[face_id*4 + vertex_id];

        if (vtx_num == vtx_num_1) {
          if (cell_vtx_quad[face_id*4 + ((vertex_id+1) % 4)] == vtx_num_2) {
            direction = 1;
            break;
          }
          else if (   cell_vtx_quad[face_id*4 + ((vertex_id-1+4) % 4)]
                   == vtx_num_2) {
            direction = -1;
            break;
          }
        }

      }

      if (direction != 0)
        break;

    }

    if (direction == -1)
      warn_orient = true;
    else if (direction == 0)
      return FVM_NODAL_FROM_DESC_FAILURE;

    cell_vtx_prism[3 + ipass]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (3 * direction))) % 4)];

    cell_vtx_prism[4 + ipass]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (2 * direction))) % 4)];

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("prism : %d %d %d %d %d %d\n",
             cell_vtx_prism[0], cell_vtx_prism[1], cell_vtx_prism[2],
             cell_vtx_prism[3], cell_vtx_prism[4], cell_vtx_prism[5]);
#endif

  if (warn_orient == true)
    return FVM_NODAL_FROM_DESC_WORIENT;
  else
    return FVM_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Hexahedron construction
 *
 * parameters:
 *   cell_vtx_quad <-- quadrangle faces connectivity
 *   cell_vtx_hexa --> hexahedron connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvm_nodal_from_desc_t
_nodal_from_desc_cnv_cel_hexa(const cs_lnum_t   cell_vtx_quad[],
                              cs_lnum_t   cell_vtx_hexa[])
{
  cs_lnum_t   vertex_id, face_id;
  cs_lnum_t   ipass, direction;
  cs_lnum_t   vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("face 1 : %d %d %d %d\n",
             cell_vtx_quad[0], cell_vtx_quad[1],
             cell_vtx_quad[2], cell_vtx_quad[3]);
  bft_printf("face 2 : %d %d %d %d\n",
             cell_vtx_quad[4], cell_vtx_quad[5],
             cell_vtx_quad[6], cell_vtx_quad[7]);
  bft_printf("face 3 : %d %d %d %d\n",
             cell_vtx_quad[8], cell_vtx_quad[9],
             cell_vtx_quad[10], cell_vtx_quad[11]);
  bft_printf("face 4 : %d %d %d %d\n",
             cell_vtx_quad[12], cell_vtx_quad[13],
             cell_vtx_quad[14], cell_vtx_quad[15]);
  bft_printf("face 5 : %d %d %d %d\n",
             cell_vtx_quad[16], cell_vtx_quad[17],
             cell_vtx_quad[18], cell_vtx_quad[19]);
  bft_printf("face 6 : %d %d %d %d\n",
             cell_vtx_quad[20], cell_vtx_quad[21],
             cell_vtx_quad[22], cell_vtx_quad[23]);
#endif

  /*
    Base : vertices of the fisrt face; we take the vertices in opposite
    order ("bottom" face numbering of the pyramid with outward
    pointing normal).

    *     8 x-------x 7
    *      /|      /|
    *     / |     / |
    *  5 x-------x6 |
    *    | 4x----|--x 3
    *    | /     | /
    *    |/      |/
    *  1 x-------x 2
  */


  cell_vtx_hexa[0] = cell_vtx_quad[3];
  cell_vtx_hexa[1] = cell_vtx_quad[2];
  cell_vtx_hexa[2] = cell_vtx_quad[1];
  cell_vtx_hexa[3] = cell_vtx_quad[0];

  /*
    We have found 4 out of 8 vertices; we first seek the quadrangle  sharing
    vertices 1 and 2, so as to determine vertices 5 and 6,
    then we seek the quadrangle sharing vertices 3 and 4, so as to
    determine vertices 7 and 8.
  */

  for (ipass = 0 ; ipass < 2 ; ipass++) {

    vtx_num_1 = cell_vtx_hexa[     ipass * 2 ];
    vtx_num_2 = cell_vtx_hexa[1 + (ipass * 2)];

    direction = 0;

    for (face_id = 1 ; face_id < 6 ; face_id++) {

      for (vertex_id = 0 ; vertex_id < 4 ; vertex_id++) {

        vtx_num = cell_vtx_quad[face_id*4 + vertex_id];

        if (vtx_num == vtx_num_1) {
          if (cell_vtx_quad[face_id*4 + ((vertex_id+1) % 4)] == vtx_num_2) {
            direction = 1;
            break;
          }
          else if (   cell_vtx_quad[face_id*4 + ((vertex_id-1+4) % 4)]
                   == vtx_num_2) {
            direction = -1;
            break;
          }
        }

      }

      if (direction != 0)
        break;

    }

    if (direction == -1)
      warn_orient = true;
    else if (direction == 0)
      return FVM_NODAL_FROM_DESC_FAILURE;

    cell_vtx_hexa[4 + (ipass * 2)]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (3 * direction))) % 4)];

    cell_vtx_hexa[5 + (ipass * 2)]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (2 * direction))) % 4)];

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("hexa : %d %d %d %d %d %d %d %d\n",
             cell_vtx_hexa[0], cell_vtx_hexa[1], cell_vtx_hexa[2],
             cell_vtx_hexa[3], cell_vtx_hexa[4], cell_vtx_hexa[5],
             cell_vtx_hexa[6], cell_vtx_hexa[7]);
#endif

  if (warn_orient == true)
    return FVM_NODAL_FROM_DESC_WORIENT;
  else
    return FVM_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Determination of the number of vertices defining a given face.
 *
 * parameters:
 *   face_id,        <-- face id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *
 * returns:
 *   number of vertices of face defined by face_id
 *----------------------------------------------------------------------------*/

inline static cs_lnum_t
_nodal_face_from_desc_size(const cs_lnum_t    face_id,
                           const int          n_face_lists,
                           const cs_lnum_t    face_list_shift[],
                           const cs_lnum_t   *face_vertex_idx[])
{
  cs_lnum_t   fl, _face_id;
  cs_lnum_t   n_face_vertices;
  cs_lnum_t   vertex_id_start, vertex_id_end;

  /* Compute number of vertices */
  /*----------------------------*/

  _face_id = face_id;
  for (fl = n_face_lists - 1 ; _face_id < face_list_shift[fl] ; fl--);
  assert(fl > -1);
  _face_id -= face_list_shift[fl];

  vertex_id_start = face_vertex_idx[fl][_face_id];
  vertex_id_end   = face_vertex_idx[fl][_face_id + 1];
  n_face_vertices = vertex_id_end - vertex_id_start;

  return n_face_vertices;

}

/*----------------------------------------------------------------------------
 * Copy the vertex ids defining a given face.
 *
 * The face_vtx[] array is filled with the vertex ids of the face,
 * and should be large enough to receive them.
 *
 * parameters:
 *   face_id,        <-- face id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   face_vtx        --> local face definition
 *----------------------------------------------------------------------------*/

inline static void
_nodal_face_from_desc_copy(const cs_lnum_t    face_id,
                           const int          n_face_lists,
                           const cs_lnum_t    face_list_shift[],
                           const cs_lnum_t   *face_vertex_idx[],
                           const cs_lnum_t   *face_vertex[],
                           cs_lnum_t         *face_vtx)
{
  cs_lnum_t   fl, _face_id;
  cs_lnum_t   vtx, vertex_id, vertex_id_start, vertex_id_end;

  /* Copy vertex ids */
  /*-----------------*/

  _face_id = face_id;
  for (fl = n_face_lists - 1 ; _face_id < face_list_shift[fl] ; fl--);
  assert(fl > -1);
  _face_id -= face_list_shift[fl];

  vertex_id_start = face_vertex_idx[fl][_face_id];
  vertex_id_end   = face_vertex_idx[fl][_face_id + 1];

  for (vtx = 0, vertex_id = vertex_id_start ;
       vertex_id < vertex_id_end ;
       vtx++, vertex_id++)
    face_vtx[vtx] = face_vertex[fl][vertex_id] + 1;
}

/*----------------------------------------------------------------------------
 * Polyhedral cells connectivity extraction.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_polys         <-- size of list_poly[]
 *   list_poly       <-- list of polyhedral cells (cell index, 1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   cell_face_list  --> numbers of faces defining polyhedra
 *----------------------------------------------------------------------------*/

static void
_fvm_nodal_extract_polyhedra(fvm_nodal_section_t  *this_section,
                             const cs_lnum_t       n_polys,
                             const cs_lnum_t       list_poly[],
                             const int             n_face_lists,
                             const cs_lnum_t       face_list_shift[],
                             const cs_lnum_t      *face_vertex_idx[],
                             const cs_lnum_t      *face_vertex[],
                             const cs_lnum_t       cell_face_idx[],
                             const cs_lnum_t       cell_face_num[],
                             cs_lnum_t            *cell_face_list[])
{
  cs_lnum_t    n_faces, n_cell_faces;
  cs_lnum_t    c_cell_face_vertex_idxs, c_cell_face_vertex_nums;
  cs_lnum_t    n_cell_face_vertex_nums;
  cs_lnum_t    face_counter, face_id, poly_id, cell_id;
  cs_lnum_t    fl, i, idx, idx_start, idx_end, num_count;

  cs_lnum_t   *local_face_num = NULL;

  int sgn;


  /* Indicators initialization */

  n_faces = face_list_shift[n_face_lists] - face_list_shift[0];

  BFT_MALLOC(local_face_num, n_faces, cs_lnum_t);

  for (face_id = 0 ; face_id < n_faces ; face_id++)
    local_face_num[face_id] = 0;

  /* Flagging of referenced faces and Cells -> Faces indexes */
  /*---------------------------------------------------------*/

  n_cell_faces = 0;

  BFT_MALLOC(this_section->_face_index, n_polys + 1, cs_lnum_t);
  this_section->face_index = this_section->_face_index;
  this_section->_face_index[0] = 0;

  for (poly_id = 0 ; poly_id < n_polys ; poly_id++) {

    cell_id = list_poly[poly_id] - 1;

    idx_start = cell_face_idx[cell_id]     - 1;
    idx_end   = cell_face_idx[cell_id  +1] - 1;

    this_section->_face_index[poly_id + 1]
      = this_section->_face_index[poly_id] + (idx_end - idx_start);

    for (idx = idx_start ; idx < idx_end ; idx++) {

      face_id = CS_ABS(cell_face_num[idx]) - 1;

      /* Mark only used values for now, local_face_num[] values
         will be computed later by looping on faces so that
         for nonzero values, local_face_num[i] > local_face_num[j]
         if i > j ; this is important for later parts of this algorithm */

      if (local_face_num[face_id] == 0)
        local_face_num[face_id] = 1;

    }

  }

  /* Counting for faces -> vertices connectivity and local face numbering */

  n_cell_face_vertex_nums = 0;
  for (face_counter = 0 ; face_counter < n_faces ; face_counter++) {

    if (local_face_num[face_counter] != 0) {

      /* Transform local_face_num[] from a marker to a renumbering array */

      n_cell_faces += 1;
      local_face_num[face_counter] = n_cell_faces;

      /* Counting for faces -> vertices connectivity */

      face_id = face_counter;
      for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
      assert(fl > -1);
      face_id -= face_list_shift[fl];

      n_cell_face_vertex_nums += (  face_vertex_idx[fl][face_id+1]
                                  - face_vertex_idx[fl][face_id]);

    }
  }

  /* Cells -> Faces Connectivity (face numbers) */
  /*--------------------------------------------*/

  BFT_MALLOC(this_section->_face_num,
             this_section->_face_index[n_polys],
             cs_lnum_t);
  this_section->face_num = this_section->_face_num;

  num_count = 0;

  for (poly_id = 0 ; poly_id < n_polys ; poly_id++) {

    cell_id = list_poly[poly_id] - 1;

    idx_start = cell_face_idx[cell_id]     - 1;
    idx_end   = cell_face_idx[cell_id  +1] - 1;

    for (idx = idx_start ; idx < idx_end ; idx++) {

      face_id = CS_ABS(cell_face_num[idx]) - 1;
      sgn = (cell_face_num[idx] > 0) ? 1 : -1;

      this_section->_face_num[num_count++]
        = sgn * local_face_num[face_id] ;

    }

  }

  assert(num_count == this_section->_face_index[n_polys]);

  /* Faces -> Vertices Connectivity */
  /*--------------------------------*/

  if (cell_face_list != NULL)
    BFT_MALLOC(*cell_face_list, n_cell_faces, cs_lnum_t);

  BFT_MALLOC(this_section->_vertex_index,
             n_cell_faces + 1,
             cs_lnum_t);
  this_section->vertex_index = this_section->_vertex_index;

  BFT_MALLOC(this_section->_vertex_num,
             n_cell_face_vertex_nums,
             cs_lnum_t);
  this_section->vertex_num = this_section->_vertex_num;

  /* Definition */

  c_cell_face_vertex_idxs = 0;
  c_cell_face_vertex_nums = 0;
  this_section->_vertex_index[0] = 0;

  for (face_counter = 0 ; face_counter < n_faces ; face_counter++) {

    if (local_face_num[face_counter] != 0) {

      if (cell_face_list != NULL)
        (*cell_face_list)[local_face_num[face_counter] -1 ] = face_counter + 1;

      face_id = face_counter;
      for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
      assert(fl < n_face_lists);
      face_id -= face_list_shift[fl];

      for (i = face_vertex_idx[fl][face_id];
           i < face_vertex_idx[fl][face_id + 1];
           i++)
        this_section->_vertex_num[c_cell_face_vertex_nums++]
          = face_vertex[fl][i] + 1;

      c_cell_face_vertex_idxs += 1;
      this_section->_vertex_index[c_cell_face_vertex_idxs]
        = c_cell_face_vertex_nums;

    }
  }

  /* No past-the-end index value counted,
     so counter = this_nodal->n_cell_faces */

  assert(c_cell_face_vertex_idxs == n_cell_faces);
  assert(c_cell_face_vertex_nums == n_cell_face_vertex_nums);

  /* Set pointers and connectivity size */

  this_section->connectivity_size = n_cell_face_vertex_nums;

  this_section->n_faces      = n_cell_faces;
  this_section->face_index   = this_section->_face_index;
  this_section->face_num     = this_section->_face_num;
  this_section->vertex_index = this_section->_vertex_index;
  this_section->vertex_num   = this_section->_vertex_num;

  /* Free memory */

  BFT_FREE(local_face_num);

}

/*----------------------------------------------------------------------------
 * Raise parent element numbering in given sections by one level, as
 * defined by the parent_element_num[] arrays.
 *
 * This is useful when a nodal mesh is extracted from a temporary subset of
 * a parent mesh (corresponding to the parent_element_num[] list, using
 * 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   n_sections         <-- size of sections array
 *   sections           <-- array of sections to reduce
 *   parent_element_num <-- element -> parent element number (1 to n) if
 *                          non-trivial (i.e. if element definitions
 *                          correspond to a subset of the parent mesh),
 *                          NULL otherwise.
 *----------------------------------------------------------------------------*/

static void
_raise_sections_parent_num(const int             n_sections,
                           fvm_nodal_section_t  *sections[],
                           const cs_lnum_t       parent_element_num[])
{
  int  section_id;
  cs_lnum_t   element_counter;

  fvm_nodal_section_t  *section;

  if (parent_element_num == NULL)
    return;

  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL) {
      if (section->_parent_element_num == NULL) {
        BFT_MALLOC(section->_parent_element_num,
                   section->n_elements,
                   cs_lnum_t);
        section->parent_element_num = section->_parent_element_num;
      }
      for (element_counter = 0 ;
           element_counter < section->n_elements ;
           element_counter++)
        section->_parent_element_num[element_counter]
          = parent_element_num[section->parent_element_num[element_counter]
                               - 1];
    }
  }

}

/*----------------------------------------------------------------------------
 * Free parent element numbering arrays when not needed.
 *
 * If the parent mesh is already complete and ordered for a given section
 * (i.e. not based on a partial extraction or then based on the n first
 * elements in the same order), the parent element number is not needed,
 * and may be freed.
 *
 * parameters:
 *   sections        <-- array of sections to reduce
 *   n_sections      <-- size of sections array
 *----------------------------------------------------------------------------*/

static void
_optimize_sections_parent_num(const int             n_sections,
                              fvm_nodal_section_t  *sections[])
{
  int  section_id;
  cs_lnum_t   element_counter;

  fvm_nodal_section_t  *section;

  /* If the parent mesh is already complete for a given element type
     (i.e. not based on a partial extraction or then based on the
     n first elements), the parent cell number is not needed */

  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL) {
      for (element_counter = 0 ;
           element_counter < section->n_elements ;
           element_counter++) {
        if (section->parent_element_num[element_counter] != element_counter + 1)
          break;
      }
      if (element_counter == section->n_elements) {
        section->parent_element_num = NULL;
        BFT_FREE(section->_parent_element_num);
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Add nodal mesh structure sections to a nodal mesh.
 *
 * Sections to add are defined by an array of section pointers,
 * which may contain NULL entries. Only the non-NULL entries are
 * added to the nodal mesh structure, and they belong to that structure
 * from this point on (i.e. their lifecycle is based upon that of
 * the nodal mesh structure).
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_sections      <-- size of sections array
 *   sections        <-- array of sections to add
 *----------------------------------------------------------------------------*/

static void
_fvm_nodal_add_sections(fvm_nodal_t          *this_nodal,
                        const int             n_sections,
                        fvm_nodal_section_t  *sections[])
{
  int  section_id, section_count;

  fvm_nodal_section_t  *section;

  /* Add sections to nodal mesh structure */

  section_count = 0;
  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL)
      section_count++;
  }

  BFT_REALLOC(this_nodal->sections,
              this_nodal->n_sections + section_count,
              fvm_nodal_section_t *);

  section_count = 0;
  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL)
      this_nodal->sections[this_nodal->n_sections + section_count++]
        = section;
  }
  this_nodal->n_sections += section_count;

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert and add cells from an descending connectivity mesh to a nodal mesh.
 *
 * If the optional filter list extr_cells[] argument is non-NULL, cells
 * {extr_cells[0], extr_cells[1], extr_cells[n_extr_cells - 1]} are converted
 * and added to the nodal mesh. If this filter is set to NULL, cells
 * {1, 2, ..., n_extr_cells} are considered.
 *
 * In addition, an optional parent_cell_num[] array may also be given, in
 * case the descending connectivity mesh definition is based on a temporary
 * subset of a parent mesh, (corresponding to the parent_cell_num[] list,
 * using 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_extr_cells    <-- count of cells to add
 *   extr_cells      <-- optional filter list of cells to extract (1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists + 1
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   cell_gc_id      <-- cell -> group class ids, or NULL
 *   parent_cell_num <-- cell -> parent cell number (1 to n) if non-trivial
 *                       (i.e. if cell definitions correspond to a subset
 *                       of the parent mesh), NULL otherwise.
 *   cell_face_list  --> numbers of faces defining polyhedra
 *----------------------------------------------------------------------------*/

void
fvm_nodal_from_desc_add_cells(fvm_nodal_t        *this_nodal,
                              const cs_lnum_t     n_extr_cells,
                              const cs_lnum_t     extr_cells[],
                              const int           n_face_lists,
                              const cs_lnum_t     face_list_shift[],
                              const cs_lnum_t    *face_vertex_idx[],
                              const cs_lnum_t    *face_vertex[],
                              const cs_lnum_t     cell_face_idx[],
                              const cs_lnum_t     cell_face_num[],
                              const int           cell_gc_id[],
                              const cs_lnum_t     parent_cell_num[],
                              cs_lnum_t          *cell_face_list[])
{
  int  type_id;
  cs_lnum_t   cell_counter, cell_id;

  fvm_element_t  cell_type;

  cs_lnum_t   cell_vtx_tria[3*4]; /* We will only seek to fill these arrays */
  cs_lnum_t   cell_vtx_quad[4*6]; /* for local faces 1-4 and 1-6 at most    */
  cs_lnum_t   *p_cell_vertex;

  cs_lnum_t   n_elements_type[FVM_N_ELEMENT_TYPES];
  cs_gnum_t   n_g_elements_type[FVM_N_ELEMENT_TYPES];

  fvm_nodal_section_t  *section;
  fvm_nodal_section_t  *sections[FVM_N_ELEMENT_TYPES];

  cs_gnum_t   n_orient_pbs = 0; /* Number of cells with potential (non-
                                   fatal) orientation problem */

  fvm_nodal_from_desc_t  retcode;

  /* Initialization */

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++) {
    n_elements_type[type_id] = 0;
    sections[type_id] = NULL;
  }

  /* Guess connectivity types */
  /*--------------------------*/

  for (cell_counter = 0 ; cell_counter < n_extr_cells ; cell_counter++) {

    if (extr_cells != NULL)
      cell_id = extr_cells[cell_counter] - 1;
    else
      cell_id = cell_counter;

    cell_type = _nodal_cell_from_desc(cell_id,
                                      n_face_lists,
                                      face_list_shift,
                                      face_vertex_idx,
                                      face_vertex,
                                      cell_face_idx,
                                      cell_face_num,
                                      NULL,
                                      NULL);

    n_elements_type[cell_type] += 1;

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("dbg : fvm_nodal_cells_from_desc\n"
             "n_tetra = %d\n"
             "n_pyram = %d\n"
             "n_prism = %d\n"
             "n_hexa  = %d\n"
             "n_poly  = %d\n",
             n_elements_type[FVM_CELL_TETRA],
             n_elements_type[FVM_CELL_PYRAM],
             n_elements_type[FVM_CELL_PRISM],
             n_elements_type[FVM_CELL_HEXA],
             n_elements_type[FVM_CELL_POLY]);
#endif

  /* Set dimensions (and reset local counters) */

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++)
    n_g_elements_type[type_id] = n_elements_type[type_id];

  cs_parall_counter(n_g_elements_type, FVM_N_ELEMENT_TYPES);

  for (cell_type = FVM_CELL_TETRA ;
       cell_type <= FVM_CELL_POLY ;
       cell_type++) {
    if (n_g_elements_type[cell_type] > 0) {
      sections[cell_type] = fvm_nodal_section_create(cell_type);
      sections[cell_type]->n_elements = n_elements_type[cell_type];
      this_nodal->n_cells += n_elements_type[cell_type];
    }
    n_elements_type[cell_type] = 0;
  }

  /* Main memory allocations */

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      if (   section->type != FVM_FACE_POLY
          && section->type != FVM_CELL_POLY) {
        section->stride = fvm_nodal_n_vertices_element[type_id];
        section->connectivity_size = section->stride * section->n_elements;
        BFT_MALLOC(section->_vertex_num, section->connectivity_size, cs_lnum_t);
        section->vertex_num = section->_vertex_num;
      }
    }
  }

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      BFT_MALLOC(section->_parent_element_num, section->n_elements, cs_lnum_t);
      section->parent_element_num = section->_parent_element_num;
    }
  }

  /* Construction of nodal connectivities */
  /*--------------------------------------*/

  for (cell_counter = 0 ; cell_counter < n_extr_cells ; cell_counter++) {

    if (extr_cells != NULL)
      cell_id = extr_cells[cell_counter] - 1;
    else
      cell_id = cell_counter;

    cell_type = _nodal_cell_from_desc(cell_id,
                                      n_face_lists,
                                      face_list_shift,
                                      face_vertex_idx,
                                      face_vertex,
                                      cell_face_idx,
                                      cell_face_num,
                                      cell_vtx_tria,
                                      cell_vtx_quad);

    section = sections[cell_type];

    p_cell_vertex =   section->_vertex_num
                    + (  n_elements_type[cell_type]
                       * fvm_nodal_n_vertices_element[cell_type]);

    switch (cell_type) {
    case FVM_CELL_TETRA:
      retcode = _nodal_from_desc_cnv_cel_tetra(cell_vtx_tria,
                                               p_cell_vertex);
      break;
    case FVM_CELL_PYRAM:
      retcode = _nodal_from_desc_cnv_cel_pyram(cell_vtx_tria,
                                               cell_vtx_quad,
                                               p_cell_vertex);
      break;
    case FVM_CELL_PRISM:
      retcode = _nodal_from_desc_cnv_cel_prism(cell_vtx_tria,
                                               cell_vtx_quad,
                                               p_cell_vertex);
      break;
    case FVM_CELL_HEXA:
      retcode = _nodal_from_desc_cnv_cel_hexa(cell_vtx_quad,
                                              p_cell_vertex);
      break;
    default:
      retcode = FVM_NODAL_FROM_DESC_SUCCESS;
      break;
    }

    /* Temporary value of parent cell num based on local cell id
       (so that the list of polyhedra given as an argument to
       _fvm_nodal_extract_polyhedra() below is based on the local
       numbering, like the cell->face connectivity */

    section->_parent_element_num[n_elements_type[cell_type]] = cell_id + 1;

    n_elements_type[cell_type] += 1;

    if (retcode == FVM_NODAL_FROM_DESC_WORIENT)
      n_orient_pbs += 1;
    else
      if (retcode == FVM_NODAL_FROM_DESC_FAILURE)
        bft_error(__FILE__, __LINE__, 0,
                  _("Incoherent connectivity for cell %d\n"),
                  cell_id);

  }

  cs_parall_counter(&n_orient_pbs, 1);

  if (n_orient_pbs > 0)
    bft_printf("Warning: Possible nodal connectivity orientation\n"
               "         problems for at least %llu cells\n",
               (unsigned long long)n_orient_pbs);

  /* Extraction of remaining polyhedra */
  /*-----------------------------------*/

  if (sections[FVM_CELL_POLY] != NULL)
    _fvm_nodal_extract_polyhedra
      (sections[FVM_CELL_POLY],
       n_elements_type[FVM_CELL_POLY],
       sections[FVM_CELL_POLY]->parent_element_num,
       n_face_lists,
       face_list_shift,
       face_vertex_idx,
       face_vertex,
       cell_face_idx,
       cell_face_num,
       cell_face_list);

  /* We can now base the final value of the parent cell number on
     the parent (and not local) numbering */

  _raise_sections_parent_num(FVM_N_ELEMENT_TYPES, sections, parent_cell_num);

  _optimize_sections_parent_num(FVM_N_ELEMENT_TYPES, sections);

  /* Add group class ids if present */

  if (cell_gc_id != NULL) {
    int section_id;
    for (section_id = 0; section_id < FVM_N_ELEMENT_TYPES; section_id++) {
      section = sections[section_id];
      if (section == NULL)
        continue;
      BFT_MALLOC(section->gc_id, section->n_elements, int);
      if (section->parent_element_num != NULL) {
        for (cell_id = 0; cell_id < section->n_elements; cell_id++)
          section->gc_id[cell_id]
            = cell_gc_id[section->parent_element_num[cell_id] - 1];
      }
      else
        memcpy(section->gc_id, cell_gc_id, section->n_elements*sizeof(int));
    }
  }

  /* Add sections to nodal mesh structure */

  _fvm_nodal_add_sections(this_nodal, FVM_N_ELEMENT_TYPES, sections);
}

/*----------------------------------------------------------------------------
 * Convert and add faces from an descending connectivity mesh to a nodal mesh.
 *
 * If the optional filter list extr_faces[] argument is non-NULL, faces
 * {extr_faces[0], extr_faces[1], extr_faces[n_extr_faces - 1]} are converted
 * and added to the nodal mesh. If this filter is set to NULL, faces
 * {1, 2, ..., n_extr_faces} are considered.
 *
 * In addition, an optional parent_face_num[] array may also be given, in
 * case the descending connectivity mesh definition is based on a temporary
 * subset of a parent mesh, (corresponding to the parent_face_num[] list,
 * using 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_extr_faces    <-- count of faces to add
 *   extr_faces      <-- optional filter list of faces to extract (1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   face_gc_id      <-- face -> group class ids, or NULL (per face list)
 *   parent_face_num <-- face -> parent face number (1 to n) if non-trivial
 *                       (i.e. if face definitions correspond to a subset
 *                       of the parent mesh), NULL otherwise.
 *----------------------------------------------------------------------------*/

void
fvm_nodal_from_desc_add_faces(fvm_nodal_t        *this_nodal,
                              const cs_lnum_t     n_extr_faces,
                              const cs_lnum_t     extr_faces[],
                              const int           n_face_lists,
                              const cs_lnum_t     face_list_shift[],
                              const cs_lnum_t    *face_vertex_idx[],
                              const cs_lnum_t    *face_vertex[],
                              const int          *face_gc_id[],
                              const cs_lnum_t     parent_face_num[])
{
  int  type_id;
  cs_lnum_t   face_counter, face_id;

  fvm_element_t  face_type;

  cs_lnum_t   *p_vertex_idx;
  cs_lnum_t   *p_vertex_num;

  cs_lnum_t   poly_connect_size;
  cs_lnum_t   n_face_vertices;
  cs_lnum_t   n_elements_type[FVM_N_ELEMENT_TYPES];
  cs_gnum_t   n_g_elements_type[FVM_N_ELEMENT_TYPES];

  fvm_nodal_section_t  *section;
  fvm_nodal_section_t  *sections[FVM_N_ELEMENT_TYPES];

  /* Initialization */

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++) {
    n_elements_type[type_id] = 0;
    sections[type_id] = NULL;
  }
  poly_connect_size = 0;

  /* Compute connectivity type */
  /*---------------------------*/

  for (face_counter = 0 ; face_counter < n_extr_faces ; face_counter++) {

    if (extr_faces != NULL)
      face_id = extr_faces[face_counter] - 1;
    else
      face_id = face_counter;

    n_face_vertices = _nodal_face_from_desc_size(face_id,
                                                 n_face_lists,
                                                 face_list_shift,
                                                 face_vertex_idx);

    switch (n_face_vertices) {
    case 3:
      face_type = FVM_FACE_TRIA;
      break;
    case 4:
      face_type = FVM_FACE_QUAD;
      break;
    default:
      face_type = FVM_FACE_POLY;
      poly_connect_size += n_face_vertices;
      break;
    }

    n_elements_type[face_type] += 1;

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("dbg : fvm_nodal_faces_from_desc\n"
             "n_tria = %d\n"
             "n_quad = %d\n"
             "n_poly = %d\n",
             n_elements_type[FVM_FACE_TRIA],
             n_elements_type[FVM_FACE_QUAD],
             n_elements_type[FVM_FACE_POLY]);
#endif

  /* Set dimensions (and reset local counters) */

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++)
    n_g_elements_type[type_id] = n_elements_type[type_id];

  cs_parall_counter(n_g_elements_type, FVM_N_ELEMENT_TYPES);

  for (face_type = FVM_FACE_TRIA ;
       face_type <= FVM_FACE_POLY ;
       face_type++) {
    if (n_g_elements_type[face_type] > 0) {
      sections[face_type] = fvm_nodal_section_create(face_type);
      sections[face_type]->n_elements = n_elements_type[face_type];
      this_nodal->n_faces += n_elements_type[face_type];
    }
    n_elements_type[face_type] = 0;
  }

  /* Main memory allocations */

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      if (section->type != FVM_FACE_POLY) {
        section->stride = fvm_nodal_n_vertices_element[type_id];
        section->connectivity_size = section->stride * section->n_elements;
        BFT_MALLOC(section->_vertex_num, section->connectivity_size, cs_lnum_t);
        section->vertex_num = section->_vertex_num;
      }
      else {
        section->stride = fvm_nodal_n_vertices_element[type_id];
        section->connectivity_size = poly_connect_size;
        BFT_MALLOC(section->_vertex_index, section->n_elements + 1, cs_lnum_t);
        BFT_MALLOC(section->_vertex_num, section->connectivity_size, cs_lnum_t);
        section->vertex_index = section->_vertex_index;
        section->vertex_num = section->_vertex_num;
        section->_vertex_index[0] = 0;
      }
    }
  }

  for (type_id = 0 ; type_id < FVM_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      BFT_MALLOC(section->_parent_element_num, section->n_elements, cs_lnum_t);
      section->parent_element_num = section->_parent_element_num;
    }
  }

  /* Construction of nodal connectivities */
  /*---------------------------------------*/

  for (face_counter = 0 ; face_counter < n_extr_faces ; face_counter++) {

    if (extr_faces != NULL)
      face_id = extr_faces[face_counter] - 1;
    else
      face_id = face_counter;

    n_face_vertices = _nodal_face_from_desc_size(face_id,
                                                 n_face_lists,
                                                 face_list_shift,
                                                 face_vertex_idx);

    switch (n_face_vertices) {
    case 3:
      face_type = FVM_FACE_TRIA;
      section = sections[face_type];
      p_vertex_num = section->_vertex_num + (n_elements_type[face_type] * 3);
      break;
    case 4:
      face_type = FVM_FACE_QUAD;
      section = sections[face_type];
      p_vertex_num = section->_vertex_num + (n_elements_type[face_type] * 4);
      break;
    default:
      face_type = FVM_FACE_POLY;
      section = sections[face_type];
      p_vertex_idx = section->_vertex_index + n_elements_type[face_type];
      *(p_vertex_idx + 1) = *p_vertex_idx + n_face_vertices;
      p_vertex_num = section->_vertex_num + (*p_vertex_idx);
      break;
    }

    _nodal_face_from_desc_copy(face_id,
                               n_face_lists,
                               face_list_shift,
                               face_vertex_idx,
                               face_vertex,
                               p_vertex_num);

    section->_parent_element_num[n_elements_type[face_type]] = face_id + 1;

    n_elements_type[face_type] += 1;

  }

  /* We can now base the final value of the parent face number on
     the parent (and not local) numbering */

  _raise_sections_parent_num(FVM_N_ELEMENT_TYPES, sections, parent_face_num);

  _optimize_sections_parent_num(FVM_N_ELEMENT_TYPES, sections);

  /* Add group class ids if present */

  if (face_gc_id != NULL) {

    int section_id;

    for (section_id = 0; section_id < FVM_N_ELEMENT_TYPES; section_id++) {

      section = sections[section_id];
      if (section == NULL)
        continue;

      BFT_MALLOC(section->gc_id, section->n_elements, int);

      if (section->parent_element_num != NULL) {
        for (face_id = 0; face_id < section->n_elements; face_id++) {
          int fl;
          cs_lnum_t _face_id = section->parent_element_num[face_id] - 1;
          for (fl = n_face_lists - 1 ; _face_id < face_list_shift[fl] ; fl--);
          assert(fl > -1);
          _face_id -= face_list_shift[fl];
          section->gc_id[face_id] = face_gc_id[fl][_face_id];
        }
      }

      else {
        for (face_id = 0; face_id < section->n_elements; face_id++) {
          int fl;
          cs_lnum_t _face_id = face_id;
          for (fl = n_face_lists - 1 ; _face_id < face_list_shift[fl] ; fl--);
          assert(fl > -1);
          _face_id -= face_list_shift[fl];
          section->gc_id[face_id] = face_gc_id[fl][_face_id];
        }
      }

    }
  }

  /* Add sections to nodal mesh structure */

  _fvm_nodal_add_sections(this_nodal, FVM_N_ELEMENT_TYPES, sections);

}

/*----------------------------------------------------------------------------
 * Determination of a given cell's type.
 *
 * If the optional cell_vtx[8] array is given, it is filled with the vertex
 * indexes of cell's vertices, unless the cell is a general polyhedron.
 *
 * parameters:
 *   cell_id         <-- cell id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   vertex_num      --> nodal connectivity of cell, if not a general
 *                       polyhedron
 *
 * returns:
 *   type of cell defined by cell_id
 *----------------------------------------------------------------------------*/

fvm_element_t
fvm_nodal_from_desc_cell(const cs_lnum_t    cell_id,
                         const int          n_face_lists,
                         const cs_lnum_t    face_list_shift[],
                         const cs_lnum_t   *face_vertex_idx[],
                         const cs_lnum_t   *face_vertex[],
                         const cs_lnum_t    cell_face_idx[],
                         const cs_lnum_t    cell_face_num[],
                         cs_lnum_t          vertex_num[8])
{
  fvm_nodal_from_desc_t  retcode;
  cs_lnum_t   cell_vtx_tria[3*4]; /* We will only seek to fill these arrays */
  cs_lnum_t   cell_vtx_quad[4*6]; /* for local faces 1-4 and 1-6 at most    */

  fvm_element_t cell_type = _nodal_cell_from_desc(cell_id,
                                                  n_face_lists,
                                                  face_list_shift,
                                                  face_vertex_idx,
                                                  face_vertex,
                                                  cell_face_idx,
                                                  cell_face_num,
                                                  cell_vtx_tria,
                                                  cell_vtx_quad);

  switch (cell_type) {
  case FVM_CELL_TETRA:
    retcode = _nodal_from_desc_cnv_cel_tetra(cell_vtx_tria,
                                             vertex_num);
    break;
  case FVM_CELL_PYRAM:
    retcode = _nodal_from_desc_cnv_cel_pyram(cell_vtx_tria,
                                             cell_vtx_quad,
                                             vertex_num);
    break;
  case FVM_CELL_PRISM:
    retcode = _nodal_from_desc_cnv_cel_prism(cell_vtx_tria,
                                             cell_vtx_quad,
                                             vertex_num);
    break;
  case FVM_CELL_HEXA:
    retcode = _nodal_from_desc_cnv_cel_hexa(cell_vtx_quad,
                                            vertex_num);
    break;
  default:
    retcode = FVM_NODAL_FROM_DESC_FAILURE;
    break;
  }

  if (retcode != FVM_NODAL_FROM_DESC_SUCCESS)
    cell_type = FVM_CELL_POLY;

  return cell_type;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
