/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Management of mesh quantities
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_perio.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

enum {X, Y, Z};

#define _CS_CROSS_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define _CS_DOT_PRODUCT(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define _CS_MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to cs_mesh_quantities_t structure for the main mesh */

cs_mesh_quantities_t  *cs_glob_mesh_quantities = NULL;

/* Choice of the algorithm for computing gravity centres of the cells */

static int cs_glob_mesh_quantities_cell_cen = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute quantities associated to faces (border or internal)
 *
 * parameters:
 *   dim             <--  dimension
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx_lst    <--  "face -> vertices" connectivity list
 *   face_normal     -->  surface normal of the face
 *
 *
 *                          Pi+1
 *              *---------*                   B  : barycentre of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_compute_face_normal(cs_int_t         dim,
                     cs_int_t         n_faces,
                     const cs_real_t  vtx_coord[],
                     const cs_int_t   face_vtx_idx[],
                     const cs_int_t   face_vtx_lst[],
                     cs_real_t        face_normal[])
{
  cs_int_t  i, face_id, tri_id, vtx_id, start_id, end_id, shift;
  cs_int_t  n_face_vertices, n_max_face_vertices;
  cs_point_t  this_face_normal, this_face_barycentre;
  cs_point_t  vect1, vect2;

  cs_point_t  *face_vtx_coord = NULL;
  cs_point_t  *triangle_normal = NULL;

  /* Return if there is not enough data (some SolCom meshes) */

  if (face_vtx_idx == NULL || face_vtx_lst == NULL)
    return;

  /* Checking */

  if (dim != 3)
    bft_error(__FILE__, __LINE__,0,
              _("Face geometric quantities computation is only\n"
                "implemented in 3D."));

  assert(face_normal != NULL || n_faces == 0);

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (face_id = 0; face_id < n_faces; face_id++) {
    n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices + 1, cs_point_t);
  BFT_MALLOC(triangle_normal, n_max_face_vertices, cs_point_t);

  /* Loop on each face */

  for (face_id = 0; face_id < n_faces; face_id++) {

    /* Initialization */

    n_face_vertices = 0;

    start_id = face_vtx_idx[face_id] - 1;
    end_id = face_vtx_idx[face_id + 1] - 1;

    for (i = 0; i < 3; i++)
      this_face_normal[i] = 0.0;

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

      shift = 3 * (face_vtx_lst[vtx_id] - 1);
      for (i = 0; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
      n_face_vertices++;

    }

    for (i = 0; i < 3; i++)
      face_vtx_coord[n_face_vertices][i] = face_vtx_coord[0][i];

    /* Compute the barycentre of the face */

    for (i = 0; i < 3; i++) {

      this_face_barycentre[i] = 0.0;
      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        this_face_barycentre[i] += face_vtx_coord[vtx_id][i];
      this_face_barycentre[i] /= n_face_vertices;

    }

    /* Loop on the triangles of the face (defined by an edge of the face
       and its barycentre) */

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      /*----------------------------------------------------------------------*/
      /* Computation of the normal of each triangle Ti :                      */
      /*                                                                      */
      /*  ->            -->   -->                                             */
      /*  N(Ti) = 1/2 ( BPi X BPi+1 )                                         */
      /*----------------------------------------------------------------------*/

      for (i = 0; i < 3; i++) {
        vect1[i] = face_vtx_coord[tri_id    ][i] - this_face_barycentre[i];
        vect2[i] = face_vtx_coord[tri_id + 1][i] - this_face_barycentre[i];
      }

      _CS_CROSS_PRODUCT(triangle_normal[tri_id], vect1, vect2);

      for (i = 0; i < 3; i++)
        triangle_normal[tri_id][i] *= 0.5;

      /*----------------------------------------------------------------------*/
      /* Computation of the normal of the polygon                             */
      /*  => vectorial sum of normals of each triangle                        */
      /*                                                                      */
      /*  ->      n-1   ->                                                    */
      /*  N(P) =  Sum ( N(Ti) )                                               */
      /*          i=0                                                         */
      /*----------------------------------------------------------------------*/

      for (i = 0; i < 3; i++)
        this_face_normal[i] += triangle_normal[tri_id][i];

    } /* End of loop on triangles of the face */

    /* Store result in appropriate structure */

    for (i = 0; i < 3; i++)
      face_normal[face_id * 3 + i] = this_face_normal[i];

  } /* End of loop on faces */

  BFT_FREE(triangle_normal);
  BFT_FREE(face_vtx_coord);

}

/*----------------------------------------------------------------------------
 * Compute quantities associated to faces (border or internal)
 *
 * parameters:
 *   dim             <--  dimension
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx_lst    <--  "face -> vertices" connectivity list
 *   face_cog        -->  coordinates of the centre of gravity of the faces
 *   face_norm       -->  face surface normals
 *   face_surf       -->  face surfaces (optional), or NULL
 *
 *                          Pi+1
 *              *---------*                   B  : barycentre of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_compute_face_quantities(const cs_int_t   dim,
                         const cs_int_t   n_faces,
                         const cs_real_t  vtx_coord[],
                         const cs_int_t   face_vtx_idx[],
                         const cs_int_t   face_vtx_lst[],
                         cs_real_t        face_cog[],
                         cs_real_t        face_norm[],
                         cs_real_t        face_surf[])
{
  cs_int_t  i, fac_id, tri_id;
  cs_int_t  vtx_id, lower_vtx_id, upper_vtx_id;
  cs_int_t  n_face_vertices, n_max_face_vertices;
  cs_int_t  lower_coord_id;
  cs_real_t  face_surface, tri_surface;
  cs_real_t  face_vol_part, tri_vol_part, rectif_cog;
  cs_point_t  face_barycentre, face_normal;
  cs_point_t  face_centre, tri_centre;
  cs_point_t  vect1, vect2;

  cs_point_t  *face_vtx_coord = NULL;
  cs_point_t  *triangle_norm = NULL;

  /* Return if there is not enough data (some SolCom meshes) */

  if (face_vtx_idx == NULL || face_vtx_lst == NULL)
    return;

  /* Checking */

  if (dim != 3)
    bft_error(__FILE__, __LINE__,0,
              _("Face geometric quantities computation is only\n"
                "implemented in 3D."));

  assert(face_cog != NULL || n_faces == 0);
  assert(face_norm != NULL || n_faces == 0);

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (fac_id = 0; fac_id < n_faces; fac_id++) {
    n_face_vertices = face_vtx_idx[fac_id + 1] - face_vtx_idx[fac_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices + 1, cs_point_t);
  BFT_MALLOC(triangle_norm, n_max_face_vertices, cs_point_t);

  /*=========================================================================*/
  /* Loop on faces                                                           */
  /*=========================================================================*/

  for (fac_id = 0; fac_id < n_faces; fac_id++) {

    tri_vol_part = 0.;
    face_surface = 0.0;

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    lower_vtx_id = face_vtx_idx[fac_id] - 1;
    upper_vtx_id = face_vtx_idx[fac_id + 1] - 1;

    n_face_vertices = 0;

    for (vtx_id = lower_vtx_id; vtx_id < upper_vtx_id; vtx_id++) {

      lower_coord_id = 3 * (face_vtx_lst[vtx_id] - 1);

      for (i = X; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[lower_coord_id + i];

      n_face_vertices++;

    }

    for (i = X; i < 3; i++)
      face_vtx_coord[n_face_vertices][i] = face_vtx_coord[0][i];

    /*------------------------------------------------------------------------
     * Compute barycentre (B) coordinates for the polygon (P)
     *
     *  -->    1   n-1  -->
     *  OB  =  -  Somme OPi
     *         n   i=0
     *------------------------------------------------------------------------*/

    for (i = X; i < 3; i++) {

      face_barycentre[i] = 0.0;

      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        face_barycentre[i] += face_vtx_coord[vtx_id][i];

      face_barycentre[i] /= n_face_vertices;

    }

    for (i = X; i < 3; i++) {
      face_normal[i] = 0.0;
      face_centre[i] = 0.0;
    }

    /* First loop on triangles of the face (computation of surface normals)   */
    /*========================================================================*/

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      /*----------------------------------------------------------------------
       * Computation of the normal of the triangle Ti :
       *
       *  ->            -->   -->
       *  N(Ti) = 1/2 ( BPi X BPi+1 )
       *----------------------------------------------------------------------*/

      for (i = X; i < 3; i++) {
        vect1[i] = face_vtx_coord[tri_id    ][i] - face_barycentre[i];
        vect2[i] = face_vtx_coord[tri_id + 1][i] - face_barycentre[i];
      }

      _CS_CROSS_PRODUCT(triangle_norm[tri_id], vect1, vect2);

      for (i = X; i < 3; i++)
        triangle_norm[tri_id][i] *= 0.5;

      /*----------------------------------------------------------------------
       * Computation of the normal of the polygon
       *  => vector sum of normals of triangles
       *
       *  ->      n-1   ->
       *  N(P) =  Sum ( N(Ti) )
       *          i=0
       *----------------------------------------------------------------------*/

      for (i = X; i < 3; i++)
        face_normal[i] += triangle_norm[tri_id][i];

    } /* End of loop on triangles of the face */

    /* Second loop on triangles of the face (for the barycentre)        */
    /*==================================================================*/

    for (tri_id = 0; tri_id < n_face_vertices; tri_id++) {

      /*----------------------------------------------------------------------
       * Compation of the gravity centre G(Ti) of each triangle (Ti)
       *
       *  -->            -->  -->   -->
       *  OG(Ti) = 1/3 ( OB + OPi + OPi+1 )
       *
       * And their part in the volume of Ti
       *
       *  -->    ->
       *  OG(Ti).N(Ti)
       *----------------------------------------------------------------------*/

      for (i = X; i < 3; i++) {

        tri_centre[i] = face_barycentre[i]
                      + face_vtx_coord[tri_id    ][i]
                      + face_vtx_coord[tri_id + 1][i];

        tri_centre[i] /= 3.0;

        tri_vol_part += (tri_centre[i] * triangle_norm[tri_id][i]);

      }

      /*----------------------------------------------------------------------
       * Computation of the area of Ti (norm of the surface normal)
       *
       *               ->
       *  Surf(Ti) = | N(Ti) |
       *----------------------------------------------------------------------*/

      tri_surface = _CS_MODULE(triangle_norm[tri_id]);

      if (_CS_DOT_PRODUCT(triangle_norm[tri_id], face_normal) < 0.0)
        tri_surface *= -1.0;

      face_surface += tri_surface;

      /*----------------------------------------------------------------------
       *   n-1
       *   Sum  Surf(Ti) G(Ti)
       *   i=0
       *----------------------------------------------------------------------*/

      for (i = X; i < 3; i++)
        face_centre[i] += tri_surface * tri_centre[i];

    } /* End of second loop  on triangles of the face */

    /*------------------------------------------------------------------------
     * Compute the centre of gravity G(P) of the polygon P :
     *
     *           n-1
     *           Sum  Surf(Ti) G(Ti)
     *           i=0
     *  G(P) = -----------------------
     *           n-1
     *           Sum  Surf(Ti)
     *           i=0
     *
     * Computation of the part of volume of the polygon (before rectification)
     *
     *  -->    ->
     *  OG(P).N(P)
     *------------------------------------------------------------------------*/

    face_vol_part = 0.0;

    for (i = X; i < 3; i++) {
      face_centre[i] = face_centre[i] / face_surface;
      face_vol_part += (face_centre[i] * face_normal[i]);
    }

    rectif_cog = (tri_vol_part - face_vol_part) / (face_surface * face_surface);

    for (i = X; i < 3; i++)
      face_centre[i] += rectif_cog * face_normal[i];

    /* Store result in appropriate structure */

    for (i = X; i < 3; i++) {
      face_cog[fac_id * 3 + i] = face_centre[i];
      face_norm[fac_id * 3 + i] = face_normal[i];
    }

  } /* End of loop on faces */

  BFT_FREE(triangle_norm);
  BFT_FREE(face_vtx_coord);

  if (face_norm == NULL || face_surf == NULL)
    return;

  if (dim != 3)
    bft_error(__FILE__, __LINE__,0,
              _("Face surface computation is only\n"
                "implemented in 3D."));

  /* Compute optional face surfaces */
  /*--------------------------------*/

  if (face_surf != NULL) {

    for (fac_id = 0; fac_id < n_faces; fac_id++) {
      double nx = face_norm[fac_id*3];
      double ny = face_norm[fac_id*3+1];
      double nz = face_norm[fac_id*3+2];
      face_surf[fac_id] = sqrt(nx*nx + ny*ny + nz*nz);
    }
  }
}

/*----------------------------------------------------------------------------*
 * Compute centre of gravity of cells C from their vertices S(i) where
 * i=0, n-1
 *
 * -->      1    n-1   -->
 * OB(C) = ---   Sum   OSi
 *          n    i=0
 *
 * parameters:
 *   mesh          <--  pointer to mesh structure
 *   cell_cen      -->  centre of gravity of cells
 *----------------------------------------------------------------------------*/

static void
_compute_cell_cen_vertex(const cs_mesh_t  *const mesh,
                         cs_real_t         cell_cen[])
{
  cs_int_t  i, j, k, cell_id, fac_id, vtx_id;
  cs_int_t  face_num, vtx_counter;

  cs_int_t  *vtx_tag = NULL;
  cs_int_t  *_face_vtx_idx = NULL, *_face_vtx_lst = NULL;
  cs_int_t  *cell_faces_idx = NULL, *cell_faces_lst = NULL;

  /* Return if there is not enough data (Solcom case except radiative module) */

  if (mesh->i_face_vtx_lst == NULL && mesh->b_face_vtx_lst == NULL)
    return;

  /* Checking */

  if (mesh->dim != 3)
    bft_error(__FILE__, __LINE__,0,
              _("Cell centre computation is only implemented in 3D."));

  assert(cell_cen != NULL);

  /* Allocation and initialization */

  BFT_MALLOC(vtx_tag, mesh->n_vertices, cs_int_t);

  for (vtx_id = 0 ; vtx_id < mesh->n_vertices ; vtx_id++)
    vtx_tag[vtx_id] = -1;

  /* Initialization */

  for (i = 0; i < 3*mesh->n_cells_with_ghosts; i++)
    cell_cen[i] = 0.0;

  /* Extract "cell -> faces" connectivity */

  cs_mesh_connect_get_cell_faces(mesh,
                                 mesh->n_cells,
                                 NULL,
                                 &cell_faces_idx,
                                 &cell_faces_lst);

  /* Loop on cells */
  /* ------------- */

  for (cell_id = 0; cell_id < mesh->n_cells; cell_id++) {

    vtx_counter = 0;

    /* Loop on faces of the cell */

    for (j = cell_faces_idx[cell_id]; j < cell_faces_idx[cell_id + 1]; j++) {

      face_num = cell_faces_lst[j - 1];

      /* Internal or border face */

      if (face_num > mesh->n_b_faces) {
        fac_id = face_num - mesh->n_b_faces - 1;
        _face_vtx_idx = mesh->i_face_vtx_idx;
        _face_vtx_lst = mesh->i_face_vtx_lst;
      }
      else {
        fac_id = face_num - 1;
        _face_vtx_idx = mesh->b_face_vtx_idx;
        _face_vtx_lst = mesh->b_face_vtx_lst;
      }

      /* Loop on vertices of the face */

      for (k = _face_vtx_idx[fac_id]-1; k < _face_vtx_idx[fac_id + 1]-1; k++) {

        vtx_id = _face_vtx_lst[k] - 1;

        if (vtx_tag[vtx_id] < cell_id) {
          for (i = 0 ; i < 3 ; i++)
            cell_cen[cell_id*3 + i] += mesh->vtx_coord[vtx_id*3 + i];
          vtx_counter += 1;
          vtx_tag[vtx_id] = cell_id;
        }

      }

    } /* End of loop on faces of the cell */

    for (i = 0; i < 3; i++)
      cell_cen[cell_id*3 + i] /= (double)vtx_counter;

  } /* End of loop on cells */

  /* Free memory */

  BFT_FREE(vtx_tag);
  BFT_FREE(cell_faces_idx);
  BFT_FREE(cell_faces_lst);

}

/*----------------------------------------------------------------------------*
 * Compute centre of gravity of cells C from their faces F(i) where i=0, n-1
 *
 *           n-1
 *           Sum  Surf(Fi) G(Fi)
 *           i=0
 *  G(C) = -----------------------
 *           n-1
 *           Sum  Surf(Fi)
 *           i=0
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   i_face_norm    <--  surface normal of internal faces
 *   i_face_cog     <--  centre of gravity of internal faces
 *   b_face_norm    <--  surface normal of border faces
 *   b_face_cog     <--  centre of gravity of border faces
 *   cell_cen       -->  centre of gravity of cells
 *----------------------------------------------------------------------------*/

static void
_compute_cell_cen_face(const cs_mesh_t  *const mesh,
                       const cs_real_t         i_face_norm[],
                       const cs_real_t         i_face_cog[],
                       const cs_real_t         b_face_norm[],
                       const cs_real_t         b_face_cog[],
                             cs_real_t         cell_cen[])
{
  cs_int_t  i, j, fac_id, cell_id, cell_id1, cell_id2;
  cs_real_t  area;
  cs_real_t  _norm[3];

  cs_real_t  *cell_area = NULL;

  /* Mesh connectivity */

  const  cs_int_t  dim = mesh->dim;
  const  cs_int_t  n_i_faces = mesh->n_i_faces;
  const  cs_int_t  n_b_faces = mesh->n_b_faces;
  const  cs_int_t  n_cells = mesh->n_cells;
  const  cs_int_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;
  const  cs_int_t  *i_face_cells = mesh->i_face_cells;
  const  cs_int_t  *b_face_cells = mesh->b_face_cells;

  /* Return if ther is not enough data (Solcom case except rediative module
     or Pre-processor 1.2.d without option "-n") */

  if (mesh->i_face_vtx_lst == NULL && mesh->b_face_vtx_lst == NULL)
    return;

  /* Checking */

  if (dim != 3)
    bft_error(__FILE__, __LINE__,0,
              _("Cell centre computation is only implemented in 3D."));

  assert(cell_cen != NULL);

  /* Initialization */

  BFT_MALLOC(cell_area, n_cells_with_ghosts, cs_real_t);

  for (j = 0; j < n_cells_with_ghosts; j++) {

    cell_area[j] = 0.;

    for (i = 0; i < dim; i++)
      cell_cen[dim*j + i] = 0. ;

  }

  /* ---------------------- */
  /* Loop on internal faces */
  /* ---------------------- */

  for (fac_id = 0; fac_id < n_i_faces; fac_id++) {

    /* ----------------------------------------------------------
     * For each cell sharing the internal face, we update
     * cell_cen and cell_area
     * ---------------------------------------------------------- */

    cell_id1 = i_face_cells[2*fac_id    ] - 1;
    cell_id2 = i_face_cells[2*fac_id + 1] - 1;

    /* Computation of the area of the face */

    for (i = 0; i < dim; i++)
      _norm[i] = i_face_norm[dim*fac_id + i];

    area = _CS_MODULE(_norm);

    cell_area[cell_id1] += area;
    cell_area[cell_id2] += area;

    /* Computation of the numerator */

    for (i = 0; i < dim; i++) {
      cell_cen[dim*cell_id1 + i] += i_face_cog[dim*fac_id + i]*area;
      cell_cen[dim*cell_id2 + i] += i_face_cog[dim*fac_id + i]*area;
    }

  } /* End of loop on internal faces */

  /* -------------------- */
  /* Loop on border faces */
  /* -------------------- */

  for (fac_id = 0; fac_id < n_b_faces; fac_id++) {

    /* -------------------------------------------------------------
     * For each cell sharing a border face, we update the numerator
     * of cell_cen and cell_area
     * ------------------------------------------------------------- */

    cell_id1 = b_face_cells[fac_id] - 1;

    /* Computation of the area of the face */

    for (i = 0; i < dim; i++)
      _norm[i] = b_face_norm[dim*fac_id + i];

    area = _CS_MODULE(_norm);
    cell_area[cell_id1] += area;

    /* Computation of the numerator */

    for (i = 0; i < dim; i++)
      cell_cen[dim*cell_id1 + i] += b_face_cog[dim*fac_id + i]*area;

  } /* End of loop on border faces */

  /* ------------------------------------------------------------------
   * Loop on cells to finalize the computation of centre of gravity
   * ------------------------------------------------------------------*/

  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    for (i = 0; i < dim; i++)
      cell_cen[cell_id*dim + i] /= cell_area[cell_id];

  } /* End of loop on cells */

  /* Free memory */

  BFT_FREE(cell_area);

}

/*----------------------------------------------------------------------------
 * Compute the volume of cells C from their n faces F(i) and their center of
 * gravity G(Fi) where i=0, n-1
 *
 *         1    n-1
 *  G(C) = - .  Sum  Surf(Fi) G(Fi)
 *         3    i=0
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   i_face_norm    <--  surface normal of internal faces
 *   i_face_cog     <--  centre of gravity of internal faces
 *   b_face_norm    <--  surface normal of border faces
 *   b_face_cog     <--  centre of gravity of border faces
 *   cell_cen       -->  centre of gravity of cells
 *----------------------------------------------------------------------------*/

static void
_compute_cell_volume(const cs_mesh_t  *const mesh,
                     const cs_real_t         i_face_norm[],
                     const cs_real_t         i_face_cog[],
                     const cs_real_t         b_face_norm[],
                     const cs_real_t         b_face_cog[],
                           cs_real_t         cell_vol[])
{
  cs_int_t  id1, id2, i, fac_id, cell_id;

  cs_real_t  flux = 0;

  const cs_real_t  a_third = 1.0/3.0;
  const cs_int_t  dim = mesh->dim;

  /* Initialization */

  for (cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
    cell_vol[cell_id] = 0;

  /* Loop on internal faces */

  for (fac_id = 0; fac_id < mesh->n_i_faces; fac_id++) {

    id1 = mesh->i_face_cells[2*fac_id] - 1;
    id2 = mesh->i_face_cells[2*fac_id + 1] - 1;

    flux = 0;
    for (i = 0; i < dim; i++)
      flux += i_face_norm[dim*fac_id + i] * i_face_cog[dim*fac_id + i];

    cell_vol[id1] += flux;
    cell_vol[id2] -= flux;

  }

  /* Loop on border faces */

  for (fac_id = 0; fac_id < mesh->n_b_faces; fac_id++) {

    id1 = mesh->b_face_cells[fac_id] - 1;

    flux = 0;
    for (i = 0; i < dim; i++)
      flux += b_face_norm[dim*fac_id + i] * b_face_cog[dim*fac_id + i];

    cell_vol[id1] += flux;

  }

  /* Computation of the volume */

  for (cell_id = 0; cell_id < mesh->n_cells; cell_id++)
    cell_vol[cell_id] *= a_third;

}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for computing cell centers.
 *
 * This function returns 1 or 2 according to the selected algorithm.
 *
 * Fortran interface :
 *
 * SUBROUTINE ALGCEN (IOPT)
 * *****************
 *
 * INTEGER          IOPT        : <-> : Choice of the algorithm
 *                                      < 0 : query
 *                                        0 : computation based
 *                                            on faces (default choice)
 *                                        1 : computation based
 *                                            on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algcen, ALGCEN) (cs_int_t  *const iopt)
{
  int  iopt_ret = cs_mesh_quantities_cell_cen_choice((int)(*iopt));

  *iopt = iopt_ret;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for computing cell centers.
 *
 *  < 0 : query
 *    0 : computation based on faces (default choice)
 *    1 : computation based on vertices
 *
 * algo_choice  <--  choice of algorithm to compute cell centers.
 *
 * returns:
 *  1 or 2 according to the selected algorithm.
 *----------------------------------------------------------------------------*/

int
cs_mesh_quantities_cell_cen_choice(const int algo_choice)
{
  if (algo_choice > 1)
    bft_error(__FILE__, __LINE__,0,
              _("The algorithm selection indicator for the cell centre of gravity computation\n"
                "can take the following values:\n"
                "  0: computation based on the face centres and surfaces\n"
                "  1: computation based on the vertices\n"
                "and not %d."), cs_glob_mesh_quantities_cell_cen);

  else if (algo_choice >= 0)
    cs_glob_mesh_quantities_cell_cen = algo_choice;

  return cs_glob_mesh_quantities_cell_cen;
}

/*----------------------------------------------------------------------------
 * Create a mesh quantities structure.
 *
 * returns:
 *   pointer to created cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

cs_mesh_quantities_t  *
cs_mesh_quantities_create(void)
{
  cs_mesh_quantities_t  *mesh_quantities = NULL;

  BFT_MALLOC(mesh_quantities, 1, cs_mesh_quantities_t);

  mesh_quantities->cell_cen = NULL;
  mesh_quantities->cell_vol = NULL;
  mesh_quantities->i_face_normal = NULL;
  mesh_quantities->b_face_normal = NULL;
  mesh_quantities->i_face_cog = NULL;
  mesh_quantities->b_face_cog = NULL;
  mesh_quantities->i_face_surf = NULL;
  mesh_quantities->b_face_surf = NULL;

  return (mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Destroy a mesh quantities structure
 *
 * mesh_quantities <-- pointer to a cs_mesh_quantities_t structure
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

cs_mesh_quantities_t *
cs_mesh_quantities_destroy(cs_mesh_quantities_t  *mesh_quantities)
{

  BFT_FREE(mesh_quantities->cell_cen);
  BFT_FREE(mesh_quantities->cell_vol);
  BFT_FREE(mesh_quantities->i_face_normal);
  BFT_FREE(mesh_quantities->b_face_normal);
  BFT_FREE(mesh_quantities->i_face_cog);
  BFT_FREE(mesh_quantities->b_face_cog);
  BFT_FREE(mesh_quantities->i_face_surf);
  BFT_FREE(mesh_quantities->b_face_surf);

  BFT_FREE(mesh_quantities);

  return (mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Compute mesh quantities
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute(const cs_mesh_t       *mesh,
                           cs_mesh_quantities_t  *mesh_quantities)
{
  cs_int_t  i, j;

  cs_int_t  dim = mesh->dim;
  cs_int_t  n_i_faces = mesh->n_i_faces;
  cs_int_t  n_b_faces = mesh->n_b_faces;
  cs_int_t  n_cells = mesh->n_cells;
  cs_int_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  /* If this is not an update, allocate members of the structure */

  if (mesh_quantities->i_face_normal == NULL)
    BFT_MALLOC(mesh_quantities->i_face_normal, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->i_face_cog == NULL)
    BFT_MALLOC(mesh_quantities->i_face_cog, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->b_face_normal == NULL)
    BFT_MALLOC(mesh_quantities->b_face_normal, n_b_faces*dim, cs_real_t);

  if (mesh_quantities->b_face_cog == NULL)
    BFT_MALLOC(mesh_quantities->b_face_cog, n_b_faces*dim, cs_real_t);

  if (mesh_quantities->cell_cen == NULL)
    BFT_MALLOC(mesh_quantities->cell_cen, n_cells_with_ghosts*dim, cs_real_t);

  if (mesh_quantities->cell_vol == NULL)
    BFT_MALLOC(mesh_quantities->cell_vol, n_cells_with_ghosts, cs_real_t);

  if (mesh_quantities->i_face_surf == NULL)
    BFT_MALLOC(mesh_quantities->i_face_surf, n_i_faces, cs_real_t);

  if (mesh_quantities->b_face_surf == NULL)
    BFT_MALLOC(mesh_quantities->b_face_surf, n_b_faces, cs_real_t);

  /* Compute centres of gravity, normals, and surfaces of interior faces */

  _compute_face_quantities(dim,
                           n_i_faces,
                           mesh->vtx_coord,
                           mesh->i_face_vtx_idx,
                           mesh->i_face_vtx_lst,
                           mesh_quantities->i_face_cog,
                           mesh_quantities->i_face_normal,
                           mesh_quantities->i_face_surf);

  /* Compute centres of gravity, normals, and surfaces of boundary faces */

  _compute_face_quantities(dim,
                           n_b_faces,
                           mesh->vtx_coord,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           mesh_quantities->b_face_cog,
                           mesh_quantities->b_face_normal,
                           mesh_quantities->b_face_surf);

  /* Compute cell centers from face barycentres or vertices */

  switch (cs_glob_mesh_quantities_cell_cen) {

  case 0:
    _compute_cell_cen_face(mesh,
                           mesh_quantities->i_face_normal,
                           mesh_quantities->i_face_cog,
                           mesh_quantities->b_face_normal,
                           mesh_quantities->b_face_cog,
                           mesh_quantities->cell_cen);
    break;

  case 1:
    _compute_cell_cen_vertex(mesh,
                             mesh_quantities->cell_cen);
    break;

  default:
    assert(0);

  }

  /* Compute the volume of cells */

  _compute_cell_volume(mesh,
                       mesh_quantities->i_face_normal,
                       mesh_quantities->i_face_cog,
                       mesh_quantities->b_face_normal,
                       mesh_quantities->b_face_cog,
                       mesh_quantities->cell_vol);

  /* Initialize values for cells in halo (standard or extended) */

  for (j = n_cells; j < n_cells_with_ghosts; j++) {

    mesh_quantities->cell_vol[j] = 0.;

    for (i = 0; i < dim; i++)
      mesh_quantities->cell_cen[j*dim + i] = 0.;

  }

  /* Parallel synchronization */

  /* Synchronize geometric quantities */

  if (mesh->halo != NULL) {

    cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                             mesh_quantities->cell_cen, 3);

    cs_halo_sync_var(mesh->halo, CS_HALO_EXTENDED, mesh_quantities->cell_vol);

    if (mesh->n_init_perio > 0) {

       cs_perio_sync_coords(mesh->halo, mesh->halo_type,
                            mesh_quantities->cell_cen);

       cs_perio_sync_var_scal(mesh->halo,
                              mesh->halo_type,
                              CS_PERIO_ROTA_COPY,
                              mesh_quantities->cell_vol);
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute internal and border face normal.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_normal <-> pointer to the internal face normal array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_face_normal(const cs_mesh_t   *mesh,
                               cs_real_t         *p_i_face_normal[],
                               cs_real_t         *p_b_face_normal[])
{
  cs_real_t  *i_face_normal = NULL, *b_face_normal = NULL;

  const cs_int_t  n_b_faces = mesh->n_b_faces;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  dim = mesh->dim;

  /* Internal face treatment */

  BFT_MALLOC(i_face_normal, n_i_faces * dim, cs_real_t);

  _compute_face_normal(dim,
                       mesh->n_i_faces,
                       mesh->vtx_coord,
                       mesh->i_face_vtx_idx,
                       mesh->i_face_vtx_lst,
                       i_face_normal);

  *p_i_face_normal = i_face_normal;

  /* Border face treatment */

  BFT_MALLOC(b_face_normal, n_b_faces * dim, cs_real_t);

  _compute_face_normal(dim,
                       mesh->n_b_faces,
                       mesh->vtx_coord,
                       mesh->b_face_vtx_idx,
                       mesh->b_face_vtx_lst,
                       b_face_normal);

  *p_b_face_normal = b_face_normal;

}

/*----------------------------------------------------------------------------
 * Compute border face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_b_face_cog    <-> pointer to the border face center array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_b_face_cog[],
                           cs_real_t         *p_b_face_normal[])
{
  cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;

  BFT_MALLOC(b_face_cog, mesh->n_b_faces * mesh->dim, cs_real_t);
  BFT_MALLOC(b_face_normal, mesh->n_b_faces * mesh->dim, cs_real_t);

  _compute_face_quantities(mesh->dim,
                           mesh->n_b_faces,
                           mesh->vtx_coord,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           b_face_cog,
                           b_face_normal,
                           NULL);

  *p_b_face_cog = b_face_cog;
  *p_b_face_normal = b_face_normal;
}

/*----------------------------------------------------------------------------
 * Check that no negative volumes are present, and exit on error otherwise.
 *
 * parameters:
 *   mesh            <-- pointer to mesh structure
 *   mesh_quantities <-- pointer to mesh quantities structure
 *   allow_error     <-- 1 if errors are allowed, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_check_vol(const cs_mesh_t             *mesh,
                             const cs_mesh_quantities_t  *mesh_quantities,
                             int                          allow_error)
{
  cs_int_t  cell_id;

  fvm_gnum_t  error_count = 0;

  for (cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
    if (mesh_quantities->cell_vol[cell_id] < 0.0)
      error_count += 1;
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    fvm_gnum_t tot_error_count = 0;
    MPI_Allreduce(&error_count, &tot_error_count, 1, FVM_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    error_count = tot_error_count;
  }
#endif

  /* Exit with error */

  if (error_count > 0) {
    const char fmt[]
      = N_("  %llu cells have a Negative volume.\n"
           " Run mesh quality check for post-processing output.\n"
           " In case of mesh joining, this may be due to overly "
           " agressive joining parameters.");

    if (allow_error) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_(fmt), (unsigned long long)error_count);
      bft_printf("\n\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(fmt), (unsigned long long)error_count);
  }
}

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_quantities_t structure
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-- pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_dump(const cs_mesh_t             *mesh,
                        const cs_mesh_quantities_t  *mesh_quantities)
{
  cs_int_t  i;

  const cs_int_t  n_cells = mesh->n_cells_with_ghosts;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  n_b_faces = mesh->n_b_faces;

  const cs_real_t  *cell_cen = mesh_quantities->cell_cen;
  const cs_real_t  *cell_vol = mesh_quantities->cell_vol;
  const cs_real_t  *i_fac_norm = mesh_quantities->i_face_normal;
  const cs_real_t  *b_fac_norm = mesh_quantities->b_face_normal;
  const cs_real_t  *i_fac_cog = mesh_quantities->i_face_cog;
  const cs_real_t  *b_fac_cog = mesh_quantities->b_face_cog;
  const cs_real_t  *i_fac_surf = mesh_quantities->i_face_surf;
  const cs_real_t  *b_fac_surf = mesh_quantities->b_face_surf;

  bft_printf("\n\nDUMP OF A MESH QUANTITIES STRUCTURE: %p\n\n",
             mesh_quantities);

  if (mesh_quantities == NULL)
    return;

  /* Cell data */

  bft_printf("\n\n"
             "    ---------------"
             "    Cell quantities"
             "    ---------------\n\n");

  bft_printf("Cell center coordinates:\n");
  for (i = 0; i < n_cells; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               cell_cen[3*i], cell_cen[3*i+1], cell_cen[3*i+2]);

  bft_printf("\nCell volume:\n");
  for (i = 0; i < n_cells; i++)
    bft_printf("    < %d >    %.3f\n", i+1, cell_vol[i]);

  /* Internal faces data */

  bft_printf("\n\n"
             "    ------------------------"
             "    Interior face quantities"
             "    ------------------------\n\n");

  bft_printf("\nInterior face normals\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               i_fac_norm[3*i], i_fac_norm[3*i+1], i_fac_norm[3*i+2]);

  bft_printf("\nInterior face centers\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               i_fac_cog[3*i], i_fac_cog[3*i+1], i_fac_cog[3*i+2]);

  bft_printf("\nInterior face surfaces\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %d >    %.3f\n", i+1, i_fac_surf[i]);

  /* Border faces data */

  bft_printf("\n\n"
             "    ------------------------"
             "    Boundary face quantities"
             "    ------------------------\n\n");

  bft_printf("\nBoundary face normals\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               b_fac_norm[3*i], b_fac_norm[3*i+1], b_fac_norm[3*i+2]);

  bft_printf("\nBoundary faces centers\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %d >    %.3f    %.3f    %.3f\n", i+1,
               b_fac_cog[3*i], b_fac_cog[3*i+1], b_fac_cog[3*i+2]);

  bft_printf("\nBoundary face surfaces\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %d >    %.3f\n", i+1, b_fac_surf[i]);

  bft_printf("\n\nEND OF DUMP OF MESH QUANTITIES STRUCTURE\n\n");
  bft_printf_flush();

}

/*----------------------------------------------------------------------------*/

#if 0 /* Test if face orientation is OK */

  cs_int_t   i, fac_id, cell_id;
  cs_real_t  cogfac[3];
  cs_real_t  cogcel[3];
  cs_real_t  normal[3];
  cs_real_t  pscal;

  for (fac_id = 0 ; fac_id < mesh->n_b_faces ; fac_id++) {

    cell_id = mesh->b_face_cells[fac_id] - 1;
    pscal = 0;

    for (i = 0 ; i < 3 ; i++) {
      cogcel[i]  = cs_glob_mesh_quantities->cell_cen[cell_id*3 + i];
      cogfac[i]  = cs_glob_mesh_quantities->b_face_cog[fac_id*3 + i];
      normal[i] = cs_glob_mesh_quantities->b_face_normal[fac_id*3 + i];
      pscal += normal[i] * (cogfac[i] - cogcel[i]);
    }

    if (pscal < 0.0)
      printf("num_fac_brd = %d, num_cel = %d, pscal = %f\n",
             fac_id + 1, cell_id + 1, pscal);
  }

#endif


/*----------------------------------------------------------------------------*/

/* Delete local macro definitions */

#undef _CS_CROSS_PRODUCT
#undef _CS_DOT_PRODUCT
#undef _CS_MODULE

/*----------------------------------------------------------------------------*/

END_C_DECLS
