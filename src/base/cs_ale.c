/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_interface.h"

#include "cs_base.h"

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ale.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell and face center of gravity, cell volume.
 *
 * Fortran Interface
 *
 * subroutine algrma
 * *****************
 *
 * min_vol           : --> : Minimum cell volume
 * max_vol           : --> : Maximum cell volume
 * tot_vol           : --> : Total mesh volume
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algrma, ALGRMA)(cs_real_t  *min_vol,
                          cs_real_t  *max_vol,
                          cs_real_t  *tot_vol)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_mesh_quantities_compute(m, mq);
  cs_mesh_bad_cells_detect(m, mq);

  *min_vol = mq->min_vol;
  *max_vol = mq->max_vol;
  *tot_vol = mq->tot_vol;
}

/*----------------------------------------------------------------------------
 * Projection on mesh vertices of the displacement (computed on cell center)
 *
 * Fortran Interface
 *
 * subroutine aledis
 * *****************
 *
 * ialtyb            : <-- : Type of boundary for ALE
 * meshv             : <-- : Mesh velocity
 * gradm             : <-- : Mesh velocity gradient (du_i/dx_j : gradv[][i][j])
 * claale            : <-- : Boundary conditions A
 * clbale            : <-- : Boundary conditions B
 * dt                : <-- : Time step
 * disp_proj         : --> : Displacement projected on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (aledis, ALEDIS)(const cs_int_t      ialtyb[],
                          const cs_real_3_t  *meshv,
                          const cs_real_33_t  gradm[],
                          const cs_real_3_t  *claale,
                          const cs_real_33_t *clbale,
                          const cs_real_t    *dt,
                          cs_real_3_t        *disp_proj)
{
  cs_int_t  j, face_id, vtx_id, cell_id, cell_id1, cell_id2;

  bool *vtx_interior_indicator = NULL;

  cs_real_t *vtx_counter = NULL;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_int_t  n_vertices = m->n_vertices;
  const cs_int_t  n_cells = m->n_cells;
  const cs_int_t  n_b_faces = m->n_b_faces;
  const cs_int_t  n_i_faces = m->n_i_faces;
  const cs_int_t  dim = m->dim;
  const cs_real_3_t *restrict vtx_coord
    = (const cs_real_3_t *restrict)m->vtx_coord;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict face_cen
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  BFT_MALLOC(vtx_counter, n_vertices, cs_real_t);
  BFT_MALLOC(vtx_interior_indicator, n_vertices, bool);

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++) {

    vtx_counter[vtx_id] = 0.;
    vtx_interior_indicator[vtx_id] = true;

    for (int i = 0; i < dim; i++)
      disp_proj[vtx_id][i] = 0.;

  }

  /* All nodes wich belongs to a boundary face where the
     displacement is imposed (that is all faces except sliding BCs)
     are boundary nodes, the others are interior nodes. */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    if (ialtyb[face_id] != 2) {

      for (j = m->b_face_vtx_idx[face_id];
           j < m->b_face_vtx_idx[face_id+1];
           j++) {


        vtx_id = m->b_face_vtx_lst[j];
        vtx_interior_indicator[vtx_id] = false;

      } /* End of loop on vertices of the face */

    }

  } /* End of loop on border faces */


  /* Interior face and nodes treatment */

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cell_id1 = m->i_face_cells[face_id][0];
    cell_id2 = m->i_face_cells[face_id][1];

    cs_real_t dvol1 = 1./fvq->cell_vol[cell_id1];
    cs_real_t dvol2 = 1./fvq->cell_vol[cell_id2];

    if (cell_id1 < n_cells) { /* Test to take into account face only once */

      for (j = m->i_face_vtx_idx[face_id];
           j < m->i_face_vtx_idx[face_id+1];
           j++) {

        /* Get the vertex number */

        vtx_id = m->i_face_vtx_lst[j];

        if (vtx_interior_indicator[vtx_id]) {

          /* Get the vector from the cell center to the node */

          cs_real_3_t cen1_node;
          cs_real_3_t cen2_node;
          for (int i = 0; i < 3; i++) {
            cen1_node[i] = vtx_coord[vtx_id][i]-cell_cen[cell_id1][i];
            cen2_node[i] = vtx_coord[vtx_id][i]-cell_cen[cell_id2][i];
          }

          for (int i = 0; i < 3; i++) {
            disp_proj[vtx_id][i] +=
              dvol1*(meshv[cell_id1][i] + gradm[cell_id1][i][0]*cen1_node[0]
                                        + gradm[cell_id1][i][1]*cen1_node[1]
                                        + gradm[cell_id1][i][2]*cen1_node[2])
              *dt[cell_id1]
             +dvol2*(meshv[cell_id2][i] + gradm[cell_id2][i][0]*cen2_node[0]
                                        + gradm[cell_id2][i][1]*cen2_node[1]
                                        + gradm[cell_id2][i][2]*cen2_node[2])
             *dt[cell_id2];
          }

          vtx_counter[vtx_id] += dvol1+dvol2;

        } /* End of Interior nodes */

      }

    }

  } /* End of loop on internal faces */

  /* Border face treatment.
     only border face contribution */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    cell_id = m->b_face_cells[face_id];

    for (j = m->b_face_vtx_idx[face_id];
         j < m->b_face_vtx_idx[face_id+1];
         j++) {

      vtx_id = m->b_face_vtx_lst[j];

      if (!vtx_interior_indicator[vtx_id]) {

        /* Get the vector from the face center to the node*/

        cs_real_3_t face_node;
        for (int i = 0; i<3; i++)
          face_node[i] = -face_cen[face_id][i] + vtx_coord[vtx_id][i];

        /* 1st order extrapolation of the mesh velocity at the face center
         * to the node */

        cs_real_3_t vel_node;
        for (int i = 0; i<3; i++)
          vel_node[i] = claale[face_id][i]
                      + gradm[cell_id][i][0]*face_node[0]
                      + gradm[cell_id][i][1]*face_node[1]
                      + gradm[cell_id][i][2]*face_node[2];

        cs_real_t dsurf = 1./fvq->b_face_surf[face_id];

        for (int i = 0; i<3; i++)
          disp_proj[vtx_id][i] += dsurf*dt[cell_id]*
            (vel_node[i] + clbale[face_id][i][0]*meshv[cell_id][0]
                         + clbale[face_id][i][1]*meshv[cell_id][1]
                         + clbale[face_id][i][2]*meshv[cell_id][2]);

        vtx_counter[vtx_id] += dsurf;

      } /* End of boundary nodes */

    } /* End of loop on vertices of the face */

  } /* End of loop on border faces */


  /* If the boundary face IS a sliding face.
     We project the displacment paralelly to the face. */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    if (ialtyb[face_id] == 2) {

      for (j = m->b_face_vtx_idx[face_id];
           j < m->b_face_vtx_idx[face_id+1];
           j++) {


        vtx_id = m->b_face_vtx_lst[j];
        disp_proj[vtx_id][0] = clbale[face_id][0][0]*disp_proj[vtx_id][0]
                             + clbale[face_id][0][1]*disp_proj[vtx_id][1]
                             + clbale[face_id][0][2]*disp_proj[vtx_id][2];
        disp_proj[vtx_id][1] = clbale[face_id][1][0]*disp_proj[vtx_id][0]
                             + clbale[face_id][1][1]*disp_proj[vtx_id][1]
                             + clbale[face_id][1][2]*disp_proj[vtx_id][2];
        disp_proj[vtx_id][2] = clbale[face_id][2][0]*disp_proj[vtx_id][0]
                             + clbale[face_id][2][1]*disp_proj[vtx_id][1]
                             + clbale[face_id][2][2]*disp_proj[vtx_id][2];

      } /* End of loop on vertices of the face */

    }

  } /* End of loop on border faces */

  if (m->vtx_interfaces != NULL) {
    cs_interface_set_sum(m->vtx_interfaces,
                         n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         disp_proj);
    cs_interface_set_sum(m->vtx_interfaces,
                         n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         vtx_counter);
  }

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++)
    for (int i = 0; i < dim; i++)
      disp_proj[vtx_id][i] /= vtx_counter[vtx_id];

  BFT_FREE(vtx_counter);
  BFT_FREE(vtx_interior_indicator);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
