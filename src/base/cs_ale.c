/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell and face centre of gravity, cell volume.
 *
 * Fortran Interface
 *
 * SUBROUTINE ALGRMA
 * *****************
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algrma, ALGRMA)(void)
{

  cs_mesh_quantities_compute(cs_glob_mesh,
                             cs_glob_mesh_quantities);
  cs_mesh_bad_cells_detect(cs_glob_mesh,
                           cs_glob_mesh_quantities);

}

/*----------------------------------------------------------------------------
 * Projection on mesh vertices of the displacement (computed on cell center)
 *
 * Fortran Interface
 *
 * SUBROUTINE ALDEPL
 * *****************
 *
 * INTEGER         IFACEL(2,NFAC)  : --> : Interior faces -> cells connectivity
 * INTEGER         IFABOR(NFABOR)  : --> : Border faces -> cells connectivity
 * INTEGER         IPNFAC(NFAC+1)  : --> : Interior faces -> vertices index
 * INTEGER         NODFAC(LNDFAC)  : --> : Interior faces -> vertices list
 * INTEGER         IPNFBR(NFABOR+1): --> : Border faces -> vertices index
 * INTEGER         NODFBR(LNDFBR)  : --> : Border faces -> vertices list
 * DOUBLE PRECISION UMA(NCELET)    : --> : Mesh velocity along X
 * DOUBLE PRECISION VMA(NCELET)    : --> : Mesh velocity along Y
 * DOUBLE PRECISION WMA(NCELET)    : --> : Mesh velocity along Z
 * DOUBLE PRECISION COEFAU(NCELET) : --> : Boundary conditions A for UMA
 * DOUBLE PRECISION COEFAV(NCELET) : --> : Boundary conditions A pour VMA
 * DOUBLE PRECISION COEFAW(NCELET) : --> : Boundary conditions A pour WMA
 * DOUBLE PRECISION COEFBU(NCELET) : --> : Boundary conditions B pour UMA
 * DOUBLE PRECISION COEFBV(NCELET) : --> : Boundary conditions B pour VMA
 * DOUBLE PRECISION COEFBW(NCELET) : --> : Boundary conditions B pour WMA
 * DOUBLE PRECISION DT(NCELET)     : --> : Time step
 * DOUBLE PRECISION DEPROJ(NNOD,3)): <-- : Displacement projected on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (aldepl, ALDEPL)(const cs_int_t    i_face_cells[],
                          const cs_int_t    b_face_cells[],
                          const cs_int_t    i_face_vtx_idx[],
                          const cs_int_t    i_face_vtx_lst[],
                          const cs_int_t    b_face_vtx_idx[],
                          const cs_int_t    b_face_vtx_lst[],
                          cs_real_t        *uma,
                          cs_real_t        *vma,
                          cs_real_t        *wma,
                          cs_real_t        *coefau,
                          cs_real_t        *coefav,
                          cs_real_t        *coefaw,
                          cs_real_t        *coefbu,
                          cs_real_t        *coefbv,
                          cs_real_t        *coefbw,
                          cs_real_t        *dt,
                          cs_real_t        *disp_proj)
{
  cs_int_t  i, j, face_id, vtx_id, cell_id, cell_id1, cell_id2;

  cs_real_t  *vtx_counter = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_b_faces = mesh->n_b_faces;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  dim = mesh->dim;

  BFT_MALLOC(vtx_counter, n_vertices, cs_real_t);

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++) {

    vtx_counter[vtx_id] = 0.;
    for (i = 0; i < dim; i++)
      disp_proj[n_vertices*i + vtx_id] = 0.;

  }

  /* Internal face treatment */

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cell_id1 = i_face_cells[2*face_id] - 1;
    cell_id2 = i_face_cells[2*face_id+1] - 1;

    if (cell_id1 < n_cells) { /* Test to take into account face only once */

      for (j = i_face_vtx_idx[face_id]; j < i_face_vtx_idx[face_id+1]; j++) {

        /* Get the vertex number */

        vtx_id = i_face_vtx_lst[j-1] - 1;

        disp_proj[vtx_id] +=
          0.5 *(dt[cell_id1]*uma[cell_id1] + dt[cell_id2]*uma[cell_id2]);

        disp_proj[vtx_id+n_vertices] +=
          0.5 *(dt[cell_id1]*vma[cell_id1] + dt[cell_id2]*vma[cell_id2]);

        disp_proj[vtx_id+2*n_vertices] +=
          0.5 *(dt[cell_id1]*wma[cell_id1] + dt[cell_id2]*wma[cell_id2]);

        vtx_counter[vtx_id] += 1.;

      }

    }

  } /* End of loop on internal faces */

  /* Border face treatment.
     We reintialize vtx_counter on border faces in order to take into account
     only border face contribution */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    for (j = b_face_vtx_idx[face_id]; j < b_face_vtx_idx[face_id+1]; j++) {

      vtx_id = b_face_vtx_lst[j-1] - 1;
      vtx_counter[vtx_id] = 0.;

      for (i = 0; i < dim; i++) disp_proj[vtx_id+i*n_vertices]=0.;

    }

  } /* End of loop on border faces */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    cell_id = b_face_cells[face_id] - 1;

    for (j = b_face_vtx_idx[face_id]; j < b_face_vtx_idx[face_id+1]; j++) {

      vtx_id = b_face_vtx_lst[j-1] - 1;

      disp_proj[vtx_id] +=
        dt[cell_id]*(coefau[face_id] + coefbu[face_id]*uma[cell_id]);

      disp_proj[vtx_id + n_vertices] +=
        dt[cell_id]*(coefav[face_id] + coefbv[face_id]*vma[cell_id]);

      disp_proj[vtx_id + 2*n_vertices] +=
        dt[cell_id]*(coefaw[face_id] + coefbw[face_id]*wma[cell_id]);

      vtx_counter[vtx_id] += 1.;

    }

  } /* End of loop on border faces */

  if (mesh->vtx_interfaces != NULL) {
    cs_interface_set_sum(mesh->vtx_interfaces,
                         n_vertices,
                         3,
                         false,
                         CS_REAL_TYPE,
                         disp_proj);
    cs_interface_set_sum(mesh->vtx_interfaces,
                         n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         vtx_counter);
  }

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++)
    for (i = 0; i < dim; i++)
      disp_proj[vtx_id + i*n_vertices] /= vtx_counter[vtx_id];

  BFT_FREE(vtx_counter);
}

/*----------------------------------------------------------------------------
 * Projection on mesh vertices of the displacement (computed on cell center)
 *
 * Fortran Interface
 *
 * subroutine aledis
 * *****************
 *
 * ifacel            : <-- : Interior faces -> cells connectivity
 * ifabor            : <-- : Border faces -> cells connectivity
 * ipnfac            : <-- : Interior faces -> vertices index
 * nodfac            : <-- : Interior faces -> vertices list
 * ipnfbr            : <-- : Border faces -> vertices index
 * nodfbr            : <-- : Border faces -> vertices list
 * ialtyb            : <-- : Type of boundary for ALE
 * meshv             : <-- : Mesh velocity
 * gradm             : <-- : Mesh velocity gradient (du_i/dx_j : gradv[][i][j])
 * claale            : <-- : Boundary conditions A
 * clbale            : <-- : Boundary conditions B
 * dt                : <-- : Time step
 * deproj            : --> : Displacement projected on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (aledis, ALEDIS)(const cs_int_t      i_face_cells[],
                          const cs_int_t      b_face_cells[],
                          const cs_int_t      i_face_vtx_idx[],
                          const cs_int_t      i_face_vtx_lst[],
                          const cs_int_t      b_face_vtx_idx[],
                          const cs_int_t      b_face_vtx_lst[],
                          const cs_int_t      ialtyb[],
                          const cs_real_t    *meshv,
                          const cs_real_33_t  gradm[],
                          const cs_real_t    *claale,
                          const cs_real_t    *clbale,
                          const cs_real_t    *dt,
                          cs_real_t          *disp_proj)
{
  cs_int_t  i, j, face_id, vtx_id, cell_id, cell_id1, cell_id2;

  cs_real_t  *vtx_counter = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_b_faces = mesh->n_b_faces;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  dim = mesh->dim;
  const cs_real_t  *vtx_coord = mesh->vtx_coord;
  const cs_real_t  *cell_cen = cs_glob_mesh_quantities->cell_cen;
  const cs_real_t  *face_cen = cs_glob_mesh_quantities->b_face_cog;

  BFT_MALLOC(vtx_counter, n_vertices, cs_real_t);

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++) {

    vtx_counter[vtx_id] = 0.;

    for (i = 0; i < dim; i++)
      disp_proj[i + dim*vtx_id] = 0.;

  }

  /* Interior face treatment */

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cell_id1 = i_face_cells[2*face_id] - 1;
    cell_id2 = i_face_cells[2*face_id+1] - 1;

    cs_real_t dvol1 = 1./cs_glob_mesh_quantities->cell_vol[cell_id1];
    cs_real_t dvol2 = 1./cs_glob_mesh_quantities->cell_vol[cell_id2];

    if (cell_id1 < n_cells) { /* Test to take into account face only once */

      for (j = i_face_vtx_idx[face_id]; j < i_face_vtx_idx[face_id+1]; j++) {

        /* Get the vertex number */

        vtx_id = i_face_vtx_lst[j-1] - 1;

        /* Get the vector from the cell center to the node*/

        cs_real_t cen1_node_x = -cell_cen[3*cell_id1]   + vtx_coord[3*vtx_id];
        cs_real_t cen2_node_x = -cell_cen[3*cell_id2]   + vtx_coord[3*vtx_id];
        cs_real_t cen1_node_y = -cell_cen[3*cell_id1+1] + vtx_coord[3*vtx_id+1];
        cs_real_t cen2_node_y = -cell_cen[3*cell_id2+1] + vtx_coord[3*vtx_id+1];
        cs_real_t cen1_node_z = -cell_cen[3*cell_id1+2] + vtx_coord[3*vtx_id+2];
        cs_real_t cen2_node_z = -cell_cen[3*cell_id2+2] + vtx_coord[3*vtx_id+2];

        disp_proj[dim*vtx_id] +=
          dvol1*dt[cell_id1]*(meshv[3*cell_id1] + gradm[cell_id1][0][0]*cen1_node_x
                                                + gradm[cell_id1][0][1]*cen1_node_y
                                                + gradm[cell_id1][0][2]*cen1_node_z)
         +dvol2*dt[cell_id2]*(meshv[3*cell_id2] + gradm[cell_id2][0][0]*cen2_node_x
                                                + gradm[cell_id2][0][1]*cen2_node_y
                                                + gradm[cell_id2][0][2]*cen2_node_z);

        disp_proj[1 + dim*vtx_id] +=
          dvol1*dt[cell_id1]*(meshv[3*cell_id1+1] + gradm[cell_id1][1][0]*cen1_node_x
                                                  + gradm[cell_id1][1][1]*cen1_node_y
                                                  + gradm[cell_id1][1][2]*cen1_node_z)
         +dvol2*dt[cell_id2]*(meshv[3*cell_id2+1] + gradm[cell_id2][1][0]*cen2_node_x
                                                  + gradm[cell_id2][1][1]*cen2_node_y
                                                  + gradm[cell_id2][1][2]*cen2_node_z);

        disp_proj[2 + dim*vtx_id] +=
          dvol1*dt[cell_id1]*(meshv[3*cell_id1+2] + gradm[cell_id1][2][0]*cen1_node_x
                                                  + gradm[cell_id1][2][1]*cen1_node_y
                                                  + gradm[cell_id1][2][2]*cen1_node_z)
         +dvol2*dt[cell_id2]*(meshv[3*cell_id2+2] + gradm[cell_id2][2][0]*cen2_node_x
                                                  + gradm[cell_id2][2][1]*cen2_node_y
                                                  + gradm[cell_id2][2][2]*cen2_node_z);

        vtx_counter[vtx_id] += dvol1+dvol2;

      }

    }

  } /* End of loop on internal faces */

  /* Border face treatment.
     We reintialize vtx_counter on border faces in order to take into account
     only border face contribution */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    /* If the boundary face is NOT a sliding face */

    if (ialtyb[face_id] != 2) {

      for (j = b_face_vtx_idx[face_id]; j < b_face_vtx_idx[face_id+1]; j++) {

        vtx_id = b_face_vtx_lst[j-1] - 1;
        vtx_counter[vtx_id] = 0.;

        for (i = 0; i < dim; i++)
          disp_proj[i + dim*vtx_id]=0.;

      }
    }

  } /* End of loop on border faces */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    cell_id = b_face_cells[face_id] - 1;

    cs_real_t dsurf = 1./cs_glob_mesh_quantities->b_face_surf[face_id];

    for (j = b_face_vtx_idx[face_id]; j < b_face_vtx_idx[face_id+1]; j++) {

      vtx_id = b_face_vtx_lst[j-1] - 1;

      /* If the boundary face is NOT a sliding face */

      if (ialtyb[face_id] != 2) {

        /* Get the vector from the face center to the node*/

        cs_real_t face_node_x = -face_cen[3*face_id]   + vtx_coord[3*vtx_id];
        cs_real_t face_node_y = -face_cen[3*face_id+1] + vtx_coord[3*vtx_id+1];
        cs_real_t face_node_z = -face_cen[3*face_id+2] + vtx_coord[3*vtx_id+2];

        cs_real_t vel_cen_x = meshv[3*cell_id];
        cs_real_t vel_cen_y = meshv[3*cell_id+1];
        cs_real_t vel_cen_z = meshv[3*cell_id+2];

        /* 1st order extrapolation of the mesh velocity at the face center to the node */

        cs_real_t vel_node_x = claale[3*face_id  ] + gradm[cell_id][0][0]*face_node_x
                                                   + gradm[cell_id][0][1]*face_node_y
                                                   + gradm[cell_id][0][2]*face_node_z;
        cs_real_t vel_node_y = claale[3*face_id+1] + gradm[cell_id][1][0]*face_node_x
                                                   + gradm[cell_id][1][1]*face_node_y
                                                   + gradm[cell_id][1][2]*face_node_z;
        cs_real_t vel_node_z = claale[3*face_id+2] + gradm[cell_id][2][0]*face_node_x
                                                   + gradm[cell_id][2][1]*face_node_y
                                                   + gradm[cell_id][2][2]*face_node_z;

        disp_proj[dim*vtx_id] += dsurf*
          dt[cell_id]*(vel_node_x + clbale[9*face_id  ]*vel_cen_x
                                  + clbale[9*face_id+1]*vel_cen_y
                                  + clbale[9*face_id+2]*vel_cen_z);

        disp_proj[1 + dim*vtx_id] += dsurf*
          dt[cell_id]*(vel_node_y + clbale[9*face_id+3]*vel_cen_x
                                  + clbale[9*face_id+4]*vel_cen_y
                                  + clbale[9*face_id+5]*vel_cen_z);

        disp_proj[2 + dim*vtx_id] += dsurf*
          dt[cell_id]*(vel_node_z + clbale[9*face_id+6]*vel_cen_x
                                  + clbale[9*face_id+7]*vel_cen_y
                                  + clbale[9*face_id+8]*vel_cen_z);

        vtx_counter[vtx_id] += dsurf;

      }

      /* If the boundary face IS a sliding face.
         We project the displacment paralelly to the face. */

      else {

        cs_real_t tempox = clbale[9*face_id  ]*disp_proj[dim*vtx_id]
                         + clbale[9*face_id+1]*disp_proj[dim*vtx_id + 1]
                         + clbale[9*face_id+2]*disp_proj[dim*vtx_id + 2];
        cs_real_t tempoy = clbale[9*face_id+3]*disp_proj[dim*vtx_id]
                         + clbale[9*face_id+4]*disp_proj[dim*vtx_id + 1]
                         + clbale[9*face_id+5]*disp_proj[dim*vtx_id + 2];
        cs_real_t tempoz = clbale[9*face_id+6]*disp_proj[dim*vtx_id]
                         + clbale[9*face_id+7]*disp_proj[dim*vtx_id + 1]
                         + clbale[9*face_id+8]*disp_proj[dim*vtx_id + 2];

        disp_proj[  dim*vtx_id] = tempox;
        disp_proj[1+dim*vtx_id] = tempoy;
        disp_proj[2+dim*vtx_id] = tempoz;

      }
    }

  } /* End of loop on border faces */

  if (mesh->vtx_interfaces != NULL) {
    cs_interface_set_sum(mesh->vtx_interfaces,
                         n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         disp_proj);
    cs_interface_set_sum(mesh->vtx_interfaces,
                         n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         vtx_counter);
  }

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++)
    for (i = 0; i < dim; i++)
      disp_proj[i + dim*vtx_id] /= vtx_counter[vtx_id];

  BFT_FREE(vtx_counter);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
