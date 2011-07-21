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
 * Functions associated to ALE formulation
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_interface.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ale.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

static fvm_interface_set_t  *_ale_interface = NULL;

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

  const cs_int_t  n_vertices = cs_glob_mesh->n_vertices;
  const cs_int_t  n_cells = cs_glob_mesh->n_cells;
  const cs_int_t  n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_int_t  n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_int_t  dim = cs_glob_mesh->dim;

  if (cs_glob_mesh->global_vtx_num != NULL && _ale_interface == NULL)
    _ale_interface
      = fvm_interface_set_create(n_vertices,
                                 NULL,
                                 cs_glob_mesh->global_vtx_num,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL,
                                 NULL);

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

    if (cell_id1 <= n_cells) { /* Test to take into account face only once */

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

  if (cs_glob_mesh->global_vtx_num != NULL) {
    cs_parall_interface_sr(_ale_interface, n_vertices, 3, disp_proj);
    cs_parall_interface_sr(_ale_interface, n_vertices, 1, vtx_counter);
  }


  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++)
    for (i = 0; i < dim; i++)
      disp_proj[vtx_id + i*n_vertices] /= vtx_counter[vtx_id];

  BFT_FREE(vtx_counter);

}

/*----------------------------------------------------------------------------
 * Destroy if necessary the associated fvm_interface_set_t structure
 *
 * Fortran Interface
 *
 * SUBROUTINE LBRALE
 * *****************
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lbrale, LBRALE)(void)
{

  if (_ale_interface != NULL)
    _ale_interface = fvm_interface_set_destroy(_ale_interface);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
