/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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
 * Optional mesh renumbering
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
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

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_prototypes.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_renumber.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Redistribute vector values in case of renubering
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   n_elts          <->  Number of elements
 *   renum           <->  Pointer to face renulbering array (1 to n)
 *   val             <->  Pointer to array of vector values
 *   tmp_val         <->  Working array (size n_elts)
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_elt_vector(cs_int_t            n_elts,
                        const cs_int_t     *renum,
                        cs_real_t          *val,
                        cs_real_t          *tmp_val)
{
  int  dim_id;
  cs_int_t  face_id, tmp_face_id;

  for (dim_id = 0; dim_id < 3; dim_id++) {

    for (face_id = 0; face_id < n_elts; face_id++) {
      tmp_face_id = renum[face_id] - 1;
      tmp_val[face_id] = val[tmp_face_id*3 + dim_id];
    }

    for (face_id = 0; face_id < n_elts; face_id++)
      val[face_id*3 + dim_id] = tmp_val[face_id];

  }

}

/*----------------------------------------------------------------------------
 * Update quantities in case they were build before renumbering.
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_quantities(cs_mesh_t             *mesh,
                               cs_mesh_quantities_t  *mesh_quantities,
                               const cs_int_t        *renum_i,
                               const cs_int_t        *renum_b)
{
  cs_real_t  *tmp_val = NULL;
  cs_int_t  n_faces_max = CS_MAX(mesh->n_i_faces, mesh->n_b_faces);

  if (mesh == NULL && mesh_quantities == NULL)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(tmp_val, n_faces_max, cs_real_t);

  /* Interior faces */

  if (renum_i != NULL) {

    if (mesh_quantities->i_face_normal != NULL)
      _cs_renumber_elt_vector(mesh->n_i_faces,
                              renum_i,
                              mesh_quantities->i_face_normal,
                              tmp_val);

    if (mesh_quantities->i_face_cog != NULL)
      _cs_renumber_elt_vector(mesh->n_i_faces,
                              renum_i,
                              mesh_quantities->i_face_cog,
                              tmp_val);

  }

  /* Boundary Faces */

  if (renum_b != NULL) {

    if (mesh_quantities->b_face_normal != NULL)
      _cs_renumber_elt_vector(mesh->n_b_faces,
                              renum_b,
                              mesh_quantities->b_face_normal,
                              tmp_val);

    if (mesh_quantities->b_face_cog != NULL)
      _cs_renumber_elt_vector(mesh->n_b_faces,
                              renum_b,
                              mesh_quantities->b_face_cog,
                              tmp_val);

  }

  /* Free Work arrays */

  BFT_FREE(tmp_val);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Try to apply renumbering of faces for vector machines.
 *
 * Renumbering can be cancelled using the IVECTI and IVECTB values in
 * Fortan common IVECTO: -1 indicates we should try to renumber,
 * 0 means we should not renumber. On exit, 0 means we have not found an
 * adequate renumbering, 1 means we have (and it was applied).
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_for_vectorizing(cs_mesh_t             *mesh,
                            cs_mesh_quantities_t  *mesh_quantities)
{
  cs_int_t   ivecti = 0, ivectb = 0;
  cs_int_t  *inumfi = NULL, *inumfb = NULL;
  cs_int_t  *ipnfaw = NULL, *ipnfbw = NULL;
  cs_int_t  *nodfaw = NULL, *nodfbw = NULL;
  cs_int_t  *iworkf = NULL, *ismbs = NULL, *ismbv = NULL;
  cs_real_t  *rworkf = NULL, *rsmbs = NULL, *rsmbv = NULL;

  cs_int_t  n_faces_max = CS_MAX(mesh->n_i_faces, mesh->n_b_faces);
  cs_int_t  n_cells_wghosts = mesh->n_cells_with_ghosts;

  /* Allocate Work arrays */

  BFT_MALLOC(inumfi, mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(inumfb, mesh->n_b_faces, cs_int_t);
  BFT_MALLOC(iworkf, n_faces_max, cs_int_t);
  BFT_MALLOC(ismbs, n_cells_wghosts, cs_int_t);
  BFT_MALLOC(ismbv, n_cells_wghosts, cs_int_t);
  BFT_MALLOC(ipnfaw, mesh->n_i_faces + 1, cs_int_t);
  BFT_MALLOC(nodfaw, mesh->i_face_vtx_connect_size, cs_int_t);
  BFT_MALLOC(ipnfbw, mesh->n_b_faces + 1, cs_int_t);
  BFT_MALLOC(nodfbw, mesh->b_face_vtx_connect_size, cs_int_t);
  BFT_MALLOC(rworkf, n_faces_max, cs_real_t);
  BFT_MALLOC(rsmbs, n_cells_wghosts, cs_real_t);
  BFT_MALLOC(rsmbv, n_cells_wghosts, cs_real_t);

  /* Try renumbering */

  CS_PROCF(numvec, NUMVEC)(&(mesh->n_cells_with_ghosts),
                           &(mesh->n_cells),
                           &(mesh->n_i_faces),
                           &(mesh->n_b_faces),
                           &(mesh->n_vertices),
                           &(mesh->i_face_vtx_connect_size),
                           &(mesh->b_face_vtx_connect_size),
                           &ivecti,
                           &ivectb,
                           mesh->i_face_cells,
                           mesh->b_face_cells,
                           mesh->i_face_family,
                           mesh->b_face_family,
                           mesh->i_face_vtx_idx,
                           mesh->i_face_vtx_lst,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           inumfi,
                           inumfb,
                           iworkf,
                           ismbs,
                           ismbv,
                           ipnfaw,
                           nodfaw,
                           ipnfbw,
                           nodfbw,
                           rworkf,
                           rsmbs,
                           rsmbv);

  /* Free Work arrays */

  BFT_FREE(rsmbv);
  BFT_FREE(rsmbs);
  BFT_FREE(rworkf);
  BFT_FREE(nodfbw);
  BFT_FREE(ipnfbw);
  BFT_FREE(nodfaw);
  BFT_FREE(ipnfaw);
  BFT_FREE(ismbv);
  BFT_FREE(ismbs);
  BFT_FREE(iworkf);

  /* Renumber mesh quantities if already computed */

  if (ivecti > 0 || ivectb > 0) {

    cs_int_t   *_inumfi = NULL;
    cs_int_t   *_inumfb = NULL;

    if (ivecti > 0)
      _inumfi = inumfi;
    if (ivectb > 0)
      _inumfb = inumfb;

    _cs_renumber_update_quantities(mesh,
                                   mesh_quantities,
                                   _inumfi,
                                   _inumfb);

  }

  /* Free final work arrays */

  BFT_FREE(inumfb);
  BFT_FREE(inumfi);

}

/*----------------------------------------------------------------------------
 * Renumber mesh elements depending on code options and target machine.
 *
 * Currently, only the legacy vectorizing renumbering is handled.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_mesh(cs_mesh_t             *mesh,
                 cs_mesh_quantities_t  *mesh_quantities)
{
  cs_renumber_for_vectorizing(mesh,
                              mesh_quantities);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
