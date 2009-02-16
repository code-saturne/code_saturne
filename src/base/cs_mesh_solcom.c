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
 * Read a mesh in "SolCom" format
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_nodal_append.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_post.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_solcom.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of vertices, tetrahedra, pyramids, prisms, and hexahedra */

static cs_int_t  cs_glob_nsom = 0;
static cs_int_t  cs_glob_ntetra = 0;
static cs_int_t  cs_glob_npyram = 0;
static cs_int_t  cs_glob_nprism = 0;
static cs_int_t  cs_glob_nhexae = 0;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate memory for a mesh in "SolCom" format.
 *
 * parameters:
 *   mesh            <-- associated mesh
 *   mesh_quantities <-- associated quantities
 *----------------------------------------------------------------------------*/

static void
_mesh_solcom_alloc_mem(cs_mesh_t             *mesh,
                       cs_mesh_quantities_t  *mesh_quantities)
{
  cs_int_t  n_elts;

  /* Allocation */
  /*------------*/

  /* Faces / cells connectivity*/
  BFT_MALLOC(mesh->i_face_cells, mesh->n_i_faces * 2, cs_int_t);
  BFT_MALLOC(mesh->b_face_cells, mesh->n_b_faces, cs_int_t);

  /* Cell centers (sized for true + ghost cells, but there should be no
     ghost cells here) */
  n_elts = mesh->dim * mesh->n_cells_with_ghosts;
  BFT_MALLOC(mesh_quantities->cell_cen, n_elts, cs_real_t);

  /* Face surfaces */
  BFT_MALLOC(mesh_quantities->i_face_normal, mesh->dim * mesh->n_i_faces,
            cs_real_t);
  BFT_MALLOC(mesh_quantities->b_face_normal, mesh->dim * mesh->n_b_faces,
            cs_real_t);

  /* Face centers */
  BFT_MALLOC(mesh_quantities->i_face_cog, mesh->dim * mesh->n_i_faces,
            cs_real_t);
  BFT_MALLOC(mesh_quantities->b_face_cog, mesh->dim * mesh->n_b_faces,
            cs_real_t);

  /* Cell and boundary face families */
  BFT_MALLOC(mesh->b_face_family, mesh->n_b_faces, cs_int_t);
  BFT_MALLOC(mesh->cell_family, mesh->n_cells_with_ghosts, cs_int_t);

  /* Family properties */
  n_elts = mesh->n_families * mesh->n_max_family_items;
  BFT_MALLOC(mesh->family_item, n_elts, cs_int_t);


  if (mesh->n_vertices > 0) {

    /* Vertex coordinates */
    n_elts = mesh->dim * mesh->n_vertices;
    BFT_MALLOC(mesh->vtx_coord, n_elts, cs_real_t);

    /* Faces / vertices connectivity */
    BFT_MALLOC(mesh->i_face_vtx_idx, mesh->n_i_faces + 1, cs_int_t);
    BFT_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_int_t);
    BFT_MALLOC(mesh->b_face_vtx_idx, mesh->n_b_faces + 1, cs_int_t);
    BFT_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_int_t);

  }

}

/*----------------------------------------------------------------------------
 * Transfer a part of the nodal connectivity to an FVM structure.
 *
 * parameters:
 *   ext_mesh      <-> mesh to complete
 *   n_elts        <-- number of elements to add
 *   type          <-- type of section to add
 *   connect       <-- connectivity to transfer
 *   tot_elt_count <-> total element counter
 *----------------------------------------------------------------------------*/

static void
_mesh_solcom_add(fvm_nodal_t    *ext_mesh,
                 cs_int_t        n_elts,
                 fvm_element_t   type,
                 cs_int_t       *connect,
                 cs_int_t       *tot_elt_count)
{
  cs_int_t   ind_elt = 0;
  cs_int_t   cpt_elt = 0;

  fvm_lnum_t   *parent_element_num = NULL;

  if (n_elts == 0)
    return;

  /* Parent element numbers if necessary */

  if (cpt_elt > 0) {
    BFT_MALLOC(parent_element_num, n_elts, fvm_lnum_t);
    for (ind_elt = 0 ; ind_elt < n_elts ; ind_elt++)
      parent_element_num[ind_elt] = ind_elt + (*tot_elt_count + 1);
  }

  /* Transfer connectivity and parent numbers */

  fvm_nodal_append_by_transfer(ext_mesh,
                               n_elts,
                               type,
                               NULL,
                               NULL,
                               NULL,
                               connect,
                               parent_element_num);

  *tot_elt_count = *tot_elt_count + n_elts;
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update mesh size information after reading a "SolCom" format file header.
 *
 * Fortran interface:
 *
 * SUBROUTINE DIMGEO (NDIM  , NCELET, NCEL  , NFAC  , NFABOR, NSOM  ,
 * *****************
 *                    LNDFAC, LNDFBR, NFML  , NPRFML,
 *                    NTETRA, NPYRAM, NPRISM, NHEXAE )
 *
 * INTEGER          NDIM        : <-- : spatial dimension (3)
 * INTEGER          NCELET      : <-- : number of extended cells
 * INTEGER          NCEL        : <-- : number of true cells
 * INTEGER          NFAC        : <-- : number of interior faces
 * INTEGER          NFABOR      : <-- : number of boundary faces
 * INTEGER          NSOM        : <-- : number of vertices (optional)
 * INTEGER          LNDFAC      : <-- : length of SOMFAC (optional)
 * INTEGER          LNDFBR      : <-- : length of SOMFBR (optional)
 * INTEGER          NFML        : <-- : number of families
 * INTEGER          NPRFML      : <-- : max. number of properties per family
 * INTEGER          NTETRA      : <-- : number of tetrahedra
 * INTEGER          NPYRAM      : <-- : number of pyramids
 * INTEGER          NPRISM      : <-- : number of prisms
 * INTEGER          NHEXAE      : <-- : number of hexahedra
 *----------------------------------------------------------------------------*/

void CS_PROCF (dimgeo, DIMGEO)
(
 const cs_int_t   *ndim,
 const cs_int_t   *ncelet,
 const cs_int_t   *ncel,
 const cs_int_t   *nfac,
 const cs_int_t   *nfabor,
 const cs_int_t   *nsom,
 const cs_int_t   *lndfac,
 const cs_int_t   *lndfbr,
 const cs_int_t   *nfml,
 const cs_int_t   *nprfml,
 const cs_int_t   *ntetra,
 const cs_int_t   *npyram,
 const cs_int_t   *nprism,
 const cs_int_t   *nhexae
)
{
  cs_mesh_t *mesh = cs_glob_mesh;

  mesh->dim = *ndim;

  mesh->n_cells = *ncel;
  mesh->n_i_faces = *nfac;
  mesh->n_b_faces = *nfabor;

  cs_glob_nsom = *nsom;

  if (*lndfac + *lndfbr > 0)
    mesh->n_vertices = *nsom;
  else
    mesh->n_vertices = 0;

  mesh->i_face_vtx_connect_size = *lndfac;
  mesh->b_face_vtx_connect_size = *lndfbr;

  mesh->n_cells_with_ghosts = *ncelet;

  assert (*ncelet == *ncel);

  mesh->n_g_cells = (fvm_gnum_t)mesh->n_cells;
  mesh->n_g_i_faces = (fvm_gnum_t)mesh->n_i_faces;
  mesh->n_g_b_faces = (fvm_gnum_t)mesh->n_b_faces;
  mesh->n_g_vertices = (fvm_gnum_t)mesh->n_vertices;

  mesh->n_max_family_items = *nprfml;
  mesh->n_families          = *nfml;

  cs_glob_ntetra = *ntetra;
  cs_glob_npyram = *npyram;
  cs_glob_nprism = *nprism;
  cs_glob_nhexae = *nhexae;

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read a mesh in "SolCom" format (prior to Code_Saturne 1.0)
 *
 * parameters:
 *   mesh            <-- associated mesh
 *   mesh_quantities <-- associated quantities
 *----------------------------------------------------------------------------*/

void
cs_mesh_solcom_read(cs_mesh_t             *mesh,
                    cs_mesh_quantities_t  *mesh_quantities)
{
  cs_int_t   indic_nodal = 0;
  cs_int_t   tot_elt_count = 0;

  cs_real_t  *vtx_coord = NULL;
  cs_int_t   *connect_tetra = NULL;
  cs_int_t   *connect_pyram = NULL;
  cs_int_t   *connect_prism = NULL;
  cs_int_t   *connect_hexae = NULL;

  fvm_nodal_t  *ext_mesh = NULL;

  /* Allocations for main mesh */

  _mesh_solcom_alloc_mem(mesh, mesh_quantities);

  /* Allocations for post-processing mesh when we do not have
     a faces -> vertices connectivity */

  if (mesh->vtx_coord != NULL)
    vtx_coord = mesh->vtx_coord;

  else {

    BFT_MALLOC(vtx_coord, cs_glob_nsom * 3, cs_real_t);
    BFT_MALLOC(connect_tetra, cs_glob_ntetra * 4, cs_int_t);
    BFT_MALLOC(connect_pyram, cs_glob_npyram * 5, cs_int_t);
    BFT_MALLOC(connect_prism, cs_glob_nprism * 6, cs_int_t);
    BFT_MALLOC(connect_hexae, cs_glob_nhexae * 8, cs_int_t);

  }

  /* Read mesh body */

  CS_PROCF (letgeo, LETGEO) (&(mesh->dim),
                             &(mesh->n_cells_with_ghosts),
                             &(mesh->n_cells),
                             &(mesh->n_i_faces),
                             &(mesh->n_b_faces),
                             &(mesh->n_families),
                             &(mesh->n_max_family_items),
                             &(cs_glob_nsom),
                             &(mesh->i_face_vtx_connect_size),
                             &(mesh->b_face_vtx_connect_size),
                             &cs_glob_ntetra,
                             &cs_glob_npyram,
                             &cs_glob_nprism,
                             &cs_glob_nhexae,
                             &indic_nodal,
                             mesh->i_face_cells,
                             mesh->b_face_cells,
                             mesh->b_face_family,
                             mesh->cell_family,
                             mesh->family_item,
                             connect_tetra,
                             connect_pyram,
                             connect_prism,
                             connect_hexae,
                             mesh->i_face_vtx_idx,
                             mesh->i_face_vtx_lst,
                             mesh->b_face_vtx_idx,
                             mesh->b_face_vtx_lst,
                             mesh_quantities->cell_cen,
                             mesh_quantities->i_face_normal,
                             mesh_quantities->b_face_normal,
                             mesh_quantities->i_face_cog,
                             mesh_quantities->b_face_cog,
                             vtx_coord);


  if (indic_nodal > 0) {

    /* Direct creation of the  post-processing mesh when we
       do not have a faces -> vertices connectivity */

    ext_mesh = fvm_nodal_create(_("Fluid volume"), 3);

    if (cs_glob_ntetra > 0)
      _mesh_solcom_add(ext_mesh,
                       cs_glob_ntetra,
                       FVM_CELL_TETRA,
                       connect_tetra,
                       &tot_elt_count);

    if (cs_glob_npyram > 0)
      _mesh_solcom_add(ext_mesh,
                       cs_glob_npyram,
                       FVM_CELL_PYRAM,
                       connect_pyram,
                       &tot_elt_count);

    if (cs_glob_nprism > 0)
      _mesh_solcom_add(ext_mesh,
                       cs_glob_nprism,
                       FVM_CELL_PRISM,
                       connect_prism,
                       &tot_elt_count);

    if (cs_glob_nhexae > 0)
      _mesh_solcom_add(ext_mesh,
                       cs_glob_nhexae,
                       FVM_CELL_HEXA,
                       connect_hexae,
                       &tot_elt_count);

    fvm_nodal_transfer_vertices(ext_mesh, vtx_coord);

    /* Transfer structure to post-processing */

    cs_post_add_existing_mesh(-1, ext_mesh, true);

  }
  else if (mesh->vtx_coord == NULL) {

    BFT_FREE(vtx_coord);
    BFT_FREE(connect_tetra);
    BFT_FREE(connect_pyram);
    BFT_FREE(connect_prism);
    BFT_FREE(connect_hexae);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
