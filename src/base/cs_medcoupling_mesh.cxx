/*============================================================================
 * Usage of MEDCoupling base components.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_timer.h"

#include "fvm_defs.h"
#include "fvm_nodal_from_desc.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#include "cs_medcoupling_mesh.hxx"

#if defined(HAVE_MEDCOUPLING)
#include <MEDCoupling_version.h>

// Medloader
#if defined(HAVE_MEDCOUPLING_LOADER)
#include <MEDFileMesh.hxx>
#endif

using namespace MEDCoupling;

#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const cs_lnum_t _perm_tri[3]  = {0, 2, 1};
static const cs_lnum_t _perm_quad[4] = {0, 3, 2, 1};
static const cs_lnum_t _perm_pent[5] = {0, 4, 3, 2, 1};

static int                     _n_sub_meshes = 0;
static cs_medcoupling_mesh_t **_sub_meshes = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get permutation array for a face given its number of vertices.
 *
 * \param[in] n_face_vertices number of vertices of the face
 *
 * \return    pointer to the permutation array
 */
/*----------------------------------------------------------------------------*/

static inline const cs_lnum_t *
_get_face_vertices_permutation(cs_lnum_t  n_face_vertices)
{
  const cs_lnum_t *perm = NULL;

  switch (n_face_vertices) {
    case 3:
      perm = _perm_tri;
      break;
    case 4:
      perm = _perm_quad;
      break;
    case 5:
      perm = _perm_pent;
      break;
    default:
      perm = NULL;
      break;
  }

  return perm;
}

#if defined(HAVE_MEDCOUPLING)

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assign vertex coordinates to a MEDCoupling mesh structure
 *
 * \param[in] mesh      pointer to cs_mesh_t structure from which data is copied
 * \param[in] n_vtx     number of vertices to assign
 * \param[in] vtx_id    pointer to vertices id's used for assigning
 * \param[in] med_mesh  pointer to MEDCouplingUMesh to which we copy the
 *                      coordinates
 */
/*----------------------------------------------------------------------------*/

static void
_assign_vertex_coords(const cs_mesh_t    *mesh,
                      cs_lnum_t           n_vtx,
                      const cs_lnum_t    *vtx_id,
                      MEDCouplingUMesh   *med_mesh)
{
  const cs_lnum_t  dim = mesh->dim;
  const cs_coord_t  *vertex_coords = mesh->vtx_coord;

  assert(med_mesh != NULL);

  /* Assign all coordinates */
  /*------------------------*/

  DataArrayDouble *med_coords = DataArrayDouble::New();
  med_coords->alloc(n_vtx, dim);

  if (vtx_id != NULL) {
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
      if (vtx_id[i] > -1) {
        for (cs_lnum_t j = 0; j < dim; j++) {
          med_coords->setIJ(vtx_id[i], j, vertex_coords[i*dim + j]);
        }
      }
    }
  }
  else {
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
      for (cs_lnum_t j = 0; j < dim; j++)
        med_coords->setIJ(i, j, vertex_coords[i*dim + j]);
    }
  }

  med_mesh->setCoords(med_coords);
  med_coords->decrRef();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assign boundary faces to a MEDCoupling mesh structure
 *
 * \param[in] mesh      pointer to cs_mesh_t structure from which data is copie
 * \param[in] n_elts    number of faces to copy
 * \param[in] elts_list list of faces to copy
 * \param[in] med_mesh  pointer to MEDCouplingUMesh to which we copy the faces
 */
/*----------------------------------------------------------------------------*/

static void
_assign_face_mesh(const cs_mesh_t   *mesh,
                  cs_lnum_t          n_elts,
                  const cs_lnum_t   *elts_list,
                  MEDCouplingUMesh  *med_mesh,
                  cs_lnum_t          new_to_old[])
{
  INTERP_KERNEL::NormalizedCellType type;

  cs_lnum_t vtx_count = -1;
  cs_lnum_t elt_buf_size = 4;
  cs_lnum_t *vtx_id = NULL;

  /* Build old->new face id indirection */

  cs_lnum_t *face_id = NULL;
  cs_lnum_t face_count = 0;
  BFT_MALLOC(face_id, mesh->n_b_faces, cs_lnum_t);
  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++)
    face_id[i] = -1;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    face_id[elts_list[i]] = face_count++;
  }

  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    new_to_old[face_id[elts_list[ii]]] = elts_list[ii];
  }
  BFT_FREE(face_id);

  /* Mark and renumber vertices */

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);

  /* Initialize the value of vtx_id */
  for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
    vtx_id[i] = -1;

  /* Case with filter list */

  if (elts_list != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t eid = elts_list[i];
      for (cs_lnum_t j = mesh->b_face_vtx_idx[eid];
           j < mesh->b_face_vtx_idx[eid+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t j = mesh->b_face_vtx_idx[i];
           j < mesh->b_face_vtx_idx[i+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }

  }

  /* Assign coordinates */

  _assign_vertex_coords(mesh, vtx_count, vtx_id, med_mesh);

  /* Assign faces */

  mcIdType *elt_buf = NULL;
  BFT_MALLOC(elt_buf, elt_buf_size, mcIdType);
  med_mesh->allocateCells(n_elts);

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    cs_lnum_t eid = (elts_list != NULL) ? elts_list[i] : i;

    assert(eid >= 0 && eid < mesh->n_b_faces);

    mcIdType n_vtx = mesh->b_face_vtx_idx[eid+1] - mesh->b_face_vtx_idx[eid];

    cs_lnum_t connect_start = mesh->b_face_vtx_idx[eid];

    if (n_vtx > elt_buf_size) { /* reallocate buffer if required */
      elt_buf_size *= 2;
      BFT_REALLOC(elt_buf, elt_buf_size, mcIdType);
    }

    const cs_lnum_t *_perm_face = _get_face_vertices_permutation(n_vtx);
    if (_perm_face != NULL) {
      for (cs_lnum_t j = 0; j < n_vtx; j++)
        elt_buf[j] = vtx_id[mesh->b_face_vtx_lst[connect_start + _perm_face[j]]];
    }
    else {
      for (cs_lnum_t j = 0; j < n_vtx; j++)
        elt_buf[j] = vtx_id[mesh->b_face_vtx_lst[connect_start + n_vtx - 1 - j]];
    }
    switch(n_vtx) {
    case 3:
      type = INTERP_KERNEL::NORM_TRI3;
      break;
    case 4:
      type = INTERP_KERNEL::NORM_QUAD4;
      break;
    default:
      type = INTERP_KERNEL::NORM_POLYGON;
      break;
    }

    med_mesh->insertNextCell(type, n_vtx, elt_buf);

  }

  med_mesh->finishInsertingCells();

  BFT_FREE(elt_buf);
  BFT_FREE(vtx_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assign cells to a MEDCoupling mesh structure
 *
 * \param[in] mesh        pointer to cs_mesh_t structure from which data is copied
 * \param[in] n_elts      number of cells to assign
 * \param[in] elts_list   list of cells to assign
 * \param[in] med_mesh    pointer to MEDCouplingUMesh to which we copy the cells
 * \param[in] new_to_old  indirection array between local mesh connectivity
 *                        and MEDCouplingUMesh connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_assign_cell_mesh(const cs_mesh_t   *mesh,
                  cs_lnum_t          n_elts,
                  const cs_lnum_t    elts_list[],
                  MEDCouplingUMesh  *med_mesh,
                  cs_lnum_t          new_to_old[])
{
  INTERP_KERNEL::NormalizedCellType type;

  cs_lnum_t vtx_count = 0, cell_count = 0;

  cs_lnum_t  elt_buf_size = 8;
  cs_lnum_t *vtx_id = NULL;
  cs_lnum_t *cell_id = NULL;
  cs_lnum_t *cell_faces_idx = NULL, *cell_faces_num = NULL;

  /* Build old->new cell id indirection */

  BFT_MALLOC(cell_id, mesh->n_cells, cs_lnum_t);
  for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
    cell_id[i] = -1;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cell_id[elts_list[i]] = cell_count++;
  }

  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    new_to_old[cell_id[elts_list[ii]]] = elts_list[ii];
  }

  /* Mark and renumber vertices */

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);
  for (cs_lnum_t vid = 0; vid < mesh->n_vertices; vid++) {
    vtx_id[vid] = -1;
  }

  for (cs_lnum_t face_id = 0; face_id < mesh->n_b_faces; face_id++) {
    cs_lnum_t c_id = cell_id[mesh->b_face_cells[face_id]];
    if (c_id > -1) {
      for (cs_lnum_t j = mesh->b_face_vtx_idx[face_id];
           j < mesh->b_face_vtx_idx[face_id+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }
  }

  for (cs_lnum_t face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    cs_lnum_t c_id1 = mesh->i_face_cells[face_id][0];
    cs_lnum_t c_id2 = mesh->i_face_cells[face_id][1];
    c_id1 = (c_id1 < mesh->n_cells) ? cell_id[c_id1] : -1;
    c_id2 = (c_id2 < mesh->n_cells) ? cell_id[c_id2] : -1;
    if (c_id1 > -1 || c_id2 > -1) {
      for (cs_lnum_t j = mesh->i_face_vtx_idx[face_id];
           j < mesh->i_face_vtx_idx[face_id+1];
           j++) {
        cs_lnum_t vid = mesh->i_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }
  }

  /* Assign coordinates */

  _assign_vertex_coords(mesh, vtx_count, vtx_id, med_mesh);

  /* Build temporary descending connectivity */

  cs_mesh_connect_get_cell_faces(mesh,
                                 mesh->n_cells, /* TODO test using n_elts */
                                 cell_id,
                                 &cell_faces_idx,
                                 &cell_faces_num);

  BFT_FREE(cell_id);

  /* Now loop on cells */

  const cs_lnum_t  face_num_shift[2] = {0, mesh->n_b_faces};

  const cs_lnum_t  *face_vertices_idx[2] = {mesh->b_face_vtx_idx,
                                            mesh->i_face_vtx_idx};
  const cs_lnum_t  *face_vertices_num[2] = {mesh->b_face_vtx_lst,
                                            mesh->i_face_vtx_lst};

  mcIdType *elt_buf = NULL;
  BFT_MALLOC(elt_buf, elt_buf_size, mcIdType);
  for (cs_lnum_t  ii = 0; ii < elt_buf_size; ii++)
    elt_buf[ii] = -1;

  /* Allocate the cells array */
  med_mesh->allocateCells(n_elts);

  for (cs_lnum_t ic = 0; ic < n_elts; ic++) {

    mcIdType  n_vtx;
    cs_lnum_t cell_vtx[8];

    cs_lnum_t i = ic;

    fvm_element_t fvm_type = fvm_nodal_from_desc_cell(i,
                                                      2,
                                                      face_num_shift,
                                                      face_vertices_idx,
                                                      face_vertices_num,
                                                      cell_faces_idx,
                                                      cell_faces_num,
                                                      cell_vtx);

    switch(fvm_type) {

    case FVM_CELL_TETRA:
      type = INTERP_KERNEL::NORM_TETRA4;
      n_vtx = 4;
      elt_buf[0] = vtx_id[cell_vtx[0]-1];
      elt_buf[1] = vtx_id[cell_vtx[2]-1];
      elt_buf[2] = vtx_id[cell_vtx[1]-1];
      elt_buf[3] = vtx_id[cell_vtx[3]-1];
      break;

    case FVM_CELL_PYRAM:
      type = INTERP_KERNEL::NORM_PYRA5;
      n_vtx = 5;
      elt_buf[0] = vtx_id[cell_vtx[0]-1];
      elt_buf[1] = vtx_id[cell_vtx[3]-1];
      elt_buf[2] = vtx_id[cell_vtx[2]-1];
      elt_buf[3] = vtx_id[cell_vtx[1]-1];
      elt_buf[4] = vtx_id[cell_vtx[4]-1];
      break;

    case FVM_CELL_PRISM:
      type = INTERP_KERNEL::NORM_PENTA6;
      n_vtx = 6;
      elt_buf[0] = vtx_id[cell_vtx[0]-1];
      elt_buf[1] = vtx_id[cell_vtx[2]-1];
      elt_buf[2] = vtx_id[cell_vtx[1]-1];
      elt_buf[3] = vtx_id[cell_vtx[3]-1];
      elt_buf[4] = vtx_id[cell_vtx[5]-1];
      elt_buf[5] = vtx_id[cell_vtx[4]-1];
      break;

    case FVM_CELL_HEXA:
      type = INTERP_KERNEL::NORM_HEXA8;
      n_vtx = 8;
      elt_buf[0] = vtx_id[cell_vtx[0]-1];
      elt_buf[1] = vtx_id[cell_vtx[3]-1];
      elt_buf[2] = vtx_id[cell_vtx[2]-1];
      elt_buf[3] = vtx_id[cell_vtx[1]-1];
      elt_buf[4] = vtx_id[cell_vtx[4]-1];
      elt_buf[5] = vtx_id[cell_vtx[7]-1];
      elt_buf[6] = vtx_id[cell_vtx[6]-1];
      elt_buf[7] = vtx_id[cell_vtx[5]-1];
      break;

    default:
      type = INTERP_KERNEL::NORM_POLYHED;

      n_vtx = 0;

      cs_lnum_t s_id = cell_faces_idx[i] - 1;
      cs_lnum_t e_id = cell_faces_idx[i+1] -1;

      for (cs_lnum_t j = s_id; j < e_id; j++) {
        int face_sgn = 0;
        cs_lnum_t face_id;
        if (cell_faces_num[j] > 0) {
          face_id  = cell_faces_num[j] - 1;
          face_sgn = 1;
        }
        else {
          face_id  = -cell_faces_num[j] - 1;
          face_sgn = -1;
        }

        int fl = 1;
        if (face_id < face_num_shift[fl])
          fl = 0;
        face_id -= face_num_shift[fl];

        cs_lnum_t v_id_start = face_vertices_idx[fl][face_id];
        cs_lnum_t v_id_end   = face_vertices_idx[fl][face_id + 1];
        cs_lnum_t n_face_vertices  = v_id_end - v_id_start;

        while (n_vtx + n_face_vertices + 1 > elt_buf_size) {
          elt_buf_size *= 2;
          BFT_REALLOC(elt_buf, elt_buf_size, mcIdType);
        }

        /* Add separator after first face */
        if (j > s_id)
          elt_buf[n_vtx++] = -1;

        const cs_lnum_t *_face_perm = NULL;
        _get_face_vertices_permutation(n_face_vertices);
        if (_face_perm != NULL) {
          for (cs_lnum_t  ik = 0; ik < n_face_vertices; ik++) {
            cs_lnum_t iik = _face_perm[ik];
            cs_lnum_t l =   v_id_start
                          + (  n_face_vertices
                             + (iik*face_sgn))%n_face_vertices;
            cs_lnum_t vid = face_vertices_num[fl][l];
            elt_buf[n_vtx++] = vtx_id[vid];
          }
        }
        else {
          for (cs_lnum_t ik = 0; ik < n_face_vertices; ik++) {
            cs_lnum_t l =   v_id_start
                          + (  n_face_vertices
                             + (ik*face_sgn))%n_face_vertices;
            cs_lnum_t vid = face_vertices_num[fl][l];
            elt_buf[n_vtx++] = vtx_id[vid];
          }
        }
      } /* Loop on j (cell faces) */
    } /* switch on cell_type */

    med_mesh->insertNextCell(type, n_vtx, elt_buf);
  } /* Loop on cells */

  med_mesh->finishInsertingCells();

  BFT_FREE(elt_buf);
  BFT_FREE(cell_faces_num);;
  BFT_FREE(cell_faces_idx);
  BFT_FREE(vtx_id);
}

#endif /* HAVE_MEDCOUPLING - BEGINNING OF PRIVATE FUNCTIONS */


/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to cs_medcoupling_mesh_t based on selection criteria
 * and element's dimension.
 *
 * \param[in] selection_criteria  selection criteria (entire mesh or part of it)
 * \param[in] elt_dim             dimension of elements. 2: faces, 3: cells
 *
 * \return pointer to found mesh instance, NULL if none found.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
_get_mesh_from_parameters(const char *selection_criteria,
                          const int   elt_dim)
{

  cs_medcoupling_mesh_t *m = NULL;

  for (int i = 0; i < _n_sub_meshes; i++) {
    cs_medcoupling_mesh_t *mt = _sub_meshes[i];
    if (elt_dim == mt->elt_dim &&
        strcmp(mt->sel_criteria, selection_criteria) == 0) {
      m = mt;
      break;
    }
  }

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a cs_medcoupling_mesh_t to the list
 *
 * \param[in] name                name of the mesh
 * \param[in] selection_criteria  selection criteria (entire mesh or part of it)
 * \param[in] elt_dim             dimension of elements. 2: faces, 3: cells
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
_add_medcoupling_mesh(const char *name,
                      const char *selection_criteria,
                      int         elt_dim)
{

  cs_medcoupling_mesh_t *m = NULL;

#if !defined(HAVE_MEDCOUPLING)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else

  BFT_REALLOC(_sub_meshes, _n_sub_meshes + 1, cs_medcoupling_mesh_t *);

  BFT_MALLOC(m, 1, cs_medcoupling_mesh_t);

  size_t len_sel_crit = strlen(selection_criteria);
  BFT_MALLOC(m->sel_criteria, len_sel_crit+1, char);
  strcpy(m->sel_criteria, selection_criteria);
  m->sel_criteria[len_sel_crit] = '\0';

  m->elt_dim  = elt_dim;
  m->n_elts   = 0;
  m->elt_list = NULL;

  m->med_mesh = MEDCouplingUMesh::New();
  m->med_mesh->setName(name);
  m->med_mesh->setTimeUnit("s");
  m->med_mesh->setMeshDimension(elt_dim);

  m->bbox = NULL;

  m->new_to_old = NULL;

  _sub_meshes[_n_sub_meshes] = m;

  _n_sub_meshes += 1;

#endif

  return m;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief copy a cs_mesh_t into a cs_medcoupling_mesh_t
 *
 * \param[in] csmesh    pointer to the cs_mesh_t struct to copy data from
 * \param[in] pmmesh    pointer to the cs_medcoupling_mesh_t for copy
 * \param[in] use_bbox  flag indicating if a reduced bounding is used. Usefull
 *                      for interpolation to reduce the matrix sice.
 *                      0: Do not use a reduced bbox
 *                      1: Use a reduced bbox
 */
/*----------------------------------------------------------------------------*/

void
_copy_mesh_from_base(cs_mesh_t              *csmesh,
                     cs_medcoupling_mesh_t  *pmmesh,
                     int                     use_bbox)
{
#if !defined(HAVE_MEDCOUPLING)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else
  if (pmmesh->elt_dim == 3) {

    /* Creation of a new nodal mesh from selected cells */

    BFT_MALLOC(pmmesh->elt_list, csmesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(pmmesh->sel_criteria,
                              &(pmmesh->n_elts),
                              pmmesh->elt_list);

    BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);

    BFT_MALLOC(pmmesh->new_to_old, pmmesh->n_elts, cs_lnum_t);

    _assign_cell_mesh(csmesh,
                      pmmesh->n_elts,
                      pmmesh->elt_list,
                      pmmesh->med_mesh,
                      pmmesh->new_to_old);

    // BBOX
    if (use_bbox) {
      if (pmmesh->bbox == NULL) {
        BFT_MALLOC(pmmesh->bbox, 6, cs_real_t);
      }
      pmmesh->med_mesh->getBoundingBox(pmmesh->bbox);
    }

  } else if (pmmesh->elt_dim == 2) {

    /* Creation of a new nodal mesh from selected border faces */

    BFT_MALLOC(pmmesh->elt_list, csmesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(pmmesh->sel_criteria,
                                &(pmmesh->n_elts),
                                pmmesh->elt_list);

    BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);

    BFT_MALLOC(pmmesh->new_to_old, pmmesh->n_elts, cs_lnum_t);

    _assign_face_mesh(csmesh,
                      pmmesh->n_elts,
                      pmmesh->elt_list,
                      pmmesh->med_mesh,
                      pmmesh->new_to_old);

  }
#endif
}


/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   create a new cs_medcoupling_mesh_t instance based on cs_mesh_t
 *
 * \param[in] csmesh              pointer to cs_mesh_t instance
 * \param[in] name                name of the mesh
 * \param[in] selection_criteria  selection criteria (entire mesh or part of it)
 * \param[in] elt_dim             dimension of elements. 2: faces, 3: cells
 * \param[in] use_bbox            Use a reduced bounding box
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_from_base(cs_mesh_t  *csmesh,
                              const char *name,
                              const char *selection_criteria,
                              const int   elt_dim,
                              const int   use_bbox)

{
  cs_medcoupling_mesh_t *m = NULL;

#if !defined(HAVE_MEDCOUPLING)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else

  const char *_sel_crit =
    (selection_criteria != NULL) ? selection_criteria : "all[]";

  m = _get_mesh_from_parameters(_sel_crit, elt_dim);

  if (m == NULL) {
    m = _add_medcoupling_mesh(name, _sel_crit, elt_dim);
    _copy_mesh_from_base(csmesh, m, use_bbox);
  }

#endif

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_medcoupling_mesh_t
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_mesh_destroy(cs_medcoupling_mesh_t  *mesh)
{
  BFT_FREE(mesh->sel_criteria);
  BFT_FREE(mesh->elt_list);
  BFT_FREE(mesh->new_to_old);
  BFT_FREE(mesh->bbox);

#if defined(HAVE_MEDCOUPLING)
  mesh->med_mesh->decrRef();
#endif

  BFT_FREE(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all cs_medcoupling_mesh_t instances.
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_mesh_destroy_all(void)
{

  for (int i = 0; i < _n_sub_meshes; i++) {
    cs_medcoupling_mesh_destroy(_sub_meshes[i]);
    _sub_meshes[i] = NULL;
  }

  BFT_FREE(_sub_meshes);

  _sub_meshes   = NULL;
  _n_sub_meshes = 0;

}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's spatial dimension
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return associated spatial dimension
 */
/*----------------------------------------------------------------------------*/

int
cs_medcoupling_mesh_get_dim(cs_medcoupling_mesh_t  *m)
{
  int retval = -1;
  if (m != NULL)
    retval = m->elt_dim;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's number of elements
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return associated number of elements
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_medcoupling_mesh_get_n_elts(cs_medcoupling_mesh_t  *m)
{
  cs_lnum_t retval = 0;
  if (m != NULL)
    retval = m->n_elts;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's (parent) elements list
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return ids of associated elements, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_medcoupling_mesh_get_elt_list(cs_medcoupling_mesh_t  *m)
{
  const cs_lnum_t *retval = NULL;

  if (m != NULL)
    retval = m->elt_list;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's (parent) elements list
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return ids of associated elements, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_medcoupling_mesh_get_connectivity(cs_medcoupling_mesh_t  *m)
{
  const cs_lnum_t *retval = NULL;

  if (m != NULL)
    retval = m->new_to_old;

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
