/*============================================================================
 * Usage of MEDCoupling base components.
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

#if defined(HAVE_MEDCOUPLING)

static const cs_lnum_t _perm_tri[3]  = {0, 2, 1};
static const cs_lnum_t _perm_quad[4] = {0, 3, 2, 1};
static const cs_lnum_t _perm_pent[5] = {0, 4, 3, 2, 1};

#endif

static int                     _n_sub_meshes = 0;
static cs_medcoupling_mesh_t **_sub_meshes = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MEDCOUPLING)

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assign vertex coordinates to a MEDCoupling mesh structure
 *
 * \param[in] mesh          pointer to mesh structure from which data is copied
 * \param[in] n_vtx         number of vertices to assign
 * \param[in, out] vtx_id   pointer to vertices id's used for assigning
 *                          (-1 if not present/0 if present in, id out)
 * \param[in, out] pmmesh   pointer to associated cs_medcoupling_mesh_t
 */
/*----------------------------------------------------------------------------*/

static void
_assign_vertex_coords(const cs_mesh_t        *mesh,
                      cs_lnum_t               n_vtx,
                      cs_lnum_t              *vtx_id,
                      cs_medcoupling_mesh_t  *pmmesh)
{
  MEDCouplingUMesh  *med_mesh = pmmesh->med_mesh;

  const cs_lnum_t  dim = mesh->dim;
  const cs_coord_t  *vertex_coords = mesh->vtx_coord;

  assert(pmmesh->n_elts > 0);

  /* Assign vertex list if not all vertices are selected
     (so as to be able to handle vertex data */

  pmmesh->n_vtx = n_vtx;

  if (n_vtx >= mesh->n_vertices) {
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
      vtx_id[i] = i;
  }
  else {
    BFT_REALLOC(pmmesh->vtx_list, n_vtx, cs_lnum_t);
    cs_lnum_t vtx_count = 0;
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
      if (vtx_id[i] > -1) {
        vtx_id[i] = vtx_count;
        pmmesh->vtx_list[vtx_count] = i;
        vtx_count++;
      }
    }
    assert(vtx_count == n_vtx);
  }

  /* Assign all coordinates */
  /*------------------------*/

  DataArrayDouble *med_coords = DataArrayDouble::New();
  med_coords->alloc(n_vtx, dim);

  if (pmmesh->vtx_list != NULL) {
    for (cs_lnum_t i = 0; i < pmmesh->n_vtx; i++) {
      cs_lnum_t v_id = pmmesh->vtx_list[i];
      for (cs_lnum_t j = 0; j < dim; j++) {
        med_coords->setIJ(i, j, vertex_coords[v_id*dim + j]);
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
 * \param[in] mesh         pointer to cs_mesh_t structure
 * \param[in, out] pmmesh  pointer to associated cs_medcoupling_mesh_t
 */
/*----------------------------------------------------------------------------*/

static void
_assign_face_mesh(const cs_mesh_t        *mesh,
                  cs_medcoupling_mesh_t  *pmmesh)
{
  const cs_lnum_t   n_elts = pmmesh->n_elts;
  const cs_lnum_t  *elts_list = pmmesh->elt_list;
  cs_lnum_t        *new_to_old = pmmesh->new_to_old;
  MEDCouplingUMesh  *med_mesh = pmmesh->med_mesh;

  INTERP_KERNEL::NormalizedCellType type;

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

  /* Mark vertices (-1 if unused, 0 if used) */

  cs_lnum_t vtx_count = 0;

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);
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
        if (vtx_id[vid] < 0) {
          vtx_id[vid] = 0;  /* mark only */
          vtx_count++;
        }
      }
    }

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t j = mesh->b_face_vtx_idx[i];
           j < mesh->b_face_vtx_idx[i+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0) {
          vtx_id[vid] = 0;  /* mark only */
          vtx_count++;
        }
      }
    }

  }

  /* Assign coordinates and renumber vertices */

  _assign_vertex_coords(mesh, vtx_count, vtx_id, pmmesh);

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
 * \param[in] mesh         pointer to cs_mesh_t structure
 * \param[in, out] pmmesh  pointer to associated cs_medcoupling_mesh_t
 */
/*----------------------------------------------------------------------------*/

static void
_assign_cell_mesh(const cs_mesh_t        *mesh,
                  cs_medcoupling_mesh_t  *pmmesh)
{
  const cs_lnum_t   n_elts = pmmesh->n_elts;
  const cs_lnum_t  *elts_list = pmmesh->elt_list;
  cs_lnum_t        *new_to_old = pmmesh->new_to_old;
  MEDCouplingUMesh  *med_mesh = pmmesh->med_mesh;

  INTERP_KERNEL::NormalizedCellType type;

  cs_lnum_t cell_count = 0;

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

  /* Mark vertices (-1 if unused, 0 if used) */

  cs_lnum_t vtx_count = 0;

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
        if (vtx_id[vid] < 0) {
          vtx_id[vid] = 0;  /* mark only */
          vtx_count++;
        }
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
        if (vtx_id[vid] < 0) {
          vtx_id[vid] = 0;  /* mark only */
          vtx_count++;
        }
      }
    }
  }

  /* Assign coordinates and renumber vertices */

  _assign_vertex_coords(mesh, vtx_count, vtx_id, pmmesh);

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assign an empty mesh
 *
 * \param[in, out] pmmesh   pointer to associated cs_medcoupling_mesh_t
 */
/*----------------------------------------------------------------------------*/

static void
_assign_empty_mesh(cs_medcoupling_mesh_t *pmmesh)
{
  MEDCouplingUMesh  *med_mesh = pmmesh->med_mesh;

  //  Vertices
  DataArrayDouble *med_coords = DataArrayDouble::New();
  med_coords->alloc(0, 3);
  med_mesh->setCoords(med_coords);
  med_coords->decrRef();

  // Elements
  med_mesh->allocateCells(0);
  med_mesh->finishInsertingCells();
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

static cs_medcoupling_mesh_t *
_get_mesh_from_criteria(const char  *selection_criteria,
                        int          elt_dim)
{
  cs_medcoupling_mesh_t *m = NULL;

  for (int i = 0; i < _n_sub_meshes; i++) {
    cs_medcoupling_mesh_t *mt = _sub_meshes[i];
    if (   elt_dim == mt->elt_dim
        && strcmp(mt->sel_criteria, selection_criteria) == 0) {
      m = mt;
      break;
    }
  }

  return m;
}

#if defined(HAVE_MEDCOUPLING)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to cs_medcoupling_mesh_t based on name
 *        and element dimension.
 *
 * \param[in]  name     name
 * \param[in]  elt_dim  dimension of elements. 2: faces, 3: cells
 *
 * \return pointer to found mesh instance, NULL if none found.
 */
/*----------------------------------------------------------------------------*/

static cs_medcoupling_mesh_t *
_get_mesh_from_name(const char  *name,
                    int          elt_dim)
{
  cs_medcoupling_mesh_t *m = NULL;

#if !defined(HAVE_MEDCOUPLING)
  CS_UNUSED(name);
  CS_UNUSED(elt_dim);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else

  for (int i = 0; i < _n_sub_meshes; i++) {
    cs_medcoupling_mesh_t *mt = _sub_meshes[i];

    if (elt_dim == mt->elt_dim) {
      std::string m_name = mt->med_mesh->getName();
      if (strcmp(m_name.c_str(), name) == 0) {
        m = mt;
        break;
      }
    }
  }

#endif

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a cs_medcoupling_mesh_t to the list
 *
 * \param[in] name      name of the mesh
 * \param[in] elt_dim   dimension of elements (2: faces, 3: cells)
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

static cs_medcoupling_mesh_t *
_add_medcoupling_mesh(const char  *name,
                      int          elt_dim)
{
  cs_medcoupling_mesh_t *m = NULL;

#if !defined(HAVE_MEDCOUPLING)
  CS_UNUSED(name);
  CS_UNUSED(elt_dim);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else

  BFT_REALLOC(_sub_meshes, _n_sub_meshes + 1, cs_medcoupling_mesh_t *);

  BFT_MALLOC(m, 1, cs_medcoupling_mesh_t);

  m->sel_criteria = NULL;
  m->elt_dim  = elt_dim;
  m->n_elts   = 0;
  m->n_vtx    = 0;
  m->elt_list = NULL;
  m->vtx_list = NULL;

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
 * \brief Build selection lists based on criteria.
 *
 * \param[in] mesh         pointer to cs_mesh_t structure
 * \param[in, out] pmmesh  pointer to associated cs_medcoupling_mesh_t
 */
/*----------------------------------------------------------------------------*/

static void
_select_from_criteria(cs_mesh_t              *csmesh,
                      cs_medcoupling_mesh_t  *pmmesh)
{
  if (pmmesh->elt_dim == 3) {

    BFT_MALLOC(pmmesh->elt_list, csmesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(pmmesh->sel_criteria,
                              &(pmmesh->n_elts),
                              pmmesh->elt_list);

    BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);

  }
  else if (pmmesh->elt_dim == 2) {

    BFT_MALLOC(pmmesh->elt_list, csmesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(pmmesh->sel_criteria,
                                &(pmmesh->n_elts),
                                pmmesh->elt_list);

    BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);

  }
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

static void
_copy_mesh_from_base(cs_mesh_t              *csmesh,
                     cs_medcoupling_mesh_t  *pmmesh,
                     int                     use_bbox)
{
#if !defined(HAVE_MEDCOUPLING)
  CS_UNUSED(csmesh);
  CS_UNUSED(pmmesh);
  CS_UNUSED(use_bbox);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else
  // Assign mesh only if non empty list of elements
  if (pmmesh->n_elts == 0) {
    _assign_empty_mesh(pmmesh);
    return;
  }
  else {
    if (pmmesh->elt_dim == 3) {
      /* Creation of a new nodal mesh from selected cells */
      BFT_MALLOC(pmmesh->new_to_old, pmmesh->n_elts, cs_lnum_t);

      _assign_cell_mesh(csmesh, pmmesh);

      // BBOX
      if (use_bbox) {
        if (pmmesh->bbox == NULL) {
          BFT_MALLOC(pmmesh->bbox, 6, cs_real_t);
        }
        pmmesh->med_mesh->getBoundingBox(pmmesh->bbox);
      }
    }
      else if (pmmesh->elt_dim == 2) {

      /* Creation of a new nodal mesh from selected border faces */

      if (pmmesh->sel_criteria != NULL) {
        BFT_MALLOC(pmmesh->elt_list, csmesh->n_b_faces, cs_lnum_t);

        cs_selector_get_b_face_list(pmmesh->sel_criteria,
                                    &(pmmesh->n_elts),
                                    pmmesh->elt_list);

        BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);
      }

      BFT_MALLOC(pmmesh->new_to_old, pmmesh->n_elts, cs_lnum_t);

      _assign_face_mesh(csmesh, pmmesh);

    }
  } /* if n_etls == 0 */
#endif
}

#endif /* defined(HAVE_MEDCOUPLING) */

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
 * \param[in] selection_criteria  selection criteria string
 * \param[in] elt_dim             dimension of elements (2: faces, 3: cells)
 * \param[in] use_bbox            Use a reduced bounding box
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_from_base(cs_mesh_t   *csmesh,
                              const char  *name,
                              const char  *selection_criteria,
                              int          elt_dim,
                              int          use_bbox)
{
  cs_medcoupling_mesh_t *m = NULL;

#if !defined(HAVE_MEDCOUPLING)
  CS_UNUSED(csmesh);
  CS_UNUSED(name);
  CS_UNUSED(selection_criteria);
  CS_UNUSED(elt_dim);
  CS_UNUSED(use_bbox);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else

  const char *_sel_crit
    = (selection_criteria != NULL) ? selection_criteria : "all[]";

  m = _get_mesh_from_criteria(_sel_crit, elt_dim);

  if (m == NULL) {
    m = _add_medcoupling_mesh(name, elt_dim);

    size_t len_sel_crit = strlen(_sel_crit);
    BFT_MALLOC(m->sel_criteria, len_sel_crit+1, char);
    strcpy(m->sel_criteria, _sel_crit);
    m->sel_criteria[len_sel_crit] = '\0';

    _select_from_criteria(csmesh, m);

    _copy_mesh_from_base(csmesh, m, use_bbox);
  }

#endif

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   create a new cs_medcoupling_mesh_t instance based on cs_mesh_t
 *
 * \param[in] csmesh      pointer to cs_mesh_t instance
 * \param[in] name        name of the mesh
 * \param[in] n_elts      local number of elements
 * \param[in] elt_ids     list of local elements
 * \param[in] elt_dim     dimension of elements (2: faces, 3: cells)
 * \param[in] use_bbox    use a reduced bounding box
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_from_ids(cs_mesh_t       *csmesh,
                             const char      *name,
                             cs_lnum_t        n_elts,
                             const cs_lnum_t  elt_ids[],
                             int              elt_dim,
                             int              use_bbox)
{
  cs_medcoupling_mesh_t *m = NULL;

#if !defined(HAVE_MEDCOUPLING)
  CS_UNUSED(csmesh);
  CS_UNUSED(name);
  CS_UNUSED(n_elts);
  CS_UNUSED(elt_ids);
  CS_UNUSED(elt_dim);
  CS_UNUSED(use_bbox);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));
#else

  m = _get_mesh_from_name(name, elt_dim);

  if (m == NULL) {
    m = _add_medcoupling_mesh(name, elt_dim);

    m->n_elts = n_elts;
    BFT_MALLOC(m->elt_list, n_elts, cs_lnum_t);
    memcpy(m->elt_list, elt_ids, n_elts*sizeof(cs_lnum_t));

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
  BFT_FREE(mesh->vtx_list);
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
 * \brief Return a cs_medcoupling_mesh_t structure's number of vertices
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return associated number of vertices
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_medcoupling_mesh_get_n_vertices(cs_medcoupling_mesh_t  *m)
{
  cs_lnum_t retval = 0;
  if (m != NULL)
    retval = m->n_vtx;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's (parent) vertices list
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return ids of associated vertices, or NULL if all or no local vertices
 *         o parent mesh are present.
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_medcoupling_mesh_get_vertex_list(cs_medcoupling_mesh_t  *m)
{
  const cs_lnum_t *retval = NULL;

  if (m != NULL)
    retval = m->vtx_list;

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
/*!
 * \brief Returns a pointer to a MEDCouplingUMesh of a plane.
 *
 * \param[in] origin   Plane origin coordinates
 * \param[in] normal   Plane normal vector
 * \param[in] length1  Plane's edge length along first axis
 * \param[in] length2  Plane's edge length along second axis
 *
 * \return pointer to the MEDCouplingUMesh structure.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
MEDCouplingUMesh *
#else
void *
#endif
cs_medcoupling_create_plane_mesh(const cs_real_t origin[],
                                 const cs_real_t normal[],
                                 const cs_real_t length1,
                                 const cs_real_t length2)
{
#if defined(HAVE_MEDCOUPLING)
  /* Compute plane orthonormal basis */
  cs_real_t basis[3][3] = {{0.}, {0.}, {0.}};
  cs_math_3_orthonormal_basis(normal, basis);

  /* Compute plane coordinates */
  cs_real_t coords[12] = {0.};
  double dx[4] = {-0.5, -0.5, 0.5, 0.5};
  double dy[4] = {-0.5, 0.5, 0.5, -0.5};

  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 3; i++)
      coords[j*3 + i] = origin[i]
                      + dx[j] * length1 * basis[0][i]
                      + dy[j] * length2 * basis[1][i];
  }

  /* Generate mesh */
  mcIdType conn[4] = {0, 1, 2, 3};

  MEDCouplingUMesh *m = MEDCouplingUMesh::New();

  m->setMeshDimension(2);
  m->allocateCells(1);
  m->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn);
  m->finishInsertingCells();

  DataArrayDouble *localCoords=DataArrayDouble::New();
  localCoords->alloc(4,3);

  std::copy(coords, coords+12, localCoords->getPointer());
  m->setCoords(localCoords);
  localCoords->decrRef();

  return m;
#else
  CS_UNUSED(origin);
  CS_UNUSED(normal);
  CS_UNUSED(length1);
  CS_UNUSED(length2);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));

  return NULL;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a pointer to a MEDCouplingUMesh of a disc.
 *
 * \param[in] origin     Disc origin coordinates
 * \param[in] normal     Disc normal vector
 * \param[in] radius     Disc radius
 * \param[in] n_sectors  Number of sectors for discretization. If negative,
 *                       default value of 36 is taken.
 *
 * \return pointer to the MEDCouplingUMesh structure.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
MEDCouplingUMesh *
#else
void *
#endif
cs_medcoupling_create_disc_mesh(const cs_real_t origin[],
                                const cs_real_t normal[],
                                const cs_real_t radius,
                                const int       n_sectors)
{
#if defined(HAVE_MEDCOUPLING)
  /* Compute number of sectors. Default is 36 if user does not provide one */
  const int _n_sectors = (n_sectors > 0) ? n_sectors : 36;

  /* Compute plane orthonormal basis */
  cs_real_t basis[3][3] = {{0.}, {0.}, {0.}};
  cs_math_3_orthonormal_basis(normal, basis);

  /* Compute coordinates for disc nodes.
   * There n_sectors+1 nodes, one for center, and one per sector.
   */
  std::vector<cs_real_t> coords(3*(_n_sectors+1), 0.);
  std::copy(origin, origin+3, coords.begin());

  const cs_real_t theta = 2. * cs_math_pi / _n_sectors;

  for (int j = 0; j < _n_sectors; j++) {
    cs_real_t cos_j = (j == 0) ? 1. : cos(theta * j);
    cs_real_t sin_j = (j == 0) ? 0. : sin(theta * j);
    for (int i = 0; i < 3; i++) {
      coords[(j+1)*3 + i] = origin[i]
                          + cos_j * radius * basis[0][i]
                          + sin_j * radius * basis[1][i];
    }
  }

  /* Generate mesh */
  MEDCouplingUMesh *m = MEDCouplingUMesh::New();
  m->setMeshDimension(2);
  m->allocateCells(_n_sectors);

  for (int i = 0; i < _n_sectors; i++) {
    /*
     * N3       N2
     * .-------.
     *  \     /
     *   \   /
     *    \./
     *     N1
     *
     * A sector has 3 nodes (TRI3).
     * N1 is the center of the disc (origin), hence node #0
     * N2 value is in [1, N_SECTORS]
     * N3 value is in [2,..., N_SECTORS, 1] -> Hence the modulo call
     */
    mcIdType conn[3] = {0, 1 + i, 1 + ((1+i) % _n_sectors)};
    m->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn);
  }
  m->finishInsertingCells();

  DataArrayDouble *localCoords=DataArrayDouble::New();
  localCoords->alloc(_n_sectors + 1, 3);
  std::copy(coords.begin(), coords.end(), localCoords->getPointer());
  m->setCoords(localCoords);
  localCoords->decrRef();

  return m;
#else
  CS_UNUSED(origin);
  CS_UNUSED(normal);
  CS_UNUSED(radius);
  CS_UNUSED(n_sectors);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));


  return NULL;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a pointer to a MEDCouplingUMesh of an annulus
 *
 * \param[in] origin     Annulus origin coordinates
 * \param[in] normal     Annulus normal vector
 * \param[in] radius1    Annulus inner radius
 * \param[in] radius2    Annulus outer radius
 * \param[in] n_sectors  Number of sectors for discretization. If negative,
 *                       default value of 36 is taken.
 *
 * \return pointer to the MEDCouplingUMesh structure.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
MEDCouplingUMesh *
#else
void *
#endif
cs_medcoupling_create_annulus_mesh(const cs_real_t origin[],
                                   const cs_real_t normal[],
                                   const cs_real_t radius1,
                                   const cs_real_t radius2,
                                   const int       n_sectors)
{
#if defined(HAVE_MEDCOUPLING)
  /* Compute number of sectors. Default is 36 if user does not provide one */
  const int _n_sectors = (n_sectors > 0) ? n_sectors : 36;

  /* Compute plane orthonormal basis */
  cs_real_t basis[3][3] = {{0.}, {0.}, {0.}};
  cs_math_3_orthonormal_basis(normal, basis);

  /* Compute coordinates for disc nodes.
   * There n_sectors*2 nodes, n_sectors for each radius.
   */
  std::vector<cs_real_t> coords(3*2*_n_sectors, 0.);

  const cs_real_t theta = 2. * cs_math_pi / _n_sectors;

  for (int j = 0; j < _n_sectors; j++) {
    cs_real_t cos_j = (j == 0) ? 1. : cos(theta * j);
    cs_real_t sin_j = (j == 0) ? 0. : sin(theta * j);
    for (int i = 0; i < 3; i++) {
      /* Inner circle, radius1 */
      coords[j*3 + i] = origin[i]
                        + cos_j * radius1 * basis[0][i]
                        + sin_j * radius1 * basis[1][i];

      /* Outer circle, radius2 */
      coords[(j+_n_sectors)*3 + i] = origin[i]
                                   + cos_j * radius2 * basis[0][i]
                                   + sin_j * radius2 * basis[1][i];
    }
  }

  /* Generate mesh */
  MEDCouplingUMesh *m = MEDCouplingUMesh::New();
  m->setMeshDimension(2);
  m->allocateCells(_n_sectors);

  for (int i = 0; i < _n_sectors; i++) {
    int ii = i;
    int jj = (i+1) % _n_sectors;

    mcIdType conn[4] = {jj, jj+_n_sectors, ii+_n_sectors, ii};
    m->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn);
  }
  m->finishInsertingCells();

  DataArrayDouble *localCoords=DataArrayDouble::New();
  localCoords->alloc(2*_n_sectors, 3);
  std::copy(coords.begin(), coords.end(), localCoords->getPointer());
  m->setCoords(localCoords);
  localCoords->decrRef();

  return m;
#else
  CS_UNUSED(origin);
  CS_UNUSED(normal);
  CS_UNUSED(radius1);
  CS_UNUSED(radius2);
  CS_UNUSED(n_sectors);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: this funnction cannot be called without "
              "MEDCoupling support\n"));

  return NULL;
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
