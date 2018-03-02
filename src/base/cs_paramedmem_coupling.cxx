/*============================================================================
 * ParaMEDMEM coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
#include <MEDCoupling_version.h>

#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>
#include "MEDCouplingRemapper.hxx"

#if defined(HAVE_MEDCOUPLING_LOADER)
#include <MEDLoader.hxx>
#endif

#if defined(HAVE_PARAMEDMEM)
#include <ParaFIELD.hxx>
#include <ParaMESH.hxx>
#include <InterpKernelDEC.hxx>
#include <OverlapDEC.hxx>
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

#include "cs_paramedmem_coupling.hxx"

/*----------------------------------------------------------------------------*/

using namespace MEDCoupling;

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * ParaMEDMED field structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int                       mesh_id;       /* Associated mesh structure id */
  int                       dim;           /* Field dimension */

  TypeOfTimeDiscretization  td;            /* NO_TIME, ONE_TIME, LINEAR_TIME,
                                              or CONST_ON_TIME_INTERVAL */

  MEDCouplingFieldDouble   *f;             /* Pointer to MED coupling field */

  ParaFIELD                *pf;            /* Pointer to ParaMEDMEM field */

} _paramedmem_field_t;

/*----------------------------------------------------------------------------
 * PAraMEDMED mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char               *sel_criteria;   /* Element selection criteria */

  int                 direction;      /* 1: send, 2: receive, 3: both */

  int                 elt_dim;        /* Element dimension */

  cs_lnum_t           n_elts;         /* Number of coupled elements */
  cs_lnum_t          *elt_list;       /* List of associated elements
                                         (0 to n-1) */

  MEDCouplingUMesh   *med_mesh;       /* MED mesh structure */
  ParaMESH           *para_mesh[2];   /* parallel MED mesh structures
                                         for send (0) and receive (1) */

  int                *new_to_old;     /* Connectivity used if only a section of
                                         the mesh is read */
  cs_real_t          *bbox;

} _paramedmem_mesh_t;

/*----------------------------------------------------------------------------
 * MEDCoupling writer/reader structure
 *----------------------------------------------------------------------------*/

struct _cs_paramedmem_coupling_t {

  char                      *name;           /* Coupling name */

  int                        n_meshes;       /* Number of meshes */
  _paramedmem_mesh_t       **meshes;         /* Array of mesh helper
                                                structures */

  int                        n_fields;       /* Number of fields */
  _paramedmem_field_t      **fields;         /* Array of field helper
                                                structures */

  InterpKernelDEC           *send_dec;       /* Send data exchange channel */
  InterpKernelDEC           *recv_dec;       /* Receive data exchange channel */
};

/*----------------------------------------------------------------------------
 * MEDCoupling parallel interpolation structure
 *----------------------------------------------------------------------------*/

struct _cs_medcoupling_remapper_t {

  char                     *name;
  char                     *medfile_path;
  char                     *mesh_name;
  char                    **field_names;

  int                       n_fields;

  _paramedmem_mesh_t       *target_mesh;

  MEDCouplingUMesh         *bbox_source_mesh;

  MEDCouplingFieldDouble  **source_fields;

  MEDCouplingRemapper      *remapper;     /* MEDCoupling remapper */

};

/*============================================================================
 * Private global variables
 *============================================================================*/

static int                          _n_remappers = 0;
static cs_medcoupling_remapper_t  **_remapper = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Assign vertex coordinates to MEDCoupling structures
 *
 * parameters:
 *   mesh     <-- pointer to parent mesh structure that should be written.
 *   n_vtx    <-- number of vertices to assign.
 *   vtx_id   <-- id of initial vertices in new (sub) mesh, or -1
 *                (size: mesh->n_vtx).
 *   med_mesh <-> pointer to associated MEDCoupling mesh structure.
 *----------------------------------------------------------------------------*/

static void
_assign_vertex_coords(const cs_mesh_t    *mesh,
                      cs_lnum_t           n_vtx,
                      const cs_lnum_t    *vtx_id,
                      MEDCouplingUMesh   *med_mesh)
{
  int  i, j;

  const int  dim = mesh->dim;
  const cs_coord_t  *vertex_coords = mesh->vtx_coord;

  assert(med_mesh != NULL);

  /* Assign all coordinates */
  /*------------------------*/

  DataArrayDouble *med_coords = DataArrayDouble::New();
  med_coords->alloc(n_vtx, dim);

  if (vtx_id != NULL) {
    for (i = 0; i < mesh->n_vertices; i++) {
      if (vtx_id[i] > -1) {
        for (j = 0; j < dim; j++) {
          med_coords->setIJ(vtx_id[i], j, vertex_coords[i*dim + j]);
        }
      }
    }
  }
  else {
    for (i = 0; i < mesh->n_vertices; i++) {
      for (j = 0; j < dim; j++)
        med_coords->setIJ(i, j, vertex_coords[i*dim + j]);
    }
  }

  med_mesh->setCoords(med_coords);
  med_coords->decrRef();

  return;
}

/*----------------------------------------------------------------------------
 * Assign boundary faces to a MEDCoupling structure
 *
 * parameters:
 *   mesh      <-- parent mesh structure
 *   n_elts    <-- number of selected elements
 *   elts_list <-- list of selected elements (1 to n numbering),
 *                 or NULL (if all are selected)
 *   med_mesh  <-> associated MEDCouplingUMesh structure
 *----------------------------------------------------------------------------*/

static void
_assign_face_mesh(const cs_mesh_t   *mesh,
                  cs_lnum_t          n_elts,
                  const cs_lnum_t   *elts_list,
                  MEDCouplingUMesh  *med_mesh)
{
  cs_lnum_t i, j;
  INTERP_KERNEL::NormalizedCellType type;

  cs_lnum_t vtx_count = -1;
  int elt_buf_size = 4;
  int *elt_buf = NULL;
  cs_lnum_t *vtx_id = NULL;

  const int perm_tri[4] = {0, 2, 1};
  const int perm_quad[4] = {0, 3, 2, 1};

  /* Mark and renumber vertices */

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);

  /* Initialize the value of vtx_id */
  for (i = 0; i < mesh->n_vertices; i++)
    vtx_id[i] = -1;

  /* Case with filter list */

  if (elts_list != NULL) {

    for (i = 0; i < n_elts; i++) {
      cs_lnum_t eid = elts_list[i];
      for (j = mesh->b_face_vtx_idx[eid];
           j < mesh->b_face_vtx_idx[eid+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }

  } else {

    for (i = 0; i < n_elts; i++) {
      for (j = mesh->b_face_vtx_idx[i];
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

  BFT_MALLOC(elt_buf, elt_buf_size, int);
  med_mesh->allocateCells(n_elts);

  for (i = 0; i < n_elts; i++) {

    cs_lnum_t eid = (elts_list != NULL) ? elts_list[i] : i;

    assert(eid >= 0 && eid < mesh->n_b_faces);

    int n_vtx = mesh->b_face_vtx_idx[eid+1] - mesh->b_face_vtx_idx[eid];

    cs_lnum_t connect_start = mesh->b_face_vtx_idx[eid];

    if (n_vtx > elt_buf_size) { /* reallocate buffer if required */
      elt_buf_size *= 2;
      BFT_REALLOC(elt_buf, elt_buf_size, int);
    }

    switch(n_vtx) {
    case 3:
      type = INTERP_KERNEL::NORM_TRI3;
      for (j = 0; j < 3; j++)
        elt_buf[j]
          = vtx_id[mesh->b_face_vtx_lst[connect_start + perm_tri[j]]];
      break;
    case 4:
      type = INTERP_KERNEL::NORM_QUAD4;
      for (j = 0; j < 4; j++)
        elt_buf[j]
          = vtx_id[mesh->b_face_vtx_lst[connect_start + perm_quad[j]]];
      break;
    default:
      type = INTERP_KERNEL::NORM_POLYGON;
      for (j = 0; j < n_vtx; j++)
        elt_buf[j]
          = vtx_id[mesh->b_face_vtx_lst[connect_start + n_vtx - 1 - j]];
      break;
    }

    med_mesh->insertNextCell(type, n_vtx, elt_buf);

  }

  med_mesh->finishInsertingCells();

  BFT_FREE(elt_buf);
  BFT_FREE(vtx_id);
}

/*----------------------------------------------------------------------------
 * Assign cells to a MEDCoupling structure
 *
 * parameters:
 *   mesh      <-- parent mesh structure
 *   n_elts    <-- number of selected elements
 *   elts_list <-- list of selected elements (1 to n numbering),
 *                 or NULL (if all are selected)
 *   med_mesh  <-> associated MEDCouplingUMesh structure
 *----------------------------------------------------------------------------*/

static void
_assign_cell_mesh(const cs_mesh_t   *mesh,
                  cs_lnum_t          n_elts,
                  const cs_lnum_t   *elts_list,
                  MEDCouplingUMesh  *med_mesh,
                  int               *new_to_old)
{
  cs_lnum_t i, j, k;
  cs_lnum_t  c_id, c_id1, c_id2, face_id;
  INTERP_KERNEL::NormalizedCellType type;

  cs_lnum_t vtx_count = 0, cell_count = 0;

  int elt_buf_size = 8;
  int *elt_buf = NULL;
  cs_lnum_t *vtx_id = NULL;
  cs_lnum_t *cell_id = NULL;
  cs_lnum_t *cell_faces_idx = NULL, *cell_faces_num = NULL;

  /* Build old->new cell id indirection */

  BFT_MALLOC(cell_id, mesh->n_cells, cs_lnum_t);
  for (i = 0; i < mesh->n_cells; i++)
    cell_id[i] = -1;

  if (elts_list != NULL) {
    for (i = 0; i < n_elts; i++){
      cell_id[elts_list[i]] = cell_count++;
    }
    for (int ii = 0; ii < n_elts; ii++) {
      new_to_old[ cell_id[elts_list[ii]] ] = elts_list[ii];
    }
  }
  else {
    for (i = 0; i < n_elts; i++)
      cell_id[i] = cell_count++;
  }

  /* Mark and renumber vertices */

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);
  for (int vid = 0; vid < mesh->n_vertices; vid++) {
    vtx_id[vid] = -1;
  }

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
    c_id = cell_id[mesh->b_face_cells[face_id]];
    if (c_id > -1) {
      for (j = mesh->b_face_vtx_idx[face_id];
           j < mesh->b_face_vtx_idx[face_id+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }
  }

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    c_id1 = mesh->i_face_cells[face_id][0];
    c_id2 = mesh->i_face_cells[face_id][1];
    c_id1 = (c_id1 < mesh->n_cells) ? cell_id[c_id1] : -1;
    c_id2 = (c_id2 < mesh->n_cells) ? cell_id[c_id2] : -1;
    if (c_id1 > -1 || c_id2 > -1) {
      for (j = mesh->i_face_vtx_idx[face_id];
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
                                 mesh->n_cells,
                                 cell_id,
                                 &cell_faces_idx,
                                 &cell_faces_num);

  BFT_FREE(cell_id);

  /* Now loop on cells */

  const cs_lnum_t  face_num_shift[2] = {0,
                                        mesh->n_b_faces};

  const cs_lnum_t  *face_vertices_idx[2] = {mesh->b_face_vtx_idx,
                                            mesh->i_face_vtx_idx};
  const cs_lnum_t  *face_vertices_num[2] = {mesh->b_face_vtx_lst,
                                            mesh->i_face_vtx_lst};

  BFT_MALLOC(elt_buf, elt_buf_size, int);

  /* Allocate the cells array */
  med_mesh->allocateCells(n_elts);

  for (i = 0; i < n_elts; i++) {

    int n_vtx;
    cs_lnum_t cell_vtx[8];

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

      int connect_size = 0;

      n_vtx = 0;

      for (j = cell_faces_idx[i]; j < cell_faces_idx[i]; j++) {

        face_id = CS_ABS(cell_faces_num[j-1]) - 1;

        int fl = 1;
        if (face_id < face_num_shift[fl])
          fl = 0;
        face_id -= face_num_shift[fl];

        cs_lnum_t v_id_start = face_vertices_idx[fl][face_id] - 1;
        cs_lnum_t v_id_end   = face_vertices_idx[fl][face_id + 1] - 1;
        int n_face_vertices = v_id_end - v_id_start;

        connect_size += n_face_vertices + 1;
        if (connect_size > elt_buf_size) { /* reallocate buffer if required */
          elt_buf_size *= 2;
          BFT_REALLOC(elt_buf, elt_buf_size, int);
        }

        /* Add separator after first face */
        if (j > cell_faces_idx[i])
          elt_buf[n_vtx++] = -1;

        /* Add face vertices */
        if (cell_faces_num[j-1] > 0) {
          for (k = v_id_start; k < v_id_end; k++)
            elt_buf[n_vtx++] = vtx_id[face_vertices_num[fl][k] - 1];
        }
        else {
          for (k = v_id_end - 1; k >= v_id_start; k--)
            elt_buf[n_vtx++] = vtx_id[face_vertices_num[fl][k] - 1];
        }

      }

    }

    med_mesh->insertNextCell(type, n_vtx, elt_buf);
  }

  med_mesh->finishInsertingCells();

  BFT_FREE(elt_buf);
  BFT_FREE(cell_faces_num);;
  BFT_FREE(cell_faces_idx);
  BFT_FREE(vtx_id);
}

/*----------------------------------------------------------------------------
 * Copy a fvm format mesh to a paramedmem mesh structure
 *
 * parameters:
 *   csmesh  <-- Code_Saturne FVM format mesh structure
 *   pmmesh  <-> partially ParaMEDMEM mesh coupling structure
 *----------------------------------------------------------------------------*/

static void
_copy_cs_mesh_to_med(cs_mesh_t          *csmesh,
                     _paramedmem_mesh_t *pmmesh)
{

  if (pmmesh->elt_dim == 3) {
  /* Creation of a new nodal mesh from selected cells */

    BFT_MALLOC(pmmesh->elt_list, csmesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(pmmesh->sel_criteria,
                              &(pmmesh->n_elts),
                              pmmesh->elt_list);

    BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);

    BFT_MALLOC(pmmesh->new_to_old, pmmesh->n_elts, int);

    _assign_cell_mesh(csmesh,
                      pmmesh->n_elts,
                      pmmesh->elt_list,
                      pmmesh->med_mesh,
                      pmmesh->new_to_old);

  } else if (pmmesh->elt_dim == 2) {
  /* Creation of a new nodal mesh from selected border faces */

    BFT_MALLOC(pmmesh->elt_list, csmesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(pmmesh->sel_criteria,
                                &(pmmesh->n_elts),
                                pmmesh->elt_list);

    BFT_REALLOC(pmmesh->elt_list, pmmesh->n_elts, cs_lnum_t);

    _assign_face_mesh(csmesh,
                      pmmesh->n_elts,
                      pmmesh->elt_list,
                      pmmesh->med_mesh);

  }

}

/*----------------------------------------------------------------------------
 * Initialize mesh for ParaMEDMEM coupling.
 *
 * parameters:
 *   coupling  <-- coupling structure.
 *   mesh      <-> partially ParaMEDMEM mesh coupling structure
 *----------------------------------------------------------------------------*/

static void
_init_mesh_coupling(cs_paramedmem_coupling_t  *coupling,
                    _paramedmem_mesh_t        *mesh)
{
  cs_mesh_t *parent_mesh = cs_glob_mesh;

  assert(mesh != NULL);

  /* Building the MED representation of the Saturne mesh */
  _copy_cs_mesh_to_med(parent_mesh, mesh);

  /* Define associated ParaMESH */
  mesh->para_mesh[0] = new ParaMESH(mesh->med_mesh,
                                    *(coupling->send_dec->getSourceGrp()),
                                    "source mesh");
  mesh->para_mesh[1] = new ParaMESH(mesh->med_mesh,
                                    *(coupling->recv_dec->getTargetGrp()),
                                    "target mesh");
}

/*----------------------------------------------------------------------------
 * Destroy coupled entity helper structure.
 *
 * parameters:
 *   coupling ent <-> pointer to structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_mesh(_paramedmem_mesh_t **mesh)
{
  _paramedmem_mesh_t *m = *mesh;

  if (m == NULL)
    return;

  for (int i = 0; i < 2; i++) {
    if (m->para_mesh[i] != NULL)
      delete m->para_mesh[i];
  }
  if (m->med_mesh != NULL)
    m->med_mesh = NULL; // delete m->med_mesh by decreasing reference counter;

  BFT_FREE(m->elt_list);

  BFT_FREE(*mesh);
}

/*----------------------------------------------------------------------------
 * Create an InterpKernelDEC object based on two lists, and their sizes,
 * of mpi ranks (within MPI_COMM_WORLD).
 *
 * parameters:
 *   grp1_global_ranks <-- array of ranks of group 1
 *   grp1_size         <-- size of grp1_global_ranks array
 *   grp2_global_ranks <-- array of ranks of group 2
 *   grp2_size         <-- size of grp2_global_ranks array
 *
 * return:
 *   new InterpKernelDEC object
 *----------------------------------------------------------------------------*/

static InterpKernelDEC *
_cs_paramedmem_create_InterpKernelDEC(int  *grp1_global_ranks,
                                      int   grp1_size,
                                      int  *grp2_global_ranks,
                                      int   grp2_size)
{
  /* Group 1 id's */
  std::set<int> grp1_ids;
  for (int i = 0; i < grp1_size; i++) {
    grp1_ids.insert(grp1_global_ranks[i]);
  }

  /* Group 2 id's */
  std::set<int> grp2_ids;
  for (int i = 0; i < grp2_size; i++) {
    grp2_ids.insert(grp2_global_ranks[i]);
  }

  /* Create the InterpKernel DEC */
  InterpKernelDEC *NewDec = new InterpKernelDEC(grp1_ids, grp2_ids);

  return NewDec;
}

#if defined(HAVE_MEDCOUPLING_LOADER)

/*----------------------------------------------------------------------------
 * Create a new cs_medcoupling_remapper_t * object.
 *
 * parameters:
 *   name            <-- new object name
 *   elt_dim         <-- element dimension
 *   select_criteria <-- selection criteria
 *   medfile_path    <-- path of associated MED file
 *   mesh_name       <-- mesh name
 *   n_fields        <-- number of fields
 *   field_names     <-- associated field names
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *
 * return:
 *   new remapper object
 *----------------------------------------------------------------------------*/

static cs_medcoupling_remapper_t *
_cs_paramedmem_create_remapper(const char   *name,
                               int           elt_dim,
                               const char   *select_criteria,
                               const char   *medfile_path,
                               const char   *mesh_name,
                               int           n_fields,
                               const char  **field_names,
                               int           iteration,
                               int           iteration_order)
{
  cs_medcoupling_remapper_t *r = NULL;
  BFT_MALLOC(r, 1, cs_medcoupling_remapper_t);

  r->n_fields = n_fields;

  BFT_MALLOC(r->name, strlen(name)+1, char);
  strcpy(r->name, name);

  // Store fields and medfile info in case updates are needed

  BFT_MALLOC(r->medfile_path, strlen(medfile_path)+1, char);
  strcpy(r->medfile_path, medfile_path);

  BFT_MALLOC(r->mesh_name, strlen(mesh_name)+1, char);
  strcpy(r->mesh_name, mesh_name);

  BFT_MALLOC(r->field_names, n_fields, char *);
  for (int i = 0; i < n_fields; i++) {
    BFT_MALLOC(r->field_names[i], strlen(field_names[i])+1, char);
    strcpy(r->field_names[i], field_names[i]);
  }

  // New MEDCoupling UMesh linked to Code_Saturne mesh

  _paramedmem_mesh_t *new_mesh = NULL;
  BFT_MALLOC(new_mesh, 1, _paramedmem_mesh_t);
  BFT_MALLOC(new_mesh->sel_criteria, strlen(select_criteria)+1, char);
  strcpy(new_mesh->sel_criteria, select_criteria);

  new_mesh->elt_dim  = elt_dim;
  new_mesh->n_elts   = 0;
  new_mesh->elt_list = NULL;

  new_mesh->bbox = NULL;

  new_mesh->med_mesh = MEDCouplingUMesh::New();
  new_mesh->med_mesh->setName("target_mesh");
  new_mesh->med_mesh->setTimeUnit("s");
  new_mesh->med_mesh->setMeshDimension(elt_dim);

  new_mesh->para_mesh[0] = NULL;
  new_mesh->para_mesh[1] = NULL;

  cs_mesh_t *parent_mesh = cs_glob_mesh;

  _copy_cs_mesh_to_med(parent_mesh, new_mesh);

  if (new_mesh->bbox == NULL) {
    BFT_MALLOC(new_mesh->bbox, 6, cs_real_t);
  }

  new_mesh->med_mesh->getBoundingBox(new_mesh->bbox);

  r->target_mesh = new_mesh;

  // MEDCoupling remapper (sequential interpolation)

  r->remapper = new MEDCouplingRemapper;
  r->remapper->setPrecision(1.0e-12);
  r->remapper->setIntersectionType(INTERP_KERNEL::Triangulation);

  // Read the fields from the medfile

  BFT_MALLOC(r->source_fields, n_fields, MEDCouplingFieldDouble *);
  for (int ii = 0; ii < n_fields; ii++) {
    r->source_fields[ii] = ReadFieldCell(medfile_path,
                                         mesh_name,
                                         0,
                                         field_names[ii],
                                         iteration,
                                         iteration_order);
  }

  // REduced file mesh: to improve the interpolation performance,
  //                    we use a reduced mesh based only on the cells which
  //                    are intersected by the local mesh bounding box

  r->bbox_source_mesh
    = dynamic_cast<MEDCouplingUMesh *>(r->source_fields[0]->getMesh());

  return r;
}

/*----------------------------------------------------------------------------
 * Add a new remapper.
 *
 * parameters:
 *   name            <-- new remapper name
 *   elt_dim         <-- element dimension
 *   select_criteria <-- selection criteria
 *   medfile_path    <-- path of associated MED file
 *   mesh_name       <-- mesh name
 *   n_fields        <-- number of fields
 *   field_names     <-- associated field names
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/

static void
_cs_paramedmem_add_remapper(const char   *name,
                            int           elt_dim,
                            const char   *select_criteria,
                            const char   *medfile_path,
                            const char   *mesh_name,
                            int           n_fields,
                            const char  **field_names,
                            int           iteration,
                            int           iteration_order)
{
  // Allocate or reallocate if needed

  if (_remapper == NULL)
    BFT_MALLOC(_remapper, 1, cs_medcoupling_remapper_t *);
  else
    BFT_REALLOC(_remapper, _n_remappers+1, cs_medcoupling_remapper_t *);

  // Initialize new remapper, and update number of remappers

  _remapper[_n_remappers] = _cs_paramedmem_create_remapper(name,
                                                           elt_dim,
                                                           select_criteria,
                                                           medfile_path,
                                                           mesh_name,
                                                           n_fields,
                                                           field_names,
                                                           iteration,
                                                           iteration_order);

  _n_remappers++;
}

#endif /* HAVE_MEDCOUPLING_LOADER */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new ParaMEDMEM coupling.
 *
 * arguments:
 *   name     <-- name of coupling
 *   send_dec <-- send Data Exchange Channel
 *   recv_dec <-- receive Data Exchange Channel
 *----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_create(const char       *name,
                     InterpKernelDEC  *send_dec,
                     InterpKernelDEC  *recv_dec)
{
  cs_paramedmem_coupling_t  *c = NULL;

  /* Add corresponding coupling to temporary ICoCo couplings array */

  BFT_MALLOC(c, 1, cs_paramedmem_coupling_t);

  BFT_MALLOC(c->name, strlen(name) + 1, char);
  strcpy(c->name, name);

  c->n_meshes = 0;
  c->meshes = NULL;

  c->n_fields = 0;
  c->fields = NULL;

  c->send_dec = send_dec;
  c->recv_dec = recv_dec;

  return c;
}

/*----------------------------------------------------------------------------
 * Create a paramedmem coupling based on an InterpKernelDEC.
 *
 * The latter is created using the the lists of ranks provided as
 * input to this function.
 *
 * parameters:
 *   name              <-- coupling name
 *   grp1_global_ranks <-- array of ranks of group 1
 *   grp1_size         <-- size of grp1_global_ranks array
 *   grp2_global_ranks <-- array of ranks of group 2
 *   grp2_size         <-- size of grp2_global_ranks array
 *
 * return:
 *   pointer to new coupling object
 *----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_interpkernel_create(const char  *name,
                                  int         *grp1_global_ranks,
                                  int          grp1_size,
                                  int         *grp2_global_ranks,
                                  int          grp2_size)
{
  cs_paramedmem_coupling_t *c = NULL;

  /* Add corresponding coupling to temporary ICoCo couplings array */

  BFT_MALLOC(c, 1, cs_paramedmem_coupling_t);

  BFT_MALLOC(c->name, strlen(name) + 1, char);
  strcpy(c->name, name);

  c->n_meshes = 0;
  c->meshes = NULL;

  c->n_fields = 0;
  c->fields = NULL;

  bool is_in_grp1 = false;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  for (int ii = 0; ii < grp1_size; ii++) {
    if (my_rank == grp1_global_ranks[ii]) {
      is_in_grp1 = true;
      break;
    }
  }

  if (is_in_grp1) {
    c->send_dec = _cs_paramedmem_create_InterpKernelDEC(grp1_global_ranks,
                                                        grp1_size,
                                                        grp2_global_ranks,
                                                        grp2_size);

    c->recv_dec = _cs_paramedmem_create_InterpKernelDEC(grp2_global_ranks,
                                                        grp2_size,
                                                        grp1_global_ranks,
                                                        grp1_size);
  } else {
    c->recv_dec = _cs_paramedmem_create_InterpKernelDEC(grp1_global_ranks,
                                                        grp1_size,
                                                        grp2_global_ranks,
                                                        grp2_size);

    c->send_dec = _cs_paramedmem_create_InterpKernelDEC(grp2_global_ranks,
                                                        grp2_size,
                                                        grp1_global_ranks,
                                                        grp1_size);
  }

  return c;
}

/*----------------------------------------------------------------------------
 * Define new ParaMEDMEM coupling.
 *
 * arguments:
 *   name     <-- name of coupling
 *   send_dec <-- send Data Exchange Channel
 *   recv_dec <-- receive Data Exchange Channel
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_destroy(cs_paramedmem_coupling_t  **coupling)
{
  cs_paramedmem_coupling_t  *c = *coupling;

  if (c != NULL) {

    BFT_FREE(c->name);

    for (int i = 0; i < c->n_fields; i++) {
      if (c->fields[i]->pf != NULL)
        delete c->fields[i]->pf;
      c->fields[i]->f = NULL;
      BFT_FREE(c->fields[i]);
    }
    BFT_FREE(c->fields);

    for (int i = 0; i < c->n_meshes; i++)
      _destroy_mesh(&(c->meshes[i]));

    c->n_meshes = 0;
    c->meshes = NULL;

    c->send_dec = NULL;
    c->recv_dec = NULL;

  }
}

/*----------------------------------------------------------------------------
 * Define mesh for ParaMEDMEM coupling from selection criteria.
 *
 * parameters:
 *   coupling        <-- partially initialized ParaMEDMEM coupling structure
 *   name            <-- name of coupling mesh
 *   select_criteria <-- selection criteria
 *   elt_dim         <-- element dimension
 *   is_source       <-- true if fields located on mesh are sent
 *   is_dest         <-- true if fields located on mesh are received
 *
 * returns:
 *   id of created mesh in coupling
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_define_mesh(cs_paramedmem_coupling_t  *coupling,
                          const char                *name,
                          const char                *select_criteria,
                          int                        elt_dim,
                          bool                       is_source,
                          bool                       is_dest)
{
  int id;

  _paramedmem_mesh_t *mesh = NULL;

  assert(coupling != NULL);

  /* Initialization */

  BFT_MALLOC(mesh, 1, _paramedmem_mesh_t);

  BFT_MALLOC(mesh->sel_criteria, strlen(select_criteria) + 1, char);
  strcpy(mesh->sel_criteria, select_criteria);

  mesh->direction = 0;

  if (is_source)
    mesh->direction += 1;
  if (is_dest)
    mesh->direction += 2;

  mesh->elt_dim = elt_dim;

  mesh->n_elts = 0;
  mesh->elt_list = NULL;

  /* Define MED mesh (connectivity will be defined later) */

  mesh->med_mesh = MEDCouplingUMesh::New();
  mesh->med_mesh->setName(name);
  mesh->med_mesh->setTimeUnit("s");
  mesh->med_mesh->setMeshDimension(elt_dim);

  mesh->para_mesh[0] = NULL;
  mesh->para_mesh[1] = NULL;

  mesh->new_to_old = NULL;

  /* Add as new MEDCoupling mesh structure */

  id = coupling->n_meshes;
  coupling->n_meshes += 1;

  BFT_REALLOC(coupling->meshes, coupling->n_meshes, _paramedmem_mesh_t *);

  coupling->meshes[id] = mesh;

  return id;
}

/*----------------------------------------------------------------------------
 * Initialize nodal coupled meshes.
 *
 * parameters:
 *   coupling  <-- partially initialized ParaMEDMEM coupling structure
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_init_meshes(cs_paramedmem_coupling_t  *coupling)
{
  for (int i = 0; i < coupling->n_meshes; i++)
    _init_mesh_coupling(coupling, coupling->meshes[i]);
}

/*----------------------------------------------------------------------------
 * Return the ParaMEDMEM mesh id associated with a given mesh name,
 * or -1 if no association found.
 *
 * parameters:
 *   coupling  <-- coupling structure
 *   mesh_name <-- mesh name
 *
 * returns:
 *    mesh id for this coupling, or -1 if mesh name is not associated
 *    with this coupling.
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_mesh_id(cs_paramedmem_coupling_t  *coupling,
                      const char                *mesh_name)
{
  int i;
  int retval = -1;

  assert(coupling != NULL);

  for (i = 0; i < coupling->n_meshes; i++) {
    if (   strcmp(mesh_name, coupling->meshes[i]->med_mesh->getName().c_str())
        == 0)
      break;
  }

  if (i < coupling->n_meshes)
    retval = i;

  return retval;
}

/*----------------------------------------------------------------------------
 * Get number of associated coupled elements in coupled mesh
 *
 * parameters:
 *   coupling <-- ParaMEDMEM coupling structure
 *   mesh_id  <-- id of coupled mesh in coupling
 *
 * returns:
 *   number of elements in coupled mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_paramedmem_mesh_get_n_elts(const cs_paramedmem_coupling_t *coupling,
                              int                             mesh_id)
{
  cs_lnum_t retval = 0;

  if (mesh_id >= 0)
    retval = coupling->meshes[mesh_id]->n_elts;

  return retval;
}

/*----------------------------------------------------------------------------
 * Get local list of coupled elements (0 to n-1 numbering) for a coupled mesh
 *
 * parameters:
 *   coupling <-- ParaMEDMEM coupling structure
 *   mesh_id  <-- id of coupled mesh in coupling
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_paramedmem_mesh_get_elt_list(const cs_paramedmem_coupling_t *coupling,
                                int                             mesh_id)
{
  const cs_lnum_t *retval = NULL;

  if (mesh_id >= 0)
    retval = coupling->meshes[mesh_id]->elt_list;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a MEDCoupling field structure.
 *
 * parameters:
 *   coupling  <-- MED coupling structure.
 *   name      <-- field name.
 *   mesh_id   <-- id of associated mesh in structure.
 *   dim       <-- number of field components.
 *   type      <-- mesh mesh (ON_NODES, ON_CELLS)
 *   td        <-- time discretization type
 *   dirflag   <-- 1: send, 2: receive
 *
 * returns
 *   field id in coupling structure
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_field_add(cs_paramedmem_coupling_t  *coupling,
                        const char                *name,
                        int                        mesh_id,
                        int                        dim,
                        TypeOfField                type,
                        TypeOfTimeDiscretization   td,
                        int                        dirflag)
{
  int f_id = -1;
  _paramedmem_mesh_t *mesh = coupling->meshes[mesh_id];

  /* Prepare coupling structure */

  f_id = coupling->n_fields;

  BFT_REALLOC(coupling->fields,
              coupling->n_fields + 1,
              _paramedmem_field_t *);

  BFT_MALLOC(coupling->fields[f_id], 1, _paramedmem_field_t);

  /* Build ParaFIELD object if required */

  MEDCouplingFieldDouble  *f = NULL;

  if (dirflag == 1 && coupling->send_dec != NULL) {
    if (mesh->para_mesh[0] == NULL) {
      mesh->para_mesh[0] = new ParaMESH(mesh->med_mesh,
                                        *(coupling->send_dec->getSourceGrp()),
                                        "source mesh");
    }
    ComponentTopology comp_topo(dim);
    coupling->fields[f_id]->pf = new ParaFIELD(type,
                                               td,
                                               mesh->para_mesh[0],
                                               comp_topo);
    f = coupling->fields[f_id]->pf->getField();
    coupling->send_dec->attachLocalField(coupling->fields[f_id]->pf);
  }
  else if (dirflag == 2 && coupling->recv_dec != NULL) {
    if (mesh->para_mesh[1] == NULL) {
      mesh->para_mesh[1] = new ParaMESH(mesh->med_mesh,
                                        *(coupling->recv_dec->getTargetGrp()),
                                        "target mesh");
    }
    ComponentTopology comp_topo(dim);
    coupling->fields[f_id]->pf = new ParaFIELD(type,
                                               td,
                                               mesh->para_mesh[1],
                                               comp_topo);

    f = coupling->fields[f_id]->pf->getField();
    coupling->recv_dec->attachLocalField(coupling->fields[f_id]->pf);
  }
  else {
    f = MEDCouplingFieldDouble::New(type, td);
  }

  coupling->fields[f_id]->td = td;
  coupling->fields[f_id]->mesh_id = mesh_id;

  /* TODO: setNature should be set by caller to allow for more options */

  f->setNature(IntensiveConservation);

  f->setName(name);

  /* Assign array to field (filled later) */

  int n_locs = 0;
  DataArrayDouble *array = DataArrayDouble::New();

  if (type == ON_NODES)
    n_locs = mesh->med_mesh->getNumberOfNodes();
  else if (type == ON_CELLS)
    n_locs = mesh->med_mesh->getNumberOfCells();

  array->alloc(n_locs, dim);
  f->setArray(array);
  f->getArray()->decrRef();

  /* Update coupling structure */

  coupling->fields[f_id]->td = td;
  coupling->fields[f_id]->dim = dim;

  coupling->fields[f_id]->f = f;

  coupling->n_fields++;

  return f_id;
}

/*----------------------------------------------------------------------------
 * Return the ParaMEDMEM field id associated with given mesh and field names,
 * or -1 if no association found.
 *
 * parameters:
 *   coupling <-- coupling structure.
 *   mesh_id  <-- id of associated mesh in structure.
 *   name     <-- field name.
 *
 * returns
 *   field id in coupling structure, or -1 if not found
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_field_get_id(cs_paramedmem_coupling_t  *coupling,
                           int                        mesh_id,
                           const char                *name)
{
  /* Loop on fields to know if field has already been created */

  for (int f_id = 0; f_id < coupling->n_fields; f_id++) {
    if (   coupling->fields[f_id]->mesh_id == mesh_id
        && strcmp(name, coupling->fields[f_id]->f->getName().c_str()) == 0)
      return f_id;
  }

  return -1;
}

/*----------------------------------------------------------------------------
 * Return ParaMEDMEM::ParaFIELD object associated with a given field id.
 *
 * parameters:
 *   coupling  <-- pointer to associated coupling
 *   field_id  <-- id of associated field structure
 *
 * returns:
 *   pointer to ParaFIELD to which values were assigned
 *----------------------------------------------------------------------------*/

MEDCoupling::ParaFIELD *
cs_paramedmem_field_get(cs_paramedmem_coupling_t  *coupling,
                        int                        field_id)
{
  ParaFIELD *pf = NULL;

  if (field_id >= 0)
    pf = coupling->fields[field_id]->pf;

  return pf;
}

/*----------------------------------------------------------------------------
 * Write field associated with a mesh to MEDCoupling.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   coupling     <-- pointer to associated coupling
 *   field_id     <-- id of associated field
 *   on_parent    <-- if true, values are defined on parent mesh
 *   field_values <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_field_export(cs_paramedmem_coupling_t  *coupling,
                           int                        field_id,
                           bool                       on_parent,
                           const double               field_values[])
{
  int mesh_id = coupling->fields[field_id]->mesh_id;
  _paramedmem_mesh_t *mesh = coupling->meshes[mesh_id];

  MEDCouplingFieldDouble *f = NULL;

  f = coupling->fields[field_id]->f;

  double  *val_ptr = f->getArray()->getPointer();
  const int dim = coupling->fields[field_id]->dim;

  /* Assign element values */
  /*-----------------------*/

  if (! on_parent) {
    for (cs_lnum_t i = 0; i < dim*mesh->n_elts; i++)
      val_ptr[i] = field_values[i];
  }
  else {
    for (cs_lnum_t i = 0; i < mesh->n_elts; i++) {
      for (int j = 0; j < dim; j++)
        val_ptr[i*dim + j] = field_values[mesh->elt_list[i]*dim + j];
    }
  }

  /* Update field status */
  /*---------------------*/

  f->getArray()->declareAsNew();
}

/*----------------------------------------------------------------------------
 * Read field associated with a mesh from MEDCoupling.
 *
 * Only double precision floating point values are considered.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   coupling     <-- pointer to associated coupling
 *   field_id     <-- id of associated field
 *   on_parent    <-- if true, values are defined on parent mesh
 *   field_values <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_field_import(cs_paramedmem_coupling_t  *coupling,
                           int                        field_id,
                           bool                       on_parent,
                           double                     field_values[])
{
  int mesh_id = coupling->fields[field_id]->mesh_id;
  _paramedmem_mesh_t *mesh = coupling->meshes[mesh_id];

  MEDCouplingFieldDouble *f = coupling->fields[field_id]->f;

  const double  *val_ptr = f->getArray()->getConstPointer();
  const int dim = coupling->fields[field_id]->dim;

  /* Import element values */
  /*-----------------------*/

  if (! on_parent) {
    for (cs_lnum_t i = 0; i < dim*mesh->n_elts; i++)
      field_values[i] = val_ptr[i];
  }
  else {
    for (cs_lnum_t i = 0; i < mesh->n_elts; i++) {
      for (int j = 0; j < dim; j++)
        field_values[mesh->elt_list[i]*dim + j] = val_ptr[i*dim + j];
    }
  }
}

/*----------------------------------------------------------------------------
 * Synchronize DEC assciated with a given coupling.
 *
 * This sync function needs to be called at least once before exchanging data.
 * dec->synchronize() creates the interpolation matrix between the two codes!
 *
 * parameters:
 *   coupling    <-- coupling structure.
 *   dec_to_sync <-- 1 for send_dec, != 1 for recv_dec
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_sync_dec(cs_paramedmem_coupling_t  *coupling,
                       int                        dec_to_sync)
{
  if (dec_to_sync == 1) {
    coupling->send_dec->synchronize();
  } else {
    coupling->recv_dec->synchronize();
  }
}

/*----------------------------------------------------------------------------
 * Send the values related to a coupling
 *
 * parameters:
 *   coupling <-> coupling structure.
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_send_data(cs_paramedmem_coupling_t  *coupling)
{
  coupling->send_dec->sendData();
}

/*----------------------------------------------------------------------------
 * Receive the values related to a coupling
 *
 * parameters:
 *   coupling <-> coupling structure.
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_recv_data(cs_paramedmem_coupling_t  *coupling)
{
  coupling->recv_dec->recvData();
}

/*----------------------------------------------------------------------------
 * Link a given field to the DEC before send/recv
 *
 * parameters:
 *   coupling <-> coupling structure.
 *   field_id <-> associated field id
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_reattach_field(cs_paramedmem_coupling_t  *coupling,
                             int                        field_id)
{
  int mesh_id = coupling->fields[field_id]->mesh_id;
  _paramedmem_mesh_t *mesh = coupling->meshes[mesh_id];

  if (mesh->direction == 1)
    coupling->send_dec->attachLocalField(coupling->fields[field_id]->pf);
  else if (mesh->direction == 2)
    coupling->recv_dec->attachLocalField(coupling->fields[field_id]->pf);
}

/*============================================================================
 * Public C functions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Map MPI ranks within cs_glob_mpi_comm to their values in MPI_COMM_WORLD.
 *
 * The caller is responsible for freeing the returned array
 *
 * return:
 *   list of ranks in MPI_COMM_WORLD
 *----------------------------------------------------------------------------*/

int *
cs_paramedmem_get_mpi_comm_world_ranks(void)
{
  /* Global rank of current rank */
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Size of the local communicator */
  int mycomm_size;
  MPI_Comm_size(cs_glob_mpi_comm, &mycomm_size);

  int *world_ranks;
  BFT_MALLOC(world_ranks, mycomm_size, int);

  MPI_Allgather(&my_rank, 1, MPI_INT, world_ranks, 1, MPI_INT, cs_glob_mpi_comm);

  return world_ranks;
}

/*----------------------------------------------------------------------------
 * Return remapper associated with a given id
 *
 * parameters:
 *   id <-- remapper id
 *
 * return:
 *   pointer to remapper
 *----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_id(int  r_id)
{
  cs_medcoupling_remapper_t *r = _remapper[r_id];

  return r;
}

/*----------------------------------------------------------------------------
 * Return remapper associated with a given name
 *
 * parameters:
 *   name <-- remapper name
 *
 * return:
 *   pointer to remapper, or NULL
 *----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_name_try(const char  *name)
{
  if (_n_remappers > 0) {
    for (int r_id = 0; r_id < _n_remappers; r_id++) {
      const char *r_name = _remapper[r_id]->name;
      if (strcmp(r_name, name) == 0) {
        return _remapper[r_id];

      }
    }
  }

  return NULL;
}

#if defined(HAVE_MEDCOUPLING_LOADER)

/*----------------------------------------------------------------------------
 * Create or update update the list of remappers in the case where
 * several remappers may be needed.
 *
 * parameters:
 *   name            <-- new remapper name
 *   elt_dim         <-- element dimension
 *   select_criteria <-- selection criteria
 *   medfile_path    <-- path of associated MED file
 *   mesh_name       <-- mesh name
 *   n_fields        <-- number of fields
 *   field_names     <-- associated field names
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *
 * return:
 *   id of the newly added remapper within the list
 *----------------------------------------------------------------------------*/

int
cs_medcoupling_remapper_initialize(const char   *name,
                                   int           elt_dim,
                                   const char   *select_criteria,
                                   const char   *medfile_path,
                                   const char   *mesh_name,
                                   int           n_fields,
                                   const char  **field_names,
                                   int           iteration,
                                   int           iteration_order)
{
  _cs_paramedmem_add_remapper(name,
                              elt_dim,
                              select_criteria,
                              medfile_path,
                              mesh_name,
                              n_fields,
                              field_names,
                              iteration,
                              iteration_order);

  int r_id = _n_remappers - 1;

  return r_id;
}

/*----------------------------------------------------------------------------
 * Update field values (if several time steps are available in the MED file).
 *
 * parameters:
 *   r               <-- remapper object
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_set_iteration(cs_medcoupling_remapper_t  *r,
                                      int                         iteration,
                                      int                         iteration_order)
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i] = ReadFieldCell(r->medfile_path,
                                        r->mesh_name,
                                        0,
                                        r->field_names[i],
                                        iteration,
                                        iteration_order);
  }
}

#endif /* HAVE_MEDCOUPLING_LOADER */

/*----------------------------------------------------------------------------
 * Create the interpolation matrix.
 *
 * This step is separated from the interpolation step since it only needs
 * to be done once per mesh, while interpolation can be done for several
 * fields.
 *
 * parameters:
 *   r               <-- remapper object
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_setup(cs_medcoupling_remapper_t  *r)
{
  int n_elts = r->target_mesh->n_elts;

  if (n_elts > 0) {
    // List of subcells intersecting the local mesh bounding box
    const DataArrayInt *subcells
      = r->bbox_source_mesh->getCellsInBoundingBox(r->target_mesh->bbox,
                                                   1.1);

    // Construction of a subfield and the submesh associated with it.
    MEDCouplingFieldDouble *source_field
      = r->source_fields[0]->buildSubPart(subcells);

    source_field->setNature(IntensiveMaximum); /* TODO options */

    // Update the remapper structure and interpolation matrix
    // TODO allow settings for precision and interpolation type
    r->remapper->setPrecision(1.e-12);
    r->remapper->setIntersectionType(INTERP_KERNEL::Triangulation);

    r->remapper->prepare(source_field->getMesh(),
                         r->target_mesh->med_mesh,
                         "P0P0");
  }
}

/*----------------------------------------------------------------------------
 * Copy interpolated values to a new array.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   field_id        <-- id of given field
 *   r               <-- pointer to remapper object
 *   default_val     <-- default value
 *
 * return:
 *   pointer to allocated values array
 *----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_remapper_copy_values(cs_medcoupling_remapper_t  *r,
                                    int                         field_id,
                                    double                      default_val)
{
  int n_elts = r->target_mesh->n_elts;
  int n_elts_loc = cs_glob_mesh->n_cells;

  cs_real_t *new_vals;
  BFT_MALLOC(new_vals, n_elts_loc, cs_real_t);
  for (int i = 0; i < n_elts_loc; i++) {
    new_vals[i] = default_val;
  }

  if (n_elts > 0) {
    // List of subcells intersecting the local mesh bounding box
    const DataArrayInt *subcells
      = r->bbox_source_mesh->getCellsInBoundingBox(r->target_mesh->bbox,
                                                   1.1);
    // Construct the subfields based on the subcells list
    MEDCouplingFieldDouble *source_field
      = r->source_fields[field_id]->buildSubPart(subcells);

    // Set the nature of the field
    source_field->setNature(IntensiveMaximum); /* TODO options */

    // Interpolate the new values
    MEDCouplingFieldDouble *target_field
      = r->remapper->transferField(source_field, default_val);

    // Generate the output array
    const double *val_ptr = target_field->getArray()->getConstPointer();
    int npts = target_field->getNumberOfValues();

    if (r->target_mesh->elt_list != NULL) {
      for (int i = 0; i < npts; i++) {
        int e_id = r->target_mesh->new_to_old[i];
        new_vals[e_id] = val_ptr[i];
      }
    } else {
      for (int i = 0; i < npts; i++) {
        new_vals[i] = val_ptr[i];
      }
    }

  }

  return new_vals;
}

/*----------------------------------------------------------------------------
 * Translate the mapped source mesh.
 *
 * Caution: cs_medcoupling_remapper_prepare() must to be called after this
 * function in order to update the interpolation matrix.
 *
 * parameters:
 *   r           <-- pointer to remapper object
 *   translation <-- translation vector
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_translate(cs_medcoupling_remapper_t  *r,
                                  cs_real_t                   translation[3])
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i]->getMesh()->translate(translation);
  }
}

/*----------------------------------------------------------------------------
 * Rotate the mapped source mesh.
 *
 * Caution: cs_medcoupling_remapper_prepare() must to be called after this
 * function in order to update the interpolation matrix.
 *
 * parameters:
 *   r         <-- pointer to remapper object
 *   invariant <-- coordinates of invariant point
 *   axis      <-- rotation axis vector
 *   angle     <-- rotation angle
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_rotate(cs_medcoupling_remapper_t  *r,
                               cs_real_t                   invariant[3],
                               cs_real_t                   axis[3],
                               cs_real_t                   angle)
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i]->getMesh()->rotate(invariant, axis, angle);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*----------------------------------------------------------------------------*/

#endif /* HAVE_MEDCOUPLING */
