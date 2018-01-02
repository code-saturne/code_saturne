/*============================================================================
 * ICoCo coupling
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
 * MED library headers
 *----------------------------------------------------------------------------*/

#include <mpi.h>
#include <ParaFIELD.hxx>
#include <ParaMESH.hxx>
#include <InterpKernelDEC.hxx>

#include <MEDCoupling_version.h>

#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>

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
 * MED field structure
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
 * MED mesh structure
 *----------------------------------------------------------------------------*/

/* Structure associated with ParaMEDMEM entity */

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

/*============================================================================
 *  Global variables
 *============================================================================*/

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
        for (j = 0; j < dim; j++)
          med_coords->setIJ(vtx_id[i], j, vertex_coords[i*dim + j]);
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

  cs_lnum_t vtx_count = 0;
  int elt_buf_size = 4;
  int *elt_buf = NULL;
  cs_lnum_t *vtx_id = NULL;

  const int perm_tri[4] = {0, 2, 1};
  const int perm_quad[4] = {0, 3, 2, 1};

  /* Mark and renumber vertices */

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);

  for (i = 0; i < mesh->n_vertices; i++)
    vtx_id[i] = -1;

  /* Case with filter list */

  if (elts_list != NULL) {

    for (i = 0; i < n_elts; i++) {
      cs_lnum_t eid = elts_list[i] - 1;
      for (j = mesh->b_face_vtx_idx[eid];
           j < mesh->b_face_vtx_idx[eid+1];
           j++) {
        cs_lnum_t vid = mesh->b_face_vtx_lst[j];
        if (vtx_id[vid] < 0)
          vtx_id[vid] = vtx_count++;
      }
    }

  }
  else {

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

    cs_lnum_t eid = (elts_list != NULL) ? elts_list[i] - 1 : i;

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
                  MEDCouplingUMesh  *med_mesh)
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
    for (i = 0; i < n_elts; i++)
      cell_id[elts_list[i]-1] = cell_count++;
  }
  else {
    for (i = 0; i < n_elts; i++)
      cell_id[elts_list[i]-1] = i;
  }

  /* Mark and renumber vertices */

  BFT_MALLOC(vtx_id, mesh->n_vertices, cs_lnum_t);

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
    c_id = cell_id[mesh->b_face_cells[face_id] - 1];
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

  const cs_lnum_t  face_num_shift[3] = {0,
                                        mesh->n_b_faces,
                                        mesh->n_b_faces + mesh->n_i_faces};
  const cs_lnum_t  *face_vertices_idx[2] = {mesh->b_face_vtx_idx,
                                            mesh->i_face_vtx_idx};
  const cs_lnum_t  *face_vertices_num[2] = {mesh->b_face_vtx_lst,
                                            mesh->i_face_vtx_lst};

  BFT_MALLOC(elt_buf, elt_buf_size, int);

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
      elt_buf[0] = vtx_id[0];
      elt_buf[1] = vtx_id[2];
      elt_buf[2] = vtx_id[1];
      elt_buf[3] = vtx_id[3];
      break;

    case FVM_CELL_PYRAM:
      type = INTERP_KERNEL::NORM_PYRA5;
      n_vtx = 5;
      elt_buf[0] = vtx_id[0];
      elt_buf[1] = vtx_id[3];
      elt_buf[2] = vtx_id[2];
      elt_buf[3] = vtx_id[1];
      elt_buf[4] = vtx_id[4];
      break;

    case FVM_CELL_PRISM:
      type = INTERP_KERNEL::NORM_PENTA6;
      n_vtx = 6;
      elt_buf[0] = vtx_id[0];
      elt_buf[1] = vtx_id[2];
      elt_buf[2] = vtx_id[1];
      elt_buf[3] = vtx_id[3];
      elt_buf[4] = vtx_id[5];
      elt_buf[5] = vtx_id[4];
      break;

    case FVM_CELL_HEXA:
      type = INTERP_KERNEL::NORM_HEXA8;
      n_vtx = 8;
      elt_buf[0] = vtx_id[0];
      elt_buf[1] = vtx_id[3];
      elt_buf[2] = vtx_id[2];
      elt_buf[3] = vtx_id[1];
      elt_buf[4] = vtx_id[4];
      elt_buf[5] = vtx_id[7];
      elt_buf[6] = vtx_id[6];
      elt_buf[7] = vtx_id[5];
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
 * Initialize mesh for ParaMEDMEM coupling.
 *
 * parameters:
 *   coupling  <-- coupling structure.
 *   mesh      <-> partially ParaMEDMEM mesh coupling structure
 *----------------------------------------------------------------------------*/

static void
_init_mesh(cs_paramedmem_coupling_t  *coupling,
           _paramedmem_mesh_t        *mesh)
{
  cs_mesh_t *parent_mesh = cs_glob_mesh;

  assert(mesh != NULL);

  /* Creation of a new nodal mesh from selected cells */

  if (mesh->elt_dim == 3) {

    BFT_MALLOC(mesh->elt_list, parent_mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(mesh->sel_criteria,
                              &(mesh->n_elts),
                              mesh->elt_list);

    BFT_REALLOC(mesh->elt_list, mesh->n_elts, cs_lnum_t);

    _assign_cell_mesh(parent_mesh,
                      mesh->n_elts,
                      mesh->elt_list,
                      mesh->med_mesh);

  }

  /* Creation of a new nodal mesh from selected border faces */

  else if (mesh->elt_dim == 2) {

    BFT_MALLOC(mesh->elt_list, parent_mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(mesh->sel_criteria,
                                &(mesh->n_elts),
                                mesh->elt_list);

    BFT_REALLOC(mesh->elt_list, mesh->n_elts, cs_lnum_t);

    _assign_face_mesh(parent_mesh,
                      mesh->n_elts,
                      mesh->elt_list,
                      mesh->med_mesh);

  }

  /* Define associated ParaMESH */

  if (mesh->direction & 1)
    mesh->para_mesh[0] = new ParaMESH(mesh->med_mesh,
                                      *(coupling->send_dec->getSourceGrp()),
                                      "source mesh");

  if (mesh->direction & 2)
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

cs_paramedmem_coupling_t  *
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

    BFT_FREE(c->fields);
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
 *   coupling        <-- partially initialized ParaMEDMEM coupling structure
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_init_meshes(cs_paramedmem_coupling_t  *coupling)
{
  for (int i = 0; i < coupling->n_meshes; i++)
    _init_mesh(coupling, coupling->meshes[i]);
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
  else
    f = MEDCouplingFieldDouble::New(type, td);

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
  //coupling->fields[f_id]->pf = NULL;

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

/*----------------------------------------------------------------------------*/
