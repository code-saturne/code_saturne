/*============================================================================
 * ParaMEDMEM coupling
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
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#include <MEDCoupling_version.h>

#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>

#if defined(HAVE_PARAMEDMEM)
#include <ParaFIELD.hxx>
#include <ParaMESH.hxx>
#include <InterpKernelDEC.hxx>

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

#include "cs_medcoupling_utils.hxx"
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
 * ParaMEDMED mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_medcoupling_mesh_t *mesh;        /* The Code_Saturne stracture englobing
                                         the medcoupling mesh structure */

  char               *sel_criteria;   /* Element selection criteria */

  int                 elt_dim;        /* Element dimension */

  cs_lnum_t           n_elts;         /* Number of coupled elements */
  cs_lnum_t          *elt_list;       /* List of associated elements
                                        (0 to n-1) */
  int                *new_to_old;     /* Connectivity used if only a section of
                                        the mesh is read */

  MEDCouplingUMesh   *med_mesh;       /* MED mesh structure */

  cs_real_t          *bbox;           /* Bounding box to optimize search */

  int                 direction;      /* 1: send, 2: receive, 3: both */

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

/*----------------------------------------------------------------------------
 * MEDCoupling parallel interpolation structure
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

  /* Building the MED representation of the internal mesh */
  cs_medcoupling_mesh_copy_from_base(parent_mesh, mesh->mesh);

  /* Linking the pointers in order to simplify the calls*/
  mesh->sel_criteria = mesh->mesh->sel_criteria;
  mesh->elt_dim      = mesh->mesh->elt_dim;
  mesh->n_elts       = mesh->mesh->n_elts;
  mesh->elt_list     = mesh->mesh->elt_list;
  mesh->new_to_old   = mesh->mesh->new_to_old;
  mesh->med_mesh     = mesh->mesh->med_mesh;
  mesh->bbox         = mesh->mesh->bbox;

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
_destroy_mesh(_paramedmem_mesh_t  **mesh)
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

  _paramedmem_mesh_t *pmmesh = NULL;
  cs_medcoupling_mesh_t *mesh = NULL;

  assert(coupling != NULL);

  /* Initialization */

  BFT_MALLOC(pmmesh, 1, _paramedmem_mesh_t);
  BFT_MALLOC(mesh, 1, cs_medcoupling_mesh_t);

  BFT_MALLOC(mesh->sel_criteria, strlen(select_criteria) + 1, char);
  strcpy(mesh->sel_criteria, select_criteria);

  pmmesh->direction = 0;

  if (is_source)
    pmmesh->direction += 1;
  if (is_dest)
    pmmesh->direction += 2;

  mesh->elt_dim = elt_dim;

  mesh->n_elts = 0;
  mesh->elt_list = NULL;

  /* Define MED mesh (connectivity will be defined later) */

  mesh->med_mesh = MEDCouplingUMesh::New();
  mesh->med_mesh->setName(name);
  mesh->med_mesh->setTimeUnit("s");
  mesh->med_mesh->setMeshDimension(elt_dim);

  pmmesh->para_mesh[0] = NULL;
  pmmesh->para_mesh[1] = NULL;

  mesh->new_to_old = NULL;

  /* Add as new MEDCoupling mesh structure */

  id = coupling->n_meshes;
  coupling->n_meshes += 1;

  BFT_REALLOC(coupling->meshes, coupling->n_meshes, _paramedmem_mesh_t *);

  pmmesh->mesh = mesh;
  coupling->meshes[id] = pmmesh;

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
    if (strcmp(mesh_name, coupling->meshes[i]->med_mesh->getName().c_str())
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
cs_paramedmem_mesh_get_n_elts(const cs_paramedmem_coupling_t  *coupling,
                              int                              mesh_id)
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
cs_paramedmem_mesh_get_elt_list(const cs_paramedmem_coupling_t  *coupling,
                                int                              mesh_id)
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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* HAVE_PARAMEDMEM */
