/*============================================================================
 * ParaMEDMEM coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 *  PLE headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

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
#include "cs_coupling.h"

#include "fvm_defs.h"
#include "fvm_nodal_from_desc.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_medcoupling_utils.hxx"
#include "cs_paramedmem_coupling.hxx"

#if defined(HAVE_PARAMEDMEM)

#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>

#include <ParaFIELD.hxx>
#include <ParaMESH.hxx>
#include <InterpKernelDEC.hxx>

using namespace MEDCoupling;
#endif

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * ParaMEDMED field structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int                       dim;           /* Field dimension */

#if defined(HAVE_PARAMEDMEM)
  TypeOfTimeDiscretization  td;            /* NO_TIME, ONE_TIME, LINEAR_TIME,
                                              or CONST_ON_TIME_INTERVAL */
  MEDCouplingFieldDouble   *f;             /* Pointer to MED coupling field */
  ParaFIELD                *pf;            /* Pointer to ParaMEDMEM field */
#else
  void                     *td;
  void                     *f;
  void                     *pf;
#endif

} _paramedmem_field_t;

/*----------------------------------------------------------------------------
 * ParaMEDMED mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_medcoupling_mesh_t *mesh;        /* The Code_Saturne stracture englobing
                                         the medcoupling mesh structure */

#if defined(HAVE_PARAMEDMEM)
  ParaMESH           *para_mesh[2];   /* parallel MED mesh structures
                                         for send (0) and receive (1) */
#else
  void               *para_mesh[2];
#endif
} _paramedmem_mesh_t;

/*----------------------------------------------------------------------------
 * MEDCoupling writer/reader structure
 *----------------------------------------------------------------------------*/

struct _cs_paramedmem_coupling_t {

  /* Coupling Name */
  char                      *name;           /* Coupling name */

  /* Current app name */
  ple_coupling_mpi_set_info_t apps[2];
  int                         loc_app;

  _paramedmem_mesh_t        *pmesh;         /* Array of mesh helper
                                                structures */

  int                        n_fields;       /* Number of fields */
  _paramedmem_field_t      **fields;         /* Array of field helper
                                                structures */

#if defined(HAVE_PARAMEDMEM)
  InterpKernelDEC           *send_dec;       /* Send data exchange channel */
  InterpKernelDEC           *recv_dec;       /* Receive data exchange channel */
#else
  void                      *send_dec;
  void                      *recv_dec;
#endif

  int                       send_synced;
  int                       recv_synced;
};

/*=============================================================================
 * Private global variables
 *============================================================================*/

static int                          _n_paramed_couplers = 0;
static cs_paramedmem_coupling_t   **_paramed_couplers = NULL;

const int cs_medcpl_cell_field = 0;
const int cs_medcpl_vertex_field = 1;

const int cs_medcpl_no_time = 0;
const int cs_medcpl_one_time = 1;
const int cs_medcpl_linear_time = 2;

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_PARAMEDMEM)

/*----------------------------------------------------------------------------
 * Initialize mesh for ParaMEDMEM coupling.
 *
 * parameters:
 *   coupling  <-- coupling structure.
 *   mesh      <-> partially ParaMEDMEM mesh coupling structure
 *----------------------------------------------------------------------------*/

static void
_init_mesh_coupling(cs_paramedmem_coupling_t  *coupling)
{
  cs_mesh_t *parent_mesh = cs_glob_mesh;

  assert(coupling->pmesh != NULL);

  _paramedmem_mesh_t *pmesh = coupling->pmesh;

  /* Building the MED representation of the internal mesh */
  cs_medcoupling_mesh_copy_from_base(parent_mesh, pmesh->mesh, 0);

  /* Define associated ParaMESH */
  pmesh->para_mesh[0] = new ParaMESH(coupling->pmesh->mesh->med_mesh,
                                     *(coupling->send_dec->getSourceGrp()),
                                     "source mesh");
  pmesh->para_mesh[1] = new ParaMESH(coupling->pmesh->mesh->med_mesh,
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
_destroy_mesh(_paramedmem_mesh_t  *pm)
{

  if (pm == NULL)
    return;

  for (int i = 0; i < 2; i++) {
    if (pm->para_mesh[i] != NULL)
      delete pm->para_mesh[i];
  }
  if (pm->mesh != NULL)
    cs_medcoupling_mesh_destroy(pm->mesh);

  BFT_FREE(pm);
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
_create_InterpKernelDEC(const int grp1_root_rank,
                        const int grp1_n_ranks,
                        const int grp2_root_rank,
                        const int grp2_n_ranks)
{
  /* Group 1 id's */
  std::set<int> grp1_ids;
  for (int i = 0; i < grp1_n_ranks; i++) {
    grp1_ids.insert(grp1_root_rank + i);
  }

  /* Group 2 id's */
  std::set<int> grp2_ids;
  for (int i = 0; i < grp2_n_ranks; i++) {
    grp2_ids.insert(grp2_root_rank + i);
  }

  /* Create the InterpKernel DEC */
  InterpKernelDEC *NewDec = new InterpKernelDEC(grp1_ids, grp2_ids);

  return NewDec;
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

static void
_add_paramedmem_interpkernel(const char                  *cpl_name,
                             ple_coupling_mpi_set_info_t  apps[2],
                             const int                    l_id)
{

  assert(l_id > -1);

  if (_paramed_couplers == NULL)
    BFT_MALLOC(_paramed_couplers,
               1,
               cs_paramedmem_coupling_t *);
  else
    BFT_REALLOC(_paramed_couplers,
                _n_paramed_couplers + 1,
                cs_paramedmem_coupling_t *);

  cs_paramedmem_coupling_t *c = NULL;

  /* Add corresponding coupling to temporary ICoCo couplings array */

  BFT_MALLOC(c, 1, cs_paramedmem_coupling_t);

  /* Apps identification */
  for (int i = 0; i < 2; i++)
    c->apps[i] = apps[i];
  c->loc_app = l_id;

  /* Set coupling name */
  BFT_MALLOC(c->name, strlen(cpl_name) + 1, char);
  strcpy(c->name, cpl_name);

  c->pmesh = NULL;

  c->n_fields = 0;
  c->fields = NULL;

  c->send_synced = 0;
  c->recv_synced = 0;

  InterpKernelDEC *dec1 = _create_InterpKernelDEC(apps[0].root_rank,
                                                  apps[0].n_ranks,
                                                  apps[1].root_rank,
                                                  apps[1].n_ranks);

  InterpKernelDEC *dec2 = _create_InterpKernelDEC(apps[1].root_rank,
                                                  apps[1].n_ranks,
                                                  apps[0].root_rank,
                                                  apps[0].n_ranks);
  if (l_id == 0) {
    c->send_dec = dec1;
    c->recv_dec = dec2;
  } else {
    c->send_dec = dec2;
    c->recv_dec = dec1;
  }

  _paramed_couplers[_n_paramed_couplers] = c;

  _n_paramed_couplers++;

}

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_by_id(int  pc_id)
{

  cs_paramedmem_coupling_t * c = _paramed_couplers[pc_id];

  return c;

}

#endif /* HAVE_PARAMEDMEM */

/*============================================================================
 * Public C functions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Define new ParaMEDMEM coupling.
 *
 * arguments:
 *   name               <-- name of coupling
 *   grp1_global_ranks  <-- First group ranks in MPI_COMM_WORLD
 *   grp1_size          <-- Number of ranks in the first group
 *   grp2_global_ranks  <-- Second group ranks in MPI_COMM_WORLD
 *   grp2_size          <-- Number of ranks in the second group
 *----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create(const char *app1_name,
                              const char *app2_name,
                              const char *cpl_name)
{
  cs_paramedmem_coupling_t *c = NULL;

  /* Check that at least on app name is provided */
  if (app1_name == NULL && app2_name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Partner application name was not provided.\n"));

  /* Find which app is which */
  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);
  const int app_id = ple_coupling_mpi_set_get_app_id(mpi_apps);

  /* Get each application ranks */
  ple_coupling_mpi_set_info_t apps[2];
  for (int i = 0; i < 2; i++)
    apps[i] = ple_coupling_mpi_set_get_info(mpi_apps, -1);

  int l_id = -1;
  if (app1_name == NULL)
    l_id = 0;
  else if (app2_name == NULL)
    l_id = 1;

  for (int i = 0; i < n_apps; i++) {
    ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (l_id > -1) {
      if (app_id == i)
        apps[l_id] = ai;
    }
    if (app1_name != NULL) {
      if (strcmp(ai.app_name, app1_name) == 0) {
        apps[0] = ai;
        if (app_id == i)
          l_id = i;
      }
    }
    if (app2_name != NULL) {
      if (strcmp(ai.app_name, app2_name) == 0) {
        apps[1] = ai;
        if (app_id == i)
          l_id = i;
      }
    }
  }

  for (int i = 0; i < 2; i++) {
    if (apps[i].root_rank < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error: Partner app was not found...\n"));
  }

  /* Set coupling name */
  char *_name = NULL;

  if (cpl_name == NULL) {
    /* string is <a_name>_<p_name>_cpl */
    size_t lname = strlen(apps[0].app_name)
                 + strlen(apps[1].app_name)
                 + 6;
    BFT_MALLOC(_name, lname, char);
    sprintf(_name, "%s_%s_cpl", apps[0].app_name, apps[1].app_name);

  } else {
    size_t lname = strlen(cpl_name);
    BFT_MALLOC(_name, lname, char);
    strcpy(_name, cpl_name);
  }


#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else
  _add_paramedmem_interpkernel(_name, apps, l_id);

  BFT_FREE(_name);

  c = _paramed_couplers[_n_paramed_couplers-1];
#endif

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

#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else
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

    _destroy_mesh(c->pmesh);

    c->pmesh = NULL;

    c->send_dec = NULL;
    c->recv_dec = NULL;

  }
#endif

  return;
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
                          int                        elt_dim)
{
  int id = -1;

#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else
  _paramedmem_mesh_t *pmmesh = NULL;
  cs_medcoupling_mesh_t *mesh = NULL;

  assert(coupling != NULL);

  if (coupling->pmesh != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: mesh was allready defined for this coupling.\n"));

  /* Initialization */

  BFT_MALLOC(pmmesh, 1, _paramedmem_mesh_t);
  BFT_MALLOC(mesh, 1, cs_medcoupling_mesh_t);

  BFT_MALLOC(mesh->sel_criteria, strlen(select_criteria) + 1, char);
  strcpy(mesh->sel_criteria, select_criteria);

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

  pmmesh->mesh = mesh;
  coupling->pmesh = pmmesh;
#endif

  return id;
}

/*----------------------------------------------------------------------------
 * Initialize nodal coupled meshes.
 *
 * parameters:
 *   coupling  <-- partially initialized ParaMEDMEM coupling structure
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_init_mesh(cs_paramedmem_coupling_t  *coupling)
{
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else
  _init_mesh_coupling(coupling);
#endif

  return;
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
cs_paramedmem_mesh_get_n_elts(const cs_paramedmem_coupling_t  *coupling)
{
  cs_lnum_t retval = 0;

#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else
  retval = coupling->pmesh->mesh->n_elts;
#endif

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
cs_paramedmem_mesh_get_elt_list(const cs_paramedmem_coupling_t  *coupling)
{
  const cs_lnum_t *retval = NULL;

#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else
  retval = coupling->pmesh->mesh->elt_list;
#endif

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
                        int                        dim,
                        int                        medcpl_field_type,
                        int                        medcpl_time_discr,
                        int                        dirflag)
{
  int f_id = -1;

#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  _paramedmem_mesh_t *pmesh = coupling->pmesh;

  /* Prepare coupling structure */
  TypeOfField type = ON_CELLS;
  if (medcpl_field_type == cs_medcpl_cell_field)
    type = ON_CELLS;
  else if (medcpl_field_type == cs_medcpl_vertex_field)
    type = ON_NODES;

  TypeOfTimeDiscretization td = NO_TIME;
  if (medcpl_time_discr == cs_medcpl_no_time)
    td = NO_TIME;
  else if (medcpl_time_discr == cs_medcpl_one_time)
    td = ONE_TIME;
  else if (medcpl_time_discr == cs_medcpl_linear_time)
    td = LINEAR_TIME;


  f_id = coupling->n_fields;

  BFT_REALLOC(coupling->fields,
              coupling->n_fields + 1,
              _paramedmem_field_t *);

  BFT_MALLOC(coupling->fields[f_id], 1, _paramedmem_field_t);

  /* Build ParaFIELD object if required */

  MEDCouplingFieldDouble  *f = NULL;

  if (dirflag == 1 && coupling->send_dec != NULL) {
    if (pmesh->para_mesh[0] == NULL) {
      pmesh->para_mesh[0] = new ParaMESH(pmesh->mesh->med_mesh,
                                        *(coupling->send_dec->getSourceGrp()),
                                        "source mesh");
    }
    ComponentTopology comp_topo(dim);
    coupling->fields[f_id]->pf = new ParaFIELD(type,
                                               td,
                                               pmesh->para_mesh[0],
                                               comp_topo);
    f = coupling->fields[f_id]->pf->getField();
    coupling->send_dec->attachLocalField(coupling->fields[f_id]->pf);
  }
  else if (dirflag == 2 && coupling->recv_dec != NULL) {
    if (pmesh->para_mesh[1] == NULL) {
      pmesh->para_mesh[1] = new ParaMESH(pmesh->mesh->med_mesh,
                                        *(coupling->recv_dec->getTargetGrp()),
                                        "target mesh");
    }
    ComponentTopology comp_topo(dim);
    coupling->fields[f_id]->pf = new ParaFIELD(type,
                                               td,
                                               pmesh->para_mesh[1],
                                               comp_topo);

    f = coupling->fields[f_id]->pf->getField();
    coupling->recv_dec->attachLocalField(coupling->fields[f_id]->pf);
  }
  else {
    f = MEDCouplingFieldDouble::New(type, td);
  }

  coupling->fields[f_id]->td = td;

  /* TODO: setNature should be set by caller to allow for more options */

  f->setNature(IntensiveConservation);

  f->setName(name);

  /* Assign array to field (filled later) */

  int n_locs = 0;
  DataArrayDouble *array = DataArrayDouble::New();

  if (type == ON_NODES)
    n_locs = pmesh->mesh->med_mesh->getNumberOfNodes();
  else if (type == ON_CELLS)
    n_locs = pmesh->mesh->med_mesh->getNumberOfCells();

  array->alloc(n_locs, dim);
  f->setArray(array);
  f->getArray()->decrRef();

  /* Update coupling structure */

  coupling->fields[f_id]->td = td;
  coupling->fields[f_id]->dim = dim;

  coupling->fields[f_id]->f = f;

  coupling->n_fields++;

#endif

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
                           const char                *name)
{

  int f_id = -1;
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  /* Loop on fields to know if field has already been created */
  for (int i = 0; i < coupling->n_fields; i++) {
    if (strcmp(name, coupling->fields[i]->f->getName().c_str()) == 0) {
      f_id = i;
      break;
    }
  }

#endif

  return f_id;
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

#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  _paramedmem_mesh_t *pmesh = coupling->pmesh;

  MEDCouplingFieldDouble *f = NULL;

  f = coupling->fields[field_id]->f;

  double  *val_ptr = f->getArray()->getPointer();
  const int dim = coupling->fields[field_id]->dim;

  /* Assign element values */
  /*-----------------------*/
  if (! on_parent) {
    cs_lnum_t n_elts = pmesh->mesh->n_elts;

    for (cs_lnum_t i = 0; i < dim*n_elts; i++)
      val_ptr[i] = field_values[i];
  }
  else {
    cs_lnum_t  n_elts   = pmesh->mesh->n_elts;
    cs_lnum_t *elt_list = pmesh->mesh->elt_list;
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (int j = 0; j < dim; j++)
        val_ptr[i*dim + j] = field_values[elt_list[i]*dim + j];
    }
  }

  /* Update field status */
  /*---------------------*/
  f->getArray()->declareAsNew();

#endif

  return;
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
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  _paramedmem_mesh_t *pmesh = coupling->pmesh;

  MEDCouplingFieldDouble *f = coupling->fields[field_id]->f;

  const double  *val_ptr = f->getArray()->getConstPointer();
  const int dim = coupling->fields[field_id]->dim;

  /* Import element values */
  /*-----------------------*/

  if (! on_parent) {
    cs_lnum_t  n_elts = pmesh->mesh->n_elts;
    cs_lnum_t *new_to_old = pmesh->mesh->new_to_old;

    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int j = 0; j < dim; j++) {
        cs_lnum_t c_id = new_to_old[i];
        field_values[dim*c_id+j] = val_ptr[i*dim + j];
      }
  }
  else {
    cs_lnum_t  n_elts   = pmesh->mesh->n_elts;
    cs_lnum_t *elt_list = pmesh->mesh->elt_list;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (int j = 0; j < dim; j++)
        field_values[elt_list[i]*dim + j] = val_ptr[i*dim + j];
    }
  }

#endif

  return;
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
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  if (dec_to_sync == 1) {
    if (coupling->send_synced == 0) {
      coupling->send_dec->synchronize();
      coupling->send_synced = 1;
    }
  } else if (coupling->recv_synced == 0) {
    coupling->recv_dec->synchronize();
    coupling->recv_synced = 1;
  }

#endif

  return;
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
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  coupling->send_dec->sendData();

#endif

  return;
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
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  coupling->recv_dec->recvData();

#endif

  return;
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
                             int                        field_id,
                             int                        dec_id)
{
#if !defined(HAVE_PARAMEDMEM)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling MPI support.\n"));
#else

  _paramedmem_mesh_t *pmesh = coupling->pmesh;

  if (dec_id == 0)
    coupling->send_dec->attachLocalField(coupling->fields[field_id]->pf);
  else if (dec_id == 1)
    coupling->recv_dec->attachLocalField(coupling->fields[field_id]->pf);
  else
    bft_error(__FILE__, __LINE__, 0, _("Error: Uknown dec id.\n"));

#endif

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
