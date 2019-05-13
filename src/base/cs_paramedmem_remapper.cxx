/*============================================================================
 * Parallel interpolation using ParaMEDMEM OverlapDEC
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_PARAMEDMEM)
#include <MEDLoader.hxx>
#include <MEDCoupling_version.h>
#include <MEDFileField.hxx>
#include <MEDFileFieldMultiTS.hxx>

#include <MEDFileMesh.hxx>
#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>

#include <ParaFIELD.hxx>
#include <ParaMESH.hxx>
#include <OverlapDEC.hxx>
#include <ParaMEDFileMesh.hxx>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_medcoupling_utils.hxx"
#include "cs_paramedmem_remapper.hxx"

/*----------------------------------------------------------------------------*/

using namespace MEDCoupling;

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * MEDCoupling/ParaMEDMEM parallel interpolation structure
 *----------------------------------------------------------------------------*/

struct _cs_paramedmem_remapper_t {

  char                             *name;           /* Coupling name */

  cs_medcoupling_mesh_t            *local_mesh;     /* Local cs_mesh in MED format */

  MEDFileUMesh                     *src_mesh;
  MEDFileFields                    *MEDFields;

  int                               ntsteps;
  int                              *iter;
  int                              *order;
  cs_real_t                        *time_steps;

  OverlapDEC                       *odec;           /* Overlap data exchange channel */
  int                               synced;
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Returns an array containing the ranks of Code_Saturne processes in
 *          MPI_COMM_WORLD.
 *
 * \return  array of ranks in MPI_COMM_WORLD
 */
/* -------------------------------------------------------------------------- */

int *
_cs_paramedmem_get_mpi_comm_world_ranks(void)
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

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Create the parallel remapper structure based on ParaMEDMEM OverlapDEC
 *
 * \param[in] name    name of the remapper
 *
 * \return  pointer to cs_paramedmem_remapper_t struct
 */
/* -------------------------------------------------------------------------- */

cs_paramedmem_remapper_t *
_cs_paramedmem_overlap_create(const char  *name)
{
  cs_paramedmem_remapper_t *r = NULL;

  /* Add corresponding coupling to temporary ICoCo couplings array */

  BFT_MALLOC(r, 1, cs_paramedmem_remapper_t);

  BFT_MALLOC(r->name, strlen(name) + 1, char);
  strcpy(r->name, name);

  r->iter       = NULL;
  r->order      = NULL;
  r->time_steps = NULL;
  r->ntsteps    = -1;

  r->synced = 0;

  /* Local id's */
  int *cs_ranks = _cs_paramedmem_get_mpi_comm_world_ranks();
  int cs_comm_size;
  MPI_Comm_size(cs_glob_mpi_comm, &cs_comm_size);

  std::set<int> grp_ids;
  for (int i = 0; i < cs_comm_size; i++) {
    grp_ids.insert(cs_ranks[i]);
  }

  /* Create the Overlap DEC */
  r->odec = new OverlapDEC(grp_ids);

  return r;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief  Set the target mesh for interpolation
 *
 * \param[in] r               pointer to cs_paramedmem_remapper_t struct
 * \param[in] name            mesh name
 * \param[in] select_criteria selection criteria for needed cells
 *
 */
/* -------------------------------------------------------------------------- */

void
_cs_paramedmem_remapper_target_mesh(cs_paramedmem_remapper_t  *r,
                                    const char                *name,
                                    const char                *select_criteria)
{

  assert(r != NULL);

  /* Initialization: Only volumes are possible for this option */
  r->local_mesh = cs_medcoupling_mesh_create(name,
                                             select_criteria,
                                             3);

  cs_mesh_t *parent_mesh = cs_glob_mesh;

  /* Building the MED representation of the internal mesh */
  cs_medcoupling_mesh_copy_from_base(parent_mesh, r->local_mesh, 0);

  return;

}

/*----------------------------------------------------------------------------*/

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Load the mesh parts on each process
 *
 * \param[in] r         pointer to cs_paramedmem_remapper_t struct
 * \param[in] fileName  name of med file containing the data
 * \param[in] meshName  name of the mesh to read in the med file
 *
 */
/* -------------------------------------------------------------------------- */

void
_cs_paramedmem_load_paramesh(cs_paramedmem_remapper_t *r,
                             char                     *fileName,
                             char                     *meshName)
{

  int myPart;
  MPI_Comm_rank(cs_glob_mpi_comm, &myPart);
  int nParts;
  MPI_Comm_size(cs_glob_mpi_comm, &nParts);

  const std::string fname = fileName;
  const std::string mname = meshName;

  // Mesh is stored with -1, -1 indices in MED files
  r->src_mesh = ParaMEDFileUMesh::New(myPart,
                                      nParts,
                                      fname,
                                      mname,
                                      -1,
                                      -1);

  r->src_mesh->forceComputationOfParts();

  MCAuto<MEDFileMeshes> ms = MEDFileMeshes::New();
  ms = MEDFileMeshes::New();
  ms->pushMesh(r->src_mesh);

  r->MEDFields = MEDFileFields::LoadPartOf(fname, true, ms);
  r->MEDFields->loadArrays();

  /* Get Time steps info */
  MCAuto<MEDFileAnyTypeFieldMultiTS> fts = r->MEDFields->getFieldAtPos(0);

  bft_printf("Got field!\n");
  bft_printf_flush();

  std::vector< std::pair<int,int> > tio = fts->getIterations();
  bft_printf("Got iterations!\n");
  bft_printf_flush();

  r->ntsteps    = tio.size();
  bft_printf("Got steps!\n");
  bft_printf_flush();

  BFT_MALLOC(r->time_steps, r->ntsteps, cs_real_t);
  BFT_MALLOC(r->iter, r->ntsteps, int);
  BFT_MALLOC(r->order, r->ntsteps, int);
  for (int i = 0; i < r->ntsteps; i++) {
    int it  = tio[i].first;
    int ord = tio[i].second;

    r->iter[i]  = it;
    r->order[i] = ord;

    r->time_steps[i] = fts->getTimeStep(it,ord)->getTime(it,ord);
  }

}

/*============================================================================
 * Public C functions
 *============================================================================*/

BEGIN_C_DECLS

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Creates a new cs_paramedmem_remapper_t instance
 *
 * \param[in] name          name of the remapper
 * \param[in] sel_criteria  cells selection criteria
 * \param[in] fileName      med file name
 * \param[in] meshName      name of the mesh in the med file
 *
 * \return  cs_paramedmem_remapper_t struct
 */
/* -------------------------------------------------------------------------- */

cs_paramedmem_remapper_t *
cs_paramedmem_remapper_create(char       *name,
                              const char *sel_criteria,
                              char       *fileName,
                              char       *meshName)
{

  cs_paramedmem_remapper_t *r = _cs_paramedmem_overlap_create(name);

  _cs_paramedmem_remapper_target_mesh(r, name, sel_criteria);

  _cs_paramedmem_load_paramesh(r, fileName, meshName);

  return r;

}

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Remaps a field from the med file to the local mesh for a given time
 *
 * \param[in] r         pointer to cs_paramedmem_remapper_t struct
 * \param[in] fieldName name of the field to remap from the file
 * \param[in] iter      time iteration to use from the file
 * \param[in] order     time order to use from the file
 *
 * \return  cs_real_t pointer containing the new values on target mesh
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
cs_paramedmem_remap_field_one_time(cs_paramedmem_remapper_t *r,
                                   char                     *fieldName,
                                   int                       dt,
                                   int                       it)
{

  /* Source Field */
  const std::string fname(fieldName);

  MCAuto<MEDFileAnyTypeField1TS> f =
    r->MEDFields->getFieldWithName(fname)->getTimeStep(dt,it);


  MCAuto<MEDFileField1TS>
    sf(MEDCoupling::DynamicCast<MEDFileAnyTypeField1TS,MEDFileField1TS>(f));

  MEDCouplingFieldDouble *tmp_field(sf->field(r->src_mesh));

  MEDCouplingFieldDouble *src_field;
  if (tmp_field->getTypeOfField() == ON_NODES)
    src_field = tmp_field->nodeToCellDiscretization();
  else
    src_field = tmp_field;

  src_field->setNature(IntensiveMaximum);
  r->odec->attachSourceLocalField(src_field);

  /* Target Field */
  MEDCouplingFieldDouble *trg_field = MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  trg_field->setMesh(r->local_mesh->med_mesh);

  DataArrayDouble *arr = DataArrayDouble::New();

  int fdim = src_field->getNumberOfComponents();
  arr->alloc(r->local_mesh->n_elts, fdim);
  trg_field->setArray(arr);
  trg_field->getArray()->decrRef();
  trg_field->setNature(IntensiveMaximum);
  r->odec->attachTargetLocalField(trg_field);

  // Sync the DEC if needed
  if (r->synced != 1) {
    r->odec->synchronize();
    r->synced = 1;
  }

  r->odec->sendData();

  // Write new values
  const double *val_ptr = trg_field->getArray()->getConstPointer();

  cs_real_t *new_vals;
  BFT_MALLOC(new_vals, r->local_mesh->n_elts, cs_real_t);

  const cs_lnum_t *new_connec = r->local_mesh->new_to_old;
  for (int ii = 0; ii < r->local_mesh->n_elts; ii++) {
    int c_id = new_connec[ii];
    new_vals[c_id] = val_ptr[ii];
  }

  return new_vals;

}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Interpolate a given field on the local mesh for a given time
 *
 * \param[in] r             pointer to cs_paramedmem_remapper_t struct
 * \param[in] fieldName     name of the field to remap from the file
 * \param[in] time_choice   Choice of the time interpolation.
 *                          0: Value of field interpolated at t=tval from the
 *                          med file.
 *                          1: Returns field values for the first time step in
 *                          the file. tval is then ignored.
 *                          2: Returns field values for the last time step in
 *                          the file. tval is then ignored.
 * \param[in] tval          requested time instant. If time choice is 0 and
 *                          tval outside of the file time bounds, return value
 *                          will be at the the first time step (if tval < tmin)
 *                          or last time step (if tval > tmax)
 *
 * \return  cs_real_t pointer containing the new values on target mesh
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
cs_paramedmem_remap_field(cs_paramedmem_remapper_t *r,
                          char                     *fieldName,
                          int                       time_choice,
                          double                    tval)
{

  cs_real_t *new_vals = NULL;

  bft_printf(" REMAPPING!\n");
  bft_printf_flush();

  if ( (time_choice == 0 && tval < r->time_steps[0]) ||
        time_choice == 1 ||
        r->ntsteps == 1) {
    /* First instance */
    int it    = r->iter[0];
    int order = r->order[0];
    new_vals = cs_paramedmem_remap_field_one_time(r, fieldName, it, order);

  } else if ( (time_choice == 0 && tval > r->time_steps[r->ntsteps-1]) ||
              time_choice == 2) {
    /* Last instance */
    int it    = r->iter[r->ntsteps-1];
    int order = r->order[r->ntsteps-1];

    new_vals = cs_paramedmem_remap_field_one_time(r, fieldName, it, order);

  } else if (time_choice == 0) {
    /* A given time within the file time bounds*/
    int id1 = -1;
    int id2 = -1;
    for (int i = 0; i < r->ntsteps-1; i++) {
      if (tval > r->time_steps[i] && tval < r->time_steps[i+1]) {
        id1 = i;
        id2 = i+1;
        break;
      }
    }

    cs_real_t t1 = r->time_steps[id1];
    cs_real_t t2 = r->time_steps[id2];

    cs_real_t *vals1 = cs_paramedmem_remap_field_one_time(r,
                                                          fieldName,
                                                          r->iter[id1],
                                                          r->order[id1]);

    cs_real_t *vals2 = cs_paramedmem_remap_field_one_time(r,
                                                          fieldName,
                                                          r->iter[id2],
                                                          r->order[id2]);

    BFT_MALLOC(new_vals, r->local_mesh->n_elts, cs_real_t);
    for (int c_id = 0; c_id < r->local_mesh->n_elts; c_id++) {
      new_vals[c_id] = vals1[c_id] +
                       (vals2[c_id]-vals1[c_id]) *
                       (tval - t1) / (t2-t1);
    }

  }



  return new_vals;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* HAVE_PARAMEDMEM */
