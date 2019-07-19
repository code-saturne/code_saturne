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
 * Header for the current file
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

  char                   *name;           /* Coupling name */

  cs_medcoupling_mesh_t  *local_mesh;     /* Local cs_mesh in MED format */

  MEDFileUMesh           *src_mesh;
  MEDFileFields          *MEDFields;

  int                     ntsteps;
  int                    *iter;
  int                    *order;
  cs_real_t              *time_steps;

  /* Bounding sphere for localization */
  cs_real_t               _sphere_cen[3] = {0., 0., 0.};
  cs_real_t               _sphere_rad    = 1.e20;

  OverlapDEC             *odec;           /* Overlap data exchange channel */
  int                     synced;
};

struct _mesh_transformation_t {

  int   type          = -1; /* 0: rotation, 1: translation */

  cs_real_t vector[3] = {0., 0., 0.};
  cs_real_t center[3] = {0., 0., 0.};
  cs_real_t angle     = 0.;
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

static int                         _n_remappers = 0;
static cs_paramedmem_remapper_t  **_remapper = NULL;

/*----------------------------------------------------------------------------*/

static int                      _n_transformations = 0;
static _mesh_transformation_t **_transformations = NULL;

static bool _transformations_applied = false;

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a mesh transformation (rotation or translation)
 *
 * \return  a mesh_transformation pointer
 */
/*----------------------------------------------------------------------------*/

static _mesh_transformation_t *
_cs_paramedmem_create_transformation(int             type,
                                     const cs_real_t center[3],
                                     const cs_real_t vector[3],
                                     const cs_real_t angle)
{

  _mesh_transformation_t *mt = NULL;
  BFT_MALLOC(mt, 1, _mesh_transformation_t);

  mt->type  = type;
  mt->angle = angle;

  for (int i = 0; i < 3; i++) {
    mt->center[i] = center[i];
    mt->vector[i] = vector[i];
  }

  return mt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Reset all mesh transformations
 *
 */
/*----------------------------------------------------------------------------*/

static void
_cs_paramedmem_reset_transformations(void)
{
  if (_transformations != NULL) {
    for (int i = 0; i < _n_transformations; i++) {
      BFT_FREE(_transformations[i]);
    }
    BFT_FREE(_transformations);
  }

  _n_transformations = 0;

  _transformations_applied = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the different mesh transformations and update the bounding
 *          sphere
 */
/*----------------------------------------------------------------------------*/

static void
_cs_paramedmem_apply_transformations(MEDCouplingFieldDouble   *field,
                                     cs_paramedmem_remapper_t *r)
{
  if (_transformations_applied == false) {
    for (int i = 0; i < _n_transformations; i++) {
      _mesh_transformation_t *mt = _transformations[i];
      if (mt->type == 0) {
        /* Rotation */
        field->getMesh()->rotate(mt->center, mt->vector, mt->angle);

        /* Bounding sphere */
        cs_real_t c = cos(mt->angle);
        cs_real_t s = sin(mt->angle);

        cs_real_t norm = cs_math_3_norm(mt->vector);
        cs_real_t vec[3] = {0., 0., 0.};
        for (int id = 0; id < 3; id++)
          vec[id] = mt->vector[id]/norm;

        cs_real_t xyz[3] = {0., 0., 0.};
        for (int id = 0; id < 3; id++)
          xyz[id] = r->_sphere_cen[id] - mt->center[id];

        r->_sphere_cen[0] = (vec[0]*vec[0]*(1. - c) + c)        * xyz[0]
                          + (vec[0]*vec[1]*(1. - c) - vec[2]*s) * xyz[1]
                          + (vec[0]*vec[2]*(1. - c) + vec[1]*s) * xyz[2]
                          + mt->center[0];

        r->_sphere_cen[1] = (vec[0]*vec[1]*(1. - c) + vec[2]*s) * xyz[0]
                          + (vec[1]*vec[1]*(1. - c) + c)        * xyz[1]
                          + (vec[1]*vec[2]*(1. - c) - vec[0]*s) * xyz[2]
                          + mt->center[1];

        r->_sphere_cen[2] = (vec[2]*vec[0]*(1. - c) - vec[1]*s) * xyz[0]
                          + (vec[2]*vec[1]*(1. - c) + vec[0]*s) * xyz[1]
                          + (vec[2]*vec[2]*(1. - c) + c)        * xyz[2]
                          + mt->center[2];

      } else if (mt->type == 1) {
        /* Translation */
        field->getMesh()->translate(mt->vector);

        /* Bounding sphere */
        for (int id = 0; id < 3; id++)
          r->_sphere_cen[id] += mt->vector[id];
      }
    }
  }

  _transformations_applied = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Returns an array containing the ranks of Code_Saturne processes in
 *          MPI_COMM_WORLD.
 *
 * \return  array of ranks in MPI_COMM_WORLD
 */
/*----------------------------------------------------------------------------*/

static int *
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create the parallel remapper structure based on ParaMEDMEM OverlapDEC
 *
 * \param[in] name    name of the remapper
 *
 * \return  pointer to cs_paramedmem_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

static cs_paramedmem_remapper_t *
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the target mesh for interpolation
 *
 * \param[in] r               pointer to cs_paramedmem_remapper_t struct
 * \param[in] name            mesh name
 * \param[in] select_criteria selection criteria for needed cells
 */
/*----------------------------------------------------------------------------*/

static void
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
/*!
 * \brief   Load the mesh parts on each process
 *
 * \param[in] r         pointer to cs_paramedmem_remapper_t struct
 * \param[in] fileName  name of med file containing the data
 * \param[in] meshName  name of the mesh to read in the med file
 *
 */
/*----------------------------------------------------------------------------*/

static void
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

  std::vector< std::pair<int,int> > tio = fts->getIterations();

  r->ntsteps    = tio.size();

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

/*----------------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------------*/

cs_paramedmem_remapper_t *
cs_paramedmem_remapper_create(char       *name,
                              const char *sel_criteria,
                              char       *fileName,
                              char       *meshName,
                              cs_real_t   center[3],
                              cs_real_t   radius)
{
  if (_remapper == NULL)
    BFT_MALLOC(_remapper, 1, cs_paramedmem_remapper_t *);
  else
    BFT_REALLOC(_remapper, _n_remappers+1, cs_paramedmem_remapper_t *);

  cs_paramedmem_remapper_t *r = _cs_paramedmem_overlap_create(name);

  _cs_paramedmem_remapper_target_mesh(r, name, sel_criteria);

  _cs_paramedmem_load_paramesh(r, fileName, meshName);

  r->_sphere_rad = radius;
  for (int i = 0; i < 3; i++)
    r->_sphere_cen[i] = center[i];

  _remapper[_n_remappers] = r;
  _n_remappers++;

  return r;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief get a remapper by its name
 *
 * \param[in] name  name of the remapper
 *
 * \return  pointer to cs_paramedmem_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_remapper_t *
cs_paramedmem_remapper_by_name_try(const char *name)
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

/*----------------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_paramedmem_remap_field_one_time(cs_paramedmem_remapper_t *r,
                                   char                     *fieldName,
                                   cs_real_t                 default_val,
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

  _cs_paramedmem_apply_transformations(src_field, r);

  r->odec->attachSourceLocalField(src_field);

  /* Target Field */
  MEDCouplingFieldDouble *trg_field
    = MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  trg_field->setMesh(r->local_mesh->med_mesh);

  DataArrayDouble *arr = DataArrayDouble::New();

  int fdim = src_field->getNumberOfComponents();
  arr->alloc(r->local_mesh->n_elts, fdim);
  trg_field->setArray(arr);
  trg_field->getArray()->decrRef();
  trg_field->setNature(IntensiveMaximum);
  r->odec->attachTargetLocalField(trg_field);

  r->odec->setDefaultValue(default_val);
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

  const int npts = trg_field->getNumberOfValues();

  const cs_real_3_t * xyzcen = (cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;

  for (int ii = 0; ii < npts; ii++) {
    int c_id = new_connec[ii];
    if (cs_math_3_distance(xyzcen[c_id], r->_sphere_cen) < r->_sphere_rad)
      new_vals[c_id] = val_ptr[ii];
    else
      new_vals[c_id] = default_val;
  }

  return new_vals;
}

/*----------------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_paramedmem_remap_field(cs_paramedmem_remapper_t *r,
                          char                     *fieldName,
                          cs_real_t                 dval,
                          int                       time_choice,
                          double                    tval)
{
  cs_real_t *new_vals = NULL;

  if ((time_choice == 0 && tval < r->time_steps[0]) ||
      time_choice == 1 ||
      r->ntsteps == 1) {
    /* First instance */
    int it    = r->iter[0];
    int order = r->order[0];
    new_vals = cs_paramedmem_remap_field_one_time(r, fieldName, dval, it, order);

  } else if ( (time_choice == 0 && tval > r->time_steps[r->ntsteps-1]) ||
              time_choice == 2) {
    /* Last instance */
    int it    = r->iter[r->ntsteps-1];
    int order = r->order[r->ntsteps-1];

    new_vals = cs_paramedmem_remap_field_one_time(r, fieldName, dval, it, order);

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
                                                          dval,
                                                          r->iter[id1],
                                                          r->order[id1]);

    cs_real_t *vals2 = cs_paramedmem_remap_field_one_time(r,
                                                          fieldName,
                                                          dval,
                                                          r->iter[id2],
                                                          r->order[id2]);

    BFT_MALLOC(new_vals, r->local_mesh->n_elts, cs_real_t);
    for (int c_id = 0; c_id < r->local_mesh->n_elts; c_id++) {
      new_vals[c_id] = vals1[c_id] +
                       (vals2[c_id]-vals1[c_id]) *
                       (tval - t1) / (t2-t1);
    }

  }

  r->synced = 0;
  _cs_paramedmem_reset_transformations();

  return new_vals;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] r            pointer to the cs_paramedmem_remapper_t struct
 * \param[in] translation  translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_remapper_translate(cs_paramedmem_remapper_t  *r,
                                 cs_real_t                  translation[3])
{
  if (_transformations == NULL)
    BFT_MALLOC(_transformations, 1, _mesh_transformation_t *);
  else
    BFT_REALLOC(_transformations, _n_transformations+1, _mesh_transformation_t *);

  cs_real_t cen[3] = {0.,0.,0.};
  _transformations[_n_transformations] =
    _cs_paramedmem_create_transformation(1, cen, translation, 0.);

  _n_transformations++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rotate the mesh using a center point, axis and angle
 *
 * \param[in] r          pointer to the cs_paramedmem_remapper_t struct
 * \param[in] invariant  coordinates of the invariant point
 * \param[in] axis       rotation axis vector
 * \param[in] angle      rotation angle in radians
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_remapper_rotate(cs_paramedmem_remapper_t  *r,
                              cs_real_t                  invariant[3],
                              cs_real_t                  axis[3],
                              cs_real_t                  angle)
{
  if (_transformations == NULL)
    BFT_MALLOC(_transformations, 1, _mesh_transformation_t *);
  else
    BFT_REALLOC(_transformations, _n_transformations+1, _mesh_transformation_t *);

  _transformations[_n_transformations] =
    _cs_paramedmem_create_transformation(0, invariant, axis, angle);

  _n_transformations++;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* HAVE_PARAMEDMEM */
