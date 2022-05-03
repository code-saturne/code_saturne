/*============================================================================
 * ParaMEDMEM coupling
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
 * Standard C/C++ library headers
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

#if defined(HAVE_PARAMEDMEM)
#include <string>
#include <vector>
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

#include "cs_medcoupling_mesh.hxx"
#include "cs_paramedmem_coupling.h"

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
 * MEDCoupling writer/reader structure
 *----------------------------------------------------------------------------*/

struct _cs_paramedmem_coupling_t {

#if defined(HAVE_PARAMEDMEM)

  /* Coupling Name */
  std::string                  name;           /* Coupling name */

  /* Current app name */
  ple_coupling_mpi_set_info_t  apps[2];

  cs_medcoupling_mesh_t       *mesh;

  /* index 0 is for send, 1 for recv */

  ParaMESH        *para_mesh;   /* Two meshes, one for send, one for recv */

  InterpKernelDEC *dec;       /* Send data exchange channel */

  std::vector<ParaFIELD *> fields;

#else

  void *para_mesh;
  void *dec;
  void *fields;

#endif

  int dec_synced;
};

/*=============================================================================
 * Private global variables
 *============================================================================*/

#if defined(HAVE_PARAMEDMEM)

static std::vector<cs_paramedmem_coupling_t *> _paramed_couplers;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_PARAMEDMEM)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generate mesh structure from user's defintion
 *
 * \param[in] c   pointer to cs_paramedmem_coupling_t struct
 * \param[in] select_criteria   selection criteria (string)
 * \param[in] elt_dim           mesh dimension (2 or 3)
 *
 */
/*----------------------------------------------------------------------------*/

static void
_generate_coupling_mesh(cs_paramedmem_coupling_t  *c,
                        const char                *select_criteria,
                        int                        elt_dim)
{
  cs_mesh_t *parent_mesh = cs_glob_mesh;

  /* Building the MED representation of the internal mesh */
  int use_bbox = (elt_dim == 3) ? 1 : 0;
  c->mesh = cs_medcoupling_mesh_from_base(parent_mesh,
                                          "CouplingMesh",
                                          select_criteria,
                                          elt_dim,
                                          use_bbox);

  /* Define associated ParaMESH */
  ProcessorGroup *Grp =
    c->dec->isInSourceSide() ? c->dec->getSourceGrp():c->dec->getTargetGrp();

  c->para_mesh = new ParaMESH(c->mesh->med_mesh, *(Grp), "CoupledMesh");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a paramedmem coupling based on an InterpKernelDEC.
 *
 * The latter is created using the the lists of ranks provided as
 * input to this function.
 *
 * \param[in] cpl_name   coupling name
 * \param[in] apps       array of size 2, containing ple_coupling_mpi_set_info_t
 *
 */
/*----------------------------------------------------------------------------*/

static int
_add_paramedmem_coupling(const std::string            cpl_name,
                         ple_coupling_mpi_set_info_t  apps[2])
{
  cs_paramedmem_coupling_t *c = new cs_paramedmem_coupling_t();

  /* Apps identification */
  for (int i = 0; i < 2; i++)
    c->apps[i] = apps[i];

  /* Set coupling name */
  c->name = cpl_name;

  c->mesh = NULL;

  std::set<int> grp1_ids;
  std::set<int> grp2_ids;
  for (int i = 0; i < apps[0].n_ranks; i++)
    grp1_ids.insert(apps[0].root_rank + i);
  for (int i = 0; i < apps[1].n_ranks; i++)
    grp2_ids.insert(apps[1].root_rank + i);

  c->dec =new InterpKernelDEC(grp1_ids, grp2_ids);

  _paramed_couplers.push_back(c);

  int retval = _paramed_couplers.size() - 1;

  return retval;
}

#endif /* HAVE_PARAMEDMEM */

/*============================================================================
 * Public C functions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve coupling struct pointer by id
 *
 * \param[in] cpl_id  index of the sought coupling
 *
 * \return pointer to cs_paramedmem_coupling_t struct. Raise an error if
 * the coupling does not exist.
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_by_id(int  cpl_id)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(cpl_id);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);
  return NULL;

#else

  if (cpl_id < 0 || cpl_id > (int)_paramed_couplers.size())
    bft_error(__FILE__, __LINE__, 0,
              _("Error: coupling with id %d does not exist\n"), cpl_id);

  cs_paramedmem_coupling_t *c = _paramed_couplers[cpl_id];

  return c;

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve coupling struct pointer by name
 *
 * \param[in] name  name of the coupling
 *
 * \return pointer to cs_paramedmem_coupling_t struct or NULL if not found.
 *
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_by_name(const char *name)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(name);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);
  return NULL;

#else

  cs_paramedmem_coupling_t *c = NULL;

  for (size_t i = 0; i < _paramed_couplers.size(); i++) {
    if (strcmp(_paramed_couplers[i]->name.c_str(), name) == 0) {
      c = _paramed_couplers[i];
      break;
    }
  }

  return c;

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new ParaMEDMEM coupling
 *
 * \param[in] app1_name  Name of app n°1 or NULL if calling app is app1
 * \param[in] app2_name  Name of app n°2 or NULL if calling app is app2
 * \param[in] cpl_name   Name of the coupling. If NULL an automatic name is generated.
 *
 * \return pointer to newly created cs_paramedmem_coupling_t structure.
 *
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create(const char  *app1_name,
                              const char  *app2_name,
                              const char  *cpl_name)
{
  cs_paramedmem_coupling_t *c = NULL;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(app1_name);
  CS_NO_WARN_IF_UNUSED(app2_name);
  CS_NO_WARN_IF_UNUSED(cpl_name);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

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
      if (strcmp(ai.app_name, app1_name) == 0)
        apps[0] = ai;
    }
    if (app2_name != NULL) {
      if (strcmp(ai.app_name, app2_name) == 0)
        apps[1] = ai;
    }
  }

  for (int i = 0; i < 2; i++) {
    if (apps[i].root_rank < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error: Partner app was not found...\n"));
  }

  /* Set coupling name */
  std::string _name = "";
  if (cpl_name == NULL) {
    /* string is <a_name>_<p_name>_cpl */
    std::stringstream ss;
    ss << apps[0].app_name << "_" << apps[1].app_name << "_cpl";
    _name = ss.str();
  } else {
    _name = cpl_name;
  }

  int cpl_id = _add_paramedmem_coupling(_name, apps);

  c = cs_paramedmem_coupling_by_id(cpl_id);

#endif

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a given ParaMEDMEM coupling structure
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_destroy(cs_paramedmem_coupling_t  *c)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  if (c != NULL) {

    /* Destroy fields */
    c->fields.clear();

    /* Deallocate mesh */
    delete c->para_mesh;

    // Mesh will deallocated afterwards since it can be shared
    c->mesh = NULL;

    /* Destroy DECs */
    delete c->dec;

    delete c;
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all coupling structures
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_finalize(void)
{
#if defined(HAVE_PARAMEDMEM)

  /* Clear vector */
  _paramed_couplers.clear();

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a selection criteria
 *
 * \param[in] c         pointer to cs_paramedmem_coupling_t struct
 * \param[in] sel_crit  geometrical selection criteria (string)
 * \param[in] elt_dim   dimension of coupled elements
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_add_mesh_from_criteria(cs_paramedmem_coupling_t  *c,
                                     const char                *sel_crit,
                                     int                        elt_dim)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(sel_crit);
  CS_NO_WARN_IF_UNUSED(elt_dim);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  _generate_coupling_mesh(c, sel_crit, elt_dim);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a cs_zone_t pointer
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t struct
 * \param[in] zone  pointer to cs_zone_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_add_mesh_from_zone(cs_paramedmem_coupling_t  *c,
                                 const cs_zone_t           *zone)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(zone);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  /* Get location id */
  int location_id = zone->location_id;

  /* Get selection criteria from location id */
  const char *sel_crit = cs_mesh_location_get_selection_string(location_id);

  /* Get location id and deduce element dimension */
  int elt_dim = -1;
  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(location_id);
  if (loc_type == CS_MESH_LOCATION_CELLS)
    elt_dim = 3;
  else if (loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
    elt_dim = 2;

  if (elt_dim < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: wrong dimension for considered zone\n"));

  _generate_coupling_mesh(c, sel_crit, elt_dim);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of defined couplings
 *
 * \return number of defined couplings (int)
 */
/*----------------------------------------------------------------------------*/
int
cs_paramedmem_get_number_of_couplings(void)
{
  int retval = 0;

#if !defined(HAVE_PARAMEDMEM)

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  retval = _paramed_couplers.size();

#endif

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements of coupled mesh
 *
 * \param[in] coupling  pointer to cs_paramedmem_coupling_t struct
 *
 * \return number of elements in mesh associated to coupling
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_paramedmem_mesh_get_n_elts(const cs_paramedmem_coupling_t  *coupling)
{
  cs_lnum_t retval = 0;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(coupling);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  retval = coupling->mesh->n_elts;

#endif

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get indirection list for elements in coupled mesh
 *
 * \param[in] coupling  pointer to cs_paramedmem_coupling_t struct
 *
 * \return cs_lnum_t pointer to indirection list
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_paramedmem_mesh_get_elt_list(const cs_paramedmem_coupling_t  *coupling)
{
  const cs_lnum_t *retval = NULL;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(coupling);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  retval = coupling->mesh->elt_list;

#endif

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a coupled field
 *
 * \param[in] c             pointer to cs_paramedmem_coupling_t struct
 * \param[in] name          name of field
 * \param[in] dim           field dimension
 * \param[in] field_nature  field nature flag
 * \param[in] space_discr   field space discretisation (nodes or cells)
 * \param[in] time_discr    field coupling time discretisation
 *
 * \return index of field within the storing vector
 *
 */
/*----------------------------------------------------------------------------*/

int
cs_paramedmem_def_coupled_field(cs_paramedmem_coupling_t  *c,
                                const char                *name,
                                int                        dim,
                                cs_medcpl_field_nature_t   field_nature,
                                cs_medcpl_space_discr_t    space_discr,
                                cs_medcpl_time_discr_t     time_discr)
{
  int f_id = -1;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(dim);
  CS_NO_WARN_IF_UNUSED(field_nature);
  CS_NO_WARN_IF_UNUSED(space_discr);
  CS_NO_WARN_IF_UNUSED(time_discr);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  /* Prepare coupling structure */
  TypeOfField type = ON_CELLS;
  switch (space_discr) {
  case CS_MEDCPL_ON_CELLS:
    type = ON_CELLS;
    break;

  case CS_MEDCPL_ON_NODES:
    type = ON_NODES;
    break;
  }

  TypeOfTimeDiscretization td = NO_TIME;
  switch (time_discr) {
  case CS_MEDCPL_NO_TIME:
    td = NO_TIME;
    break;

  case CS_MEDCPL_ONE_TIME:
    td = ONE_TIME;
    break;

  case CS_MEDCPL_LINEAR_TIME:
    td = LINEAR_TIME;
    break;
  }


  /* Build ParaFIELD object if required */
  ComponentTopology comp_topo(dim);
  ParaFIELD *pf = new ParaFIELD(type, td, c->para_mesh, comp_topo);

  if (type == ON_CELLS)
    c->dec->setMethod("P0");
  else
    c->dec->setMethod("P1");

  c->fields.push_back(pf);
  /* TODO: setNature should be set by caller to allow for more options */

  /* Set nature of field
   * As default we use intensive conservation since code_saturne is
   * a CFD code :)
   */
  NatureOfField nature = IntensiveConservation;
  switch (field_nature) {
  case CS_MEDCPL_FIELD_EXT_CONSERVATION:
    nature = ExtensiveConservation;
    break;

  case CS_MEDCPL_FIELD_EXT_MAXIMUM:
    nature = ExtensiveMaximum;
    break;

  case CS_MEDCPL_FIELD_INT_CONSERVATION:
    nature = IntensiveConservation;
    break;

  case CS_MEDCPL_FIELD_INT_MAXIMUM:
    nature = IntensiveMaximum;
    break;

  default:
    assert(0);
  }

  pf->getField()->setNature(nature);

  pf->getField()->setName(name);

#endif

  return f_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a coupled field based on a cs_field_t pointer
 *
 * \param[in] c           pointer to cs_paramedmem_coupling_t struct
 * \param[in] f           pointer to cs_field_t struct
 * \param[in] fn          field nature flag
 * \param[in] time_discr  field coupling time discretisation
 *
 * \return index of field within the storing vector
 */
/*----------------------------------------------------------------------------*/

int
cs_paramedmem_def_coupled_field_from_cs_field(cs_paramedmem_coupling_t *c,
                                              cs_field_t               *f,
                                              cs_medcpl_field_nature_t  fn,
                                              cs_medcpl_time_discr_t    td)
{
  int f_id = -1;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(f);
  CS_NO_WARN_IF_UNUSED(fn);
  CS_NO_WARN_IF_UNUSED(td);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(f->location_id);
  cs_medcpl_space_discr_t sd       = CS_MEDCPL_ON_CELLS;

  if (loc_type == CS_MESH_LOCATION_CELLS ||
      loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
    sd = CS_MEDCPL_ON_CELLS;
  else if (loc_type == CS_MESH_LOCATION_VERTICES)
    sd = CS_MEDCPL_ON_NODES;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Non-compatible field location for '%s'\n"), f->name);

  f_id = cs_paramedmem_def_coupled_field(c,
                                         f->name,
                                         f->dim,
                                         fn,
                                         sd,
                                         td);
#endif

  return f_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write values before sending operation
 *
 * \param[in] c      pointer to cs_paramedmem_coupling_t structure
 * \param[in] name   name of field
 * \param[in] values array of values to write
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_field_export(cs_paramedmem_coupling_t  *c,
                           const char                *name,
                           const double               values[])
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(values);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  MEDCouplingFieldDouble *f = NULL;
  for (size_t i = 0; i < c->fields.size(); i++) {
    if (strcmp(name, c->fields[i]->getField()->getName().c_str()) == 0) {
      f = c->fields[i]->getField();
      break;
    }
  }

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Could not find field '%s'\n"), name);

  double  *val_ptr = f->getArray()->getPointer();
  const int dim = f->getNumberOfComponents();

  /* Assign element values */
  /*-----------------------*/
  cs_lnum_t n_elts = c->mesh->n_elts;
  cs_lnum_t *elt_list = c->mesh->elt_list;
  if (elt_list == NULL) {
    for (cs_lnum_t i = 0; i < dim*n_elts; i++)
      val_ptr[i] = values[i];
  } else {
    for (int j = 0; j < dim; j++) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = elt_list[i];
        val_ptr[i + j*n_elts] = values[c_id + j*n_elts];
      }
    }
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read values before sending operation
 *
 * \param[in] c      pointer to cs_paramedmem_coupling_t structure
 * \param[in] name   name of field
 * \param[in] values array in which values will be stored
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_field_import(cs_paramedmem_coupling_t  *c,
                           const char                *name,
                           double                     values[])
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(values);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  MEDCouplingFieldDouble *f = NULL;
  for (size_t i = 0; i < c->fields.size(); i++) {
    if (strcmp(name, c->fields[i]->getField()->getName().c_str()) == 0) {
      f = c->fields[i]->getField();
      break;
    }
  }

  const double  *val_ptr = f->getArray()->getConstPointer();
  const int dim = f->getNumberOfComponents();

  /* Import element values */
  /*-----------------------*/

  cs_lnum_t *connec = c->mesh->new_to_old;
  cs_lnum_t  n_elts = c->mesh->n_elts;
  if (connec != NULL) {
    for (int j = 0; j < dim; j++) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = connec[i];
        values[c_id + j*n_elts] = val_ptr[i + j*n_elts];
      }
    }
  } else {
    for (int j = 0; j < dim; j++) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        values[i + j*n_elts] = val_ptr[i + j*n_elts];
      }
    }
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sync the coupling's InterpKernelDEC
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_sync_dec(cs_paramedmem_coupling_t  *c)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  if (c->dec_synced == 0) {
    c->dec->synchronize();
    c->dec_synced = 1;
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of field attached to DEC
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_send_data(cs_paramedmem_coupling_t  *c)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  c->dec->sendData();

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of field attached to DEC
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_recv_data(cs_paramedmem_coupling_t  *c)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  c->dec->recvData();

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a field to InterpKernelDEC for send operation using its index
 *
 * \param[in] c         pointer to cs_paramedmem_coupling_t structure
 * \param[in] field_id  index of field in storing vector
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_attach_field_by_id(cs_paramedmem_coupling_t  *c,
                                 int                        field_id)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(field_id);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  c->dec->attachLocalField(c->fields[field_id]);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a field to InterpKernelDEC for send operation using its name
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field (string)
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_attach_field_by_name(cs_paramedmem_coupling_t  *c,
                                   const char                *name)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(name);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  ParaFIELD *pf = NULL;

  for (size_t i = 0; i < c->fields.size(); i++) {
    if (strcmp(name, c->fields[i]->getField()->getName().c_str()) == 0) {
      pf = c->fields[i];
      break;
    }
  }

  if (pf == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Could not find field '%s'\n"), name);

  c->dec->attachLocalField(pf);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of a field. If vals pointer is non-null,
 * values are updated before send
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field
 * \param[in] vals  array of values to write
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_send_field_vals(cs_paramedmem_coupling_t *c,
                              const char               *name,
                              const double             *vals)
{
#if defined(HAVE_PARAMEDMEM)

  /* If provided, export data to DEC */
  if (vals != NULL)
    cs_paramedmem_field_export(c, name, vals);

  /* Attach field to DEC for sending */
  cs_paramedmem_attach_field_by_name(c, name);

  /* Send data */
  cs_paramedmem_send_data(c);

#else

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(vals);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of a field.
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field
 * \param[in] vals  array of values to write
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_recv_field_vals(cs_paramedmem_coupling_t *c,
                              const char               *name,
                              double                   *vals)
{
#if defined(HAVE_PARAMEDMEM)

  /* Attach field to DEC for sending */
  cs_paramedmem_attach_field_by_name(c, name);

  /* Recieve data */
  cs_paramedmem_recv_data(c);

  /* Read values */
  cs_paramedmem_field_import(c, name, vals);

#else

  CS_NO_WARN_IF_UNUSED(c);
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(vals);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log ParaMEDMEM coupling setup information
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_log_setup(void)
{
#if defined(HAVE_PARAMEDMEM)

  /* Get number of couplings */
  int ncpls = cs_paramedmem_get_number_of_couplings();

  if (ncpls > 0) {

    cs_log_printf(CS_LOG_SETUP,
                  _("\nParaMEDMEM coupling\n"
                    "-------------------\n\n"
                    "  number of couplings: %d\n"),
                  ncpls);

    for (int cpl_id = 0; cpl_id < ncpls; cpl_id++) {
      cs_paramedmem_coupling_t *c = _paramed_couplers[cpl_id];

      cs_log_printf(CS_LOG_SETUP,
                    _("\n"
                      "  Coupling: %s\n"
                      "    App1: %s\n"
                      "    App2: %s\n"),
                    c->name.c_str(), c->apps[0].app_name, c->apps[1].app_name);

      if (c->mesh->elt_dim == 3)
        cs_log_printf(CS_LOG_SETUP, _("    Type: Volume coupling\n"));
      else if (c->mesh->elt_dim == 2)
        cs_log_printf(CS_LOG_SETUP, _("    Type: Boundary coupling\n"));

    }

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "  Couplings not usable in non-MPI builds.n"));
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize couplings based on user functions
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_init(void)
{
  cs_user_paramedmem_define_couplings();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize coupled mesh and fields based on user functions
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_define_mesh_fields(void)
{
#if defined(HAVE_PARAMEDMEM)

  cs_user_paramedmem_define_meshes();
  cs_user_paramedmem_define_fields();

  /* Sync dec after declaration using first field */
  int ncpls = cs_paramedmem_get_number_of_couplings();
  for (int cpl_id = 0; cpl_id < ncpls; cpl_id++) {
    cs_paramedmem_coupling_t *c = cs_paramedmem_coupling_by_id(cpl_id);
    if (c->fields.size() < 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error: You did not define fields to couple with "
                  "coupling '%s'\n"),
                c->name.c_str());

    cs_paramedmem_attach_field_by_id(c, 0);
    cs_paramedmem_sync_dec(c);
  }

  cs_paramedmem_coupling_log_setup();

#endif /* defined(HAVE_PARAMEDMEM) */
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
