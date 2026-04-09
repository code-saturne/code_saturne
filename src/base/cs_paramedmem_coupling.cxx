/*============================================================================
 * MEDCoupling ParaMESH/ParaFIELD wrapper functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_coupling.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"
#include "base/cs_prototypes.h"
#include "base/cs_selector.h"
#include "base/cs_timer.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_connect.h"

#include "fvm/fvm_defs.h"
#include "fvm/fvm_nodal_from_desc.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_paramedmem_coupling_priv.h"

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*=============================================================================
 * Private global variables
 *============================================================================*/

#if defined(HAVE_PARAMEDMEM)

static std::vector<cs_paramedmem_coupling_t *> _paramed_couplers;

#endif

#if defined(HAVE_PARAMEDMEM)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a paramedmem coupling based on an InterpKernelDEC.
 *
 * The latter is created using the the lists of ranks provided as
 * input to this function.
 *
 * \param[in]  cpl_name  coupling name
 * \param[in]  apps      array containing ple_coupling_mpi_set_info_t
 */
/*----------------------------------------------------------------------------*/

static int
_add_paramedmem_coupling(const std::string            cpl_name,
                         ple_coupling_mpi_set_info_t  apps[2])
{
#if defined(HAVE_MPI)
  cs_paramedmem_coupling_t *c = new cs_paramedmem_coupling_t();

  c->dec_synced = false;
  c->_curr_field = nullptr;
  c->para_mesh = nullptr;

  /* Apps identification */
  for (int i = 0; i < 2; i++)
    c->apps[i] = apps[i];

  // Create intracomm if more than 2 apps, to avoid deadlock inside MEDCoupling
  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();
  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  const MPI_Comm _set_comm = ple_coupling_mpi_set_get_base_comm(mpi_apps);
  c->comm = _set_comm;

  std::set<int> grp1_ids;
  std::set<int> grp2_ids;

  if (n_apps > 2) {
    const int app_id = ple_coupling_mpi_set_get_app_id(mpi_apps);
    ple_coupling_mpi_set_info_t _app =
      ple_coupling_mpi_set_get_info(mpi_apps, app_id);

    int dist_root = -1;
    int loc_id = -1;
    if (_app.root_rank == apps[0].root_rank) {
      dist_root = apps[1].root_rank;
      loc_id = 0;
    }
    else if (_app.root_rank == apps[1].root_rank) {
      dist_root = apps[0].root_rank;
      loc_id = 1;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Could not find partner app...\n", __func__);

    int local_range[2] = {-1, -1};
    int distant_range[2] = {-1, -1};
    ple_coupling_mpi_intracomm_create(_set_comm,
                                      cs_glob_mpi_comm,
                                      dist_root,
                                      &(c->comm),
                                      local_range,
                                      distant_range);
    if (loc_id == 0) {
      for (int i = local_range[0]; i < local_range[1]; i++)
        grp1_ids.insert(i);
      for (int i = distant_range[0]; i < distant_range[1]; i++)
        grp2_ids.insert(i);
    }
    else {
      for (int i = local_range[0]; i < local_range[1]; i++)
        grp2_ids.insert(i);
      for (int i = distant_range[0]; i < distant_range[1]; i++)
        grp1_ids.insert(i);
    }
  }
  else {
    for (int i = 0; i < apps[0].n_ranks; i++)
      grp1_ids.insert(apps[0].root_rank + i);
    for (int i = 0; i < apps[1].n_ranks; i++)
      grp2_ids.insert(apps[1].root_rank + i);
  }
  /* Set coupling name */
  c->_name = cpl_name;

  c->mesh = nullptr;

  c->dec = new InterpKernelDEC(grp1_ids, grp2_ids, c->comm);

  _paramed_couplers.push_back(c);

  int retval = _paramed_couplers.size() - 1;

  return retval;
#else
  return -1;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a paramedmem coupling for "dry run" mode.
 *
 * \param[in]  cpl_name  coupling name
 */
/*----------------------------------------------------------------------------*/

static int
_add_paramedmem_coupling_dry_run(const std::string cpl_name)
{
  cs_paramedmem_coupling_t *c = new cs_paramedmem_coupling_t();

  c->dec_synced  = false;
  c->_curr_field = nullptr;
  c->para_mesh  = nullptr;

  memset(c->apps, 0, 2 * sizeof(ple_coupling_mpi_set_info_t));

  /* Set coupling name */
  c->_name = cpl_name;

  c->mesh = nullptr;

  c->dec = nullptr;

  _paramed_couplers.push_back(c);

  int retval = _paramed_couplers.size() - 1;

  return retval;
}

#endif /* HAVE_PARAMEDMEM */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public functions
 *============================================================================*/

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
  return nullptr;

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
 * \return pointer to cs_paramedmem_coupling_t struct or null if not found.
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
  return nullptr;

#else

  cs_paramedmem_coupling_t *c = nullptr;

  for (size_t i = 0; i < _paramed_couplers.size(); i++) {
    if (strcmp(_paramed_couplers[i]->getName().c_str(), name) == 0) {
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
 * \param[in] app1_name  Name of app n°1 or null if calling app is app1
 * \param[in] app2_name  Name of app n°2 or null if calling app is app2
 * \param[in] cpl_name   Name of the coupling.
 *                       If null an automatic name is generated.
 *
 * \return pointer to newly created cs_paramedmem_coupling_t structure.
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create(const char  *app1_name,
                              const char  *app2_name,
                              const char  *cpl_name)
{
  cs_paramedmem_coupling_t *c = nullptr;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(app1_name);
  CS_NO_WARN_IF_UNUSED(app2_name);
  CS_NO_WARN_IF_UNUSED(cpl_name);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  /* Check that at least on app name is provided */
  if (app1_name == nullptr && app2_name == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: coupled application name not provided.\n"),
              __func__);

  /* Find which app is which */
  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);
  const int app_id = ple_coupling_mpi_set_get_app_id(mpi_apps);

  /* Get each application's ranks */
  ple_coupling_mpi_set_info_t apps[2];
  for (int i = 0; i < 2; i++)
    apps[i] = ple_coupling_mpi_set_get_info(mpi_apps, -1);

  int l_id = -1;
  if (app1_name == nullptr)
    l_id = 0;
  else if (app2_name == nullptr)
    l_id = 1;

  for (int i = 0; i < n_apps; i++) {
    ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (l_id > -1) {
      if (app_id == i)
        apps[l_id] = ai;
    }
    if (app1_name != nullptr) {
      if (strcmp(ai.app_name, app1_name) == 0)
        apps[0] = ai;
    }
    if (app2_name != nullptr) {
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
  if (cpl_name == nullptr) {
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
 * \brief Create a new ParaMEDMEM handler structure with no actual coupling.
 *
 * This can be useful for a "dry run" when setting up a coupling, so as to
 * first debug local commands before actually running in coupled mode.
 *
 * In this case, data "received" matches the initialized values.
 *
 * \param[in] app1_name  Name of app n°1 or null if calling app is app1
 * \param[in] app2_name  Name of app n°2 or null if calling app is app2
 * \param[in] cpl_name   Name of the coupling.
 *                       If null an automatic name is generated.
 *
 * \return pointer to newly created cs_paramedmem_coupling_t structure.
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create_uncoupled(const char  *cpl_name)
{
  cs_paramedmem_coupling_t *c = nullptr;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(cpl_name);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: %s cannot be called without "
              "MEDCoupling MPI support."), __func__);

#else

  /* Set coupling name */
  std::string _name = cpl_name;

  int cpl_id = _add_paramedmem_coupling_dry_run(_name);

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

  if (c != nullptr) {

    /* Destroy fields */
    c->fields.clear();

    c->_curr_field = nullptr;

    /* Deallocate mesh */
    delete c->para_mesh;

    // Mesh will deallocated afterwards since it can be shared
    c->mesh = nullptr;

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

      c->log();
    }

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "  Couplings not usable in non-MPI builds.n"));
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize couplings based on user functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_init(void)
{
  cs_user_paramedmem_define_couplings();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize coupled mesh and fields based on user functions.
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
      bft_error(__FILE__,
                __LINE__,
                0,
                _("Error: You did not define fields to couple with "
                  "coupling '%s'\n"),
                c->getName().c_str());

    c->attach_field_by_id(0);
    c->sync_dec();
  }

  cs_paramedmem_coupling_log_setup();

#endif /* defined(HAVE_PARAMEDMEM) */
}

#ifdef __cplusplus

#if defined(HAVE_PARAMEDMEM)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a local field
 *
 */
/*----------------------------------------------------------------------------*/

void
#if USE_PARAFIELD == 1
_cs_paramedmem_coupling_t::_attachLocalField(const ParaFIELD *field)
#else
_cs_paramedmem_coupling_t::_attachLocalField(
  const MEDCouplingFieldDouble *field)
#endif
{
  if (field == nullptr) {
    return;
  }

  if (this->_curr_field) {
    TypeOfField   type, type_curr;
    NatureOfField nature, nature_curr;

#if USE_PARAFIELD == 1
    type      = field->getField()->getTypeOfField();
    type_curr = this->_curr_field->getField()->getTypeOfField();

    nature      = field->getField()->getNature();
    nature_curr = this->_curr_field->getField()->getNature();
#else
    type      = field->getTypeOfField();
    type_curr = this->_curr_field->getTypeOfField();

    nature      = field->getNature();
    nature_curr = this->_curr_field->getNature();
#endif

    if (type != type_curr) {
      bft_error(__FILE__,
                __LINE__,
                0,
                _("%s: Impossible to attach a field with a different nature."),
                __func__);
    }

    if (nature != nature_curr) {
      this->dec_synced = false;
    }
  }
  this->_curr_field = field;

  if (this->dec != nullptr) {
    this->dec->attachLocalField(field);
  }
};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a field by its name
 *
 * \param[in] name          name of field
 *
 * \return medcoupling field
 */
/*----------------------------------------------------------------------------*/

MEDCouplingFieldDouble *
_cs_paramedmem_coupling_t::_get_field(const char *name) const
{
  MEDCouplingFieldDouble *f = nullptr;
  for (size_t i = 0; i < this->fields.size(); i++) {
#if USE_PARAFIELD == 1
    if (strcmp(name, this->fields[i]->getField()->getName().c_str()) == 0) {
      f = this->fields[i]->getField();
      break;
    }
#else
    if (strcmp(name, this->fields[i]->getName().c_str()) == 0) {
      f = this->fields[i];
      break;
    }
#endif
  }

  if (f == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: Could not find field '%s'."),
              name);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generate mesh structure from user's defintion
 *
 * \param[in] select_criteria   selection criteria (string)
 * \param[in] elt_dim           mesh dimension (2 or 3)
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::_generate_coupling_mesh(const char *select_criteria,
                                                   int         elt_dim)
{
  cs_mesh_t *parent_mesh = cs_glob_mesh;

  if (this->mesh != nullptr) {
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: mesh is already added '%s'."),
              this->mesh->med_mesh->getName().c_str());
  }

  /* Building the MED representation of the internal mesh */
  int use_bbox = (elt_dim == 3) ? 1 : 0;
  this->mesh   = cs_medcoupling_mesh_from_base(parent_mesh,
                                               "CouplingMesh",
                                               select_criteria,
                                               elt_dim,
                                               use_bbox);

  /* Define associated ParaMESH */

  if (this->dec != nullptr) {
    ProcessorGroup *Grp = this->dec->isInSourceSide()
                            ? this->dec->getSourceGrp()
                            : this->dec->getTargetGrp();

    this->para_mesh = new ParaMESH(this->mesh->med_mesh, *(Grp), "CoupledMesh");
  }

  this->dec_synced = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generate mesh structure from user's defintion
 *
 * \param[in] n_elts   local number of elements
 * \param[in] elt_ids  list of local elements
 * \param[in] elt_dim  dimension of elements (2: faces, 3: cells)
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::_generate_coupling_mesh_from_ids(
  cs_lnum_t       n_elts,
  const cs_lnum_t elt_ids[],
  int             elt_dim)
{
  cs_mesh_t *parent_mesh = cs_glob_mesh;

  if (this->mesh != nullptr) {
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: mesh is already added '%s'."),
              this->mesh->med_mesh->getName().c_str());
  }

  /* Building the MED representation of the internal mesh */
  int use_bbox = (elt_dim == 3) ? 1 : 0;
  this->mesh   = cs_medcoupling_mesh_from_ids(parent_mesh,
                                              this->_name.c_str(),
                                              n_elts,
                                              elt_ids,
                                              elt_dim,
                                              use_bbox);

  /* Define associated ParaMESH */
  if (this->dec != nullptr) {
    ProcessorGroup *Grp = this->dec->isInSourceSide()
                            ? this->dec->getSourceGrp()
                            : this->dec->getTargetGrp();

    this->para_mesh = new ParaMESH(this->mesh->med_mesh, *(Grp), "CoupledMesh");
  }

  this->dec_synced = false;
}

#endif /* HAVE_PARAMEDMEM */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a selection criteria
 *
 * \param[in] sel_crit  geometrical selection criteria (string)
 * \param[in] elt_dim   dimension of coupled elements
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::add_mesh_from_criteria(const char *sel_crit,
                                                  int         elt_dim)
{
#if !defined(HAVE_PARAMEDMEM)

  this->_error_without_paramedmem();

#else

  this->_generate_coupling_mesh(sel_crit, elt_dim);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a cs_zone_t pointer
 *
 * \param[in] zone  pointer to cs_zone_t struct
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::add_mesh_from_zone(const cs_zone_t *zone)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(zone);
  this->_error_without_paramedmem();

#else

  /* Get location id */
  int location_id = zone->location_id;

  /* Get selection criteria from location id */
  const char *sel_crit = cs_mesh_location_get_selection_string(location_id);

  /* Get location id and deduce element dimension */
  int                     elt_dim  = -1;
  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(location_id);
  if (loc_type == CS_MESH_LOCATION_CELLS)
    elt_dim = 3;
  else if (loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
    elt_dim = 2;

  if (elt_dim < 1)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: wrong dimension for considered zone\n"));

  this->_generate_coupling_mesh(sel_crit, elt_dim);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a cs_zone_t pointer
 *
 * \param[in] n_elts   local number of elements
 * \param[in] elt_ids  list of local elements
 * \param[in] elt_dim  dimension of elements (2: faces, 3: cells)
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::add_mesh_from_ids(cs_lnum_t       n_elts,
                                             const cs_lnum_t elt_ids[],
                                             int             elt_dim)
{
#if !defined(HAVE_PARAMEDMEM)

  this->_error_without_paramedmem();

#else

  if (elt_dim < 2 || elt_dim > 3)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: wrong dimension for considered zone\n"));

  this->_generate_coupling_mesh_from_ids(n_elts, elt_ids, elt_dim);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a coupled field
 *
 * \param[in] name          name of field
 * \param[in] dim           field dimension
 * \param[in] field_nature  field nature flag
 * \param[in] space_discr   field space discretisation (nodes or cells)
 * \param[in] time_discr    field coupling time discretisation
 *
 * \return index of field within the storing vector
 */
/*----------------------------------------------------------------------------*/

int
_cs_paramedmem_coupling_t::add_field(const char              *name,
                                     int                      dim,
                                     cs_medcpl_field_nature_t field_nature,
                                     cs_medcpl_space_discr_t  space_discr,
                                     cs_medcpl_time_discr_t   time_discr)
{
  int f_id = -1;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(dim);
  CS_NO_WARN_IF_UNUSED(field_nature);
  CS_NO_WARN_IF_UNUSED(space_discr);
  CS_NO_WARN_IF_UNUSED(time_discr);

  this->_error_without_paramedmem();

#else

  /* Prepare coupling structure */
  TypeOfField type;
  switch (space_discr) {
    case CS_MEDCPL_ON_CELLS:
      type = ON_CELLS;
      break;

    case CS_MEDCPL_ON_NODES:
      type = ON_NODES;
      break;

    default:
      bft_error(__FILE__,
                __LINE__,
                0,
                _("%s: The type of field is not supported"),
                __func__);
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

    default:
      bft_error(__FILE__,
                __LINE__,
                0,
                _("%s: The type of time is not supported"),
                __func__);
  }

  /* Build ParaFIELD object if required */
  ComponentTopology comp_topo(dim);
#if USE_PARAFIELD == 1
  ParaFIELD *pf = new ParaFIELD(type, td, this->para_mesh, comp_topo);
#else
  MEDCouplingFieldDouble *pf = MEDCouplingFieldDouble::New(type, td);
  pf->setMesh(this->mesh->med_mesh);
  this->mesh->med_mesh->decrRef();
  DataArrayDouble *arr = DataArrayDouble::New();
  int n_elts = (type == ON_CELLS) ? this->mesh->n_elts : get_n_vertices(c);
  arr->alloc(n_elts, dim);
  pf->setArray(arr);
  pf->getArray()->decrRef();
#endif

  if (this->dec != nullptr) {
    if (type == ON_CELLS)
      this->dec->setMethod("P0");
    else
      this->dec->setMethod("P1");
  }

  this->fields.push_back(pf);

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
      bft_error(__FILE__,
                __LINE__,
                0,
                _("%s: The nature of field is not supported"),
                __func__);
  }

#if USE_PARAFIELD == 1
  pf->getField()->setNature(nature);
  pf->getField()->setName(name);
  pf->getField()->getArray()->fillWithZero();
#else
  pf->setNature(nature);
  pf->setName(name);
  pf->getArray()->fillWithZero();
#endif

#endif

  return f_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a coupled field based on a cs_field_t pointer
 *
 * \param[in] f           pointer to cs_field_t struct
 * \param[in] fn          field nature flag
 * \param[in] time_discr  field coupling time discretisation
 *
 * \return index of field within the storing vector
 */
/*----------------------------------------------------------------------------*/

int
_cs_paramedmem_coupling_t::add_field(const cs_field_t        *f,
                                     cs_medcpl_field_nature_t fn,
                                     cs_medcpl_time_discr_t   td)
{
  int f_id = -1;

#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(f);
  CS_NO_WARN_IF_UNUSED(fn);
  CS_NO_WARN_IF_UNUSED(td);

  this->_error_without_paramedmem();

#else

  cs_mesh_location_type_t loc_type = cs_mesh_location_get_type(f->location_id);
  cs_medcpl_space_discr_t sd       = CS_MEDCPL_ON_CELLS;

  if (loc_type == CS_MESH_LOCATION_CELLS ||
      loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
    sd = CS_MEDCPL_ON_CELLS;
  else if (loc_type == CS_MESH_LOCATION_VERTICES)
    sd = CS_MEDCPL_ON_NODES;
  else
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: Non-compatible field location for '%s'\n"),
              f->name);

  f_id = this->add_field(f->name, f->dim, fn, sd, td);
#endif

  return f_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values based on parent mesh location to associated
 *        ParaFIELD objects.
 *
 * \param[in]  name    name of field
 * \param[in]  values  array of values to write
 *                     (defined on parent mesh location)
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::set_values(const char *name, const double values[])
{
#if !defined(HAVE_PARAMEDMEM)

  this->_error_without_paramedmem();

#else
  MEDCouplingFieldDouble *f = this->_get_field(name);

  double *val_ptr = f->getArray()->getPointer();

  /* Assign element values */
  cs_lnum_t *elt_list = nullptr;
  cs_lnum_t  n_elts   = 0;
  if (f->getTypeOfField() == ON_NODES) {
    elt_list = this->mesh->vtx_list;
    n_elts   = this->get_n_vertices();
  }
  else {
    elt_list = this->mesh->elt_list;
    n_elts   = this->get_n_elts();
  }

  if (elt_list == nullptr) {
    const cs_lnum_t n_vals = (cs_lnum_t)f->getNumberOfValues();
    cs_array_copy(n_vals, values, val_ptr);
  }
  else {
    const cs_lnum_t _dim = (cs_lnum_t)f->getNumberOfComponents();
    assert(n_elts * _dim == f->getNumberOfValues());
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t c_id = elt_list[i];
      for (cs_lnum_t j = 0; j < _dim; j++) {
        val_ptr[i * _dim + j] = values[c_id * _dim + j];
      }
    }
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values based on mesh location corresponding to coupled
 *        elements (and associated ParaMESH) to associated ParaFIELD objects.
 *
 * If the whole mesh is coupled, the behavior is the sames as that of
 * \ref set_values.
 *
 * \param[in]  name    name of field
 * \param[in]  values  array of values to write
 *                     (defined on selected mesh subset)
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::set_values_l(const char *name, const double values[])
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(values);

  this->_error_without_paramedmem();

#else

  MEDCouplingFieldDouble *f = this->_get_field(name);

  double         *val_ptr = f->getArray()->getPointer();
  const cs_lnum_t n_vals  = (cs_lnum_t)f->getNumberOfValues();

  /* Assign element values */
  cs_array_copy(n_vals, values, val_ptr);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy values from associated ParaFIELD object to array defined
 *        parent mesh location.
 *
 * \param[in]  name    name of field
 * \param[in]  values  array in which values will be stored
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::get_values(const char *name, double values[]) const
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(values);

  this->_error_without_paramedmem();

#else

  MEDCouplingFieldDouble *f = this->_get_field(name);

  const double *val_ptr = f->getArray()->getConstPointer();

  /* Import element values */

  cs_lnum_t *connec = this->mesh->new_to_old;

  if (connec != nullptr) {
    assert(f->getTypeOfField() == ON_CELLS);
    cs_lnum_t _dim   = (cs_lnum_t)f->getNumberOfComponents();
    cs_lnum_t n_elts = this->mesh->n_elts;
    assert(n_elts * _dim == f->getNumberOfValues());
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t c_id = connec[i];
      for (cs_lnum_t j = 0; j < _dim; j++) {
        values[c_id * _dim + j] = val_ptr[i * _dim + j];
      }
    }
  }
  else {
    const cs_lnum_t n_vals = (cs_lnum_t)f->getNumberOfValues();
    cs_array_copy(n_vals, val_ptr, values);
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy values from associated ParaFIELD structure to array defined
 *        on mesh location corresponding to coupled elements
 *        (and associated ParaMESH).
 *
 * If the whole mesh is coupled, the behavior is the sames as that of
 * \ref get_values.
 *
 * \param[in]  name    name of field
 * \param[in]  values  array in which values will be stored
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::get_values_l(const char *name, double values[]) const
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(values);

  this->_error_without_paramedmem();

#else

  MEDCouplingFieldDouble *f = this->_get_field(name);

  const double   *val_ptr = f->getArray()->getConstPointer();
  const cs_lnum_t n_vals  = (cs_lnum_t)f->getNumberOfValues();

  /* Import element values */

  cs_array_copy(n_vals, val_ptr, values);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sync the coupling's InterpKernelDEC
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::sync_dec()
{
#if !defined(HAVE_PARAMEDMEM)

  this->_error_without_paramedmem();

#else

  if (!this->dec_synced) {
    if (this->dec != nullptr) {
      this->dec->synchronize();
    }
    this->dec_synced = true;
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of field attached to DEC
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::send_data() const
{
#if !defined(HAVE_PARAMEDMEM)

  this->_error_without_paramedmem();

#else

  if (this->dec != nullptr) {
    this->dec->sendData();
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of a field. If vals pointer is non-null,
 * values are updated before send
 *
 * \param[in] name  name of field
 * \param[in] vals  array of values to write
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::send_data(const char *name, const double *vals)
{
#if defined(HAVE_PARAMEDMEM)

  /* If provided, export data to DEC */
  if (vals != nullptr) {
    this->set_values(name, vals);
  }

  /* Attach field to DEC for sending */
  this->attach_field_by_name(name);

  /* Send data */
  this->sync_dec();
  this->send_data();

#else

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(vals);

  this->_error_without_paramedmem();

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of a field. If vals pointer is non-null,
 * values are updated before send
 *
 * \param[in] name  name of field
 * \param[in] vals  array of values to write
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::send_data_l(const char *name, const double *vals)
{
#if defined(HAVE_PARAMEDMEM)

  /* If provided, export data to DEC */
  if (vals != nullptr) {
    this->set_values_l(name, vals);
  }

  /* Attach field to DEC for sending */
  this->attach_field_by_name(name);

  /* Send data */
  this->sync_dec();
  this->send_data();

#else

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(vals);

  this->_error_without_paramedmem();

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of field attached to DEC
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::recv_data()
{
#if !defined(HAVE_PARAMEDMEM)

  this->_error_without_paramedmem();

#else

  if (this->dec != nullptr) {
    this->dec->recvData();
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of a field.
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field
 * \param[in] vals  array of values to read
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::recv_data(const char *name, double *vals)
{
#if defined(HAVE_PARAMEDMEM)

  /* Attach field to DEC for receiving */
  this->attach_field_by_name(name);

  /* Recieve data */
  this->sync_dec();
  this->recv_data();

  /* Read values */
  this->get_values(name, vals);

#else

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(vals);

  this->_error_without_paramedmem();
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of a field.
 *
 * \param[in] name  name of field
 * \param[in] vals  array of values to read
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::recv_data_l(const char *name, double *vals)
{
#if defined(HAVE_PARAMEDMEM)

  /* Attach field to DEC for receiving */
  this->attach_field_by_name(name);

  /* Recieve data */
  this->sync_dec();
  this->recv_data();

  /* Read values */
  this->get_values_l(name, vals);

#else

  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(vals);

  this->_error_without_paramedmem();

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a field to the DEC for send operation using its index
 *
 * \param[in] field_id  index of field in storing vector
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::attach_field_by_id(int field_id)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(field_id);
  this->_error_without_paramedmem();

#else

  this->_attachLocalField(this->fields[field_id]);

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a field to the DEC for send operation using its name
 *
 * \param[in] name  name of field (string)
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::attach_field_by_name(const char *name)
{
#if !defined(HAVE_PARAMEDMEM)

  CS_NO_WARN_IF_UNUSED(name);
  this->_error_without_paramedmem();

#else

#if USE_PARAFIELD == 1

  ParaFIELD *pf = nullptr;

  for (size_t i = 0; i < this->fields.size(); i++) {
    if (strcmp(name, this->fields[i]->getField()->getName().c_str()) == 0) {
      pf = this->fields[i];
      break;
    }
  }

  if (pf == nullptr) {
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: Could not find field '%s'\n"),
              name);
  }

  this->_attachLocalField(pf);

#else

  MEDCouplingFieldDouble *f = nullptr;

  for (size_t i = 0; i < this->fields.size(); i++) {
    if (strcmp(name, this->fields[i]->getName().c_str()) == 0) {
      f = this->fields[i];
      break;
    }
  }

  if (f == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: Could not find field '%s'\n"),
              name);

  this->_attachLocalField(f);

#endif

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log ParaMEDMEM coupling setup information
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_paramedmem_coupling_t::log() const
{
#if defined(HAVE_PARAMEDMEM)

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "  Coupling: %s\n"
                  "    App1: %s\n"
                  "    App2: %s\n"),
                this->getName().c_str(),
                this->apps[0].app_name,
                this->apps[1].app_name);

  if (this->mesh->elt_dim == 3)
    cs_log_printf(CS_LOG_SETUP, _("    Type: Volume coupling\n"));
  else if (this->mesh->elt_dim == 2)
    cs_log_printf(CS_LOG_SETUP, _("    Type: Boundary coupling\n"));

#endif
};

#endif // cplusplus

/*----------------------------------------------------------------------------*/
