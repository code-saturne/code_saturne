/*============================================================================
 * SYRTHES coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#if defined(HAVE_MPI)
#include "cs_coupling.h"
#endif

#include "cs_cf_thermo.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_physical_model.h"
#include "cs_thermal_model.h"
#include "cs_syr4_coupling.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_syr_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print information on yet unmatched SYRTHES couplings.
 *
 * parameters:
 *   n_unmatched    <--  number of unmatched couplings
 *   unmatched_ids  <--  array of unmatched couplings
 *----------------------------------------------------------------------------*/

static void
_print_all_unmatched_syr(int        n_unmatched,
                         const int  unmatched_ids[])
{
  /* Loop on defined SYRTHES instances */

  for (int i = 0; i < n_unmatched; i++) {

    cs_syr4_coupling_t *syr_coupling
      = cs_syr4_coupling_by_id(unmatched_ids[i]);
    const char *local_name = cs_syr4_coupling_get_name(syr_coupling);

    bft_printf(_(" SYRTHES coupling:\n"
                 "   coupling id:              %d\n"
                 "   local name:               \"%s\"\n\n"),
               i, local_name);

  }

  bft_printf_flush();
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize MPI SYRTHES couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   n_unmatched    <->  pointer to number of unmatched couplings
 *   unmatched_ids  <->  pointer to array of unmatched couplings
 *----------------------------------------------------------------------------*/

static void
_init_all_mpi_syr(int  *n_unmatched,
                  int  **unmatched_ids)
{
  int _n_unmatched = *n_unmatched;
  int *_unmatched_ids = *unmatched_ids;

  const int n_couplings = cs_syr4_coupling_n_couplings();

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* Loop on applications */

  for (int i = 0; i < n_apps; i++) {

    ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(mpi_apps, i);

    if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0) {

      int  match_queue_id = -1;
      int  coupling_id = -1;

      if (n_apps == 2 && n_couplings == 1 && _n_unmatched == 1) {
        match_queue_id = 0;
        coupling_id = 0;
      }
      else if (ai.app_name != NULL) {
        for (int j = 0; j < _n_unmatched; j++) {
          int k = _unmatched_ids[j];
          cs_syr4_coupling_t *scpl = cs_syr4_coupling_by_id(k);
          if (strcmp(ai.app_name, cs_syr4_coupling_get_name(scpl)) == 0) {
            coupling_id = k;
            match_queue_id = j;
            break;
          }
        }
      }

      if (coupling_id > -1) {

        /* Remove from unmatched queue */
        _n_unmatched -= 1;
        for (int l = match_queue_id; l < _n_unmatched; l++)
          _unmatched_ids[l] = _unmatched_ids[l+1];
        if (_n_unmatched == 0)
          BFT_FREE(_unmatched_ids);

        /* Set communicator */
        cs_syr4_coupling_init_comm(cs_syr4_coupling_by_id(coupling_id),
                                   coupling_id,
                                   ai.root_rank,
                                   ai.n_ranks);

        /* Print matching info */

        const char *syr_version = cs_empty_string;
        const char *local_name = cs_empty_string;
        const char *distant_name = cs_empty_string;

        if (ai.app_name != NULL)
          local_name = ai.app_name;
        if (ai.app_type != NULL)
          syr_version = ai.app_type;
        if (ai.app_name != NULL)
          distant_name = ai.app_name;

        bft_printf(_(" SYRTHES coupling:\n"
                     "   coupling id:              %d\n"
                     "   version:                  \"%s\"\n"
                     "   local name:               \"%s\"\n"
                     "   distant application name: \"%s\"\n"
                     "   MPI application id:       %d\n"
                     "   MPI root rank:            %d\n"
                     "   number of MPI ranks:      %d\n\n"),
                   coupling_id, syr_version, local_name, distant_name,
                   i, ai.root_rank, ai.n_ranks);
      }

      /* Note that a SYRTHES app may be present in the coupling set, but
         not coupled to the current code_saturne instance, so
         coupling_id < 0 here should not be reported as an error or
         complained about here. In case if missing matches, only the
         codes having defined and missing couplings should complain. */

    }

  } /* End of loop on applications */

  bft_printf_flush();

  /* Set return values */

  *n_unmatched = _n_unmatched;
  *unmatched_ids = _unmatched_ids;
}

/*----------------------------------------------------------------------------
 * Find name of single SYRTHES coupling using MPI.
 *
 * If no coupling or multiple couplings are present, the default cannot be
 * determined, so NULL is returned.
 *----------------------------------------------------------------------------*/

static const char *
_mpi_syr_default_name(void)
{
  const char *retval = NULL;

  int n_syr4_apps = 0;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return NULL;

  int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* First pass to count available SYRTHES couplings */

  for (int i = 0; i < n_apps; i++) {
    const ple_coupling_mpi_set_info_t
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0) {
      if (n_syr4_apps == 0)
        retval = ai.app_name;
      else
        retval = NULL;
      n_syr4_apps += 1;
    }
  }

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define new SYRTHES coupling.
 *
 * \param[in] syrthes_name      matching SYRTHES application name
 * \param[in] boundary_criteria surface selection criteria, or NULL
 * \param[in] volume_criteria   volume selection criteria, or NULL
 * \param[in] projection_axis   x', 'y', or 'y' for 2D projection axis (case
 *                              independent), or ' ' for standard 3D coupling
 * \param[in] allow_nonmatching allow nearest-neighbor mapping where matching
 *                              within tolerance is not available (useful
 *                              when meshes have a different level of detail)
 * \param[in] tolerance         addition to local extents of each element
 *                              extent = base_extent * (1 + tolerance)
 * \param[in] verbosity         verbosity level
 * \param[in] visualization     visualization output level (0 or 1)
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * 'syrthes_name' argument is ignored, as there is only one matching
 * possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances based on the 'syrthes_name' argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_define(const char  *syrthes_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       bool         allow_nonmatching,
                       float        tolerance,
                       int          verbosity,
                       int          visualization)
{
  int dim = 3;
  int ref_axis = -1;

  switch (projection_axis) {
  case 'x':
  case 'X':
    dim = 2;
    ref_axis = 0;
    break;
  case 'y':
  case 'Y':
    dim = 2;
    ref_axis = 1;
    break;
  case 'z':
  case 'Z':
    dim = 2;
    ref_axis = 2;
    break;
  default:
    break;
  }

  /* Ensure name is available */

#if defined(HAVE_MPI)
  if (syrthes_name == NULL)
    syrthes_name = _mpi_syr_default_name();
#endif

  if (syrthes_name == NULL)
    syrthes_name = cs_empty_string;

  /* Define additional coupling */

  cs_syr4_coupling_t  *syr_coupling = cs_syr4_coupling_define(dim,
                                                              ref_axis,
                                                              syrthes_name,
                                                              allow_nonmatching,
                                                              tolerance,
                                                              verbosity,
                                                              visualization);

  /* Add locations if done at that stage (deprecated) */

  int n_locations = cs_mesh_location_n_locations();

  const char *sel_criteria[2] = {boundary_criteria, volume_criteria};
  const char *type_name[2] = {"faces", "cells"};
  cs_mesh_location_type_t type_filter[2] = {CS_MESH_LOCATION_BOUNDARY_FACES,
                                            CS_MESH_LOCATION_CELLS};

  for (int i = 0; i < 2; i++) {

    if (sel_criteria[i] != NULL) {
      for (int j = 0; j < n_locations && sel_criteria[i] != NULL; j++) {
        cs_mesh_location_type_t l_type = cs_mesh_location_get_type(j);
        if (l_type & type_filter[i]) {
          const char *c = cs_mesh_location_get_selection_string(j);
          if (c != NULL) {
            if (strcmp(c, sel_criteria[i]) == 0) {
              cs_syr4_coupling_add_location(syr_coupling, j);
              sel_criteria[i] = NULL;
            }
          }
        }
      }
    }

    if (sel_criteria[i] != NULL) {

      char *name;
      size_t l = strlen(syrthes_name) + strlen(type_name[i]) + 2;
      BFT_MALLOC(name, l, char);
      snprintf(name, l, "%s_%s", syrthes_name, type_name[i]);

      int j = cs_mesh_location_add(name, type_filter[i], sel_criteria[i]);

      BFT_FREE(name);

      cs_syr4_coupling_add_location(syr_coupling, j);

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associated a zone to a defined SYRTHES coupling.
 *
 * \param[in] syrthes_name  matching SYRTHES application name
 * \param[in] z             pointer to matching zone
 *                          (boundary or volume)
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_add_zone(const char       *syrthes_name,
                         const cs_zone_t  *z)
{
  /* Ensure name is available */

#if defined(HAVE_MPI)
  if (syrthes_name == NULL)
    syrthes_name = _mpi_syr_default_name();
#endif

  if (syrthes_name == NULL)
    syrthes_name = cs_empty_string;

  /* Search for matching name in existing couplings */

  int n_couplings = cs_syr4_coupling_n_couplings();
  bool match = false;

  for (int i = 0; i < n_couplings; i++) {

    cs_syr4_coupling_t  *syr_coupling = cs_syr4_coupling_by_id(i);
    const char *cmp_name = cs_syr4_coupling_get_name(syr_coupling);

    if (strcmp(syrthes_name, cmp_name) == 0) {
      cs_syr4_coupling_add_location(syr_coupling, z->location_id);
      match = true;
      break;
    }

  }

  if (match == false)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: no defined SYRTHES coupling named \"%s\"."),
              __func__, syrthes_name);
}

/*----------------------------------------------------------------------------
 * Initialize SYRTHES couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_init(void)
{
  int n_unmatched = cs_syr4_coupling_n_couplings();

  int *unmatched_ids;
  BFT_MALLOC(unmatched_ids, n_unmatched, int);

  for (int i = 0; i < n_unmatched; i++)
    unmatched_ids[i] = i;

  /* First try using MPI */

#if defined(HAVE_MPI)

  if (n_unmatched > 0)
    _init_all_mpi_syr(&n_unmatched, &unmatched_ids);

#endif

  if (n_unmatched > 0) {

    bft_printf("Unmatched SYRTHES couplings:\n"
               "----------------------------\n\n");

    _print_all_unmatched_syr(n_unmatched, unmatched_ids);

    BFT_FREE(unmatched_ids);

    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 SYRTHES coupling was defined for which\n"
                "no communication with a SYRTHES instance is possible."));
  }
}

/*----------------------------------------------------------------------------
 * Finalize all SYRTHES couplings.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_finalize(void)
{
  cs_syr4_coupling_all_destroy();
}

/*----------------------------------------------------------------------------
 * Return number of SYRTHES couplings.
 *
 * return:
 *   number of SYRTHES couplings defined
 *----------------------------------------------------------------------------*/

int
cs_syr_coupling_n_couplings(void)
{
  return cs_syr4_coupling_n_couplings();
}

/*----------------------------------------------------------------------------
 * Set conservativity forcing flag to True (1) or False (0) for all defined
 * SYRTHES couplings
 *
 * parameter:
 *   flag     <--  Conservativity forcing flag to set
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_conservativity(int  flag)
{
  assert(flag == 0 || flag == 1);
  cs_syr4_coupling_set_conservativity(flag);
}

/*----------------------------------------------------------------------------
 * Set explicit treatment for the source terms in SYRTHES volume couplings
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_explicit_treatment(void)
{
  cs_syr4_coupling_set_explicit_treatment();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log SYRTHES coupling setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_log_setup(void)
{
  /* Get the number of SYRTHES couplings */
  int n_coupl = cs_syr4_coupling_n_couplings();
  const int keysca = cs_field_key_id("scalar_id");
  const int kcpsyr = cs_field_key_id("syrthes_coupling");
  int icpsyr;

  if (n_coupl >= 1) {

    cs_log_printf
      (CS_LOG_SETUP,
       _("SYRTHES coupling\n"
         "----------------\n\n"
         "    number of couplings: %d\n"),
         n_coupl);

    int n_surf_coupl = 0, n_vol_coupl = 0, issurf, isvol;

    for (int coupl_id = 0; coupl_id < n_coupl; coupl_id++) {
      cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(coupl_id);

      /* Add a new surface coupling if detected */
      issurf = cs_syr4_coupling_is_surf(syr_coupling);
      n_surf_coupl += issurf;

      /* Add a new volume coupling if detected */
      isvol = cs_syr4_coupling_is_vol(syr_coupling);
      n_vol_coupl += isvol;
    }

    cs_log_printf
      (CS_LOG_SETUP,
       _("    with             %d surface coupling(s)\n"
         "    with             %d volume coupling(s)\n"),
         n_surf_coupl, n_vol_coupl);

    cs_log_printf
      (CS_LOG_SETUP,
       _("\n"
         "   Coupled scalars\n"
         "------------------------\n"
         " Scalar    Number icpsyr\n"
         "------------------------\n"));


    for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
      cs_field_t  *f = cs_field_by_id(f_id);
      if ((f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_USER)) {
        int ii = cs_field_get_key_int(f, keysca);
        if (ii > 0) {
          icpsyr = cs_field_get_key_int(f, kcpsyr);
          cs_log_printf
            (CS_LOG_SETUP,
             _(" %s %7d %7d\n"),cs_field_get_label(f),ii, icpsyr);
        }
      }
    }
    cs_log_printf
      (CS_LOG_SETUP,
       _("------------------------\n\n"
         "    icpsyr = 0 or 1         (1: scalar coupled to SYRTHES)\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create coupled meshes and setup PLE locator for Syrthes couplings.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_init_meshes(void)
{
  int n_coupl = cs_syr4_coupling_n_couplings();

  for (int coupl_id = 0; coupl_id < n_coupl; coupl_id++) {
    cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(coupl_id);
    cs_syr4_coupling_init_mesh(syr_coupling);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given SYRTHES coupling number is a surface couplings.
 *
 * \param[in] cpl_id   matching SYRTHES coupling id
 *
 * \return 1 if the coupling includes the surface, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_syr_coupling_is_surf(int  cpl_id)
{
  int retval = 0;  /* Default initialization */

  cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id, cs_syr4_coupling_n_couplings());

  retval = cs_syr4_coupling_is_surf(syr_coupling);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read boundary field/variable values relative to a SYRTHES coupling.
 *
 * \param[in]       nvar     number of variables
 * \param[in]       bc_type  boundary condition type
 * \param[in, out]  icodcl   boundary condition codes
 * \param[in, out]  rcodcl   boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_recv_boundary(int        nvar,
                              int        bc_type[],
                              int        icodcl[],
                              cs_real_t  rcodcl[])
{
  /* SYRTHES coupling: get wall temperature
     ====================================== */

  const int kcpsyr = cs_field_key_id("syrthes_coupling");

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_vars = nvar;  /* cast to cs_lnum_t because this
                                     is used in address computations here */

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Loop on fields, handling only those coupled with Syrthes */

  int n_fields = cs_field_n_fields();
  for (int field_id = 0 ; field_id <  n_fields; field_id++) {

    cs_field_t  *f = cs_field_by_id(field_id);

    int icpsyr = 0;
    if (f->type & CS_FIELD_VARIABLE)
      icpsyr = cs_field_get_key_int(f, kcpsyr);

    if (icpsyr < 1)
      continue;

    /* Loop on couplings: get wall temperature array for each coupling
       and apply matching boundary condition. */

    for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

      cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

      if (! cs_syr4_coupling_is_surf(syr_coupling))  /* ignore if volume-only */
        continue;

      cs_lnum_t n_cpl_faces = cs_syr4_coupling_get_n_elts(syr_coupling, 0);

      /* Get list of coupled faces */

      cs_lnum_t  *f_ids;
      BFT_MALLOC(f_ids, n_cpl_faces, cs_lnum_t);
      cs_syr4_coupling_get_elt_ids(syr_coupling, f_ids, 0);

      /* Read wall temperature and interpolate if necessary */

      cs_real_t *t_solid;
      BFT_MALLOC(t_solid, n_cpl_faces, cs_real_t);
      cs_syr4_coupling_recv_tsolid(syr_coupling, t_solid, 0);

      /*  For scalars coupled with SYRTHES, prescribe a Dirichlet
          condition at coupled faces.
          For the time being, pass here only once, as only one scalar is
          coupled with SYRTHES.
          For the compressible module, solve in energy, but save the
          temperature separately, for BC's to be clearer. */

      const int k_var_id = cs_field_key_id("variable_id");
      int var_id = cs_field_get_key_int(f, k_var_id) - 1;

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {
        if (f == CS_F_(e_tot)) {
          const cs_field_t *f_t_kelvin = CS_F_(t_kelvin);
          var_id = cs_field_get_key_int(f_t_kelvin, k_var_id);
        }
        else
          bft_error
            (__FILE__, __LINE__, 0,
             _("With the compressible module, only the \"total energy\"\n"
               "scalar field may be coupled with SYRTHES.\n"
               "Here, one tries to couple with the field \"%s\"."),
             f->name);
      }

      int  *_icodcl = icodcl + (var_id*n_b_faces);
      cs_real_t  *_rcodcl1 = rcodcl + (var_id*n_b_faces);
      cs_real_t  *_rcodcl2 = rcodcl + (n_b_faces*n_vars + var_id*n_b_faces);
      cs_real_t  *_rcodcl3 = rcodcl + (2*n_b_faces*n_vars + var_id*n_b_faces);

      for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {

        cs_lnum_t face_id = f_ids[i];

        if (   _icodcl[face_id] != CS_INDEF
            && _icodcl[face_id] != CS_SMOOTHWALL
            && _icodcl[face_id] != CS_ROUGHWALL) {
          if (bc_type[face_id] == CS_SMOOTHWALL)
            _icodcl[face_id] = CS_SMOOTHWALL;
          else if (bc_type[face_id] == CS_ROUGHWALL)
            _icodcl[face_id] = CS_ROUGHWALL;
        }

        _rcodcl1[face_id] = t_solid[i];
        _rcodcl2[face_id] = cs_math_infinite_r;
        _rcodcl3[face_id] = 0.;

      }

      /* Require temperature -> enthalpy conversion */

      if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

        if (f == cs_thermal_model_field()) {
          for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
            cs_lnum_t face_id = f_ids[i];
            _icodcl[face_id] *= -1;
          }
        }

      } /* End case for enthalpy */

      BFT_FREE(f_ids);
      BFT_FREE(t_solid);

    } /* End loop on couplings */

  } /* End loop on fields */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Send field/variable values relative to a SYRTHES coupling.
 *
 * \param[in]  h_wall   wall thermal exchange coefficient
 * \param[in]  v_fluid  near-wall fluid thermal variable
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_send_boundary(const cs_real_t  h_wall[],
                              cs_real_t        v_fluid[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Check if we have a boundary coupling. */

  bool have_boundary_cpl = false;
  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {
    cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);
    if (cs_syr4_coupling_is_surf(syr_coupling)) {
      have_boundary_cpl = true;
      break;
    }
  }

  if (! have_boundary_cpl)
    return;

  /* Build arrays large enough for all cases */

  cs_lnum_t  *f_ids;
  cs_real_t  *t_fluid, *h_cpl;
  BFT_MALLOC(f_ids, n_b_faces, cs_lnum_t);
  BFT_MALLOC(t_fluid, n_b_faces, cs_real_t);
  BFT_MALLOC(h_cpl, n_b_faces, cs_real_t);

  /* Prepare conversion to temperature for enthalpy or energy
     (check for surface couplings to make sure it is needed,
     exit earlier otherwise) */

  cs_real_t  *wa = NULL;
  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {
    BFT_MALLOC(wa, n_b_faces, cs_real_t);
    CS_PROCF(b_h_to_t, B_H_TO_T)(v_fluid, wa);
  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TOTAL_ENERGY) {
    /* Epsilon sup for perfect gas at cells */
    BFT_MALLOC(wa, n_cells, cs_real_t);
    cs_cf_thermo_eps_sup(CS_F_(rho)->val, wa, n_cells);
  }

  /* Loop on couplings */

  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

    cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

    if (! cs_syr4_coupling_is_surf(syr_coupling))  /* ignore if volume-only */
      continue;

    cs_lnum_t n_cpl_faces = cs_syr4_coupling_get_n_elts(syr_coupling, 0);

    /* Get list of coupled faces */

    cs_syr4_coupling_get_elt_ids(syr_coupling, f_ids, 0);

    switch (cs_glob_thermal_model->itherm) {

    case CS_THERMAL_MODEL_TEMPERATURE:
      {
        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];

          /* Saved fluid temperatures and exchange coefficients */
          t_fluid[i] = v_fluid[face_id];
          h_cpl[i] = h_wall[face_id];
        }
      }
      break;

    case CS_THERMAL_MODEL_ENTHALPY:
      {
        /* In enthalpy formulation, transform to temperatures for SYRTHES
         *  To conserve flux Phi = (lambda/d     ) Delta T
         *                 or Phi = (lambda/(d Cp)) Delta H
         * recall      hbord = lambda/d.
         *  Conservation is not guaranteed, so we add a warning. */

        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];

          t_fluid[i] = wa[face_id];
          h_cpl[i] = h_wall[face_id];
        }
      }
      break;

    case CS_THERMAL_MODEL_TOTAL_ENERGY:
      {
        /* In energy formulation, transform to temperatures for SYRTHES
         *  To conserve flux Phi = (lambda/d     ) Delta T
         *                or Phi = (lambda/(d Cp)) Delta H
         *  Recall      hbord = lambda/ d
         *  Note that Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
         *  and  that Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
         *    (the difference is thus Cv Delta T)

         * Modify temperature and exchange coefficient

         * Compute e - CvT */

        const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

        cs_real_t   cv0 = cs_glob_fluid_properties->cv0;
        const cs_real_t  *cv = NULL;
        cs_lnum_t   cv_step = 0;

        const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

        if (CS_F_(cv) != NULL) {
          cv = (const cs_real_t *)CS_F_(cv)->val;
          cv_step = 1;
        }
        else
          cv = &cv0;

        const cs_real_t *cvar_e_tot = (const cs_real_t *)CS_F_(e_tot)->val;

        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];

          cs_real_t cvt =   cvar_e_tot[face_id]
                          - 0.5*cs_math_3_square_norm(cvar_vel[cell_id])
                          + wa[cell_id];

          t_fluid[i] = cvt / cv[cell_id * cv_step];
          h_cpl[i] = h_wall[face_id];
        }
      }
      break;

    default:
      break;
    }

    /* Fluxes are multiplied by porosity if present.
     * Here as the flux is expressed as h.(Tw-Tf), the exchange coefficient
     * is multipled by the porosity. */

    const cs_field_t *f_poro = CS_F_(poro);

    if (f_poro != NULL) {

      const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

      if (f_poro->dim == 1) {
        const cs_real_t * cpro_poro = (const cs_real_t *)f_poro->val;
        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];
          h_cpl[i] *= cpro_poro[cell_id];
        }
      }
      else if (f_poro->dim == 6) {
        const cs_real_6_t * cpro_poro = (const cs_real_6_t *)f_poro->val;
        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];
          /* TODO: using a product between the porosity
             and the boundary face normal would be more precise. */
          h_cpl[i] *= 1./3. * (  cpro_poro[cell_id][0]
                               + cpro_poro[cell_id][1]
                               + cpro_poro[cell_id][2]);
        }
      }

    }

    cs_syr4_coupling_send_tf_hf(syr_coupling, f_ids, t_fluid, h_cpl, 0);

  } /* End loop on couplings */

  BFT_FREE(wa);

  BFT_FREE(f_ids);
  BFT_FREE(t_fluid);
  BFT_FREE(h_cpl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Exchange volume values relative to a SYRTHES coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_exchange_volume(void)
{
  const int kcpsyr = cs_field_key_id("syrthes_coupling");

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Loop on fields, handling only those coupled with Syrthes */

  int n_fields = cs_field_n_fields();
  for (int field_id = 0 ; field_id <  n_fields; field_id++) {

    cs_field_t  *f = cs_field_by_id(field_id);

    int icpsyr = 0;
    if (f->type & CS_FIELD_VARIABLE)
      icpsyr = cs_field_get_key_int(f, kcpsyr);

    if (icpsyr < 1)
      continue;

    /* Sanity check : only temperature is possible when doing a
       volume coupling with SYRTHES */
    if (f != cs_thermal_model_field())
      bft_error
        (__FILE__, __LINE__, 0,
         _("SYRTHES volume coupling possible only with temperature variable,\n"
           "not \"%s\"."),
         f->name);

    /* Loop on couplings: get wall temperature array for each coupling
       and apply matching boundary condition. */

    for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

      cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

      if (! cs_syr4_coupling_is_vol(syr_coupling))  /* ignore if surface-only */
        continue;

      cs_lnum_t n_cpl_cells = cs_syr4_coupling_get_n_elts(syr_coupling, 1);

      /* Get list of coupled cells */

      cs_lnum_t  *c_ids;
      cs_real_t *t_fluid, *h_vol;
      BFT_MALLOC(c_ids, n_cpl_cells, cs_lnum_t);
      BFT_MALLOC(t_fluid, n_cpl_cells, cs_real_t);
      BFT_MALLOC(h_vol, n_cpl_cells, cs_real_t);

      cs_syr4_coupling_get_elt_ids(syr_coupling, c_ids, 1);

      for (cs_lnum_t i = 0; i < n_cpl_cells; i++) {
        h_vol[i] = 0.;
      }

      /* Receive solid temperature.
       * This temperature is stored in a C structure for a future
       * use in source term definition. */

      cs_syr4_coupling_recv_tsolid(syr_coupling, t_fluid, 1);

      const cs_real_t  *cvar_t = (const cs_real_t *)f->val;

      const char  *syrthes_name = cs_syr4_coupling_get_name(syr_coupling);


      cs_user_syrthes_coupling_volume_h(cpl_id,
                                        syrthes_name,
                                        n_cpl_cells,
                                        c_ids,
                                        h_vol);

      for (cs_lnum_t i = 0; i < n_cpl_cells; i++)
        t_fluid[i] = cvar_t[c_ids[i]];

      cs_syr4_coupling_send_tf_hf(syr_coupling, c_ids, t_fluid, h_vol, 1);

      BFT_FREE(c_ids);
      BFT_FREE(t_fluid);
      BFT_FREE(h_vol);

    } /* End loop on couplings */

  } /* End loop on fields */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the source term (implicit and/or explicit part) for a
 *         volume coupling with SYRTHES.
 *
 * \param[in]       field_id  field id
 * \param[in, out]  st_exp    explicit source term
 * \param[in, out]  st_imp    implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_volume_source_terms(int        field_id,
                                    cs_real_t  st_exp[],
                                    cs_real_t  st_imp[])
{
  cs_field_t  *f = cs_field_by_id(field_id);

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Sanity check : only temperature is possible when doing a
     volume coupling with SYRTHES */
  if (f != cs_thermal_model_field())
    bft_error
      (__FILE__, __LINE__, 0,
       _("SYRTHES volume coupling possible only with temperature variable,\n"
         "not \"%s\"."),
       f->name);

  /* Loop on couplings: get wall temperature array for each coupling
     and apply matching boundary condition. */

  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

    cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

    if (! cs_syr4_coupling_is_vol(syr_coupling))  /* ignore if surface-only */
      continue;

    cs_lnum_t n_cpl_cells = cs_syr4_coupling_get_n_elts(syr_coupling, 1);

    /* Get list of coupled cells */

    cs_lnum_t  *c_ids;
    cs_real_t *t_fluid, *ctbimp, *ctbexp;
    BFT_MALLOC(c_ids, n_cpl_cells, cs_lnum_t);
    BFT_MALLOC(t_fluid, n_cpl_cells, cs_real_t);
    BFT_MALLOC(ctbimp, n_cpl_cells, cs_real_t);
    BFT_MALLOC(ctbexp, n_cpl_cells, cs_real_t);

    cs_syr4_coupling_get_elt_ids(syr_coupling, c_ids, 1);

    /* Loop on coupled cells to initialize arrays */

    const cs_real_t *cvara_vart = (const cs_real_t *)f->vals[1];

    for (cs_lnum_t i = 0; i < n_cpl_cells; i++) {
      t_fluid[i] = cvara_vart[c_ids[i]];
      ctbimp[i] = 0.;
      ctbexp[i] = 0.;
    }

    /* Loop on coupled cells to compute crvexp and crvimp */

    for (cs_lnum_t i = 0; i < n_cpl_cells; i++) {

      cs_lnum_t c_id = c_ids[i];

      cs_real_t tsexp = (ctbexp[i] - ctbimp[i]*t_fluid[i]) * cell_f_vol[c_id];
      cs_real_t tsimp =  ctbimp[i] * cell_f_vol[c_id];

      st_exp[c_id] += tsexp;
      st_imp[c_id] += tsimp;

    }

    BFT_FREE(c_ids);
    BFT_FREE(t_fluid);
    BFT_FREE(ctbimp);
    BFT_FREE(ctbexp);

  } /* End loop on couplings */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get number of coupled elements with SYRTHES.
 *
 * \param[in]   cpl_id  coupling id
 * \param[in]   mode    0 for boundary, 1 for volume
 *
 * \return  number of coupled elements for this coupling
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_syr_coupling_n_elts(int  cpl_id,
                       int  mode)
{
  cs_lnum_t retval = 0;

  cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id,cs_syr4_coupling_n_couplings());

  else
    retval = cs_syr4_coupling_get_n_elts(syr_coupling, mode);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get local ids of elements coupled with SYRTHES
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[out]   elt_ids  ids of coupled elements (preallocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_elt_ids(int        cpl_id,
                        int        mode,
                        cs_lnum_t  elt_ids[])
{
  cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id,cs_syr4_coupling_n_couplings());

  else
    cs_syr4_coupling_get_elt_ids(syr_coupling, elt_ids, mode);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Receive coupling variables from SYRTHES.
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[out]   t_solid  solid temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_recv_tsolid(int        cpl_id,
                            int        mode,
                            cs_real_t  t_solid[])
{
  cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id,cs_syr4_coupling_n_couplings());

  else
    cs_syr4_coupling_recv_tsolid(syr_coupling, t_solid, mode);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Send coupling variables to SYRTHES.
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[in]    elt_ids  ids of coupled elements
 * \param[in]    t_fluid  fluid temperature
 * \param[in]    h_fluid  fluid exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_send_tf_hf(int              cpl_id,
                           int              mode,
                           const cs_lnum_t  elt_ids[],
                           cs_real_t        t_fluid[],
                           cs_real_t        h_fluid[])
{
  cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id,cs_syr4_coupling_n_couplings());

  else
    cs_syr4_coupling_send_tf_hf(syr_coupling, elt_ids,
                                t_fluid, h_fluid, mode);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
