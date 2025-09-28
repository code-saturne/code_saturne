/*============================================================================
 * Management of the GUI parameters file: conjugate heat transfer
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_mem.h"
#include "gui/cs_gui_util.h"
#include "gui/cs_gui_boundary_conditions.h"
#include "mesh/cs_mesh.h"
#include "base/cs_parameters.h"
#include "base/cs_syr_coupling.h"
#include "base/cs_tree.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "gui/cs_gui_conjugate_heat_transfer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define SYRTHES couplings based on the GUI-generated setup.
 *
 * This variant applies to code_saturne v7.0 and later setup files.
 */
/*----------------------------------------------------------------------------*/

static void
_cs_gui_syrthes_coupling_legacy(void)
{
  const int *v_i = nullptr;
  const cs_real_t *v_r = nullptr;

  const char path_c[] = "conjugate_heat_transfer/external_coupling";
  cs_tree_node_t *tn_c = cs_tree_find_node(cs_glob_tree, path_c);

  for (cs_tree_node_t *tn = cs_tree_get_node(tn_c, "syrthes");
       tn != nullptr;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *syrthes_name
      = cs_tree_node_get_child_value_str(tn, "syrthes_name");

    v_r = cs_tree_node_get_child_values_real(tn, "tolerance");
    double tolerance = (v_r != nullptr) ? v_r[0] : 0.1;

    v_i = cs_tree_node_get_child_values_int(tn, "verbosity");
    int verbosity = (v_i != nullptr) ? v_i[0] : 0;

    v_i = cs_tree_node_get_child_values_int(tn, "visualization");
    int visualization = (v_i != nullptr) ? v_i[0] : 1;

    char projection_axis = ' ';
    const char *_projection_axis
      = cs_tree_node_get_child_value_str(tn, "projection_axis");
    if (_projection_axis != nullptr) {
      char c = _projection_axis[0];
      if (   c == 'x' || c == 'X'
          || c == 'y' || c == 'Y'
          || c == 'z' || c == 'Z')
        projection_axis = c;
    }

    bool allow_nonmatching = false;
    v_i = cs_tree_node_get_child_values_int(tn, "allow_nonmatching");
    if (v_i != nullptr) {
      if (v_i[0] > 0) allow_nonmatching = true;
    }

    const char *boundary_criteria
      = cs_tree_node_get_child_value_str(tn, "selection_criteria");
    const char *volume_criteria
      = cs_tree_node_get_child_value_str(tn, "volume_criteria");

    cs_syr_coupling_define(syrthes_name,
                           boundary_criteria,
                           volume_criteria,
                           projection_axis,
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           visualization);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define SYRTHES couplings based on the GUI-generated setup.
 *
 * This variant applies to code_saturne v7.0 and later setup files.
 */
/*----------------------------------------------------------------------------*/

static void
_cs_gui_syrthes_coupling_bc(void)
{
  /* First get the CHT parameters node */
  const char path_c[] = "conjugate_heat_transfer/external_coupling";
  cs_tree_node_t *tn_c = cs_tree_find_node(cs_glob_tree, path_c);

  const cs_real_t *v_r
    = cs_tree_node_get_child_values_real(tn_c, "tolerance");
  double tolerance = (v_r != nullptr) ? v_r[0] : 0.1;

  const int *v_e
    = cs_tree_node_get_child_values_int(tn_c, "verbosity");
  int verbosity = (v_e != nullptr) ? v_e[0] : 0;

  const int *v_i
    = cs_tree_node_get_child_values_int(tn_c, "visualization");
  int visualization = (v_i != nullptr) ? v_i[0] : 1;

  char projection_axis = ' ';
  const char *_projection_axis
    = cs_tree_node_get_child_value_str(tn_c, "projection_axis");
  if (_projection_axis != nullptr) {
    char c = _projection_axis[0];
    if (   c == 'x' || c == 'X'
        || c == 'y' || c == 'Y'
        || c == 'z' || c == 'Z')
      projection_axis = c;
  }

  bool allow_nonmatching = false;
  const int *v_n =
    cs_tree_node_get_child_values_int(tn_c, "allow_nonmatching");
  if (v_n != nullptr) {
    if (v_n[0] > 0) allow_nonmatching = true;
  }

  /* Loop on the list of defined coupled Syrthes instances */
  cs_tree_node_t *tn_syr_inst = cs_tree_find_node(tn_c, "syrthes_instances");

  for (cs_tree_node_t *tn_cpl = cs_tree_find_node(tn_syr_inst, "instance");
      tn_cpl != nullptr;
      tn_cpl = cs_tree_node_get_next_of_name(tn_cpl)) {

    const char *syrthes_name = cs_tree_node_get_tag(tn_cpl, "name");

    cs_syr_coupling_define(syrthes_name,
                           nullptr, /* No bnd criteria */
                           nullptr, /* No vol criteria */
                           projection_axis,
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           visualization);

    /* Loop on coupled boundray zones */
    for (cs_tree_node_t *tn_bnd = cs_tree_find_node(tn_cpl, "coupled_boundary");
         tn_bnd != nullptr;
         tn_bnd = cs_tree_node_get_next_of_name(tn_bnd)) {
      const char *zone_name = cs_tree_node_get_tag(tn_bnd, "label");
      assert(zone_name != nullptr);

      const cs_zone_t *z = cs_boundary_zone_by_name(zone_name);
      cs_syr_coupling_add_zone(syrthes_name, z);

    } /* loop on coupled boundaries */
  } /* loop on SYRTHES instances */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define SYRTHES couplings based on the GUI-generated setup.
 *
 * In the case of a single code_saturne and single SYRTHES instance, the
 * syrthes_name argument is ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_syrthes_coupling(void)
{
  /* Check if new style of coupling is available */
  const char path_c[]
    = "conjugate_heat_transfer/external_coupling/syrthes_instances/instance";
  cs_tree_node_t *tn_c = cs_tree_find_node(cs_glob_tree, path_c);

  /* Prioritize BC declaration of SYRTHES */
  if (tn_c != nullptr) {
    _cs_gui_syrthes_coupling_bc();
  }
  else {
    _cs_gui_syrthes_coupling_legacy();
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
