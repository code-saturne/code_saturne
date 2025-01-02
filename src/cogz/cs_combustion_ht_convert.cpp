/*============================================================================
 * Gas combustion model: enthaly to and from temperature conversion.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions.h"
#include "cogz/cs_combustion_boundary_conditions.h"
#include "cogz/cs_combustion_gas.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_ht_convert.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal ht_convert.cpp
        Enthalpy to and from temperature conversion for coal combustion.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for Fortran subroutines
 *============================================================================*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Convert enthalpy to temperature at boundary for gas combustion.
 *
 * \param[in]   h    enthalpy at boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_ht_convert_h_to_t_faces(const cs_real_t  h[],
                                      cs_real_t        t[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_t  n_b_faces = m->n_b_faces;

  const int *pm_flag = cs_glob_physical_model_flag;

  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0) {

    const cs_real_t *bym1 = cs_field_by_name("boundary_ym_fuel")->val;
    const cs_real_t *bym2 = cs_field_by_name("boundary_ym_oxydizer")->val;
    const cs_real_t *bym3 = cs_field_by_name("boundary_ym_product")->val;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_real_t hbl = h[face_id];

      cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
      coefg[0] = bym1[face_id];
      coefg[1] = bym2[face_id];
      coefg[2] = bym3[face_id];
      for (int i = 3; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
        coefg[0] = 0;

      t[face_id] = cs_gas_combustion_h_to_t(coefg, hbl);
    }

  }
  else if (pm_flag[CS_COMBUSTION_SLFM] >= 0) {

    const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
    const cs_real_t *cpro_temp = CS_F_(t)->val;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_lnum_t cell_id = b_face_cells[face_id];
      t[face_id] = cpro_temp[cell_id];
    }

    /* Specific values at inlets */

    assert(cs_glob_boundaries != NULL);
    const cs_boundary_t *bdy = cs_glob_boundaries;

    for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

      if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
        continue;

      const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

      auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>
                  (cs_boundary_conditions_get_model_inlet(z));

      if (ci->ientfu == 1 || ci->ientox == 1) {
        const cs_real_t t_b_z = (ci->ientfu == 1) ? cm->tinfue : cm->tinoxy;

        const cs_lnum_t n_elts = z->n_elts;
        const cs_lnum_t *elt_ids = z->elt_ids;

        for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
          cs_lnum_t face_id = elt_ids[elt_idx];
          t[face_id] = t_b_z;
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at selected boundary faces
 *        for gas combustion.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   n_faces   number of selected boundary faces
 * \param[in]   face_ids  list of associated face ids
 * \param[in]   t         temperature values
 * \param[out]  h         enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_ht_convert_t_to_h_faces_l(cs_lnum_t        n_faces,
                                        const cs_lnum_t  face_ids[],
                                        const cs_real_t  t[],
                                        cs_real_t        h[])
{
  const int *pm_flag = cs_glob_physical_model_flag;

  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0) {

    const cs_real_t *bym1 = cs_field_by_name("boundary_ym_fuel")->val;
    const cs_real_t *bym2 = cs_field_by_name("boundary_ym_oxydizer")->val;
    const cs_real_t *bym3 = cs_field_by_name("boundary_ym_product")->val;

    for (cs_lnum_t face_idx = 0; face_idx < n_faces; face_idx++) {
      cs_lnum_t  face_id = face_ids[face_idx];
      cs_real_t tbl = t[face_id];

      cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
      coefg[0] = bym1[face_id];
      coefg[1] = bym2[face_id];
      coefg[2] = bym3[face_id];
      for (int i = 3; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
        coefg[0] = 0;

      h[face_id] = cs_gas_combustion_t_to_h(coefg, tbl);
    }

  }
  else if (pm_flag[CS_COMBUSTION_SLFM] >= 0) {

    const cs_mesh_t *m = cs_glob_mesh;
    const cs_lnum_t *b_face_cells = m->b_face_cells;
    const cs_lnum_t  n_b_faces = m->n_b_faces;

    const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
    const cs_real_t *cvar_h = CS_F_(h)->val;

    cs_real_t *h_in;
    CS_MALLOC(h_in, n_b_faces, cs_real_t);
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      h_in[face_id] = -HUGE_VALF;
    }

    /* Specific values at inlets */

    assert(cs_glob_boundaries != NULL);
    const cs_boundary_t *bdy = cs_glob_boundaries;

    for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

      if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
        continue;

      const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

      auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>
                  (cs_boundary_conditions_get_model_inlet(z));

      if (ci->ientfu == 1 || ci->ientox == 1) {
        const cs_real_t h_b_z = (ci->ientfu == 1) ? cm->hinfue : cm->hinoxy;

        const cs_lnum_t n_elts = z->n_elts;
        const cs_lnum_t *elt_ids = z->elt_ids;

        for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
          cs_lnum_t face_id = elt_ids[elt_idx];
          h_in[face_id] = h_b_z;
        }
      }

    }

    /* Now get values at selected faces */

    for (cs_lnum_t face_idx = 0; face_idx < n_faces; face_idx++) {
      cs_lnum_t  face_id = face_ids[face_idx];

      if (h_in[face_id] > -HUGE_VALF)
        h[face_id] = h_in[face_id];
      else {
        cs_lnum_t cell_id = b_face_cells[face_id];
        h[face_id] = cvar_h[cell_id];
      }
    }

    CS_FREE(h_in);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy for a given boundary zone with
 *        a gas combustio model, using dense storage for temperature
 *        and enthalpy arrays.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   z  pointer to selected zone.
 * \param[in]   t  temperature values
 * \param[out]  h  enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_ht_convert_t_to_h_faces_z(const cs_zone_t *z,
                                        const cs_real_t  t[],
                                        cs_real_t        h[])
{
  const int *pm_flag = cs_glob_physical_model_flag;

  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0) {

    const cs_real_t *bym1 = cs_field_by_name("boundary_ym_fuel")->val;
    const cs_real_t *bym2 = cs_field_by_name("boundary_ym_oxydizer")->val;
    const cs_real_t *bym3 = cs_field_by_name("boundary_ym_product")->val;

    const cs_lnum_t n_elts = z->n_elts;
    const cs_lnum_t *elt_ids = z->elt_ids;

    for (cs_lnum_t face_idx = 0; face_idx < n_elts; face_idx++) {
      cs_lnum_t  face_id = elt_ids[face_idx];
      cs_real_t tbl = t[face_id];

      cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
      coefg[0] = bym1[face_id];
      coefg[1] = bym2[face_id];
      coefg[2] = bym3[face_id];
      for (int i = 3; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
        coefg[0] = 0;

      h[face_id] = cs_gas_combustion_t_to_h(coefg, tbl);
    }

  }
  else if (pm_flag[CS_COMBUSTION_SLFM] >= 0) {

    const cs_mesh_t *m = cs_glob_mesh;
    const cs_lnum_t *b_face_cells = m->b_face_cells;

    const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    /* Specific values at inlets */

    void *ci_p = cs_boundary_conditions_get_model_inlet_try(z);

    if (ci_p != nullptr) {
      auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>(ci_p);

      if (ci->ientfu == 1 || ci->ientox == 1) {
        const cs_real_t h_b_z = (ci->ientfu == 1) ? cm->hinfue : cm->hinoxy;

        const cs_lnum_t n_elts = z->n_elts;
        const cs_lnum_t *elt_ids = z->elt_ids;

        for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
          cs_lnum_t face_id = elt_ids[elt_idx];
          h[face_id] = h_b_z;
        }

        return;
      }

    }

    /* Other zones */

    {
      const cs_real_t *cvar_h = CS_F_(h)->val;

      const cs_lnum_t n_elts = z->n_elts;
      const cs_lnum_t *elt_ids = z->elt_ids;

      for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
        cs_lnum_t face_id = elt_ids[elt_idx];
        cs_lnum_t cell_id = b_face_cells[face_id];
        h[face_id] = cvar_h[cell_id];
      }
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
