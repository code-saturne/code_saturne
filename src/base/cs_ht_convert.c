/*============================================================================
 * Enthaly to and from temperature conversion.
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

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_elec_model.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ht_convert.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_ht_convert.c
        Enthalpy to and from temperature conversion.

        Other fields may be involved in the conversion.

        TODO: when possible (based on calling functions's conversion to C)
              a function pointer-based logic would allow migrating this
              functionnality to cs_physical_properties, without adding
              high level physical model dependencies (i.e. cs_physical_model.h)
              to that lower-level API.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for Fortran subroutines
 *============================================================================*/

extern void CS_PROCF(coh2tb, COH2TB)
(
  const cs_real_t  *h,
  cs_real_t        *t
);

extern void CS_PROCF(cot2hb, COT2HB)
(
  const cs_lnum_t  *n_faces,
  const cs_lnum_t  *face_ids,
  const cs_real_t  *t,
  cs_real_t        *h
);

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
/*!
 * \brief Convert enthalpy to temperature at all cells.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   h   enthalpy values
 * \param[out]  t   temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_h_to_t_cells(const cs_real_t  h[],
                           cs_real_t        t[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const int *pm_flag = cs_glob_physical_model_flag;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  bool need_solid_default = (mq->has_disable_flag) ? true : false;

  cs_real_t *cpro_t = NULL;

  /* Gas combustion: premixed or diffusion flame */
  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0
      || pm_flag[CS_COMBUSTION_SLFM]>= 0)
    cpro_t = CS_F_(t)->val;

  /* Pulverized coal or fuel combustion */
  else if (   pm_flag[CS_COMBUSTION_COAL] >= 0
           || pm_flag[CS_COMBUSTION_FUEL] >= 0)
    cpro_t = CS_F_(t)->val;

  /* Electric arcs */
  else if (   pm_flag[CS_JOULE_EFFECT] >= 1
           || pm_flag[CS_ELECTRIC_ARCS] >= 1)
    cpro_t = CS_F_(t)->val;

  /* When temperature maintained by model is available
     ------------------------------------------------- */

  if (cpro_t != NULL) {

    for (cs_lnum_t i = 0; i < n_cells; i++)
      t[i] = cpro_t[i];

  }

  /* Default for other cases
     ----------------------- */

  else {

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_cells; i++)
        t[i] = h[i] / cpro_cp[i];
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_cells; i++)
        t[i] = h[i] / cp0;
    }

    need_solid_default = false;

  }

  /* Handle solid zones */

  if (need_solid_default) {

    const int *disable_flag = mq->c_disable_flag;

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        if (disable_flag[i])
          t[i] = h[i] / cpro_cp[i];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        if (disable_flag[i])
          t[i] = h[i] / cp0;
      }
    }

  }

  /* Allow user functions */

  int n_zones = cs_volume_zone_n_zones();
  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    /* Note: we could restrict this to
       z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES
       but the user can also easily handle this */

    cs_user_physical_properties_h_to_t(cs_glob_domain,
                                       z,
                                       false,  /* z_local */
                                       h,
                                       t);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at solid cells only.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_h_to_t_cells_solid(void)
{
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  if (mq->has_disable_flag == 0 || CS_F_(h) == NULL || CS_F_(t) == NULL)
    return;

  const cs_real_t *h = CS_F_(h)->val;
  cs_real_t *t = CS_F_(t)->val;

  int n_zones = cs_volume_zone_n_zones();
  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (   z->type & CS_VOLUME_ZONE_SOLID
        && z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES) {

      const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
      if (f_cp != NULL) {
        const cs_real_t *cpro_cp = f_cp->val;
        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t c_id = z->elt_ids[i];
          t[c_id] = h[c_id] / cpro_cp[c_id];
        }
      }
      else {
        const double cp0 = cs_glob_fluid_properties->cp0;
        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t c_id = z->elt_ids[i];
          t[c_id] = h[c_id] / cp0;
        }
      }

      cs_user_physical_properties_h_to_t(cs_glob_domain,
                                         z,
                                         false,  /* z_local */
                                         h,
                                         t);

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at all boundary faces.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   h   enthalpy values
 * \param[out]  t   temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_ht_convert_h_to_t_faces(const cs_real_t  h[],
                           cs_real_t        t[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t has_dc = mq->has_disable_flag;

  bool need_solid_default = (has_dc) ? true : false;

  const int *pm_flag = cs_glob_physical_model_flag;

  /* Gas combustion: premixed or diffusion flame */
  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0
      || pm_flag[CS_COMBUSTION_SLFM]>= 0)
    CS_PROCF(coh2tb, COH2TB)(h, t);

  /* Pulverized coal combustion */
  else if (pm_flag[CS_COMBUSTION_COAL] >= 0)
    cs_coal_thfieldconv1(CS_MESH_LOCATION_BOUNDARY_FACES, h, t);

  /* Pulverized fuel combustion */
  else if (pm_flag[CS_COMBUSTION_FUEL] >= 0)
    cs_fuel_thfieldconv1(CS_MESH_LOCATION_BOUNDARY_FACES, h, t);

  /* Electric arcs */
  else if (pm_flag[CS_JOULE_EFFECT] < 1 && pm_flag[CS_ELECTRIC_ARCS] >= 1)
    cs_elec_convert_h_to_t_faces(h, t);

  /* Default for other cases
     ----------------------- */

  else {

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        cs_lnum_t c_id = b_face_cells[i];
        t[i] = h[i] / cpro_cp[c_id];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_b_faces; i++)
        t[i] = h[i] / cp0;
    }

    need_solid_default = false;

  }

  /* Default for solid zones
     ----------------------- */

  if (need_solid_default) {

    assert(has_dc == 1);
    const int *disable_flag = mq->c_disable_flag;

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        cs_lnum_t c_id = b_face_cells[i];
        if (disable_flag[c_id])
          t[i] = h[i] / cpro_cp[c_id];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        cs_lnum_t c_id = b_face_cells[i];
        if (disable_flag[c_id])
          t[i] = h[i] / cp0;
      }
    }

  }

  /* Allow user functions */

  int n_zones = cs_boundary_zone_n_zones();
  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);

    cs_user_physical_properties_h_to_t(cs_glob_domain,
                                       z,
                                       false,  /* z_local */
                                       h,
                                       t);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at selected boundary faces.
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
cs_ht_convert_t_to_h_faces_l(cs_lnum_t        n_faces,
                             const cs_lnum_t  face_ids[],
                             const cs_real_t  t[],
                             cs_real_t        h[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t has_dc = mq->has_disable_flag;

  const int *pm_flag = cs_glob_physical_model_flag;

  bool need_solid_default = (has_dc) ? true : false;

  /* Gas combustion: premixed or diffusion flame */
  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0
      || pm_flag[CS_COMBUSTION_SLFM] >= 0)
    CS_PROCF(cot2hb, COT2HB)(&n_faces, face_ids, t, h);

  /* Pulverized coal combustion */
  else if (pm_flag[CS_COMBUSTION_COAL] >= 0)
    cs_coal_bt2h(n_faces, face_ids, t, h);

  /* Pulverized fuel combustion */
  else if (pm_flag[CS_COMBUSTION_FUEL] >= 0)
    cs_fuel_bt2h(n_faces, face_ids, t, h);

  /* Electric arcs */
  else if (pm_flag[CS_JOULE_EFFECT] < 1 && pm_flag[CS_ELECTRIC_ARCS] >= 1)
    cs_elec_convert_t_to_h_faces(n_faces,  face_ids, t, h);

  /* Default for other cases
     ----------------------- */

  else {

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        h[f_id] = t[f_id] * cpro_cp[c_id];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        h[f_id] = t[f_id] * cp0;
      }
    }

    need_solid_default = false;

  }

  /* Default for solid zones
     ----------------------- */

  if (need_solid_default) {

    assert(has_dc == 1);
    const int *disable_flag = mq->c_disable_flag;

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        if (disable_flag[c_id])
          h[f_id] = t[f_id] * cpro_cp[c_id];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        if (disable_flag[c_id])
          h[f_id] = t[f_id] * cp0;
      }
    }

  }

  /* Allow user functions */

  cs_real_t *h_f;
  BFT_MALLOC(h_f, n_b_faces, cs_real_t);
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    h_f[i] = h[i];

  int n_zones = cs_boundary_zone_n_zones();
  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);

    cs_user_physical_properties_t_to_h(cs_glob_domain,
                                       z,
                                       false,  /* z_local */
                                       t,
                                       h_f);
  }

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    cs_lnum_t f_id = face_ids[i];
    h[f_id] = h_f[f_id];
  }

  BFT_FREE(h_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy for a given boundary zone,
 *        using dense storage for temperature and enthalpy arrays.
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
cs_ht_convert_t_to_h_faces_z(const cs_zone_t *z,
                             const cs_real_t  t[],
                             cs_real_t        h[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_t has_dc = mq->has_disable_flag;

  const cs_lnum_t  n_faces = z->n_elts;
  const cs_lnum_t  *face_ids = z->elt_ids;

  const int *pm_flag = cs_glob_physical_model_flag;

  bool need_solid_default = (has_dc) ? true : false;

  cs_real_t *t_b = NULL, *h_b = NULL;

  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0
      || pm_flag[CS_COMBUSTION_SLFM] >= 0
      || pm_flag[CS_COMBUSTION_COAL] >= 0
      || pm_flag[CS_COMBUSTION_FUEL] >= 0
      || (pm_flag[CS_JOULE_EFFECT] < 1 && pm_flag[CS_ELECTRIC_ARCS] >= 1)) {

    BFT_MALLOC(t_b, m->n_b_faces, cs_real_t);
    BFT_MALLOC(h_b, m->n_b_faces, cs_real_t);

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t f_id = face_ids[i];
      t_b[f_id] = t[i];
    }
  }

  /* Gas combustion: premixed or diffusion flame */
  if (   pm_flag[CS_COMBUSTION_EBU] >= 0
      || pm_flag[CS_COMBUSTION_3PT] >= 0
      || pm_flag[CS_COMBUSTION_SLFM] >= 0)
    CS_PROCF(cot2hb, COT2HB)(&n_faces, face_ids, t_b, h_b);

  /* Pulverized coal combustion */
  else if (pm_flag[CS_COMBUSTION_COAL] >= 0)
    cs_coal_bt2h(n_faces, face_ids, t_b, h_b);

  /* Pulverized fuel combustion */
  else if (pm_flag[CS_COMBUSTION_FUEL] >= 0)
    cs_fuel_bt2h(n_faces, face_ids, t_b, h_b);

  /* Electric arcs */
  else if (pm_flag[CS_JOULE_EFFECT] < 1 && pm_flag[CS_ELECTRIC_ARCS] >= 1)
    cs_elec_convert_t_to_h_faces(n_faces,  face_ids, t_b, h_b);

  /* Default for other cases
     ----------------------- */

  else {

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        h[i] = t[i] * cpro_cp[c_id];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        h[i] = t[i] * cp0;
      }
    }

    need_solid_default = false;

  }

  /* Gather values if scattered */

  if (h_b != NULL) {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t f_id = face_ids[i];
      h[i] = h_b[f_id];
    }

    BFT_FREE(t_b);
    BFT_FREE(h_b);
  }

  /* Default for solid zones
     ----------------------- */

  if (need_solid_default) {

    assert(has_dc == 1);
    const int *disable_flag = mq->c_disable_flag;

    const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
    if (f_cp != NULL) {
      const cs_real_t *cpro_cp = f_cp->val;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        if (disable_flag[c_id])
          h[i] = t[i] * cpro_cp[c_id];
      }
    }
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t f_id = face_ids[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        if (disable_flag[c_id])
          h[i] = t[i] * cp0;
      }
    }

  }

  /* Allow user functions */

  cs_user_physical_properties_t_to_h(cs_glob_domain,
                                     z,
                                     true,  /* z_local */
                                     t,
                                     h);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
