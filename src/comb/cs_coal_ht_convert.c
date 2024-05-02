/*============================================================================
 * Coal combustion model: enthaly to and from temperature conversion.
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

#include "cs_coal.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal_ht_convert.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal ht_convert.c
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
/*!
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        concentrations at cells for coal combustion.
 *
 * \param[in]   location_id     mesh location id (cells or boundary faces)
 * \param[in]   eh              gas enthalpy (j/kg of gaseous mixture)
 * \param[out]  tp              gas temperature in kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_ht_convert_h_to_t_gas(int              location_id,
                              const cs_real_t  eh[],
                              cs_real_t        tp[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  cs_lnum_t n_elts = 0;
  if (location_id == CS_MESH_LOCATION_CELLS)
    n_elts = mesh->n_cells;
  else if (location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
    n_elts = mesh->n_b_faces;
  else {
    bft_error
      (__FILE__, __LINE__, 0,
       _(" %s: called for mesh location %d but only handles locations:\n"
         "   CS_MESH_LOCATION_CELLS:\n"
         "   CS_MESH_LOCATION_BOUNDARY_FACES."), __func__, location_id);
  }

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int ih2s = cm->ih2s -1;
  int ihy =cm->ihy -1;
  int ihcn = cm->ihcn -1;
  int inh3 = cm->inh3 -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int iso2 = cm->iso2 -1;
  int in2 = cm->in2 -1;

  const cs_real_t *fuel1 = cs_field_by_id(cm->iym1[ichx1])->val;
  const cs_real_t *fuel2 = cs_field_by_id(cm->iym1[ichx2])->val;
  const cs_real_t *fuel3 = cs_field_by_id(cm->iym1[ico])->val;
  const cs_real_t *fuel4 = cs_field_by_id(cm->iym1[ih2s])->val;
  const cs_real_t *fuel5 = cs_field_by_id(cm->iym1[ihy])->val;
  const cs_real_t *fuel6 = cs_field_by_id(cm->iym1[ihcn])->val;
  const cs_real_t *fuel7 = cs_field_by_id(cm->iym1[inh3])->val;
  const cs_real_t *oxyd = cs_field_by_id(cm->iym1[io2])->val;
  const cs_real_t *prod1 = cs_field_by_id(cm->iym1[ico2])->val;
  const cs_real_t *prod2 = cs_field_by_id(cm->iym1[ih2o])->val;
  const cs_real_t *prod3 = cs_field_by_id(cm->iym1[iso2])->val;
  const cs_real_t *xiner = cs_field_by_id(cm->iym1[in2])->val;

  // Mass fraction of gas
  const cs_real_t *x1 = cs_field_by_name("x_c")->val;

  const cs_real_t *cvar_f1m[CS_COMBUSTION_MAX_COALS];
  const cs_real_t *cvar_f2m[CS_COMBUSTION_MAX_COALS];

  for (int icha = 0; icha < cm->n_coals; icha++) {
    cvar_f1m[icha] = cs_field_by_id(cm->if1m[icha])->val;
    cvar_f2m[icha] = cs_field_by_id(cm->if2m[icha])->val;
  }

  # pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {

    cs_real_t f1mc[CS_COMBUSTION_MAX_COALS];
    cs_real_t f2mc[CS_COMBUSTION_MAX_COALS];
    cs_real_t xesp[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

    cs_lnum_t elt_id = elt_idx;
    if (location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
      elt_id = b_face_cells[elt_idx];

    /* Precompute quantities independent of interpolation point */

    xesp[ichx1] = fuel1[elt_id];
    xesp[ichx2] = fuel2[elt_id];
    xesp[ico]   = fuel3[elt_id];
    xesp[ih2s]  = fuel4[elt_id];
    xesp[ihy]   = fuel5[elt_id];
    xesp[ihcn]  = fuel6[elt_id];
    xesp[inh3]  = fuel7[elt_id];
    xesp[io2]   = oxyd[elt_id];
    xesp[ico2]  = prod1[elt_id];
    xesp[ih2o]  = prod2[elt_id];
    xesp[iso2]  = prod3[elt_id];
    xesp[in2]   = xiner[elt_id];

    for (int icha = 0; icha < cm->n_coals; icha++) {
      f1mc[icha] = cvar_f1m[icha][elt_id] / x1[elt_id];
      f2mc[icha] = cvar_f2m[icha][elt_id] / x1[elt_id];
    }

    /* Now interpolate values */

    tp[elt_idx]
      = cs_coal_ht_convert_h_to_t_gas_by_yi_f1f2(eh[elt_idx], xesp, f1mc, f2mc);

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        given mass fractions and average f1/f2 for coal combustion.
 *
 * \param[in]  eh            gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 * \param[in]  xesp          mass fraction (yi) of species
 * \param[in]  f1mc          average f1 per coal
 * \param[in]  f2mc          average f2 per coal
 *
 * \return  gas temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_gas_by_yi_f1f2(cs_real_t        eh,
                                         const cs_real_t  xesp[],
                                         const cs_real_t  f1mc[],
                                         const cs_real_t  f2mc[])
{
  cs_real_t  tp = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int ih2s = cm->ih2s -1;
  int ihy = cm->ihy -1;
  int ihcn = cm->ihcn -1;
  int inh3 = cm->inh3 -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int iso2 = cm->iso2 -1;
  int in2 = cm->in2 -1;

  cs_real_t den1[CS_COMBUSTION_MAX_COALS];
  cs_real_t den2[CS_COMBUSTION_MAX_COALS];

  cs_real_t ychx10 = 0, ychx20 = 0;

  /* Precompute quantities independent of interpolation point */

  for (int icha = 0; icha < cm->n_coals; icha++) {

    int ichx1c_icha = cm->ichx1c[icha] -1;
    int ichx2c_icha = cm->ichx2c[icha] -1;

    den1[icha] = 1. / (  cm->a1[icha]*cm->wmole[ichx1c_icha]
                       + cm->b1[icha]*cm->wmole[ico]
                       + cm->c1[icha]*cm->wmole[ih2o]
                       + cm->d1[icha]*cm->wmole[ih2s]
                       + cm->e1[icha]*cm->wmole[ihcn]
                       + cm->f1[icha]*cm->wmole[inh3]);

    ychx10 += den1[icha]*(f1mc[icha]*cm->a1[icha]*cm->wmole[ichx1c_icha]);

    den2[icha] = 1. / (  cm->a2[icha]*cm->wmole[ichx2c_icha]
                       + cm->b2[icha]*cm->wmole[ico]
                       + cm->c2[icha]*cm->wmole[ih2o]
                       + cm->d2[icha]*cm->wmole[ih2s]
                       + cm->e2[icha]*cm->wmole[ihcn]
                       + cm->f2[icha]*cm->wmole[inh3]);

    ychx20 += den2[icha]*(f2mc[icha]*cm->a2[icha]*cm->wmole[ichx2c_icha]);

  }

  /* Calculation of enthalpy of the gaseous species CHx1m
   *                                            and CHx2m at */

  cs_real_t eh0 = -HUGE_VAL;

  for (int i = 0; i < cm->n_tab_points && tp <= -HUGE_VAL; i++) {

    cs_real_t ehchx1 = 0, ehchx2 = 0;

    if (ychx10 > cs_math_epzero) {
      for (int icha = 0; icha < cm->n_coals; icha++) {
        int ichx1c_icha = cm->ichx1c[icha] -1;
        ehchx1 +=   den1[icha]
                  * (  cm->ehgaze[i][ichx1c_icha]
                     * f1mc[icha]
                     * cm->a1[icha]
                     * cm->wmole[ichx1c_icha]);
      }
      ehchx1 /= ychx10;
    }
    else
      ehchx1 = cm->ehgaze[i][ichx1];

    if (ychx20 > cs_math_epzero) {
      for (int icha = 0; icha < cm->n_coals; icha++) {
        int ichx2c_icha = cm->ichx2c[icha] -1;
        ehchx2 +=   den2[icha]
                  * (  cm->ehgaze[i][ichx2c_icha]
                     * f2mc[icha]
                     * cm->a2[icha]
                     * cm->wmole[ichx2c_icha]);
      }
      ehchx2 /= ychx20;
    }
    else
      ehchx2 = cm->ehgaze[i][ichx2];

    cs_real_t eh1 =   xesp[ichx1]*ehchx1
                    + xesp[ichx2]*ehchx2
                    + xesp[ico]  *cm->ehgaze[i][ico]
                    + xesp[ih2s] *cm->ehgaze[i][ih2s]
                    + xesp[ihy]  *cm->ehgaze[i][ihy]
                    + xesp[ihcn] *cm->ehgaze[i][ihcn]
                    + xesp[inh3] *cm->ehgaze[i][inh3]
                    + xesp[io2]  *cm->ehgaze[i][io2]
                    + xesp[ico2] *cm->ehgaze[i][ico2]
                    + xesp[ih2o] *cm->ehgaze[i][ih2o]
                    + xesp[iso2] *cm->ehgaze[i][iso2]
                    + xesp[in2]  *cm->ehgaze[i][in2];

    /* Interpolate, with clipping at bounds */

    if (eh <= eh1) {
      if (i == 0)
        tp = cm->th[0];
      else {
        assert(eh >= eh0);
        tp = cm->th[i-1] + (eh-eh0) * (cm->th[i]-cm->th[i-1]) / (eh1-eh0);
      }
    }
    else if (i == cm->n_tab_points-1) {
      tp = cm->th[i];
    }

    eh0 = eh1;

  } /* loop on interpolation points */

  return tp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas enthalpy from gas temperature and
 *        given mass fractions and average f1/f2 for coal combustion.
 *
 * \param[in]  tp            gas temperature (in kelvin)
 * \param[in]  xesp          mass fraction (yi) of species
 * \param[in]  f1mc          average f1 per coal
 * \param[in]  f2mc          average f2 per coal
 *
 * \return  gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_gas_by_yi_f1f2(cs_real_t        tp,
                                         const cs_real_t  xesp[],
                                         const cs_real_t  f1mc[],
                                         const cs_real_t  f2mc[])
{
  cs_real_t  eh = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int ih2s = cm->ih2s -1;
  int ihy = cm->ihy -1;
  int ihcn = cm->ihcn -1;
  int inh3 = cm->inh3 -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int iso2 = cm->iso2 -1;
  int in2 = cm->in2 -1;

  cs_real_t den1[CS_COMBUSTION_MAX_COALS];
  cs_real_t den2[CS_COMBUSTION_MAX_COALS];

  cs_real_t ychx10 = 0, ychx20 = 0;

  /* Precompute quantities independent of interpolation point */

  for (int icha = 0; icha < cm->n_coals; icha++) {

    int ichx1c_icha = cm->ichx1c[icha] -1;
    int ichx2c_icha = cm->ichx2c[icha] -1;

    den1[icha] = 1. / (  cm->a1[icha]*cm->wmole[ichx1c_icha]
                       + cm->b1[icha]*cm->wmole[ico]
                       + cm->c1[icha]*cm->wmole[ih2o]
                       + cm->d1[icha]*cm->wmole[ih2s]
                       + cm->e1[icha]*cm->wmole[ihcn]
                       + cm->f1[icha]*cm->wmole[inh3]);

    ychx10 += den1[icha]*(f1mc[icha]*cm->a1[icha]*cm->wmole[ichx1c_icha]);

    den2[icha] = 1. / (  cm->a2[icha]*cm->wmole[ichx2c_icha]
                       + cm->b2[icha]*cm->wmole[ico]
                       + cm->c2[icha]*cm->wmole[ih2o]
                       + cm->d2[icha]*cm->wmole[ih2s]
                       + cm->e2[icha]*cm->wmole[ihcn]
                       + cm->f2[icha]*cm->wmole[inh3]);

    ychx20 += den2[icha]*(f2mc[icha]*cm->a2[icha]*cm->wmole[ichx2c_icha]);

  }

  /* Calculation of enthalpy of the gaseous species CHx1m
   *                                            and CHx2m at */

  cs_real_t eh0 = -HUGE_VAL;

  int s_id = 0, e_id = cm->n_tab_points;

  if (tp <= cm->th[0])
    e_id = 1;
  else if (tp >= cm->th[cm->n_tab_points - 1])
    s_id = cm->n_tab_points - 1;
  else {
    for (int i = 1; i < cm->n_tab_points; i++) {
      if (tp <= cm->th[i]) {
        s_id = i-1;
        e_id = i+1;
        break;
      }
    }
  }

  for (int i = s_id; i < e_id && eh <= -HUGE_VAL; i++) {

    cs_real_t ehchx1 = 0, ehchx2 = 0;

    if (ychx10 > cs_math_epzero) {
      for (int icha = 0; icha < cm->n_coals; icha++) {
        int ichx1c_icha = cm->ichx1c[icha] -1;
        ehchx1 +=   den1[icha]
                  * (  cm->ehgaze[i][ichx1c_icha]
                     * f1mc[icha]
                     * cm->a1[icha]
                     * cm->wmole[ichx1c_icha]);
      }
      ehchx1 /= ychx10;
    }
    else
      ehchx1 = cm->ehgaze[i][ichx1];

    if (ychx20 > cs_math_epzero) {
      for (int icha = 0; icha < cm->n_coals; icha++) {
        int ichx2c_icha = cm->ichx2c[icha] -1;
        ehchx2 +=   den2[icha]
                  * (  cm->ehgaze[i][ichx2c_icha]
                     * f2mc[icha]
                     * cm->a2[icha]
                     * cm->wmole[ichx2c_icha]);
      }
      ehchx2 /= ychx20;
    }
    else
      ehchx2 = cm->ehgaze[i][ichx2];

    cs_real_t eh1 =   xesp[ichx1]*ehchx1
                    + xesp[ichx2]*ehchx2
                    + xesp[ico]  *cm->ehgaze[i][ico]
                    + xesp[ih2s] *cm->ehgaze[i][ih2s]
                    + xesp[ihy]  *cm->ehgaze[i][ihy]
                    + xesp[ihcn] *cm->ehgaze[i][ihcn]
                    + xesp[inh3] *cm->ehgaze[i][inh3]
                    + xesp[io2]  *cm->ehgaze[i][io2]
                    + xesp[ico2] *cm->ehgaze[i][ico2]
                    + xesp[ih2o] *cm->ehgaze[i][ih2o]
                    + xesp[iso2] *cm->ehgaze[i][iso2]
                    + xesp[in2]  *cm->ehgaze[i][in2];

    /* Interpolate, with clipping at bounds */

    /* Linear interpolation */
    if (e_id - s_id == 2) {
      if (i == s_id) {
        /* First pass: prepare for second */
        eh0 = eh1;
      }
      else {
        /* Second pass: compute value */
        eh = eh0 + (eh1-eh0) * (tp-cm->th[i-1]) / (cm->th[i]-cm->th[i-1]);
      }
    }

    /* Clipping at lower or upper bound */
    else {
      assert(e_id - s_id == 1);
      eh = eh1;
    }

  } /* loop on interpolation points */

  return eh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        given mass fractions for coal combustion.
 *
 * \param[in]  eh            gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_gas_by_yi(cs_real_t        eh,
                                    const cs_real_t  xesp[])
{
  cs_real_t  tp = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int ih2s = cm->ih2s -1;
  int ihy = cm->ihy -1;
  int ihcn = cm->ihcn -1;
  int inh3 = cm->inh3 -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int iso2 = cm->iso2 -1;
  int in2 = cm->in2 -1;

  /* Calculation of enthalpy of the gaseous species CHx1m
   *                                            and CHx2m at */

  cs_real_t eh0 = -HUGE_VAL;

  for (int i = 0; i < cm->n_tab_points && tp <= -HUGE_VAL; i++) {

    cs_real_t ehchx1 = cm->ehgaze[i][ichx1];
    cs_real_t ehchx2 = cm->ehgaze[i][ichx2];

    cs_real_t eh1 =   xesp[ichx1]*ehchx1
                    + xesp[ichx2]*ehchx2
                    + xesp[ico]  *cm->ehgaze[i][ico]
                    + xesp[ih2s] *cm->ehgaze[i][ih2s]
                    + xesp[ihy]  *cm->ehgaze[i][ihy]
                    + xesp[ihcn] *cm->ehgaze[i][ihcn]
                    + xesp[inh3] *cm->ehgaze[i][inh3]
                    + xesp[io2]  *cm->ehgaze[i][io2]
                    + xesp[ico2] *cm->ehgaze[i][ico2]
                    + xesp[ih2o] *cm->ehgaze[i][ih2o]
                    + xesp[iso2] *cm->ehgaze[i][iso2]
                    + xesp[in2]  *cm->ehgaze[i][in2];

    /* Interpolate, with clipping at bounds */

    if (eh <= eh1) {
      if (i == 0)
        tp = cm->th[0];
      else {
        assert(eh >= eh0);
        tp = cm->th[i-1] + (eh-eh0) * (cm->th[i]-cm->th[i-1]) / (eh1-eh0);
      }
    }
    else if (i == cm->n_tab_points-1) {
      tp = cm->th[i];
    }

    eh0 = eh1;

  } /* loop on interpolation points */

  return tp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas enthalpy from gas temperature and
 *        given mass fractions for coal combustion.
 *
 * \param[in]  tp            gas temperature (in kelvin)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_gas_by_yi(cs_real_t        tp,
                                    const cs_real_t  xesp[])
{
  cs_real_t  eh = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int ih2s = cm->ih2s -1;
  int ihy = cm->ihy -1;
  int ihcn = cm->ihcn -1;
  int inh3 = cm->inh3 -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int iso2 = cm->iso2 -1;
  int in2 = cm->in2 -1;

  /* Calculation of enthalpy of the gaseous species CHx1m
   *                                            and CHx2m at */

  cs_real_t eh0 = -HUGE_VAL;

  int s_id = 0, e_id = cm->n_tab_points;

  if (tp <= cm->th[0])
    e_id = 1;
  else if (tp >= cm->th[cm->n_tab_points - 1])
    s_id = cm->n_tab_points - 1;
  else {
    for (int i = 1; i < cm->n_tab_points; i++) {
      if (tp <= cm->th[i]) {
        s_id = i-1;
        e_id = i+1;
        break;
      }
    }
  }

  for (int i = s_id; i < e_id && eh <= -HUGE_VAL; i++) {

    cs_real_t ehchx1 = cm->ehgaze[i][ichx1];
    cs_real_t ehchx2 = cm->ehgaze[i][ichx2];

    cs_real_t eh1 =   xesp[ichx1]*ehchx1
                    + xesp[ichx2]*ehchx2
                    + xesp[ico]  *cm->ehgaze[i][ico]
                    + xesp[ih2s] *cm->ehgaze[i][ih2s]
                    + xesp[ihy]  *cm->ehgaze[i][ihy]
                    + xesp[ihcn] *cm->ehgaze[i][ihcn]
                    + xesp[inh3] *cm->ehgaze[i][inh3]
                    + xesp[io2]  *cm->ehgaze[i][io2]
                    + xesp[ico2] *cm->ehgaze[i][ico2]
                    + xesp[ih2o] *cm->ehgaze[i][ih2o]
                    + xesp[iso2] *cm->ehgaze[i][iso2]
                    + xesp[in2]  *cm->ehgaze[i][in2];

    /* Interpolate, with clipping at bounds */

    /* Linear interpolation */
    if (e_id - s_id == 2) {
      if (i == s_id) {
        /* First pass: prepare for second */
        eh0 = eh1;
      }
      else {
        /* Second pass: compute value */
        eh = eh0 + (eh1-eh0) * (tp-cm->th[i-1]) / (cm->th[i]-cm->th[i-1]);
      }
    }

    /* Clipping at lower or upper bound */
    else {
      assert(e_id - s_id == 1);
      eh = eh1;
    }

  } /* loop on interpolation points */

  return eh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        given mass fractions for coal combustion with drying.
 *
 * \param[in]  eh            gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_gas_by_yi_with_drying(cs_real_t        eh,
                                                const cs_real_t  xesp[])
{
  /* Remark: this function is very similar to the main (non-drying)
     variant. With the correct combination of values in xesp, and
     zero values for the c, d, e, and f coefficients, the general
     function should provide the same results. So we should check
     if this is not always the case using the drying model,
     in which case we could simply use the general function and remove
     this one. */

  cs_real_t  tp = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int in2 = cm->in2 -1;

  /* Calculation of enthalpy of the gaseous species CHx1m
   *                                            and CHx2m at */

  cs_real_t eh0 = -HUGE_VAL;

  for (int i = 0; i < cm->n_tab_points && tp <= -HUGE_VAL; i++) {

    cs_real_t ehchx1 = cm->ehgaze[i][ichx1];
    cs_real_t ehchx2 = cm->ehgaze[i][ichx2];

    cs_real_t eh1 =   xesp[ichx1]*ehchx1
                    + xesp[ichx2]*ehchx2
                    + xesp[ico]  *cm->ehgaze[i][ico]
                    + xesp[io2]  *cm->ehgaze[i][io2]
                    + xesp[ico2] *cm->ehgaze[i][ico2]
                    + xesp[ih2o] *cm->ehgaze[i][ih2o]
                    + xesp[in2]  *cm->ehgaze[i][in2];

    /* Interpolate, with clipping at bounds */

    if (eh <= eh1) {
      if (i == 0)
        tp = cm->th[0];
      else {
        assert(eh >= eh0);
        tp = cm->th[i-1] + (eh-eh0) * (cm->th[i]-cm->th[i-1]) / (eh1-eh0);
      }
    }
    else if (i == cm->n_tab_points-1) {
      tp = cm->th[i];
    }

    eh0 = eh1;

  } /* loop on interpolation points */

  return tp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas enthalpy from gas temperature and
 *        given mass fractions for coal combustion with drying.
 *
 * \param[in]  tp            gas temperature (in kelvin)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_gas_by_yi_with_drying(cs_real_t        tp,
                                                const cs_real_t  xesp[])
{
  /* Remark: this function is very similar to the main (non-drying)
     variant. With the correct combination of values in xesp, and
     zero values for the c, d, e, and f coefficients, the general
     function should provide the same results. So we should check
     if this is not always the case using the drying model,
     in which case we could simply use the general function and remove
     this one. */

  cs_real_t  eh = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  int ichx1 = cm->ichx1 -1;
  int ichx2 = cm->ichx2 -1;
  int ico = cm->ico -1;
  int io2 = cm->io2 -1;
  int ico2 = cm->ico2 -1;
  int ih2o = cm->ih2o -1;
  int in2 = cm->in2 -1;

  /* Calculation of enthalpy of the gaseous species CHx1m
   *                                            and CHx2m at */

  cs_real_t eh0 = -HUGE_VAL;

  int s_id = 0, e_id = cm->n_tab_points;

  if (tp <= cm->th[0])
    e_id = 1;
  else if (tp >= cm->th[cm->n_tab_points - 1])
    s_id = cm->n_tab_points - 1;
  else {
    for (int i = 1; i < cm->n_tab_points; i++) {
      if (tp <= cm->th[i]) {
        s_id = i-1;
        e_id = i+1;
        break;
      }
    }
  }

  for (int i = s_id; i < e_id && eh <= -HUGE_VAL; i++) {

    cs_real_t ehchx1 = cm->ehgaze[i][ichx1];
    cs_real_t ehchx2 = cm->ehgaze[i][ichx2];

    cs_real_t eh1 =   xesp[ichx1]*ehchx1
                    + xesp[ichx2]*ehchx2
                    + xesp[ico]  *cm->ehgaze[i][ico ]
                    + xesp[io2]  *cm->ehgaze[i][io2]
                    + xesp[ico2] *cm->ehgaze[i][ico2]
                    + xesp[ih2o] *cm->ehgaze[i][ih2o]
                    + xesp[in2]  *cm->ehgaze[i][in2];

    /* Interpolate, with clipping at bounds */

    /* Linear interpolation */
    if (e_id - s_id == 2) {
      if (i == s_id) {
        /* First pass: prepare for second */
        eh0 = eh1;
      }
      else {
        /* Second pass: compute value */
        eh = eh0 + (eh1-eh0) * (tp-cm->th[i-1]) / (cm->th[i]-cm->th[i-1]);
      }
    }

    /* Clipping at lower or upper bound */
    else {
      assert(e_id - s_id == 1);
      eh = eh1;
    }

  } /* loop on interpolation points */

  return eh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the particles temperature from particles enthalpy and
 *        concentrations at cells for coal combustion.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_ht_convert_h_to_t_particles(void)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_coal_model_t  *cm = cs_glob_coal_model;

  cs_lnum_t n_cells = mesh->n_cells;

  const cs_real_t *cpro_temp = cs_field_by_name("temperature")->val;

  const int ihflt2 = 1; // Conversion mode

  /* H2 linear function of T2
     ------------------------ */

  if (ihflt2 == 0) {

    for (int icla = 0; icla < cm->nclacp; icla++) {

      const int icha = cm->ichcor[icla] - 1;
      const cs_real_t h02ch_icha = cm->h02ch[icha];
      const cs_real_t cp2ch_icha = cm->cp2ch[icha];
      const cs_real_t *cvar_h2cl = cs_field_by_id(cm->ih2[icla])->val;
      cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[icla])->val;

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        // FIXME divide by x2
        cpro_temp2[cell_id] =   (cvar_h2cl[cell_id] - h02ch_icha) / cp2ch_icha
                              + cs_coal_trefth;
      } /* Loop on cells */

    } /* Loop on coal classes */

    return;
  }

  /* H2 tabulated
     ------------ */

  for (int icla = 0; icla < cm->nclacp; icla++) {

    const int icha = cm->ichcor[icla] - 1;
    const int ich_icha = cm->ich[icha] - 1;
    const int ick_icha = cm->ick[icha] - 1;
    const int iash_icha = cm->iash[icha] - 1;
    const int iwat_icha = cm->iwat[icha] - 1;
    const cs_real_t xmash_icla = cm->xmash[icla];
    const cs_real_t xmp0_icla = cm->xmp0[icla];

    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[icla])->val;
    const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[icla])->val;
    const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[icla])->val;
    const cs_real_t *cvar_xwtcl = NULL;
    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
      cvar_xwtcl = cs_field_by_id(cm->ixwt[icla])->val;
    const cs_real_t *cvar_h2cl = cs_field_by_id(cm->ih2[icla])->val;
    cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[icla])->val;

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

      const cs_real_t xch  = cvar_xchcl[cell_id];
      const cs_real_t xck  = cvar_xckcl[cell_id];
      const cs_real_t xnp  = cvar_xnpcl[cell_id];
      const cs_real_t xash = xmash_icla*xnp;
      const cs_real_t xwat = (cvar_xwtcl != NULL) ? cvar_xwtcl[cell_id] : 0.;

      const cs_real_t x2 = xch + xck + xash + xwat;
      const cs_real_t xtes = xmp0_icla * xnp;

      if (xtes > cs_coal_epsilon && x2 > cs_coal_epsilon*100) {

        const cs_real_t xch_d_x2  = xch / x2;
        const cs_real_t xck_d_x2  = xck / x2;
        const cs_real_t xash_d_x2 = xash / x2;;
        const cs_real_t xwat_d_x2 = xwat / x2;;
        const cs_real_t h2 = cvar_h2cl[cell_id] / x2;

        cs_real_t eh0 = -HUGE_VAL, t2 = -HUGE_VAL;

        for (int i = 0; i < cm->npoc && t2 <= -HUGE_VAL; i++) {

          cs_real_t eh1 =   xch_d_x2  * cm->ehsoli[i][ich_icha]
                          + xck_d_x2  * cm->ehsoli[i][ick_icha]
                          + xash_d_x2 * cm->ehsoli[i][iash_icha]
                          + xwat_d_x2 * cm->ehsoli[i][iwat_icha];

          /* Interpolate, with clipping at bounds */

          if (h2 <= eh1) {
            if (i == 0)
              t2 = cm->thc[0];
            else {
              assert(h2 >= eh0);
              t2 = cm->thc[i-1] + (h2-eh0) *  (cm->thc[i]-cm->thc[i-1])
                                             / (eh1-eh0);
            }
          }
          else if (i == cm->npoc-1) {
            t2 = cm->thc[i];
          }

          eh0 = eh1;

        }

        cpro_temp2[cell_id] = t2;

      }

      else {
        cpro_temp2[cell_id] = cpro_temp[cell_id];  /* gas mix temperature */
      }

    } /* Loop on cells */

  } /* Loop on coal classes */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the particles temperature from particles enthalpy and
 *        given mass fractions for coal combustion.
 *
 * \remark  Function not called in code, so should probably be removed,
 *          unless useful for advanced postprocessing.
 *
 * \param[in]  enthal     mass enthalpy (\f$ j . kg^{-1} \f$)
 * \param[in]  class_id   class id (0 to n-1)
 * \param[in]  xesp       mass fraction of components
 *                        (size: cm->nsolid)
 * \param[in]  t1         coal inlet/boundary temperature
 *
 * \return   temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_particles_by_yi(cs_real_t        enthal,
                                          int              class_id,
                                          const cs_real_t  xsolid[],
                                          cs_real_t        t1)
{
  cs_real_t  temper = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  const int ihflt2 = 1; // Conversion mode

  /* H2 linear function
     ------------------ */

  if (ihflt2 == 0) {

    temper = (enthal - cm->h02ch[class_id]) / cm->cp2ch[class_id] + cs_coal_trefth;
    return temper;

  }

  /* H2 tabulated
     ------------ */

  cs_real_t x2 = 0;
  for (int i = 0; i < cm->nsolid; i++)
    x2 += xsolid[i];

  if (x2 > cs_coal_epsilon) {

    cs_real_t eh0 = -HUGE_VAL;

    for (int i = 0; i < cm->npoc && temper <= -HUGE_VAL; i++) {

      cs_real_t eh1 = 0.;
      for (int j = 0; j < cm->nsolid; j++)
        eh1 += xsolid[j]*cm->ehsoli[i][j];

      /* Interpolate, with clipping at bounds */

      if (enthal <= eh1) {
        if (i == 0)
          temper = cm->thc[0];
        else {
          assert(enthal >= eh0);
          temper = cm->thc[i-1] + (enthal-eh0) *   (cm->thc[i]-cm->thc[i-1])
                                                 / (eh1-eh0);
        }
      }
      else if (i == cm->npoc-1) {
        temper = cm->thc[i];
      }

      eh0 = eh1;

    } /* loop on interpolation points */

  }
  else
    temper = t1;

  return temper;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the particles enthalpy from particles temperature and
 *        given mass fractions for coal combustion.
 *
 * \param[in]  temper        temperature (in kelvin)
 * \param[in]  class_id      class id (0 to n-1)
 * \param[in]  xesp          mass fraction of components
 *                           (size: cm->nsolid)
 *
 * \return  mass enthalpy (\f$ j . kg^{-1} \f$)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_particles_by_yi(cs_real_t        temper,
                                          int              class_id,
                                          const cs_real_t  xsolid[])
{
  cs_real_t  enthal = -HUGE_VAL;

  const cs_coal_model_t  *cm = cs_glob_coal_model;

  const int ihflt2 = 1; // Conversion mode

  /* H2 linear function
     ------------------ */

  if (ihflt2 == 0) {

    enthal = cm->h02ch[class_id] + cm->cp2ch[class_id]*(temper-cs_coal_trefth);
    return enthal;

  }

  /* H2 tabulated
     ------------ */

  cs_real_t eh0 = -HUGE_VAL;

  int s_id = 0, e_id = cm->npoc;

  if (temper <= cm->thc[0])
    e_id = 1;
  else if (temper >= cm->thc[cm->npoc - 1])
    s_id = cm->npoc - 1;
  else {
    for (int i = 1; i < cm->npoc; i++) {
      if (temper <= cm->thc[i]) {
        s_id = i-1;
        e_id = i+1;
        break;
      }
    }
  }

  for (int i = s_id; i < e_id && enthal <= -HUGE_VAL; i++) {

    cs_real_t eh1 = 0.;
    for (int j = 0; j < cm->nsolid; j++)
      eh1 += xsolid[j]*cm->ehsoli[i][j];

    /* Interpolate, with clipping at bounds */

    /* Linear interpolation */
    if (e_id - s_id == 2) {
      if (i == s_id) {
        /* First pass: prepare for second */
        eh0 = eh1;
      }
      else {
        /* Second pass: compute value */
        enthal = eh0 + (eh1-eh0) * (temper-cm->th[i-1])
                                 / (cm->th[i]-cm->th[i-1]);
      }
    }

    /* Clipping at lower or upper bound */
    else {
      assert(e_id - s_id == 1);
      enthal = eh1;
    }

  } /* loop on interpolation points */

  return enthal;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
