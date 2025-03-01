/*============================================================================
 * Infinitely fast 3-point combustion model
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_physical_constants.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_properties.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"
#include "mesh/cs_mesh.h"
#include "rayt/cs_rad_transfer.h"
#include "turb/cs_turbulence_model.h"

#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_boundary_conditions.h"
#include "cogz/cs_combustion_ht_convert.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_bsh.h"
#include "cogz/cs_combustion_d3p.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_d3p.cpp
        Eddy-Break-Up gas combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of local stoechiometric enthalpy.
 *
 * \param[in]       n_cells    number of cells
 * \param[in]       indpdf     indicator for pdf integration or mean value
 * \param[in]       dirmin     Dirac's peak value at \f$ f_{min} \f$
 * \param[in]       dirmax     Dirac's peak value at \f$ f_{max} \f$
 * \param[in]       fdeb       abscissa of rectangle low boundary
 * \param[in]       ffin       abscissa of rectangle high boundary
 * \param[in]       hrec       rectangle height
 * \param[in]       fm         mean mixture fraction at cell centers
 * \param[in]       hm         mean mixture enthalpy at cell centers
 * \param[in, out]  hstoe      local stoechiometric enthalpy
 */
/*----------------------------------------------------------------------------*/

static void
_d3phst(cs_lnum_t  n_cells,
        int        indpdf[],
        cs_real_t  dirmin[],
        cs_real_t  dirmax[],
        cs_real_t  fdeb[],
        cs_real_t  ffin[],
        cs_real_t  hrec[],
        cs_real_t  fm[],
        cs_real_t  hm[],
        cs_real_t  hstoe[])
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const cs_real_t fsir = cm->fs[0];
  constexpr cs_real_t epsi = 1.e-6;

  cs_gnum_t n1 = 0, n2 = 0;
  cs_real_t hsmin = HUGE_VALF;
  cs_real_t hsmax = -HUGE_VALF;

  bool log_active = cs_log_default_is_active();

  const cs_real_t hinoxy = cm->hinoxy;
  const cs_real_t hinfue = cm->hinfue;

  const int nmaxh = 9;
  const double *hh = cm->hh;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (indpdf[c_id] == 0) {

      /* Determine hstoe without integration
         ----------------------------------- */

      if (fm[c_id] <= fsir && fm[c_id] > epsi) {
        hstoe[c_id] =   (fsir*hm[c_id]+(fm[c_id]-fsir)*hinoxy)
                      / fm[c_id];
      }
      else if (fm[c_id] < (1.-epsi)) {
        hstoe[c_id] =  ((1.-fsir)*hm[c_id]+(fsir-fm[c_id])*hinfue)
                      / (1.-fm[c_id]);
      }

    }
    else {

      /* Determine hstoe with integration
         -------------------------------- */

      cs_real_t hct = dirmin[c_id]*hinoxy+dirmax[c_id]*hinfue;
      cs_real_t hhh = 0.0;
      if (hrec[c_id] > epsi) {
        if (fdeb[c_id] <= fsir) {
          cs_real_t f1 = fdeb[c_id];
          cs_real_t f2 = std::min(fsir, ffin[c_id]);
          hct += hrec[c_id]*(f2-f1)*hinoxy*(2.*fsir-f1-f2)/(2.*fsir);
          hhh += hrec[c_id]*(f2*f2-f1*f1)/(2.*fsir);
        }
        if (ffin[c_id] > fsir) {
          cs_real_t f1 = std::max(fsir, fdeb[c_id]);
          cs_real_t f2 = ffin[c_id];
          hct += hrec[c_id]*(f2-f1)*hinfue*(f2+f1-2.*fsir)/(2.*(1.-fsir));
          hhh += hrec[c_id]*(f2-f1)*(2.-f1-f2)/(2.*(1.-fsir));
        }
        hstoe[c_id] = (hm[c_id]-hct) / hhh;
      }
    }

    // Clipping to hstoea = hh[0] at max
    // Clipping to hstoea = hh[nmaxh-1] at min

    if (hstoe[c_id] > hh[0]) {
      n1 += 1;
      hsmax = std::max(hstoe[c_id], hsmax);
      hstoe[c_id] = hh[0];
    }

    if (hstoe[c_id] < hh[nmaxh-1]) {
      n2 += 1;
      hsmin = std::min(hstoe[c_id], hsmin);
      hstoe[c_id] = hh[nmaxh-1];
    }
  }

  if (log_active) {

    cs_parall_sum_scalars(n1, n2);
    cs_parall_max_scalars(hsmax);
    cs_parall_min_scalars(hsmin);

    if (n1 > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _(" Clip hstoe to max at %lu points\n"
                      "   Max value:  %15.7g\n"
                      "   Clip value: %15.7g\n\n"),
                    (unsigned long)n1, hsmax, hh[0]);
    }
    if (n2 > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _(" Clip hstoe to min at %lu points\n"
                      "   Min value:  %15.7g\n"
                      "   Clip value: %15.7g\n\n"),
                    (unsigned long)n2, hsmin, hh[nmaxh-1]);
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for infinitely fast 3-point chemistry
 *        combustion model, stage 0
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_d3p_fields_init0(void)
{
  // Only when not a restart
  if (cs_restart_present())
    return;

  // Local variables

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  cs_real_t t0 = cs_glob_fluid_properties->t0;

  assert(cm->type / 100 == 1);

  cs_real_t *cvar_scalt = nullptr;
  if (cm->type % 2 == 1)
    cvar_scalt = CS_F_(h)->val;

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int igg = 0; igg < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; igg++)
    coefg[igg] = 0;

  // Initialize fields
  // -----------------

  // No need to set fm or fp2m to 0 here,
  // as this is the default for all fields.

  if (cvar_scalt != nullptr) {
    cs_real_t h_air;

    if (cm->type % 100 > 1) {  // Burke-Schumann model
      cs_real_t yg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
      cs_real_t ye[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
      for (int igg = 0; igg < cm->n_gas_species; igg++)
        yg[igg] = 0;

      yg[0] = 0.;
      yg[1] = 1.;
      yg[2] = 0.;

      cs_real_t ytot = 0.;
      for (int ii = 0; ii < cm->n_gas_el_comp - 1; ii++) {
        ye[ii] = 0.;
        for (int igg = 0; igg < cm->n_gas_species; igg++) {
          ye[ii] += cm->coefeg[igg][ii] * yg[igg];
        }
        ytot += ye[ii];
      }

      ye[cm->n_gas_el_comp - 1] = 1. - ytot;

      h_air = cs_compute_burke_schumann_enthalpy(t0, ye);
    }
    else {
      coefg[0] = 0.;
      coefg[1] = 1.;
      coefg[2] = 0.;
      h_air = cs_gas_combustion_t_to_h(coefg, t0);
    }

    cs_array_real_set_scalar(n_cells_ext, h_air, cvar_scalt);
  }

  if (cm->type % 100 < 2) { // Legacy, not Burke-Schumann model
    // User initialization, HINFUE and HINOXY are needed

    // Oxydant enthalpy hinoxy at tinoxy.
    coefg[0] = 0.;
    coefg[1] = 1.;
    coefg[2] = 0.;
    cm->hinoxy = cs_gas_combustion_t_to_h(coefg, cm->tinoxy);

    // Fuel enthalpy hinfue at tinfue.
    coefg[0] = 1.;
    coefg[1] = 0.;
    coefg[2] = 0.;
    cm->hinfue = cs_gas_combustion_t_to_h(coefg, cm->tinfue);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for infinitely fast 3-point chemistry
 *        combustion model, stage 1.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_d3p_fields_init1(void)
{
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  assert(cm->type / 100 == 1);

  // Determine thermochemical quantities
  // -----------------------------------

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int igg = 0; igg < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; igg++)
    coefg[igg] = 0;

  const double *fs = cm->fs;

  cm->hstoea = fs[0]*cm->hinfue + (1.-fs[0])*cm->hinoxy;

  // Construction of an enthalpy table based on richness and
  // stoichiometrc enthalpy, of dimension 9x9.

  const int nmaxh = 9, nmaxf = 9;

  double fsir = fs[0];

  // Compute ff array
  for (int i_f = 0; i_f < nmaxf/2+1; i_f++) {
    cm->ff[i_f] = fsir * (double)(2*i_f)/(double)(nmaxf-1);
  }
  for (int i_f = nmaxf/2+1; i_f < nmaxf; i_f++) {
    cm->ff[i_f] = fsir +   (double)(2*i_f+2-nmaxf-1)
                         / (double)(nmaxf-1)*(1.-fsir);
  }

  // Compute hh array
  coefg[0] = 0.;
  coefg[1] = 0.;
  coefg[2] = 1.;
  cs_real_t tin = std::min(cm->tinfue, cm->tinoxy);
  cm->hh[nmaxh-1] = cs_gas_combustion_t_to_h(coefg, tin);
  cm->hh[0] = cm->hstoea;
  for (int i_h = 1; i_h < nmaxh-1; i_h++) {
    cm->hh[i_h] = cm->hh[0] +   (cm->hh[nmaxh-1] - cm->hh[0])
                              * (double)(i_h)/double(nmaxh-1);
  }

  // Compute tfh[i_h][i_f]

  for (int i_h = 0; i_h < nmaxh; i_h++) {
    for (int i_f = 0; i_f < nmaxf/2+1; i_f++) {
      // Poor mixture
      coefg[0] = 0.;
      coefg[1] = (fsir - cm->ff[i_f]) / fsir;
      coefg[2] = cm->ff[i_f] / fsir;
      cs_real_t hhloc = cm->hinoxy +   (double)(2*i_f)/(double)(nmaxf-1)
                                     * (cm->hh[i_h] - cm->hinoxy);
      cm->tfh[i_h][i_f] = cs_gas_combustion_h_to_t(coefg, hhloc);
    }
    for (int i_f = nmaxf/2+1; i_f < nmaxf; i_f++) {
      // Rich mixture
      coefg[0] = (cm->ff[i_f] - fsir) / (1.-fsir);
      coefg[1] = 0.;
      coefg[2] = (1.-cm->ff[i_f]) / (1.-fsir);
      cs_real_t hhloc =   (  (double)(2*i_f+2)*(cm->hinfue-cm->hh[i_h])
                           + (double)(2*nmaxf)*cm->hh[i_h]
                           - cm->hinfue*(double)(nmaxf+1))
                        / (double)(nmaxf-1);
      cm->tfh[i_h][i_f] = cs_gas_combustion_h_to_t(coefg, hhloc);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Integration of thermodynamic variablesfunction of mixture fraction.
 *
 * Remark: temperature integration could be weighted by Cp.
 *
 * \param[in]     indpdf        indicator for pdf integration or mean value
 * \param[in]     dirmin        Dirac's peak value at \f$ f_{min} \f$
 * \param[in]     dirmax        Dirac's peak value at \f$ f_{max} \f$
 * \param[in]     fdeb          abscissa of rectangle low boundary
 * \param[in]     ffin          abscissa of rectangle high boundary
 * \param[in]     hrec          rectangle height
 * \param[in]     w1            work array
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_d3p_integration(int        indpdf[],
                              cs_real_t  dirmin[],
                              cs_real_t  dirmax[],
                              cs_real_t  fdeb[],
                              cs_real_t  ffin[],
                              cs_real_t  hrec[],
                              cs_real_t  w1[])
{
  static int n_calls = 0;  // Count calls
  n_calls++;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const int sub_type = cm->type % 100;

  const int n_gas_species = cm->n_gas_species;

  const cs_real_t tinfue = cm->tinfue;
  const cs_real_t tinfue3 = cs_math_pow3(cm->tinfue);
  const cs_real_t tinfue4 = cs_math_pow4(cm->tinfue);
  const cs_real_t tinoxy = cm->tinoxy;
  const cs_real_t tinoxy3 = cs_math_pow3(cm->tinoxy);
  const cs_real_t tinoxy4 = cs_math_pow4(cm->tinoxy);
  const double *wmolg = cm->wmolg;

  cs_real_t *cpro_csca[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
    cpro_csca[i] = nullptr;

  for (int icg = 0; icg < n_gas_species; icg++)
    cpro_csca[icg] = cm->ym[icg]->val;

  cs_real_t *cvar_fm = cm->fm->val;
  cs_real_t *cvar_scalt = nullptr;
  if (cm->type %2 == 1)
    cvar_scalt = CS_F_(h)->val;

  bool update_rho = false;
  if (n_calls > 1 || cs_restart_get_field_read_status(CS_F_(rho)->id) == 1)
    update_rho = true;

  const cs_real_t pther = cs_glob_fluid_properties->pther;
  const cs_real_t srrom = cm->srrom;

  const cs_real_t epzero = cs_math_epzero;
  const int nmaxf = 9, nmaxh = 9;
  const double *ff = cm->ff;
  const double *hh = cm->hh;
  double tfh[9][9];
  memcpy(tfh, cm->tfh, sizeof(cm->tfh));
  const double *ckabsg = cm->ckabsg;

  /* Integrate n_gas_species mass fractions of global species
   * -------------------------------------------------------- */

  // For 3-point chemistry diffusion flame:
  // - There is only one global reaction (ir).
  // - The mixture fraction varies between 0 and 1.
  cs_real_t fsir = cm->fs[0];
  cs_real_t fmini = 0.;
  cs_real_t fmaxi = 1.;

  cs_host_context ctx;

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    cs_real_t aa1 = 0., aa2 = 0., bb1 = 0., bb2 = 0.;

    cs_real_t fm   = cvar_fm[c_id];

    for (int icg = 0; icg < n_gas_species; icg++) {

      // Determine parameters of variable lines (f)
      //  Y = A + B F
      // With a diffusion flame, there is a single global reaction (ir=1).
      // By definition the mass fractions of global species are then

      //   F       0  FS  1
      //   Yfuel   0  0   1
      //   Yoxyd   1  0   0
      //   Yprod   0  1   0

      if (icg == 0) {
        // Fuel
        aa1 =  0.;
        bb1 =  0.;
        aa2 = -fsir/(1.-fsir);
        bb2 =  1./(1.-fsir);
      }
      else if (icg == 1) {
        // Oxydant
        aa1 =  1.;
        bb1 = -1./fsir;
        aa2 =  0.;
        bb2 =  0.;
      }
      else if (icg == 2) {
        // Products
        aa1 =  0.;
        bb1 =  1./fsir;
        aa2 =  1./(1.-fsir);
        bb2 = -1./(1.-fsir);
      }

      if (indpdf[c_id] == 1) {
        // Integrate PDF

        cpro_csca[icg][c_id] =   dirmin[c_id] * (aa1 + bb1*fmini)
                               + dirmax[c_id] * (aa2 + bb2*fmaxi);
        if (fdeb[c_id] < fsir) {
          cs_real_t f1 = fdeb[c_id];
          cs_real_t f2 = std::min(fsir, ffin[c_id]);
          cpro_csca[icg][c_id] += hrec[c_id]*(f2-f1)*(aa1+bb1*5.e-1*(f2+f1));
        }
        if (ffin[c_id] > fsir) {
          cs_real_t f1 = std::max(fsir, fdeb[c_id]);
          cs_real_t f2 = ffin[c_id];
          cpro_csca[icg][c_id] += hrec[c_id]*(f2-f1)*(aa2+bb2*5.e-1*(f2+f1));
        }
      }
      else {
        // Degenerate to mean value

        if (fm <= fsir) {
          cpro_csca[icg][c_id] = aa1+bb1*fm;
        }
        else {
          cpro_csca[icg][c_id] = aa2+bb2*fm;
        }
      }
    }
  });  // Loop on cells

  /* Local determination of enthalpy of stoichiometric burned gases
     in permeatic (non adiabatic), stored in work array w1.
     -------------------------------------------------------------- */

  // Initialization

  cs_array_real_set_scalar(n_cells, cm->hstoea, w1);

  if (sub_type == 1) {
    _d3phst(n_cells, indpdf, dirmin, dirmax, fdeb, ffin, hrec,
            cvar_fm, cvar_scalt, w1);
  }

  /* Integration a) of temperature
   *             b) of the radiative absorption coefficient
   *             c) of terms T^4 and T^3 in case of radiation
   *             d) of density
   *--------------------------------------------------------- */

  // Positions of variables, coefficients

  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_rho = CS_F_(rho)->val;

  cs_real_t *cpro_ckabs = nullptr, *cpro_t4m = nullptr, *cpro_t3m = nullptr;
  if (cs_glob_rad_transfer_params->type != CS_RAD_TRANSFER_NONE) {
    cpro_ckabs = cm->ckabs->val;
    cpro_t4m = cm->t4m->val;
    cpro_t3m = cm->t3m->val;
  }

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    cs_real_t fm   = cvar_fm[c_id];

    cs_real_t temsmm = 0.;

    if (indpdf[c_id] == 1) {
      // Integrate PDF

      int i_h = 0;
      for (int jh = 0; jh < (nmaxh-1); jh++) {
        if (w1[c_id] > hh[jh+1] && w1[c_id] <= hh[jh]) {
          i_h = jh;
          break;
        }
      }
      if (w1[c_id] >= hh[0])
        i_h = 0;
      if (w1[c_id] <= hh[nmaxh-1])
        i_h = nmaxh-2;

      cpro_temp[c_id] = dirmin[c_id]*tinoxy + dirmax[c_id]*tinfue;
      temsmm = dirmin[c_id]/wmolg[1]*tinoxy + dirmax[c_id]/wmolg[0]*tinfue;

      if (cpro_ckabs != nullptr) {
        cpro_ckabs[c_id] = dirmin[c_id]*ckabsg[1]  + dirmax[c_id]*ckabsg[0];
        cpro_t4m[c_id] =   dirmin[c_id]*tinoxy4 + dirmax[c_id]*tinfue4;
        cpro_t3m[c_id] =   dirmin[c_id]*tinoxy3 + dirmax[c_id]*tinfue3;
      }

      int i_f = 0;
      for (int jf = 0; jf < (nmaxf-1); jf++) {
        if (fdeb[c_id] >= ff[jf] && fdeb[c_id] < ff[jf+1]) {
          i_f = jf;
          break;
        }
      }
      if (fdeb[c_id] <= ff[0])
        i_f = 0;
      if (fdeb[c_id] >= ff[nmaxf-1])
        i_f = nmaxf-2;

      cs_real_t f2 = 0.;
      cs_real_t f1 = fdeb[c_id];

      while ((ffin[c_id]-f2) > epzero) {
        f2 = std::min(ff[i_f+1], ffin[c_id]);
        // In array TFH,
        // we extract on each line i : T = Ai+Bi*F
        // and we build for the current value of hstoe (w1)
        // T = A+B*F
        cs_real_t aa1 = tfh[i_h][i_f];
        cs_real_t bb1 = (tfh[i_h][i_f+1]-tfh[i_h][i_f]) / (ff[i_f+1]-ff[i_f]);
        cs_real_t aa2 = tfh[i_h+1][i_f];
        cs_real_t bb2 =   (tfh[i_h+1][i_f+1]-tfh[i_h+1][i_f])
                        / (ff[i_f+1]-ff[i_f]);
        cs_real_t a = aa1 + (w1[c_id]-hh[i_h])*(aa2-aa1) / (hh[i_h+1]-hh[i_h]);
        cs_real_t b = bb1 + (w1[c_id]-hh[i_h])*(bb2-bb1) / (hh[i_h+1]-hh[i_h]);
        a = a - b*ff[i_f];

        // Compute temperature by integration

        cpro_temp[c_id] += hrec[c_id]*(f2-f1)*(a+b*(f1+f2)/2.);

        // Prepare computation of absorption coefficient, of T^4 and T^3
        //   Poor side
        //     unsmm = (fs-f)/fs / wmolg[1]+ f/fs / wmolg[2]
        //     ckabs = (fs-f)/fs * ckabsg[1] + f/fs * ckabsg[2]
        //   Rich side
        //     unsmm = (f-fs)/(1-fs)/wmolg[0] + (1-f)/(1-fs)/wmolg[2]
        //     ckabs = (f-fs)/(1-fs)*ckabsg[0]  + (1-f)/(1-fs)*ckabsg[2]
        //   Everywhere
        //     unsmm = c + df
        //     ckabs = u + vF
        //     temsmm = t*unsmm = (c+df)*(a+bf) = ca +(cb+ad)f + bd f^2
        //     T^4 = (a+bf)^4
        //         = a4 + 4a3b f + 6a2b2 f^2 + 4ab3 f^3 + b4 f^4
        //     T^3 = (a+bf)^3
        //         = a3 + 3a2b f + 3ab2 f^2 + b3 f^3

        cs_real_t c, d;
        if (f1 < fsir) {
          // We started on poor side
          c =   1./wmolg[1];
          d = (-1./wmolg[1] + 1./wmolg[2]) / fsir;
        }
        else {
          // We finish on rich side (starting with f1=fs)
          c = (  -fsir/wmolg[0] + 1./wmolg[2]) / (1.-fsir);
          d = (   1./wmolg[0]   - 1./wmolg[2]) / (1.-fsir);
        }

        cs_real_t f1_2 = cs_math_pow2(f1), f2_2 = cs_math_pow2(f2);
        cs_real_t f1_3 = cs_math_pow3(f1), f2_3 = cs_math_pow3(f2);

        if (cpro_ckabs != nullptr) { // If radiative model is active
          cs_real_t u, v;

          if (f1 < fsir) {
            // We started on poor side
            u =   ckabsg[1];
            v = (-ckabsg[1] + ckabsg[2]) / fsir;
          }
          else {
            // We finish on rich side (starting with f1=fs)
            u = (-fsir*ckabsg[0] + ckabsg[2]) / (1.-fsir);
            v = (      ckabsg[0] - ckabsg[2]) / (1.-fsir);
          }

          // Compute absorption coefficien and terms T^4 and T^3
          cpro_ckabs[c_id] += hrec[c_id]*(u*(f2-f1) + v*(f2*f2-f1*f1)*0.5);

          cs_real_t a_2 = cs_math_pow2(a), b_2 = cs_math_pow2(b);
          cs_real_t a_3 = cs_math_pow3(a), b_3 = cs_math_pow3(b);
          cs_real_t a_4 = cs_math_pow4(a), b_4 = cs_math_pow4(b);
          cs_real_t f1_4 = cs_math_pow4(f1), f2_4 = cs_math_pow4(f2);
          cs_real_t f1_5 = cs_math_pow5(f1), f2_5 = cs_math_pow5(f2);

          cpro_t4m[c_id]
            += hrec[c_id] * (  (   a_4         ) * (f2-f1)
                             + (4.*a_3  *b     ) * (f2_2-f1_2)/2.
                             + (6.*(a_2)*(b_2) ) * (f2_3-f1_3)/3.
                             + (4.*a    *(b_3) ) * (f2_4-f1_4)/4.
                             + (         (b_4) ) * (f2_5-f1_5)/5.);

          cpro_t3m[c_id]
            += hrec[c_id] * (  (    a_3       ) * (f2-f1)
                             + (3.*(a_2)*b    ) * (f2_2-f1_2)/2.
                             + (3.*a    *(b_2)) * (f2_3-f1_3)/3.
                             + (         (b_3)) * (f2_4-f1_4)/4.);
        }

        // Compute temperature / molar mass term

        temsmm += hrec[c_id] * (   a*c      * (f2-f1)
                                + (c*b+a*d) * (f2_2-f1_2)/2.
                                +  b*d      * (f2_3-f1_3)/3.);

        i_f++;
        f1 = f2;
      }

    }
    else {
      // Degenerate to mean value

      int i_h = 0;
      for (int jh = 0; jh < (nmaxh-1); jh++) {
        if (w1[c_id] > hh[jh+1] && w1[c_id] <= hh[jh]) {
          i_h = jh;
          break;
        }
      }
      if (w1[c_id] >= hh[0])
        i_h = 0;
      if (w1[c_id] <= hh[nmaxh-1])
        i_h = nmaxh-2;

      int i_f = 0;
      for (int jf = 0; jf < (nmaxf-1); jf++) {
        if (fdeb[c_id] >= ff[jf] && fdeb[c_id] < ff[jf+1]) {
          i_f = jf;
          break;
        }
      }
      if (fm <= ff[0])
        i_f = 0;
      if (fm >= ff[nmaxf-1])
        i_f = nmaxf-2;

      cs_real_t aa1 = tfh[i_h][i_f];
      cs_real_t bb1 = (tfh[i_h][i_f+1]-tfh[i_h][i_f])/(ff[i_f+1]-ff[i_f]);
      cs_real_t aa2 = tfh[i_h+1][i_f];
      cs_real_t bb2 = (tfh[i_h+1][i_f+1]-tfh[i_h+1][i_f])/(ff[i_f+1]-ff[i_f]);
      cs_real_t a  = aa1 + (w1[c_id]-hh[i_h])*(aa2-aa1)/(hh[i_h+1]-hh[i_h]);
      cs_real_t b  = bb1 + (w1[c_id]-hh[i_h])*(bb2-bb1)/(hh[i_h+1]-hh[i_h]);
      a  = a - b*ff[i_f];

      // Compute temperature from mean value

      cpro_temp[c_id] = a + b*fm;

      cs_real_t c, d;
      if (fm < fsir) {
        // We started on poor side
        c =   1./wmolg[1];
        d = (-1./wmolg[1] + 1./wmolg[2]) / fsir;
      }
      else {
        // We finish on rich side (starting with f1=fs)
        c = (  -fsir/wmolg[0] + 1./wmolg[2]) / (1.-fsir);
        d = (     1./wmolg[0] - 1./wmolg[2]) / (1.-fsir);
      }

      cs_real_t fm_2 = cs_math_pow2(fm);

      if (cpro_ckabs != nullptr) { // If radiative model is active
        cs_real_t u, v;

        if (fm < fsir) {
          // We started on poor side
          u =   ckabsg[1];
          v = (-ckabsg[1]+ ckabsg[2]) / fsir;
        }
        else {
          // We finish on rich side (starting with f1=fs)
          u = (-fsir*ckabsg[0]+ ckabsg[2]) / (1.-fsir);
          v = (      ckabsg[0]- ckabsg[2]) / (1.-fsir);
        }

        // Compute absorption coefficient and terms T^4 and T^3
        // based on the mean value, in case of radiation.

        cs_real_t fm_3 = cs_math_pow3(fm);
        cs_real_t fm_4 = cs_math_pow4(fm);
        cs_real_t a_2 = cs_math_pow2(a), b_2 = cs_math_pow2(b);
        cs_real_t a_3 = cs_math_pow3(a), b_3 = cs_math_pow3(b);
        cs_real_t a_4 = cs_math_pow4(a), b_4 = cs_math_pow4(b);

        cpro_ckabs[c_id] = u + v*fm;
        cpro_t4m[c_id] =   a_4
                         + (4.*(a_3)*b    ) * fm
                         + (6.*(a_2)*(b_2)) * fm_2
                         + (4.*a    *(b_3)) * fm_3
                         + (         (b_4)) * fm_4;

        cpro_t3m[c_id] =   a_3
                         + ( 3.*(a_2)*b    ) * fm
                         + ( 3.*a    *(b_2)) * fm_2
                         + (          (b_3)) * fm_3;
      }

      // Compute temperature / molar mass term

      temsmm = a*c +(c*b+a*d)*fm + b*d*fm_2;
    }

    // Density

    if (update_rho) {
      cpro_rho[c_id] =   srrom*cpro_rho[c_id]
                       + (1.-srrom)*(pther/(cs_physical_constants_r*temsmm));
    }
  });
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
