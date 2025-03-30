/*============================================================================
 * Coal combustion model: physical properties
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
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary_conditions.h"
#include "comb/cs_coal.h"
#include "comb/cs_coal_boundary_conditions.h"
#include "comb/cs_coal_ht_convert.h"
#include "comb/cs_coal_source_terms.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_halo.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "comb/cs_coal_physical_properties.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal_noxst.cpp

  \brief Coal combustion model: NOx source terms computation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute mean gaseous concentrations.
 *
 * \param[in]     indpdf  use pdf
 * \param[in]     f1m     mean of tracer 1 mvl [chx1m+co]
 * \param[in]     f2m     mean of tracer 2 mvl [chx2m+co]
 * \param[in]     f3m     mean of tracer 3 (oxydant 1)
 * \param[in]     f4m     mean of tracer 4 (oxydant 2)
 * \param[in]     f5m     mean of tracer 5 (oxydant 3)
 * \param[in]     f6m     mean of tracer 6 (humidity)
 * \param[in]     f7m     mean of tracer 7 (C + O2)
 * \param[in]     f8m     mean of tracer 8 (C + CO2)
 * \param[in]     f9m     mean of tracer 9 (C + H2O)
 * \param[in]     pdfm1   lower bound of pdf
 * \param[in]     pdfm2   upper bound of pdf
 * \param[in]     dfuel   amplitude of Dirac at 0
 * \param[in]     doxud   amplitude of Dirac at 0
 * \param[in]     hrec    height of pdf rectangle
 * \param[in]     fuel1   mass fraction of chx1m
 * \param[in]     fuel2   mass fraction of chx2m
 * \param[in]     fuel3   mass fraction of co
 * \param[in]     fuel4   mass fraction of h2s
 * \param[in]     fuel5   mass fraction of h2
 * \param[in]     fuel6   mass fraction of hcn
 * \param[in]     oxyd    mass fraction of o2
 * \param[in]     prod1   mass fraction of co2
 * \param[in]     prod2   mass fraction of h2o
 * \param[in]     prod3   mass fraction of so2
 * \param[in]     prod4   mass fraction of nh3
 * \param[in]     xiner   mass fraction of n2
 */
/*----------------------------------------------------------------------------*/

static void
_gas_comb(cs_lnum_t        n_cells,
          int              icb1,
          int              icb2,
          int              indpdf[],
          const cs_real_t  f1m[],
          const cs_real_t  f2m[],
          const cs_real_t  f3m[],
          const cs_real_t  f4m[],
          const cs_real_t  f5m[],
          const cs_real_t  f6m[],
          const cs_real_t  f7m[],
          const cs_real_t  f8m[],
          const cs_real_t  f9m[],
          const cs_real_t  pdfm1[],
          const cs_real_t  pdfm2[],
          const cs_real_t  doxyd[],
          const cs_real_t  dfuel[],
          const cs_real_t  hrec[],
          const cs_real_t  af1[],
          const cs_real_t  af2[],
          const cs_real_t  cx1m[],
          const cs_real_t  cx2m[],
          const cs_real_t  wmf1[],
          const cs_real_t  wmf2[],
          cs_real_t        fuel1[],
          cs_real_t        fuel2[],
          cs_real_t        fuel3[],
          cs_real_t        fuel4[],
          cs_real_t        fuel5[],
          cs_real_t        fuel6[],
          cs_real_t        fuel7[],
          cs_real_t        oxyd[],
          cs_real_t        prod1[],
          cs_real_t        prod2[],
          cs_real_t        prod3[],
          cs_real_t        xiner[],
          cs_real_t        fs3no[],
          cs_real_t        fs4no[],
          cs_real_t        yfs4no[])
{
  const cs_coal_model_t *cm = cs_glob_coal_model;

  /* Initialization
   * -------------- */

  const int ico = cm->ico -1;
  const int ih2s = cm->ih2s -1;
  const int ihy =cm->ihy -1;
  const int ihcn = cm->ihcn -1;
  const int inh3 = cm->inh3 -1;
  const int io2 = cm->io2 -1;
  const int ico2 = cm->ico2 -1;
  const int ih2o = cm->ih2o -1;
  const int iso2 = cm->iso2 -1;
  const int in2 = cm->in2 -1;

  /* Aliases for simpler syntax */

  const cs_real_t *wmole = cm->wmole;
  const int n_gas_sp = cm->n_gas_species;

  const int nox_model = cm->ieqnox;

  /* Preliminary computations
     ------------------------ */

  cs_real_t *cvar_yco2 = nullptr, *x1 = nullptr;

  const auto af3 = cm->af3, af4 = cm->af4;
  const auto af5 = cm->af5, af6 = cm->af6;
  const auto af7 = cm->af7, af8 = cm->af8;
  const auto af9 = cm->af9;

  if (cm->ieqco2 == 1)
    cvar_yco2 = cs_field_by_id(cm->iyco2)->val;

  // Massic fraction of gas
  x1 = cs_field_by_name("x_c")->val;

  // Molar masses
  double wmh2s = cm->wmole[ih2s];
  double wmh2  = cm->wmole[ihy];
  double wmhcn = wmole[ihcn];
  double wmnh3 = wmole[inh3];
  double wmco  = cm->wmole[ico];
  double wmo2  = cm->wmole[io2];
  double wmco2 = cm->wmole[ico2];
  double wmh2o = cm->wmole[ih2o];
  double wmso2 = cm->wmole[iso2];
  double wmn2  = cm->wmole[in2];

  // Initialization of mass fractions with oxydant 1

  #pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {
    fuel1[c_id]  = wmf1[c_id]*af3[icb1];
    fuel2[c_id]  = wmf2[c_id]*af3[icb2];
    fuel3[c_id]  = wmco     *af3[ico];
    fuel4[c_id]  = wmh2s    *af3[ih2s];
    fuel5[c_id]  = wmh2     *af3[ihy];
    fuel6[c_id]  = wmhcn    *af3[ihcn];
    fuel7[c_id]  = wmnh3    *af3[inh3];
    oxyd[c_id]   = wmo2     *af3[io2];
    prod1[c_id]  = wmco2    *af3[ico2];
    prod2[c_id]  = wmh2o    *af3[ih2o];
    prod3[c_id]  = wmso2    *af3[iso2];
    xiner[c_id]  = wmn2     *af3[in2];
  }

  /* Compute mixture composition without the PDF
     if the fluctuations are too small
     ------------------------------------------- */

  #pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {

    if (indpdf[c_id] != 0)
      continue;

    cs_real_t zz[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

    // preliminary calculations
    cs_real_t zco2t = 0.0;
    if (cvar_yco2 != nullptr)
      zco2t = (cvar_yco2[c_id] / x1[c_id]) / wmco2;

    // Composition of gaseous phase before combustion

    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zz[g_id] =   af1[n_cells*g_id + c_id] * f1m[c_id]
                 + af2[n_cells*g_id + c_id] * f2m[c_id]
                 + af3[g_id]*f3m[c_id] + af4[g_id]*f4m[c_id]
                 + af5[g_id]*f5m[c_id] + af6[g_id]*f6m[c_id]
                 + af7[g_id]*f7m[c_id] + af8[g_id]*f8m[c_id]
                 + af9[g_id]*f9m[c_id];
    }

    // Compute mixture composition

    // 1st reaction:

    cs_real_t anu1   = 0.5;
    cs_real_t reac1  = std::min(zz[ihy],(zz[io2]/anu1));
    zz[ihy] = zz[ihy]  -      reac1;
    zz[io2] = zz[io2]  - anu1*reac1;
    zz[ih2o]= zz[ih2o] +      reac1;

    // 2nd reaction:

    cs_real_t anu2   = 0.25*abs(cx1m[c_id]-cx2m[c_id]);
    cs_real_t reac2  = std::min(zz[icb1], (zz[io2]/anu2));
    zz[icb1]= zz[icb1] -           reac2;
    zz[icb2]= zz[icb2] +           reac2;
    zz[io2]  = zz[io2]   -      anu2*reac2;
    zz[ih2o] = zz[ih2o]  + 2.0 *anu2*reac2;

    // 3rd reaction:

    cs_real_t anu3  = 0.25*(2.0 + cx2m[c_id]);
    cs_real_t reac3 = std::min(zz[icb2], (zz[io2]/anu3));
    zz[icb2]  = zz[icb2] -                reac3;
    zz[ico]   = zz[ico]  +                reac3;
    zz[io2]   = zz[io2]  -           anu3*reac3;
    zz[ih2o]  = zz[ih2o] + 0.5*cx2m[c_id]*reac3;

    // 4th reaction:

    cs_real_t anu4  = 1.5;
    cs_real_t reac4 = std::min(zz[ih2s], (zz[io2]/anu4));
    zz[ih2s] = zz[ih2s]   -      reac4;
    zz[io2]  = zz[io2]    - anu4*reac4;
    zz[ih2o] = zz[ih2o]   +      reac4;
    zz[iso2] = zz[iso2]   +      reac4;

    // 5th reaction:

    cs_real_t anu5   = 0.5;
    cs_real_t reac5  = std::min(zz[ico], (zz[io2]/anu5));
    if (cvar_yco2 != nullptr)
      reac5  = std::min(std::max(zco2t-zz[ico2], 0.0), reac5);

    zz[ico] = zz[ico]    -      reac5;
    zz[io2] = zz[io2]    - anu5*reac5;
    zz[ico2]= zz[ico2]   +      reac5;

    fuel1[c_id] = zz[icb1] * wmf1[c_id];
    fuel2[c_id] = zz[icb2] * wmf2[c_id];
    fuel3[c_id] = zz[ico]  * wmco;
    fuel4[c_id] = zz[ih2s] * wmh2s;
    fuel5[c_id] = zz[ihy]  * wmh2;
    fuel6[c_id] = zz[ihcn] * wmhcn;
    fuel7[c_id] = zz[inh3] * wmnh3;
    oxyd [c_id] = zz[io2]  * wmo2;
    prod1[c_id] = zz[ico2] * wmco2;
    prod2[c_id] = zz[ih2o] * wmh2o;
    prod3[c_id] = zz[iso2] * wmso2;
    xiner[c_id] = zz[in2]  * wmn2;

  }

  /* Compute mixture composition with the PDF
     ---------------------------------------- */

  { // Extra brace to work around nvcc 11 compiler bug

  #pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {

    if (indpdf[c_id] == 0)
      continue;

    cs_real_t zz[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t zzox[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t zzcl[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t zzs1[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t zzs2[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t zzs3[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t zzs4[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

    // preliminary calculations
    cs_real_t zco2t = 0.0;
    if (cvar_yco2 != nullptr)
      zco2t = (cvar_yco2[c_id] / x1[c_id]) / wmco2;

    // local oxydizer
    cs_real_t soxy =   f3m[c_id] + f4m[c_id] + f5m[c_id] + f6m[c_id]
                     + f7m[c_id] + f8m[c_id] + f9m[c_id];

    zzox[icb1] = 0.0;
    zzox[icb2] = 0.0;
    for (cs_lnum_t g_id = 2; g_id < n_gas_sp; g_id++) {
      zzox[g_id] = (  af3[g_id]*f3m[c_id] + af4[g_id]*f4m[c_id]
                    + af5[g_id]*f5m[c_id] + af6[g_id]*f6m[c_id]
                    + af7[g_id]*f7m[c_id] + af8[g_id]*f8m[c_id]
                    + af9[g_id]*f9m[c_id]) / soxy;
    }

    // Local fuel
    cs_real_t scomb = f1m[c_id] + f2m[c_id];
    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zzcl[g_id] = (  af1[n_cells*g_id + c_id] * f1m[c_id]
                    + af2[n_cells*g_id + c_id] * f2m[c_id]) / scomb;
    }

    // 1st reaction: hydrogen recombination

    cs_real_t reac1 = std::min(zzox[ihy], (2.0*zzox[io2]));
    zzox[ihy]  = zzox[ihy]  -     reac1;
    zzox[ih2o] = zzox[ih2o] +     reac1;
    zzox[io2]  = zzox[io2]  - 0.5*reac1;

    // 2nd reaction: CHx1 + (x1-x2)/4 O2  => CHx2 + (x1-x2)/2 H2O

    cs_real_t fs1 = 1.0;
    if (zzcl[icb1] > 0.0  &&  zzox[io2] > 0.0)
      fs1 =   zzox[io2]/(  std::abs(cx1m[c_id]-cx2m[c_id])*0.25 * zzcl[icb1]
                         + zzox[io2]);

    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zzs1[g_id] = fs1*zzcl[g_id] + (1.0-fs1)*zzox[g_id];
    }
    zzs1[icb2] = zzs1[icb1] + zzs1[icb2];
    zzs1[ih2o] = zzs1[ih2o] + 0.5 * std::abs(cx1m[c_id]-cx2m[c_id]) * zzs1[icb1];
    zzs1[icb1] = 0.0;
    zzs1[io2]  = 0.0;

    // 3rd reaction: CHx2 + (2+x2)/4 O2 => CO + x2/2 H2O

    cs_real_t fs2 = fs1;
    if (zzs1[icb2] > 0.0 && zzox[io2] > 0.0)
      fs2 = fs1 * zzox[io2] / (  (2.0+cx2m[c_id])*0.25 * zzs1[icb2]
                               + zzox[io2]);

    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zzs2[g_id] = (fs2/fs1)*zzs1[g_id]+(1.0-(fs2/fs1))*zzox[g_id];
    }
    zzs2[ico] = zzs2[ico] + zzs2[icb2];
    zzs2[ih2o]= zzs2[ih2o] + cx2m[c_id] * 0.5 * zzs2[icb2];
    zzs2[icb2]= 0.0;
    zzs2[io2] = 0.0;

    // 4eme reaction:  H2S + 3/2 O2 => H2O + SO2

    cs_real_t fs3 = fs2;
    if (zzs2[ih2s] > 0.0 && zzox[io2] > 0.0)
      fs3 = fs2 * zzox[io2] / (1.5*zzs2[ih2s]+zzox[io2]);

    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zzs3[g_id] = (fs3/fs2)*zzs2[g_id]+(1.0-(fs3/fs2))*zzox[g_id];
    }
    zzs3[iso2] = zzs3[iso2] + zzs3[ih2s];
    zzs3[ih2o] = zzs3[ih2o] + zzs3[ih2s];
    zzs3[ih2s] = 0.0;
    zzs3[io2]  = 0.0;

    // 5th reaction CO+1/2 O2 => CO2

    cs_real_t fs4 = fs3;
    if (   (zzs3[ico] > 0.0 && zzox[io2] > 0.0)
        && (cvar_yco2 == nullptr)) {
      fs4 = fs3 * zzox[io2] / (zzs3[ico] + zzox[io2]);
    }

    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zzs4[g_id] = (fs4/fs3)*zzs3[g_id] + (1.0-(fs4/fs3))*zzox[g_id];
    }

    if (cvar_yco2 == nullptr) {
      zzs4[ico2] = zzs4[ico2] + zzs4[ico];
      zzs4[io2]  = zzs4[io2] * 0.5;
      zzs4[ico]  = 0.0;
    }

    // Store fs3, fs4 and concentrations in fs4 for the NOx model

    if (nox_model == 1) {

      fs3no[c_id] = fs3;

      if (zzs3[ico] > 0.0 && zzox[io2] > 0.0)
        fs4no[c_id] = fs3 * zzox[io2] / (zzs3[ico] + zzox[io2]);
      else
        fs4no[c_id] = fs3;

      for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
        yfs4no[n_cells*g_id + c_id]
          =   (fs4no[c_id]/fs3) * zzs3[g_id]
            + (1.0-(fs4no[c_id]/fs3)) * zzox[g_id];
      }
      yfs4no[n_cells*ico2 + c_id] += yfs4no[n_cells*ico + c_id];
      yfs4no[n_cells*io2  + c_id] *= 0.5;
      yfs4no[n_cells*ico + c_id] = 0.0;

      yfs4no[n_cells*icb1 + c_id] *= wmf1[c_id];
      yfs4no[n_cells*icb2 + c_id] *= wmf2[c_id];
      for (cs_lnum_t g_id = ico; g_id < n_gas_sp; g_id++) {
        yfs4no[n_cells*g_id + c_id] *= wmole[g_id];
      }

    }

    // We now know the concentrations
    // cl,s1,s2,s3,s4,Ox
    // The intermediate concentrations are piecewise linear
    // and the parameters of the pdf dfuel, doxyd, pdfm1, pdfm2, hrec

    // Initialize by concentrations at extremities
    for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
      zz[g_id] = dfuel[c_id]*zzcl[g_id] + doxyd[c_id]*zzox[g_id];
    }

    // Integration on first richness interval (between s1 and 1)
    cs_real_t bb1 = std::max(pdfm1[c_id], fs1);
    cs_real_t bb2 = std::min(pdfm2[c_id], 1.0);
    if (bb2 > bb1) {
      for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
        zz[g_id] =  zz[g_id] + hrec[c_id]*(bb2-bb1)/(1.0-fs1)
                  * (   zzs1[g_id]-fs1*zzcl[g_id]
                     + (zzcl[g_id]-zzs1[g_id])*(bb1+bb2)*0.5);
      }
    }

    // Integration on second richness interval (between s2 and s1)
    bb1 = std::max(pdfm1[c_id], fs2);
    bb2 = std::min(pdfm2[c_id], fs1);
    if (bb2 > bb1) {
      for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
        zz[g_id] =   zz[g_id] + hrec[c_id]*(bb2-bb1)/(fs1-fs2)
                   * (   fs1*zzs2[g_id]-fs2*zzs1[g_id]
                      + (zzs1[g_id]-zzs2[g_id])*(bb1+bb2)*0.5);
      }
    }

    // Integration on third richness interval (between s3 and s2)
    bb1 = std::max(pdfm1[c_id], fs3);
    bb2 = std::min(pdfm2[c_id], fs2);
    if (bb2 > bb1) {
      for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
        zz[g_id] =   zz[g_id] + hrec[c_id]*(bb2-bb1)/(fs2-fs3)
                   * (  fs2*zzs3[g_id]-fs3*zzs2[g_id]
                      + (zzs2[g_id]-zzs3[g_id])*(bb1+bb2)*0.5);
      }
    }

    // Integration on fourth richness interval (between s4 and s3)
    bb1 = std::max(pdfm1[c_id], fs4);
    bb2 = std::min(pdfm2[c_id], fs3);
    if (bb2 > bb1) {
      for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
        zz[g_id] =   zz[g_id] + hrec[c_id]*(bb2-bb1)/(fs3-fs4)
                   * (  fs3*zzs4[g_id]-fs4*zzs3[g_id]
                      + (zzs3[g_id]-zzs4[g_id])*(bb1+bb2)*0.5);
      }
    }

    // Integration on fifth richness interval (between 0and s4)
    bb1 = std::max(pdfm1[c_id], 0.0);
    bb2 = std::min(pdfm2[c_id], fs4);
    if (bb2 > bb1) {
      for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {
        zz[g_id] =   zz[g_id] + hrec[c_id]*(bb2-bb1)/(fs4-0.0)
                   * (  fs4*zzox[g_id]-0.0*zzs4[g_id]
                      + (zzs4[g_id]-zzox[g_id])*(bb1+bb2)*0.5);
      }
    }

    // CO2 transport
    if (cvar_yco2 != nullptr) {
      cs_real_t anu5   = 0.5;
      cs_real_t reac5  = std::min(zz[ico], zz[io2]/anu5);
      reac5  = std::min(std::max(zco2t-zz[ico2], 0.0),  reac5);

      zz[ico]  = zz[ico]  -      reac5;
      zz[io2]  = zz[io2]  - anu5*reac5;
      zz[ico2] = zz[ico2] +      reac5;
    }

    fuel1[c_id] = zz[icb1] * wmf1[c_id];
    fuel2[c_id] = zz[icb2] * wmf2[c_id];
    fuel3[c_id] = zz[ico]  * wmco;
    fuel4[c_id] = zz[ih2s] * wmh2s;
    fuel5[c_id] = zz[ihy]  * wmh2;
    fuel6[c_id] = zz[ihcn] * wmhcn;
    fuel7[c_id] = zz[inh3] * wmnh3;
    oxyd[c_id]  = zz[io2]  * wmo2;
    prod1[c_id] = zz[ico2] * wmco2;
    prod2[c_id] = zz[ih2o] * wmh2o;
    prod3[c_id] = zz[iso2] * wmso2;
    xiner[c_id] = zz[in2]  * wmn2;

  }}

  /* Logging
     ------- */

  if (cs_log_default_is_active() == false)
    return;

  cs_gnum_t n1 = (cs_gnum_t)n_cells;
  cs_gnum_t n2 = 0, n3 = 0, n4 = 0, n5 = 0, n6 = 0, n7 = 0, n8 = 0;
  cs_gnum_t n9 = 0, n10 = 0, n11 = 0, n12 = 0, n13 = 0, n14 = 0, n15 = 0;

  // Control pdf parameters and different values of mass fractions

  cs_real_t sommin =  HUGE_VALF;
  cs_real_t sommax = -HUGE_VALF;

  const cs_real_t epzero = cs_coal_epsilon;

  for (auto c_id = 0; c_id < n_cells; c_id++) {

    if (indpdf[c_id] != 0)
      n2 += 1;

    cs_real_t somm =   fuel1[c_id] + fuel2[c_id] + fuel3[c_id]
                     + fuel4[c_id] + fuel5[c_id] + fuel6[c_id] + fuel7[c_id]
                     + oxyd[c_id]
                     + prod1[c_id] + prod2[c_id] + prod3[c_id]
                     + xiner[c_id];

    sommin = std::min(sommin, somm);
    sommax = std::max(sommax, somm);

    if (std::abs(somm-1.0) < epzero)
      n3 += 1;

    if (fuel1[c_id] < -epzero || fuel1[c_id] > (1.+epzero)) n4  += 1;
    if (fuel2[c_id] < -epzero || fuel2[c_id] > (1.+epzero)) n5  += 1;
    if (fuel3[c_id] < -epzero || fuel3[c_id] > (1.+epzero)) n6  += 1;
    if (fuel4[c_id] < -epzero || fuel4[c_id] > (1.+epzero)) n7  += 1;
    if (fuel5[c_id] < -epzero || fuel5[c_id] > (1.+epzero)) n8  += 1;
    if (fuel6[c_id] < -epzero || fuel6[c_id] > (1.+epzero)) n9  += 1;
    if (fuel7[c_id] < -epzero || fuel7[c_id] > (1.+epzero)) n10 += 1;
    if (oxyd[c_id]  < -epzero || oxyd[c_id]  > (1.+epzero)) n11 += 1;
    if (xiner[c_id] < -epzero || xiner[c_id] > (1.+epzero)) n12 += 1;
    if (prod1[c_id] < -epzero || prod1[c_id] > (1.+epzero)) n13 += 1;
    if (prod2[c_id] < -epzero || prod2[c_id] > (1.+epzero)) n14 += 1;
    if (prod3[c_id] < -epzero || prod3[c_id] > (1.+epzero)) n15 += 1;

  }

  cs_parall_sum_scalars(n1, n2, n3, n4, n5, n6, n7, n8, n9,
                        n10, n11, n12, n13, n14, n15);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       "Combustion modeling with turbulent diffusion model (CPCYM2)\n"
       "Fast 3-point chemistry - extension to 3 fuels\n"
       "===========================================================\n"
       "  Nb. of computation points                       : %llu\n"
       "  Nb. of turbulent points (using PDFs)            : %llu\n\n"),
     (unsigned long long)n1, (unsigned long long)n2);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("Control of mass fraction values\n"
       "  Nb. computation points verifying sum of Yi = 1  : %llu\n"
       "  Nb. points YCHX1, YCHX2  < 0 or > 1             : %llu, %llu\n"
       "  Nb. points YC0, YH2S     < 0 or > 1             : %llu, %llu\n"
       "  Nb. points YH2, YHCN     < 0 or > 1             : %llu, %llu\n"
       "  Nb. points YNH3          < 0 or > 1             : %llu\n"
       "  Nb. points YO2, YN2      < 0 or > 1             : %llu, %llu\n"
       "  Nb. points YCO2, YH2O    < 0 or > 1             : %llu, %llu\n"
       "  Nb. points YSO2          < 0 or > 1             : %llu\n"),
     (unsigned long long)n3, (unsigned long long)n4,
     (unsigned long long)n5, (unsigned long long)n6,
     (unsigned long long)n7, (unsigned long long)n8,
     (unsigned long long)n9, (unsigned long long)n10,
     (unsigned long long)n11, (unsigned long long)n12,
     (unsigned long long)n13, (unsigned long long)n14,
     (unsigned long long)n15);

  cs_parall_min(1, CS_REAL_TYPE, &sommin);
  cs_parall_max(1, CS_REAL_TYPE, &sommax);

  cs_log_printf(CS_LOG_DEFAULT,
                _(" Sum Min Max : %g %g\n"),
                sommin, sommax);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute physical properties in gaseous phase.
 *
 * Temperature, mass density and average concentrations
 * (use of a pdf rectangle-dirac)
 * ==> fast chemistry model in 3 points
 * Extension to 3 combustibles for pulverized coal
 *
 * Heterogeneous reactions
 *  - Pyrolysis
 *    Elementary composition of the mol of volatile materials
 *    The reactive coal is written C(1)H(ALPHA)O(BETA)
 *      -(k1)-> alpha/4 CH4  + beta CO + (1-alpha/4-beta)    Coke
 *    Reactive coal
 *      -(k2)-> alpha/Y CXHY + beta CO + (1-alpha/RYSX-beta) Coke
 *      With RYSX = Y/X
 *  - Heterogeneous combustion
 *    Coke + 1/2 (O2 + XSI N2) -> CO + XSI/2 N2
 *  - Gas-phase reaction
 * (4/(4-RYSX)) CH4 + (O2 + XSI N2)   -(1)->  4/X/(4-RYSX)*CXHY + 2 H2O
 *                                            + XSI N2
 * CXHY + X/4*(2+RYSX) (O2 + XSI N2)  -(2)->  X CO + Y/2 H2O
 *                                           + X/4*(2+RYSX)*XSI N2
 *            CO + 1/2 (O2 + XSI N2)  -(3)->  CO2 + XSI/2 N2
 * Variable choice
 * - F1 is the mass fractions of volatile materials: CH4  + CO
 * - F2 is the mass fractions of volatile materials: CXHY + CO
 * - F3 is the mass fraction of carbon coming from the homogeneous combustion
 * Let Y be the mass fractions and Z be the concentrations (moles/kg)
 * of index f before reaction, b final
 *
 * Joint PDF degenerated into a 1D PDF of type RECTANGLE - DIRAC
 *
 * \param[in]       f1m     mean of tracer 1 mvl [chx1m+co]
 * \param[in]       f2m     mean of tracer 2 mvl [chx2m+co]
 * \param[in]       f3m     mean of tracer 3 (oxydant 1)
 * \param[in]       f4m     mean of tracer 4 (oxydant 2)
 * \param[in]       f5m     mean of tracer 5 (oxydant 3)
 * \param[in]       f6m     mean of tracer 6 (humidity)
 * \param[in]       f7m     mean of tracer 7 (C + O2)
 * \param[in]       f8m     mean of tracer 8 (C + CO2)
 * \param[in]       f9m     mean of tracer 9 (C + H2O)
 * \param[in, out]  fvp2m   f1f2 variance
 * \param[in]       enth    enthalpy in \f$ j . kg^{-1} \f$  either of
 *                          the gas or of the mixture
 * \param[in]       enthox  oxydant enthalpy
 * \param[out]      rom1    gas density
 */
/*----------------------------------------------------------------------------*/

static void
_physprop1(const  cs_real_t  f1m[],
           const  cs_real_t  f2m[],
           const  cs_real_t  f3m[],
           const  cs_real_t  f4m[],
           const  cs_real_t  f5m[],
           const  cs_real_t  f6m[],
           const  cs_real_t  f7m[],
           const  cs_real_t  f8m[],
           const  cs_real_t  f9m[],
           cs_real_t         fvp2m[],
           const  cs_real_t  enth[],
           const  cs_real_t  enthox[],
           cs_real_t         rom1[])
{
  static int ipass = 0;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;

  const cs_coal_model_t *cm = cs_glob_coal_model;

  /* Initialization
   * -------------- */

  const int ichx1 = cm->ichx1 -1;
  const int ichx2 = cm->ichx2 -1;
  const int ico = cm->ico -1;
  const int ih2s = cm->ih2s -1;
  const int ihy =cm->ihy -1;
  const int ihcn = cm->ihcn -1;
  const int inh3 = cm->inh3 -1;
  const int io2 = cm->io2 -1;
  const int ico2 = cm->ico2 -1;
  const int ih2o = cm->ih2o -1;
  const int iso2 = cm->iso2 -1;
  const int in2 = cm->in2 -1;

  /* Aliases for simpler syntax */

  const cs_real_t *wmole = cm->wmole;
  const int n_gas_sp = cm->n_gas_species;

  // Mass fraction of gas
  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;

  /* Memory allocation
     ----------------- */

  int *intpdf;
  CS_MALLOC(intpdf, n_cells, int);

  cs_real_t *ffuel, *dfuel, *doxyd, *pdfm1;
  cs_real_t *pdfm2, *hrec, *cx1m, *cx2m, *wmchx1, *wmchx2;
  cs_real_t *af1, *af2;
  CS_MALLOC(ffuel, n_cells, cs_real_t);
  CS_MALLOC(dfuel, n_cells, cs_real_t);
  CS_MALLOC(doxyd, n_cells, cs_real_t);
  CS_MALLOC(pdfm1, n_cells, cs_real_t);
  CS_MALLOC(pdfm2, n_cells, cs_real_t);
  CS_MALLOC(hrec, n_cells, cs_real_t);
  CS_MALLOC(cx1m, n_cells, cs_real_t);
  CS_MALLOC(cx2m, n_cells, cs_real_t);
  CS_MALLOC(wmchx1, n_cells, cs_real_t);
  CS_MALLOC(wmchx2, n_cells, cs_real_t);
  CS_MALLOC(af1, n_cells*n_gas_sp, cs_real_t);
  CS_MALLOC(af2, n_cells*n_gas_sp, cs_real_t);

  cs_real_t *fs3no = nullptr, *fs4no = nullptr, *yfs4no = nullptr;

  if (cm->ieqnox == 1) {
    CS_MALLOC(fs3no, n_cells, cs_real_t);
    CS_MALLOC(fs4no, n_cells, cs_real_t);
    CS_MALLOC(yfs4no, n_cells*n_gas_sp, cs_real_t);
  }

  // Pointer to CHx1 and CHx2
  cs_real_t *cpro_cyf1 = cs_field_by_id(cm->iym1[ichx1])->val;
  cs_real_t *cpro_cyf2 = cs_field_by_id(cm->iym1[ichx2])->val;
  cs_real_t *cpro_cyf3 = cs_field_by_id(cm->iym1[ico])->val;
  cs_real_t *cpro_cyf4 = cs_field_by_id(cm->iym1[ih2s])->val;
  cs_real_t *cpro_cyf5 = cs_field_by_id(cm->iym1[ihy])->val;
  cs_real_t *cpro_cyf6 = cs_field_by_id(cm->iym1[ihcn])->val;
  cs_real_t *cpro_cyf7 = cs_field_by_id(cm->iym1[inh3])->val;
  cs_real_t *cpro_cyox = cs_field_by_id(cm->iym1[io2])->val;
  cs_real_t *cpro_cyp1 = cs_field_by_id(cm->iym1[ico2])->val;
  cs_real_t *cpro_cyp2 = cs_field_by_id(cm->iym1[ih2o])->val;
  cs_real_t *cpro_cyp3 = cs_field_by_id(cm->iym1[iso2])->val;
  cs_real_t *cpro_cyin = cs_field_by_id(cm->iym1[in2])->val;

  /* Determine the type of pdf
     ------------------------- */

  cs_real_t *fmini, *fmaxi, *tpdf;
  CS_MALLOC(fmini, n_cells, cs_real_t);
  CS_MALLOC(fmaxi, n_cells, cs_real_t);
  CS_MALLOC(tpdf, n_cells_ext, cs_real_t);

  # pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {
    // min et max boundary of the pdf
    fmini[c_id] = 0;;
    fmaxi[c_id] = 1.;
    // Sum F1+F2
    ffuel[c_id] = f1m[c_id] + f2m[c_id];
  }

  cs_combustion_dirac_pdf(n_cells, intpdf, tpdf, ffuel, fvp2m, fmini, fmaxi,
                          doxyd, dfuel, pdfm1, pdfm2, hrec);

  CS_FREE(tpdf);
  CS_FREE(fmini);
  CS_FREE(fmaxi);

  ipass++;

  /* Calculation of average concentrations
     ------------------------------------- */

  // Calculation of F8MC = sum(F8M(coal_id))
  // Calculation of F9MC = sum(F9M(coal_id))

  cs_arrays_set_value<cs_real_t, 1>(n_cells*n_gas_sp, 0., af1, af2);
  cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., cx1m, cx2m);

  const cs_real_t epzero = cs_math_epzero;

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

    /* Precompute quantities independent of interpolation point
       and set array pointers */

    int ichx1c_coal = cm->ichx1c[coal_id] -1;
    int ichx2c_coal = cm->ichx2c[coal_id] -1;

    cs_real_t den1_coal = 1. / (  cm->a1[coal_id]*cm->wmole[ichx1c_coal]
                                + cm->b1[coal_id]*cm->wmole[ico]
                                + cm->c1[coal_id]*cm->wmole[ih2o]
                                + cm->d1[coal_id]*cm->wmole[ih2s]
                                + cm->e1[coal_id]*cm->wmole[ihcn]
                                + cm->f1[coal_id]*cm->wmole[inh3]);

    cs_real_t den2_coal = 1. / (  cm->a2[coal_id]*cm->wmole[ichx2c_coal]
                                + cm->b2[coal_id]*cm->wmole[ico]
                                + cm->c2[coal_id]*cm->wmole[ih2o]
                                + cm->d2[coal_id]*cm->wmole[ih2s]
                                + cm->e2[coal_id]*cm->wmole[ihcn]
                                + cm->f2[coal_id]*cm->wmole[inh3]);

    const cs_real_t *cvar_f1m_coal = cs_field_by_id(cm->if1m[coal_id])->val;
    const cs_real_t *cvar_f2m_coal = cs_field_by_id(cm->if2m[coal_id])->val;

    const cs_real_t da1w = den1_coal * cm->a1[coal_id] * wmole[ichx1c_coal];
    const cs_real_t da1 = den1_coal * cm->a1[coal_id];
    const cs_real_t db1 = den1_coal * cm->b1[coal_id];
    const cs_real_t dc1 = den1_coal * cm->c1[coal_id];
    const cs_real_t dd1 = den1_coal * cm->d1[coal_id];
    const cs_real_t de1 = den1_coal * cm->e1[coal_id];
    const cs_real_t df1 = den1_coal * cm->f1[coal_id];

    const cs_real_t da2w = den2_coal * cm->a2[coal_id] * wmole[ichx2c_coal];
    const cs_real_t da2 = den2_coal * cm->a2[coal_id];
    const cs_real_t db2 = den2_coal * cm->b2[coal_id];
    const cs_real_t dc2 = den2_coal * cm->c2[coal_id];
    const cs_real_t dd2 = den2_coal * cm->d2[coal_id];
    const cs_real_t de2 = den2_coal * cm->e2[coal_id];
    const cs_real_t df2 = den2_coal * cm->f2[coal_id];

    # pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (auto c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t f1mc = cvar_f1m_coal[c_id] / cpro_x1[c_id];
      cx1m[c_id] += f1mc * da1w;

      af1[n_cells*ichx1 + c_id] += f1mc * da1;
      af1[n_cells*ico   + c_id] += f1mc * db1;
      af1[n_cells*ih2o  + c_id] += f1mc * dc1;
      af1[n_cells*ih2s  + c_id] += f1mc * dd1;
      af1[n_cells*ihcn  + c_id] += f1mc * de1;
      af1[n_cells*inh3  + c_id] += f1mc * df1;
    }

    # pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (auto c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t f2mc = cvar_f2m_coal[c_id] / cpro_x1[c_id];
      cx2m[c_id] += f2mc * da2w;

      af2[n_cells*ichx2 + c_id] += f2mc * da2;
      af2[n_cells*ico   + c_id] += f2mc * db2;
      af2[n_cells*ih2o  + c_id] += f2mc * dc2;
      af2[n_cells*ih2s  + c_id] += f2mc * dd2;
      af2[n_cells*ihcn  + c_id] += f2mc * de2;
      af2[n_cells*inh3  + c_id] += f2mc * df2;

    }

  }

  const cs_real_t wmolat_iatc = cm->wmolat[cs_coal_atom_id_c];
  const cs_real_t wmolat_iath = cm->wmolat[cs_coal_atom_id_h];
  const cs_real_t wmolat_iato = cm->wmolat[cs_coal_atom_id_o];

  # pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t c_af1 = af1[n_cells*ichx1 + c_id];
    cs_real_t c_af2 = af2[n_cells*ichx2 + c_id];

    if (c_af1 > epzero)
      cx1m[c_id] = (cx1m[c_id]/c_af1 - wmolat_iatc) / wmolat_iath;
    else
      cx1m[c_id] = 4.0;

    if (c_af2 > epzero)
      cx2m[c_id] = (cx2m[c_id]/c_af2 - wmolat_iatc) / wmolat_iath;
    else
      cx2m[c_id] = 2.0;

    wmchx1[c_id] = wmolat_iatc + cx1m[c_id]*wmolat_iath;
    wmchx2[c_id] = wmolat_iatc + cx2m[c_id]*wmolat_iath;

  }

  for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {

    # pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (auto c_id = 0; c_id < n_cells; c_id++) {
      if (f1m[c_id] > 0)
        af1[n_cells*g_id + c_id] /= f1m[c_id];
      else
        af1[n_cells*g_id + c_id] = 0.;

      if (f2m[c_id] > 0)
        af2[n_cells*g_id + c_id] /= f2m[c_id];
      else
        af2[n_cells*g_id + c_id] = 0.;
    }

  }

  _gas_comb(n_cells, ichx1, ichx2,
            intpdf, f1m, f2m, f3m, f4m, f5m, f6m, f7m, f8m, f9m,
            pdfm1, pdfm2, doxyd, dfuel, hrec,
            af1, af2, cx1m, cx2m, wmchx1, wmchx2,
            cpro_cyf1, cpro_cyf2, cpro_cyf3,
            cpro_cyf4, cpro_cyf5, cpro_cyf6,
            cpro_cyf7, cpro_cyox, cpro_cyp1,
            cpro_cyp2, cpro_cyp3, cpro_cyin,
            fs3no, fs4no, yfs4no);

  // Eventual clipping of mass fractions.

  for (cs_lnum_t g_id = 0; g_id < n_gas_sp; g_id++) {

    cs_real_t *cpro_ym1_sp = cs_field_by_id(cm->iym1[g_id])->val;

    # pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (auto c_id = 0; c_id < n_cells; c_id++) {

      if (cpro_ym1_sp[c_id] < cs_coal_epsilon)
        cpro_ym1_sp[c_id] = 0;
      else if (cpro_ym1_sp[c_id] > 1.)
        cpro_ym1_sp[c_id] = 1.;

    }

  }

  /* Calculation of temperature and density
     -------------------------------------- */

  cs_real_t *cpro_temp = cs_field_by_name("temperature")->val;
  cs_real_t *cpro_mmel = cs_field_by_id(cm->immel)->val;

  cs_coal_ht_convert_h_to_t_gas(CS_MESH_LOCATION_CELLS, enth, cpro_temp);

  const cs_real_t p0 = cs_glob_fluid_properties->p0;
  const cs_real_t r = cs_physical_constants_r;

  # pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t wmolme =   cpro_cyf1[c_id]/wmchx1[c_id]
                       + cpro_cyf2[c_id]/wmchx2[c_id]
                       + cpro_cyf3[c_id]/wmole[ico]
                       + cpro_cyf4[c_id]/wmole[ih2s]
                       + cpro_cyf5[c_id]/wmole[ihy]
                       + cpro_cyf6[c_id]/wmole[ihcn]
                       + cpro_cyf7[c_id]/wmole[inh3]
                       + cpro_cyox[c_id]/wmole[io2]
                       + cpro_cyp1[c_id]/wmole[ico2]
                       + cpro_cyp2[c_id]/wmole[ih2o]
                       + cpro_cyp3[c_id]/wmole[iso2]
                       + cpro_cyin[c_id]/wmole[in2];

    // storage of the mixture molar mass

    cpro_mmel[c_id] = 1.0 / wmolme;

    // We do not include the mecanical pressure IPR but P0

    rom1[c_id] = p0 / (wmolme*r*cpro_temp[c_id]);

  }

  /* Nox model
     --------- */

  if (cm->ieqnox == 1) {

    if (ipass > 1) {
      // Nox model not used at the first relative iteration

      cs_coal_noxst(intpdf,
                    pdfm1, pdfm2, doxyd, dfuel, hrec,
                    f3m, f4m, f5m, f6m, f7m, f8m, f9m,
                    fs3no, fs4no, yfs4no, enthox);
    }
    else {

      cs_real_t *cpro_ghcn1 = cs_field_by_id(cm->ighcn1)->val;
      cs_real_t *cpro_ghcn2 = cs_field_by_id(cm->ighcn2)->val;
      cs_real_t *cpro_gnoth = cs_field_by_id(cm->ignoth)->val;
      cs_real_t *cpro_gnh31 = cs_field_by_id(cm->ignh31)->val;
      cs_real_t *cpro_gnh32 = cs_field_by_id(cm->ignh32)->val;
      cs_real_t *cpro_fhcnr = cs_field_by_id(cm->ifhcnr)->val;
      cs_real_t *cpro_fhcnd = cs_field_by_id(cm->ifhcnd)->val;
      cs_real_t *cpro_fhcnc = cs_field_by_id(cm->ifhcnc)->val;
      cs_real_t *cpro_fnh3d = cs_field_by_id(cm->ifnh3d)->val;
      cs_real_t *cpro_fnh3c = cs_field_by_id(cm->ifnh3c)->val;
      cs_real_t *cpro_fnohc = cs_field_by_id(cm->ifnohc)->val;
      cs_real_t *cpro_fnonh = cs_field_by_id(cm->ifnonh)->val;
      cs_real_t *cpro_fnoch = cs_field_by_id(cm->ifnoch)->val;
      cs_real_t *cpro_fnoth = cs_field_by_id(cm->ifnoth)->val;
      cs_real_t *cpro_cnohc = cs_field_by_id(cm->icnohc)->val;
      cs_real_t *cpro_cnonh = cs_field_by_id(cm->icnonh)->val;
      cs_real_t *cpro_cnorb = cs_field_by_id(cm->icnorb)->val;
      cs_real_t *cpro_grb = cs_field_by_id(cm->igrb)->val;

      cs_arrays_set_value<cs_real_t, 1>(n_cells, 0.,
                                        cpro_ghcn1, cpro_ghcn2,
                                        cpro_gnoth,
                                        cpro_gnh31, cpro_gnh32,
                                        cpro_fhcnr, cpro_fhcnd,
                                        cpro_fhcnc,
                                        cpro_fnh3d, cpro_fnh3c,
                                        cpro_fnohc, cpro_fnonh,
                                        cpro_fnoch, cpro_fnoth,
                                        cpro_cnohc, cpro_cnonh,
                                        cpro_cnorb,  cpro_grb);
    }
  }

  /* Calculation of balances in C, O et H
     ------------------------------------ */

  cs_real_t *cpro_bcarbone = cs_field_by_id(cm->ibcarbone)->val;
  cs_real_t *cpro_boxygen = cs_field_by_id(cm->iboxygen)->val;
  cs_real_t *cpro_bhydrogen = cs_field_by_id(cm->ibhydrogen)->val;

  # pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {

    cpro_bcarbone[c_id] = cpro_x1[c_id]
                          * (  cpro_cyf1[c_id]*wmolat_iatc/wmchx1[c_id]
                             + cpro_cyf2[c_id]*wmolat_iatc/wmchx2[c_id]
                             + cpro_cyf3[c_id]*wmolat_iatc/wmole[ico]
                             + cpro_cyf6[c_id]*wmolat_iatc/wmole[ihcn]
                             + cpro_cyp1[c_id]*wmolat_iatc/wmole[ico2]);

    cpro_boxygen[c_id] = cpro_x1[c_id]
                         * (  cpro_cyf3[c_id]*    wmolat_iato/wmole[ico]
                            + cpro_cyox[c_id]*2.0*wmolat_iato/wmole[io2]
                            + cpro_cyp1[c_id]*2.0*wmolat_iato/wmole[ico2]
                            + cpro_cyp2[c_id]*    wmolat_iato/wmole[ih2o]
                            + cpro_cyp3[c_id]*2.0*wmolat_iato/wmole[iso2]);

    cpro_bhydrogen[c_id] = cpro_x1[c_id]
                           *(  cpro_cyf1[c_id]*cx1m[c_id]*wmolat_iath/wmchx1[c_id]
                             + cpro_cyf2[c_id]*cx2m[c_id]*wmolat_iath/wmchx2[c_id]
                             + cpro_cyf4[c_id]*2.0       *wmolat_iath/wmole[ih2s]
                             + cpro_cyf5[c_id]*2.0       *wmolat_iath/wmole[ihy]
                             + cpro_cyf6[c_id]*           wmolat_iath/wmole[ihcn]
                             + cpro_cyf7[c_id]*3.0       *wmolat_iath/wmole[inh3]
                             + cpro_cyp2[c_id]*2.0       *wmolat_iath/wmole[ih2o]);

  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const int coal_id = cm->ichcor[class_id] -1;

    cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
    cs_real_t *cvar_xwtcl = nullptr;
    if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] == 1)
      cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;

    // Calculation of total X2 per zone

    cs_real_t somch =   cm->cch[coal_id]
                      + cm->hch[coal_id]
                      + cm->och[coal_id]
                      + cm->sch[coal_id]
                      + cm->nch[coal_id];

    cs_real_t somck =   cm->cck[coal_id]
                      + cm->hck[coal_id]
                      + cm->ock[coal_id]
                      + cm->sck[coal_id]
                      + cm->nck[coal_id];

    cs_real_t chxc  = cm->cch[coal_id]/somch;
    cs_real_t chxh  = cm->hch[coal_id]/somch;
    cs_real_t chxo  = cm->och[coal_id]/somch;

    cs_real_t ckxc = 1., ckxh = 0., ckxo = 0.;
    if (somck > 0.) {
      ckxc  = cm->cck[coal_id]/somck;
      ckxh  = cm->hck[coal_id]/somck;
      ckxo  = cm->ock[coal_id]/somck;
    }

    # pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (auto c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t xch  = cvar_xchcl[c_id];
      cs_real_t xck  = cvar_xckcl[c_id];
      cs_real_t xwat = (cvar_xwtcl != nullptr) ? cvar_xwtcl[c_id] : 0.0;

      cpro_bcarbone[c_id] += xch*chxc + xck*ckxc;

      cpro_boxygen[c_id] +=   xch*chxo
                            + xck*ckxo
                            + xwat*wmolat_iato/wmole[ih2o];

      cpro_bhydrogen[c_id] +=   xch*chxh
                              + xck*ckxh
                              + xwat*2.0 *wmolat_iath/wmole[ih2o];

    }
  }

  // Free work arrays
  CS_FREE(intpdf);
  CS_FREE(ffuel);
  CS_FREE(dfuel);
  CS_FREE(doxyd);
  CS_FREE(pdfm1);
  CS_FREE(pdfm2);
  CS_FREE(hrec);
  CS_FREE(cx1m);
  CS_FREE(cx2m);
  CS_FREE(wmchx1);
  CS_FREE(wmchx2);
  CS_FREE(af1);
  CS_FREE(af2);

  if (cm->ieqnox == 1) {
    CS_FREE(fs3no);
    CS_FREE(fs4no);
    CS_FREE(yfs4no);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the physical properties of the dispersed phase
 *        (classes of particules)
 *
 * Cell values
 * -----------
 * - Mass fraction of solid
 *   and eventual clippings
 * - Diameter
 * - Mass density
 *   and eventual clippings
 *
 * \param[in]  n_cells   number of cells in mesh
 */
/*----------------------------------------------------------------------------*/

static void
_physprop2(cs_lnum_t  n_cells)
{
  constexpr cs_real_t c_1ov3 = 1./3.;

  // Coefficient relating to the coke diameter
  const cs_real_t coedmi = 1.2;

  const cs_coal_model_t *cm = cs_glob_coal_model;

  bool log_is_active = cs_log_default_is_active();

  /* Calculation for each class
     --------------------------
      - of the solid mass fraction
      - of the coke diameter
      - of the coal mass density */

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    cs_gnum_t n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0,  n6 = 0, n7 = 0, n8 = 0;
    cs_real_t x2min  = HUGE_VAL, x2max  = -HUGE_VAL;
    cs_real_t dchmin = HUGE_VAL, dchmax = -HUGE_VAL;
    cs_real_t dckmin = HUGE_VAL, dckmax = -HUGE_VAL;
    cs_real_t romin  = HUGE_VAL, romax  = -HUGE_VAL;

    const int coal_id = cm->ichcor[class_id] - 1;

    char f_name[32];

    const cs_real_t *nagcpi = nullptr;
    cs_real_t *agecpi = nullptr;
    if (cm->idrift >= 1) {
      snprintf(f_name, 31, "n_p_age_%02d", class_id+1); f_name[31] = '\0';
      nagcpi = cs_field_by_name(f_name)->val;

      snprintf(f_name, 31, "age_p_%02d", class_id+1); f_name[31] = '\0';
      agecpi = cs_field_by_name(f_name)->val;
    }

    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
    const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[class_id])->val;
    const cs_real_t *cvar_xwtcl = nullptr;
    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
      cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;

    const cs_real_t diam20  = cm->diam20[class_id];
    const cs_real_t rho20  = cm->rho20[class_id];
    const cs_real_t xmash = cm->xmash[class_id];
    const cs_real_t xmp0  = cm->xmp0[class_id];

    const cs_real_t rhock  = cm->rhock[coal_id];
    const cs_real_t xashcl = cm->xashch[coal_id];
    const cs_real_t pi = cs_math_pi;

    cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
    cs_real_t *cpro_rom2 = cs_field_by_id(cm->irom2[class_id])->val;

    cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      const cs_real_t xck    = cvar_xckcl[c_id];
      const cs_real_t xch    = cvar_xchcl[c_id];
      const cs_real_t xnp    = cvar_xnpcl[c_id];
      const cs_real_t xuash = xnp*(1.-xashcl)*xmp0;

      //  Calculation of the solid mass fraction
      cpro_x2[c_id] = xch + xck + xnp*xmash;
      // Taking into account the humidity
      if (cvar_xwtcl != nullptr)
        cpro_x2[c_id] += cvar_xwtcl[c_id];

      // Eventual clipping for the solid mass fraction
      if (cpro_x2[c_id] > (1. + cs_coal_epsilon)) {
        n1++;
        x2max = fmax(cpro_x2[c_id], x2max);
        cpro_x2[c_id] = 1.;
      }
      else if (cpro_x2[c_id] < -cs_coal_epsilon) {
        n2++;
        x2min = fmin(cpro_x2[c_id], x2min);
        cpro_x2[c_id] = 0.;
      }

      // Initialization

      cpro_rom2[c_id] = rho20;
      cpro_diam2[c_id] = diam20;

      cs_real_t dch;

      if (xuash > cs_coal_epsilon) {

        // Calculation of the reactive coal diameter: Dch
        dch = diam20 * pow(xch/xuash, c_1ov3);

        // Eventual clipping for the reactive coal diameter
        if (dch > (diam20 + cs_coal_epsilon)) {
          n3++;
          dchmax = fmax(dch, dchmax);
          dch = diam20;
        }
        else if (dch < - cs_coal_epsilon) {
          n4++;
          dchmin = fmin(dch, dchmin);
          dch = 0.;
        }

        // Calculation of the coke diameter: Dck stores in cpro_diam2[c_id]
        cs_real_t dck = pow(  (xch/rho20 + xck/rhock)
                            / ((1.-xashcl)*pi/6.0*xnp), c_1ov3);

        // Eventual clipping for the coke diameter
        if (dck >  coedmi*diam20) {
          n5++;
          dckmax = fmax(dck, dckmax);
          dck = diam20*coedmi;
        }
        else if (dck < -cs_coal_epsilon) {
          n6++;
          dckmin = fmin(dck, dckmin);
          dck = 0.;
        }
        cpro_diam2[c_id] = dck;

        // Density

        cs_real_t ro2ini = rho20;
        // Taking into account the humidity
        if (cvar_xwtcl != nullptr) {
          // at the moment we asume that ROH2O is constant
          constexpr cs_real_t roh2o = 998.203;
          ro2ini += cvar_xwtcl[c_id] * roh2o;
        }

        const cs_real_t diam20_3 = cs_math_pow3(diam20);
        const cs_real_t dch_3 = cs_math_pow3(dch);
        const cs_real_t dck_3 = cs_math_pow3(dck);

        cpro_rom2[c_id]
          =   (      xashcl  * diam20_3 * rho20
               + (1.-xashcl) * (dck_3 - dch_3) * rhock
               + (1.-xashcl) * dch_3 * ro2ini)
            / (      xashcl  * diam20_3
               + (1.-xashcl) * dck_3);

        // Clipping for density

        if (cpro_rom2[c_id] > (ro2ini + cs_coal_epsilon)) {
          n7++;
          romax = fmax(cpro_rom2[c_id], romax);
          cpro_rom2[c_id] = rho20;
        }
        if (cpro_rom2[c_id] < (rhock - cs_coal_epsilon)) {
          n8++;
          romin = fmin(cpro_rom2[c_id], romin);
          cpro_rom2[c_id] = rhock;
        }

      } // (xuash > cs_coal_epsilon)

      // Particles' age of each particle class
      if (agecpi != nullptr) {
        if (xnp >= cs_coal_epsilon)
          agecpi[c_id] = nagcpi[c_id]/xnp;
        else
          agecpi[c_id] = 0.0;
      }

    } // Loop on cells

    if (log_is_active) {

      cs_parall_sum_scalars(n1, n2, n3, n4, n5, n6, n7, n8);
      cs_parall_max_scalars(x2max, dchmax, dckmax, romax);
      cs_parall_min_scalars(x2min, dchmin, dckmin, romin);

      if (n1 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on max for solid mass frac. for class %d\n"
                        "    number of points: %lu\n"
                        "    max. value:       %g\n"),
                      class_id+1, (unsigned long)n1, x2max);
      if (n2 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on min for solid mass frac. for class %d\n"
                        "    number of points: %lu\n"
                        "    min. value:       %g\n"),
                      class_id+1, (unsigned long)n2, x2min);
      if (n3 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on max on coal diameter for class %d\n"
                        "    number of points: %lu\n"
                        "    max. value:       %g\n"),
                      class_id+1, (unsigned long)n3, dchmax);
      if (n4 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on min on coal diameter for class %d\n"
                        "    number of points: %lu\n"
                        "    min. value:       %g\n"),
                      class_id+1, (unsigned long)n4, dchmin);
      if (n5 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on max on coke diameter for class %d\n"
                        "    number of points: %lu\n"
                        "    max. value:       %g\n"),
                      class_id+1, (unsigned long)n5, dckmax);
      if (n6 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on min on coke diameter for class %d\n"
                        "    number of points: %lu\n"
                        "    min. value:       %g\n"),
                      class_id+1, (unsigned long)n6, dckmin);
      if (n7 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on max on mass density for class %d\n"
                        "    number of points: %lu\n"
                        "    max. value:       %g\n"),
                      class_id+1, (unsigned long)n7, romax);
      if (n8 > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" clipping on min on mass density for class %d\n"
                        "    number of points: %lu\n"
                        "    max. value:       %g\n"),
                      class_id+1, (unsigned long)n8, romin);

    } // Log is active

  } // Loop on classes
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definition of physical variable laws for combustion with a drift.
 *
 * \param[in]  n_cells   number of cells in mesh
 */
/*----------------------------------------------------------------------------*/

static void
_physical_properties_combustion_drift(cs_lnum_t  n_cells)
{
  /* Initializations to keep
     ----------------------- */

  const cs_coal_model_t *cm = cs_glob_coal_model;

  cs_real_t *visco;
  CS_MALLOC(visco, n_cells, cs_real_t);

  const cs_real_t *cpro_ym1_3 = cs_field_by_id(cm->iym1[3-1])->val;
  const cs_real_t *cpro_ym1_5 = cs_field_by_id(cm->iym1[5-1])->val;
  const cs_real_t *cpro_ym1_7 = cs_field_by_id(cm->iym1[7-1])->val;
  const cs_real_t *cpro_ym1_8 = cs_field_by_id(cm->iym1[8-1])->val;
  const cs_real_t *cpro_ym1_9 = cs_field_by_id(cm->iym1[9-1])->val;
  const cs_real_t *cpro_ym1_11 = cs_field_by_id(cm->iym1[11-1])->val;
  const cs_real_t *cpro_ym1_12 = cs_field_by_id(cm->iym1[12-1])->val;

  // Key id for drift scalar and coal scalar class
  const int keydri = cs_field_key_id("drift_scalar_model");
  const int keyccl = cs_field_key_id("scalar_class");

  const int n_fields = cs_field_n_fields();

  const cs_real_t *cpro_temp = cs_field_by_name("temperature")->val;

  // gas density
  cs_real_t *cpro_rom1 = cs_field_by_id(cm->irom1)->val;

  cs_host_context ctx;

  // First initialization
  if (cs_glob_time_step->nt_cur <= 1) {
    const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;
    const cs_real_t ro0 = cs_glob_fluid_properties->ro0;

    cs_arrays_set_value<cs_real_t, 1>(n_cells, viscl0, visco);
    cs_arrays_set_value<cs_real_t, 1>(n_cells, ro0, cpro_rom1);

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
      cs_real_t *cpro_rom2 = cs_field_by_id(cm->irom2[class_id])->val;

      cs_arrays_set_value<cs_real_t, 1>
        (n_cells, cm->rho20[class_id], cpro_rom2);
      cs_arrays_set_value<cs_real_t, 1>
        (n_cells, cm->diam20[class_id], cpro_diam2);
    }
  }

  /* Gas viscosity function of temperature
     ------------------------------------- */

  // 1-O2 2-CO 3-H2 4-N2 5-SO2 6-NH3 7-CO2

  double aa1 = 4.0495e-6;
  double bb1 = 6.22e-8;
  double cc1 = -2.3032e-11;
  double dd1 = 4.4077e-15;

  double aa2 = 9.9987e-6;
  double bb2 = 5.1578e-8;
  double cc2 = -1.8383e-11;
  double dd2 = 3.33307e-15;

  double aa3 = 2.894e-6;
  double bb3 = 2.22508e-8;
  double cc3 = -8.041e-12;
  double dd3 = 1.4619e-15;

  double aa4 = 4.3093e-6;
  double bb4 = 5.0516e-8;
  double cc4 = -1.7869e-11;
  double dd4 = 3.2136e-15;

  double aa5 = -1.9889e-6;
  double bb5 = 5.365e-8;
  double cc5 = -1.4286e-11;
  double dd5 = 2.1639e-15;

  double aa6 = -1.293e-6;
  double bb6 = 4.1194e-8;
  double cc6 = -1.772e-11;
  double dd6 = 1.8699e-15;

  double aa7 = 4.4822e-7;
  double bb7 = 5.4327e-8;
  double cc7 = -1.7581e-11;
  double dd7 = 2.9979e-15;

  /*
   *   law                 mu   = a + b T + c T**2 + d T**3
   *   so      cpro_viscl[c_id] = a +b*xvart+c*xvart**2 + d*xvart**3
   */

  if (cs_glob_time_step->nt_cur <= 1) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

      cs_real_t xvart = cpro_temp[c_id];
      cs_real_t xvart_2 = cs_math_pow2(xvart);
      cs_real_t xvart_3 = cs_math_pow3(xvart);

      cs_real_t visco_o2  = aa1 + xvart*bb1 + cc1*xvart_2 + dd1*xvart_3;
      cs_real_t visco_co  = aa2 + xvart*bb2 + cc2*xvart_2 + dd2*xvart_3;
      cs_real_t visco_h2  = aa3 + xvart*bb3 + cc3*xvart_2 + dd3*xvart_3;
      cs_real_t visco_n2  = aa4 + xvart*bb4 + cc4*xvart_2 + dd4*xvart_3;
      cs_real_t visco_so2 = aa5 + xvart*bb5 + cc5*xvart_2 + dd5*xvart_3;
      cs_real_t visco_nh3 = aa6 + xvart*bb6 + cc6*xvart_2 + dd6*xvart_3;
      cs_real_t visco_co2 = aa7 + xvart*bb7 + cc7*xvart_2 + dd7*xvart_3;

      // Viscosity of the mixing
      visco[c_id] =   (  cpro_ym1_8[c_id] * visco_o2
                       + cpro_ym1_3[c_id] * visco_co
                       + cpro_ym1_5[c_id] * visco_h2
                       + cpro_ym1_12[c_id]* visco_n2
                       + cpro_ym1_11[c_id]* visco_so2
                       + cpro_ym1_7[c_id] * visco_nh3
                       + cpro_ym1_9[c_id] * visco_co2)
                    / (  cpro_ym1_8[c_id] + cpro_ym1_3[c_id]
                       + cpro_ym1_5[c_id] + cpro_ym1_12[c_id]
                       + cpro_ym1_11[c_id]+ cpro_ym1_7[c_id]
                       + cpro_ym1_9[c_id] + 1.e-12);

    });
  }

  /* Loop over coal particle classes
   * We only handle here coal class with a drift
   * ------------------------------------------- */

  for (int fld_id = 0; fld_id < n_fields; fld_id++) {

    const cs_field_t *fld = cs_field_by_id(fld_id);

    // Index of the scalar class (<0 if the scalar belongs to the gas phase)
    int class_id = cs_field_get_key_int(fld, keyccl) - 1;
    int isccdri = cs_field_get_key_int(fld, keydri);

    // We only handle here one scalar with a drift per particle class
    if (class_id < 0 || (isccdri & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX) == 0)
      continue;

    const cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
    const cs_real_t *cpro_rom2 = cs_field_by_id(cm->irom2[class_id])->val;

    cs_real_t *cpro_taup
      = cs_field_by_composite_name(fld->name, "drift_tau")->val;

    // Computation of the relaxation time of the particles
    // the drift is therefore v_g = tau_p * g
    // Here, it is possible to include differential thermo (or electro) phoresis
    //----------------------------------------------------

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cpro_taup[c_id] = 0.;
      if (visco[c_id] > 1e-17) {
        // Simple model for Low Reynolds Numbers
        cpro_taup[c_id] =   cpro_rom2[c_id] * cs_math_pow2(cpro_diam2[c_id])
                          / (18.0 * visco[c_id]);
      }
    });

  } // Loop on fields

  CS_FREE(visco);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute \f$ \rho \f$ of the pulverized coal combustion mixture.
 *
 * \param[in, out]   mbrom    filling indicator of romb
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_physprop(int  *mbrom)
{
  static int ipass = 0; // call counter

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;

  const cs_coal_model_t *cm = cs_glob_coal_model;

  // Additional properties for drift

  if (cm->idrift > 1)
    _physical_properties_combustion_drift(n_cells);

  ipass++;

  /* Initialization
   * -------------- */

  /* Aliases for simpler syntax */

  const int n_c_classes = cm->nclacp;

  const cs_real_t *cvar_hgas = cs_field_by_id(cm->ihgas)->val;

  const cs_real_t *cpro_x2b[CS_COMBUSTION_COAL_MAX_CLASSES];
  const cs_real_t *cpro_ro2[CS_COMBUSTION_COAL_MAX_CLASSES];
  for (int class_id = 0; class_id < CS_COMBUSTION_COAL_MAX_CLASSES; class_id++) {
    if (class_id < cm->nclacp) {
      cpro_x2b[class_id] = cs_field_by_id(cm->ix2[class_id])->val;
      cpro_ro2[class_id] = cs_field_by_id(cm->irom2[class_id])->val;
    }
    else {
      cpro_x2b[class_id] = nullptr;
      cpro_ro2[class_id] = nullptr;
    }
  }

  /* Allocate arrays */

  cs_real_t *f1m, *f2m;
  CS_MALLOC(f1m, n_cells_ext, cs_real_t);
  CS_MALLOC(f2m, n_cells_ext, cs_real_t);
  cs_real_t *f3m, *f4m, *f5m, *f6m, *f7m, *f8m, *f9m;

  CS_MALLOC(f3m, n_cells_ext, cs_real_t);
  CS_MALLOC(f4m, n_cells_ext, cs_real_t);
  CS_MALLOC(f5m, n_cells_ext, cs_real_t);
  CS_MALLOC(f6m, n_cells_ext, cs_real_t);
  CS_MALLOC(f7m, n_cells_ext, cs_real_t);
  CS_MALLOC(f8m, n_cells_ext, cs_real_t);
  CS_MALLOC(f9m, n_cells_ext, cs_real_t);

  cs_real_t *enth1, *fvp2m;
  CS_MALLOC(enth1, n_cells, cs_real_t);
  CS_MALLOC(fvp2m, n_cells, cs_real_t);

  /* Calculation of the physical properties of the dispersed phase
   * -------------------------------------------------------------
   *
   * cell values
   * -----------
   *   Mass fraction of solid
   *   Diameter
   *   Mass density */

  _physprop2(n_cells);

  /* Calculation of the physical properties of the gaseous phase
   * -----------------------------------------------------------
   *
   * cell values
   * -----------
   *   Temperature
   *   Mass density
   *   Concentrations of the gaseous species */

  cs_host_context ctx;

  // Mass fraction of gas

  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;
  cs_arrays_set_value<cs_real_t, 1>(n_cells, 1., cpro_x1);
  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cs_axpy(n_cells, -1, cpro_x2b[class_id], cpro_x1);
  }
  cs_halo_sync(mesh->halo, false, cpro_x1);

  // Calculation of the gas enthalpy  enth1
  //   of F1M
  //   of F2M
  //   of F3M      in W3=1-F1M-F2M-F4M-F5M-F6M-F7M-F8M-F9M
  //   of F4M
  //   of F5M
  //   of F6M
  //   of F7M
  //   of F8M
  //   of F9M
  //   of FVP2M

  cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., f1m, f2m);

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
    const cs_real_t *cvar_f1m = cs_field_by_id(cm->if1m[coal_id])->val;
    const cs_real_t *cvar_f2m = cs_field_by_id(cm->if2m[coal_id])->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      f1m[c_id] += cvar_f1m[c_id];
      f2m[c_id] += cvar_f2m[c_id];
    });
  }

  cs_real_t *enthox = nullptr;
  if (cm->ieqnox == 1) {
    const cs_real_t *cvar_hox = cs_field_by_id(cm->ihox)->val;
    CS_MALLOC(enthox, n_cells, cs_real_t);

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      const cs_real_t xoxyd = cpro_x1[c_id] - f1m[c_id] - f2m[c_id];
      enthox[c_id] = cvar_hox[c_id] / xoxyd;
    });
  }

  cs_array_real_copy(n_cells, cs_field_by_id(cm->if7m)->val, f7m);
  cs_array_real_copy(n_cells, cs_field_by_id(cm->ifvp2m)->val, fvp2m);

  cs_arrays_set_value<cs_real_t, 1>(n_cells, 0.,
                                    f4m, f5m, f6m, f8m, f9m);

  if (cm->noxyd >= 2) {
    const cs_real_t *cvar_f4m = cs_field_by_id(cm->if4m)->val;
    cs_array_real_copy(n_cells, cvar_f4m, f4m);
    if (cm->noxyd == 3) {
      const cs_real_t *cvar_f5m = cs_field_by_id(cm->if5m)->val;
      cs_array_real_copy(n_cells, cvar_f5m, f5m);
    }
  }
  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
    const cs_real_t *cvar_f6m = cs_field_by_id(cm->if6m)->val;
    cs_array_real_copy(n_cells, cvar_f6m, f6m);
  }

  if (cm->ihtco2 == 1) {
    const cs_real_t *cvar_f8m = cs_field_by_id(cm->if8m)->val;
    cs_array_real_copy(n_cells, cvar_f8m, f8m);
  }
  if (cm->ihth2o == 1) {
    const cs_real_t *cvar_f9m = cs_field_by_id(cm->if9m)->val;
    cs_array_real_copy(n_cells, cvar_f9m, f9m);
  }

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    const cs_real_t uns1pw = 1. / cpro_x1[c_id];

    // Units: [kg scalars / kg gas]
    f1m[c_id]  *= uns1pw;
    f2m[c_id]  *= uns1pw;
    f4m[c_id]  *= uns1pw;
    f5m[c_id]  *= uns1pw;
    f6m[c_id]  *= uns1pw;
    f7m[c_id]  *= uns1pw;
    f8m[c_id]  *= uns1pw;
    f9m[c_id]  *= uns1pw;

    fvp2m[c_id ] *= uns1pw;

    f3m[c_id] = 1. - (  f1m[c_id] + f2m[c_id] + f4m[c_id] + f5m[c_id]
                      + f6m[c_id] + f7m[c_id] + f8m[c_id] + f9m[c_id]);
  });

  // Gas Enthalpy h1 (cpro_x1 h1 is transported)
  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    enth1[c_id] = cvar_hgas[c_id] / cpro_x1[c_id];
  });

  // ctx.wait();

  /* Logging */

  if (cs_log_default_is_active()) {
    cs_real_t ff3min = HUGE_VAL, ff3max = -HUGE_VAL;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      ff3max = cs::max(ff3max, f3m[c_id]);
      ff3min = cs::min(ff3min, f3m[c_id]);
    }

    cs_log_printf(CS_LOG_DEFAULT,
                  _("\n Values of F3 min and max: %g %g\n"),
                  ff3min, ff3max);
  }

  cs_real_t *cpro_rom1 = cs_field_by_id(cm->irom1)->val;

  _physprop1(f1m, f2m, f3m, f4m, f5m, f6m, f7m, f8m, f9m,
             fvp2m, enth1, enthox, cpro_rom1);

  CS_FREE(f1m);
  CS_FREE(f2m);
  CS_FREE(f3m);
  CS_FREE(f4m);
  CS_FREE(f5m);
  CS_FREE(f6m);
  CS_FREE(f7m);
  CS_FREE(f8m);
  CS_FREE(f9m);

  CS_FREE(fvp2m);
  CS_FREE(enth1);
  CS_FREE(enthox);

  /* Calculation of the physical properties of the dispersed phase
   * -------------------------------------------------------------
   *
   * cell values
   * -----------
   *   Temperature */

  //Transport of H2

  cs_coal_ht_convert_h_to_t_particles();

  /* Calculation of the physical properties of the mixture
   * -----------------------------------------------------
   *
   * cell values
   * -----------
   *   Mass density */

  // Calculation of Rho of the mixture: 1/Rho = X1/Rho1 + Sum(X2/Rho2)
  // We relax when we have a rho n available, ie
  // from the second pass or
  // from the first pass if we are in continuation of the calculation and
  // that we have re-read the mass density in the restart file.

  cs_real_t *crom = CS_F_(rho)->val;

  cs_real_t srrom1 = cm->srrom;
  if (ipass == 1) {
    if (cs_restart_present()) {
      if (   cs_restart_get_field_read_status(CS_F_(rho)->id) == 0
          || cs_restart_get_field_read_status(CS_F_(rho_b)->id) == 0)
        srrom1 = 0.;
    }
    else
      srrom1 = 0.;
  }

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    cs_real_t x2sro2 = 0.;
    for (int class_id = 0; class_id < n_c_classes; class_id++) {
      x2sro2 += cpro_x2b[class_id][c_id] / cpro_ro2[class_id][c_id];
    }
    cs_real_t x1sro1 = cpro_x1[c_id] / cpro_rom1[c_id];
    // Eventual relaxation
    crom[c_id] = srrom1*crom[c_id] + (1.-srrom1)/(x1sro1+x2sro2);
  });

  /* Calculation of the physical properties of the mixture
   * -----------------------------------------------------
   *
   * face values
   * ----------- */

  *mbrom = 1;

  cs_coal_boundary_conditions_density();

  /* Compute the drift velocity if needed
   * ------------------------------------ */

  cs_real_3_t *vdc = nullptr;

  if (cm->idrift >= 1) {

    const int n_fields = cs_field_n_fields();

    // Key id for drift scalar and coal scalar class
    const int keydri = cs_field_key_id("drift_scalar_model");
    const int keyccl = cs_field_key_id("scalar_class");

    const cs_real_t gravity[3] = {cs_glob_physical_constants->gravity[0],
                                  cs_glob_physical_constants->gravity[1],
                                  cs_glob_physical_constants->gravity[2]};

    /* Compute the limit velocity
       -------------------------- */

    // Loop over coal particle classes
    // We only handle here coal class with a drift

    for (int fld_id = 0; fld_id < n_fields; fld_id++) {

      const cs_field_t *fld = cs_field_by_id(fld_id);

      // Index of the scalar class (<0 if the scalar belongs to the gas phase)
      int class_id = cs_field_get_key_int(fld, keyccl) - 1;
      int isccdri = cs_field_get_key_int(fld, keydri);

      // We only handle here one scalar with a drift per particle class
      if (class_id < 0 || (isccdri & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX) == 0)
        continue;

      // Position of variables, coefficients

      // Corresponding relaxtion time
      cs_real_t *cpro_taup
        = cs_field_by_composite_name(fld->name, "drift_tau")->val;

      char fld_name[32];
      snprintf(fld_name, 31, "vg_lim_p_%02d", class_id+1);
      fld_name[31] = '\0';
      cs_real_3_t *vg_lim_pi = (cs_real_3_t *)(cs_field_by_name(fld_name)->val);

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        vg_lim_pi[c_id][0] = cpro_taup[c_id] * gravity[0];
        vg_lim_pi[c_id][1] = cpro_taup[c_id] * gravity[1];
        vg_lim_pi[c_id][2] = cpro_taup[c_id] * gravity[2];
      });

    } // loop on fields

    /* Init of the drift velocity of the continuous phase (gas)
       -------------------------------------------------------- */

    vdc = (cs_real_3_t *)(cs_field_by_name("vd_c")->val);

    cs_array_real_fill_zero(n_cells*3, (cs_real_t *)vdc);
  }

  /* Transported particle velocity
     ----------------------------- */

  if (cm->idrift == 1) {

    // Get all needed fields
    const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      char class_name[8];
      snprintf(class_name, 8, "%02d", class_id+1); class_name[7] = '\0';

      cs_real_t *v_x_pi
        = cs_field_by_composite_name("v_x_p", class_name)->val;
      cs_real_t *v_y_pi
        = cs_field_by_composite_name("v_y_p", class_name)->val;
      cs_real_t *v_z_pi
        = cs_field_by_composite_name("v_z_p", class_name)->val;
      cs_real_3_t *vdp_i
        = (cs_real_3_t *)(cs_field_by_composite_name("vd_p", class_name)->val);

      cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        // Vdi = Vpi-Vs
        if (cpro_x2[c_id] > 1.e-7) {
          vdp_i[c_id][0] = v_x_pi[c_id] - cvar_vel[c_id][0];
          vdp_i[c_id][1] = v_y_pi[c_id] - cvar_vel[c_id][1];
          vdp_i[c_id][2] = v_z_pi[c_id] - cvar_vel[c_id][2];
        }
        else {
          vdp_i[c_id][0] = 0.;
          vdp_i[c_id][1] = 0.;
          vdp_i[c_id][2] = 0.;
        }
        vdc[c_id][0] -= cpro_x2[c_id] * vdp_i[c_id][0];
        vdc[c_id][1] -= cpro_x2[c_id] * vdp_i[c_id][1];
        vdc[c_id][2] -= cpro_x2[c_id] * vdp_i[c_id][2];
      });
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t x1_inv = 1. / cpro_x1[c_id];
      vdc[c_id][0] *= x1_inv;
      vdc[c_id][1] *= x1_inv;
      vdc[c_id][2] *= x1_inv;
    });

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      char class_name[32];
      snprintf(class_name, 31, "%02d", class_id+1); class_name[31] = '\0';

      cs_real_3_t *vdp_i
        = (cs_real_3_t *)(cs_field_by_composite_name("vd_p", class_name)->val);
      cs_real_3_t *vg_pi
        = (cs_real_3_t *)(cs_field_by_composite_name("vg_p", class_name)->val);

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        vg_pi[c_id][0] = vdp_i[c_id][0] - vdc[c_id][0];
        vg_pi[c_id][1] = vdp_i[c_id][1] - vdc[c_id][1];
        vg_pi[c_id][2] = vdp_i[c_id][2] - vdc[c_id][2];
      });
    }

  }

  /* Prescribed drift
     ---------------- */

  else if (cm->idrift > 1) {

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      char class_name[8];
      snprintf(class_name, 8, "%02d", class_id+1); class_name[7] = '\0';

      cs_real_3_t *vg_lim_pi
        = (cs_real_3_t *)(cs_field_by_composite_name("vg_lim_p",
                                                     class_name)->val);
      cs_real_3_t *vg_pi
        = (cs_real_3_t *)(cs_field_by_composite_name("vg_p",
                                                     class_name)->val);
      cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        // FIXME vg_ is useless!
        vg_pi[c_id][0] = vg_lim_pi[c_id][0];
        vg_pi[c_id][1] = vg_lim_pi[c_id][1];
        vg_pi[c_id][2] = vg_lim_pi[c_id][2];
        vdc[c_id][0] = vdc[c_id][0] -cpro_x2[c_id] * vg_pi[c_id][0];
        vdc[c_id][1] = vdc[c_id][1] -cpro_x2[c_id] * vg_pi[c_id][1];
        vdc[c_id][2] = vdc[c_id][2] -cpro_x2[c_id] * vg_pi[c_id][2];
      });
    }

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      char class_name[32];
      snprintf(class_name, 31, "%02d", class_id+1); class_name[31] = '\0';

      cs_real_3_t *vg_pi
        = (cs_real_3_t *)(cs_field_by_composite_name("vg_p",
                                                     class_name)->val);
      cs_real_3_t *vdp_i
        = (cs_real_3_t *)(cs_field_by_composite_name("vd_p",
                                                     class_name)->val);

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        vdp_i[c_id][0] = vdc[c_id][0] + vg_pi[c_id][0];
        vdp_i[c_id][1] = vdc[c_id][1] + vg_pi[c_id][1];
        vdp_i[c_id][2] = vdc[c_id][2] + vg_pi[c_id][2];
      });
    }

  } // End for drift

  // ctx.wait();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
