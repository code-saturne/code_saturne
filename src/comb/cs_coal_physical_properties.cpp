/*============================================================================
 * Coal combustion model: physical properties
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

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_coal.h"
#include "cs_coal_ht_convert.h"
#include "cs_coal_source_terms.h"
#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal_physical_properties.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal_noxst.c

  \brief Coal combustion model: NOx source terms computation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Prototypes for Fortran functions
 *============================================================================*/

void
CS_PROCF(pppdfr, PPPDFR)(cs_lnum_t        *ncelet,
                         cs_lnum_t        *ncel,
                         int              *intpdf,
                         cs_real_t        *tpdf,
                         const cs_real_t  *ffuel,
                         const cs_real_t  *fvp2m,
                         const cs_real_t  *fmini,
                         const cs_real_t  *fmaxi,
                         const cs_real_t  *doxyd,
                         const cs_real_t  *dfuel,
                         const cs_real_t  *pdfm1,
                         const cs_real_t  *pdfm2,
                         const cs_real_t  *hrec);

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

    // Store de fs3, fs4 and concentrations in fs4 for the NOx model

    if (cvar_yco2 != nullptr) {

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
    else
      assert(0);

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

/*============================================================================
 * Public function definitions
 *============================================================================*/

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
 * \param[in]     f1m     mean of tracer 1 mvl [chx1m+co]
 * \param[in]     f2m     mean of tracer 2 mvl [chx2m+co]
 * \param[in]     f3m     mean of tracer 3 (oxydant 1)
 * \param[in]     f4m     mean of tracer 4 (oxydant 2)
 * \param[in]     f5m     mean of tracer 5 (oxydant 3)
 * \param[in]     f6m     mean of tracer 6 (humidity)
 * \param[in]     f7m     mean of tracer 7 (C + O2)
 * \param[in]     f8m     mean of tracer 8 (C + CO2)
 * \param[in]     f9m     mean of tracer 9 (C + H2O)
 * \param[in]     fvp2m   f1f2 variance
 * \param[in]     enth    enthalpy in \f$ j . kg^{-1} \f$  either of
 *                        the gas or of the mixture
 * \param[in]     enthox  oxydant enthalpy
 * \param[out]    rom1    gas density
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_physprop1(const  cs_real_t  f1m[],
                  const  cs_real_t  f2m[],
                  const  cs_real_t  f3m[],
                  const  cs_real_t  f4m[],
                  const  cs_real_t  f5m[],
                  const  cs_real_t  f6m[],
                  const  cs_real_t  f7m[],
                  const  cs_real_t  f8m[],
                  const  cs_real_t  f9m[],
                  const  cs_real_t  fvp2m[],
                  const  cs_real_t  enth[],
                  const  cs_real_t  enthox[],
                  cs_real_t         rom1[])
{
  static int ipass = 0;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

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
  BFT_MALLOC(intpdf, n_cells, int);

  cs_real_t *ffuel, *dfuel, *doxyd, *pdfm1;
  cs_real_t *pdfm2, *hrec, *cx1m, *cx2m, *wmchx1, *wmchx2;
  cs_real_t *af1, *af2;
  BFT_MALLOC(ffuel, n_cells, cs_real_t);
  BFT_MALLOC(dfuel, n_cells, cs_real_t);
  BFT_MALLOC(doxyd, n_cells, cs_real_t);
  BFT_MALLOC(pdfm1, n_cells, cs_real_t);
  BFT_MALLOC(pdfm2, n_cells, cs_real_t);
  BFT_MALLOC(hrec, n_cells, cs_real_t);
  BFT_MALLOC(cx1m, n_cells, cs_real_t);
  BFT_MALLOC(cx2m, n_cells, cs_real_t);
  BFT_MALLOC(wmchx1, n_cells, cs_real_t);
  BFT_MALLOC(wmchx2, n_cells, cs_real_t);
  BFT_MALLOC(af1, n_cells*n_gas_sp, cs_real_t);
  BFT_MALLOC(af2, n_cells*n_gas_sp, cs_real_t);

  cs_real_t *fs3no = nullptr, *fs4no = nullptr, *yfs4no = nullptr;

  if (cm->ieqnox == 1) {
    BFT_MALLOC(fs3no, n_cells, cs_real_t);
    BFT_MALLOC(fs4no, n_cells, cs_real_t);
    BFT_MALLOC(yfs4no, n_cells*n_gas_sp, cs_real_t);
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
  BFT_MALLOC(fmini, n_cells, cs_real_t);
  BFT_MALLOC(fmaxi, n_cells, cs_real_t);
  BFT_MALLOC(tpdf, n_cells_ext, cs_real_t);

  # pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (auto c_id = 0; c_id < n_cells; c_id++) {
    // min et max boundary of the pdf
    fmini[c_id] = 0;;
    fmaxi[c_id] = 1.;
    // Sum F1+F2
    ffuel[c_id] = f1m[c_id] + f2m[c_id];
  }

  cs_lnum_t ncelet = n_cells_ext;
  cs_lnum_t ncel = n_cells;

  CS_PROCF(pppdfr, PPPDFR)(&ncelet, &ncel, intpdf, tpdf, ffuel, fvp2m,
                           fmini, fmaxi, doxyd, dfuel, pdfm1, pdfm2, hrec);

  BFT_FREE(tpdf);
  BFT_FREE(fmini);
  BFT_FREE(fmaxi);

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
  BFT_FREE(intpdf);
  BFT_FREE(ffuel);
  BFT_FREE(dfuel);
  BFT_FREE(doxyd);
  BFT_FREE(pdfm1);
  BFT_FREE(pdfm2);
  BFT_FREE(hrec);
  BFT_FREE(cx1m);
  BFT_FREE(cx2m);
  BFT_FREE(wmchx1);
  BFT_FREE(wmchx2);
  BFT_FREE(af1);
  BFT_FREE(af2);

  if (cm->ieqnox == 1) {
    BFT_FREE(fs3no);
    BFT_FREE(fs4no);
    BFT_FREE(yfs4no);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
