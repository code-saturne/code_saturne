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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

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

void
cs_gascomb(cs_lnum_t        n_cells,
           int              ichx1,
           int              ichx2,
           int              intpdf[],
           const cs_real_t  f1m[],
           const cs_real_t  f2m[],
           const cs_real_t  f3m[],
           const cs_real_t  f4m[],
           const cs_real_t  f5m[],
           const cs_real_t  f6m[],
           const cs_real_t  f7m[],
           const cs_real_t  f8m[],
           const cs_real_t  f9m[],
           cs_real_t        pdfm1[],
           cs_real_t        pdfm2[],
           cs_real_t        doxyd[],
           cs_real_t        dfuel[],
           cs_real_t        hrec[],
           cs_real_t        af1[],
           cs_real_t        af2[],
           cs_real_t        cx1m[],
           cs_real_t        cx2m[],
           cs_real_t        wmchx1[],
           cs_real_t        wmchx2[],
           const cs_real_t  cpro_cyf1[],
           const cs_real_t  cpro_cyf2[],
           const cs_real_t  cpro_cyf3[],
           const cs_real_t  cpro_cyf4[],
           const cs_real_t  cpro_cyf5[],
           const cs_real_t  cpro_cyf6[],
           const cs_real_t  cpro_cyf7[],
           const cs_real_t  cpro_cyox[],
           const cs_real_t  cpro_cyp1[],
           const cs_real_t  cpro_cyp2[],
           const cs_real_t  cpro_cyp3[],
           const cs_real_t  cpro_cyin[],
           cs_real_t  fs3no[],
           cs_real_t  fs4no[],
           cs_real_t  yfs4no[]);

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
  const cs_real_t *cpro_cyf1 = cs_field_by_id(cm->iym1[ichx1])->val;
  const cs_real_t *cpro_cyf2 = cs_field_by_id(cm->iym1[ichx2])->val;
  const cs_real_t *cpro_cyf3 = cs_field_by_id(cm->iym1[ico])->val;
  const cs_real_t *cpro_cyf4 = cs_field_by_id(cm->iym1[ih2s])->val;
  const cs_real_t *cpro_cyf5 = cs_field_by_id(cm->iym1[ihy])->val;
  const cs_real_t *cpro_cyf6 = cs_field_by_id(cm->iym1[ihcn])->val;
  const cs_real_t *cpro_cyf7 = cs_field_by_id(cm->iym1[inh3])->val;
  const cs_real_t *cpro_cyox = cs_field_by_id(cm->iym1[io2])->val;
  const cs_real_t *cpro_cyp1 = cs_field_by_id(cm->iym1[ico2])->val;
  const cs_real_t *cpro_cyp2 = cs_field_by_id(cm->iym1[ih2o])->val;
  const cs_real_t *cpro_cyp3 = cs_field_by_id(cm->iym1[iso2])->val;
  const cs_real_t *cpro_cyin = cs_field_by_id(cm->iym1[in2])->val;

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

  int icb1 = ichx1 + 1; // For Fortran index
  int icb2 = ichx2 + 1;

  cs_gascomb(n_cells, icb1, icb2,
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
