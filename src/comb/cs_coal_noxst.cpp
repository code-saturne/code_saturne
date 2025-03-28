/*============================================================================
 * Coal combustion model: NOx source term computation
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
#include "comb/cs_coal.h"
#include "comb/cs_coal_ht_convert.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "comb/cs_coal_source_terms.h"

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

#define CS_COAL_NOXST_NPART  200

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute NOx source terms for pulverized coal flame.
 *
 * K1 exp(-E1/RT)   (conversion HCN to N2)
 * K2 exp(-E2/RT)   (conversion HCN to NO)
 * K3 exp(-E3/RT)   (themal NO)
 *
 * \param[in]      indpdf  use pdf
 * \param[in]      pdfm1   lower bound of pdf
 * \param[in]      pdfm2   upper bound of pdf
 * \param[in]      doxyd   amplitude of Dirac peak in 1
 * \param[in]      dfuel   amplitude of Dirac peak in 0
 * \param[in]      hrec    height of pdf rectangle
 * \param[in]      f1m     mean of tracer 1 mvl [chx1m+co]
 * \param[in]      f2m     mean of tracer 2 mvl [chx2m+co]
 * \param[in]      f3m     mean of tracer 3 (oxydant 1)
 * \param[in]      f4m     mean of tracer 4 (oxydant 2)
 * \param[in]      f5m     mean of tracer 5 (oxydant 3)
 * \param[in]      f6m     mean of tracer 6 (humidity)
 * \param[in]      f7m     mean of tracer 7 (C + O2)
 * \param[in]      f8m     mean of tracer 8 (C + CO2)
 * \param[in]      f9m     mean of tracer 9 (C + H2O)
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_noxst(const  int        indpdf[],
              const  cs_real_t  pdfm1[],
              const  cs_real_t  pdfm2[],
              const  cs_real_t  doxyd[],
              const  cs_real_t  dfuel[],
              const  cs_real_t  hrec[],
              const  cs_real_t  f3m[],
              const  cs_real_t  f4m[],
              const  cs_real_t  f5m[],
              const  cs_real_t  f6m[],
              const  cs_real_t  f7m[],
              const  cs_real_t  f8m[],
              const  cs_real_t  f9m[],
              const  cs_real_t  fs3no[],
              const  cs_real_t  fs4no[],
              const  cs_real_t  yfs4no[],
              const  cs_real_t  enthox[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

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
  const int in2 = cm->in2 -1;

  /* Aliases for simpler syntax */

  const cs_real_t *wmole = cm->wmole;
  const int ngazem = CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS;
  const int n_gas_sp = cm->n_gas_species;

  /* Preliminary computations
     ======================== */

  /* Molar masses */

  const cs_real_t wmo2  = wmole[io2];

  /* Pointers */

  cs_real_t *cpro_exp1 = cs_field_by_id(cm->ighcn1)->val;
  cs_real_t *cpro_exp2 = cs_field_by_id(cm->ighcn2)->val;
  cs_real_t *cpro_exp3 = cs_field_by_id(cm->ignoth)->val;
  cs_real_t *cpro_exp4 = cs_field_by_id(cm->ignh31)->val;
  cs_real_t *cpro_exp5 = cs_field_by_id(cm->ignh32)->val;
  cs_real_t *cpro_exprb = cs_field_by_id(cm->igrb)->val;

  /* Parameters of Arrhenius laws */

  const cs_real_t kk1 = 3.0e12;
  const cs_real_t ee1 = 3.0e4;
  const cs_real_t kk2 = 1.2e10;
  const cs_real_t ee2 = 3.35e4;
  const cs_real_t kk3 = 3.4e12;
  const cs_real_t ee3 = 6.69e4;

  const cs_real_t kk4 = 4.1e6;
  const cs_real_t ee4 = 1.6e4;
  const cs_real_t kk5 = 1.8e8;
  const cs_real_t ee5 = 1.35e4;
  const cs_real_t kkrb= 2.7e12;
  const cs_real_t eerb= 9.467e3;

  /* For terms, indicator for computation by pdf
   *  = 1 --> use pdf
   *  = 0 --> do not use pdf */

  int ipdf1 = 0;
  int ipdf2 = 0;
  int ipdf3 = 1;

  int ipdf4 = 0;
  int ipdf5 = 0;

  // Initialization

  cs_array_real_fill_zero(n_cells, cpro_exp1);
  cs_array_real_fill_zero(n_cells, cpro_exp2);
  cs_array_real_fill_zero(n_cells, cpro_exp3);
  cs_array_real_fill_zero(n_cells, cpro_exp4);
  cs_array_real_fill_zero(n_cells, cpro_exp5);
  cs_array_real_fill_zero(n_cells, cpro_exprb);

  /* Computation without the PDFs
     ============================ */

  const cs_real_t *cpro_temp = CS_F_(t)->val;
  const cs_real_t *cpro_yo2 = cs_field_by_id(cm->iym1[io2])->val;
  const cs_real_t *cpro_mmel = cs_field_by_id(cm->immel)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t tg  = cpro_temp[c_id];
    cs_real_t yo2 = cpro_yo2[c_id];
    cs_real_t xo2 = yo2*cpro_mmel[c_id]/wmo2;

    // Reaction: HCN + NO + 1/4 O2 ---> N2 + 1/2 H2O + CO
    cpro_exp1[c_id]  = kk1*exp(-ee1/tg);

    // Reaction: HCN + 5/4 O2 --> NO + 1/2 H2O  + CO
    cs_real_t xo2_pow_bb;
    if (xo2 > 0.018)
      xo2_pow_bb = 1.0; // bb = 0.
    else if (xo2 < 0.0025)
      xo2_pow_bb = xo2; // bb= 1.
    else
      xo2_pow_bb = pow(xo2, ((0.018-xo2) / (0.018 - 0.0025)));

    cpro_exp2[c_id] = kk2 * exp(-ee2/tg) * xo2_pow_bb;

    // Thermal No: Zeldovich
    cpro_exp3[c_id] = kk3 * exp(-ee3/tg) * sqrt(xo2);

    // Reaction: NH3 + O2 --> NO + ...
    cpro_exp4[c_id] = kk4 * exp(-ee4/tg) * xo2_pow_bb;

    // Reaction NH3 + NO --> N2 + ...
    cpro_exp5[c_id] = kk5 *exp(-ee5/tg);

    // Reburning (Chen's model)
    cpro_exprb[c_id] = kkrb * exp(-eerb/tg);

  }

  if (ipdf1 + ipdf2 + ipdf3 + ipdf4 + ipdf5 == 0)
    return;

  /* Computation with the PDFs
     ========================= */

  // Arrays of pointers containing the field values for each class
  // (loop on cells outside loop on classes)

  /* Mark cells for which computation is needed */

  cs_lnum_t *compute_id;
  CS_MALLOC(compute_id, n_cells, cs_lnum_t);

  cs_lnum_t n_compute = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (indpdf[c_id] == 1) {
      if ((fs3no[c_id] > fs4no[c_id]) &&  (fs4no [c_id] < 1.0)) {
        compute_id[n_compute] = c_id;
        n_compute++;
      }
    }
  }

  cs_real_t *xmx2, *tfuel;
  CS_MALLOC(xmx2, n_compute, cs_real_t);
  CS_MALLOC(tfuel, n_compute, cs_real_t);

  cs_array_real_fill_zero(n_compute, xmx2);
  cs_array_real_fill_zero(n_compute, tfuel);

  /* Compute mean tfuel */

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const cs_real_t xmash = cm->xmash[class_id];
    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
    const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[class_id])->val;

    for (cs_lnum_t c_idx = 0; c_idx < n_compute; c_idx++) {
      cs_lnum_t c_id = compute_id[c_idx];
      xmx2[c_idx] += cvar_xchcl[c_id] + cvar_xckcl[c_id] + cvar_xnpcl[c_id]*xmash;
    }

    // Acount for humidity
    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
      const cs_real_t *cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;
      for (cs_lnum_t c_idx = 0; c_idx < n_compute; c_idx++) {
        cs_lnum_t c_id = compute_id[c_idx];
        xmx2[c_idx] += cvar_xwtcl[c_id];
      }
    }

  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const cs_real_t xmash = cm->xmash[class_id];
    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
    const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[class_id])->val;
    const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;

    const cs_real_t *cvar_xwtcl = NULL;
    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
      cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;

    for (cs_lnum_t c_idx = 0; c_idx < n_compute; c_idx++) {
      cs_lnum_t c_id = compute_id[c_idx];
      if (xmx2[c_idx] > 0) {
        tfuel[c_idx] += (  cvar_xchcl[c_id] + cvar_xckcl[c_id]
                         + (cvar_xnpcl[c_id] * xmash)) * cpro_temp2[c_id];

        if (cvar_xwtcl != NULL)
          tfuel[c_idx] += cvar_xwtcl[c_id] * cpro_temp2[c_id];
      }
    }

  }

  for (cs_lnum_t c_idx = 0; c_idx < n_compute; c_idx++) {
    if (xmx2[c_idx] > 0) {
      tfuel[c_idx] /= xmx2[c_idx];
    }
    else {
      cs_lnum_t c_id = compute_id[c_idx];
      tfuel[c_idx] = cpro_temp[c_id];
    }
  }

  CS_FREE(xmx2);

  cs_gnum_t inok = 0;
  cs_gnum_t i300 = 0, i000 = 0, imini= 0;
  cs_gnum_t i2500= 0, i2600= 0, i2700= 0, i2800= 0;
  cs_gnum_t i3000= 0, i3500= 0, i4000= 0, i5000= 0;
  cs_gnum_t imaxi= 0, nbpt = 0, nbclip1 = 0, nbclip2 = 0;
  cs_gnum_t nbclip30 = 0, nbclip31 = 0;

  cs_real_t ts4min = 1.e+20, ts4max = -1.e+20;
  cs_real_t summin = 1.e+20, summax = -1.e+20;
  cs_real_t ts4admin = 1.e+20, ts4admax = -1.e+20;
  cs_real_t toxmin = 1.e+20, toxmax = -1.e+20;
  cs_real_t yo2oxmin = 1.e+20, yo2oxmax = -1.e20;

  cs_real_t yo2min = 1.e+20, yo2max  = -1.e+20;
  cs_real_t yo2min1 = 1.e+20, yo2max1 = -1.e+20;

  /* Compute coal enthalpy at Tfuel */

  cs_real_t *xxf, *hhf;
  CS_MALLOC(xxf, n_compute, cs_real_t);
  CS_MALLOC(hhf, n_compute, cs_real_t);

  cs_array_real_fill_zero(n_compute, xxf);
  cs_array_real_fill_zero(n_compute, hhf);

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

    const cs_real_t *cvar_f1m_c = cs_field_by_id(cm->if1m[coal_id])->val;
    const cs_real_t *cvar_f2m_c = cs_field_by_id(cm->if2m[coal_id])->val;

    cs_real_t coefe_1[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    cs_real_t coefe_2[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    for (int ige = 0; ige < ngazem; ige++) {
      coefe_1[ige] = 0;
      coefe_2[ige] = 0;
    }

    int ichx1c_coal = cm->ichx1c[coal_id] - 1;
    int ichx2c_coal = cm->ichx2c[coal_id] - 1;

    cs_real_t den;
    den =   cm->a1[coal_id]*wmole[ichx1c_coal]
          + cm->b1[coal_id]*wmole[ico]
          + cm->c1[coal_id]*wmole[ih2o]
          + cm->d1[coal_id]*wmole[ih2s]
          + cm->e1[coal_id]*wmole[ihcn]
          + cm->f1[coal_id]*wmole[inh3];

    coefe_1[ichx1] = cm->a1[coal_id]*wmole[ichx1c_coal] / den;
    coefe_1[ico]   = cm->b1[coal_id]*wmole[ico]         / den;
    coefe_1[ih2o]  = cm->c1[coal_id]*wmole[ih2o]        / den;
    coefe_1[ih2s]  = cm->d1[coal_id]*wmole[ih2s]        / den;
    coefe_1[ihcn]  = cm->e1[coal_id]*wmole[ihcn]        / den;
    coefe_1[inh3]  = cm->f1[coal_id]*wmole[inh3]        / den;

    den =   cm->a2[coal_id]*wmole[ichx2c_coal]
          + cm->b2[coal_id]*wmole[ico]
          + cm->c2[coal_id]*wmole[ih2o]
          + cm->d2[coal_id]*wmole[ih2s]
          + cm->e2[coal_id]*wmole[ihcn]
          + cm->f2[coal_id]*wmole[inh3];

    coefe_2[ichx2] = cm->a2[coal_id]*wmole[ichx2c_coal] / den;
    coefe_2[ico]   = cm->b2[coal_id]*wmole[ico]         / den;
    coefe_2[ih2o]  = cm->c2[coal_id]*wmole[ih2o]        / den;
    coefe_2[ih2s]  = cm->d2[coal_id]*wmole[ih2s]        / den;
    coefe_2[ihcn]  = cm->e2[coal_id]*wmole[ihcn]        / den;
    coefe_2[inh3]  = cm->f2[coal_id]*wmole[inh3]        / den;

    cs_real_t f1mc[CS_COMBUSTION_MAX_COALS], f2mc[CS_COMBUSTION_MAX_COALS];
    for (int icha = 0; icha < CS_COMBUSTION_MAX_COALS; icha++) {
      f1mc[icha] = 0;
      f2mc[icha] = 0;
    }
    f1mc[coal_id] = 1.;
    f2mc[coal_id] = 1.;

    // H(mv1, TFUEL) and H(mv2, TFUEL)

    for (cs_lnum_t c_idx = 0; c_idx < n_compute; c_idx++) {
      cs_lnum_t c_id = compute_id[c_idx];

      const cs_real_t _tfuel = tfuel[c_idx];

      cs_real_t xhf1
        = cs_coal_ht_convert_t_to_h_gas_by_yi_f1f2(_tfuel, coefe_1,
                                                   f1mc, f2mc);

      cs_real_t xhf2
        = cs_coal_ht_convert_t_to_h_gas_by_yi_f1f2(_tfuel, coefe_2,
                                                   f1mc, f2mc);

      xxf[c_idx] += cvar_f1m_c[c_id] + cvar_f2m_c[c_id];
      hhf[c_idx] +=   cvar_f1m_c[c_id]*xhf1
                    + cvar_f2m_c[c_id]*xhf2;

    }

  } // loop on coals

  for (cs_lnum_t c_idx = 0; c_idx < n_compute; c_idx++) {
    cs_lnum_t c_id = compute_id[c_idx];

    const cs_real_t _tfuel = tfuel[c_idx];

    // Compute Yo2 in oxydant
    //         Yo2 in fs4

    cs_real_t bb1 = fmax(0., pdfm1[c_id]);
    cs_real_t bb2 = fmin(fs3no[c_id], pdfm2[c_id]);
    cs_real_t lro = (bb2 > bb1) ? hrec[c_id] : 0.;
    cs_real_t qqq = bb2*bb2 - bb1*bb1;
    cs_real_t rrr = bb2     - bb1;
    cs_real_t gt1 = lro*qqq / (2.*fs3no[c_id]);
    cs_real_t gt2 = lro*rrr;
    cs_real_t gt3 = doxyd[c_id];
    cs_real_t yo2cb = 0.;
    cs_real_t yo2ox = 0.;
    cs_real_t yo2s4 = 0.;

    if (cpro_yo2[c_id] > 0.) {

      yo2ox = cpro_yo2[c_id] / (-gt1+gt2+gt3);

      yo2oxmin = fmin(yo2oxmin, yo2ox);
      yo2oxmax = fmax(yo2oxmax, yo2ox);

      cs_real_t yo2moy = cpro_yo2[c_id];
      cs_real_t dirac  = dfuel[c_id]*yo2cb + doxyd[c_id]*yo2ox;

      bb1 = fmax(0., pdfm1[c_id]);
      bb2 = fmin(fs4no[c_id], pdfm2[c_id]);
      cs_real_t bb3 = fmax(fs4no[c_id], pdfm1[c_id]);
      cs_real_t bb4 = fmin(fs3no[c_id], pdfm2[c_id]);
      if (bb2 > bb1)
        lro = hrec[c_id];
      else
        lro = 0.;
      cs_real_t lrf = (bb4 > bb3) ? hrec[c_id] : 0.;

      qqq = bb2*bb2 - bb1*bb1;
      rrr = bb2     - bb1;
      cs_real_t sss = bb4*bb4 - bb3*bb3;
      cs_real_t ttt = bb4     - bb3;
      cs_real_t uuu = fs4no[c_id] - fs3no[c_id];

      gt1 = lro*qqq / (2.0*fs4no[c_id]);
      gt2 = lro*rrr;

      cs_real_t gt10 = lrf*sss / (2.0*uuu);
      cs_real_t gt20 = lrf*ttt*fs3no[c_id] / uuu;

      cs_real_t yo24num = yo2moy - dirac + yo2ox*(gt1 -gt2);
      cs_real_t yo24den = gt1+gt10-gt20;

      yo2s4  = yo24num/yo24den;

    }

    // Logging and clipping

    yo2min = fmin(yo2s4, yo2min);
    yo2max = fmax(yo2s4, yo2max);
    if (yo2s4 < 0) {
      nbclip30++;
      yo2s4 = 0.;
    }
    if (yo2s4 > yo2ox) {
      nbclip31++;
      yo2s4 = yo2ox;
    }
    yo2min1 = fmin(yo2s4, yo2min1);
    yo2max1 = fmax(yo2s4, yo2max1);

    // Get value of Toxyd from hoxyd

    cs_real_t hoxyd = enthox[c_id];

    cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    for (int ige = 0; ige < ngazem; ige++)
      coefe[ige] = 0;

    const cs_real_t _f3m = f3m[c_id];
    const cs_real_t _f4m = f4m[c_id];
    const cs_real_t _f5m = f5m[c_id];
    const cs_real_t _f6m = f6m[c_id];
    const cs_real_t _f7m = f7m[c_id];
    const cs_real_t _f8m = f8m[c_id];
    const cs_real_t _f9m = f9m[c_id];

    int el[] = {io2, in2, ico2, ih2o, ico, ihy};

    for (int el_id = 0; el_id < 6; el_id++) {
      int e = el[el_id];
      coefe[e] = (  cm->af3[e]*_f3m + cm->af4[e]*_f4m
                  + cm->af5[e]*_f5m + cm->af6[e]*_f6m
                  + cm->af7[e]*_f7m + cm->af8[e]*_f8m
                  + cm->af9[e]*_f9m) * wmole[e];
    }

    // Nb moles which reacted with CO to generate CO2

    cs_real_t deltamol   = (coefe[io2]-yo2ox) / wmole[io2];
    cs_real_t react      = fmin(2.0*deltamol, coefe[ico]/wmole[ico]);
    coefe[io2] = coefe[io2] - react/2.0*wmole[io2];
    coefe[ico] = coefe[ico]-react     *wmole[ico];
    coefe[ico2] = coefe[ico2]+react   *wmole[ico2];

    cs_real_t sum =  coefe[io2] +coefe[in2]+coefe[ico2]
                    +coefe[ih2o]+coefe[ico]+coefe[ihy];
    coefe[io2] = coefe[io2] / sum;
    coefe[in2] = coefe[in2] /sum;
    coefe[ico2] = coefe[ico2] /sum;
    coefe[ih2o] = coefe[ih2o] /sum;
    coefe[ico] = coefe[ico] /sum;
    coefe[ihy] = coefe[ihy] /sum;

    cs_real_t toxyd = cs_coal_ht_convert_h_to_t_gas_by_yi(hoxyd, coefe);

    toxmin = fmin(toxmin, toxyd);
    toxmax = fmax(toxmax, toxyd);

    if (toxyd > cpro_temp[c_id])
      toxyd = cpro_temp[c_id];

    // Initialize with temperatures Toxy et Tfuel at bounds.

    cs_real_t dirac = dfuel[c_id]*_tfuel + doxyd[c_id]*toxyd;

    // Get value of mean temperature.

    cs_real_t tmpgaz = cpro_temp[c_id];

    bb1 = fmax(0., pdfm1[c_id]);
    bb2 = fmin(fs4no[c_id], pdfm2[c_id]);
    cs_real_t bb3 = fmax(fs4no[c_id], pdfm1[c_id]);
    cs_real_t bb4 = fmin(1., pdfm2[c_id]);

    if (bb2 > bb1)
      lro = hrec[c_id];
    else
      lro = 0.;
    cs_real_t lrf = (bb4 > bb3) ? hrec[c_id] : 0.;

    qqq = bb2*bb2 - bb1*bb1;
    rrr = bb2     - bb1;
    cs_real_t sss = bb4*bb4 - bb3*bb3;
    cs_real_t ttt = bb4    - bb3;
    cs_real_t uuu = 1.0    - fs4no[c_id];

    gt1 = lro*qqq / (2.0*fs4no[c_id]);
    gt2 = lro*rrr;

    cs_real_t gt10 = lrf*sss / (2.0*uuu);
    cs_real_t gt20 = lrf*ttt;
    cs_real_t gt30 = lrf*ttt / uuu;

    cs_real_t ts4num = tmpgaz - dirac + toxyd*(gt1 -gt2)
                              - _tfuel*(gt10+gt20-gt30);
    cs_real_t ts4den = gt1-gt10+gt30;

    cs_real_t ts4 = ts4num / ts4den;

    if (xxf[c_idx] > cs_coal_epsilon) {

      cs_real_t hfuel = hhf[c_idx] / xxf[c_idx];

      cs_real_t hfs4ad = fs4no[c_id]*hfuel + (1.-fs4no[c_id])*hoxyd;

      for (int ige = 0; ige < ngazem; ige++)
        coefe[ige] = 0;

      // TODO check if we should not have an additional indirection below
      // (not needed if gas species come first).

      for (int ige = 0; ige < n_gas_sp; ige++)
        coefe[ige] = yfs4no[ige*n_cells + c_id];

      cs_real_t tfs4ad = cs_coal_ht_convert_h_to_t_gas_by_yi(hfs4ad, coefe);

      // Compute for logging

      ts4min = fmin(ts4min, ts4);
      ts4max = fmax(ts4max, ts4);

      ts4admin = fmin(ts4admin, tfs4ad);
      ts4admax = fmax(ts4admax, tfs4ad);

      cs_real_t summ = 0.;
      for (int ige = 0; ige < n_gas_sp; ige++)
        summ += yfs4no[ige*n_cells + c_id];

      summin = fmin(summin, summ);
      summax = fmax(summax, summ);

      cs_real_t  min_togf = fmin(toxyd, fmin(tmpgaz, _tfuel));

      if ((ts4 > min_togf) && (ts4 < 2400.0))
        inok++;

      if (ts4 < min_togf) {
        if (ts4 >= 300.0)
          i300++;
        else if (ts4 > 0)
          i000++;
        else
          imini++;
      }

      if (ts4 > 2400.) {
        if (ts4 < 2500)
          i2500++;
        else if (ts4 < 2600)
          i2600++;
        else if (ts4 < 2700)
          i2700++;
        else if (ts4 < 2800)
          i2800++;
        else if (ts4 < 3000)
          i3000++;
        else if (ts4 < 3500)
          i3500++;
        else if (ts4 < 4000)
          i4000++;
        else if (ts4 < 5000)
          i5000++;
        else
          imaxi++;
      }

      // End of computation for logging

      // Clipping of Ts4: to min(toxyd, tmpgaz, tfuel) on min
      //                  to ts4ad                     on max

      nbpt++;
      if (ts4 < min_togf) {
        nbclip1++;
        ts4 = min_togf;
      }

      if (ts4 > tfs4ad) {
        nbclip2++;
        ts4 = tfs4ad;
      }

      // Oxygen concentration

      cs_real_t xo2 = cpro_yo2[c_id] * cpro_mmel[c_id] / wmole[io2];

      // Integration

      cs_real_t val[CS_COAL_NOXST_NPART+1];
      cs_real_t tt[CS_COAL_NOXST_NPART+1];
      cs_real_t gs[CS_COAL_NOXST_NPART+1];
      cs_real_t yyo2[CS_COAL_NOXST_NPART+1];

      cs_real_t f_npart = CS_COAL_NOXST_NPART;

      for (int i = 0; i < CS_COAL_NOXST_NPART+1; i++) {
        gs[i] = pdfm1[c_id] +   (cs_real_t)i / f_npart
                              * (pdfm2[c_id]-pdfm1[c_id]);
        // compute T
        if (gs[i] < fs4no[c_id])
          tt[i] = (ts4-toxyd) / fs4no[c_id] * gs[i] + toxyd;
        else
          tt[i] =   (_tfuel-ts4) / (1.0-fs4no[c_id]) * gs[i]
                  + _tfuel - (_tfuel-ts4) / (1.0-fs4no[c_id]);

        // compute yo2
        if (gs[i] < fs4no[c_id])
          yyo2[i] = (yo2s4-yo2ox) / fs4no[c_id] * gs[i] + yo2ox;
        else if (gs[i] < fs3no[c_id]) {
          cs_real_t aa = yo2s4 / (fs4no[c_id]-fs3no[c_id]);
          yyo2[i] = aa * (gs[i] -fs3no[c_id]);
        }
        else
          yyo2[i] = 0.;
      }

      // No integration

      cs_real_t dgs = (pdfm2[c_id] - pdfm1[c_id]) / f_npart;

      cs_real_t bb;

      // Compute K1*EXP(-E1/T)

      if (ipdf1 == 1) {

        cpro_exp1[c_id] =   kk1 * exp(-ee1/toxyd) * doxyd[c_id]
                          + kk1 * exp(-ee1/_tfuel) * dfuel[c_id];

        for (int i = 0; i < CS_COAL_NOXST_NPART+1; i++)
          val[i] = kk1 * exp(-ee1/tt[i]) * hrec[c_id];

        for (int i = 0; i < CS_COAL_NOXST_NPART; i++)
          cpro_exp1[c_id] += 0.5 * dgs * (val[i] + val[i+1]);

      }

      // Compute K2*EXP(-E2/T)

      if (ipdf2 == 1) {

        if (xo2 > 0.) {

          if (xo2 > 0.018)
            bb = 0.;
          else if (xo2 < 0.0025)
            bb = 1.;
          else
            bb = (0.018-xo2) / (0.018-0.0025);

          cs_real_t xo2_p_bb = pow(xo2, bb);

          cpro_exp2[c_id] =   kk2
                            * (  exp(-ee2/toxyd)*doxyd[c_id]
                               + exp(-ee2/_tfuel)*dfuel[c_id]) * xo2_p_bb;

          for (int i = 0; i < CS_COAL_NOXST_NPART+1; i++)
            val[i] = kk2 * exp(-ee2/tt[i]) * hrec[c_id];

          for (int i = 0; i < CS_COAL_NOXST_NPART; i++)
            cpro_exp2[c_id] += 0.5 * dgs * (val[i]+val[i+1]) * xo2_p_bb;

        }
        else
          cpro_exp2[c_id] = 0;

      }

      // Compute K3*EXP(-E3/T)

      if (ipdf3 == 1) {

        if (xo2 > 0.) {
          cpro_exp3[c_id] =   kk3
                            * (  exp(-ee3/toxyd)*doxyd[c_id]*sqrt(yo2ox)
                               + exp(-ee3/_tfuel)*dfuel[c_id]*sqrt(yo2cb));

          for (int i = 0; i < CS_COAL_NOXST_NPART+1; i++) {
            if (yyo2[i] > 0) {
              if (gs[i] <= fs3no[c_id])
                val[i] = kk3*exp(-ee3/tt[i]) * hrec[c_id] * sqrt(yyo2[i]);
              else
                val[i] = 0.;
            }
            else
              val[i] = 0.;
          }

          for (int i = 0; i < CS_COAL_NOXST_NPART; i++)
            cpro_exp3[c_id] += 0.5 * dgs * (val[i] + val[i+1]);
        }
        else
          cpro_exp3[c_id] = 0.;

      }

      // Compute K4*EXP(-E4/T)

      if (ipdf4 == 1) {

        if (xo2 > 0.) {

          if (xo2 > 0.018)
            bb = 0.;
          else if (xo2 < 0.0025)
            bb = 1.;
          else
            bb = (0.018-xo2) / (0.018-0.0025);

          cs_real_t xo2_p_bb = pow(xo2, bb);

          cpro_exp4[c_id] =   kk4
                            * (  exp(-ee4/toxyd) * doxyd[c_id]
                               + exp(-ee4/_tfuel) * dfuel[c_id]) * xo2_p_bb;

          for (int i = 0; i < CS_COAL_NOXST_NPART+1; i++)
            val[i] = kk4 * exp(-ee4/tt[i]) * hrec[c_id];

          for (int i = 0; i < CS_COAL_NOXST_NPART; i++)
            cpro_exp4[c_id] += 0.5 * dgs * (val[i] + val[i+1]) * xo2_p_bb;

        }
        else
          cpro_exp4[c_id] = 0.;

      }

      // Compute K5*EXP(-E5/T)

      if (ipdf5 == 1) {

        cpro_exp5[c_id] =  kk5
                          * (  exp(-ee5/toxyd)*doxyd[c_id]
                             + exp(-ee5/_tfuel)*dfuel[c_id]);

        for (int i = 0; i < CS_COAL_NOXST_NPART+1; i++)
          val[i] = kk5 * exp(-ee5/tt[i]) * hrec[c_id];

        for (int i = 0; i < CS_COAL_NOXST_NPART; i++)
          cpro_exp5[c_id] += 0.5 * dgs * (val[i] + val[i+1]);

      }

    } // if (xxf[c_idx] > cs_coal_epsilon)

  } // loop on cells

  CS_FREE(xxf);
  CS_FREE(hhf);

  CS_FREE(tfuel);
  CS_FREE(compute_id);

  if (cs_log_default_is_active()) {

    cs_gnum_t cpt[] = {
      inok,  imini, i000, i300, i2500, i2600,
      i2700, i2800, i3000, i3500, i4000, i5000,
      imaxi, nbpt, nbclip1, nbclip2, nbclip30, nbclip31
    };

    cs_parall_counter(cpt, 18);

    cs_real_t minmax[] = {
      ts4min, -ts4max, ts4admin, -ts4admax, summin, -summax,
      yo2min, -yo2max, yo2min1, -yo2max1,
      toxmin, -toxmax, yo2oxmin, -yo2oxmax
    };

    cs_parall_min(14, CS_REAL_TYPE, minmax);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n"
         " Min max of TSox  %g %g\n"
         " Min max of TS4  %g %g\n"
         "    number of points at Tmini < ts4 < 2400   %lu\n"
         "    number of points at         ts4 < 0      %lu\n"
         "    number of points at 0     < ts4 < 300    %lu\n"
         "    number of points at 300   < ts4 < Tmini  %lu\n"
         "    number of points at 2400  < ts4 < 2500   %lu\n"
         "    number of points at 2500  < ts4 < 2600   %lu\n"
         "    number of points at 2600  < ts4 < 2700   %lu\n"
         "    number of points at 2700  < ts4 < 2800   %lu\n"
         "    number of points at 2800  < ts4 < 3000   %lu\n"
         "    number of points at 3000  < ts4 < 3500   %lu\n"
         "    number of points at 3500  < ts4 < 4000   %lu\n"
         "    number of points at 4000  < ts4 < 5000   %lu\n"
         "    number of points at 5000  < ts4          %lu\n"),
       minmax[10], -minmax[11], minmax[0], -minmax[1],
       (unsigned long)cpt[0], (unsigned long)cpt[1],
       (unsigned long)cpt[2], (unsigned long)cpt[3],
       (unsigned long)cpt[4], (unsigned long)cpt[5],
       (unsigned long)cpt[6], (unsigned long)cpt[7],
       (unsigned long)cpt[8], (unsigned long)cpt[9],
       (unsigned long)cpt[10], (unsigned long)cpt[11],
       (unsigned long)cpt[12]);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _(" Min max of TS4ad  %g %g\n"
         " Min max concentration in fs4  %g %g\n"
         " Clip to min at %lu of %lu points\n"
         " Clip to max at %lu of %lu points\n"),
       minmax[2], -minmax[3], minmax[4], -minmax[5],
       (unsigned long)cpt[14], (unsigned long)cpt[13],
       (unsigned long)cpt[15], (unsigned long)cpt[13]);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n"
         " Min max of Yo2ox at 0  %g %g\n"
         " Min max of Yo2 in fs4 before clipping  %g %g\n"
         " Clip to min over Yo2 in fs4            %lu\n"
         " Clip to max over Yo2 in fs4            %lu\n"
         " Min max of Yo2 in fs4 after clipping   %g %g\n"),
       minmax[12], -minmax[13], minmax[6], -minmax[7],
       (unsigned long)cpt[16], (unsigned long)cpt[17],
       minmax[8], -minmax[9]);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
