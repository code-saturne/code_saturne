/*============================================================================
 * Coal combustion model: mass transfer between continous and dispersed phase
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
#include "comb/cs_coal.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_thermal_model.h"

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
  \file cs_coal_masstransfer.cpp

  \brief Coal combustion model: Calculation of the terms of mass transfer
         between the continous phase and the dispersed phase.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute mass transfer terms between the continous and dispersed phase.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_mass_transfer(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t pi = cs_math_pi;
  const cs_real_t r = cs_physical_constants_r;
  const cs_real_t prefth = cs_coal_prefth;

  cs_coal_model_t *cm = cs_glob_coal_model;

  /* Initialization
   * -------------- */

  const int io2 = cm->io2 -1;
  const int ico2 = cm->ico2 -1;
  const int ih2o = cm->ih2o -1;

  // Molar masses
  double wmo2  = cm->wmole[io2];
  double wmco2 = cm->wmole[ico2];
  double wmh2o = cm->wmole[ih2o];

  /* Aliases for simpler syntax */

  const cs_real_t eps_cp = cs_coal_epsilon;

  cs_host_context ctx;

  /* Initialization and preliminary computations
     ------------------------------------------- */

  const cs_real_t *cpro_temp = CS_F_(t)->val;
  const cs_real_t *cpro_rom1 = cs_field_by_id(cm->irom1)->val;

  cs_real_t *x2, *x2srho2, *rho1, *w1;
  CS_MALLOC(x2, n_cells, cs_real_t);
  CS_MALLOC(x2srho2, n_cells, cs_real_t);
  CS_MALLOC(rho1, n_cells, cs_real_t);
  CS_MALLOC(w1, n_cells, cs_real_t);

  // Initialization of mass transfer terms

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
    cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;
    cs_real_t *cpro_cgch = cs_field_by_id(cm->igmdch[class_id])->val;
    cs_real_t *cpro_cght = cs_field_by_id(cm->igmhet[class_id])->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cpro_cgd1[c_id] = 0.;
      cpro_cgd2[c_id] = 0.;
      cpro_cgch[c_id] = 0.;
      cpro_cght[c_id] = 0.;
    });
  }

  // Calculation of mass density of the gas mixture

  const cs_real_t *crom = CS_F_(rho)->val;

  // Calculation of x2=sum(X2i) and x2sro2 = sum(X2i/Rho2i)

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    x2[c_id] = 0.;
    x2srho2[c_id] = 0.;
  });

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const cs_real_t xmash = cm->xmash[class_id];
    const cs_real_t *cpro_rom2 = cs_field_by_id(cm->irom2[class_id])->val;
    const cs_real_t *cvara_xchcl = cs_field_by_id(cm->ixch[class_id])->val_pre;
    const cs_real_t *cvara_xckcl = cs_field_by_id(cm->ixck[class_id])->val_pre;
    const cs_real_t *cvara_xnpcl = cs_field_by_id(cm->inp[class_id])->val_pre;
    const cs_real_t *cvara_xwtcl = nullptr;
    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
      cvara_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val_pre;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t xch  = cvara_xchcl[c_id];
      cs_real_t xck  = cvara_xckcl[c_id];
      cs_real_t xash = cvara_xnpcl[c_id]*xmash;
      cs_real_t xx2  = xch + xck + xash;
      if (cvara_xwtcl != nullptr)  // humidity
        xx2 += cvara_xwtcl[c_id];

      x2[c_id] += xx2;
      x2srho2[c_id] += xx2 / cpro_rom2[c_id];
    });

  }

  // Rho 1

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    rho1[c_id] = (1.-x2[c_id]) / (1./crom[c_id]-x2srho2[c_id]);
  });

  /* Mass transfer by devolatilization
     --------------------------------- */

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const int coal_id = cm->ichcor[class_id] - 1;

    cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
    cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;
    cs_real_t *cpro_cgch = cs_field_by_id(cm->igmdch[class_id])->val;
    const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;

    const cs_real_t *y1ch = cm->y1ch;
    const cs_real_t *y2ch = cm->y2ch;
    const cs_real_t *a1ch = cm->a1ch;
    const cs_real_t *a2ch = cm->a2ch;
    const cs_real_t *e1ch = cm->e1ch;
    const cs_real_t *e2ch = cm->e2ch;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      // Mass transfer due to light mass density (s-1) < 0
      cpro_cgd1[c_id] =   -y1ch[coal_id]*a1ch[coal_id]
                        * exp(-e1ch[coal_id]/(r*cpro_temp2[c_id]));

      // Mass transfer due to heavy mass density (s-1) < 0
      cpro_cgd2[c_id] =   -y2ch[coal_id]*a2ch[coal_id]
                        * exp(-e2ch[coal_id]/(r*cpro_temp2[c_id]));

      // Rate of disappearance of reactive coal (s-1) < 0

      cpro_cgch[c_id] = -a1ch[coal_id] * exp(-e1ch[coal_id]/(r*cpro_temp2[c_id]))
                       - a2ch[coal_id] * exp(-e2ch[coal_id]/(r*cpro_temp2[c_id]));
    });

  }

  /* Calculation of average RHO_COKE for each coal
     ---------------------------------------------
   * It is assumed for the calculation of the mass density of the coke
   * that the devolatization takes place at constant volume */

  cs_real_t devto1[CS_COMBUSTION_MAX_COALS];
  cs_real_t devto2[CS_COMBUSTION_MAX_COALS];

  // Initialization
  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
    devto1[coal_id] = 0;
    devto2[coal_id] = 0;
    cm->rhock[coal_id] = cm->rho0ch[coal_id];
  }

  // Calculation of the integral of GMDEV1 and GMDEV2 for each coal
  // TODO: use more generic and accurate reduction instead of naive one here.

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const int coal_id = cm->ichcor[class_id] - 1;

    const cs_real_t *cvara_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
    const cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;

    # pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (auto c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t xch = cvara_xchcl[c_id];
      devto1[coal_id] -= cpro_cgd1[c_id]*xch*crom[c_id]*volume[c_id];
      devto2[coal_id] -= cpro_cgd2[c_id]*xch*crom[c_id]*volume[c_id];
    }

  }

  cs_parall_sum(cm->n_coals, CS_REAL_TYPE, devto1);
  cs_parall_sum(cm->n_coals, CS_REAL_TYPE, devto2);

  // Calculation of average mass density of coke

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
    cs_real_t den =   cm->y2ch[coal_id]*devto1[coal_id]
                    + cm->y1ch[coal_id]*devto2[coal_id];
    if (den > eps_cp) {
      cm->rhock[coal_id] =   cm->rho0ch[coal_id]
                           * (1.0 - (   cm->y1ch[coal_id]*cm->y2ch[coal_id]
                                     * (devto1[coal_id]+devto2[coal_id]) / den));
    }
  }

  /* Mass transfer by heterogeneous combustion with O2
   * ------------------------------------------------- */

  const cs_real_t *cpro_yox = cs_field_by_id(cm->iym1[io2])->val;

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {

    const cs_real_t *cvara_xnpcl = cs_field_by_id(cm->inp[class_id])->val_pre;
    const cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
    const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;
    cs_real_t *cpro_cght = cs_field_by_id(cm->igmhet[class_id])->val;

    const int coal_id = cm->ichcor[class_id] - 1;

    const cs_real_t diam20 = cm->diam20[class_id];
    const cs_real_t rho20 = cm->rho20[class_id];
    const cs_real_t xmp0  = cm->xmp0[class_id];

    const int iochet = cm->iochet[coal_id];
    const cs_real_t rhock = cm->rhock[coal_id];
    const cs_real_t xashch = cm->xashch[coal_id];
    const cs_real_t ahetch = cm->ahetch[coal_id];
    const cs_real_t ehetch = cm->ehetch[coal_id];

    const cs_real_t diam20_2 = diam20*diam20;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

      const cs_real_t xnp   = cvara_xnpcl[c_id];
      const cs_real_t xuash = xnp*(1.-xashch)*xmp0;

      // Compute partial pressure of oxygen [atm]

      // PO2 = Rho1*cs_physical_constants_r*T*YO2/MO2

      cs_real_t pparo2 = rho1[c_id] * r * cpro_temp[c_id]
                                        * cpro_yox[c_id] / wmo2;
      pparo2 /= prefth;

      // Coefficient of chemical kinetics of formation of CO
      // in [kg.m-2.s-1.atm(-n)]

      cs_real_t xdfchi =   ahetch* exp(-ehetch*4185.
                         / (r * cpro_temp2[c_id]));

      // Diffusion coefficient  in kg/m2/s/[atm] : XDFEXT
      // Global coefficient for n=0.5 in kg/m2/s : XDFTOT0
      // Global coefficient for n=1   in Kg/m2/s : XDFTOT1

      cs_real_t xdftot0, xdftot1;

      cs_real_t diacka = cpro_diam2[c_id]/diam20;
      if (diacka > eps_cp) {
        cs_real_t xdfext =   2.53e-7*(pow(cpro_temp2[c_id], 0.750))
                           / cpro_diam2[c_id] * 2.;
        cs_real_t xdfext_2 = xdfext*xdfext;
        cs_real_t xdfchi_2 = xdfchi*xdfchi;
        cs_real_t xdfchi_4 = xdfchi_2*xdfchi_2;
        xdftot1 = pparo2 / (1./xdfchi + 1./xdfext);
        xdftot0 = - xdfchi_2 / (2. * xdfext)
                  + sqrt(pparo2*xdfchi_2 + xdfchi_4/(4.*xdfext_2));
      }
      else {
        xdftot1 = xdfchi * pparo2;
        xdftot0 = xdfchi * sqrt(pparo2);
      }

      // Remark AE: For now, we limit ourselves to this test
      //         The introduction of a blowing correlation will come later

      // Calculation of coxck such as: Se = coxck * Xck**(2/3)

      constexpr cs_real_t c_2ov3 = 2./3.;

      cs_real_t coxck = 0.;
      if (xuash > eps_cp) {
        coxck = pi * diam20_2 * pow(rho20/(rhock*xuash), c_2ov3);
      }

      // Calculation of cpro_cght[c_id] = - coxck*Xdftoto*PPARO2*Xnp < 0
      // or  cpro_cght[c_id] = - coxck*XDFTOT1*PPARO2*Xnp < 0

      if (iochet == 1)
        cpro_cght[c_id] = - xdftot1*coxck*xnp;
      else
        cpro_cght[c_id] = - xdftot0*coxck*xnp;
    });

  }

  /* Mass transfer by heterogeneous combustion with CO2
   * -------------------------------------------------- */

  if (cm->ihtco2 == 1) {

    const cs_real_t *cpro_yco2 = cs_field_by_id(cm->iym1[ico2])->val;

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {

      const cs_real_t *cvara_xnpcl = cs_field_by_id(cm->inp[class_id])->val_pre;
      const cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
      const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;
      cs_real_t *cpro_cght = cs_field_by_id(cm->ighco2[class_id])->val;

      const int coal_id = cm->ichcor[class_id] - 1;

      const cs_real_t diam20 = cm->diam20[class_id];
      const cs_real_t rho20 = cm->rho20[class_id];
      const cs_real_t xmp0  = cm->xmp0[class_id];

      const int iochet = cm->iochet[coal_id];
      const cs_real_t rhock = cm->rhock[coal_id];
      const cs_real_t xashch = cm->xashch[coal_id];
      const cs_real_t ahetc2 = cm->ahetc2[coal_id];
      const cs_real_t ehetc2 = cm->ehetc2[coal_id];

      const cs_real_t diam20_2 = diam20*diam20;

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

        const cs_real_t xnp   = cvara_xnpcl[c_id];
        const cs_real_t xuash = xnp*(1.-xashch)*xmp0;

        // Compute partial pressure of CO2 [atm]

        // PCO2 = Rho1*cs_physical_constants_r*T*YCO2/MCO2

        cs_real_t pprco2 = rho1[c_id] * r * cpro_temp[c_id]
                                          * cpro_yco2[c_id] / wmco2;
        pprco2 /= prefth;

        // Coefficient of chemical kinetics of formation of CO
        // in [kg.m-2.s-1.atm(-n)]

        cs_real_t xdfchi =   ahetc2* exp(-ehetc2*4185.
                           / (r * cpro_temp2[c_id]));

        // Diffusion coefficient  in kg/m2/s/[atm] : XDFEXT
        // Global coefficient for n=0.5 in kg/m2/s : XDFTOT0
        // Global coefficient for n=1   in Kg/m2/s : XDFTOT1

        cs_real_t xdftot0, xdftot1;

        cs_real_t diacka = cpro_diam2[c_id]/diam20;
        if (diacka > eps_cp) {
          cs_real_t xdfext =   2.53e-7*(pow(cpro_temp2[c_id], 0.750))
                             / cpro_diam2[c_id] * 2.;
          cs_real_t xdfext_2 = xdfext*xdfext;
          cs_real_t xdfchi_2 = xdfchi*xdfchi;
          cs_real_t xdfchi_4 = xdfchi_2*xdfchi_2;
          xdftot1 = pprco2 / (1./xdfchi + 1./xdfext);
          xdftot0 = - xdfchi_2 / (2. * xdfext)
                    + sqrt(pprco2*xdfchi_2 + xdfchi_4/(4.*xdfext_2));
        }
        else {
          xdftot1 = xdfchi * pprco2;
          xdftot0 = xdfchi * sqrt(pprco2);
        }

        // Remark AE: For now, we limit ourselves to this test
        //         The introduction of a blowing correlation will come later

        // Calculation of coxck such as: Se = coxck * Xck**(2/3)

        constexpr cs_real_t c_2ov3 = 2./3.;

        cs_real_t coxck = 0.;
        if (xuash > eps_cp) {
          coxck = pi * diam20_2 * pow(rho20/(rhock*xuash), c_2ov3);
        }

        // Calculation of cpro_cght[c_id] = - coxck*Xdftoto*PPRCO2*Xnp < 0
        // or  cpro_cght[c_id] = - coxck*XDFTOT1*PPRCO2*Xnp < 0

        if (iochet == 1)
          cpro_cght[c_id] = - xdftot1*coxck*xnp;
        else
          cpro_cght[c_id] = - xdftot0*coxck*xnp;
      });

    }

  }

  /* Mass transfer by heterogeneous combustion with H20
   * -------------------------------------------------- */

  if (cm->ihth2o == 1) {

    const cs_real_t *cpro_yh2o = cs_field_by_id(cm->iym1[ih2o])->val;

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {

      const cs_real_t *cvara_xnpcl = cs_field_by_id(cm->inp[class_id])->val_pre;
      const cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
      const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;
      cs_real_t *cpro_cght = cs_field_by_id(cm->ighh2o[class_id])->val;

      const int coal_id = cm->ichcor[class_id] - 1;

      const cs_real_t diam20 = cm->diam20[class_id];
      const cs_real_t rho20 = cm->rho20[class_id];
      const cs_real_t xmp0  = cm->xmp0[class_id];

      const int iochet = cm->iochet[coal_id];
      const cs_real_t rhock = cm->rhock[coal_id];
      const cs_real_t xashch = cm->xashch[coal_id];
      const cs_real_t ahetc2 = cm->ahetc2[coal_id];
      const cs_real_t ehetc2 = cm->ehetc2[coal_id];

      const cs_real_t diam20_2 = diam20*diam20;

      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

        const cs_real_t xnp   = cvara_xnpcl[c_id];
        const cs_real_t xuash = xnp*(1.-xashch)*xmp0;

        // Compute partial pressure of H2O [atm]

        // PH2O = Rho1*cs_physical_constants_r*T*YH2O/MH2O

        cs_real_t pprh2o = rho1[c_id] * r * cpro_temp[c_id]
                                          * cpro_yh2o[c_id] / wmh2o;
        pprh2o /= prefth;

        // Coefficient of chemical kinetics of formation of CO
        // in [kg.m-2.s-1.atm(-n)]

        cs_real_t xdfchi =   ahetc2* exp(-ehetc2*4185.
                           / (r * cpro_temp2[c_id]));

        // Diffusion coefficient  in kg/m2/s/[atm] : XDFEXT
        // Global coefficient for n=0.5 in kg/m2/s : XDFTOT0
        // Global coefficient for n=1   in Kg/m2/s : XDFTOT1

        cs_real_t xdftot0, xdftot1;

        cs_real_t diacka = cpro_diam2[c_id]/diam20;
        if (diacka > eps_cp) {
          cs_real_t xdfext =   2.53e-7*(pow(cpro_temp2[c_id], 0.750))
                             / cpro_diam2[c_id] * 2.;
          cs_real_t xdfext_2 = xdfext*xdfext;
          cs_real_t xdfchi_2 = xdfchi*xdfchi;
          cs_real_t xdfchi_4 = xdfchi_2*xdfchi_2;
          xdftot1 = pprh2o / (1./xdfchi + 1./xdfext);
          xdftot0 = - xdfchi_2 / (2. * xdfext)
                    + sqrt(pprh2o*xdfchi_2 + xdfchi_4/(4.*xdfext_2));
        }
        else {
          xdftot1 = xdfchi * pprh2o;
          xdftot0 = xdfchi * sqrt(pprh2o);
        }

        // Remark AE: For now, we limit ourselves to this test
        //         The introduction of a blowing correlation will come later

        // Calculation of coxck such as: Se = coxck * Xck**(2/3)

        constexpr cs_real_t c_2ov3 = 2./3.;

        cs_real_t coxck = 0.;
        if (xuash > eps_cp) {
          coxck = pi * diam20_2 * pow(rho20/(rhock*xuash), c_2ov3);
        }

        // Calculation of cpro_cght[c_id] = - coxck*Xdftoto*PPRCO2*Xnp < 0
        // or  cpro_cght[c_id] = - coxck*XDFTOT1*PPRCO2*Xnp < 0

        if (iochet == 1)
          cpro_cght[c_id] = - xdftot1*coxck*xnp;
        else
          cpro_cght[c_id] = - xdftot0*coxck*xnp;
      });

    }

  }

  /* Mass transfer during the dryer phase
     ------------------------------------ */

  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {

    // Latent heat in J/kg
    cs_real_t lv = 2.263e+6;
    cs_real_t tebl = 100. + cs_physical_constants_celsius_to_kelvin;

    cs_real_t shrd  = 2.;
    cs_real_t xmeau = 0.018;

    cs_real_t tlimit = 302.24;
    cs_real_t tmini = tlimit * (1.-tlimit/(lv*xmeau));
    const cs_real_t *cpro_viscls = nullptr, *cpro_cp = nullptr;

    const int kivisl = cs_field_key_id("diffusivity_id");
    const int kvisl0 = cs_field_key_id("diffusivity_ref");
    const cs_field_t *fld_th = cs_thermal_model_field();
    cs_real_t cp0 = cs_glob_fluid_properties->cp0;
    double visls_0 = cs_field_get_key_double(fld_th, kvisl0);
    int ifcvsl = cs_field_get_key_int(fld_th, kivisl);
    if (ifcvsl >= 0)
      cpro_viscls = cs_field_by_id(ifcvsl)->val;
    int icp = cs_glob_fluid_properties->icp;
    if (icp >= 0)
      cpro_cp = cs_field_by_id(icp)->val;

    if (ifcvsl >= 0) {
      if (icp >= 0) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
          w1[c_id] = cpro_viscls[c_id] * cpro_cp[c_id];
        });
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
          w1[c_id] = cpro_viscls[c_id] * cp0;
        });
      }
    }
    else {
      if (icp >= 0) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
          w1[c_id] = visls_0 * cpro_cp[c_id];
        });
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
          w1[c_id] = visls_0 * cp0;
        });
      }
    }

    /* Logging
       ------- */

    const cs_real_t *cpro_mmel = cs_field_by_id(cm->immel)->val;

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {

      const int coal_id = cm->ichcor[class_id] - 1;

      const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;
      const cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
      const cs_real_t *cpro_yh2o = cs_field_by_id(cm->iym1[ih2o])->val;

      const cs_real_t *cvara_xnpcl
        = cs_field_by_id(cm->inp[class_id])->val_pre;
      const cs_real_t *cvara_xwtcl
        = cs_field_by_id(cm->ixwt[class_id])->val_pre;

      cs_real_t *cpro_csec = cs_field_by_id(cm->igmsec[class_id])->val;

      const cs_real_t diam20 = cm->diam20[class_id];

      const cs_real_t xashch = cm->xashch[coal_id];

      const cs_real_t diam20_2 = diam20*diam20;

      cs_gnum_t npoin1 = 0, npoin2 = 0, npoin3 = 0, npoin4 = 0, npoin63 = 0;
      cs_gnum_t npyv = 0;
      cs_real_t yymax = 0;

      // Calculation of the diameter of particles in W2
      //   d20 = sqrt(A0^2 + (1-A0)*DCK^2)

      cs_real_t tmax = -HUGE_VAL;
      cs_real_t tmin = HUGE_VAL;
      cs_real_t yvmax = -HUGE_VAL;
      cs_real_t yvmin = HUGE_VAL;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        cpro_csec[c_id] = 0;

        if (cvara_xwtcl[c_id] > eps_cp) {

          npoin1++;

          cs_real_t xnp = cvara_xnpcl[c_id];

          // Calculation of the diameter of the particules in W2
          // d20 = sqrt((A0.^2+(1-A0)*DCK^2)

          cs_real_t dp
            = sqrt(      xashch*diam20_2
                   + (1.-xashch)*cs_math_pow2(cpro_diam2[c_id]));

          npoin2++;

          cs_real_t xmgaz = cpro_mmel[c_id];

          cs_real_t yvs;
          if (cpro_temp2[c_id] > tlimit) {
            yvs   =   xmeau/xmgaz
                    * exp(lv*xmeau * (1./tebl-1./cpro_temp2[c_id])/r);
          }
          else {
            yvs   =   xmeau/xmgaz
                    * exp(lv*xmeau *(1./tebl-1./tlimit)/r)
                    * (lv*xmeau*(cpro_temp2[c_id]-tmini))
                    / (tlimit*tlimit);
          }

          cs_real_t yv = yvs;

          if (yv < 0.) {
            bft_error
              (__FILE__, __LINE__, 0,
               _("yv = %g, yvs = %g yh2o = %g\n"
                 "temp2 = %g tmini = %g"),
               yv, yvs, cpro_yh2o[c_id], cpro_temp2[c_id], tmini);
          }

          if (yv >= 1.) {
            yymax = cs::max(yymax, yv);
            npyv++;
            yv= 0.99;
          }

          yvmin = cs::min(yvmin, yv);
          yvmax = cs::max(yvmax, yv);

          if (yv > eps_cp && yv < 1.) {
            npoin3++;
            cpro_csec[c_id] =   pi*dp*cpro_rom1[c_id]*visls_0*shrd*xnp
                              * log((1.-cpro_yh2o[c_id])/(1.-yv));

            if (cpro_csec[c_id] < 0.) {
              cpro_csec[c_id] = 0.;
              npoin63++;
            }
          }
          else {
            npoin4++;
            cpro_csec[c_id] = 0.;
          }

        }
        else {
          cpro_csec[c_id] = 0.;
        }

        tmax = cs::max(cpro_csec[c_id], tmax);
        tmin = cs::min(cpro_csec[c_id], tmin);

      } /* Loop on cells */

      if (cs_log_default_is_active()) {

        cs_gnum_t npoint = (cs_gnum_t)n_cells;

        cs_parall_sum_scalars(npoint, npoin1, npoin2, npoin3, npoin4,
                              npoin63, npyv);
        cs_parall_min_scalars(tmin, yvmin);
        cs_parall_max_scalars(tmax, yvmax, yymax);

        cs_log_printf(CS_LOG_DEFAULT,
                      _("\n For the class %d, n_points = %lu\n"
                        "   Number of pts Xwat > 0 = %lu\n"
                        "                 T < Tebl = %lu\n"
                        "                 Yv   > 0 = %lu\n"
                        "                 Yv   < 0 = %lu\n"
                        "     Cancel G    T < Tebl = %lu\n"
                        "   Min Max ST = %g %g\n"
                        "   Min Max Yv = %g %g\n"
                        "   Clipping of YV at Max %lu, %g\n"),
                      class_id+1, (unsigned long)npoint,
                      (unsigned long)npoin1, (unsigned long)npoin2,
                      (unsigned long)npoin3, (unsigned long)npoin4,
                      (unsigned long)npoin63,
                      tmin, tmax, yvmin, yvmax,
                      (unsigned long)npyv, yymax);

      } /* Logging */

    } /* Loop on coal classes */

  } /* Case with drying */


  /* Free memory */

  CS_FREE(w1);
  CS_FREE(rho1);
  CS_FREE(x2srho2);
  CS_FREE(x2);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
