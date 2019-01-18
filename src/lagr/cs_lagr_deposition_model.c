/*============================================================================
 * Methods for particle deposition
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*============================================================================
 * Functions dealing with particle deposition
 *============================================================================*/

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"

#include "cs_physical_constants.h"
#include "cs_random.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_deposition_model.h"

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Boltzmann constant */
static const double _k_boltz = 1.38e-23;

/*=============================================================================
 * Prototypes for private functions used in a cross-recursive manner (ugly)
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Management of the diffusion in the internal zone (y^+ < dintrf)
 *
 * parameters:
 *   dx        -->    wall-normal displacement
 *   vvue      <->    wall-normal velocity of the flow seen
 *   vpart     <->    particle wall-normal velocity
 *   marko     <->    state of the jump process
 *   tempf     <--    temperature of the fluid
 *   dtl       -->    Lagrangian time step
 *   tstruc    <->    coherent structure mean duration
 *   tdiffu    <->    diffusion phase mean duration
 *   ttotal    <->    tdiffu + tstruc
 *   vstruc    <--    coherent structure velocity
 *   romp      <--    particle density
 *   taup      <->    particle relaxation time
 *   kdif      <->    diffusion phase diffusion coefficient
 *   tlag2     <->    diffusion relaxation timescale
 *   lvisq     <--    wall-unit lengthscale
 *   yplus     <--    wall-normal velocity of the flow seen
 *   unif1     <->    random number (uniform law)
 *   unif2     <->    random number (uniform law)
 *   dintrf    <--    extern-intern interface location
 *   rpart     <--    particle radius
 *   kdifcl    <->    internal zone diffusion coefficient
 *   indint    <->    interface indicator
 *   gnorm     <--    wall-normal gravity component
 *   vnorm     <--    wall-normal fluid (Eulerian) velocity
 *   grpn      <--    wall-normal pressure gradient
 *   piiln     <--    SDE integration auxiliary term
 *----------------------------------------------------------------------------*/

static void
_dep_inner_zone_diffusion(cs_real_t *dx,
                          cs_real_t *vvue,
                          cs_real_t *vpart,
                          cs_lnum_t *marko,
                          cs_real_t  tempf,
                          cs_real_t *depint,
                          cs_real_t  dtl,
                          cs_real_t *tstruc,
                          cs_real_t *tdiffu,
                          cs_real_t *ttotal,
                          cs_real_t *vstruc,
                          cs_real_t  romp,
                          cs_real_t  taup,
                          cs_real_t *kdif,
                          cs_real_t *tlag2,
                          cs_real_t *yplus,
                          cs_real_t  lvisq,
                          cs_real_t *unif1,
                          cs_real_t *unif2,
                          cs_real_t *dintrf,
                          cs_real_t *rpart,
                          cs_real_t *kdifcl,
                          cs_lnum_t *indint,
                          cs_real_t *gnorm,
                          cs_real_t *vnorm,
                          cs_real_t *grpn,
                          cs_real_t *piiln);

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Management of the ejection coherent structure (marko = CS_LAGR_COHERENCE_STRUCT_EJECTION)
 *
 * parameters:
 *   marko     -->    state of the jump process
 *   depint    <--    interface location near-wall/core-flow
 *   rpart     <--    particle radius
 *   kdifcl    <--    internal zone diffusion coefficient
 *   dtp       <--    Lagrangian timestep
 *   tstruc    <--    coherent structure mean duration
 *   vstruc    <--    coherent structure velocity
 *   lvisq     <--    wall-unit lengthscale
 *   dx        <->    wall-normal displacement
 *   vpart     <->    particle wall-normal velocity
 *   vvue      <->    wall-normal velocity of the flow seen
 *   taup      <--    particle relaxation time
 *   yplus     <--    particle wall-normal normalized distance
 *   unif1     <--    random number (uniform law)
 *   unif2     <--    random number (uniform law)
 *   dintrf    <--    extern-intern interface location
 *   gnorm     <--    wall-normal gravity component
 *   vnorm     <--    wall-normal fluid (Eulerian) velocity
 *----------------------------------------------------------------------------*/

static void
_dep_ejection(cs_lnum_t *marko,
              cs_real_t *depint,
              cs_real_t  dtp,
              cs_real_t *tstruc,
              cs_real_t *vstruc,
              cs_real_t  lvisq,
              cs_real_t *dx,
              cs_real_t *vvue,
              cs_real_t *vpart,
              cs_real_t  taup,
              cs_real_t *yplus,
              cs_real_t *unif1,
              cs_real_t *dintrf,
              cs_real_t *gnorm,
              cs_real_t *vnorm)
{
  cs_real_t  vpart0, vvue0, ypaux;

  vvue0  = *vvue;
  vpart0 = *vpart;

  /* Gravity and ormal fluid velocity added   */

  *vvue  =  -*vstruc + *gnorm * taup + *vnorm;
  *vpart = vpart0 * exp ( -dtp / taup) + (1 - exp ( -dtp / taup)) * vvue0;
  *dx    =   vvue0 * dtp
           + vvue0 * taup * (exp ( -dtp / taup) - 1.)
           + vpart0 * taup * (1. - exp ( -dtp / taup));
  ypaux  = *yplus - *dx / lvisq;

  /* --------------------------------------------------------
   *    Dissociation of cases by the arrival position
   * --------------------------------------------------------*/

  if (ypaux > *depint)
    *marko    =  -2;
  else if (ypaux < *dintrf)
    *marko    = CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF;
  else {

    if (*unif1 < (dtp / *tstruc))
      *marko  = CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION;
    else
      *marko  = CS_LAGR_COHERENCE_STRUCT_EJECTION;

  }
}

/*-----------------------------------------------------------------------------
 * Management of the sweep coherent structure (marko = CS_LAGR_COHERENCE_STRUCT_SWEEP)
 *
 * parameters:
 *   dx        <->    wall-normal displacement
 *   vvue      <->    wall-normal velocity of the flow seen
 *   vpart     <->    particle wall-normal velocity
 *   marko     -->    state of the jump process
 *   tempf     <--    temperature of the fluid
 *   dtp       <--    Lagrangian time step
 *   tstruc    <--    coherent structure mean duration
 *   tdiffu    <--    diffusion phase mean duration
 *   ttotal    <--    tdiffu + tstruc
 *   vstruc    <--    coherent structure velocity
 *   romp      <--    particle density
 *   taup      <--    particle relaxation time
 *   kdif      <--    diffusion phase diffusion coefficient
 *   tlag2     <--    diffusion relaxation timescale
 *   lvisq     <--    wall-unit lengthscale
 *   yplus     <--    wall-normal velocity of the flow seen
 *   unif1     <--    random number (uniform law)
 *   unif2     <--    random number (uniform law)
 *   dintrf    <--    extern-intern interface location
 *   rpart     <--    particle radius
 *   kdifcl    <--    internal zone diffusion coefficient
 *   gnorm     <--    wall-normal gravity component
 *   vnorm     <--    wall-normal fluid (Eulerian) velocity
 *   grpn      <--    wall-normal pressure gradient
 *   depint    <--    interface location near-wall/core-flow
 *   piiln     <--    SDE integration auxiliary term
 *----------------------------------------------------------------------------*/

static void
_dep_sweep(cs_real_t *dx,
           cs_real_t *vvue,
           cs_real_t *vpart,
           cs_lnum_t *marko,
           cs_real_t  tempf,
           cs_real_t *depint,
           cs_real_t  dtp,
           cs_real_t *tstruc,
           cs_real_t *tdiffu,
           cs_real_t *ttotal,
           cs_real_t *vstruc,
           cs_real_t  romp,
           cs_real_t  taup,
           cs_real_t *kdif,
           cs_real_t *tlag2,
           cs_real_t  lvisq,
           cs_real_t *yplus,
           cs_real_t *unif1,
           cs_real_t *unif2,
           cs_real_t *dintrf,
           cs_real_t *rpart,
           cs_real_t *kdifcl,
           cs_real_t *gnorm,
           cs_real_t *vnorm,
           cs_real_t *grpn,
           cs_real_t *piiln)
{
  /* -------------------------------------------------------
   *  Computation phase
   * ------------------------------------------------------- */

  cs_real_t vvue0  = *vvue;
  cs_real_t vpart0 = *vpart;
  *vvue  = *vstruc + *gnorm * taup + *vnorm;

  /*  Deposition submodel */

  *vpart = vpart0 * exp ( -dtp / taup) + (1. - exp ( -dtp / taup)) * vvue0;
  *dx    =   vvue0 * dtp + vvue0 * taup * (exp ( -dtp / taup) - 1.)
           + vpart0 * taup * (1. - exp ( -dtp / taup));

  /* Y+ at arrival */
  cs_real_t yplusa = *yplus - *dx / lvisq;

  /* --------------------------------------------------------
   *    Dissociation of cases by the arrival position
   * --------------------------------------------------------*/

  if (yplusa > *depint)
    *marko    =  -2;

  else if (yplusa < *dintrf) {

    cs_real_t dtp1 = (*dintrf - yplusa) * lvisq / CS_ABS(*vpart);
    *dx *= (*dintrf - *yplus) / (yplusa - *yplus);
    cs_real_t dxaux     = *dx;
    cs_real_t ypluss    = *yplus;
    *yplus    = *dintrf;
    *vvue     =  -*vstruc + *gnorm * taup + *vnorm;
    *marko    = CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF;
    cs_lnum_t indint    = 1;

    _dep_inner_zone_diffusion(dx,
                              vvue,
                              vpart,
                              marko,
                              tempf,
                              depint,
                              dtp1,
                              tstruc,
                              tdiffu,
                              ttotal,
                              vstruc,
                              romp,
                              taup,
                              kdif,
                              tlag2,
                              yplus,
                              lvisq,
                              unif1,
                              unif2,
                              dintrf,
                              rpart,
                              kdifcl,
                              &indint,
                              gnorm,
                              vnorm,
                              grpn,
                              piiln);

    indint  = 0;
    *dx    += dxaux;
    yplusa  = ypluss - *dx / lvisq;

    if (yplusa > *dintrf) {

      *marko  = CS_LAGR_COHERENCE_STRUCT_EJECTION;
      *vvue   =  -*vstruc + *gnorm * taup + *vnorm;
      _dep_ejection(marko,
                    depint,
                    dtp1,
                    tstruc,
                    vstruc,
                    lvisq,
                    dx,
                    vvue,
                    vpart,
                    taup,
                    yplus,
                    unif1,
                    dintrf,
                    gnorm,
                    vnorm);

      *dx     = *dx + dxaux;

    }

  }
  else {

    if (*unif1 < (dtp / *tstruc))
      *marko  = CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION;
    else
      *marko  = CS_LAGR_COHERENCE_STRUCT_SWEEP;

  }
}

/*----------------------------------------------------------------------------
 * Management of the diffusion phases (marko = CS_LAGR_COHERENCE_STRUCT_DIFFUSION)
 *
 * parameters:
 *   dx        -->    wall-normal displacement
 *   vvue      <->    wall-normal velocity of the flow seen
 *   vpart     <->    particle wall-normal velocity
 *   marko     <->    state of the jump process
 *   tempf     <--    temperature of the fluid
 *   depint    <--    interface location near-wall/core-flow
 *   dtl       -->    Lagrangian time step
 *   tstruc    <->    coherent structure mean duration
 *   tdiffu    <->    diffusion phase mean duration
 *   ttotal    <->    tdiffu + tstruc
 *   vstruc    <--    coherent structure velocity
 *   romp      <--    particle density
 *   taup      -->    particle relaxation time
 *   kdif      <->    diffusion phase diffusion coefficient
 *   tlag2     <->    diffusion relaxation timescale
 *   lvisq     <--    wall-unit lengthscale
 *   yplus     <--    wall-normal velocity of the flow seen
 *   unif1     <->    random number (uniform law)
 *   unif2     <->    random number (uniform law)
 *   dintrf    <--    extern-intern interface location
 *   rpart     <--    particle radius
 *   kdifcl    <->    internal zone diffusion coefficient
 *   indint    <->    interface indicator
 *   gnorm     <--    wall-normal gravity component
 *   vnorm     <--    wall-normal fluid (Eulerian) velocity
 *   grpn      <--    wall-normal pressure gradient
 *   piiln     <--    SDE integration auxiliary term
 *----------------------------------------------------------------------------*/

static void
_dep_diffusion_phases(cs_real_t *dx,
                      cs_real_t *vvue,
                      cs_real_t *vpart,
                      cs_lnum_t *marko,
                      cs_real_t  tempf,
                      cs_real_t *depint,
                      cs_real_t  dtl,
                      cs_real_t *tstruc,
                      cs_real_t *tdiffu,
                      cs_real_t *ttotal,
                      cs_real_t *vstruc,
                      cs_real_t  romp,
                      cs_real_t  taup,
                      cs_real_t *kdif,
                      cs_real_t *tlag2,
                      cs_real_t  lvisq,
                      cs_real_t *yplus,
                      cs_real_t *unif1,
                      cs_real_t *unif2,
                      cs_real_t *dintrf,
                      cs_real_t *rpart,
                      cs_real_t *kdifcl,
                      cs_lnum_t *indint,
                      cs_real_t *gnorm,
                      cs_real_t *vnorm,
                      cs_real_t *grpn,
                      cs_real_t *piiln)
{
  cs_real_t  tci, force, aux1, aux2, aux3, aux4, aux5, aux6;
  cs_real_t  aux7, aux8, aux9, aux10, aux11, aa, bb, cc, dd, ee, ter1x, ter2x;
  cs_real_t  ter3x, ter4x, ter5x, ter1f, ter2f, ter3f, ter1p, ter2p, ter3p;
  cs_real_t  ter4p, ter5p, gama2, omegam, omega2, p11, p21, p22, p31, p32;
  cs_real_t  p33, grga2, gagam, gaome, yplusa, dxaux, dtp1;

  cs_real_t vagaus[4];
  cs_random_normal(4, vagaus);

  cs_real_t vpart0 = *vpart;

  cs_real_t vvue0;
  if (*marko == CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION)
    vvue0 = vagaus[3] * sqrt (cs_math_pow2(*kdif) * *tlag2 / 2.0);
  else
    vvue0 = *vvue;

  tci = *piiln * *tlag2 + *vnorm;

  force = (- *grpn / romp + *gnorm) * taup;

  /* --> Coefficients and deterministic terms computation
   *     ------------------------------------------------    */

  aux1 = exp ( -dtl / taup);
  aux2 = exp ( -dtl / *tlag2);
  aux3 = *tlag2 / (*tlag2 - taup);
  aux4 = *tlag2 / (*tlag2 + taup);
  aux5 = *tlag2 * (1.0 - aux2);
  aux6 = cs_math_pow2(*kdif) * *tlag2;
  aux7 = *tlag2 - taup;
  aux8 = cs_math_pow2(*kdif) * cs_math_pow2(aux3);

  /* --> terms for the trajectory   */
  aa     = taup * (1.0 - aux1);
  bb     = (aux5 - aa) * aux3;
  cc     = dtl - aa - bb;
  ter1x  = aa * vpart0;
  ter2x  = bb * vvue0;
  ter3x  = cc * tci;
  ter4x  = (dtl - aa) * force;

  /* --> vu fluid terms   */
  ter1f  = vvue0 * aux2;
  ter2f  = tci * (1.0 - aux2);

  /* --> particles velocity terms   */
  dd     = aux3 * (aux2 - aux1);
  ee     = 1.0 - aux1;
  ter1p  = vpart0 * aux1;
  ter2p  = vvue0 * dd;
  ter3p  = tci * (ee - dd);
  ter4p  = force * ee;

  /* --> Coefficients computation for the stochastic integrals:
   * --> Integral on the particles position   */
  gama2  = 0.5 * (1.0 - aux2 * aux2);
  omegam = 0.5 * aux4 * (aux5 - aux2 * aa) - 0.5 * aux2 * bb;
  omegam = omegam * sqrt (aux6);
  omega2 =   aux7 * (aux7 * dtl - 2.0 * (*tlag2 * aux5 - taup * aa))
           + 0.5 * *tlag2 * *tlag2 * aux5 * (1.0 + aux2)
           + 0.5 * taup * taup * aa * (1.0 + aux1)
           - 2.0 * aux4 * *tlag2 * taup * taup * (1.0 - aux1 * aux2);
  omega2      = aux8 * omega2;

  if (CS_ABS(gama2) > cs_math_epzero) {

    p21  = omegam / sqrt (gama2);
    p22  = omega2 - cs_math_pow2(p21);
    p22  = sqrt (CS_MAX(0.0, p22));

  }
  else {

    p21  = 0.0;
    p22  = 0.0;

  }

  ter5x  = p21 * vagaus[0] + p22 * vagaus[1];

  /* --> vu fluid integral     */
  p11    = sqrt (gama2 * aux6);
  ter3f  = p11 * vagaus[0];

  /* --> particles velocity Integral     */
  aux9   = 0.5 * *tlag2 * (1.0 - aux2 * aux2);
  aux10  = 0.5 * taup * (1.0 - aux1 * aux1);
  aux11  = taup * *tlag2 * (1.0 - aux1 * aux2) / (taup + *tlag2);

  grga2  = (aux9 - 2.0 * aux11 + aux10) * aux8;

  gagam  = (aux9 - aux11) * (aux8 / aux3);

  gaome  = ((*tlag2 - taup) * (aux5 - aa) - *tlag2 * aux9 - taup * aux10 + (*tlag2 + taup) * aux11) * aux8;

  if (p11 > cs_math_epzero)
    p31  = gagam / p11;
  else
    p31  = 0.0;

  if (p22 > cs_math_epzero)
    p32  = (gaome - p31 * p21) / p22;
  else
    p32  = 0.0;

  p33    = grga2 - cs_math_pow2(p31) - cs_math_pow2(p32);
  p33    = sqrt (CS_MAX(0.0, p33));
  ter5p  = p31 * vagaus[0] + p32 * vagaus[1] + p33 * vagaus[2];

  /* ====================================================================
   * 3. Writings finalization
   * ====================================================================*/

  /* --> trajectory  */
  *dx = ter1x + ter2x + ter3x + ter4x + ter5x;

  /* --> seen fluid velocity     */
  *vvue = ter1f + ter2f + ter3f;

  /* --> particles velocity    */
  *vpart = ter1p + ter2p + ter3p + ter4p + ter5p;
  yplusa = *yplus - *dx / lvisq;

  /* --------------------------------------------------
   *  Dissociation of cases by the arrival position
   * -------------------------------------------------- */

  if (yplusa > *depint)
    *marko    =  -2;
  else if (yplusa < *dintrf) {

    *marko = 0;
    *vvue  = sqrt ( cs_math_pow2((*kdifcl)) * *tlag2 / 2.0) * sqrt (2.0 * cs_math_pi) * 0.5;
    *dx   *= (*dintrf - *yplus) / (yplusa - *yplus);
    dxaux  = *dx;
    *vpart = (*yplus - yplusa) * lvisq / dtl;
    dtp1   = dtl * (*dintrf - yplusa) / (*yplus - yplusa);
    *yplus = *dintrf;
    _dep_inner_zone_diffusion(dx,
                              vvue,
                              vpart,
                              marko,
                              tempf,
                              depint,
                              dtp1,
                              tstruc,
                              tdiffu,
                              ttotal,
                              vstruc,
                              romp,
                              taup,
                              kdif,
                              tlag2,
                              yplus,
                              lvisq,
                              unif1,
                              unif2,
                              dintrf,
                              rpart,
                              kdifcl,
                              indint,
                              gnorm,
                              vnorm,
                              grpn,
                              piiln);

    *dx  = dxaux + *dx;

  }
  else {

    if (*unif1 < (dtl / *tdiffu)) {

      if (*unif2 < 0.5) {

        *marko = CS_LAGR_COHERENCE_STRUCT_SWEEP;
        *vvue  = *vstruc + *gnorm * taup + *vnorm;

      }
      else {

        *marko = CS_LAGR_COHERENCE_STRUCT_EJECTION;
        *vvue  =  -*vstruc + *gnorm * taup + *vnorm;

      }
    }
    else
      *marko = CS_LAGR_COHERENCE_STRUCT_DIFFUSION;

  }
}

/*----------------------------------------------------------------------------
 * Management of the diffusion in the internal zone (y^+ < dintrf)
 *
 * parameters:
 *   dx        -->    wall-normal displacement
 *   vvue      <->    wall-normal velocity of the flow seen
 *   vpart     <->    particle wall-normal velocity
 *   marko     <->    state of the jump process
 *   tempf     <--    temperature of the fluid
 *   dtl       -->    Lagrangian time step
 *   tstruc    <->    coherent structure mean duration
 *   tdiffu    <->    diffusion phase mean duration
 *   ttotal    <->    tdiffu + tstruc
 *   vstruc    <--    coherent structure velocity
 *   romp      -->    particle density
 *   taup      <->    particle relaxation time
 *   kdif      <->    diffusion phase diffusion coefficient
 *   tlag2     <->    diffusion relaxation timescale
 *   lvisq     <--    wall-unit lengthscale
 *   yplus     <--    wall-normal velocity of the flow seen
 *   unif1     <->    random number (uniform law)
 *   unif2     <->    random number (uniform law)
 *   dintrf    <--    extern-intern interface location
 *   rpart     <--    particle radius
 *   kdifcl    <->    internal zone diffusion coefficient
 *   indint    <->    interface indicator
 *   gnorm     <--    wall-normal gravity component
 *   vnorm     <--    wall-normal fluid (Eulerian) velocity
 *   grpn      <--    wall-normal pressure gradient
 *   piiln     <--    SDE integration auxiliary term
 *----------------------------------------------------------------------------*/

static void
_dep_inner_zone_diffusion(cs_real_t *dx,
                          cs_real_t *vvue,
                          cs_real_t *vpart,
                          cs_lnum_t *marko,
                          cs_real_t  tempf,
                          cs_real_t *depint,
                          cs_real_t  dtl,
                          cs_real_t *tstruc,
                          cs_real_t *tdiffu,
                          cs_real_t *ttotal,
                          cs_real_t *vstruc,
                          cs_real_t  romp,
                          cs_real_t  taup,
                          cs_real_t *kdif,
                          cs_real_t *tlag2,
                          cs_real_t *yplus,
                          cs_real_t  lvisq,
                          cs_real_t *unif1,
                          cs_real_t *unif2,
                          cs_real_t *dintrf,
                          cs_real_t *rpart,
                          cs_real_t *kdifcl,
                          cs_lnum_t *indint,
                          cs_real_t *gnorm,
                          cs_real_t *vnorm,
                          cs_real_t *grpn,
                          cs_real_t *piiln)
{
  cs_real_t vagaus[3], vagausbr[2];
  cs_random_normal(3, vagaus);
  cs_random_normal(2, vagausbr);

  cs_real_t force  = *gnorm * taup;
  cs_real_t vvue0  = *vvue;
  cs_real_t vpart0 = *vpart;

  cs_real_t  argt, kaux, tci;
  if (*yplus < 5.0) {

    argt = cs_math_pi * *yplus / 5.0;
    kaux = *kdifcl * 0.5 * (1.0 - cos (argt));
    tci  =  -cs_math_pow2(*tlag2) * 0.5 * cs_math_pow2(*kdifcl) * cs_math_pi * sin (argt)
           * (1.0 - cos (argt)) / (2.0 * 5.0) / lvisq;

  }
  else {

    kaux      = *kdifcl;
    /* Interpolation of the decreasing normal fluid velocity around zero:     */
    tci  = *vnorm * *yplus / *dintrf;

  }

  /* -----------------------------------------------
   *  Brownian motion
   * -----------------------------------------------    */

  cs_real_t mpart    = 4.0 / 3.0 * cs_math_pi * cs_math_pow3(*rpart) * romp;
  cs_real_t kdifbr   = sqrt (2.0 * _k_boltz * tempf / (mpart * taup));
  cs_real_t kdifbrtp = kdifbr * taup;

  /* -----------------------------------------------    */

  cs_real_t dtstl  = dtl / *tlag2;
  cs_real_t dtstp  = dtl / taup;
  cs_real_t tlmtp  = *tlag2 - taup;
  cs_real_t tlptp  = *tlag2 + taup;
  cs_real_t tltp   = *tlag2 * taup;
  cs_real_t tl2    = cs_math_pow2(*tlag2);
  cs_real_t tp2    = cs_math_pow2(taup);
  cs_real_t thet   = *tlag2 / tlmtp;
  cs_real_t the2   = cs_math_pow2(thet);
  cs_real_t etl    = exp ( -dtstl);
  cs_real_t etp    = exp ( -dtstp);
  cs_real_t l1l    = 1.0 - etl;
  cs_real_t l1p    = 1.0 - etp;
  cs_real_t l2l    = 1.0 - etl * etl;
  cs_real_t l2p    = 1.0 - etp * etp;
  cs_real_t l3     = 1.0 - etl * etp;
  cs_real_t kaux2  = cs_math_pow2(kaux);
  cs_real_t k2the2 = kaux2 * the2;
  cs_real_t aa1    = taup * l1p;
  cs_real_t bb1    = thet * (*tlag2 * l1l - aa1);
  cs_real_t cc1    = dtl - aa1 - bb1;
  cs_real_t dd1    = thet * (etl - etp);
  cs_real_t ee1    = l1p;

  /* --------------------------------------------------
   * Auxiliary terms for Brownian motion
   * -------------------------------------------------- */

  cs_real_t xiubr  = 0.5 * cs_math_pow2((kdifbrtp * l1p));
  cs_real_t ucarbr = kdifbrtp * kdifbr * 0.5 * l2p;
  cs_real_t xcarbr = cs_math_pow2(kdifbrtp) * (dtl - l1p * (2.0 + l1p) * 0.5 * taup);
  cs_real_t ubr    = sqrt (CS_MAX(ucarbr, 0.0));

  /* ---------------------------------------------------
   * Deterministic terms computation
   * ---------------------------------------------------*/

  *vvue  = vvue0 * etl + tci * l1l;
  *vpart = vpart0 * etp + dd1 * vvue0 + tci * (ee1 - dd1) + force * ee1;
  *dx    = aa1 * vpart0 + bb1 * vvue0 + cc1 * tci + (dtl - aa1) * force;

  /* ---------------------------------------------------
   * Correlation matrix
   * ---------------------------------------------------*/

  cs_real_t pgam2  = 0.5 * kaux2 * *tlag2 * l2l;
  cs_real_t ggam2  = the2 * pgam2 + k2the2 * (  l3 * ( -2 * tltp / tlptp)
                                              + l2p * (taup * 0.5));
  cs_real_t ome2   = k2the2 * ( dtl * cs_math_pow2(tlmtp) + l2l * (tl2 * *tlag2 * 0.5)
                               + l2p * (tp2 * taup * 0.5)
                               + l1l * ( -2.0 * tl2 * tlmtp)
                               + l1p * (2.0 * tp2 * tlmtp)
                               + l3 * ( -2.0 * (cs_math_pow2(tltp)) / tlptp));

  cs_real_t pgagga = thet * (pgam2 - kaux2 * tltp / tlptp * l3);
  cs_real_t pgaome = thet * *tlag2 * ( -pgam2 + kaux2 * (  l1l * tlmtp
                                                         + l3 * tp2 /  tlptp));
  cs_real_t ggaome = k2the2 * (  tlmtp * ( *tlag2 * l1l + l1p * ( -taup))
                               + l2l * ( -tl2 * 0.5)
                               + l2p * ( -tp2 * 0.5) + l3 * tltp);

  /* -----------------------------------------------------
   *  Choleski decomposition
   * -----------------------------------------------------*/

  /*  P11 Computation     */
  cs_real_t p11 = sqrt (CS_MAX(0.0, pgam2));
  cs_real_t p21;

  /*  P21 and P22 computations */
  if (CS_ABS(p11) > cs_math_epzero)
    p21  = pgagga / p11;
  else
    p21  = 0.0;

  cs_real_t p22 = sqrt (CS_MAX (0.0, ggam2 - cs_math_pow2(p21)));

  /*  P31, P32 and P33 computations */
  cs_real_t p31;
  if (CS_ABS(p11) > cs_math_epzero) {
    p31  = pgaome / p11;
  }
  else {
    p31  = 0.0;
  }

  cs_real_t p32;
  if (CS_ABS(p22) > cs_math_epzero)
    p32  = (ggaome - p21 * p31) / p22;
  else
    p32  = 0.0;

  cs_real_t p33 = sqrt(CS_MAX(0.0, ome2 - p31*p31 - p32*p32));

  /* ----------------------------------------------------------
   *  Brownian motion term
   * ----------------------------------------------------------*/

  cs_real_t p11br  = ubr;

  cs_real_t p21br;
  if (CS_ABS(p11br) > cs_math_epzero)
    p21br = xiubr / p11br;
  else
    p21br = 0.0;

  cs_real_t p22br  = sqrt (CS_MAX (xcarbr - cs_math_pow2(p21br), 0.0));

  /* ----------------------------------------------------------
   *  The random terms are consequently:
   * ----------------------------------------------------------*/

  cs_real_t terf = p11 * vagaus[0];
  cs_real_t terp = p21 * vagaus[0] + p22 * vagaus[1];
  cs_real_t terx = p31 * vagaus[0] + p32 * vagaus[1] + p33 * vagaus[2];

  /* ----------------------------------------------------------
   *  Brownian motion
   * ----------------------------------------------------------*/

  cs_real_t terpbr = p11br * vagausbr[0];
  cs_real_t terxbr = p21br * vagausbr[0] + p22br * vagausbr[1];

  /* -----------------------------------------------------------
   *  Finalization (First order)
   * ----------------------------------------------------------*/

  *vvue  = *vvue + terf;
  *vpart = *vpart + terp + terpbr;
  *dx    = *dx + terx + terxbr;
  cs_real_t yplusa = *yplus - *dx / lvisq;

  if ((yplusa * lvisq) < *rpart) {

    *dx  = *dx + 2 * *rpart;
    return;

  }

  if ((yplusa > *dintrf) && (*indint != 1)) {

    *marko = CS_LAGR_COHERENCE_STRUCT_DIFFUSION;
    *vvue  =  -sqrt(cs_math_pow2((*kdifcl * (*ttotal / *tdiffu))) * *tlag2 / 2.0)
             * sqrt (2.0 * cs_math_pi) * 0.5;
    *dx   *= (*dintrf - *yplus) / (yplusa - *yplus);
    *vpart = (*yplus - yplusa) * lvisq / dtl;
    cs_real_t dxaux  = *dx;
    cs_real_t dtp1   = dtl * (*dintrf - yplusa) / (*yplus - yplusa);
    *yplus = *dintrf;
    _dep_diffusion_phases(dx,
                          vvue,
                          vpart,
                          marko,
                          tempf,
                          depint,
                          dtp1,
                          tstruc,
                          tdiffu,
                          ttotal,
                          vstruc,
                          romp,
                          taup,
                          kdif,
                          tlag2,
                          lvisq,
                          yplus,
                          unif1,
                          unif2,
                          dintrf,
                          rpart,
                          kdifcl,
                          indint,
                          gnorm,
                          vnorm,
                          grpn,
                          piiln);

    *dx  = dxaux + *dx;

  }
  else {

    if (yplusa > 0.0) {

      cs_real_t tcin1;
      cs_real_t kauxn1;

      if (yplusa < 5.0) {

        cs_real_t argtn1  = cs_math_pi * yplusa / 5.0;
        kauxn1  = *kdifcl * 0.5 * (1.0 - cos (argtn1));
        tcin1   =  cs_math_pow2(*tlag2) * 0.5 * cs_math_pow2(*kdifcl)
                   * cs_math_pi * sin (argtn1) * (1.0 - cos (argtn1))
                   / (2.0 * 5.0) / lvisq;
      }
      else {
        kauxn1     = *kdifcl;
        tcin1      = 0.0;
      }

      /* ---------------------------------------------------------------
       *   Auxiliary computation
       * ---------------------------------------------------------------*/

      cs_real_t pox1 = l1l / dtstl;
      cs_real_t pox2 = tlptp / dtl * l1p;
      cs_real_t aa2  = -etl + pox1;
      cs_real_t bb2  = 1.0 - pox1;
      cs_real_t c2c  = *tlag2 / tlmtp * (etl - etp);
      cs_real_t a2c  = -etp + pox2 - (1.0 + *tlag2 / dtl) * c2c;
      cs_real_t b2c  = 1.0 - pox2 + (*tlag2 / dtl) * c2c;
      cs_real_t a22  = l2l + l2l / (2 * dtstl) - 1.0;
      cs_real_t b22  = 1.0 - l2l / (2 * dtstl);

      /* ---------------------------------------------------------------
       *   Deterministic terms computation
       * ---------------------------------------------------------------*/

      *vvue   = vvue0 * etl + aa2 * tci + bb2 * tcin1;
      *vpart  = vpart0 * etp + vvue0 * c2c + a2c * tci + b2c * tcin1
                             + force * (1.0 - (etp - 1.0) / ( -dtstp));

      /* -------------------------------------------------
       *  Diffusion coefficient computation
       * --------------------------------------------------*/

      cs_real_t ketoi   = (a22 * kaux + b22 * kauxn1) / l2l;
      cs_real_t ketoi2  = cs_math_pow2(ketoi);

      /* --------------------------------------------------
       *  Correlation matrix computation
       * --------------------------------------------------*/

      pgam2   = 0.5 * ketoi2 * *tlag2 * l2l;
      ggam2   = the2 * (pgam2 + ketoi2 * (l3 * ( -2.0 * tltp / tlptp) + l2p *
                                          taup * 0.5));
      pgagga  = thet * (pgam2 - ketoi2 * tltp / tlptp * l3);

      /* -----------------------------------------------------
       *  Choleski decomposition
       * -----------------------------------------------------*/

      /*  P11 computation     */
      p11 = sqrt (CS_MAX (0.0, pgam2));

      /*  P21 and P22 computation  */
      if (CS_ABS(p11) > cs_math_epzero)
        p21 = pgagga / p11;
      else
        p21 = 0.0;

      p22 = sqrt (CS_MAX (0.0, ggam2 - cs_math_pow2(p21)));

      /* ----------------------------------------------------------
       *  The random terms are:
       * ----------------------------------------------------------*/

      terf = p11 * vagaus[0];
      terp = p21 * vagaus[0] + p22 * vagaus[1];

      /* -----------------------------------------------------------
       *  Finalization (Second order)
       * -----------------------------------------------------------*/

      *vvue   = *vvue + terf;
      *vpart  = *vpart + terp + terpbr;

    }

  }
}

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Deposition submodel:
 *   1/ Parameter initialization
 *   2/ Call of the different subroutines with respect to the marko indicator
 *
 * parameters:
 *   marko     <->    state of the jump process
 *   tempf     <--    temperature of the fluid
 *   lvisq     <--    wall-unit lenghtscale
 *   tvisq     <--    wall-unit timescale
 *   vpart     <--    particle wall-normal velocity
 *   vvue      <--    wall-normal velocity of the flow seen
 *   dx        <--    wall-normal displacement
 *   diamp     <--    particle diameter
 *   romp      <--    particle density
 *   taup      <--    particle relaxation time
 *   yplus     <--    particle wall-normal normalized distance
 *   dintrf    <--    extern-intern interface location
 *   enertur   <--    turbulent kinetic energy
 *   gnorm     <--    wall-normal gravity component
 *   vnorm     <--    wall-normal fluid (Eulerian) velocity
 *   grpn      <--    wall-normal pressure gradient
 *   piiln     <--    SDE integration auxiliary term
 *   depint    <--    interface location near-wall/core-flow
 *----------------------------------------------------------------------------*/

void
cs_lagr_deposition(cs_real_t  dtp,
                   cs_lnum_t *marko,
                   cs_real_t  tempf,
                   cs_real_t  lvisq,
                   cs_real_t  tvisq,
                   cs_real_t *vpart,
                   cs_real_t *vvue,
                   cs_real_t *dx,
                   cs_real_t *diamp,
                   cs_real_t  romp,
                   cs_real_t  taup,
                   cs_real_t *yplus,
                   cs_real_t *dintrf,
                   cs_real_t *enertur,
                   cs_real_t *gnorm,
                   cs_real_t *vnorm,
                   cs_real_t *grpn,
                   cs_real_t *piiln,
                   cs_real_t *depint)
{
  /* ====================================================================   */
  /* 1. Initializations   */
  /* ====================================================================   */
  /* Ratio between k and v'    */

  cs_real_t rapkvp = 0.39;

  /* The temporal parameters estimated from the DNS computations  */
  /* and written in adimensional form    */

  cs_real_t tlag2  = 3.0 * tvisq;
  cs_real_t tstruc = 30.0 * tvisq;
  cs_real_t tdiffu = 10.0 * tvisq;
  cs_real_t ttotal = tstruc + tdiffu;

  /* The velocity Vstruc is estimated as the square-root of turbulent kinetic
   * energy times rapkvp, which corresponds roughly to v' in most of the
   * turbulent boundary layer */

  cs_real_t vstruc = sqrt(*enertur * rapkvp);

  /* From Vstruc we estimate the diffusion coefficient Kdif to balance the fluxes.
     * Kdif is roughly equal to sqrt(k/(4*pi)) in the core flow (which is the
     * theoretical value of the standard Langevin model with a C0 = 2.1) such as:
     *    flux_langevin = sig / sqrt(2*pi)  = v' / sqrt(2*pi)
     * and (v') = k * C0 /( 1 + 3*C0/2 ) = approx. k/2 (see Minier & Pozorski, 1999) */

  cs_real_t tt = (sqrt(cs_math_pi * rapkvp) * tstruc);

  cs_real_t kdif = 0;

  if (ttotal > tt)
    kdif =   sqrt(*enertur / tlag2) * (ttotal - tt) / tdiffu;

  else
    bft_error(__FILE__, __LINE__, 0,
              "Incorrect parameter values in %s", __func__);

  /* Ratios computation of the flux to determine the kdifcl value */

  cs_real_t ectype = sqrt (cs_math_pow2(kdif) * tlag2 / 2.0);
  cs_real_t paux   = sqrt (cs_math_pi / 2.0) * tstruc * vstruc / (ectype * tdiffu);
  paux   = paux / (1.0 + paux);
  cs_real_t kdifcl = kdif * (tdiffu / ttotal);

  cs_real_t unif[2];
  cs_random_uniform(2, unif);

  cs_lnum_t indint = 0;

  /* ====================================================================   */
  /* 2. Treatment of the 'degenerated' cases
   * (marko =
   *   CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF,
   *   CS_LAGR_COHERENCE_STRUCT_DEGEN_SWEEP,
   *   CS_LAGR_COHERENCE_STRUCT_DEGEN_EJECTION) */
  /* ====================================================================   */

  cs_real_t unif1;
  if (*marko == CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF) {

    *marko    = CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF;
    *vvue     = 0.0;

  }
  else if (*marko == CS_LAGR_COHERENCE_STRUCT_DEGEN_SWEEP) {

    cs_random_uniform(1, &unif1);

    if (unif1 < paux)
      *marko  = CS_LAGR_COHERENCE_STRUCT_SWEEP;

    else
      *marko  = CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION;

  }
  else if (*marko == CS_LAGR_COHERENCE_STRUCT_DEGEN_EJECTION) {

    cs_random_uniform(1, &unif1);

    if (unif1 < 0.5)
      *marko  = CS_LAGR_COHERENCE_STRUCT_SWEEP;

    else
      *marko  = CS_LAGR_COHERENCE_STRUCT_EJECTION;

  }

  /* ====================================================================
   * 2. Call of different subroutines given the value of marko
   *   marko = CS_LAGR_COHERENCE_STRUCT_SWEEP --> sweep phase
   *   marko = CS_LAGR_COHERENCE_STRUCT_DIFFUSION
   *         or CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION --> diffusion phase
   *   marko = CS_LAGR_COHERENCE_STRUCT_EJECTION --> ejection phase
   *   marko = CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF --> inner-zone diffusion
   * =====================================================================  */

  cs_real_t rpart = *diamp * 0.5;

  if (*marko == CS_LAGR_COHERENCE_STRUCT_SWEEP)
    _dep_sweep(dx,
               vvue,
               vpart,
               marko,
               tempf,
               depint,
               dtp,
               &tstruc,
               &tdiffu,
               &ttotal,
               &vstruc,
               romp,
               taup,
               &kdif,
               &tlag2,
               lvisq,
               yplus,
               &(unif[0]),
               &(unif[1]),
               dintrf,
               &rpart,
               &kdifcl,
               gnorm,
               vnorm,
               grpn,
               piiln);

  else if (*marko == CS_LAGR_COHERENCE_STRUCT_DIFFUSION
      || *marko == CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION)
    _dep_diffusion_phases(dx,
                          vvue,
                          vpart,
                          marko,
                          tempf,
                          depint,
                          dtp,
                          &tstruc,
                          &tdiffu,
                          &ttotal,
                          &vstruc,
                          romp,
                          taup,
                          &kdif,
                          &tlag2,
                          lvisq,
                          yplus,
                          &(unif[0]),
                          &(unif[1]),
                          dintrf,
                          &rpart,
                          &kdifcl,
                          &indint,
                          gnorm,
                          vnorm,
                          grpn,
                          piiln);

  else if (*marko == CS_LAGR_COHERENCE_STRUCT_EJECTION)
    _dep_ejection(marko,
                  depint,
                  dtp,
                  &tstruc,
                  &vstruc,
                  lvisq,
                  dx,
                  vvue,
                  vpart,
                  taup,
                  yplus,
                  &(unif[0]),
                  dintrf,
                  gnorm,
                  vnorm);

  else if (*marko == CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF)
    _dep_inner_zone_diffusion(dx,
                              vvue,
                              vpart,
                              marko,
                              tempf,
                              depint,
                              dtp,
                              &tstruc,
                              &tdiffu,
                              &ttotal,
                              &vstruc,
                              romp,
                              taup,
                              &kdif,
                              &tlag2,
                              yplus,
                              lvisq,
                              &(unif[0]),
                              &(unif[1]),
                              dintrf,
                              &rpart,
                              &kdifcl,
                              &indint,
                              gnorm,
                              vnorm,
                              grpn,
                              piiln);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
