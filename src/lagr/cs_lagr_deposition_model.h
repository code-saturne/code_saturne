#ifndef __CS_LAGR_DEPOSITION_MODEL_H__
#define __CS_LAGR_DEPOSITION_MODEL_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
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
                   cs_real_t *depint);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_DEPOSITION_MODEL_H__ */
