#ifndef __COMMUNICATION_H__
#define __COMMUNICATION_H__

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

#ifdef __cplusplus
extern "C" {
#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive geometric data and sets variables nb_for, nb_dyn and lref
 *----------------------------------------------------------------------------*/

int
recv_geom(void *component);

/*----------------------------------------------------------------------------
 * Send geometric data to Code_Aster
 * To be removed in initialization stage
 *----------------------------------------------------------------------------*/

int
send_geom(void* component);

/*----------------------------------------------------------------------------
 * Sends time step computed by middle component to codes
 *----------------------------------------------------------------------------*/

int
send_pdt(void *component,
         double dt,
         int numpdt);

/*----------------------------------------------------------------------------
 * Receives time steps from Code_Aster and Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_pdt(void *component,
         double *dt_ast,
         double *dt_sat,
         int numpdt);

/*----------------------------------------------------------------------------
 * Sends the following parameters:
 *                                 nbpdtm
 *                                 nbssit
 *                                 epsilo
 *                                 isyncp
 *                                 ntchr
 *                                 ttpabs
 *----------------------------------------------------------------------------*/

int
send_param(void *component);

/*----------------------------------------------------------------------------
 * Receives displacements and velocities from Code_Aster at current time step
 *----------------------------------------------------------------------------*/

int
recv_dyn(void *component);

/*----------------------------------------------------------------------------
 * Send predicted displacements to Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_dyn(void *component);

/*----------------------------------------------------------------------------
 * Receive efforts from Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_for(void *component);

/*----------------------------------------------------------------------------
 * Send predicted efforts to Code_Aster
 *----------------------------------------------------------------------------*/

int
send_for(void *component);

/*----------------------------------------------------------------------------
 * Send convergence indicator to Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_icv1(void *component,
          int icv);

/*----------------------------------------------------------------------------
 * Receive convergence indicator from Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_icv(void *component,
         int *icv);

/*----------------------------------------------------------------------------
 * Send convergence indicator to Code_Aster
 *----------------------------------------------------------------------------*/

int
send_icv2(void *component,
          int icv);

/*----------------------------------------------------------------------------
 * Initialize communication with Calcium
 *----------------------------------------------------------------------------*/

int
inicom(void *component);

/*----------------------------------------------------------------------------
 * End communication with Calcium and stop calculation
 *----------------------------------------------------------------------------*/

int
calfin(void *component);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* __COMMUNICATION_H__ */
