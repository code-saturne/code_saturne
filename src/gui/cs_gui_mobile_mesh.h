#ifndef __CS_GUI_MOBILE_MESH_H__
#define __CS_GUI_MOBILE_MESH_H__

/*============================================================================
 * Management of the GUI parameters file: mobile mesh
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * ALE method.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALIN()
 * *****************
 *
 * nalinf  <->   number of subiterations of initialization of
 *               fluid
 * nalimx  <->   max number of iterations of implicitation of
 *               the displacement of the structures
 * epalim  <->   realtive precision of implicitation of
 *               the displacement of the structures
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialin, UIALIN) (int    *nalinf,
                                int    *nalimx,
                                double *epalim);

/*----------------------------------------------------------------------------
 * ALE diffusion type
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALVM
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialvm, UIALVM) (void);

/*-----------------------------------------------------------------------------
 * uialcl
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALCL
 * *****************
 *
 * parameters:
 *   ialtyb       --> ialtyb
 *   impale       --> uialcl_fixed_displacement
 *   disale       --> See uialcl_fixed_displacement
 *   iuma         <-- See uialcl_fixed_velocity
 *   ivma         <-- See uialcl_fixed_velocity
 *   iwma         <-- See uialcl_fixed_velocity
 *   rcodcl       --> See uialcl_fixed_velocity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialcl, UIALCL) (const int *const    ibfixe,
                                const int *const    igliss,
                                const int *const    ivimpo,
                                const int *const    ifresf,
                                int       *const    ialtyb,
                                int       *const    impale,
                                cs_real_3_t        *disale,
                                const int *const    iuma,
                                const int *const    ivma,
                                const int *const    iwma,
                                double *const       rcodcl);

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling. Called once at initialization
 *
 * Fortran Interface:
 *
 * subroutine uistr1
 * *****************
 *
 * parameters:
 *   idfstr   --> Structure definition
 *   mbstru   <-- number of previous structures (-999 or by restart)
 *   aexxst   --> Displacement prediction alpha
 *   bexxst   --> Displacement prediction beta
 *   cfopre   --> Stress prediction alpha
 *   ihistr   --> Monitor point synchronisation
 *   xstr0    <-> Values of the initial displacement
 *   xstreq   <-> Values of the equilibrium displacement
 *   vstr0    <-> Values of the initial velocity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uistr1, UISTR1) (cs_lnum_t        *idfstr,
                                const int        *mbstru,
                                double           *aexxst,
                                double           *bexxst,
                                double           *cfopre,
                                int              *ihistr,
                                double           *xstr0,
                                double           *xstreq,
                                double           *vstr0);

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling. Called at each step
 *
 * Fortran Interface:
 *
 * SUBROUTINE UISTR2
 * *****************
 *
 * parameters:
 * xmstru       --> Mass matrix
 * xcstr        --> Damping matrix
 * xkstru       --> Stiffness matrix
 * forstr       --> Fluid force matrix
 * dtref        <--  time step
 * ttcabs       <-- current time
 * ntcabs       <-- current iteration number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uistr2, UISTR2) ( double *const  xmstru,
                                 double *const  xcstru,
                                 double *const  xkstru,
                                 double *const  forstr,
                                 double *const  dtref,
                                 double *const  ttcabs,
                                 int *const     ntcabs);

/*-----------------------------------------------------------------------------
 * Retrieve data for external coupling
 *
 * parameters:
 *   nfabor    <-- Number of boundary faces
 *   idfstr    <-- Structure definition
 *   asddlf    <-- Block of the DDL forces
 *----------------------------------------------------------------------------*/

void
CS_PROCF(uiaste, UIASTE)(int       *idfstr,
                         cs_int_t  *asddlf);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the viscosity's type of ALE method
 *
 * parameters:
 *   type        <--  type of viscosity's type
 *----------------------------------------------------------------------------*/

void
cs_gui_get_ale_viscosity_type(int  *type);

/*----------------------------------------------------------------------------
 * Mesh viscosity setting.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_viscosity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the fixed velocity for a boundary
 *
 * \param[in]  label boundary condition label
 * \param[out] vel   imposed mesh velocity
 *----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_get_fixed_velocity(const char*    label,
                                      cs_real_t     *vel);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_MOBILE_MESH_H__ */
