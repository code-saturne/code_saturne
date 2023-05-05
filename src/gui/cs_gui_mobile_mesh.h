#ifndef __CS_GUI_MOBILE_MESH_H__
#define __CS_GUI_MOBILE_MESH_H__

/*============================================================================
 * Management of the GUI parameters file: mobile mesh
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include "cs_domain.h"

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
 * epalim  <->   relative precision of implicitation of
 *               the displacement of the structures
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialin, UIALIN) (int    *nalinf,
                                int    *nalimx,
                                double *epalim);

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
 *----------------------------------------------------------------------------*/

void cs_gui_mobile_mesh_boundary_conditions(int         *const  ialtyb,
                                            int         *const  impale,
                                            cs_real_3_t        *disale);


/*-----------------------------------------------------------------------------
 * Retrieve data for external coupling
 *
 * parameters:
 *   nfabor    <-- Number of boundary faces
 *   idfstr    <-- Structure definition
 *----------------------------------------------------------------------------*/

void
CS_PROCF(uiaste, UIASTE)(int  *idfstr);

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
 * Set ALE diffusion type from GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_ale_diffusion_type(void);

/*----------------------------------------------------------------------------
 * Mesh viscosity setting.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_viscosity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Translate the user settings for the domain boundaries into a
 *         structure storing the ALE boundaries (New mechanism used in CDO)
 *
 * \param[in, out]  domain   pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_get_boundaries(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the fixed velocity for a boundary
 *
 * \param[in]  label boundary condition label
 *
 * \return a pointer to an array of cs_real_t values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_gui_mobile_mesh_get_fixed_velocity(const char  *label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add mobile structures based on GUI BC definitions.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_structures_add(void);

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling. Called once at initialization
 *
 * parameters:
 *   mbstru   <-- number of previous structures (-999 or by restart)
 *   aexxst   --> Displacement prediction alpha
 *   bexxst   --> Displacement prediction beta
 *   cfopre   --> Stress prediction alpha
 *   ihistr   --> Monitor point synchronisation
 *   xstr0    <-> Values of the initial displacement
 *   xstreq   <-> Values of the equilibrium displacement
 *   vstr0    <-> Values of the initial velocity
 *----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_init_structures(int         mbstru,
                                   double     *aexxst,
                                   double     *bexxst,
                                   double     *cfopre,
                                   int        *ihistr,
                                   double     *xstr0,
                                   double     *xstreq,
                                   double     *vstr0);

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling. Called at each step
 *
 * parameters:
 *   xmstru       --> Mass matrix
 *   xcstr        --> Damping matrix
 *   xkstru       --> Stiffness matrix
 *   forstr       --> Fluid force matrix
 *----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_internal_structures(cs_real_t  xmstru[][3][3],
                                       cs_real_t  xcstru[][3][3],
                                       cs_real_t  xkstru[][3][3],
                                       cs_real_t  forstr[][3]);

/*-----------------------------------------------------------------------------
 * Retrieve structure id associated to faces for structure coupling
 *
 * parameters:
 *   idfstr    <-- structure number associated to each boundary face.
 *----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_bc_structures(int  *idfstr);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_MOBILE_MESH_H__ */
