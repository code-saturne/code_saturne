/*============================================================================
 * User definitions for fluid-structure interaction using ALE.
 *============================================================================*/

/* VERS */

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
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_fluid_structure_interaction.c
 *
 * \brief User-defined functions dedicated to Fluid-Structure interaction
 *        modeling.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definition of internal mobile structures and corresponding initial
 *        conditions (initial displacement and velocity ).
 *
 * \param[in]       is_restart         indicate if computation is restarted
 * \param[in]       n_structs          number of mobile structures
 * \param[in, out]  plot;              monitoring format mask
 *                                       0: no plot
 *                                       1: plot to text (.dat) format
 *                                       2: plot to .csv format
 *                                       3: plot to both formats
 * \param[in, out]  plot_time_control  plot time output frequency control
 * \param[in, out]  aexxst             coefficient for predicted displacement
 * \param[in, out]  bexxst             coefficient for predicted displacement
 * \param[in, out]  cfopre             coefficient for predicted force
 * \param[in, out]  xstr0              initial displacement per structure
 * \param[in, out]  vstr0              initial velocity per structure
 * \param[in, out]  xstreq             displacement of initial mesh relative to
 *                                     structures position at equilibrium
 *
 * \param[in, out]  plot_time_control  time control associated to plotting
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_fsi_structure_define
void
cs_user_fsi_structure_define(int                 is_restart,
                             int                 n_structs,
                             int                *plot,
                             cs_time_control_t  *plot_time_control,
                             cs_real_t          *aexxst,
                             cs_real_t          *bexxst,
                             cs_real_t          *cfopre,
                             cs_real_t           xstr0[][3],
                             cs_real_t           vstr0[][3],
                             cs_real_t           xstreq[][3])
{
  CS_NO_WARN_IF_UNUSED(is_restart);
  CS_NO_WARN_IF_UNUSED(n_structs);
  CS_NO_WARN_IF_UNUSED(plot);
  CS_NO_WARN_IF_UNUSED(plot_time_control);
  CS_NO_WARN_IF_UNUSED(aexxst);
  CS_NO_WARN_IF_UNUSED(bexxst);
  CS_NO_WARN_IF_UNUSED(cfopre);
  CS_NO_WARN_IF_UNUSED(xstr0);
  CS_NO_WARN_IF_UNUSED(vstr0);
  CS_NO_WARN_IF_UNUSED(xstreq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Time-based settings for internal mobile structures.
 *
 * \param[in]       n_structs  number of mobile structures
 * \param[in]       ts         time step structure
 * \param[in]       xstreq     displacement of initial mesh rel. to equilibrium
 * \param[in]       xstr       structural displacement
 * \param[in]       vstr       structural velocity
 * \param[in, out]  xmstru     matrix of structural mass
 * \param[in, out]  xcstru     matrix of structural friction
 * \param[in, out]  xkstru     matrix of structural stiffness
 * \param[in, out]  forstr     forces acting on structures (take forces)
 * \param[in, out]  dtstr      structural time step
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_fsi_structure_values
void
cs_user_fsi_structure_values(int                    n_structs,
                             const cs_time_step_t  *ts,
                             const cs_real_t        xstreq[][3],
                             const cs_real_t        xstr[][3],
                             const cs_real_t        vstr[][3],
                             cs_real_t              xmstru[][3][3],
                             cs_real_t              xcstru[][3][3],
                             cs_real_t              xkstru[][3][3],
                             cs_real_t              forstr[][3],
                             cs_real_t              dtstr[])
{
  CS_NO_WARN_IF_UNUSED(n_structs);
  CS_NO_WARN_IF_UNUSED(ts);
  CS_NO_WARN_IF_UNUSED(xstreq);
  CS_NO_WARN_IF_UNUSED(xstr);
  CS_NO_WARN_IF_UNUSED(vstr);
  CS_NO_WARN_IF_UNUSED(xmstru);
  CS_NO_WARN_IF_UNUSED(xcstru);
  CS_NO_WARN_IF_UNUSED(xkstru);
  CS_NO_WARN_IF_UNUSED(forstr);
  CS_NO_WARN_IF_UNUSED(dtstr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define structure numbers for faces associated with internal
 *        or external (code_aster) structures.
 *
 * Structure numbers associated to a given face have the following values:
 * - -i where coupled to  i-th (1-to n) external (code_aster) structure.
 * - 0 where not coupled with an internal or external structure.
 * - i  where coupled to  i-th (1-to n) internal (mass-spring) structure.
 *
 * \param[in, out]  domain         pointer to a cs_domain_t structure
 * \param[in, out]  structure_num  structure id associated to each face
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_fsi_structure_num
void
cs_user_fsi_structure_num(cs_domain_t  *domain,
                          int           structure_num[])
{
  CS_UNUSED(domain);
  CS_UNUSED(structure_num);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
