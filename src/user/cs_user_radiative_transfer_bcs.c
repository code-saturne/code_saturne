/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <string.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_gui_util.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"
#include "cs_rad_transfer.h"
#include "cs_thermal_model.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_user_radiative_transfer_bcs.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of radiative transfer boundary conditions.
 *
 * See \subpage cs_user_radiative_transfer for examples.
 *
 * \warning the temperature unit here is the Kelvin
 *
 * \section cs_user_radiative_transfer_bcs_zones  Zone definitions
 *
 *   We define zones of wall boundaries, and we assign a type.
 *     This allows to apply the boundary conditions and realize
 *     balance sheets by treating them separately for each zone.
 *   For each boundary face face_id (not just wall faces) a zone number
 *     IZFRDP(face_id) must be assigned.
 *   Warning: it is essential that ALL boundary faces
 *     have been assigned to a zone.
 *   The number of zones (the value of IZFRDP(face_id)) is
 *     arbitrarily chosen by the user, but must be a positive integer
 *     less than or equal to cs_glob_rad_transfer_params->nbzrdm
 *     (value set in parameter cs_user_radiation_parameters.h).
 *
 \section cs_user_radiative_transfer_bcs_wall  Wall characteristics
 *
 * The following face characteristics must be set:
 *  - isothp(face_id) boundary face type
 *               = itpimp -> Gray wall with fixed inside temperature
 *               = ipgrno -> Gray wall with fixed outside temperature
 *               = iprefl -> Reflecting wall with fixed outside temperature
 *               = ifgrno -> Gray wall with fixed conduction flux
 *               = ifrefl -> Reflecting wall with fixed conduction flux
 *  - tintp(face_id) inside wall temperature (Kelvin)
 *               initialize thwall at the first time step.
 *               If isothp = itpimp, the value of thwall is fixed to tintp
 *               In the other case, tintp is only for initialization.
 *
 * Depending on the value of isothp, other values may also need to be set:
 *  - rcodcl = conduction flux
 *  - epsp   = emissivity
 *  - xlamp  = conductivity (W/m/K)
 *  - epap   = thickness (m)
 *  - textp  = outside temperature (K)
 *
 * \param[in]     nvarcl        total number of variable BC's
 * \param[in]     bc_type       boundary face types
 * \param[in]     icodcl        boundary face code
 *                                - 1  -> Dirichlet
 *                                - 2  -> convective outlet
 *                                - 3  -> flux density
 *                                - 4  -> sliding wall and u.n=0 (velocity)
 *                                - 5  -> friction and u.n=0 (velocity)
 *                                - 6  -> roughness and u.n=0 (velocity)
 *                                - 9  -> free inlet/outlet (velocity)
 *                                inflowing possibly blocked
 * \param[in]     izfrdp        boundary faces -> zone number
 * \param[in]     isothp        boundary face type for radative transfer
 *                                - itpimp -> Gray wall with fixed inside temp
 *                                - ipgrno -> Gray wall with fixed outside temp
 *                                - iprefl -> Reflecting wall with fixed
 *                                         outside temp
 *                                - ifgrno -> Gray wall with fixed
 *                                      conduction flux
 *                                - ifrefl -> Reflecting wall with fixed
 *                                      conduction flux
 * \param[out]    tmin          min allowed value of the wall temperature
 * \param[out]    tmax          max allowed value of the wall temperature
 * \param[in]     tx            relaxation coefficient (0 < tx < 1)
 * \param[in]     dt            time step (per cell)
 * \param[in]     rcodcl        boundary condition values
 *                                rcodcl(3) = flux density value
 *                                (negative for gain) in W/m2
 * \param[in]     thwall        inside current wall temperature (K)
 * \param[in]     qincid        radiative incident flux  (W/m2)
 * \param[in]     hfcnvp        convective exchange coefficient (W/m2/K)
 * \param[in]     flcnvp        convective flux (W/m2)
 * \param[out]    xlamp         conductivity (W/m/K)
 * \param[out]    epap          thickness (m)
 * \param[out]    epsp          emissivity (>0)
 * \param[out]    textp         outside temperature (K)
 * \param[out]    tintp         initial inside temperature (K)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_radiative_transfer_bcs(int               nvarcl,
                               const int         bc_type[],
                               int               icodcl[],
                               int               izfrdp[],
                               int               isothp[],
                               cs_real_t        *tmin,
                               cs_real_t        *tmax,
                               cs_real_t        *tx,
                               const cs_real_t   dt[],
                               cs_real_t         rcodcl[],
                               const cs_real_t   thwall[],
                               const cs_real_t   qincid[],
                               cs_real_t         hfcnvp[],
                               cs_real_t         flcnvp[],
                               cs_real_t         xlamp[],
                               cs_real_t         epap[],
                               cs_real_t         epsp[],
                               cs_real_t         textp[],
                               cs_real_t         tintp[])
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS

