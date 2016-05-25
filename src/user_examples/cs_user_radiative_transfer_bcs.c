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
 *     izfrdp[face_id]) must be assigned.
 *   Warning: it is essential that ALL boundary faces
 *     have been assigned to a zone.
 *   The number of zones (the value of izfrdp[face_id]) is
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
  /*< [loc_var]*/
  cs_real_t tkelvi = 273.15;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  /*< [loc_var]*/

  /*< [allocate]*/
  /* Allocate a temporary array for boundary faces selection */

  cs_lnum_t nlelt;
  cs_lnum_t  *lstelt;
  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);
  /*< [allocate]*/

  /*< [ivar]*/
  cs_field_t *fth;

  switch (cs_glob_thermal_model->itherm) {
  case 1:
    fth = CS_F_(t);
    break;
  case 2:
    fth = CS_F_(h);
    break;
  default:
    fth = NULL;
  }

  const cs_lnum_t ivart
    = cs_field_get_key_int(fth, cs_field_key_id("variable_id")) - 1;
  /*< [ivar]*/

  /* Min and max values for the wall temperatures (clipping otherwise)
   * TMIN and TMAX are given in Kelvin. */

  /*< [temp]*/
  *tmin = 0.0;
  *tmax = cs_math_big_r + tkelvi;
  /*<[temp]*/

  /* Zone definitions */
  /*------------------*/

   /* Example: for wall boundary faces, selection criteria: color 1;
    * gray or black wall with profile of fixed inside temperature
    * ------------------------------------------------------------*/

  /*< [example_1]*/
  cs_selector_get_b_face_list("1",
                              &nlelt,
                              lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {
      /* zone number */
      izfrdp[face_id] = 51;

      /* Type of condition: gray or black wall with fixed inside temperature */
      isothp[face_id] = cs_glob_rad_transfer_params->itpimp;

      /* Emissivity */
      epsp[face_id] = 0.1;

      /* Fixed inside temperature */
      tintp[face_id] = 200 + tkelvi;
    }

  }
  /*< [example_1]*/

  /* Example: for wall boundary faces, selection criteria: color 2;
   * gray or black wall with fixed outside temperature TEXTP
   * --------------------------------------------------------*/

  /*< [example_2]*/
  cs_selector_get_b_face_list("2",
                              &nlelt,
                              lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_ROUGHWALL) {
      /* zone number */
      izfrdp[face_id] = 52;

      /* Type of condition: gray or black wall with fixed
         outside temperature TEXTP */
      isothp[face_id] = cs_glob_rad_transfer_params->ipgrno;

      /* Emissivity */
      epsp[face_id] = 0.9;

      /* Conductivity (W/m/K)*/
      xlamp[face_id] = 3.0;

      /* Thickness (m)*/
      epap[face_id] = 0.1;

      /* Fixed outside temperature: 473.16 K */
      textp[face_id] = 200. + tkelvi;

      /* Initial inside temperature: 473.16 K */
      tintp[face_id] = 200. + tkelvi;
    }

  }
  /*< [example_2]*/

  /* Example: for wall boundary faces, selection criteria: color 3
   * reflecting wall (EPSP = 0) with fixed outside temperature TEXTP
   * --------------------------------------------------------------- */

  /*< [example_3]*/
  cs_selector_get_b_face_list("3",
                                  &nlelt,
                                  lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {
      /* zone number */
      izfrdp[face_id] = 53;

      /* Type of condition: reflecting wall with fixed outside temperature TEXTP */
      isothp[face_id] = cs_glob_rad_transfer_params->iprefl;

      /* Conductivity (W/m/K) */
      xlamp[face_id] = 3.0;

      /* Thickness    (m)*/
      epap[face_id] = 0.10;

      /* Fixed outside temperature: 473.16 K */
      textp[face_id] = 200.0 + tkelvi;

      /* Initial inside temperature: 473.16 K */
      tintp[face_id] = 200.0 + tkelvi;
    }

  }
  /*< [example_3]*/

  /* Example: for wall boundary faces which have the color 4:
   * gray or black wall and fixed conduction flux through the wall
   *
   *        XLAMP
   *        -----(Tparop-Textp) = fixed conduction flux     (W/m2)
   *        EPAP
   *                         = RODCL(FACE_ID,IVAR,3)

   *       If the conduction flux is zero then the wall is adiabatic.
   *       The array RCODCL(FACE_ID,IVAR,3) has the value of the flux.
   *       Flux density (< 0 if gain for the fluid)
   *         For temperatures T,    in Watt/m2:
   *            RCODCL(FACE_ID,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
   *         For enthalpies H,      in Watt/m2:
   *            RCODCL(FACE_ID,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
   */

  /*< [example_4]*/
  cs_selector_get_b_face_list("4",
                              &nlelt,
                              lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {
      /* zone number */
      izfrdp[face_id] = 54;

      /* Type of condition: gray or black wall with fixed conduction
         flux through the wall */
      isothp[face_id] = cs_glob_rad_transfer_params->ifgrno;

      /* Emissivity */
      epsp[face_id] = 0.9;

      /* Conduction flux (W/m2) */
      rcodcl[face_id + ivart * n_b_faces + 2 * nvarcl * n_b_faces ] = 0.0;

      /* Initial inside temperature: 473.16 K */
      tintp[face_id] = 200.0 + tkelvi;
    }

  }
  /*< [example_4]*/

  /* Example: for wall boundary faces which have the color 5:
   * reflecting wall and fixed conduction flux through the wall
   *
   *      Equivalent to impose a Neumann to the fluid
   *
   *        XLAMP
   *        -----(Tparop-Textp) = fixed conduction flux and EPSP = 0
   *        EPAP
   *                         = RODCL(FACE_ID,IVAR,3)

   *       If the conduction flux is zero then the wall is adiabatic.
   *       Flux density (< 0 if gain for the fluid)
   *         For temperatures T,    in Watt/m2:
   *            RCODCL(FACE_ID,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
   *         For enthalpies H,      in Watt/m2:
   *            RCODCL(FACE_ID,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
   */

  /*< [example_5]*/
  cs_selector_get_b_face_list("5",
                              &nlelt,
                              lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {
      /* zone number */
      izfrdp[face_id] = 55;

      /* Type of condition: reflecting wall with fixed conduction
         flux through the wall */
      isothp[face_id] = cs_glob_rad_transfer_params->ifrefl;

      /* Conduction flux (W/m2)*/
      rcodcl[face_id + ivart * n_b_faces + 2 * nvarcl * n_b_faces ] = 0.0;

      /* Initial inside temperature: 473.16 K */
      tintp[face_id] = 200.0 + tkelvi;
    }

  }
  /*< [example_5]*/

  /* WARNING: for all boundary faces, even when not a wall, it is MANDATORY to
   *          assign a zone number in the array izfrdp. */

  /* Example for assigning zones */

  for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_b_faces; face_id++) {

    if (bc_type[face_id] == CS_OUTLET)
      izfrdp[face_id] = 60;
    else if (bc_type[face_id] == CS_FREE_INLET)
      izfrdp[face_id] = 61;
    else if (bc_type[face_id] == CS_INLET)
      izfrdp[face_id] = 62;
    else if (bc_type[face_id] == CS_CONVECTIVE_INLET)
      izfrdp[face_id] = 63;
    else if (bc_type[face_id] == CS_SYMMETRY)
      izfrdp[face_id] = 64;

  }

  /* Deallocate the temporary array */
  BFT_FREE(lstelt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
