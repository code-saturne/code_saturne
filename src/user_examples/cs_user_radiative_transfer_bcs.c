/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* VERS */

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

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
 * See \ref cs_user_radiative_transfer for examples.
 *
 * \warning the temperature unit here is the Kelvin
 *
 * \section cs_user_radiative_transfer_bcs_zones  Zone definitions
 *
 * For each boundary face face_id, a specific output (logging and
 * postprocessing) class id may be assigned. This allows realizing balance
 * sheets by treating them separately for each zone. By default, the
 * output class id is set to the general (input) zone id associated to a face.
 *
 * To access output class ids (both for reading and modifying), use the
 * \ref cs_boundary_zone_face_class_id function.
 * The zone id values are arbitrarily chosen by the user, but must be
 * positive integers; very high numbers may also lead to higher memory
 * consumption.
 *
 * \section cs_user_radiative_transfer_bcs_wall  Wall characteristics
 *
 * The following face characteristics must be set:
 *  - isothp(face_id) boundary face type
 *    * CS_BOUNDARY_RAD_WALL_GRAY:
 *      Gray wall with temperature based on fluid BCs
 *    * CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T:
 *      Gray wall with fixed outside temperature
 *    * CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T:
 *      Reflecting wall with fixed outside temperature
 *    * CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX:
 *      Gray wall with fixed conduction flux
 *    * CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX:
 *      Reflecting wall with fixed conduction flux
 *
 * Depending on the value of isothp, other values may also need to be set:
 *  - rcodcl = conduction flux
 *  - epsp   = emissivity
 *  - xlamp  = conductivity (W/m/K)
 *  - epap   = thickness (m)
 *  - textp  = outside temperature (K)
 *
 * \param[in, out]  domain        pointer to a cs_domain_t structure
 * \param[in]       bc_type       boundary face types
 * \param[out]      isothp        boundary face type for radiative transfer
 * \param[out]      tmin          min allowed value of the wall temperature
 * \param[out]      tmax          max allowed value of the wall temperature
 * \param[in]       tx            relaxation coefficient (0 < tx < 1)
 * \param[in]       dt            time step (per cell)
 * \param[in]       thwall        inside current wall temperature (K)
 * \param[in]       qincid        radiative incident flux  (W/m2)
 * \param[in]       hfcnvp        convective exchange coefficient (W/m2/K)
 * \param[in]       flcnvp        convective flux (W/m2)
 * \param[out]      xlamp         conductivity (W/m/K)
 * \param[out]      epap          thickness (m)
 * \param[out]      epsp          emissivity (>0)
 * \param[out]      textp         outside temperature (K)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_radiative_transfer_bcs(cs_domain_t      *domain,
                               const int         bc_type[],
                               int               isothp[],
                               cs_real_t        *tmin,
                               cs_real_t        *tmax,
                               cs_real_t        *tx,
                               const cs_real_t   dt[],
                               const cs_real_t   thwall[],
                               const cs_real_t   qincid[],
                               cs_real_t         hfcnvp[],
                               cs_real_t         flcnvp[],
                               cs_real_t         xlamp[],
                               cs_real_t         epap[],
                               cs_real_t         epsp[],
                               cs_real_t         textp[])
{
  /*< [loc_var]*/
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  const cs_zone_t *zone = NULL;
  int *izfrdp = cs_boundary_zone_face_class_id();
  /*< [loc_var]*/

  /*< [ivar]*/
  cs_field_t *f_th = cs_thermal_model_field();
  /*< [ivar]*/

  /* Min and max values for the wall temperatures (clipping otherwise)
   * TMIN and TMAX are given in Kelvin. */

  /*< [temp]*/
  *tmin = 0.0;
  *tmax = cs_math_big_r + tkelvi;
  /*<[temp]*/

   /* Example: for wall boundary faces, in zone named "wall_1";
    * gray or black wall with profile of fixed inside temperature
    * ------------------------------------------------------------*/

  /*< [example_1]*/
  zone = cs_boundary_zone_by_name("wall_1");

  for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

    cs_lnum_t face_id = zone->elt_ids[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {

      /* logging zone number */
      izfrdp[face_id] = 51;

      /* Type of condition: gray or black wall with fixed interior
         temperature (based on main temperature BC's) */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_GRAY;

      /* Emissivity */
      epsp[face_id] = 0.1;

    }

  }
  /*< [example_1]*/

  /* Example: for wall boundary faces, zone named "wall_2";
   * gray or black wall with fixed outside temperature TEXTP
   * --------------------------------------------------------*/

  /*< [example_2]*/
  zone = cs_boundary_zone_by_name("wall_2");

  for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

    cs_lnum_t face_id = zone->elt_ids[ilelt];

    if (bc_type[face_id] == CS_ROUGHWALL) {

      /* logging zone number */

      izfrdp[face_id] = 52;

      /* Gray or black wall with fixed exterior temperature */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T;

      /* Emissivity */
      epsp[face_id] = 0.9;

      /* Conductivity (W/m/K)*/
      xlamp[face_id] = 3.0;

      /* Thickness (m)*/
      epap[face_id] = 0.1;

      /* Fixed exterior temperature: 473.15 K */
      textp[face_id] = 200. + tkelvi;

    }

  }
  /*< [example_2]*/

  /* Example: for wall boundary faces, zone named "wall_3";
   * reflecting wall (EPSP = 0) with fixed outside temperature TEXTP
   * --------------------------------------------------------------- */

  /*< [example_3]*/
  zone = cs_boundary_zone_by_name("wall_3");

  for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

    cs_lnum_t face_id = zone->elt_ids[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {

      /* log zone number */
      izfrdp[face_id] = 53;

      /* Type of condition: reflecting wall with fixed outside temperature */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T;

      /* Conductivity (W/m/K) */
      xlamp[face_id] = 3.0;

      /* Thickness    (m)*/
      epap[face_id] = 0.10;

      /* Fixed outside temperature: 473.15 K */
      textp[face_id] = 200.0 + tkelvi;

    }

  }
  /*< [example_3]*/

  /* Example: for wall boundary faces of zone named "wall_4":
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
  zone = cs_boundary_zone_by_name("wall_4");

  for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

    cs_lnum_t face_id = zone->elt_ids[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {

      /* log zone number */
      izfrdp[face_id] = 54;

      /* Type of condition: gray or black wall with fixed conduction
         flux through the wall */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX;

      /* Emissivity */
      epsp[face_id] = 0.9;

      /* Conduction flux (W/m2) */
      f_th->bc_coeffs->rcodcl3[face_id] = 0.0;

    }

  }
  /*< [example_4]*/

  /* Example: for wall boundary faces of zone named "wall_5":
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
  zone = cs_boundary_zone_by_name("wall_5");

  for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

    cs_lnum_t face_id = zone->elt_ids[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {

      /* log zone number */
      izfrdp[face_id] = 55;

      /* Type of condition: reflecting wall with fixed conduction
         flux through the wall */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX;

      /* Conduction flux (W/m2)*/
      f_th->bc_coeffs->rcodcl3[face_id] = 0.0;

    }

  }
  /*< [example_5]*/

  /*< [example_6]*/

  /* Example:
   * For wall boundary faces which were marked with cs_user_1d_wall_thermal.c
   * Heat transfer solved in a gray wall exposed to a radiative and convective
   * flux with the 1D wall thermal module */

  if (cs_glob_1d_wall_thermal->nfpt1d > 0) {
    for (cs_lnum_t ii = 0 ; ii < cs_glob_1d_wall_thermal->nfpt1d; ii++) {
      cs_lnum_t face_id = cs_glob_1d_wall_thermal->ifpt1d[ii]-1;

      /* Zone number */
      izfrdp[face_id] = 56;

      /* Type of condition: heat transfer equation solved */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_GRAY_1D_T;

      /* Emissivity */
      epsp[face_id] = 0.9;

    }
  }
  /*< [example_6]*/

  /* Example for assigning specific output zones to non-wall faces */

  for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_b_faces; face_id++) {

    if (bc_type[face_id] == CS_OUTLET)
      izfrdp[face_id] = 60;
    else if (bc_type[face_id] == CS_FREE_INLET)
      izfrdp[face_id] = 61;
    else if (   bc_type[face_id] == CS_INLET
             || bc_type[face_id] == CS_CONVECTIVE_INLET)
      izfrdp[face_id] = 62;
    else if (bc_type[face_id] == CS_SYMMETRY)
      izfrdp[face_id] = 63;

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
