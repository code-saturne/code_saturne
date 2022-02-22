/*============================================================================
 * Phase thermodynamic for the compressible homogeneous two-phase model
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cf_thermo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hgn_phase_thermo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_hgn_phase_thermo.c
 *
 *  \brief Phase thermodynamic for compressible homogeneous two-phase model.
 *
 *  Wrappers of thermodynamic relation function are implemented here.
 */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static variables
 *============================================================================*/

/* Stiffened gas parameters and associated pointer */

static cs_stiffened_gas_t _stiffened_gas[2] = {
  {
    .cv    = -1.,
    .gamma = -1.,
    .pinf  = -1.,
    .qprim = 0.,
    .q     = 0.
  },
  {
    .cv    = -1.,
    .gamma = -1.,
    .pinf  = -1.,
    .qprim = 0.,
    .q     = 0.
  }
};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define stiffened gas parameters for phase iph.
 *
 * \param[in]  iph        index of phase (0 or 1)
 * \param[in]  cv         heat capacity (in J/kg/K)
 * \param[in]  gamma      polytropic coefficient (-)
 * \param[in]  pinf       minimum pressure (in Pa)
 * \param[in]  qprim      entropy parameter (in J/kg/K)
 * \param[in]  q          reference for specific enthalpy (in J/kg)
 */
/*----------------------------------------------------------------------------*/

void
cs_hgn_thermo_define_stiffened_gas(int        iph,
                                   cs_real_t  cv,
                                   cs_real_t  gamma,
                                   cs_real_t  pinf,
                                   cs_real_t  qprim,
                                   cs_real_t  q)
{
  if (iph > 1) {
    bft_error(__FILE__, __LINE__, 0,
              "Error while defining a stiffened gas with homogeneous two-phase"
              " flow model,\n two phases allowed.");
  }

  _stiffened_gas[iph].cv    = cv;
  _stiffened_gas[iph].gamma = gamma;
  _stiffened_gas[iph].pinf  = pinf;
  _stiffened_gas[iph].qprim = qprim;
  _stiffened_gas[iph].q     = q;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the entropy (in J/kg/K) of phase iph (0 or 1).
 *
 *        This function returns the entropy of phase iph with respect to the
 *        specific volume of phase iph and to the specific energy of phase iph.
 *
 * \param[in]     vol   specific volume of phase iph (in m^3/kg)
 * \param[in]     energ specific energy of phase iph (in J/kg)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_entropy_ve(cs_real_t vol,
                               cs_real_t energ,
                               int       iph)
{
  return cs_cf_thermo_entropy_sg_ve(vol, energ, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the temperature (in K) of phase iph.
 *
 *        This function returns the temperature of phase iph with respect to the
 *        specific volume of phase iph and to the specific energy of phase iph.
 *
 * \param[in]     vol   the specific volume of phase iph  (in m^3/kg)
 * \param[in]     energ the specific energy of phase iph (in J/kg)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_temperature_ve(cs_real_t vol,
                                   cs_real_t energ,
                                   int       iph)
{
  return cs_cf_thermo_temperature_sg_ve(vol, energ, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the pressure (in Pa) of phase iph.
 *
 *        This function returns the pressure of phase iph with respect to the
 *        specific volume of phase iph and to the specific energy of phase iph.
 *
 * \param[in]     vol   the specific volume of phase iph (in m^3/kg)
 * \param[in]     energ the specific energy of phase iph (in J/kg)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_pressure_ve(cs_real_t vol,
                                cs_real_t energ,
                                int       iph)
{
  return cs_cf_thermo_pressure_sg_ve(vol, energ, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of internal energy (in J/kg) of phase iph in plane (T,P).
 *
 *        This function returns the internal energy of phase iph with respect
 *        to the temperature of phase iph and to the pressure energy of phase iph.
 *
 * \param[in]     T    temperature (in K)
 * \param[in]     P    pressure    (in Pa)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_internal_energy_tp(cs_real_t T,
                                       cs_real_t P,
                                       int       iph)
{
  return cs_cf_thermo_internal_energy_sg_tp(T, P, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of specific volume (in m^3/kg) of phase iph in plane (T,P).
 *
 *        This function returns the specific volume of phase iph with respect
 *        to the temperature of phase iph and to the pressure energy of phase iph.
 *
 * \param[in]     T    temperature (in K)
 * \param[in]     P    pressure    (in Pa)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_specific_volume_tp(cs_real_t T,
                                       cs_real_t P,
                                       int       iph)
{
  return cs_cf_thermo_specific_volume_sg_tp(T, P, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of specific entropy (in J/kg/K) of phase iph in plane (T,P).
 *
 *        This function returns the specific entropy of phase iph with respect
 *        to the temperature of phase iph and to the pressure energy of phase iph.
 *
 * \param[in]     T    temperature (in K)
 * \param[in]     P    pressure    (in Pa)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_entropy_tp(cs_real_t T,
                               cs_real_t P,
                               int       iph)
{
 return cs_cf_thermo_entropy_sg_tp(T, P, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of internal energy (in J/kg) of phase iph in plane (s,v).
 *
 *        This function returns the internal energy of phase iph with respect
 *        to the specific entropy of phase iph and to the specific volume of
 *        phase iph.
 *
 * \param[in]     s    specific entropy          (in J/kg/K)
 * \param[in]     v    specific volume  (in m^3/kg)
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_internal_energy_sv(cs_real_t s,
                                       cs_real_t v,
                                       int       iph)
{
  return cs_cf_thermo_internal_energy_sg_sv(s, v, _stiffened_gas[iph]);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
