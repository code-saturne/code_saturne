#ifndef __CS_HGN_PHASE_THERMO_H__
#define __CS_HGN_PHASE_THERMO_H__

/*============================================================================
 * Phase thermodynamic for the compressible homogeneous two-phase model
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define stiffened gas parameters for phase iph.
 *
 * \param[in]  iph        index of phase (0 or 1)
 * \param[in]  cv         heat capacity
 * \param[in]  gamma      polytropic coefficient
 * \param[in]  pinf       minimum pressure
 * \param[in]  qprim      entropy parameter
 * \param[in]  q
 */
/*----------------------------------------------------------------------------*/

void
cs_hgn_thermo_define_stiffened_gas(int        iph,
                                   cs_real_t  cv,
                                   cs_real_t  gamma,
                                   cs_real_t  pinf,
                                   cs_real_t  qprim,
                                   cs_real_t  q);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the entropy of phase iph (0 or 1).
 *
 *        This function returns the entropy of phase iph with respect to the
 *        specific volume of phase iph and to the specific energy of phase iph.
 *
 * \param[in]     vol   specific volume of phase iph
 * \param[in]     energ specific energy of phase iph
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_entropy_ve(cs_real_t vol,
                               cs_real_t energ,
                               int       iph);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the temperature of phase iph.
 *
 *        This function returns the temperature of phase iph with respect to the
 *        specific volume of phase iph and to the specific energy of phase iph.
 *
 * \param[in]     vol   the specific volume of phase iph
 * \param[in]     energ the specific energy of phase iph
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_temperature_ve(cs_real_t vol,
                                   cs_real_t energ,
                                   int       iph);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the pressure of phase iph.
 *
 *        This function returns the pressure of phase iph with respect to the
 *        specific volume of phase iph and to the specific energy of phase iph.
 *
 * \param[in]     vol   the specific volume of phase iph
 * \param[in]     energ the specific energy of phase iph
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_pressure_ve(cs_real_t vol,
                                cs_real_t energ,
                                int       iph);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of internal energy of phase iph in plane (T,P).
 *
 * \param[in]     T    temperature
 * \param[in]     P    pressure
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_internal_energy_tp(cs_real_t T,
                                       cs_real_t P,
                                       int       iph);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of specific volume of phase iph in plane (T,P).
 *
 * \param[in]     T    temperature
 * \param[in]     P    pressure
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_specific_volume_tp(cs_real_t T,
                                       cs_real_t P,
                                       int       iph);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of entropy of phase iph in plane (T,P).
 *
 * \param[in]     T    temperature
 * \param[in]     P    pressure
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_entropy_tp(cs_real_t T,
                               cs_real_t P,
                               int       iph);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of internal energy of phase iph in plane (s,v).
 *
 * \param[in]     s    entropy
 * \param[in]     v    specific volume
 * \param[in]     iph   index of phase
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_phase_thermo_internal_energy_sv(cs_real_t s,
                                       cs_real_t v,
                                       int       iph);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif  /* __CS_HGN_PHASE_THERMO_H__ */
