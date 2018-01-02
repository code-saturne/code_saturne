#ifndef __CS_LAGR_LAGNEW_H__
#define __CS_LAGR_LAGNEW_H__

/*============================================================================
 * Handling of new particles.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inject a series of particles at random positions on a given zone.
 *
 * \warning Currently works only for tri and quadrangular faces.
 *
 * The fluid velocity and other variables and attributes are computed here.
 *
 * \param[inout]   npt      current number of particles
 * \param[in]      nznew    number of added particles
 * \param[in]      zone_id  zone number (zone index + 1 in arrays)
 * \param[in]      ifrlag   boundary zone number for lagrangian
 * \param[in,out]  iworkp   array containing injection face number
 *----------------------------------------------------------------------------*/

void
cs_lagr_new(cs_lnum_t  *npt,
            cs_lnum_t   nznew,
            int         zone_id,
            const int   ifrlag[],
            cs_lnum_t   iworkp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialization for new particles.
 *
 * The fluid velocity seen is computed here.
 *
 * \param[in]   p_id_l     lower particle id bound (included)
 * \param[in]   p_id_u     uppder particle id bound (excluded)
 * \param[in]   time_id    associated time id (0: current, 1: previous)
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new_particle_init(cs_lnum_t   p_id_l,
                          cs_lnum_t   p_id_u,
                          cs_lnum_t   time_id,
                          cs_real_t   vislen[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_LAGNEW_H__ */
