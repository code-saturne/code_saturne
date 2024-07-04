#ifndef __CS_WALL_DISTANCE_H__
#define __CS_WALL_DISTANCE_H__

/*============================================================================
 * Compute distance to wall.
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {
  int need_compute;  /* Is the wall distance needs to be computed */

  int is_up_to_date; /* Is the wall distance up to date */

  int method;        /* Computation method for the wall distance */

} cs_wall_distance_options_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main time step structure */

extern const cs_wall_distance_options_t *cs_glob_wall_distance_options;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_wall_distance.c
 *
 * \brief Compute distance to wall.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \param[in]  iterns        iteration number on Navier-Stokes equations
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_distance(int  iterns);

/*----------------------------------------------------------------------------*/
/*
 * \param[in]     visvdr        dynamic viscosity in edge cells after
 *                              driest velocity amortization
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_distance_yplus(cs_real_t  visvdr[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes distance to wall by a brute force geometric approach
 *        (serial only)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_distance_geometric(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide read/write access to cs_glob_wall_distance
 *
 * \return pointer to global wall distance structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_distance_options_t *
cs_get_glob_wall_distance_options(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALL_DISTANCE_H__ */
