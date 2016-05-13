#ifndef __CS_LAGR_POST_H__
#define __CS_LAGR_POST_H__

/*============================================================================
 * Lagrangian module postprocessing
 *============================================================================*/

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

#include "assert.h"
#include "cs_base.h"
#include "cs_field.h"

#include "cs_lagr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  /*! \anchor iensi3
    activation (=1) or not (=0) of the recording of the particle/boundary
    interactions in  \ref bound_stat, and of the calculation of the
    statistics at the corresponding boundaries.
    By default, the statistics are unsteady (reset to zero at every
    time step). They may be steady if \ref isttio=1 (i.e.
    calculation of a cumulated value over time, and then calculation of an
    average over time or over the number of interactions with the boundary).*/
  int  iensi3;

  /*! associates (=1) or not (=0) the variable "velocity of the locally
    undisturbed fluid flow field" with the output of particles or
    trajectories. */
  int  ivisv1;

  /*! associates (=1) or not (=0) the variable "particle velocity"
    with the output of particles or trajectories. */
  int  ivisv2;

  /*! associates (=1) or not (=0) the variable "residence time"
    with the output of particles or trajectories. */
  int  ivistp;

  /*! associates (=1) or not (=0) the variable "particle diameter"
    with the output of particles or trajectories. */
  int  ivisdm;

  /*! associates (=1) or not (=0) the variable "particle temperature"
     with the output of particles or trajectories. */
  int  iviste;

  /*! associates (=1) or not (=0) the variable "particle mass"
    with the output of particles or trajectories. */
  int  ivismp;

  /*! associates (=1) or not (=0) the variable "shrinking core diameter of
    the coal particles" with the output of particles or trajectories.
    useful only if \ref iphyla = 2 */
  int  ivisdk;

  /*! associates (=1) or not (=0) the variable "mass of reactive coal of the
    coal particles" with the output of particles or trajectories.
    useful only if \ref iphyla = 2 */
  int  ivisch;

  /*! associates (=1) or not (=0) the variable "mass of coal of the
    coal particles" with the output of particles or trajectories.
    useful only if \ref iphyla = 2 */
  int  ivisck;

  /*! associates (=1) or not (=0) the variable "mass of water of the
    coal particles" with the output of particles or trajectories.
    useful only if \ref iphyla = 2 */
  int  iviswat;

} cs_lagr_post_options_t;

/*============================================================================
  Global variables
  ============================================================================*/

/*! read-only pointer to Lagrangian postprocessing options */

extern const cs_lagr_post_options_t *cs_glob_lagr_post_options;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize Lagrangian postprocessing.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_post_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the global cs_lagr_post_options_t structure.
 *
 * This pointer allows write access to the structure.
 *
 * \return pointer to cs_glob_lagr_post_options
 */
/*----------------------------------------------------------------------------*/

cs_lagr_post_options_t *
cs_lagr_post_get_options(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_POST_H__ */
