#ifndef __CS_STOKES_MODEL_H__
#define __CS_STOKES_MODEL_H__

/*============================================================================
 * Stokes equation model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Stokes equation model descriptor */
/*---------------------------*/

typedef struct {

  int           ivisse;       /* take viscous term of transposed velocity
                                 gradient into account in momentum equation
                                 - 1: true (default)
                                 - 0: false */
  int           irevmc;       /* reconstruction of the velocity field with the
                                 updated pressure option
                                 - 0: default */
  int           iprco;        /* compute the pressure step thanks to the
                                 continuity equation
                                 - 1: true (default)
                                 - 0: false */
  double        rnormp;       /* normed residual for the pressure step */
  double        arak;         /* Arakawa multiplicator for the Rhie and Chow
                                 filter (1 by default) */
  int           ipucou;       /* pseudo coupled pressure-velocity solver
                                 - 1: true (default)
                                 - 0: false */
  int           iccvfg;       /* calculation with a fixed velocity field
                                 - 1: true (default)
                                 - 0: false */
  int           idilat;       /* algorithm to take into account the density
                                 variation in time
                                 - 1: dilatable steady algorithm (default)
                                 - 2: dilatable unsteady algorithm
                                 - 3: low-Mach algorithm
                                 - 4: algorithm for fire
                                 - 0: boussinesq algorithm with constant
                                 density */
  double        epsdp;        /* parameter of diagonal pressure strengthening */
  int           itbrrb;       /* accurate treatment of the wall temperature
                                 - 1: true
                                 - 0: false (default) */
  int           iphydr;       /* improve hydrostatic pressure algorithm
                                 - 1: impose the equilibrium of the hydrostaic
                                   part of the pressure with any external force,
                                   even head losses
                                 - 2: compute an hydrostatic pressure due to
                                   buoyancy forces before the prediction step
                                 - 0: no treatment (default) */
  int           igprij;       /* improve static pressure algorithm
                                 - 1: take -div(rho R) in the static pressure
                                   treatment IF iphydr=1
                                 - 0: no treatment (default) */
  int           igpust;       /* improve static pressure algorithm
                                 - 1: take user momemtum source terms in the
                                   static pressure treatment IF iphydr=1 (default)
                                 - 0: no treatment */
  int           iifren;       /* indicates the presence of a Bernoulli boundary
                                 face (automatically computed)
                                 - 0: no face
                                 - 1: at least one face */
  int           icalhy;       /* compute the hydrostatic pressure in order to
                                 compute the Dirichlet conditions on the
                                 pressure at outlets
                                 - 1: true
                                 - 0: false (default) */
  int           irecmf;       /* use interpolated face diffusion coefficient
                                 instead of cell diffusion coefficient for the
                                 mass flux reconstruction for the
                                 non-orthogonalities
                                 - 1: true
                                 - 0: false (default) */

} cs_stokes_model_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main Stokes model structure */

extern const cs_stokes_model_t  *cs_glob_stokes_model;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide acces to cs_glob_stokes_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_stokes_model_t *
cs_get_glob_stokes_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the stokes model parameters to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_stokes_model_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_STOKES_MODEL_H__ */
