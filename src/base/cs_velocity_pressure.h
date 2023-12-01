#ifndef __CS_VELOCITY_PRESSURE_H__
#define __CS_VELOCITY_PRESSURE_H__

/*============================================================================
 * Velocity-pressure coupling model and parameters.
 *============================================================================*/

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

/* Velocity-pressure coupling model descriptor */
/*---------------------------------------------*/

typedef struct {

  int         ivisse;         /* take viscous term of transposed velocity
                                 gradient into account in momentum equation
                                 - 1: true (default)
                                 - 0: false */

  int         idilat;         /* algorithm to take into account the density
                                 variation in time
                                 - 1: dilatable steady algorithm (default)
                                 - 2: dilatable unsteady algorithm
                                 - 3: low-Mach algorithm
                                 - 4: algorithm for fire
                                 - 0: boussinesq algorithm with constant
                                 density */

  bool        fluid_solid;    /* Has a solid zone where dynamics must be killed?
                                 - false (default)
                                 - true */

  int         n_buoyant_scal; /* number of buoyant scalars,
                                 zero if there is no buoyant scalar */

  int         iprcdo;         /* Discretization method for pressure:
                                  - 0: Legacy finite Volume method.
                                  - 1: CDO method, if CDO/FV is coupled. */

} cs_velocity_pressure_model_t;

/*----------------------------------------------------------------------------
 * Velocity-pressure coupling parameters
 *----------------------------------------------------------------------------*/

typedef struct {

  int         iphydr;         /* improve hydrostatic pressure algorithm
                                 - 1: impose the equilibrium of the hydrostaic
                                   part of the pressure with any external force,
                                   even head losses
                                 - 2: compute an hydrostatic pressure due to
                                   buoyancy forces before the prediction step
                                 - 0: no treatment (default) */

  int         icalhy;         /* compute the hydrostatic pressure in order to
                                 compute the Dirichlet conditions on the
                                 pressure at outlets
                                 - 1: true
                                 - 0: false (default) */

  int         iprco;          /* compute the pressure step thanks to the
                                 continuity equation
                                 - 1: true (default)
                                 - 0: false */

  int         ipredfl;        /* deprecated:
                                 switch on mass flux prediction before momentum
                                 solving to be fully conservative in momentum
                                 over time for variable density flows. */

  int         irevmc;         /* reconstruction of the velocity field with the
                                 updated pressure option
                                 - 0: default */

  int         iifren;         /* indicates the presence of a Bernoulli boundary
                                 face (automatically computed)
                                 - 0: no face
                                 - 1: at least one face */
  int         irecmf;         /* use interpolated face diffusion coefficient
                                 instead of cell diffusion coefficient for the
                                 mass flux reconstruction for the
                                 non-orthogonalities
                                 - 1: true
                                 - 0: false (default) */

  int         igprij;         /* improve static pressure algorithm
                                 - 1: take -div(rho R) in the static pressure
                                      treatment IF iphydr=1
                                 - 0: no treatment (default) */

  int         igpust;         /* improve static pressure algorithm
                                 - 1: take user momemtum source terms in the
                                      static pressure treatment IF iphydr=1
                                      (default)
                                 - 0: no treatment */

  int         igrdpp;         /* For the compressible algorithm, indicate whether
                                 the pressure should be updated after solution
                                 of the acoustic equation.
                                 - 1: true (default)
                                 - 0: false */

  int         ipucou;         /* pseudo coupled pressure-velocity solver
                                 - 1: true (default)
                                 - 0: false */

  int         itpcol;         /* time scheme option:
                                  - 0: staggered.
                                  - 1: colocated time scheme. */

  double      arak;           /* Arakawa multiplier for the Rhie and Chow
                                 filter (1 by default) */

  int         rcfact;         /* Indicates the factor of the Rhie and Chow
                                 filter
                                 - 1: dt (default)
                                 - 0: 1/A_u */

  int         staggered;      /* indicates if one works with the 1D staggered
                                 scheme
                                 - 0: colocated (default)
                                 - 1: staggered */

  int         nterup;         /* number of iterations on the pressure-velocity
                                 coupling on Navier-Stokes */

  double      epsup;          /* relative precision for the convergence test of
                                 the iterative process on pressure-velocity
                                 coupling */

  double      xnrmu;          /* norm  of the increment
                                 \f$ \vect{u}^{k+1} - \vect{u}^k \f$
                                 of the iterative process on pressure-velocity
                                 coupling */

  double      xnrmu0;         /* norm of \f$ \vect{u}^0 \f$ */

  double      epsdp;          /* parameter of diagonal pressure strengthening */

} cs_velocity_pressure_param_t;

/* Deprecated structures (partial compatibilty mode) */

typedef cs_velocity_pressure_model_t cs_stokes_model_t;
typedef cs_velocity_pressure_param_t cs_piso_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main velocity pressure coupling model structure */
extern const cs_velocity_pressure_model_t  *cs_glob_velocity_pressure_model;

/* Pointer to main velocity pressure coupling parameters structure */
extern const cs_velocity_pressure_param_t  *cs_glob_velocity_pressure_param;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide read/write access to cs_glob_velocity_pressure_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_velocity_pressure_model_t *
cs_get_glob_velocity_pressure_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_velocity_pressure_param
 *
 * needed to initialize structure with GUI and user C functions.
 *
 * \return  velocity_pressure_param information structure
 */
/*----------------------------------------------------------------------------*/

cs_velocity_pressure_param_t *
cs_get_glob_velocity_pressure_param(void);

/*----------------------------------------------------------------------------
 *!
 * \brief Count and set number of buoyant scalars.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_set_n_buoyant_scalars(void);

/*----------------------------------------------------------------------------
 *!
 * \brief Set `fluid_solid` flag if solid zones are present.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_set_solid(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the stokes model parameters to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_model_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print Velocity-pressure parameters to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_param_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VELOCITY_PRESSURE_H__ */
