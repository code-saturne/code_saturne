#ifndef __CS_PARAMETERS_H__
#define __CS_PARAMETERS_H__

/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Structure of variable calculation options
 *----------------------------------------------------------------------------*/

typedef struct {
  int     iwarni;
  int     iconv;
  int     istat;
  int     idiff;
  int     idifft;
  int     idften;
  int     iswdyn;
  int     ischcv;
  int     ibdtso;
  int     isstpc;
  int     nswrgr;
  int     nswrsm;
  int     imrgra;
  int     imligr;
  int     ircflu;
  int     iwgrec;       /* gradient calculation
                           - 0: standard (default)
                           - 1: weighted (could be used with imvisf = 1) */
  double  thetav;
  double  blencv;
  double  epsilo;
  double  epsrsm;
  double  epsrgr;
  double  climgr;
  double  extrag;
  double  relaxv;
} cs_var_cal_opt_t;

/*----------------------------------------------------------------------------
 * Structure of the solving info
 *----------------------------------------------------------------------------*/

typedef struct {
  int     n_it;
  double  rhs_norm;
  double  res_norm;
  double  derive;
  double  l2residual;
} cs_solving_info_t;

/*----------------------------------------------------------------------------
 * Structure of condensation modelling physical properties
 *----------------------------------------------------------------------------*/

typedef struct {
  double  mol_mas;
  double  cp;
  double  vol_dif;
  double  mu_a;
  double  mu_b;
  double  lambda_a;
  double  lambda_b;
} cs_gas_mix_species_prop_t;


/*----------------------------------------------------------------------------
 * Boundary condition types
 *----------------------------------------------------------------------------*/

enum {
  CS_INDEF = 1,
  CS_INLET = 2,
  CS_OUTLET = 3,
  CS_SYMMETRY = 4,
  CS_SMOOTHWALL = 5,
  CS_ROUGHWALL = 6,
  CS_ESICF = 7,
  CS_SSPCF = 8,
  CS_SOPCF = 9,
  CS_ERUCF = 10,
  CS_EPHCF = 11,
  CS_EQHCF = 12,
  CS_COUPLED = 13,
  CS_FREE_INLET = 14,
  CS_CONVECTIVE_INLET = 15
};

/*----------------------------------------------------------------------------
 * Space discretisation options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  int           imvisf;       /* face viscosity field interpolation
                                 - 1: harmonic
                                 - 0: arithmetic (default) */

  int           imrgra;       /* type of gradient reconstruction
                                 - 0: iterative process
                                 - 1: standard least square method
                                 - 2: least square method with extended
                                      neighbourhood
                                 - 3: least square method with reduced extended
                                      neighbourhood
                                 - 4: iterative precess initialized by the least
                                      square method */

  double        anomax;       /* non orthogonality angle of the faces, in radians.
                                 For larger angle values, cells with one node
                                 on the wall are kept in the extended support of
                                 the neighbouring cells. */

  int           iflxmw;       /* method to compute interior mass flux due to ALE
                                 mesh velocity
                                 - 1: based on cell center mesh velocity
                                 - 0: based on nodes displacement */

} cs_space_disc_t;

/*----------------------------------------------------------------------------
 * PISO descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  int           nterup;       /* number of interations on the pressure-velocity
                                 coupling on Navier-Stokes */

  double        epsup;        /* relative precision for the convergence test of
                                 the iterative process on pressure-velocity
                                 coupling */

  double        xnrmu;        /* norm  of the increment
                                 \f$ \vect{u}^{k+1} - \vect{u}^k \f$
                                 of the iterative process on pressure-velocity
                                 coupling */

  double        xnrmu0;       /* norm of \f$ \vect{u}^0 \f$ */

} cs_piso_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to space discretisation options structure */

extern const cs_space_disc_t  *cs_glob_space_disc;

/* Pointer to PISO structure */

extern const cs_piso_t  *cs_glob_piso;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define general field keys.
 *
 * A recommended practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 *----------------------------------------------------------------------------*/

void
cs_parameters_define_field_keys(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define field key for condensation.
 *
 * Note: this should be moved in the future to a condensation-specific file.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_key_gas_mix(void);

/*----------------------------------------------------------------------------
 * Read general restart info.
 *
 * This updates the previous time step info.
 *----------------------------------------------------------------------------*/

void
cs_parameters_read_restart_info(void);

/*----------------------------------------------------------------------------
 * Define a user variable.
 *
 * Solved variables are always defined on cells.
 *
 * parameters:
 *   name <-- name of variable and associated field
 *   dim  <-- variable dimension
 *----------------------------------------------------------------------------*/

void
cs_parameters_add_variable(const char  *name,
                           int          dim);

/*----------------------------------------------------------------------------
 * Define a user variable which is a variance of another variable.
 *
 * Only variances of thermal or user-defined variables are currently handled.
 *
 * parameters:
 *   name          <-- name of variance and associated field
 *   variable_name <-- name of associated variable
 *----------------------------------------------------------------------------*/

void
cs_parameters_add_variable_variance(const char  *name,
                                    const char  *variable_name);

/*----------------------------------------------------------------------------
 * Define a user property.
 *
 * parameters:
 *   name        <-- name of property and associated field
 *   dim         <-- property dimension
 *   location_id <-- id of associated mesh location
 *----------------------------------------------------------------------------*/

void
cs_parameters_add_property(const char  *name,
                           int          dim,
                           int          location_id);

/*----------------------------------------------------------------------------
 * Return the number of defined user variables not added yet.
 *
 * This number is reset to 0 when cs_parameters_create_added_variables()
 * is called.
 *
 * returns:
 *   number of defined user variables
 *----------------------------------------------------------------------------*/

int
cs_parameters_n_added_variables(void);

/*----------------------------------------------------------------------------
 * Return the number of defined user properties not added yet.
 *
 * This number is reset to 0 when cs_parameters_create_added_properties()
 * is called.
 *
 * returns:
 *   number of defined user properties
 *----------------------------------------------------------------------------*/

int
cs_parameters_n_added_properties(void);

/*----------------------------------------------------------------------------
 * Create previously added user variables.
 *----------------------------------------------------------------------------*/

void
cs_parameters_create_added_variables(void);

/*----------------------------------------------------------------------------
 * Create previously added user properties.
 *----------------------------------------------------------------------------*/

void
cs_parameters_create_added_properties(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMETERS_H__ */
