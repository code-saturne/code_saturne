#ifndef __CS_PARAMETERS_H__
#define __CS_PARAMETERS_H__

/*============================================================================
 * General parameters management.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_tree.h"

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
 * (now an alias of cs_equation_param_t)
 *----------------------------------------------------------------------------*/

typedef cs_equation_param_t cs_var_cal_opt_t;

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
  CS_EPHCF = 10,
  CS_EQHCF = 11,
  CS_COUPLED = 12,           /* coupled face */
  CS_COUPLED_FD = 13,        /* coupled face with decentered flux */
  CS_FREE_INLET = 14,
  CS_FREE_SURFACE = 15,
  CS_CONVECTIVE_INLET = 16
};

/*----------------------------------------------------------------------------
 * flag for computing the drift mass flux:
 * (for coal classes for instance, only the first
 *  scalar of a class compute the drift flux of the class
 *  and the other scalars use it without recomputing it)
 *----------------------------------------------------------------------------*/

enum {
  CS_DRIFT_SCALAR_ON = (1 << 0),
  CS_DRIFT_SCALAR_ADD_DRIFT_FLUX = (1 << 1),
  CS_DRIFT_SCALAR_THERMOPHORESIS = (1 << 2),
  CS_DRIFT_SCALAR_TURBOPHORESIS = (1 << 3),
  CS_DRIFT_SCALAR_ELECTROPHORESIS = (1 << 4),
  CS_DRIFT_SCALAR_CENTRIFUGALFORCE = (1 << 5),
  CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX = (1 << 6),
  CS_DRIFT_SCALAR_ZERO_BNDY_FLUX = (1 << 7),
  CS_DRIFT_SCALAR_ZERO_BNDY_FLUX_AT_WALLS = (1 << 8)
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
                                      neighborhood
                                 - 3: least square method with reduced extended
                                      neighborhood
                                 - 4: Green-Gauss using least squares face
                                      values interpolation */

  int           iflxmw;       /* method to compute interior mass flux due to ALE
                                 mesh velocity
                                 - 1: based on cell center mesh velocity
                                 - 0: based on nodes displacement */

  int           itbrrb;       /* accurate treatment of the wall temperature
                                 - 1: true
                                 - 0: false (default) */

} cs_space_disc_t;

/*----------------------------------------------------------------------------
 * Time scheme descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  int           time_order;   /* Global time order of the time stepping */

  int           isto2t;       /* time scheme activated for the source
                                 terms of turbulent equations */

  double        thetst;       /* value of \f$theta\f$ for turbulence */

  int           iccvfg;       /* calculation with a fixed velocity field
                                 - 1: true (default)
                                 - 0: false */
} cs_time_scheme_t;

/*----------------------------------------------------------------------------
 * Auxiliary checkpoint/restart file parameters
 *----------------------------------------------------------------------------*/

typedef struct {

  int     read_auxiliary;   /* Activate reading of auxiliary restart file */
  int     write_auxiliary;  /* Activate output of auxiliary restart file */

} cs_restart_auxiliary_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to space discretisation options structure */

extern const cs_space_disc_t  *cs_glob_space_disc;

/* Pointer to time scheme  options structure */

extern const cs_time_scheme_t  *cs_glob_time_scheme;

/* Pointer to auxiliary checkpoint/restart file parameters */

extern cs_restart_auxiliary_t  *cs_glob_restart_auxiliary;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Global parameters tree structure */

extern cs_tree_node_t  *cs_glob_tree;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief For a given field, returns the scalar number of the fluctuating field
 * if given field is a variance.
 *
 * \param[in]  f  field
 *
 * \return        if f is a variance: scalar number of fluctuating field
 *                else if f is not a variance: 0
 *                else if f is the variance of a field that is not a scalar: -1
 */
/*----------------------------------------------------------------------------*/

static inline int
cs_parameters_iscavr(cs_field_t *f)
{
  int iscvr = 0, f_id = 0;
  int kscavr = cs_field_key_id("first_moment_id");
  int keysca = cs_field_key_id("scalar_id");

  if (kscavr >= 0) {
    f_id = cs_field_get_key_int(f, kscavr);
    if (f_id >= 0)
      iscvr = cs_field_get_key_int(cs_field_by_id(f_id), keysca);
  }

  return iscvr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_space_disc
 *
 * needed to initialize structure in GUI and user C functions.
 *
 * \return  space discretization description structure
 */
/*----------------------------------------------------------------------------*/

cs_space_disc_t *
cs_get_glob_space_disc(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_time_scheme
 *
 * needed to initialize structure with GUI and user C functions.
 *
 * \return  time scheme information structure
 */
/*----------------------------------------------------------------------------*/

cs_time_scheme_t *
cs_get_glob_time_scheme(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define general field keys.
 *
 * A recommended practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_keys(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read general restart info.
 *
 * This updates the previous time step info.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_read_restart_info(void);

/*----------------------------------------------------------------------------*/
/*!
 * Define a user variable.
 *
 * \brief Solved variables are always defined on cells.
 *
 * \param[in]  name  name of variable and associated field
 * \param[in]  dim   variable dimension
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_variable(const char  *name,
                           int          dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user variable which is a variance of another variable.
 *
 * Only variances of thermal or user-defined variables are currently handled.
 *
 * \param[in]  name           name of variance and associated field
 * \param[in]  variable_name  name of associated variable
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_variable_variance(const char  *name,
                                    const char  *variable_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user property.
 *
 * \param[in]  name         name of property and associated field
 * \param[in]  dim          property dimension
 * \param[in]  location_id  id of associated mesh location
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_property(const char  *name,
                           int          dim,
                           int          location_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined user variables not added yet.
 *
 * This number is reset to 0 when cs_parameters_create_added_variables()
 * is called.
 *
 * \return  number of defined user variables
 */
/*----------------------------------------------------------------------------*/

int
cs_parameters_n_added_variables(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined user properties not added yet.
 *
 * This number is reset to 0 when cs_parameters_create_added_properties()
 * is called.
 *
 * \return   number of defined user properties
 */
/*----------------------------------------------------------------------------*/

int
cs_parameters_n_added_properties(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create previously added user variables.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_create_added_variables(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create previously added user properties.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_create_added_properties(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a boundary values field for a variable field.
 *
 * \param[in]  f  pointer to field structure
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_values(cs_field_t  *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a boundary values field for temperature, if applicable.
 *
 * When a volume temperature variable field already exists, this amounts
 * to calling \ref cs_parameters_add_boundary_values for that field.
 * When such a variable does not exist but we have an Enthalpy variables,
 * an associated temperature boundary field is returned.
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_temperature(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Complete general equation parameter definitions.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_eqp_complete(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Complete general output options definitions.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_output_complete(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a local variable calculation options structure,
 *        with default options.
 *
 * \return  variable calculations options structure
 */
/*----------------------------------------------------------------------------*/

cs_var_cal_opt_t
cs_parameters_var_cal_opt_default(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the time scheme structure to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_scheme_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the space discretization structure to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_space_disc_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMETERS_H__ */
