#ifndef __CS_PARAMETERS_H__
#define __CS_PARAMETERS_H__

/*============================================================================
 * General parameters management.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
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
 *----------------------------------------------------------------------------*/

typedef struct {
  int     iwarni;
  int     iconv;
  int     istat;
  int     idircl;
  int     ndircl;
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
  int     icoupl;       /* internal coupling
                           - -1: not coupled (default)
                           -  1: coupled                                 */

  double  thetav;
  double  blencv;
  double  blend_st;
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
  double  muref;     /* ref. viscosity for Sutherland law                */
  double  lamref;    /* ref. thermal conductivity for Sutherland law     */
  double  trefmu;    /* ref. temperature for viscosity in Sutherland law */
  double  treflam;   /* ref. temperature for conductivity Sutherland law */
  double  smu;       /* Sutherland temperature for viscosity             */
  double  slam;      /* Sutherland temperature for conductivity          */
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
  CS_DRIFT_SCALAR_ZERO_BNDY_FLUX = (1 << 7)
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
                                 - 4: iterative process initialized by the least
                                      square method */

  double        anomax;       /* non orthogonality angle of the faces, in radians.
                                 For larger angle values, cells with one node
                                 on the wall are kept in the extended support of
                                 the neighboring cells. */

  int           iflxmw;       /* method to compute interior mass flux due to ALE
                                 mesh velocity
                                 - 1: based on cell center mesh velocity
                                 - 0: based on nodes displacement */

} cs_space_disc_t;

/*----------------------------------------------------------------------------
 * PISO descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  int           nterup;         /* number of interations on the pressure-velocity
                                   coupling on Navier-Stokes */

  double        epsup;          /* relative precision for the convergence test of
                                   the iterative process on pressure-velocity
                                   coupling */

  double        xnrmu;          /* norm  of the increment
                                   \f$ \vect{u}^{k+1} - \vect{u}^k \f$
                                   of the iterative process on pressure-velocity
                                   coupling */

  double        xnrmu0;         /* norm of \f$ \vect{u}^0 \f$ */

  int           n_buoyant_scal; /* number of buoyant scalars,
                                   zero if there is no buoyant scalar */

} cs_piso_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to space discretisation options structure */

extern const cs_space_disc_t  *cs_glob_space_disc;

/* Pointer to PISO options structure */

extern const cs_piso_t        *cs_glob_piso;

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
 * \return   space discretization description structure
 */
/*----------------------------------------------------------------------------*/

cs_space_disc_t *
cs_get_glob_space_disc(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_piso
 *
 * needed to initialize structure with GUI and user C functions.
 *
 * \return   piso information structure
 */
/*----------------------------------------------------------------------------*/

cs_piso_t *
cs_get_glob_piso(void);

/*----------------------------------------------------------------------------
 *!
 * \brief Count and set number of buoyant scalars.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_set_n_buoyant_scalars(void);

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
 * \brief Define field key for condensation.
 *
 * Note: this should be moved in the future to a condensation-specific file.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_key_gas_mix(void);

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
 * When such a variblae does not exist but we have an Enthalpy variables,
 * an associated temperature boundary field is returned.
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_temperature(void);

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
 * \brief Print the space discretization structure to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_space_disc_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMETERS_H__ */
