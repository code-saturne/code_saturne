#ifndef __CS_CDO_H__
#define __CS_CDO_H__

/*============================================================================
 * General functions or variables for the INNOV module
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



/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_CDO_LEN_NAME 64

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef unsigned short int cs_flag_t;

/* Type of numerical scheme for the discretization in space */
typedef enum {

  CS_SPACE_SCHEME_CDOVB,   /* CDO scheme with vertex-based positionning */
  CS_SPACE_SCHEME_CDOFB,   /* CDO cell-based scheme with hybridization */
  CS_SPACE_N_SCHEMES

} cs_space_scheme_t;

/* Vector-valued quantity stored using its measure (i.e. length) and
   its direction given by a unitary vector */
typedef struct {

  double  meas;
  double  unitv[3];

} cs_nvec3_t;

/* Values associated to the different ways to retrieve data */
typedef union {

  cs_flag_t           flag;       // flag
  int                 id;         // identification number
  cs_lnum_t           num;        // local number
  cs_real_t           val;        // value
  cs_real_2_t         couple;     // two values
  cs_real_3_t         vect;       // vector: 3 values
  cs_nvec3_t          nvec3;      // meas + unit vector
  cs_real_6_t         twovects;   // two vectors
  cs_real_33_t        tens;       // tensor: 9 values

} cs_get_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic analytic function
 *
 * \param[in]      time       when ?
 * \param[in]      xyz        where ?
 * \param[in, out] retval     result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_analytic_func_t) (cs_real_t           time,
                      const cs_real_3_t   xyz,
                      cs_get_t           *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Simple function to define the time step according to the number of
 *         iteration already done
 *
 * \param[in]      time_iter  current number of iterations
 *
 * \return the value of the time step
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_timestep_func_t) (int    time_iter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a quantity according to a law depending only
 *         on one variable.
 *         This law is described by a set of parameters stored in a structure.
 *         result = law(var_value)
 *
 * \param[in]      var_value  value of the variable attached to this law
 * \param[in]      law_param  set of paramters related to the current law
 * \param[in, out] retval     result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_onevar_law_func_t) (double         var_value,
                        const void    *law_param,
                        cs_get_t      *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a quantity according to a law depending only
 *         on two variables.
 *         This law is described by a set of parameters stored in a structure.
 *         result = law(var1_value, var2_value)
 *
 * \param[in]      var1_value  value of the first variable attached to this law
 * \param[in]      var2_value  value of the second variable attached to this law
 * \param[in]      law_param   set of paramters related to the current law
 * \param[in, out] retval      result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_twovar_law_func_t) (double         var1_value,
                        double         var2_value,
                        const void    *law_param,
                        cs_get_t      *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a quantity according to a law depending only
 *         on two variables (the first one is a scalar and the second one a
 *         vector)
 *         This law is described by a set of parameters stored in a structure.
 *         result = law(var1_value, var2_value)
 *
 * \param[in]      var1_value  value of the first variable attached to this law
 * \param[in]      var2_value  value of the second variable attached to this law
 * \param[in]      law_param   set of paramters related to the current law
 * \param[in, out] retval      result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_scavec_law_func_t) (double          var1_value,
                        const double    var2_vect[],
                        const void     *law_param,
                        cs_get_t       *retval);

/*============================================================================
 * Global variables
 *============================================================================*/

/* Separation lines: long, medium, short */
extern const char lsepline[];
extern const char msepline[];
extern const char ssepline[];

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string "true" or "false" according to the boolean
 *
 * \param[in]  boolean     bool  type
 *
 * \return a string "true" or "false"
 */
/*----------------------------------------------------------------------------*/

const char *
cs_base_strtf(bool  boolean);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute epsilon which is the machine precision
 */
/*----------------------------------------------------------------------------*/

void
cs_set_eps_machine(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the machine precision
 */
/*----------------------------------------------------------------------------*/

double
cs_get_eps_machine(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the threshold under which one considers it's zero
 */
/*----------------------------------------------------------------------------*/

double
cs_get_zero_threshold(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_nvec3_t structure from a cs_real_3_t
 *
 * \param[in]  v     vector of size 3
 * \param[out] qv    pointer to a cs_nvec3_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_nvec3(const cs_real_3_t    v,
         cs_nvec3_t          *qv);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_H__ */
