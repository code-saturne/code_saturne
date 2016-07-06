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
#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Macros used for building a CDO system */
#define CDO_DIFFUSION   0
#define CDO_ADVECTION   1
#define CDO_REACTION    2
#define CDO_TIME        3
#define CDO_SOURCETERM  4
#define N_CDO_TERMS     5

/* Flag related to CDO system building */
#define CS_CDO_BUILD_HCONF      (1 <<  0)  // 1: build conforming Hodge op.
#define CS_CDO_BUILD_LOC_HCONF  (1 <<  1)  // 2: build the local conf. Hodge
#define CS_CDO_PRIMAL_SOURCE    (1 <<  2)  // 4: primal source term requested
#define CS_CDO_DUAL_SOURCE      (1 <<  3)  // 8: dual source term requested

/* Flags use to identify the nature/status of an object (variable, property) */
#define CS_FLAG_STATE_UNIFORM     (1 << 0) //    1: uniform (in space)
#define CS_FLAG_STATE_CELLWISE    (1 << 1) //    2: cellwise uniform
#define CS_FLAG_STATE_UNSTEADY    (1 << 2) //    4: unsteady
#define CS_FLAG_STATE_POTENTIAL   (1 << 3) //    8: potential
#define CS_FLAG_STATE_CIRCULATION (1 << 4) //   16: circulation
#define CS_FLAG_STATE_FLUX        (1 << 5) //   32: flux
#define CS_FLAG_STATE_DENSITY     (1 << 6) //   64: density
#define CS_FLAG_STATE_OWNER       (1 << 7) //  128: owner

/* Flags use to identify where is located a variable and how to access to
   its values */
#define CS_FLAG_PRIMAL       (1 <<  0) //    1: on primal mesh
#define CS_FLAG_DUAL         (1 <<  1) //    2: on dual mesh
#define CS_FLAG_VERTEX       (1 <<  2) //    4: on vertices
#define CS_FLAG_EDGE         (1 <<  3) //    8: on edges
#define CS_FLAG_FACE         (1 <<  4) //   16: on faces
#define CS_FLAG_CELL         (1 <<  5) //   32: on cells
#define CS_FLAG_BORDER       (1 <<  6) //   64: located on the boundary
#define CS_FLAG_SCAL         (1 <<  7) //  128: scalar-valued (stride = 1)
#define CS_FLAG_VECT         (1 <<  8) //  256: vector-valued (stride = 3)
#define CS_FLAG_TENS         (1 <<  9) //  512: tensor-valued (stride = 9)
#define CS_FLAG_SCAN_BY_CELL (1 << 10) // 1024: by cell (c2e, c2f, c2v)

/* Flags use to identify the type of numerical schemes requested for computing
   the different equations attached to the computational domain. If flag is
   activated, then at least one equation solved is discretized using thiq kind
   of numerical scheme. */
#define CS_SCHEME_FLAG_CDOVB     (1 << 0) //  1: CDO vertex-based scheme
#define CS_SCHEME_FLAG_CDOVCB    (1 << 1) //  2: CDO vertex+cell-based scheme
#define CS_SCHEME_FLAG_CDOFB     (1 << 2) //  4: CDO face-based scheme
#define CS_SCHEME_FLAG_HHO       (1 << 3) //  8: Hybrid-High Order scheme
#define CS_SCHEME_FLAG_SCALAR    (1 << 4) // 16: scheme for scalar eq.
#define CS_SCHEME_FLAG_VECTOR    (1 << 5) // 32: scheme for a vector eq.

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of numerical scheme for the discretization in space */
typedef enum {

  CS_SPACE_SCHEME_CDOVB,   /* CDO scheme with vertex-based positionning */
  CS_SPACE_SCHEME_CDOVCB,  /* CDO scheme with vertex+cell-based positionning */
  CS_SPACE_SCHEME_CDOFB,   /* CDO cell-based scheme with hybridization */
  CS_SPACE_SCHEME_HHO,     /* Hybrid High Order scheme (CDO-FB + high-order) */
  CS_SPACE_N_SCHEMES

} cs_space_scheme_t;

/* Description of an object (property, advection field, array..) using
   mask of bits (i.e flag) */

typedef struct {

  cs_flag_t  location;  /* where is defined this object */
  cs_flag_t  state;     /* nature and additional information on this object */

} cs_desc_t;

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
 * \param[in]   time_iter   current number of iterations
 * \param[in]   time        value of the time at the end of the last iteration
 *
 * \return the value of the time step
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_timestep_func_t) (int      time_iter,
                      double   time);

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

/* Default locations */
extern const cs_flag_t  cs_cdo_primal_vtx;
extern const cs_flag_t  cs_cdo_primal_face;
extern const cs_flag_t  cs_cdo_primal_cell;
extern const cs_flag_t  cs_cdo_dual_vtx;
extern const cs_flag_t  cs_cdo_dual_face;
extern const cs_flag_t  cs_cdo_dual_cell;
extern const cs_flag_t  cs_cdo_dual_face_byc;

/*============================================================================
 * Static inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a location matches a referenced support
 *
 * \param[in]  location      flag corresponding to the location to check
 * \param[in]  reference     flag corresponding to the referenced support
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_cdo_same_support(cs_flag_t   location,
                    cs_flag_t   reference)
{
  if ((location & reference) == reference)
    return true;
  else
    return false;
}

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
