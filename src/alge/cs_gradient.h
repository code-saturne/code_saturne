#ifndef __CS_GRADIENT_H__
#define __CS_GRADIENT_H__

/*============================================================================
 * Gradient reconstruction.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Gradient types
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_GRADIENT_ITER,        /* Iterative method */
  CS_GRADIENT_LSQ_STD,     /* Least-square method */
  CS_GRADIENT_LSQ_EXT,     /* Least-square method with extended neighborhood  */
  CS_GRADIENT_LSQ_EXT_RED, /* Least-square method with reduced extended neig. */
  CS_GRADIENT_LSQ_ITER,    /* LSQ followed with iterative */
  CS_GRADIENT_N_TYPES

} cs_gradient_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for gradient types */

extern const char *cs_gradient_type_name[];

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdcel, CGDCEL)
(
 const cs_int_t   *const ivar,        /* <-- variable number                  */
 const cs_int_t   *const imrgra,      /* <-- gradient computation mode        */
 const cs_int_t   *const inc,         /* <-- 0 or 1: increment or not         */
 const cs_int_t   *const iccocg,      /* <-- 1 or 0: recompute COCG or not    */
 const cs_int_t   *const imobil,      /* <-- 1 for mobile mesh, 0 otherwise   */
 const cs_int_t   *const iale,        /* <-- 1 for ALE, 0 otherwise           */
 const cs_int_t   *const nswrgp,      /* <-- >1: with reconstruction          */
 const cs_int_t   *const idimtr,      /* <-- 0, 1, 2: scalar, vector, tensor
                                             in case of rotation              */
 const cs_int_t   *const iphydp,      /* <-- use hydrosatatic pressure        */
 const cs_int_t   *const iwarnp,      /* <-- verbosity level                  */
 const cs_int_t   *const imligp,      /* <-- type of clipping                 */
 const cs_real_t  *const epsrgp,      /* <-- precision for iterative gradient
                                             calculation                      */
 const cs_real_t  *const extrap,      /* <-- extrapolate gradient at boundary */
 const cs_real_t  *const climgp,      /* <-- clipping coefficient             */
 const cs_int_t          isympa[],    /* <-- indicator for symmetry faces     */
       cs_real_t         fextx[],     /* <-- components of the exterior force */
       cs_real_t         fexty[],     /*     generating the hydrostatic       */
       cs_real_t         fextz[],     /*     pressure                         */
 const cs_real_t         coefap[],    /* <-- boundary condition term          */
 const cs_real_t         coefbp[],    /* <-- boundary condition term          */
       cs_real_t         pvar[],      /* <-- gradient's base variable         */
       cs_real_t         grad[]       /* <-> gradient                         */
);

/*----------------------------------------------------------------------------
 * Compute cell gradient of vector field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdvec, CGDVEC)
(
 const cs_int_t         *const ivar,
 const cs_int_t         *const imrgra,  /* <-- gradient computation mode      */
 const cs_int_t         *const inc,     /* <-- 0 or 1: increment or not       */
 const cs_int_t         *const nswrgp,  /* <-- >1: with reconstruction        */
 const cs_int_t         *const iwarnp,  /* <-- verbosity level                */
 const cs_int_t         *const imligp,  /* <-- type of clipping               */
 const cs_real_t        *const epsrgp,  /* <-- precision for iterative
                                               gradient calculation           */
 const cs_real_t        *const climgp,  /* <-- clipping coefficient           */
 const cs_real_3_t  (*restrict coefav), /* <-- boundary condition term        */
 const cs_real_33_t (*restrict coefbv), /* <-- boundary condition term        */
 const cs_real_3_t  (*restrict pvar),   /* <-- gradient's base variable       */
       cs_real_33_t (*restrict gradv)   /* <-> gradient of the variable       */
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GRADIENT__ */
