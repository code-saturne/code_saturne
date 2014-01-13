#ifndef __CS_WALL_FUNCTIONS_H__
#define __CS_WALL_FUNCTIONS_H__

/*============================================================================
 * Wall functions
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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the friction velocity and y+/u+
 *----------------------------------------------------------------------------*/

void CS_PROCF (wallfunctions, WALLFUNCTIONS)
(
 const cs_int_t   *const iwallf,      /* <-- wall function type               */
 const cs_int_t   *const ifac,        /* <-- face number                      */
 const cs_real_t  *const xkappa,      /* <-- Von Karman constant              */
 const cs_real_t  *const cstlog,      /* <-- Log law constant                 */
 const cs_real_t  *const cmu025,      /* <-- C_mu^1/4                         */
 const cs_real_t  *const ypluli,      /* <-- y+ limit                         */
 const cs_real_t  *const apow,        /* <-- Coef of the Wener law            */
 const cs_real_t  *const bpow,        /* <-- Coef of the Wener law            */
 const cs_real_t  *const cpow,        /* <-- Coef of the Wener law            */
 const cs_real_t  *const viscosity,   /* <-- kinematic viscosity              */
 const cs_real_t  *const t_visc,      /* <-- turbulent kinematic viscosity    */
 const cs_real_t  *const vel,         /* <-- wall projected
                                             cell center velocity             */
 const cs_real_t  *const y,           /* <-- wall distance                    */
 const cs_real_t  *const kinetic_en,  /* <-- turbulente kinetic energy        */
       cs_int_t         *iuntur,      /* <-- indicator:
                                              0 in the viscous sublayer       */
       cs_int_t         *nsubla,      /* <-- counter of cell in the viscous
                                             sublayer                         */
       cs_int_t         *nlogla,      /* <-- counter of cell in the log-layer */
       cs_real_t        *ustar,       /* --> friction velocity                */
       cs_real_t        *uk,          /* --> friction velocity                */
       cs_real_t        *yplus,       /* --> non-dimension wall distance      */
       cs_real_t        *ypup,        /* --> yplus projected vel ratio        */
       cs_real_t        *cofimp,      /* --> |U_F|/|U_I^p| to ensure a good
                                             turbulence production            */
       cs_real_t        *dplus        /* --> dimensionless shift to the wall
                                             for scalable wall functions      */
);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALL_FUNCTIONS_H__ */
