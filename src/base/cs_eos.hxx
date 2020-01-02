#ifndef __CS_EOS_HXX__
#define __CS_EOS_HXX__

/*============================================================================
 * Equation of state
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_physical_properties.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new EOS instances
 *
 * parameters:
 *   EOSMethod  <-- table for EOS thermo properties (CATHARE, THETIS, ...)
 *   EOSRef     <-- reference table for EOS thermo properties
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_eos_create(const char *EOSMethod,
              char       *EOSRef);

/*----------------------------------------------------------------------------
 * Delete EOS instances
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_eos_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * Computes physical properties in (P,h,Yi) for compressible flow.
 *
 * parameters:
 * thermo_plane <--  type of thermal plane
 * property     <--  type of property to compute
 * n_vals       <--  size of properties arrays
 * var1         <--  array of pressure
 * var2         <--  array of thermal properties
 * val          -->  array of property
 */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_phys_prop_eos(cs_phys_prop_thermo_plane_type_t   thermo_plane,
                 cs_phys_prop_type_t                property,
                 const cs_lnum_t                    n_vals,
                 double                             var1[],
                 double                             var2[],
                 cs_real_t                          val[]);

/*----------------------------------------------------------------------------*/

#endif /* __CS_EOS_HXX__ */
