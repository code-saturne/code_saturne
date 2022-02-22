#ifndef __CS_COOLPROP_HXX__
#define __CS_COOLPROP_HXX__

/*============================================================================
 * Equation of state
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

/*----------------------------------------------------------------------------*/
/*!
 * Computes physical properties in (P,h,Yi) for compressible flow.
 *
 * \param[in]   coolprop_material  CoolProp material
 * \param[in]   coolprop_backend   CoolProp backend ("HEOS" by default,
 *                                 "SRK" for cubic, "TTSE&HEOS" or
 *                                 "BICUBIC&HEOS" for tabulated)
 * \param[in]   thermo_plane       type of thermal plane
 * \param[in]   property           type of property to compute
 * \param(in]   n_vals             size of variable and property arrays
 * \param[in]   var1               first variable of thermodynamic plane
 *                                 (pressure)
 * \param[in]   var2               second variable of thermodynamic plane
 * \param[out]  val                computed property values
 */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_phys_prop_coolprop(char                              *coolprop_material,
                      const char                        *coolprop_backend,
                      cs_phys_prop_thermo_plane_type_t   thermo_plane,
                      cs_phys_prop_type_t                property,
                      const cs_lnum_t                    n_vals,
                      const cs_real_t                    var1[],
                      const cs_real_t                    var2[],
                      cs_real_t                          val[]);

/*----------------------------------------------------------------------------*/

#endif /* __CS_COOLPROP_HXX__ */

