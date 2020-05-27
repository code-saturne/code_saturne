#ifndef __CS_ETHEQ_H__
#define __CS_ETHEQ_H__

/*============================================================================
 * Compute etheta and eq variable knowing the saturation.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute etheta and eq variable knowing the saturation.
 *
 * \param[in]       pphy        pressure [Pa]
 * \param[in]       thetal      Liquid potential temperature
 * \param[in]       qw          total water amount
 * \param[in]       qldia       mean liquid water content
 * \param[in]       xnebdia     nebulosity after the diagnostic
 * \param[in]       xnn         second order moment "n" <s'ql'>/(2*sigma_s**2)
 * \param[out]      etheta      sensible heat part of buoyancy flux
 * \param[out]      eq          latent heat part of buoyancy flux
 */
/*----------------------------------------------------------------------------*/

void
cs_etheq(cs_real_t   pphy,
         cs_real_t   thetal,
         cs_real_t   qw,
         cs_real_t   qldia,
         cs_real_t   xnebdia,
         cs_real_t   xnn,
         cs_real_t  *etheta,
         cs_real_t  *eq);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ETHEQ_H__ */
