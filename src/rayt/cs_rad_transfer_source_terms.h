#ifndef __CS_RAD_TRANSFER_SOURCE_TERMS_H__
#define __CS_RAD_TRANSFER_SOURCE_TERMS_H__

/*============================================================================
 * Radiation solver operations.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

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
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Implicit and explicit radiative source terms for thermal scalar.
 *
 * \param[in,out]  smbrs   work array for right hand side
 * \param[in,out]  rovsdt  work array for unsteady term
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_source_terms(cs_real_t  smbrs[],
                             cs_real_t  rovsdt[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_SOURCE_TERMS_H__ */
