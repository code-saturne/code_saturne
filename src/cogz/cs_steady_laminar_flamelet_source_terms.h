#ifndef __CS_STEADY_LAMINAR_FLAMELET_SOURCE_TERMS__
#define __CS_STEADY_LAMINAR_FLAMELET_SOURCE_TERMS__

/*============================================================================
 * Defines the source terms for the soot mass fraction and the precursor
 * number for soot model of Moss et al for one time step.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_steady_laminar_flamelet_source_terms.c
 *
 * \brief Specific physic routine: STE/VTE and progress variable equations.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \param[in]      fld_id        field id
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_steady_laminar_flamelet_source_terms(int        fld_id,
                                        cs_real_t  smbrs[],
                                        cs_real_t  rovsdt[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_STEADY_LAMINAR_FLAMELET_SOURCE_TERMS__ */
