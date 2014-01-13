#ifndef __CS_EOS_HXX__
#define __CS_EOS_HXX__

/*============================================================================
 * Equation of state
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * EOS library headers
 *----------------------------------------------------------------------------*/

#include <EOS/API/EOS.hxx>
#include <EOS/API/EOS_Field.hxx>
#include <EOS/API/EOS_Fields.hxx>
#include <EOS/API/EOS_Error_Field.hxx>
#include <EOS/API/EOS_enums.hxx>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

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

void
cs_eos_create(char *EOSMethod,
              char *EOSRef);

/*----------------------------------------------------------------------------
 * Delete EOS instances
 *----------------------------------------------------------------------------*/

void
cs_eos_destroy(void);

#endif /* __CS_EOS_HXX__ */
