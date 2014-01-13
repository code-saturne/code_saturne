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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * EOS library headers
 *----------------------------------------------------------------------------*/

#include <EOS/API/EOS.hxx>
#include <EOS/API/EOS_Field.hxx>
#include <EOS/API/EOS_Fields.hxx>
#include <EOS/API/EOS_Error_Field.hxx>
#include <EOS/API/EOS_enums.hxx>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_eos.hxx"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 *  Global variables
 *============================================================================*/

/* Pointer on a EOS instances array */

NEPTUNE::EOS *eos;

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*============================================================================
 * Public function definitions
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
              char *EOSRef)
{
    eos = new NEPTUNE::EOS(EOSMethod, EOSRef);
}

/*----------------------------------------------------------------------------
 * Delete EOS instances
 *----------------------------------------------------------------------------*/

void
cs_eos_destroy(void)
{
    delete eos;
}

/*----------------------------------------------------------------------------*/
