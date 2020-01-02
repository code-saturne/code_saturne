#ifndef __CS_GUI_VARIABLES_H__
#define __CS_GUI_VARIABLES_H__

/*============================================================================
 * Management of the GUI parameters file: variables
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Variables and scalars management structure
 *----------------------------------------------------------------------------*/

typedef struct {
  char  *model;            /* predifined physics model                        */
  char  *model_value;      /* predifined physics model value                  */
} cs_var_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_var_t    *cs_glob_var;   /* Pointer to main variables structure */

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_VARIABLES_H__ */
