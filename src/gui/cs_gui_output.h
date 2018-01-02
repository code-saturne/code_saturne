#ifndef __CS_GUI_OUTPUT_H__
#define __CS_GUI_OUTPUT_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Activation of a writer depending of a formula
 *
 * Fortran Interface:
 *
 * subroutine uinpst (ttcabs, ntcabs)
 * *****************
 *
 * integer          uref  <-- reference velocity
 * double          almax  <-- reference length
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinpst, UINPST) (const cs_int_t  *ntcabs,
                                const cs_real_t *ttcabs);


/*----------------------------------------------------------------------------
 * Determine output boundary fields
 *----------------------------------------------------------------------------*/

void CS_PROCF (cspstb, CSPSTB) (cs_int_t        *ipstdv);

/*----------------------------------------------------------------------------
 * Determine output options
 *----------------------------------------------------------------------------*/

void CS_PROCF (csenso, CSENSO) (cs_int_t        *iecaux);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define postprocessing meshes using an XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_postprocess_meshes(void);

/*----------------------------------------------------------------------------
 * Define postprocessing writers using an XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_postprocess_writers(void);

/*-----------------------------------------------------------------------------
 * Post-processing options for fields
 *
 * These options are used for fields not mapped to variables or properties.
 *----------------------------------------------------------------------------*/

void
cs_gui_postprocess_fields(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_OUTPUT_H__ */
