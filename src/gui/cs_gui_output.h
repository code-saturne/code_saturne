#ifndef __CS_GUI_OUTPUT_H__
#define __CS_GUI_OUTPUT_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void cs_gui_postprocess_writers(void);
void cs_gui_postprocess_meshes(void);

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Activation of a writer depending of a formula
 *
 * Fortran Interface:
 *
 * SUBROUTINE uinpst (ttcabs, ntcabs)
 * *****************
 *
 * INTEGER          UREF   <--   reference velocity
 * DOUBLE          ALMAX  <--   reference length
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinpst, UINPST) (const cs_int_t  *ntcabs,
                                const cs_real_t *ttcabs);


/*----------------------------------------------------------------------------
 * Input/output treatment
 *----------------------------------------------------------------------------*/

void CS_PROCF (csenso, CSENSO) (const cs_int_t  *nvppmx,
                                cs_int_t        *ncapt,
                                cs_int_t        *nthist,
                                cs_real_t       *frhist,
                                cs_int_t        *ntlist,
                                cs_int_t        *iecaux,
                                cs_int_t        *ipstdv,
                                cs_int_t        *ichrvr,
                                cs_int_t        *ilisvr,
                                cs_int_t        *ihisvr,
                                cs_int_t        *tplfmt,
                                const cs_int_t  *isca,
                                const cs_int_t  *iscapp,
                                const cs_int_t  *ipprtp,
                                cs_real_t       *xyzcap);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_OUTPUT_H__ */
