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


void CS_PROCF (csenso, CSENSO) (const    int *const nvppmx,
                                         int *const ncapt,
                                         int *const nthist,
                                      double *const frhist,
                                         int *const ntlist,
                                         int *const iecaux,
                                         int *const ipstdv,
                                         int *const ipstyp,
                                         int *const ipstcl,
                                         int *const ipstft,
                                         int *const ipstfo,
                                         int *const ichrvr,
                                         int *const ilisvr,
                                         int *const ihisvr,
                                         int *const tplfmt,
                                const    int *const isca,
                                const    int *const iscapp,
                                const    int *const ipprtp,
                                      double *const xyzcap);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_OUTPUT_H__ */
