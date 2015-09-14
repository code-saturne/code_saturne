#ifndef __CS_GUI_RADIATIVE_TRANSFER_H__
#define __CS_GUI_RADIATIVE_TRANSFER_H__

/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray1, UIRAY1) (int  *iirayo,
                                int  *isuird,
                                int  *i_quad,
                                int  *ndirec,
                                int  *nfreqr,
                                int  *idiver,
                                int  *iimpar,
                                int  *iimlum);

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray2, UIRAY2) (const    int  *itypfb,
                                const    int  *iparoi,
                                const    int  *iparug,
                                const    int  *ivart,
                                         int  *izfrdp,
                                         int  *isothp,
                                const    int  *itpimp,
                                const    int  *ipgrno,
                                const    int  *iprefl,
                                const    int  *ifgrno,
                                const    int  *ifrefl,
                                const    int  *nzoppm,
                                const    int  *nfabor,
                                const    int  *nvar,
                                      double  *epsp,
                                      double  *epap,
                                      double  *tintp,
                                      double  *textp,
                                      double  *xlamp,
                                      double  *rcodcl);

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/


void CS_PROCF (uiray3, UIRAY3) (      double  *ck,
                                const    int  *ncel,
                                         int  *imodak);
/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray4, UIRAY4) (int  *iirayo);

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *
 * Fortran Interface:
 *
 * subroutine memui2
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui2, MEMUI2) (void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_RADIATIVE_TRANSFER_H__ */
