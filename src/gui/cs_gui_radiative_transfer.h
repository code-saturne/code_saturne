#ifndef __CS_GUI_RADIATIVE_TRANSFER_H__
#define __CS_GUI_RADIATIVE_TRANSFER_H__

/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

void CS_PROCF (uiray1, UIRAY1) (int *const iirayo,
                                int *const isuird,
                                int *const i_quad,
                                int *const ndirec,
                                int *const nfreqr,
                                int *const idiver,
                                int *const iimpar,
                                int *const iimlum);


/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(fcnmra, FCNMRA)
(
 const char      *const fstr,    /* --> Fortran string */
 int             *const len,     /* --> String Length  */
 int             *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(cfnmra, CFNMRA)
(
 char          *const fstr,    /* --> Fortran string */
 int           *const len,     /* --> String Length  */
 int           *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/


void CS_PROCF (uiray2, UIRAY2) (const    int *const itypfb,
                                const    int *const iparoi,
                                const    int *const iparug,
                                const    int *const ivart,
                                         int *const izfrdp,
                                         int *const isothp,
                                const    int *const itpimp,
                                const    int *const ipgrno,
                                const    int *const iprefl,
                                const    int *const ifgrno,
                                const    int *const ifrefl,
                                const    int *const nzoppm,
                                const    int *const nfabor,
                                const    int *const nvar,
                                      double *const epsp,
                                      double *const epap,
                                      double *const tintp,
                                      double *const textp,
                                      double *const xlamp,
                                      double *const rcodcl);

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/


void CS_PROCF (uiray3, UIRAY3) (      double *const ck,
                                const    int *const ncel,
                                         int *const imodak);
/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray4, UIRAY4) (int *const nbrayf,
                                int *const iirayo,
                                int *const irayvf);

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of radiative transfer
 *----------------------------------------------------------------------------*/

void CS_PROCF (uirapr, UIRAPR) (const int *const nprayc,
                                const int *const nrphas,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ilumin,
                                const int *const iqx,
                                const int *const iqy,
                                const int *const iqz,
                                const int *const itsre,
                                const int *const itsri,
                                const int *const iabs,
                                const int *const iemi,
                                const int *const icak);

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *
 * Fortran Interface:
 *
 * SUBROUTINE MEMUI2
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui2, MEMUI2) (void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_RADIATIVE_TRANSFER_H__ */
