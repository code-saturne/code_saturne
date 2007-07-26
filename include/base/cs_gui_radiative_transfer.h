/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
*
*     contact: saturne-support@edf.fr
*
*     The Code_Saturne Kernel is free software; you can redistribute it
*     and/or modify it under the terms of the GNU General Public License
*     as published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*
*     The Code_Saturne Kernel is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
*     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public License
*     along with the Code_Saturne Kernel; if not, write to the
*     Free Software Foundation, Inc.,
*     51 Franklin St, Fifth Floor,
*     Boston, MA  02110-1301  USA
*
*============================================================================*/

#ifndef __CS_GUI_RADIATIVE_TRANSFER_H__
#define __CS_GUI_RADIATIVE_TRANSFER_H__

/*============================================================================
 * Reader of the parameters file: radiative transfer
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/


#include "cs_base.h"


/*----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 * C API public functions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *----------------------------------------------------------------------------*/

void cs_gui_clean_memory_rayt(void);

/*============================================================================
 * Fortran API public functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray1, UIRAY1) (int *const nbrayb,
                                int *const nbrayf,
                                int *const nphas,
                                int *const irayon,
                                int *const isuird,
                                int *const ndirec,
                                int *const nfreqr,
                                int *const iimpar,
                                int *const idiver,
                                int *const iimlum,
                                int *const irayvp,
                                int *const irayvf);

void CS_PROCF(fcnmra, FCNMRA)
(
 const char      *const fstr,    /* --> Fortran string */
 int             *const len,     /* --> String Length  */
 int             *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

void CS_PROCF(cfnmra, CFNMRA)
(
 char          *const fstr,    /* --> Fortran string */
 int           *const len,     /* --> String Length  */
 int           *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);


void CS_PROCF (uiray2, UIRAY2) (const    int *const itypfb,
                                const    int *const iparoi,
                                const    int *const ivart,
                                const    int *const iph,
                                const    int *const nphast,
                                         int *const izfrdp,
                                         int *const isothp,
                                const    int *const itpimp,
                                const    int *const ipgrno,
                                const    int *const iprefl,
                                const    int *const ifgrno,
                                const    int *const ifrefl,
                                const    int *const nfabor,
                                const    int *const nfml,
                                const    int *const ifmfbr,
                                const    int *const iprfml,
                                const    int *const nvar,
                                      double *const epsp,
                                      double *const epap,
                                      double *const tintp,
                                      double *const textp,
                                      double *const xlamp,
                                      double *const rcodcl);


void CS_PROCF (uiray3, UIRAY3)
(      double *const ck,
 const    int *const iph,
 const    int *const ncelet,
 const    int *const ncel);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_GUI_RADIATIVE_TRANSFER_H__ */
