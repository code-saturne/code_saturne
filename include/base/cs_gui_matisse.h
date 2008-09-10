/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_GUI_MATISSE_H__
#define __CS_GUI_MATISSE_H__

/*============================================================================
 * Reader of the parameters file: matisse
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/


#include "cs_base.h"


/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public functions prototypes
 *============================================================================*/


void CS_PROCF (csgein, CSGEIN) (int    *const nptran,
                                int    *const nplgrs,
                                int    *const nelgrs,
                                int    *const nchest,
                                int    *const netran,
                                int    *const itypen);

void CS_PROCF (csgedb, CSGEDB) (double *const epregi,
                                double *const epchem,
                                double *const hconve,
                                double *const rconve,
                                double *const hchali,
                                double *const hcheva,
                                double *const hfttoi,
                                double *const ptrres,
                                double *const frdtra,
                                double *const plgres,
                                double *const epchel,
                                double *const dmcont);

void CS_PROCF (csphdb, CSPHDB) (double *const dtdtmx,
                                double *const puicon,
                                double *const tinit,
                                double *const tcrit,
                                double *const emicon,
                                double *const emimur,
                                double *const hepcnt,
                                double *const dhpcnt,
                                double *const debmas,
                                double *const pdccha,
                                double *const pdcfch,
                                double *const dhchea,
                                double *const sdchea,
                                double *const pdcche,
                                double *const pdccch,
                                double *const dhches,
                                double *const sdches,
                                double *const pdcalg,
                                double *const pdcatv,
                                double *const argamt,
                                double *const pdcslg,
                                double *const pdcstv,
                                double *const argavl,
                                double *const amppdc,
                                double *const dhalve,
                                double *const hreso,
                                double *const hplen,
                                double *const dpvent);

void CS_PROCF (csphat, CSPHAT) (int    *const imdcnt,
                                int    *const icofor,
                                int    *const iconlg,
                                int    *const ialveo);

void CS_PROCF (csmtpr, CSMTPR) (int          *imatis);

void CS_PROCF (csnbmp, CSNBMP) (int    *const direction,
                                int    *const carte,
                                int    *const nb);

void CS_PROCF (csdfmp, CSDFMP) (int    *const zone,
                                int    *const direction,
                                int    *const carte,
                                double *const min,
                                double *const max,
                                double *const value);

void CS_PROCF (csmhdb, CSMHDB) (double *const jeuchr,
                                double *const jeurcl,
                                double *const jeuclr,
                                double *const jeurch,
                                int    *const nechrg,
                                int    *const nergrs,
                                int    *const neclrg,
                                int    *const nergch,
                                double *const hbdtoi,
                                int    *const neciel);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_MATISSE_H__ */


