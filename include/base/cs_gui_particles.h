/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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

#ifndef __CS_GUI_PARTICLES_H__
#define __CS_GUI_PARTICLES_H__

/*============================================================================
 * Reader of the parameters file: lagrangian particles
 *============================================================================*/

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
 * Copy variable name from C to Fortran
 *----------------------------------------------------------------------------*/

void CS_PROCF(cfname, CFNAME)
(
 int           *const flag,    /* --> flag for array = 1, 2, or 3 */
 char          *const fstr,    /* --> Fortran string */
 int           *const len,     /* --> String Length  */
 int           *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
 );

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fclag1, FCLAG1)
(
 const char          *const fstr,    /* --> Fortran string */
 int                 *const len,     /* --> String Length  */
 int                 *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fclag2, FCLAG2)
(
 const char          *const fstr,    /* --> Fortran string */
 int                 *const len,     /* --> String Length  */
 int                 *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fclag3, FCLAG3)
(
 const char          *const fstr,    /* --> Fortran string */
 int                 *const len,     /* --> String Length  */
 int                 *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

/*-----------------------------------------------------------------------------
 * Lagrangian: global settings, particles model, 2 way coupling, numerical ordering.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UILAG1
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uilag1, UILAG1) (int *const iilagr,
                                int *const isuila,
                                int *const isuist,
                                int *const nbpmax,
                                int *const isttio,
                                int *const injcon,
                                int *const iphyla,
                                int *const idpvar,
                                int *const itpvar,
                                int *const impvar,
                                int *const iencra,
                                double tprenc[],
                                double visref[],
                                double enc1[],
                                double enc2[],
                                int *const nstits,
                                int *const lstdyn,
                                int *const ltsmas,
                                int *const ltsthe,
                                int *const nordre,
                                int *const idistu,
                                int *const idiffl,
                                int *const modcpl,
                                int *const idirla,
                                int *const iensi1,
                                int *const iensi2,
                                int *const ntlal,
                                int *const nbvis,
                                int *const nvisla,
                                int *const ivisv1,
                                int *const ivisv2,
                                int *const ivistp,
                                int *const ivisdm,
                                int *const iviste,
                                int *const ivismp,
                                int *const ivishp,
                                int *const ivisdk,
                                int *const ivisch,
                                int *const ivisck,
                                int *const istala,
                                int *const nbclst,
                                double *const seuil,
                                int *const idstnt,
                                int ihslag[],
                                int *const iensi3,
                                double *const seuilf,
                                int *const nstbor,
                                int *const inbrbd, 
                                int *const iflmbd,
                                int *const iangbd,
                                int *const ivitbd,
                                int *const iencbd,
                                int imoybr[]);

/*-----------------------------------------------------------------------------
 * Fortran Interface:
 *
 * subroutine uilag2
 * *****************
 *
 * integer          iphyla     -->   physica model associated to the particles
 * integer          nclagm     <--   max number of classes for particles
 * integer          nflagm     <--   max number of boundaries
 * integer          iusncl     <--   array for particles class(es) number
 * integer          iusclb     <--   array for particles boundary conditions
 * integer          iuslag     <--   array for integer variables
 * double precision ruslag     <--   array for real variables
 *----------------------------------------------------------------------------*/

void CS_PROCF (uilag2, UILAG2) (const int *const iphyla,
                                const int *const nclagm,
                                const int *const nflagm,
                                int     iusncl[],
                                int     iusclb[],
                                int     iuslag[],
                                double  ruslag[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_PARTICLES_H__ */
