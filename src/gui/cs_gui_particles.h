#ifndef __CS_GUI_PARTICLES_H__
#define __CS_GUI_PARTICLES_H__

/*============================================================================
 * Reader of the parameters file: lagrangian particles
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

void CS_PROCF (uilag1, UILAG1) (int *const nlayer,
                                int *const iilagr,
                                int *const isuila,
                                int *const isuist,
                                int *const isttio,
                                int *const injcon,
                                int *const idepst,
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
                                int *const ltsdyn,
                                int *const ltsmas,
                                int *const ltsthe,
                                int *const nordre,
                                int *const idistu,
                                int *const idiffl,
                                int *const modcpl,
                                int *const idirla,
                                int *const ntlal,
                                int *const ivisv1,
                                int *const ivisv2,
                                int *const ivistp,
                                int *const ivisdm,
                                int *const iviste,
                                int *const ivismp,
                                int *const ivisdk,
                                int *const iviswat,
                                int *const ivisch,
                                int *const ivisck,
                                int *const istala,
                                int *const nbclst,
                                double *const seuil,
                                int *const idstnt,
                                int *const nstist,
                                int ihslag[],
                                int *const iensi3,
                                double *const seuilf,
                                int *const nstbor,
                                int *const inbrbd,
                                int *const iflmbd,
                                int *const iangbd,
                                int *const ivitbd,
                                int *const iencnbbd,
                                int *const iencmabd,
                                int *const iencdibd,
                                int *const iencckbd,
                                int imoybr[],
                                int *const iactfv,
                                int *const iactvx,
                                int *const iactvy,
                                int *const iactvz,
                                int *const iactts);

/*-----------------------------------------------------------------------------
 * Fortran Interface:
 *
 * subroutine uilag2
 * *****************
 *
 * integer          nfabor  -->  number of boundary faces
 * integer          nozppm  -->  max number of boundary conditions zone
 * integer          nclagm  -->  max number of classes
 * integer          nflagm  -->  max number of boundaries
 * integer          iphyla  -->  physica model associated to the particles
 * ..
 * integer          nlayer  <--  number of layer for coal
 * integer          inuchl  <--  particle coal number
 * integer          irawcl  <--  coal particle composition injection condition
 * integer          ihpt    <--  coal temperature in K (for each layer)
 * integer          ifrlag  -->  type of boundary face
 * integer          iusncl  <--  array for particles class(es) number
 * integer          iusclb  <--  array for particles boundary conditions
 *----------------------------------------------------------------------------*/

void CS_PROCF (uilag2, UILAG2) (const int *const nfabor,
                                const int *const nozppm,
                                const int *const ientrl,
                                const int *const isortl,
                                const int *const idepo1,
                                const int *const idepo2,
                                const int *const idepfa,
                                const int *const iencrl,
                                const int *const irebol,
                                const int *const isymtl,
                                const int *const iphyla,
                                const int *const ijnbp,
                                const int *const ijfre,
                                const int *const iclst,
                                const int *const ijuvw,
                                const int *const iuno,
                                const int *const iupt,
                                const int *const ivpt,
                                const int *const iwpt,
                                const int *const ijprpd,
                                const int *const ipoit,
                                const int *const idebt,
                                const int *const ijprdp,
                                const int *const idpt,
                                const int *const ivdpt,
                                const int *const iropt,
                                const int *const ijprtp,
                                const int *const itpt,
                                const int *const icpt,
                                const int *const iepsi,
                                const int *const nlayer,
                                const int *const inuchl,
                                const int *const irawcl,
                                const int  const ihpt[],
                                int     ifrlag[],
                                int     iusncl[],
                                int     iusclb[]);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Free global GUI structures related to particles.
 *----------------------------------------------------------------------------*/

void
cs_gui_particles_free(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_PARTICLES_H__ */
