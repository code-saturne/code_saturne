#ifndef __CS_GUI_SPECIFIC_PHYSICS_H__
#define __CS_GUI_SPECIFIC_PHYSICS_H__

/*============================================================================
 * Management of the GUI parameters file: specific physics
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
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public Fortran function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Predefined physics indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPPMO
 * *****************
 *
 * INTEGER          IPPMOD <--  specific physics indicator array
 * INTEGER          ICOD3P  --> diffusion flame in fast complete chemistry
 * INTEGER          ICODEQ  --> diffusion flame in fast chemistry to equilibrium
 * INTEGER          ICOEBU  --> Eddy Break Up premixing flame
 * INTEGER          ICOBML  --> Bray - Moss - Libby premixing flame
 * INTEGER          ICOLWC  --> Libby Williams premixing flame
 * INTEGER          ICP3PL  --> Coal combustion. Combustible moyen local
 * INTEGER          ICPL3C  --> Coal combustion coupled with lagrangien approach
 * INTEGER          ICFUEL  --> Fuel combustion
 * INTEGER          IELJOU  --> Joule effect
 * INTEGER          IELARC  --> electrical arc
 * INTEGER          IELION  --> ionique mobility
 * INTEGER          ICOMPF  --> compressible without shock
 * INTEGER          IATMOS  --> atmospheric flows
 * INTEGER          IAEROS  --> cooling tower
 * INTEGER          INDJON  --> INDJON=1: a JANAF enthalpy-temperature
 *                              tabulation is used. INDJON=1: users tabulation
 * INTEGER          IEOS    --> compressible
 * INTEGER          IEQCO2  --> CO2 massic fraction transport
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uippmo, UIPPMO) (int *const ippmod,
                                int *const icod3p,
                                int *const icodeq,
                                int *const icoebu,
                                int *const icobml,
                                int *const icolwc,
                                int *const iccoal,
                                int *const icpl3c,
                                int *const icfuel,
                                int *const ieljou,
                                int *const ielarc,
                                int *const ielion,
                                int *const icompf,
                                int *const iatmos,
                                int *const iaeros,
                                int *const ieos,
                                int *const ieqco2);

/*----------------------------------------------------------------------------
 * Density under relaxation
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI1 (SRROM)
 * *****************
 * DOUBLE PRECISION SRROM   <--   density relaxation
 * DOUBLE PRECISION DIFTL0  <--   dynamic diffusion
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi1, UICPI1) (double *const srrom,
                                double *const diftl0);

/*----------------------------------------------------------------------------
 * Temperature for D3P Gas Combustion
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI2 (SRROM)
 * *****************
 * DOUBLE PRECISION Toxy   <--   Oxydant temperature
 * DOUBLE PRECISION Tfuel  <--   Fuel temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi2, UICPI2) (double *const toxy,
                                double *const tfuel);

/*----------------------------------------------------------------------------
 * Pointers definition for scalars and coal combustion
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpsc, UICPSC) (const int *const ncharb,
                                const int *const nclass,
                                const int *const noxyd,
                                const int *const ippmod,
                                const int *const iccoal,
                                const int *const ieqnox,
                                const int *const ieqco2,
                                const int *const ihtco2,
                                const int *const ihth2o,
                                const int *const ihm,
                                const int *const inp,
                                const int *const ixch,
                                const int *const ixck,
                                const int *const ixwt,
                                const int *const ih2,
                                const int *const if1m,
                                const int *const if2m,
                                const int *const if4m,
                                const int *const if5m,
                                const int *const if6m,
                                const int *const if7m,
                                const int *const if8m,
                                const int *const ifvp2m,
                                const int *const iyco2,
                                const int *const if9m,
                                const int *const iyhcn,
                                const int *const iyno,
                                const int *const ihox,
                                const int *const iynh3);

/*----------------------------------------------------------------------------
 * Defintion des pointeurs des proprietes pour la combustion gaz
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicppr, UICPPR) (const int *const nclass,
                                const int *const nsalpp,
                                const int *const ippmod,
                                const int *const iccoal,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ieqnox,
                                const int *const ihtco2,
                                const int *const ihth2o,
                                const int *const itemp1,
                                const int *const irom1,
                                const int *const ym1,
                                const int *const ighcn1,
                                const int *const ighcn2,
                                const int *const ignoth,
                                const int *const ignh31,
                                const int *const ignh32,
                                const int *const ifhcnd,
                                const int *const ifhcnc,
                                const int *const ifnh3d,
                                const int *const ifnh3c,
                                const int *const ifnohc,
                                const int *const ifnonh,
                                const int *const ifnoch,
                                const int *const ifnoth,
                                const int *const icnohc,
                                const int *const icnonh,
                                const int *const ifhcnr,
                                const int *const icnorb,
                                const int *const igrb,
                                const int *const immel,
                                const int *const itemp2,
                                const int *const ix2,
                                const int *const irom2,
                                const int *const idiam2,
                                const int *const igmdch,
                                const int *const igmdv1,
                                const int *const igmdv2,
                                const int *const igmhet,
                                const int *const ighco2,
                                const int *const ighh2o,
                                const int *const igmsec,
                                const int *const ibcarbone,
                                const int *const iboxygen,
                                const int *const ibhydrogen);

/*----------------------------------------------------------------------------
 * Pointers definition for scalars for compressible model
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicfsc, UICFSC) (const int *const ienerg,
                                const int *const itempk);

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics (electrical model)
 *----------------------------------------------------------------------------*/
void CS_PROCF (uielpr, UIELPR) (const int *const nsalpp,
                                const int *const ippmod,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ieljou,
                                const int *const ielarc,
                                const int *const itemp,
                                const int *const iefjou,
                                const int *const idjr,
                                const int *const idji,
                                const int *const ilapla,
                                const int *const idrad,
                                const int *const ivisls,
                                const int *const ipotr,
                                const int *const ixkabe);

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar (electrical model)
 *----------------------------------------------------------------------------*/
void CS_PROCF (uielsc, UIELSC) (const int *const ippmod,
                                const int *const ieljou,
                                const int *const ielarc,
                                const int *const ngazg,
                                const int *const ihm,
                                const int *const ipotr,
                                const int *const iycoel,
                                const int *const ipoti,
                                const int *const ipotva);

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics (gaz combustion)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicopr, UICOPR) (const int *const nsalpp,
                                const int *const ippmod,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icolwc,
                                const int *const iirayo,
                                const int *const itemp,
                                const int *const imam,
                                const int *const iym,
                                const int *const ickabs,
                                const int *const it4m,
                                const int *const it3m,
                                const int *const itsc,
                                const int *const irhol,
                                const int *const iteml,
                                const int *const ifmel,
                                const int *const ifmal,
                                const int *const iampl,
                                const int *const itscl,
                                const int *const imaml);

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar (gas combustion)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicosc, UICOSC) (const int *const ippmod,
                                const int *const icolwc,
                                const int *const icoebu,
                                const int *const icod3p,
                                const int *const ihm,
                                const int *const ifm,
                                const int *const ifp2m,
                                const int *const iygfm,
                                const int *const iyfm,
                                const int *const iyfp2m,
                                const int *const icoyfp);

/*----------------------------------------------------------------------------
 * Electrical model : read parameters
 *
 * Fortran Interface:
 *
 * subroutine uieli1
 * *****************
 * integer         ieljou    -->   joule model
 * integer         ielarc    -->   arc model
 * integer         ielcor    <--   scaling electrical variables
 * double          couimp    <--   imposed current intensity
 * double          puisim    <--   imposed power
 * integer         modrec    <--   scaling type for electric arc
 * integer         idrecal   <--   current density component used to scaling
 *                                 (modrec ==2)
 * char            crit_reca <--   define criteria for plane used to scaling (modrec ==2)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uieli1, UIELI1) (const int    *const ncelet,
                                const int    *const ieljou,
                                const int    *const ielarc,
                                      int    *const ielcor,
                                      double *const couimp,
                                      double *const puisim,
                                      int    *const modrec,
                                      int    *const idreca,
                                      double *const crit_reca);

/*----------------------------------------------------------------------------
 * Electrical model : define plane for elreca
 *
 * Fortran Interface:
 *
 * subroutine uielrc
 * *****************
 * integer         izreca    <--   define plane used to scaling (modrec ==2)
 * char            crit_reca <--   define criteria for plane used to scaling (modrec ==2)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uielrc, UIELRC) (const int    *const ncelet,
                                      int    *const izreca,
                                      double *const crit_reca);

/*----------------------------------------------------------------------------
 * Atmospheric flows: read of meteorological file of data
 *
 * Fortran Interface:
 *
 * subroutine uiati1
 * *****************
 * integer         imeteo   <--   on/off index
 * char(*)         fmeteo   <--   meteo file name
 * int             len      <--   meteo file name destination string length
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiati1, UIATI1) (int           *imeteo,
                                char          *fmeteo,
                                int           *len
                                CS_ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 * Atmospheric flows: indirection between the solver numbering and the XML one
 * for physical properties
 *
 * Fortran Interface:
 *
 * subroutine uiatpr
 * *****************
 * integer         nsalpp   -->
 * integer         nsalto   -->
 * integer         ippmod   -->   specific physics indicator array
 * integer         iatmos   -->   index for atmospheric flow
 * integer         ipppro   -->
 * integer         ipproc   -->
 * integer         itempc   -->   index for real temperature
 * integer         iliqwt   -->   index for liquid water
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiatpr, UIATPR) (const int *const nsalpp,
                                const int *const nsalto,
                                const int *const ippmod,
                                const int *const iatmos,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const itempc,
                                const int *const iliqwt);

/*----------------------------------------------------------------------------
 * Atmospheric flows: indirection between the solver numbering and the XML one
 * for models scalars.
 *
 * Fortran Interface:
 *
 * subroutine uiatsc
 * *****************
 * integer         ippmod   -->   specific physics indicator array
 * integer         iatmos   -->   index for atmospheric flow
 * integer         itempp   -->   index for potential temperature
 * integer         itempl   -->   index for liquid potential temperature
 * integer         itotwt   -->   index for total water content
 * integer         intdrp   -->   index for total number of droplets
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiatsc, UIATSC) (const int *const ippmod,
                                const int *const iatmos,
                                const int *const itempp,
                                const int *const itempl,
                                const int *const itotwt,
                                const int *const intdrp);

/*----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics (pulverized solid fuels)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uisofu, UISOFU) (const int    *const iirayo,
                                const int    *const iihmpr,
                                const int    *const ncharm,
                                      int    *const ncharb,
                                      int    *const nclpch,
                                      int    *const nclacp,
                                const int    *const ncpcmx,
                                      int    *const ichcor,
                                      double *const diam20,
                                      double *const cch,
                                      double *const hch,
                                      double *const och,
                                      double *const nch,
                                      double *const sch,
                                      double *const ipci,
                                      double *const pcich,
                                      double *const cp2ch,
                                      double *const rho0ch,
                                      double *const cck,
                                      double *const hck,
                                      double *const ock,
                                      double *const nck,
                                      double *const sck,
                                      double *const pcick,
                                      double *const xashch,
                                      double *const xashsec,
                                      double *const xwatch,
                                      double *const h0ashc,
                                      double *const cpashc,
                                      int    *const iy1ch,
                                      double *const y1ch,
                                      int    *const iy2ch,
                                      double *const y2ch,
                                      double *const a1ch,
                                      double *const a2ch,
                                      double *const e1ch,
                                      double *const e2ch,
                                      double *const crepn1,
                                      double *const crepn2,
                                      double *const ahetch,
                                      double *const ehetch,
                                      int    *const iochet,
                                      double *const ahetc2,
                                      double *const ehetc2,
                                      int    *const ioetc2,
                                      double *const ahetwt,
                                      double *const ehetwt,
                                      int    *const ioetwt,
                                      int    *const ieqnox,
                                      int    *const imdnox,
                                      int    *const irb,
                                      int    *const ihtco2,
                                      int    *const ihth2o,
                                      double *const qpr,
                                      double *const fn,
                                      double *const ckabs1,
                                      int    *const noxyd,
                                      double *const oxyo2,
                                      double *const oxyn2,
                                      double *const oxyh2o,
                                      double *const oxyco2,
                                      double *const repnck,
                                      double *const repnle,
                                      double *const repnlo);

/*----------------------------------------------------------------------------
 * Copy name of thermophysical data file from C to Fortran
 *----------------------------------------------------------------------------*/

void CS_PROCF(cfnmtd, CFNMTD) (char          *fstr,    /* --> Fortran string */
                               int           *len      /* --> String Length  */
                               CS_ARGF_SUPP_CHAINE);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the name of a thermophysical model.
 *
 * parameter:
 *   model_thermo          -->  thermophysical model
 *----------------------------------------------------------------------------*/

char *
cs_gui_get_thermophysical_model(const char *const model_thermo);

/*-----------------------------------------------------------------------------
 * Modify double numerical parameters.
 *
 * parameters:
 *   param               -->  label of the numerical parameter
 *   keyword            <-->  value of the numerical parameter
 *----------------------------------------------------------------------------*/

void
cs_gui_numerical_double_parameters(const char   *const param,
                                         double *const keyword);

/*-----------------------------------------------------------------------------
 * Return if a predifined physics model is activated.
 *----------------------------------------------------------------------------*/

int
cs_gui_get_activ_thermophysical_model(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_SPECIFIC_PHYSICS_H__ */
