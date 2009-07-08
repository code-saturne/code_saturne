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

#ifndef __CS_GUI_PREDEFINED_PHYSICS_H__
#define __CS_GUI_PREDEFINED_PHYSICS_H__

/*============================================================================
 * Management of the GUI parameters file: specific physics
 *============================================================================*/

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
 * INTEGER          IEQCO2  --> CO2 massic fraction transport
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uippmo, UIPPMO) (int *const ippmod,
                                int *const icod3p,
                                int *const icodeq,
                                int *const icoebu,
                                int *const icobml,
                                int *const icolwc,
                                int *const icp3pl,
                                int *const icpl3c,
                                int *const icfuel,
                                int *const ieljou,
                                int *const ielarc,
                                int *const ielion,
                                int *const icompf,
                                int *const iatmos,
                                int *const iaeros,
                                int *const indjon,
                                int *const ieqco2);

/*----------------------------------------------------------------------------
 * Density under relaxation
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI1 (SRROM)
 * *****************
 * DOUBLE PRECISION SRROM   <--   density relaxation
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi1, UICPI1) (double *const srrom);

/*----------------------------------------------------------------------------
 * Pointers definition for scalars and coal combustion
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpsc, UICPSC) (const int *const ncharb,
                                const int *const nclass,
                                const int *const noxyd,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ieqco2,
                                const int *const ihtco2,
                                const int *const ihm,
                                const int *const inp,
                                const int *const ixch,
                                const int *const ixck,
                                const int *const ixwt,
                                const int *const ih2,
                                const int *const if1m,
                                const int *const if2m,
                                const int *const if3m,
                                const int *const if3mc2,
                                const int *const if4p2m,
                                const int *const if5m,
                                const int *const if6m,
                                const int *const if7m,
                                const int *const iyco2);

/*----------------------------------------------------------------------------
 * Defintion des pointeurs des proprietes pour le charbon
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicppr, UICPPR) (const int *const nclass,
                                const int *const nsalpp,
                                const int *const nsalto,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ihtco2,
                                const int *const itemp1,
                                const int *const irom1,
                                const int *const ym1,
                                const int *const imel,
                                const int *const itemp2,
                                const int *const ix2,
                                const int *const irom2,
                                const int *const idiam2,
                                const int *const igmdch,
                                const int *const igmdv1,
                                const int *const igmdv2,
                                const int *const igmhet,
                                const int *const ighco2,
                                const int *const igmsec);

/*----------------------------------------------------------------------------
 * Atmospheric flows: read of meteorological file of data
 *
 * Fortran Interface:
 *
 * subroutine uiati1
 * *****************
 * integer         imeteo   <--   on/off index
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiati1, UIATI1) (int *const imeteo);

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

#endif /* __CS_GUI_PREDEFINED_PHYSICS_H__ */
