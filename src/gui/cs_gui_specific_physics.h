#ifndef __CS_GUI_SPECIFIC_PHYSICS_H__
#define __CS_GUI_SPECIFIC_PHYSICS_H__

/*============================================================================
 * Management of the GUI parameters file: specific physics
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
                                      int    *const ipci,
                                      double *const pcich,
                                      double *const cp2ch,
                                      double *const rho0ch,
                                      double *const thcdch,
                                      double *const cck,
                                      double *const hck,
                                      double *const ock,
                                      double *const nck,
                                      double *const sck,
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


/*----------------------------------------------------------------------------
 * groundwater model : read parameters
 *
 * Fortran Interface:
 *
 * subroutine uidai1
 * *****************
 * integer         permeability    <--   permeability type
 * integer         dispersion      <--   dispersion type
 * integer         unsteady        <--   steady flow
 * integer         gravity         <--   check if gravity is taken into account
 * integer         unsaturated     <--   take into account unsaturated zone
 *----------------------------------------------------------------------------*/

void CS_PROCF (uidai1, UIDAI1) (int    *const permeability,
                                int    *const dispersion,
                                int    *const unsteady,
                                int    *const gravity,
                                int    *const unsaturated);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Activate specific physical models based on XML settings.
 *
 * parameters:
 *   ieos    --> compressible
 *   ieqco2  --> CO2 massic fraction transport (for combustion only)
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_model_select(cs_int_t  *ieos,
                             cs_int_t  *ieqco2);

/*----------------------------------------------------------------------------
 * Electrical model: read parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_elec_model(void);

/*----------------------------------------------------------------------------
 * Electrical model: define plane for elreca
 *----------------------------------------------------------------------------*/

void
cs_gui_elec_model_rec(void);

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
