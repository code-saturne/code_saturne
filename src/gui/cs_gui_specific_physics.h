#ifndef __CS_GUI_SPECIFIC_PHYSICS_H__
#define __CS_GUI_SPECIFIC_PHYSICS_H__

/*============================================================================
 * Management of the GUI parameters file: specific physics
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics
 * (pulverized solid fuels)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uisofu, UISOFU) (const int    *iirayo,
                                const int    *ncharm,
                                int          *ncharb,
                                int          *nclpch,
                                int          *nclacp,
                                const int    *ncpcmx,
                                int          *ichcor,
                                double       *diam20,
                                double       *cch,
                                double       *hch,
                                double       *och,
                                double       *nch,
                                double       *sch,
                                int          *ipci,
                                double       *pcich,
                                double       *cp2ch,
                                double       *rho0ch,
                                double       *thcdch,
                                double       *cck,
                                double       *hck,
                                double       *ock,
                                double       *nck,
                                double       *sck,
                                double       *xashch,
                                double       *xashsec,
                                double       *xwatch,
                                double       *h0ashc,
                                double       *cpashc,
                                int          *iy1ch,
                                double       *y1ch,
                                int          *iy2ch,
                                double       *y2ch,
                                double       *a1ch,
                                double       *a2ch,
                                double       *e1ch,
                                double       *e2ch,
                                double       *crepn1,
                                double       *crepn2,
                                double       *ahetch,
                                double       *ehetch,
                                int          *iochet,
                                double       *ahetc2,
                                double       *ehetc2,
                                int          *ioetc2,
                                double       *ahetwt,
                                double       *ehetwt,
                                int          *ioetwt,
                                int          *ieqnox,
                                int          *ieqco2,
                                int          *imdnox,
                                int          *irb,
                                int          *ihtco2,
                                int          *ihth2o,
                                double       *qpr,
                                double       *fn,
                                double       *ckabs1,
                                int          *noxyd,
                                double       *oxyo2,
                                double       *oxyn2,
                                double       *oxyh2o,
                                double       *oxyco2,
                                double       *repnck,
                                double       *repnle,
                                double       *repnlo);

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
 * Activate specific physical models based on XML settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_model_select(void);

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
 *   model_thermo -->  thermophysical model category
 *----------------------------------------------------------------------------*/

const char *
cs_gui_get_thermophysical_model(const char  *model_thermo);

/*----------------------------------------------------------------------------
 * groundwater model : read parameters
 *
 * parameters:
 *   permeability    <--   permeability type
 *   unsteady        <--   steady flow
 *   gravity         <--   check if gravity is taken into account
 *   unsaturated     <--   take into account unsaturated zone
 *----------------------------------------------------------------------------*/

void
cs_gui_gwf_model(int  *permeability,
                 int  *unsteady,
                 int  *gravity,
                 int  *unsaturated);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_SPECIFIC_PHYSICS_H__ */
