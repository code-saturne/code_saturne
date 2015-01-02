#ifndef __CS_GUI_H__
#define __CS_GUI_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
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
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialise the global 'vars' structure.
 *
 * Fortran Interface:
 *
 * subroutine uiinit
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiinit, UIINIT) (void);

/*----------------------------------------------------------------------------
 * Thermal model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTHER (ITHERM)
 * *****************
 *
 * INTEGER          ITHERM  --> thermal model
 * integer          itpscl  --> temperature scale if itherm = 1
 *----------------------------------------------------------------------------*/


void CS_PROCF (csther, CSTHER) (int  *itherm,
                                int  *itpscl);

/*----------------------------------------------------------------------------
 * Turbulence model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTURB
 * *****************
 *
 * INTEGER          ITURB   -->   turbulence model
 * INTEGER          IWALLF  -->   wall law treatment
 * INTEGER          IGRAKE  -->   k-eps gravity effects
 * INTEGER          IGRAKI  -->   Rij-eps gravity effects
 * DOUBLE PRECISION XLOMLG  -->   mixing_length_scale
 *----------------------------------------------------------------------------*/

void CS_PROCF (csturb, CSTURB) (int    *iturb,
                                int    *iwallf,
                                int    *igrake,
                                int    *igrari,
                                double *xlomlg);

/*----------------------------------------------------------------------------
 * Specific heat variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCPVA
 * *****************
 *
 * INTEGER          ICP     -->   Specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscpva, CSCPVA) (int *icp);

/*----------------------------------------------------------------------------
 * Volumic viscosity variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCVVVA (ICP)
 * *****************
 *
 * INTEGER          IVISCV     -->   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvvva, CSVVVA) (int *iviscv);

/*----------------------------------------------------------------------------
 * User thermal scalar.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UITHSC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uithsc, UITHSC) (void);

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * subroutine csivis
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (void);

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIDTV (IDTVAR)
 * *****************
 *
 * INTEGER          IDTVAR  -->   fixed or variable time step
 *----------------------------------------------------------------------------*/

void CS_PROCF(csidtv, CSIDTV) (int *idtvar);

/*----------------------------------------------------------------------------
 * Hydrostatic pressure parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIPHY (IPHYDR)
 * *****************
 *
 * INTEGER          IPHYDR  -->   hydrostatic pressure
 *----------------------------------------------------------------------------*/

void CS_PROCF (csiphy, CSIPHY) (int *iphydr);

/*----------------------------------------------------------------------------
 * Hydrostatic equilibrium parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCFGP (ICFGRP)
 * *****************
 *
 * INTEGER          ICFGRP  -->   hydrostatic equilibrium
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscfgp, CSCFGP) (int *icfgrp);

/*----------------------------------------------------------------------------
 * Restart parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISUI
 * *****************
 *
 * INTEGER          NTSUIT  -->   checkpoint frequency
 * INTEGER          ILEAUX  -->   restart with auxiliary
 * INTEGER          ICCFVG  -->   restart with frozen field
 *----------------------------------------------------------------------------*/


void CS_PROCF (csisui, CSISUI) (int *ntsuit,
                                int *ileaux,
                                int *iccvfg);

/*----------------------------------------------------------------------------
 * Time passing parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTIME
 * *****************
 *
 * INTEGER          INPDT0  -->   zero tim step
 * INTEGER          IPTLTO  -->   thermal time step control
 * INTEGER          NTMABS  -->   iterations numbers
 * INTEGER          IDTVAR  -->   time step's options
 * DOUBLE PRECISION DTREF   -->   time step
 * DOUBLE PRECISION DTMIN   -->   minimal time step
 * DOUBLE PRECISION DTMAX   -->   maximal time step
 * DOUBLE PRECISION COUMAX  -->   maximal courant number
 * DOUBLE PRECISION FOUMAX  -->   maximal fournier number
 * DOUBLE PRECISION VARRDT  -->   max time step variation between 2 iterations
 * DOUBLE PRECISION RELXST  -->   relaxation coefficient if idtvar = -1
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstime, CSTIME) (int    *inpdt0,
                                int    *iptlro,
                                int    *ntmabs,
                                int    *idtvar,
                                double *dtref,
                                double *dtmin,
                                double *dtmax,
                                double *coumax,
                                double *foumax,
                                double *varrdt,
                                double *relxst);

/*----------------------------------------------------------------------------
 *
 * Fortran Interface:
 *
 * SUBROUTINE UINUM1
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (double *blencv,
                                   int *ischcv,
                                   int *isstpc,
                                   int *ircflu,
                                double *cdtvar,
                                double *epsilo,
                                   int *nswrsm);

/*----------------------------------------------------------------------------
 * Global numerical parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNUM2
 * *****************
 *
 * INTEGER          IVISSE  -->   gradient transpose
 * INTEGER          RELAXP  -->   pressure relaxation
 * INTEGER          IPUCOU  -->   velocity pressure coupling
 * INTEGER          EXTRAG  -->   wall pressure extrapolation
 * INTEGER          IMRGRA  -->   gradient reconstruction
 * INTEGER          NTERUP  -->   piso sweep number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2) (   int *ivisse,
                                double *relaxp,
                                   int *ipucou,
                                double *extrag,
                                   int *imrgra,
                                   int *nterup);

void CS_PROCF (csphys, CSPHYS) (const    int *nmodpp,
                                         int *irovar,
                                         int *ivivar,
                                         int *icorio,
                                      double *gx,
                                      double *gy,
                                      double *gz,
                                      double *ro0,
                                      double *viscl0,
                                      double *viscv0,
                                      double *visls0,
                                      double *cp0,
                                      double *t0,
                                      double *p0,
                                      double *xmasmr,
                                         int *itempk);

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * subroutine cssca2
 * *****************
 *
 * integer          iturb    <--  turbulence model
 * integer          iturt    -->  turbulent flux model
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (const int  *iturb,
                                int        *iturt);

void CS_PROCF (cssca3, CSSCA3) (double     *visls0,
                                double     *t0,
                                double     *p0,
                                double     *cp0);

/*----------------------------------------------------------------------------
 * Turbulence initialization parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTINI
 * *****************
 *
 * INTEGER          UREF   -->   reference velocity
 * INTEGER          ALMAX  -->   reference length
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstini, CSTINI) (double  *uref,
                                double  *almax);

/*----------------------------------------------------------------------------
 * Solver taking a scalar porosity into account
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIIPSU
 * *****************
 *
 * INTEGER          IPOROS     -->   porosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiipsu, UIIPSU) (int *iporos);

/*----------------------------------------------------------------------------
 * Define porosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPORO
 * *****************
 *
 * INTEGER          IPOROS     <--   porosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiporo, UIPORO) (const int *ncelet,
                                const int *iporos);

/*----------------------------------------------------------------------------
 * User momentum source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsnv (ncelet, vel, tsexp, tsimp)
 * *****************
 *
 * integer          ncelet   <--  number of cells with halo
 * double precision vel      <--  fluid velocity
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsnv, UITSNV)(const cs_real_3_t *restrict vel,
                              cs_real_3_t       *restrict tsexp,
                              cs_real_33_t      *restrict tsimp);


/*----------------------------------------------------------------------------
 * User scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitssc (f_id, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitssc, UITSSC)(const int                  *f_id,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp);


/*----------------------------------------------------------------------------
 * Thermal scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsth (f_id, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsth, UITSTH)(const int                  *f_id,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp);

/*----------------------------------------------------------------------------
 * Variables and user scalars initialization.
 *
 * Fortran Interface:
 *
 * subroutine uiiniv
 * *****************
 *
 * integer          ncelet   <--  number of cells with halo
 * integer          isuite   <--  restart indicator
 * integer          idarcy   <--  darcy module activate or not
 * integer          iccfth   <--  type of initialisation(compressible model)
 * double precision ro0      <--  value of density if IROVAR=0
 * double precision cp0      <--  value of specific heat if ICP=0
 * double precision viscl0   <--  value of viscosity if IVIVAR=0
 * double precision uref     <--  value of reference velocity
 * double precision almax    <--  value of reference length
 * double precision xyzcen   <--  cell's gravity center
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int          *ncelet,
                              const int          *isuite,
                              const int          *idarcy,
                                    int          *iccfth,
                              const cs_real_t    *ro0,
                              const cs_real_t    *cp0,
                              const cs_real_t    *viscl0,
                              const cs_real_t    *uref,
                              const cs_real_t    *almax,
                              const double       *xyzcen);

/*----------------------------------------------------------------------------
 * User law for material Properties
 *
 * Fortran Interface:
 *
 * subroutine uiphyv
 * *****************
 *
 * integer          ncel     <--  number of cells whithout halo
 * integer          ncelet   <--  number of cells whith halo
 * integer          icp      <--  pointer for specific heat Cp
 * integer          irovar   <--  =1 if rho variable, =0 if rho constant
 * integer          ivivar   <--  =1 if mu variable, =0 if mu constant
 * integer          iviscv   <--  pointer for volumic viscosity viscv
 * integer          itempk   <--  pointer for temperature (in K)
 * double precision p0       <--  pressure reference value
 * double precision t0       <--  temperature reference value
 * double precision ro0      <--  density reference value
 * double precision cp0      <--  specific heat reference value
 * double precision viscl0   <--  dynamic viscosity reference value
 * double precision visls0   <--  diffusion coefficient of the scalars
 * double precision viscv0   <--  volumic viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const cs_int_t  *ncel,
                              const cs_int_t  *ncelet,
                              const cs_int_t  *icp,
                              const cs_int_t  *irovar,
                              const cs_int_t  *ivivar,
                              const cs_int_t  *iviscv,
                              const cs_int_t  *itempk,
                              const cs_real_t *p0,
                              const cs_real_t *t0,
                              const cs_real_t *ro0,
                              const cs_real_t *cp0,
                              const cs_real_t *viscl0,
                              const cs_real_t *visls0,
                              const cs_real_t *viscv0);

/*----------------------------------------------------------------------------
 * Head losses definition
 *
 * Fortran Interface:
 *
 * subroutine uikpdc
 * *****************
 *
 * integer          iappel   <--  number of calls during a time step
 * integer          ncelet   <--  number of cells with halo
 * integer          ncepdp  -->   number of cells with head losses
 * integer          icepdc  -->   ncepdp cells number with head losses
 * double precision ckupdc  -->   head losses matrix
 *----------------------------------------------------------------------------*/

void CS_PROCF(uikpdc, UIKPDC)(const int*   iappel,
                              const int*   ncelet,
                                    int*   ncepdp,
                                    int    icepdc[],
                                    double ckupdc[]);

/*----------------------------------------------------------------------------
 * 1D profile postprocessing
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPROF
 * *****************
 *
 * INTEGER          NCELET   <--  number of cells with halo
 * INTEGER          NCEL     <--  number of cells without halo
 * INTEGER          NTMABS   <--  max iterations numbers
 * INTEGER          NTCABS   <--  current iteration number
 * DOUBLE PRECISION TTCABS   <--  current physical time
 * DOUBLE PRECISION TTMABS   <--  max physical time
 * DOUBLE PRECISION TTPABS   <--  physical time at calculation beginning
 * DOUBLE PRECISION XYZCEN   <--  cell's gravity center
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprof, UIPROF)(const int    *ncelet,
                               const int    *ncel,
                               const int    *ntmabs,
                               const int    *ntcabs,
                               const double *ttcabs,
                               const double *ttmabs,
                               const double *ttpabs,
                               const double *xyzcen);


/*----------------------------------------------------------------------------
 * darcy model : read laws for capacity, saturation and permeability
 *
 * Fortran Interface:
 *
 * subroutine uidapp
 * *****************
 * integer         permeability    <--  permeability type
 * integer         diffusion       <--  diffusion type
 * integer         gravity         <--  check if gravity is taken into account
 * double          gravity_x       <--   x component for gravity vector
 * double          gravity_y       <--   y component for gravity vector
 * double          gravity_z       <--   z component for gravity vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (uidapp, UIDAPP) (const cs_int_t  *permeability,
                                const cs_int_t  *diffusion,
                                const cs_int_t  *gravity,
                                const double    *gravity_x,
                                const double    *gravity_y,
                                const double    *gravity_z);

/*----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables.
 *
 * Fortran Interface:
 *
 * SUBROUTINE MEMUI1
 * *****************
 *
 * INTEGER          NCHARB  <-- number of coal
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui1, MEMUI1) (const int *ncharb);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get thermal scalar model.
 *
 * return:
 *   value of itherm
 *----------------------------------------------------------------------------*/

int
cs_gui_thermal_model(void);

/*----------------------------------------------------------------------------
 * Define user variables through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_user_variables(void);

/*-----------------------------------------------------------------------------
 * Get initial value from property markup.
 *
 * parameters:
 *   property_name      <--  name of the property
 *   value              -->  new initial value of the property
 *----------------------------------------------------------------------------*/

void
cs_gui_properties_value(const char  *property_name,
                        double      *value);

/*-----------------------------------------------------------------------------
 * Initialization choice of the reference variables parameters.
 *
 * parameters:
 *   name            <--   parameter name
 *   value           -->   parameter value
 *----------------------------------------------------------------------------*/

void
cs_gui_reference_initialization(const char  *param,
                                double      *value);

/*-----------------------------------------------------------------------------
 * Selection of linear solvers.
 *----------------------------------------------------------------------------*/

void
cs_gui_linear_solvers(void);

/*-----------------------------------------------------------------------------
 * Set turbomachinery model
 *----------------------------------------------------------------------------*/

void
cs_gui_turbomachinery(void);

/*-----------------------------------------------------------------------------
 * Set turbomachinery options.
 *----------------------------------------------------------------------------*/

void
cs_gui_turbomachinery_rotor(void);

/*-----------------------------------------------------------------------------
 * Set partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_gui_partition(void);

/*-----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_parallel_io(void);

/*----------------------------------------------------------------------------
 * Time moments definition
 *----------------------------------------------------------------------------*/

void
cs_gui_time_moments(void);

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *----------------------------------------------------------------------------*/

void
cs_gui_clean_memory(void);

/*----------------------------------------------------------------------------
 * Logging output for MEI usage.
 *----------------------------------------------------------------------------*/

void
cs_gui_usage_log(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_H__ */
