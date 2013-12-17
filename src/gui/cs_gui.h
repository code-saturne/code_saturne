#ifndef __CS_GUI_H__
#define __CS_GUI_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
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
 * INTEGER          IDEUCH  -->   wall law treatment
 * INTEGER          IGRAKE  -->   k-eps gravity effects
 * INTEGER          IGRAKI  -->   Rij-eps gravity effects
 * DOUBLE PRECISION XLOMLG  -->   mixing_length_scale
 *----------------------------------------------------------------------------*/

void CS_PROCF (csturb, CSTURB) (int *const iturb,
                                int *const ideuch,
                                int *const igrake,
                                int *const igrari,
                                double *const xlomlg);

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

void CS_PROCF (cscpva, CSCPVA) (int *const icp);

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

void CS_PROCF (csvvva, CSVVVA) (int *const iviscv);

/*----------------------------------------------------------------------------
 * User scalars number.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNSCA
 * *****************
 *
 * INTEGER          NSCAUS     -->   user scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnsca, CSNSCA) (int *const nscaus);

/*----------------------------------------------------------------------------
 * User scalars labels
 *
 * Fortran Interface:
 *
 * SUBROUTINE UISCAU
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiscau, UISCAU) (void);

/*----------------------------------------------------------------------------
 * User thermal scalar.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UITHSC
 * *****************
 *
 * INTEGER          ISCALT     -->   thermal scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uithsc, UITHSC) (int *const iscalt);

/*----------------------------------------------------------------------------
 * User scalars which are variance.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISCA (ISCAVR)
 * *****************
 *
 * INTEGER          ISCAVR     -->   user scalars variance array
 * integer          itherm     <--  type of thermal model
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisca, CSISCA) (      int *const iscavr,
                                      int *const itherm,
                                const int *const iscapp);

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIVIS
 * *****************
 *
 * INTEGER          ISCAVR  <->  number of the related variance if any
 * INTEGER          IVISLS  -->  indicator for the user scalar viscosity
 * INTEGER          ISCALT  <->  number of the user thermal scalar if any
 * INTEGER          ISCSTH  <->  type of the user thermal scalar
 * INTEGER          ITEMPK  -->  rtp index for temperature (in K)
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * subroutine csivis (iscavr, ivisls, iscalt, itherm, itempk)
 * *****************
 *
 * integer          iscavr  <-->  number of the related variance if any
 * integer          ivisls  <--   indicator for the user scalar viscosity
 * integer          iscalt  <-->  number of the user thermal scalar if any
 * integer          itherm  <-->  type of thermal model
 * integer          itempk   -->  rtp index for temperature (in K)
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (int *const iscavr,
                                int *const ivisls,
                                int *const iscalt,
                                int *const itherm,
                                int *const itempk);

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

void CS_PROCF(csidtv, CSIDTV) (int *const idtvar);

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

void CS_PROCF (csiphy, CSIPHY) (int *const iphydr);

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

void CS_PROCF (cscfgp, CSCFGP) (int *const icfgrp);

/*----------------------------------------------------------------------------
 *
 * SUBROUTINE CSVNUM()
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvnum, CSVNUM) (const int *const nvar,
                                const int *const iu,
                                const int *const iv,
                                const int *const iw,
                                const int *const ipr,
                                const int *const iturb,
                                const int *const ik,
                                const int *const iep,
                                const int *const ir11,
                                const int *const ir22,
                                const int *const ir33,
                                const int *const ir12,
                                const int *const ir13,
                                const int *const ir23,
                                const int *const iomg,
                                const int *const iphi,
                                const int *const ifb,
                                const int *const ial,
                                const int *const inusa,
                                const int *const iale,
                                const int *const iuma,
                                const int *const ivma,
                                const int *const iwma,
                                const int *const isca,
                                const int *const iscapp,
                                const int *const itherm);

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


void CS_PROCF (csisui, CSISUI) (int *const ntsuit,
                                int *const ileaux,
                                int *const iccvfg);

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

void CS_PROCF (cstime, CSTIME) (int    *const inpdt0,
                                int    *const iptlro,
                                int    *const ntmabs,
                                int    *const idtvar,
                                double *const dtref,
                                double *const dtmin,
                                double *const dtmax,
                                double *const coumax,
                                double *const foumax,
                                double *const varrdt,
                                double *const relxst);

/*----------------------------------------------------------------------------
 *
 * Fortran Interface:
 *
 * SUBROUTINE UINUM1
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (const    int *const isca,
                                const    int *const iscapp,
                                      double *const blencv,
                                         int *const ischcv,
                                         int *const isstpc,
                                         int *const ircflu,
                                      double *const cdtvar,
                                         int *const nitmax,
                                      double *const epsilo,
                                         int *const iresol,
                                         int *const imgrpr,
                                         int *const nswrsm);

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

void CS_PROCF (csnum2, CSNUM2) (   int *const ivisse,
                                double *const relaxp,
                                   int *const ipucou,
                                double *const extrag,
                                   int *const imrgra,
                                   int *const nterup);

void CS_PROCF (csphys, CSPHYS) (const    int *const nmodpp,
                                         int *const irovar,
                                         int *const ivivar,
                                         int *const icorio,
                                      double *const gx,
                                      double *const gy,
                                      double *const gz,
                                      double *const omegax,
                                      double *const omegay,
                                      double *const omegaz,
                                      double *const ro0,
                                      double *const viscl0,
                                      double *const viscv0,
                                      double *const visls0,
                                      double *const cp0,
                                      double *const t0,
                                      double *const p0,
                                      double *const xmasmr,
                                         int *const itempk);

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * subroutine cssca2
 * *****************
 *
 * integer          iscavr   <--  number of the related variance if any
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (const    int *const iscavr);

void CS_PROCF (cssca3, CSSCA3) (const    int *const itherm,
                                const    int *const iscalt,
                                const    int *const iscavr,
                                double *const visls0,
                                double *const t0,
                                double *const p0);

/*----------------------------------------------------------------------------
 * Array of properties used in the calculation
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprop, UIPROP) (const int *const irom,
                                const int *const iviscl,
                                const int *const ivisct,
                                const int *const ivisls,
                                const int *const icour,
                                const int *const ifour,
                                const int *const ismago,
                                const int *const iale,
                                const int *const icp,
                                const int *const iscalt,
                                const int *const iscavr,
                                const int *const iprtot,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icmome,
                                const int *const ipptx,
                                const int *const ippty,
                                const int *const ipptz,
                                const int *const ippdt,
                                const int *const ivisma,
                                const int *const idtvar,
                                const int *const ipucou,
                                const int *const iappel);

/*----------------------------------------------------------------------------
 * Temporal averaging treatment
 *----------------------------------------------------------------------------*/

void CS_PROCF (uimoyt, UIMOYT) (const int *const ndgmox,
                                      int *const ntdmom,
                                      int *const imoold,
                                      int *const idfmom);

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

void CS_PROCF (cstini, CSTINI) (double *const uref,
                                double *const almax);

void CS_PROCF(fcnmva, FCNMVA)
(
 const char      *const fstr,    /* Fortran string */
 int             *const len,     /* String Length  */
 int             *const ipp      /* Post Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

void CS_PROCF(cfnmva, CFNMVA)
(
 char          *const fstr,    /* Fortran string */
 int           *const len,     /* String Length  */
 int           *const ipp      /* Post Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

void CS_PROCF(nvamem, NVAMEM) (void);

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
 * subroutine uitssc (iscal, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          iscal    <-- index of the corresponding scalar
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitssc, UITSSC)(const int                  *iscal,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp);


/*----------------------------------------------------------------------------
 * Thermal scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsth (iscal, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          iscal    <-- index of the corresponding scalar
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsth, UITSTH)(const int                  *iscal,
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
 * integer          isca     <--  indirection array for scalar number
 * integer          iscold   <--  scalar number for restart
 * integer          iccfth   <--  type of initialisation(compressible model)
 * integer          ipr      <--  rtp index for pressure
 * integer          iscalt   <--  index of the thermal scalar
 * integer          itempk   <--  rtp index for temperature (in K)
 * integer          ienerg   <--  rtp index for energy total
 * double precision ro0      <--  value of density if IROVAR=0
 * double precision cp0      <--  value of specific heat if ICP=0
 * double precision viscl0   <--  value of viscosity if IVIVAR=0
 * double precision uref     <--  value of reference velocity
 * double precision almax    <--  value of reference length
 * double precision xyzcen   <--  cell's gravity center
 * double precision rtp      -->  variables and scalars array
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int          *ncelet,
                              const int          *isuite,
                              const int          *isca,
                              const int          *iscold,
                                    int          *iccfth,
                              const int *const    ipr,
                              const int *const    iscalt,
                              const int *const    itempk,
                              const int *const    ienerg,
                              const cs_real_t    *ro0,
                              const cs_real_t    *cp0,
                              const cs_real_t    *viscl0,
                              const cs_real_t    *uref,
                              const cs_real_t    *almax,
                              const double *const xyzcen,
                                    double        rtp[]);

/*----------------------------------------------------------------------------
 * User law for material Properties
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPHYV (NCELET, ISCA, RTP)
 * *****************
 *
 * INTEGER          NCEL     <--  number of cells whithout halo
 * INTEGER          NCELET   <--  number of cells whith halo
 * INTEGER          NSCAUS   <--  number of user scalar including thermal scalar
 * INTEGER          IVISCL   <--  pointer for mulecular viscosity mu
 * INTEGER          ICP      <--  pointer for predifined heat Cp
 * INTEGER          IVISLS   <--  pointer for Lambda/Cp
 * INTEGER          IROVAR   <--  =1 if rho variable, =0 if rho constant
 * INTEGER          IVIVAR   <--  =1 if mu variable, =0 if mu constant
 * INTEGER          ISCA     <--  indirection array for scalar number
 * INTEGER          ISCALT   <--  pointer for the thermal scalar in ISCA
 * INTEGER          ISCAVR   <--  scalars that are variance
 * INTEGER          IPPROC   <--  indirection array for cell properties
 * INTEGER          IVISCV   <--  pointer for volumic viscosity viscv
 * INTEGER          ITEMPK   <--  pointer for temperature (in K)
 * DOUBLE PRECISION P0       <--  pressure reference value
 * DOUBLE PRECISION T0       <--  temperature reference value
 * DOUBLE PRECISION RO0      <--  density reference value
 * DOUBLE PRECISION CP0      <--  specific heat reference value
 * DOUBLE PRECISION VISCL0   <--  dynamic viscosity reference value
 * DOUBLE PRECISION VISLS0   <--  diffusion coefficient of the scalars
 * DOUBLE PRECISION VISCV0   <--  volumic viscosity
 * DOUBLE PRECISION RTP      <--  variables and scalars array
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const cs_int_t  *const ncel,
                              const cs_int_t  *const ncelet,
                              const cs_int_t  *const nscaus,
                              const cs_int_t  *const itherm,
                              const cs_int_t         iviscl[],
                              const cs_int_t         icp[],
                              const cs_int_t         ivisls[],
                              const cs_int_t         irovar[],
                              const cs_int_t         ivivar[],
                              const cs_int_t         isca[],
                              const cs_int_t         iscalt[],
                              const cs_int_t         iscavr[],
                              const cs_int_t         ipproc[],
                              const cs_int_t         iviscv[],
                              const cs_int_t         itempk[],
                              const cs_real_t        p0[],
                              const cs_real_t        t0[],
                              const cs_real_t        ro0[],
                              const cs_real_t        cp0[],
                              const cs_real_t        viscl0[],
                              const cs_real_t        visls0[],
                              const cs_real_t        viscv0[],
                              const cs_real_t        rtp[],
                                    double           propce[]);

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
 * double precision rtpa     <--  variables array at previous time step
 *----------------------------------------------------------------------------*/

void CS_PROCF(uikpdc, UIKPDC)(const int*   iappel,
                              const int*   ncelet,
                                    int*   ncepdp,
                                    int    icepdc[],
                                    double ckupdc[],
                              const double rtpa[] );

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
 * DOUBLE PRECISION RTP      <--  variables and scalars array
 * DOUBLE PRECISION PROPCE   <--  property array
 * INTEGER          IPPROC   <--  indirection array for cell properties
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprof, UIPROF)(const int    *const ncelet,
                               const int    *const ncel,
                               const int    *const ntmabs,
                               const int    *const ntcabs,
                               const double *const ttcabs,
                               const double *const ttmabs,
                               const double *const ttpabs,
                               const double *const xyzcen,
                               const double *const rtp,
                               const double *const propce,
                               const int    *const ipproc);

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

void CS_PROCF (memui1, MEMUI1) (const int *const ncharb);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

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
 * Set partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_gui_partition(void);

/*-----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_parallel_io(void);

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
