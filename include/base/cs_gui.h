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

#ifndef __CS_GUI_H__
#define __CS_GUI_H__

/*============================================================================
 * Reader of the parameters file: main parameters, boundary conditions
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
 * Type definitions
 *============================================================================*/

typedef enum {
  DIRICHLET,
  FLOW1,
  HYDRAULIC_DIAMETER,
  TURBULENT_INTENSITY,
  NEUMANN,
  COEF_ECHANGE,
  COALFLOW,
  WALL_FUNCTION
} cs_boundary_value_t;


/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *----------------------------------------------------------------------------*/

void cs_gui_clean_memory(void);

/*-----------------------------------------------------------------------------
 * Return the name of a thermophysical model.
 *
 * parameter:
 *   model_thermo          -->  thermophysical model
 *----------------------------------------------------------------------------*/

char *cs_gui_get_thermophysical_model(const char *const model_thermo);

/*-----------------------------------------------------------------------------
 * Return if a particular physics model is activated.
 *----------------------------------------------------------------------------*/

int cs_gui_get_activ_thermophysical_model(void);

/*-----------------------------------------------------------------------------
 * Return number of boundary regions definition
 *----------------------------------------------------------------------------*/

int cs_gui_boundary_zones_number(void);

/*-----------------------------------------------------------------------------
 * Return the nature of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_nature(const int ith_zone);

/*-----------------------------------------------------------------------------
 * Return the label of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_label(const int ith_zone);


/*-----------------------------------------------------------------------------
 * Return the zone number of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

int cs_gui_boundary_zone_number(const int ith_zone);

/*-----------------------------------------------------------------------------
 * Return the description of a boundary zone
 *
 * parameters:
 *   nature                -->  nature of boundary zone (inlet, wall,...)
 *   label                 -->  label of boundary zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_localization(const char *const nature,
                                        const char *const label);

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Turbulence model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTURB (ITURB, IDEUCH, IGRAKE, IGRAKI, XLOMLG)
 * *****************
 *
 * INTEGER          ITURB   <--   turbulence model
 * INTEGER          IDEUCH  <--   wall law treatment
 * INTEGER          IGRAKE  <--   k-eps gravity effects
 * INTEGER          IGRAKI  <--   Rij-eps gravity effects
 * DOUBLE PRECISION XLOMLG  <--   mixing_length_scale
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
 * SUBROUTINE CSCPVA (ICP)
 * *****************
 *
 * INTEGER          ICP     <--   Specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscpva, CSCPVA) (int *const icp);

/*----------------------------------------------------------------------------
 * User scalars number.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNSCA (NSCAUS)
 * *****************
 *
 * INTEGER          NSCAUS     <--   user scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnsca, CSNSCA) (int *const nscaus);

/*-----------------------------------------------------------------------------
 * Predefined physics indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPPMO
 * *****************
 *
 * INTEGER          IPPMOD <--  predefined physics indicator array
 * INTEGER          ICOD3P  --> diffusion flame en chimie complete rapide
 * INTEGER          ICODEQ  --> diffusion flame en chimie rapide vers l'equilibre
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
                                int *const indjon,
                                int *const ieqco2);

/*----------------------------------------------------------------------------
 * User scalars which are variance.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISCA (ISCAVR)
 * *****************
 *
 * INTEGER          ISCAVR     <--   user scalars variance array
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisca, CSISCA) (int *const iscavr);

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIVIS (IDTVAR)
 * *****************
 *
 * INTEGER          ISCAVR  <-->  number of the related variance if any
 * INTEGER          IVISLS  <--   indicator for the user scalar viscosity
 * INTEGER          ISCALT  <-->  number of the user thermal scalar if any
 * INTEGER          ISCSTH  <-->  type of the user thermal scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (int *const iscavr,
                                int *const ivisls,
                                int *const iscalt,
                                int *const iscsth);

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIDTV (IDTVAR)
 * *****************
 *
 * INTEGER          IDTVAR  <--   fixed or variable time step
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
 * INTEGER          IPHYDR  <--   hydrostatic pressure
 *----------------------------------------------------------------------------*/

void CS_PROCF (csiphy, CSIPHY) (int *const iphydr);


/*----------------------------------------------------------------------------
 * Est appele juste avant le 3eme appel a VARPOS
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALIN()
 * *****************
 *
 * INTEGER          IALE    <--   iale method activation
 * INTEGER          NALINF  <--   number of sub iteration of initialization of fluid
 * INTEGER          NALIMX  <--   max number of iterations of implicitation of
 *                                the displacement of the structures
 * DOUBLE           EPALIM  <--   realtive precision of implicitation of
 *                                the displacement of the structures
 * INTEGER          IORTVM  <--   type of viscosity of mesh
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialin, UIALIN) (int    *const iale,
                                int    *const nalinf,
                                int    *const nalimx,
                                double *const epalim,
                                int    *const iortvm);

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
                                const int *const iale,
                                const int *const iuma,
                                const int *const ivma,
                                const int *const iwma,
                                const int *const isca,
                                const int *const iscapp);

/*----------------------------------------------------------------------------
 * Restart files format.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIFOA (IFOAVA, IFOAVX)
 * *****************
 *
 * INTEGER          IFOAVA  <--   main restart file format
 * INTEGER          IFOAVX  <--   auxiliary restart file format
 *----------------------------------------------------------------------------*/


void CS_PROCF (csifoa, CSIFOA) (int *const ifoava,
                                int *const ifoavx);

/*----------------------------------------------------------------------------
 * Restart parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISUI (ISUITE, ILEAUX, ICCVFG)
 * *****************
 *
 * INTEGER          ISUITE  <--   restart
 * INTEGER          ILEAUX  <--   restart with auxiliary
 * INTEGER          ICCFVG  <--   restart with frozen field
 *----------------------------------------------------------------------------*/


void CS_PROCF (csisui, CSISUI) (int *const isuite,
                                int *const ileaux,
                                int *const iccvfg);

/*----------------------------------------------------------------------------
 * Time passing parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTIME (INPDT0, IPTLTO, NTMABS, DTREF,
 * *****************  DTMIN, DTMAX, COUMAX, FOUMAX, VARRDT)
 *
 * INTEGER          INPDT0  <--   zero tim step
 * INTEGER          IPTLTO  <--   thermal time step control
 * INTEGER          NTMABS  <--   iterations numbers
 * INTEGER          IDTVAR  <--   time step's options
 * DOUBLE PRECISION DTREF   <--   time step
 * DOUBLE PRECISION DTMIN   <--   minimal time step
 * DOUBLE PRECISION DTMAX   <--   maximal time step
 * DOUBLE PRECISION COUMAX  <--   maximal courant number
 * DOUBLE PRECISION FOUMAX  <--   maximal fournier number
 * DOUBLE PRECISION VARRDT  <--   max time step variation between 2 iterations
 * DOUBLE PRECISION RELXST  <--   relaxation coefficient if idtvar = -1
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
 * Check if a users thermal scalar is defined.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSSCA1 (ISCALT, ISCSTH)
 * *****************
 *
 * INTEGER          ISCALT  <--   number of the user thermal scalar if any
 * INTEGER          ISCSTH  <--   type of the user thermal scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca1, CSSCA1) (int *const iscalt,
                                int *const iscsth);


void CS_PROCF (uinum1, UINUM1) (const    int *const isca,
                                const    int *const iscapp,
                                      double *const blencv,
                                         int *const ischcv,
                                         int *const isstpc,
                                         int *const ircflu,
                                      double *const cdtvar,
                                         int *const nitmax,
                                      double *const epsilo);

/*----------------------------------------------------------------------------
 * Global numerical parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNUM2 (IVISSE, RELAXP, IPUCOU, EXTRAG, IMRGRA)
 * *****************
 *
 * INTEGER          IVISSE  <--   gradient transpose
 * INTEGER          RELAXP  <--   pressure relaxation
 * INTEGER          IPUCOU  <--   velocity pressure coupling
 * INTEGER          EXTRAG  <--   wall pressure extrapolation
 * INTEGER          IMRGRA  <--   gradient reconstruction
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2) (   int *const ivisse,
                                double *const relaxp,
                                   int *const ipucou,
                                double *const extrag,
                                   int *const imrgra);

void CS_PROCF (csphys, CSPHYS) (const    int *const nmodpp,
                                         int *const irovar,
                                         int *const ivivar,
                                      double *const gx,
                                      double *const gy,
                                      double *const gz,
                                      double *const ro0,
                                      double *const viscl0,
                                      double *const cp0,
                                      double *const t0,
                                      double *const p0);

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSSCA2 (ISCAVR, SCAMIN, SCAMAX)
 * *****************
 *
 * INTEGER          ISCAVR   -->  number of the related variance if any
 * DOUBLE PRECISION SCAMIN  <--   user scalar min array
 * DOUBLE PRECISION SCAMAX  <--   user scalar max array
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) ( const    int *const iscavr,
                                       double *const scamin,
                                       double *const scamax);

void CS_PROCF (cssca3, CSSCA3) (const    int *const iscalt,
                                const    int *const iscavr,
                                      double *const visls0,
                                      double *const t0,
                                      double *const p0);


/*----------------------------------------------------------------------------
 * Tableau des propriétés utilisées dans le calcul
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
                                const int *const iappel);

/*----------------------------------------------------------------------------
 * Traitement des moyennes temporelles
 *----------------------------------------------------------------------------*/

void CS_PROCF (uimoyt, UIMOYT) (const int *const ndgmox,
                                const int *const isca,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icmome,
                                      int *const ntdmom,
                                      int *const imoold,
                                      int *const idfmom);

/*----------------------------------------------------------------------------
 * Turbulence initialization parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTINI (UREF, ALMAX)
 * *****************
 *
 * INTEGER          UREF   <--   reference velocity
 * INTEGER          ALMAX  <--   reference length
 *----------------------------------------------------------------------------*/


void CS_PROCF (cstini, CSTINI) (double *const uref,
                                double *const almax);

void CS_PROCF (csenso, CSENSO) (const    int *const nvppmx,
                                         int *const ncapt,
                                         int *const nthist,
                                         int *const ntlist,
                                         int *const ichrvl,
                                         int *const ichrbo,
                                         int *const ichrsy,
                                         int *const ichrmd,
                                        char *const fmtchr,
                                         int *const size_fmt,
                                        char *const optchr,
                                         int *const size_opt,
                                         int *const ntchr,
                                         int *const iecaux,
                                         int *const ichrvr,
                                         int *const ilisvr,
                                         int *const ihisvr,
                                const    int *const isca,
                                const    int *const iscapp,
                                const    int *const ipprtp,
                                const    int *const ipppro,
                                const    int *const ipproc,
                                      double *const xyzcap);

void CS_PROCF(fcnmva, FCNMVA)
(
 const char      *const fstr,    /* --> Fortran string */
 int             *const len,     /* --> String Length  */
 int             *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

void CS_PROCF(cfnmva, CFNMVA)
(
 char          *const fstr,    /* --> Fortran string */
 int           *const len,     /* --> String Length  */
 int           *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Users arrays
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIUSAR (ICOFTU)
 * *****************
 *
 * INTEGER          ICOFTU   -->  Dimension coef for user arrays
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiusar, UIUSAR) (int *const icoftu);

/*----------------------------------------------------------------------------
 * Variables and user scalars initialization
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIINIV (NCELET, ISCA, RTP)
 * *****************
 *
 * INTEGER          ISCAVR   -->  number of cells
 * INTEGER          ISCA     -->  indirection array for scalar number
 * DOUBLE PRECISION RTP     <--   variables and scalars array
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV) (const int    *const ncelet,
                               const int    *const isca,
                                     double *const rtp);

/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLIM
 * *****************
 *
 * INTEGER          NOZPPM  --> max number of boundary conditions zone
 * INTEGER          NFABOR  --> number of boundary faces
 * INTEGER          IINDEF  --> type of boundary: not defined
 * INTEGER          IENTRE  --> type of boundary: inlet
 * INTEGER          IPAROI  --> type of boundary: smooth wall
 * INTEGER          IPARUG  --> type of boundary: rough wall
 * INTEGER          ISYMET  --> type of boundary: symmetry
 * INTEGER          ISOLIB  --> type of boundary: outlet
 * INTEGER          IQIMP   --> 1 if flow rate is applied
 * INTEGER          ICALKE  --> 1 for automatic turbulent boundary conditions
 * INTEGER          ITYPFB  --> type of boundary for each face
 * INTEGER          IZFPPP  --> zone number
 * INTEGER          ICODCL  --> boundary conditions array type
 * DOUBLE PRECISION QIMP    --> flow rate value if applied
 * DOUBLE PRECISION DH      --> hydraulic diameter
 * DOUBLE PRECISION XINTUR  --> turbulent intensity
 * DOUBLE PRECISION RCODCL  --> boundary conditions array value
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM) (const    int *const nozppm,
                                const    int *const nfabor,
                                const    int *const iindef,
                                const    int *const ientre,
                                const    int *const iparoi,
                                const    int *const iparug,
                                const    int *const isymet,
                                const    int *const isolib,
                                         int *const iqimp,
                                         int *const icalke,
                                         int *const itypfb,
                                         int *const izfppp,
                                         int *const icodcl,
                                      double *const qimp,
                                      double *const dh,
                                      double *const xintur,
                                      double *const rcodcl);


void CS_PROCF (uicpcl, UICPCL) (const    int *const nozppm,
                                const    int *const ncharm,
                                const    int *const ncharb,
                                const    int *const nclpch,
                                const    int *const nfabor,
                                const    int *const iindef,
                                const    int *const ientre,
                                const    int *const iparoi,
                                const    int *const iparug,
                                const    int *const isymet,
                                const    int *const isolib,
                                         int *const itypfb,
                                         int *const icodcl,
                                      double *const rcodcl,
                                      double *const surfbo,
                                         int *const ientat,
                                         int *const iqimp,
                                      double *const qimpat,
                                      double *const timpat,
                                         int *const ientcp,
                                      double *const qimpcp,
                                      double *const timpcp,
                                      double *const distch,
                                         int *const icalke,
                                      double *const dh,
                                      double *const xintur,
                                         int *const izfppp);


/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *
 * INTEGER          NFABOR  --> number of boundary faces
 * INTEGER          IINDEF  --> type of boundary: not defined
 * INTEGER          IENTRE  --> type of boundary: inlet
 * INTEGER          IPAROI  --> type of boundary: wall
 * INTEGER          IPARUG  --> type of boundary: wall with rugosity
 * INTEGER          ISYMET  --> type of boundary: symmetry
 * INTEGER          ISOLIB  --> type of boundary: outlet
 * INTEGER          ITYPFB  --> type of boundary for each face
 * INTEGER          IZFPPP  --> zone number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE) (const int *const nfabor,
                                const int *const iindef,
                                const int *const ientre,
                                const int *const iparoi,
                                const int *const iparug,
                                const int *const isymet,
                                const int *const isolib,
                                      int *const itypfb,
                                      int *const izfppp);

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
 * Defintion des pointeurs des scalaires model pour le charbon
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpsc, UICPSC) (const int *const ncharb,
                                const int *const nclass,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ieqco2,
                                const int *const ihm,
                                const int *const inp,
                                const int *const ixch,
                                const int *const ixck,
                                const int *const ixwt,
                                const int *const ih2,
                                const int *const if1m,
                                const int *const if2m,
                                const int *const if3m,
                                const int *const if4p2m,
                                const int *const if5m,
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
                                const int *const igmsec,
                                const int *const ilumi);

/*----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables.
 *
 * Fortran Interface:
 *
 * SUBROUTINE MEMUI1
 * *****************
 *
 * INTEGER          NCHARB  --> number of coal
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui1, MEMUI1) (const int *const ncharb);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_GUI_H__ */
