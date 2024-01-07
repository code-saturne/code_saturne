#ifndef __CS_COAL_H__
#define __CS_COAL_H__

/*============================================================================
 * Coal combustion model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*! Maximum number of coals */
#define  CS_COMBUSTION_MAX_COALS  5

/*! Maximum number of oxydants */
#define CS_COMBUSTION_COAL_MAX_OXYDANTS 3

/*! Maximum number of coal classes per coal */
#define  CS_COMBUSTION_MAX_CLASSES_PER_COAL  20

/*! Maximum total number of coal classes */
#define  CS_COMBUSTION_COAL_MAX_CLASSES    CS_COMBUSTION_MAX_COALS \
                                         * CS_COMBUSTION_MAX_CLASSES_PER_COAL

/*! Maximum number of elementary gas components */
#define  CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS  20

/*! Maximum number of tabulation points */
#define  CS_COMBUSTION_COAL_MAX_TABULATION_POINTS  500

/*! Maximum number of solid constituants */
#define  CS_COMBUSTION_COAL_MAX_SOLIDS  CS_COMBUSTION_MAX_COALS * 4

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Coal combustion model type */
/*------------------------------*/

typedef enum {

  CS_COMBUSTION_COAL_NONE = -1,

  CS_COMBUSTION_COAL_STANDARD = 0,
  CS_COMBUSTION_COAL_WITH_DRYING = 1

} cs_coal_model_type_t;

/*! Coal combustion model parameters structure */
/*---------------------------------------------*/

typedef struct {

  /* Generic members
     ---------------
     (keep aligned with gas combustion, so that we can
     move to an object inheritance model in the future) */

  int     n_gas_el_comp;             /*!< number of elementary gas components */
  int     n_gas_species;             /*!< number of global species */
  int     n_atomic_species;          /*!< number of atomic species */

  int     n_reactions;               /*!< number of global reactions
                                      *   in gas phase */

  int     n_tab_points;              /*!< number of tabulation points */

  double  pcigas;                    /*!< combustible reaction enthalpy
                                       (Lower Calorific Value)*/
  double  xco2;                      /*!< molar coefficient of CO2 */
  double  xh2o;                      /*!< molar coefficient of H2O */

  /*! molar mass of an elementary gas component */
  double  wmole[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

  /*! composition of oxidants in O2 */
   double oxyo2[CS_COMBUSTION_COAL_MAX_OXYDANTS];

  /*! composition of N2 oxidants */
  double oxyn2[CS_COMBUSTION_COAL_MAX_OXYDANTS];

  /*! composition of H2O oxidants */
  double oxyh2o[CS_COMBUSTION_COAL_MAX_OXYDANTS];

  /*! composition of CO2 oxidants */
  double oxyco2[CS_COMBUSTION_COAL_MAX_OXYDANTS];

  /*! temperature in K */
  double th[CS_COMBUSTION_COAL_MAX_TABULATION_POINTS];

  /* Members specific to the coal combustion model
     --------------------------------------------- */

  cs_coal_model_type_t  type;  /*!< combustion model type */

  int     n_coals;     /*!< number of coal types */
  int     nclacp;      /*!< number of coal classes */

  int     nsolim;      /*!< number of solid constituents
                         (reactive coal, coke, ash) */

  int     idrift;      /*!< drift (0: off, 1: on) */

  int     ieqnox;      /*!< NOx model (0: off; 1: on) */

  int     ieqco2;      /*!< kinetic model for CO <=> CO2
                         - 0  unused (maximal conversion
                                      in turbulent model)
                         - 1  transport of CO2 mass fraction
                         - 2  transport of CO mass fraction  */

  double  ckabs0;      /*!< absorption coefficient of gas mix */

  /*! number of classes per coal */
  int     n_classes_per_coal[CS_COMBUSTION_MAX_COALS];

  /* Properties of dry coal
     ---------------------- */

  /*! elementary composition of coal in C over dry (%) */
  double cch[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coal in H over dry (%) */
  double hch[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coal in O over dry (%) */
  double och[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coal in S over dry (%) */
  double sch[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coal in N over dry (%) */
  double nch[CS_COMBUSTION_MAX_COALS];

  /*! composition of reactive coal: alpha(c) = hch(c)/cch(c) */
  double alpha[CS_COMBUSTION_MAX_COALS];

  /*! composition of reactive coal: beta (c) = och(c)/cch(c) */
  double beta[CS_COMBUSTION_MAX_COALS];

  /*! composition of reactive coal: teta (c) = sch(c)/cch(c) */
  double teta[CS_COMBUSTION_MAX_COALS];

  /*! composition of reactive coal: omega (c) = nch(c)/cch(c) */
  double omega[CS_COMBUSTION_MAX_COALS];

  /*! coal pci (J/kg) */
  double pcich[CS_COMBUSTION_MAX_COALS];

  /*! initial density (kg/m3) */
  double rho0ch[CS_COMBUSTION_MAX_COALS];

  /*! coal thermal conductivity (W/m/K) */
  double thcdch[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coke in C over dry (%) */
  double cck[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coke in H over dry (%) */
  double hck[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coke in O over dry (%) */
  double ock[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coke in S over dry (%) */
  double sck[CS_COMBUSTION_MAX_COALS];

  /*! elementary composition of coke in N over dry (%) */
  double nck[CS_COMBUSTION_MAX_COALS];

  /*! composition of coke: gamma(c) = hck(c)/cck(c) */
  double gamma[CS_COMBUSTION_MAX_COALS];

  /*! composition of coke: delta(c) = ock(c)/cck(c) */
  double delta[CS_COMBUSTION_MAX_COALS];

  /*! composition of coke: kappa(c) = sck(c)/cck(c) */
  double kappa[CS_COMBUSTION_MAX_COALS];

  /*! composition of coke: zeta(c) = nck(c)/cck(c) */
  double zeta[CS_COMBUSTION_MAX_COALS];

  /*! coke pci (J/kg) */
  double pcick[CS_COMBUSTION_MAX_COALS];

  /*! coke density (kg/m3) */
  double rhock[CS_COMBUSTION_MAX_COALS];

  /*! Cp of ash (J/kg/K) */
  double cpashc[CS_COMBUSTION_MAX_COALS];

  /*! enthalpy of ash formation (J/kg) */
  double h0ashc[CS_COMBUSTION_MAX_COALS];

  /*! H0 of coal */
  double h02ch[CS_COMBUSTION_MAX_COALS];

  /*! Cp of water */
  double cp2wat[CS_COMBUSTION_MAX_COALS];

  /*! distribution of N2 in HCN and No reaction 1 */
  double  crepn1[CS_COMBUSTION_MAX_COALS][2];

  /*! distribution of N2 in HCN and No reaction 2 */
  double  crepn2[CS_COMBUSTION_MAX_COALS][2];

  /*! coal specific heat */
  double  cp2ch[CS_COMBUSTION_MAX_COALS];

  /*! Ash fraction (kg/kg) in percentage */
  double xashsec[CS_COMBUSTION_MAX_COALS];

  /*! ashes concentration (kg/kg) */
  double  xashch[CS_COMBUSTION_MAX_COALS];

  /*! humidity (kg/kg) */
  double  xwatch[CS_COMBUSTION_MAX_COALS];

  /* Kinetic parameters for devolatilization (Kobayashi's model)
     ----------------------------------------------------------- */

  /*! Indicator: 0 if MVl = {CH4;CO};  1 if MVl = {CHz;CO} */
  int iy1ch [CS_COMBUSTION_MAX_COALS];

  /*! Indicator 0 if MVL = {C2H4;CO};  1 if MVL = {CxHy;CO} */
  int iy2ch [CS_COMBUSTION_MAX_COALS];

  /*! Order of the reaction for heterogeneous coke/O2 combustion
    (0.5 if 0, 1 if 1) */
  int iochet [CS_COMBUSTION_MAX_COALS];

  /*! Order of the reaction for heterogeneous coke/CO2 combustion
    (0.5 if 0, 1 if 1) */
  int ioetc2[CS_COMBUSTION_MAX_COALS];

  /*! Order of the reaction for heterogeneous coke/H2O combustion
    (0.5 if 0, 1 if 1) */
  int ioetwt[CS_COMBUSTION_MAX_COALS];

  /*! stoechiometric coeffficient; computed if iy1ch = 0; given if iy1ch = 1 */
  double  y1ch[CS_COMBUSTION_MAX_COALS];

  /*! pre-exponetial factor (1/s) */
  double  a1ch[CS_COMBUSTION_MAX_COALS];

  /*! activation energy (J/mol) */
  double  e1ch[CS_COMBUSTION_MAX_COALS];

  /*! stoechiometric coeffficient; computed if iy2ch = 0; given if iy2ch = 1 */
  double  y2ch[CS_COMBUSTION_MAX_COALS];

  /*! pre-exponetial factor (1/s) */
  double  a2ch[CS_COMBUSTION_MAX_COALS];

  /*! activation energy (J/mol) */
  double  e2ch[CS_COMBUSTION_MAX_COALS];

  /* Kinetics  of heterogeneous coke combustion (shrinking sphere model)
     ------------------------------------------------------------------ */

  /*! pre-exponential constant for combustion of coke with O2 (kg/m2/s/atm) */
  double  ahetch[CS_COMBUSTION_MAX_COALS];

  /*! activation energy for combustion of coke with O2 (kcal/mol) */
  double  ehetch[CS_COMBUSTION_MAX_COALS];

  /*! pre-exponential constant for combustion of coke with CO2 (kg/m2/s/atm) */
  double  ahetc2[CS_COMBUSTION_MAX_COALS];

  /*! activation energy for combustion of coke with CO2 (kcal/mol) */
  double  ehetc2[CS_COMBUSTION_MAX_COALS];

  /*! pre-exponential constant for combustion of coke with H2O (kg/m2/s/atm) */
  double  ahetwt[CS_COMBUSTION_MAX_COALS];

  /*! activation energy for combustion of coke with H2O (kcal/mol) */
  double  ehetwt[CS_COMBUSTION_MAX_COALS];

  /* Enthalpy of reactive coal, coke, and ash
     ---------------------------------------- */

  /*! position in ehsoli array for reactive coal */
  int     ich[CS_COMBUSTION_MAX_COALS];

  /*! position in ehsoli array for coke */
  int     ick[CS_COMBUSTION_MAX_COALS];

  /*! position in ehsoli array for ash */
  int     iash[CS_COMBUSTION_MAX_COALS];

  /*! position in ehsoli array for humidity */
  int     iwat[CS_COMBUSTION_MAX_COALS];

  /*! mass enthalpy (J/kg) at temperature T of solid constituent S */
  double  ehsoli[CS_COMBUSTION_COAL_MAX_TABULATION_POINTS]
                [CS_COMBUSTION_COAL_MAX_SOLIDS];

  /*! molar mass of solid constituents */
  double  wmols[CS_COMBUSTION_COAL_MAX_SOLIDS];

  /*! formation enthalpy (J/kg) of solid constituents */
  double  eh0sol[CS_COMBUSTION_COAL_MAX_SOLIDS];

  /* By class (deduced properties)
     ----------------------------- */

  /*! coal id if considered class belongs to coal ich[1, 2, ...] */
  int     ichcor[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! initial diameter (m) */
  double  diam20[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! minimum diameter (m) */
  double  dia2mn[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! initial density (kg/m^3) */
  double  rho20[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! minimal density (kg/m^3) */
  double  rho2mn[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! initial particle mass (kg) */
  double  xmp0[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! particle ashes mass (kg) */
  double  xmasch[CS_COMBUSTION_COAL_MAX_CLASSES];

  /* Data relative to combustion of gaseous species
     ---------------------------------------------- */

  int     ico;         /*!< index of co in wmole */
  int     io2;         /*!< index of o2 in wmole */
  int     ih2o;        /*!< index of h2o in wmole */
  int     in2;         /*!< index of n2 in wmole */
  int     ico2;        /*!< index of co2 in wmole */

  /*! index of CHx1 in ehgaze and wmole */
  int ichx1c[CS_COMBUSTION_MAX_COALS];

  /*! index of CHx2 in ehgaze and wmole */
  int ichx2c[CS_COMBUSTION_MAX_COALS];

  /*! index of CHx1m in ehgaze and wmole */
  int ichx1;

  /*! index of CHx2m in ehgaze and wmole */
  int ichx2;

  /*! Composition of hydrocarbon relative to MVl: CH(X1) */
  double chx1[CS_COMBUSTION_MAX_COALS];

  /*! Composition of hydrocarbon relative to MVl: CH(X2) */
  double chx2[CS_COMBUSTION_MAX_COALS];

  //@{
  /*! low T devolatilization molar stoechiometric coefficients */
  double a1[CS_COMBUSTION_MAX_COALS];
  double b1[CS_COMBUSTION_MAX_COALS];
  double c1[CS_COMBUSTION_MAX_COALS];
  double d1[CS_COMBUSTION_MAX_COALS];
  double e1[CS_COMBUSTION_MAX_COALS];
  double f1[CS_COMBUSTION_MAX_COALS];
  double a2[CS_COMBUSTION_MAX_COALS];
  double b2[CS_COMBUSTION_MAX_COALS];
  double c2[CS_COMBUSTION_MAX_COALS];
  double d2[CS_COMBUSTION_MAX_COALS];
  double e2[CS_COMBUSTION_MAX_COALS];
  double f2[CS_COMBUSTION_MAX_COALS];
  //@}

  /* Complement table
     ---------------- */

  /*! temperature values in enthalpy/temperature law tabulation */
  double thc[CS_COMBUSTION_COAL_MAX_TABULATION_POINTS];

  /*! number of tabulation points for enthalpy/temperature law */
  int npoc;

  /* Combustion of coke by H2O
     ------------------------- */

  /*! mass transfer by heterogeneous combustion with H2O */
  int ihth2o;

  /*! id of field "het_ts_h2o_p_<class>" */
  int ighh2o[CS_COMBUSTION_COAL_MAX_CLASSES];

  /*! PCI computation mode:
      - 1: dry -> pure (schaff's formula)
      - 2: raw -> pure
      - 3: pure -> pure
      - 4: dry -> pure
      - 5: raw -> pure
      - 6: IGT correclation */
  int ipci[CS_COMBUSTION_MAX_COALS];

  /* NOx model
     --------- */

  /*! percentage of Nitrogen freed in devolatilization. */
  double qpr[CS_COMBUSTION_MAX_COALS];

  /*! concentration in Nitrogen relative to pure */
  double fn[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of HCN in light volatile matters. */
  double yhcnle[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of HCN in heavy volatile matters. */
  double yhcnlo[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of NH3 in light volatile matters. */
  double ynh3le[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of NH3 in heavy volatile matters. */
  double ynh3lo[CS_COMBUSTION_MAX_COALS];

  /*! Percentage of total N in coal of char1 */
  double repnle[CS_COMBUSTION_MAX_COALS];

  /*! Percentage of total N in coal of char2 */
  double repnlo[CS_COMBUSTION_MAX_COALS];

  /*! Percentage of HCN produced by heteorgeneous combustion */
  double repnck[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of HCN in products of heteorgeneous combustion of char 1 */
  double yhcnc1[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of NO in products of heteorgeneous combustion of char 1 */
  double ynoch1[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of HCN in products of heteorgeneous combustion of char 2 */
  double yhcnc2[CS_COMBUSTION_MAX_COALS];

  /*! mass fraction of NO in products of heteorgeneous combustion of char 2 */
  double ynoch2[CS_COMBUSTION_MAX_COALS];

  /*! molar mass of CHx1 */
  double wmchx1;

  /*! molar mass of CHx2 */
  double wmchx2;

} cs_coal_model_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*! Combustion model parameters structure */

extern cs_coal_model_t  *cs_glob_coal_model;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Return coal combustion model type.
 *
 * \return type of active coal combustion model
 *         (CS_COMBUSTION_COAL_NONE if model not active)
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_type_t
cs_coal_model_get_type(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Activate coal combustion model.
 *
 * \return  pointer to coal combustion model structure.
 *
 * \param[in]  type  coal combustion model type
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_t *
cs_coal_model_set_model(cs_coal_model_type_t  type);

/*----------------------------------------------------------------------------*/
/*
 * \brief Take in account the radiative source terms in the particle equation
 *        of a given class for pulverized coal flame.
 *
 * \param[in]      f       pointer to scalar field
 * \param[in, out] smbrs   right and side (explicit ST contribution)
 * \param[in, out] rovsdt  system diagonal (implicit ST contribution)
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_rad_transfer_st(const cs_field_t  *f,
                        cs_real_t         *smbrs,
                        cs_real_t         *rovsdt);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute gas enthalpy
 *        Function with gas temperature and concentrations
 *
 * \param[in] xesp      mass fraction of species
 * \param[in] f1mc      average f1
 * \param[in] f2mc      average f2
 * \param[in] tp        gas temperature (in kelvin)
 *
 * \return   gas enthalpy (in \f$ j . kg^{-1}) \f$ of mixed gas
 */
/*----------------------------------------------------------------------------*/

double
cs_coal_thconvers1(const double  xesp[],
                   const double  f1mc[],
                   const double  f2mc[],
                   cs_real_t     tp);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute particles enthalpy
 *         Function with temperature and concentrations

 * \param[in]     class_id      class id
 * \param[in]     xsolid        mass fraction of components
 * \param[in,out] temper        temperature (in kelvin)

 * \return   mass enthalpy (in \f$ j . kg^{-1}) \f$
 */
/*----------------------------------------------------------------------------*/

double
cs_coal_thconvers2(int           class_id,
                   const double  xsolid[],
                   cs_real_t     temper);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COAL__ */
