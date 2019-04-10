#ifndef __CS_LAGR_H__
#define __CS_LAGR_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "assert.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_field.h"

#include "cs_lagr_injection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for computation of particle injection profile.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * \param[in]   zone_id      id of associated mesh zone
 * \param[in]   location_id  id of associated mesh location
 * \param[in]   input        pointer to optional (untyped) value or structure.
 * \param[in]   n_elts       number of zone elements
 * \param[in]   elt_ids      ids of zone elements
 * \param[out]  profile      weight of a given zone element (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_lagr_injection_profile_compute_t) (int               zone_id,
                                       int               location_id,
                                       const void       *input,
                                       cs_lnum_t         n_elts,
                                       const cs_lnum_t   elt_ids[],
                                       cs_real_t         profile[]);

/*! Lagrangian boundary condition types */
/*--------------------------------------*/

typedef enum {

  CS_LAGR_BC_UNDEFINED,
  CS_LAGR_INLET,
  CS_LAGR_OUTLET,
  CS_LAGR_REBOUND,
  CS_LAGR_DEPO1,
  CS_LAGR_DEPO2,
  CS_LAGR_FOULING,
  CS_LAGR_JBORD1,
  CS_LAGR_JBORD2,
  CS_LAGR_JBORD3,
  CS_LAGR_JBORD4,
  CS_LAGR_JBORD5,
  CS_LAGR_DEPO_DLVO,
  CS_LAGR_SYM

} cs_lagr_bc_type_t;

/*! Lagrangian deposition state */

typedef enum {

  CS_LAGR_PART_IN_FLOW        = 0,
  CS_LAGR_PART_DEPOSITED      = 1,
  CS_LAGR_PART_ROLLING        = 2,
  CS_LAGR_PART_TO_DELETE      = 3,
  CS_LAGR_PART_NO_MOTION      = 10,
  CS_LAGR_PART_IMPOSED_MOTION = 11

} cs_lagr_deposition_state_t;

  /*! Lagrangian module status.
     the different values correspond to the following coupling:
     - CS_LAGR_OFF: Lagrangian module off
     - CS_LAGR_ONEWAY_COUPLING: Lagrangian two-phase flow in one-way coupling
         (no influence of the particles on the continuous phase)
     - CS_LAGR_TWOWAY_COUPLING: Lagrangian two-phase flow with two-way coupling
         (influence of the particles on the dynamics of the continuous phase).
         Dynamics, temperature and mass may be coupled independently.
     - CS_LAGR_FROZEN_CONTINUOUS_PHASE: Lagrangian two-phase flow on frozen i
         continuous phase. This option
         only only be used in case of a calculation restart. All the
         Eulerian fields are frozen (including the scalar fields).
         This option automatically implies \ref iccvfg = 1 */

typedef enum {
  CS_LAGR_OFF = 0,
  CS_LAGR_ONEWAY_COUPLING = 1,
  CS_LAGR_TWOWAY_COUPLING = 2,
  CS_LAGR_FROZEN_CONTINUOUS_PHASE = 3
} cs_lagr_module_status_t;


/*! Fixed maximum sizes */
/*----------------------*/

typedef struct {

  int nusbrd;  /*!< maximum number of additional user
                    particle/boundary interactions */

  int ndlaim;  /*!< maximum number of particle integer data */

  int ncharm2; /*!< maximum number of coal classes */
  int nlayer;  /*!< maximum number of coal layers */

} cs_lagr_const_dim_t;

/*! General dimensions */
/*---------------------*/

typedef struct {

  int   ntersl;  /*!< number of source terms for return coupling */
  int   n_boundary_stats;  /*!< number of boundary statistics */

} cs_lagr_dim_t;

/*! Time and coupling scheme for the Lagrangian module */
/*-----------------------------------------------------*/

typedef struct {

  /*! Lagrangian module status.
     - CS_LAGR_OFF: Lagrangian module off
     - CS_LAGR_ONEWAY_COUPLING: Lagrangian two-phase flow in one-way coupling
         (no influence of the particles on the continuous phase)
     - CS_LAGR_TWOWAY_COUPLING: Lagrangian two-phase flow with two-way coupling
         (influence of the particles on the dynamics of the continuous phase).
         Dynamics, temperature and mass may be coupled independently.
     - CS_LAGR_FROZEN_CONTINUOUS_PHASE: Lagrangian two-phase flow on frozen i
         continuous phase. This option
         only only be used in case of a calculation restart. All the
         Eulerian fields are frozen (including the scalar fields).
         This option automatically implies \ref iccvfg = 1 */
  int  iilagr;

  /*  indicates the steady (=1) or unsteady (=0) state of the
      continuous phase flow
      in particular, \ref isttio = 1 is needed in order to:
      calculate steady statistics in the volume or at the boundaries
      (starting respectively from the iterations \ref nstist)
      and calculate time-averaged two-way coupling source terms (from the
      time step \ref nstits).
      Useful if \ref iilagr = CS_LAGR_ONEWAY_COUPLING
      or \ref iilagr = CS_LAGR_TWOWAY_COUPLING
      (if \ref iilagr = CS_LAGR_FROZEN_CONTINUOUS_PHASE,
      then \ref isttio=1 automatically) */
  int  isttio;

  /*! activation (=1) or not (=0) of a Lagrangian calculation restart.
    The calculation restart file read when this option is activated
    only contains the data related to the particles (see also \ref isuist)
    the global calculation must also be a restart calculation
  */
  int  isuila;

  /*! trajectory algorithm order in time */
  int  t_order;

  /*! activates (>0) or not (=0) the complete turbulent dispersion model.
    When \ref modcpl is strictly positive, its value is interpreted as the
    absolute Lagrangian time step number (including restarts) after which the
    complete model is applied.
    Since the complete model uses volume statistics, \ref modcpl must
    either be 0 or be larger than \ref idstnt. */
  int     modcpl;

  /*!  direction (1=x, 2=y, 3=z) of the complete model.
    it corresponds to the main directions of the flow.
    Useful if \ref modcpl > 0 */
  int     idirla;

  /*!  activation (=1) or not (=0) of the particle turbulent dispersion.
    The turbulent dispersion is compatible only with the RANS turbulent models
    (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$, v2f or \f$k-\omega\f$).
    (\ref iturb=20, 21, 30, 31, 50 or 60). */
  int     idistu;

  /*! \ref idiffl=1 suppresses the crossing trajectory effect, making
    turbulent dispersion for the particles identical to the turbulent
    diffusion of fluid particles.
    Useful if \ref idistu=1 */
  int     idiffl;

  /*! activation (=1) or not (=0) of the solution of a Poisson's equation for
    the correction of the particle instantaneous velocities
    (in order to obtain a null divergence).
    this option is not validated and reserved to the development team.
    Do not change the default value */
  int     ilapoi;

  /*! activation (=1) or not (=0) of the added-mass term.
   \f[ \DP{u_p} = - \dfrac{1}{\rho_p} \grad P + \dfrac{u_s-u_p}{\tau_p}
         + g
         +1/2 C_A \dfrac{\rho_f}{\rho_p} \left( \dfrac{Du}{Dt}-\DP{u_p} \right)
   \f]
   and
    \f[ \rho_f \dfrac{Du}{Dt} \simeq  - \grad P + \rho_f g \f]
   with \f$ C_A = 1\f$. Then
    \f[ \DP{u_p} = - \dfrac{1}{\rho_p} \dfrac{1+C_A/2}
                                        {1+C_A/2\dfrac{\rho_f}{\rho_p}} \grad P
            + \dfrac{u_s-u_p}{\widetilde{\tau}_p}
            + g
    \f]
   with
   \f[ \widetilde{\tau_p} = (1 + C_A /2 \dfrac{\rho_f}{\rho_p}) \tau_p \f] */
  int     iadded_mass;

  /*! Added-mass constant (\f$ C_A = 1\f$) */
  cs_real_t        added_mass_const;

} cs_lagr_time_scheme_t;

/*! Main physical model parameters for the Lagrangian module */
/*-----------------------------------------------------------*/

typedef struct {

  /*! activates (>0) or deactivates (=0) the physical models associated to the
    particles:
    - 1: allows to associate with the particles evolution equations on
         their temperature (in degrees Celsius), their diameter and
         their mass
    - = 2: the particles are pulverised coal particles.
    Evolution equations on temperature (in degree Celsius), mass of
    reactive coal, mass of char and diameter of the shrinking core are
    associated with the particles. This option is available only if the
    continuous phase represents a pulverised coal flame. */
  int  physical_model;  /* FIXME: => enum: CS_LAGR_PHYS_STD,
                                           CS_LAGR_PHYS_COAL,
                                           CS_LAGR_PHYS_HEAT... */
  int  n_temperature_layers;

  int  deposition;
  int  dlvo;

  /*! - 0: no DLVO conditions with roughness surface
      - 1: DLVO conditions with roughness surface */
  int  roughness;

  /*!- 0: no resuspension model
     - 1: resuspension model */
  int  resuspension;

  /* - 0: no clogging model
     - 1: clogging model */
  int  clogging;

  /* - 0: no consolidation model
     - 1: consolidation model */
  int  consolidation;

  int  precipitation;
  int  fouling;

  /*! - 0: no agglomeration model
      - 1: agglomeration model used */
  int agglomeration;

  /*! - 0: no fragmentation model
      - 1: fragmentation model used */
  int fragmentation;

  int  n_stat_classes;

  /*! number of aggregates that form the particle */
  int n_particle_aggregates;

  int  n_user_variables;

} cs_lagr_model_t;

/*! Particle counters for the Lagrangian module */
/*----------------------------------------------*/

typedef struct {

  /*! total number of injected particles, since the beginning,
    including calculation restarts */
  cs_gnum_t   n_g_cumulative_total;

  /*! total number of failed particles, since the beginning,
    including calculation restarts */
  cs_gnum_t   n_g_cumulative_failed;

  /*! total number of particles */
  cs_gnum_t   n_g_total;

  /*! total number of particles*/
  cs_gnum_t   n_g_new;

  /*! number of merged particles*/
  cs_gnum_t   n_g_merged;

  /*! number of exited particles*/
  cs_gnum_t   n_g_exit;

  /*! number of deposited particles */
  cs_gnum_t   n_g_deposited;

  /*! number of fouling particles */
  cs_gnum_t   n_g_fouling;

  /*! number of re-entrained particles*/
  cs_gnum_t   n_g_resuspended;

  /*! total number of failed particles */
  cs_gnum_t   n_g_failed;

  /*! total weight of particles*/
  cs_real_t   w_total;

  /*! weight of new particles*/
  cs_real_t   w_new;

  /*! weight of exited particles*/
  cs_real_t   w_exit;

  /*! weight of deposited particles */
  cs_real_t   w_deposited;

  /*! number of fouling particles */
  cs_real_t   w_fouling;

  /*! weight of resuspended particles */
  cs_real_t   w_resuspended;

} cs_lagr_particle_counter_t;

/*! Specific physical model options for the Lagrangian module */
/* ---------------------------------------------------------- */

typedef struct {

  /*  activation (=1) or not (=0) of an evolution equation on the particle
      temperature (in degrees Celsius).
      Useful if \ref physical_model=1 and if there is a thermal scalar
      associated with the continuous phase.
  */
  int   itpvar;

  /*  activation (=1) or not (=0) of an evolution equation on the particle
      diameter. Useful if \ref physical_model = 1.
  */
  int   idpvar;

  /*  activation (=1) or not (=0) of an evolution equation on the particle mass
      Useful if \ref physical_model = 1
  */
  int   impvar;

  /*  initialization temperature (in degree Celsius) for the particles already
      present in the calculation domain when an evolution equation on
      the particle temperature is activated during a calculation
      (\ref physical_model = 1 and \ref itpvar = 1).
      Useful if \ref isuila = 1 and \ref itpvar = 0 in the previous calculation.
  */
  cs_real_t          tpart;

  /* initialization value for the specific heat (\f$ J.kg^{-1}.K^{-1} \f$)
     of the particles already present
     in the calculation domain when an evolution equation
     on the particle temperature is activated during a calculation
     (\ref physical_model = 1 and \ref itpvar = 1).
     Useful if \ref isuila = 1 and \ref itpvar = 0 in the previous calculation
  */
  cs_real_t          cppart;

} cs_lagr_specific_physics_t;

/*! Parameters of the reentrainment model */
/* -------------------------------------- */

typedef struct {

  /* - 0: no resuspension model
     - 1: resuspension model */
  int   ireent;

  /*  - 0: no head losses calculation for influence of the deposit on the flow
      - 1: head losses calculation for influence of the deposit on the flow */
  int   iflow;

  /* Parameters of the particle resuspension model*/
  cs_real_t          espasg;
  cs_real_t          denasp;
  cs_real_t          modyeq;
  cs_real_t          rayasp;
  cs_real_t          rayasg;

} cs_lagr_reentrained_model_t;

/*! Parameters of the precipitation model */
/* -------------------------------------- */

typedef struct {

  /* number of particle classes*/
  int   nbrclas;
  /* diameter of particles formed by precipitation*/
  cs_real_t          diameter;
  /* density of particles formed by precipitation*/
  cs_real_t          rho;
  /* number of precipitated particles */
  int   *nbprec;
  /*  */
  cs_real_t          *solub;
  /* number of precipitated particles */
  cs_real_t          *mp_diss;

} cs_lagr_precipitation_model_t;

/*! Parameters of the particle clogging model */
/* ------------------------------------------ */

typedef struct {

  cs_real_t          jamlim;
  cs_real_t          mporos;
  cs_real_t          csthpp;
  cs_real_t          diam_mean;

} cs_lagr_clogging_model_t;

/*! Parameters of the particle agglomeration model */
/* ------------------------------------------ */

typedef struct {

  cs_real_t          scalar_kernel;
  cs_real_t          base_diameter;
  cs_real_t          (*function_kernel)(cs_lnum_t, cs_lnum_t);

} cs_lagr_agglomeration_model_t;

/*! Parameters of the particle fragmentation model */
/* ------------------------------------------ */

typedef struct {

  cs_real_t          scalar_kernel;
  cs_real_t          base_diameter;
  cs_real_t          (*function_kernel)(cs_lnum_t);

} cs_lagr_fragmentation_model_t;

/*! Parameters of the particle consolidation model */
/* ----------------------------------------------- */

typedef struct {

  cs_lnum_t          iconsol;
  cs_real_t          rate_consol;
  cs_real_t          slope_consol;
  cs_real_t          force_consol;

} cs_lagr_consolidation_model_t;

/*! Lagrangian time stepping status */
/*----------------------------------*/

typedef struct {

  /* current step id (for 2nd order scheme) */
  int    nor;

  /* duration of a Lagrangian iteration */
  cs_real_t          dtp;

  /* physical time of the Lagrangian simulation */
  cs_real_t          ttclag;

} cs_lagr_time_step_t;

/*! Particle injection parameters for a given zone and particle set */
/*------------------------------------------------------------------*/

typedef struct {

  int         zone_id;               /*!< associated zone id */
  int         set_id;                /*!< associated set id */
  int         location_id;           /*!< associated mesh location id */

  cs_gnum_t   n_inject;              /*!< number of particles injected
                                          at a time for this class and zone */

  int         injection_frequency;   /*!< injection frequency
                                          (if =< 0, only at first iteration) */

  /*! optional injection profile computation function, or NULL */
  cs_lagr_injection_profile_compute_t  *injection_profile_func;

  /*! optional injection profile input data, or NULL */
  void                                 *injection_profile_input;

  /*! velocity condition type:
    - -1 imposed fluid velocity (from cell velocity)
    -  0 imposed velocity along the normal of the boundary face
    -  1 imposed velocity: \ref velocity must be set. */
  int         velocity_profile;

  /*! temperature condition type:
    - 0 fluid temperature
    - 1 imposed temperature */
  int         temperature_profile;

  int         coal_number;          /*!< particle coal number
                                      (if \ref physical_model=2) */

  int         cluster;              /*!< statistical cluster id */

  int         particle_aggregate;

  cs_real_t   velocity_magnitude;   /*!< particle velocity magnitude */
  cs_real_t   velocity[3];          /*!< particle velocity components */

  cs_real_t   temperature;          /*!< particle temperature */

  cs_real_t   diameter;             /*!< particle diameter */
  cs_real_t   diameter_variance;    /*!< particle diameter variance */

  cs_real_t   density;              /*!< particle density */

  cs_real_t   fouling_index;        /*!< fouling index */

  cs_real_t   cp;                   /*!< particle specific heat */

  cs_real_t   stat_weight;          /*!< particle statitistical weight */

  cs_real_t   flow_rate;            /*!< flow rate */

  cs_real_t   emissivity;           /*!< particle emissivity */

} cs_lagr_injection_set_t;

/*! 2-way coupling and source term information. */
/*----------------------------------------------*/

typedef struct {

  /*! activation (=1) or not (=0) of the two-way coupling on the dynamics
    of the continuous phase.
    Useful if \ref iilagr = CS_LAGR_TWOWAY_COUPLING and \ref iccvfg = 0 */
  int  ltsdyn;

  /*! activation (=1) or not (=0) of the two-way coupling on the mass.
    Useful if \ref iilagr = CS_LAGR_TWOWAY_COUPLING, \ref physical_model = 1 and \ref impvar = 1 */
  int  ltsmas;

  /*  if \ref physical_model = 1 and \ref itpvar = 1, \ref ltsthe
   activates (=1) or not (=0) the two-way coupling on temperature.
   if \ref physical_model = 2, \ref ltsthe activates (=1) or not (=0) the
   two-way coupling on the eulerian variables related to pulverised
   coal combustion.
   Useful if \ref iilagr = CS_LAGR_TWOWAY_COUPLING */
  int  ltsthe;

  /*! implicit source term for the continuous phase velocity and
    for the turbulent energy if the \f$k-\varepsilon\f$ model is used */
  int  itsli;

  /*  explicit source term for the turbulent dissipation and the
   turbulent energy if the \f$k-\varepsilon\f$ turbulence model is used
   for the continuous phase */
  int  itske;

  /*! explicit thermal source term for the thermal scalar of
    the continuous phase */
  int  itste;

  /*! implicit thermal source term for the thermal scalar of
    the continuous phase */
  int  itsti;

  /*! mass source term */
  int  itsmas;

  /*  source term for the light volatile matters */
  int  *itsmv1;  //ncharm2

  /*  source term for the heavy volatile matters */
  int  *itsmv2;  //ncharm2

  /*! source term for the carbon released during heterogeneous combustion */
  int  itsco;

  /*! variance of the air scalar */
  int  itsfp4;

  /*! number of absolute time steps (including the restarts)
    after which a time-average of the two-way coupling source terms is
    calculated.
    indeed, if the flow is steady (\ref isttio=1), the average quantities
    that appear in the two-way coupling source terms can be calculated over
    different time steps, in order to get a better precision.
    if the number of absolute time steps is strictly inferior to
    \ref nstits, the code considers that the flow has not yet reached its
    steady state (transition period) and the averages appearing in the source
    terms are reinitialized at each time step, as it is the case for unsteady
    flows (\ref isttio=0).
    Useful if \ref iilagr = CS_LAGR_TWOWAY_COUPLING and \ref isttio = 1 */
  int  nstits;

  /*! number of time steps for source terms accumulations */
  int  npts;

  /*! number of cells, whose vulumetric rate DODO
      (concentration ?)is greather than 0.8 */
  int  ntxerr;

  /*! maximum volumetric concentration reached */
  cs_real_t      vmax;

  /*! maximum massic concentration reached */
  cs_real_t      tmamax;

  /*! source term values */
  cs_real_t     *st_val;

} cs_lagr_source_terms_t;

/*! Boundary or volume condition definitions and data */
/*----------------------------------------------------*/

typedef struct {

  int                         location_id;         /*!< mesh location id */

  int                         n_zones;             /*!< number of zones */
  int                        *zone_type;           /*!< zone type */

  int                        *n_injection_sets;    /*!< number of injection
                                                        sets per zone */
  cs_lagr_injection_set_t   **injection_set;       /*!< injection data per
                                                        set per zone */

  char                       *elt_type;            /*! zone type per
                                                       element, or NULL */

  cs_real_t                  *particle_flow_rate;  /*!< particle flow rate
                                                        per zone per
                                                        statistical class */

} cs_lagr_zone_data_t;

/*! Internal face condition definitions */
/*---------------------------------------*/

typedef struct {

  int  *i_face_zone_id;

} cs_lagr_internal_condition_t;

/*! Encrustation model parameters */
/*--------------------------------*/

typedef struct {

  /*  activates (=1) or not (=0) the option of coal particle fouling.
   It then is necessary to specify the domain boundaries
   on which fouling may take place. Useful if \ref physical_model = 2*/
  int  iencra;

  /*  encrustation data*/
  int  npencr;
  // TODO cf particles->n_part_fou in cs_lagr_tracking.c

  /*  encrustation data*/
//TODO
  cs_real_t  *enc1;//ncharm2
  /*  encrustation data*/
//TODO
  cs_real_t  *enc2;//ncharm2

  /*  limit temperature (in degree Celsius) below which the coal particles do
   not cause any fouling (if the fouling model is activated).
   Useful if \ref physical_model = 2 and \ref iencra = 1*/
//TODO
  cs_real_t  *tprenc;//ncharm2

  /*  ash critical viscosity in \f$ kg.m^{-1}.s^{-1} \f$, in the fouling model
   cf J.D. Watt et T. Fereday (J.Inst.Fuel, Vol.42-p99).
   Useful if \ref physical_model = 2 and \ref iencra = 1*/
//TODO
  cs_real_t  *visref;//ncharm2

  /*  encrustation data */
  cs_real_t  dnpenc;

} cs_lagr_encrustation_t;

/*! Physical and chemical model parameters */
/*-----------------------------------------*/

typedef struct {

  /*! Hamaker constant for the particle/fluid/substrate system */
  cs_real_t  cstham;

  /*! Retardation wavelength for VDW forces
      for the particle/fluid/substrate system */
  cs_real_t  lambda_vdw;

  /*! Dielectric constant of the fluid */
  cs_real_t  epseau;

  /*! Electrokinetic potential of the first solid - particle */
  cs_real_t  phi_p;

  /*! Electrokinetic potential of the second solid - surface */
  cs_real_t  phi_s;

  /*! Valence of ions in the solution (used for EDL forces) */
  cs_real_t  valen;

  /*! Ionic force */
  cs_real_t  fion;

} cs_lagr_physico_chemical_t;

/*! Brownian movement parameters */
/*-------------------------------*/

typedef struct {

  int  lamvbr;  /*!< brownian motion activation */

} cs_lagr_brownian_t;

/*! Boundary interactions statistics parameters */
/*----------------------------------------------*/

typedef struct {

  /*! number of iterations during which steady boundary statistics have
    been accumulated.
    Useful if \ref isttio=1 and \ref nstist inferior
    or equal to the current time step.
    \ref npstf is initialized and updated automatically by the code,
    its value is not to be modified by the user */
  int  npstf;

  /*! number of iterations during which boundary statistics have
    been calculated
    (the potential iterations during which unsteady
    statistics have been calculated are counted in \ref npstft).
    \ref npstft is initialized and updated automatically by the code,
    its value is not to be modified by the user */
  int  npstft;

  /*! activation (=1) or not (=0) of the recording of the number of
    particle/boundary interactions, and of the calculation of the associated
    boundary statistics.
    \ref has_part_impact_nbr is set to 1 when using the particulate average
    \ref imoybr = 2. */
  int  has_part_impact_nbr;

  /*!  activation (=1) or not (=0) of the recording of the particulate mass flow
    related to the particle/boundary interactions, and of the calculation of
    the associated boundary statistics.
    If set to 2 or 3 respectively, time averages and variances are also
    computed.
    \ref has_part_impact_nbr is set to 1 when using \ref iflmbd=1.
   */
  int  iflmbd;

  /*!  activation (=1) or not (=0) of the recording of the angle between a
    particle trajectory and a boundary face involved in a particle/boundary
    interaction, and of the calculation of the associated boundary statistics. */
  int  iangbd;

  /*!  activation (=1) or not (=0) of the recording of the velocity of a particle
    involved in a particle/boundary interaction, and of the calculation of
    the associated boundary statistics. */
  int  ivitbd;

  /*!  activation (=1) or not (=0) of the recording of clogging parameters
    involved in a particle/boundary interaction, and of the calculation of
    the associated boundary statistics. */
  int  iclgst;

  /*! flag for number of recorded particle/boundary interactions with fouling */
  int  iencnbbd;

  /*! flag for mass of fouled coal particles */
  int  iencmabd;

  /*! flag for diameter of fouled coal particles */
  int  iencdibd;

  /*! flag for coke fraction of fouled coal particles */
  int  iencckbd;

  /*!  id for number of particle/boundary interactions */
  int  inbr;

  /*!  id for particle mass flow at the boundary faces */
  int  iflm;

  /*!  id for mean interaction angle with the boundary faces */
  int  iang;

  /*!  id for mean interaction velocity with the boundary faces */
  int  ivit;

  /*!  id for number of resuspended particles */
  int  ires;

  /*!  id for mass flow of resuspended particles at the boundary faces */
  int  iflres;

  /*! flag for number of recorded particle/boundary interactions with fouling */
  int  iencnb;

  /*! id for mass of fouled coal particles */
  int  iencma;

  /*! id for diameter of fouled coal particles */
  int  iencdi;

  /*! id for coke fraction of fouled coal particles */
  int  iencck;

  /*! the recordings in \ref bound_stat at every particle/boundary interaction are
      cumulated values (possibly reset to zero at every iteration in the
      unsteady case). They must therefore be divided by a quantity to
      get boundary statistics. The user can choose between two average types:
      - = 0: no average is applied to the recorded cumulated values.
      - = 1: a time-average is calculated. The cumulated value
             is divided by the physical duration in the case of steady
             averages (\ref isttio=1). The cumulated value is divided by the
             value of the last time step in the case of unsteady averages
             (\ref isttio=0), and also in the case of steady averages while the
             absolute iteration number is inferior to \ref nstist.
      - = 2: a particulate average is calculated. The cumulated value is divided
             by the number of particle/boundary interactions (in terms of
             statistical weight) recorded in \ref bound_stat "bound_stat"(nfabor,inbr).
             This average can only be calculated when \ref has_part_impact_nbr = 1.
             Only the cumulated value is recorded in the restart file. */
  int  *imoybr;

  /*! id for number of deposited particles */
  int  inclg;

  /*! id for particle deposition part */
  int  inclgt;

  /*! id for particle deposition time */
  int  iclogt;

  /*! id for particle consolidation height */
  int  iclogh;

  /*! id for particle surface coverage */
  int  iscovc;

  /* id for mean of particle deposition height */
  int  ihdepm;

  /* id for variance of particle deposition height */
  int  ihdepv;

  /* id for mean diameter of deposited particles */
  int  ihdiam;

  /* id for sum of deposited particle diameters */
  int  ihsum;

  /*! If the recording of the boundary statistics is steady, \ref tstatp
    contains the cumulated physical duration of the recording of the boundary
    statistics.
    If the recording of the boundary statistics is unsteady, then
    \ref tstatp = dtp (it is the Lagrangian time step, because the
    statistics are reset to zero at every time step). */
  cs_real_t tstatp;

  /*!  name of the boundary statistics, displayed in the log
    and the post-processing files.
    Warning: this name is also used to reference information in the restart
    file (\ref isuist =1). If the name of a variable is changed between two
    calculations, it will not be possible to read its value from the restart
    file */
  char           **nombrd;

} cs_lagr_boundary_interactions_t;

/*! Pointers to external (Eulerian solver) data. */
/*-----------------------------------------------*/

typedef struct {

  /* Turbulence model */
  int iturb;
  int itytur;

  /* cpincl */
  int ncharb;

  /* ppppar */
  int ncharm;

  /* radiation */
  int radiative_model;

  /* icp */
  int icp;

  /* diftl0 */
  cs_real_t diftl0;

  /* cmu */
  cs_real_t cmu;

  /* visls0 */
  cs_real_t visls0;

  /* Referenced fields
     ----------------- */

  /* wall ustar */
  cs_field_t *ustar;

  /* Fluid density */
  cs_field_t *cromf;

  /* Fluid pressure */
  cs_field_t *pressure;

  /* Fluid temparature */
  cs_field_t *scal_t;
  cs_field_t *temperature;
  cs_field_t *t_gaz;

  /* Fluid velocity */
  cs_field_t *vel;

  /* Fluid viscosity */
  cs_field_t *viscl;

  /* Fluid viscosity */
  cs_field_t *cpro_viscls;

  /* Fluid specific heat capacity */
  cs_field_t *cpro_cp;

  /* Radiat.       */
  cs_field_t *luminance;

  /* Combustion    */
  cs_field_t *x_oxyd;
  cs_field_t *x_eau;
  cs_field_t *x_m;

  /* Turbulence */
  /* Turbulent intensity */
  cs_field_t *cvar_k;

  /* Turbulent dissipation */
  cs_field_t *cvar_ep;

  /* Omega from k-omega SST model*/
  cs_field_t *cvar_omg;

  /* Reynolds stress component Rxx */
  cs_field_t *cvar_r11;
  /* Reynolds stress component Ryy */
  cs_field_t *cvar_r22;
  /* Reynolds stress component Rzz */
  cs_field_t *cvar_r33;

  /* Reynolds Stress Tensor */
  cs_field_t *cvar_rij;

  /* Total pressure gradient */
  cs_real_3_t *grad_pr;

  /* velocity gradient */
  cs_real_33_t *grad_vel;

} cs_lagr_extra_module_t;

/*! External data relative to coal combustion */
/*--------------------------------------------*/

typedef struct {

  int         ih2o;   // cpincl
  int         io2;    // cpincl
  int         ico;    // cpincl

  int         iatc;   // ppthch
  cs_real_t   prefth; // ppthch
  cs_real_t   trefth; // ppthch

  int         natom;  // = 5;
  cs_real_t  *wmolat; // dim = natom

  int         ngazem; // = 20;
  cs_real_t  *wmole;  // ngazem
  int        *iym1;

  int         ncharm;  // cpincl
  cs_real_t  *a1ch;   // ncharm
  cs_real_t  *h02ch;
  cs_real_t  *e1ch;   //
  cs_real_t  *a2ch;   //
  cs_real_t  *e2ch;   //
  cs_real_t  *y1ch;   //
  cs_real_t  *y2ch;   //
  cs_real_t  *cp2ch;  //
  cs_real_t  *ahetch; //
  cs_real_t  *ehetch; //
  cs_real_t  *rho0ch; //
  cs_real_t  *xwatch; //
  cs_real_t  *xashch; //
  cs_real_t  *thcdch; //

} cs_lagr_coal_comb_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Fixed constants */

extern const cs_lagr_const_dim_t          *cs_glob_lagr_const_dim;

/*! General dimensions */

extern cs_lagr_dim_t                      *cs_glob_lagr_dim;

/*! Time and Lagrangian-Eulerian coupling scheme */
extern cs_lagr_time_scheme_t              *cs_glob_lagr_time_scheme;

/*! Main Lagragian physical model parameters */
extern cs_lagr_model_t                    *cs_glob_lagr_model;

/*! Read-only pointer to global particle counter */
extern const cs_lagr_particle_counter_t   *cs_glob_lagr_particle_counter;

/* Lagrangian log output every frequency_n time steps */

extern int cs_glob_lagr_log_frequency_n;

/* Statisics on boundaries */

extern cs_real_t *bound_stat;

extern cs_lagr_specific_physics_t            *cs_glob_lagr_specific_physics;
extern cs_lagr_reentrained_model_t           *cs_glob_lagr_reentrained_model;
extern cs_lagr_precipitation_model_t         *cs_glob_lagr_precipitation_model;
extern cs_lagr_clogging_model_t              *cs_glob_lagr_clogging_model;

extern cs_lagr_agglomeration_model_t         *cs_glob_lagr_agglomeration_model;
extern cs_lagr_fragmentation_model_t         *cs_glob_lagr_fragmentation_model;

extern cs_lagr_consolidation_model_t         *cs_glob_lagr_consolidation_model;
extern cs_lagr_time_step_t                   *cs_glob_lagr_time_step;
extern cs_lagr_source_terms_t                *cs_glob_lagr_source_terms;
extern cs_lagr_encrustation_t                *cs_glob_lagr_encrustation;
extern cs_lagr_physico_chemical_t            *cs_glob_lagr_physico_chemical;
extern cs_lagr_brownian_t                    *cs_glob_lagr_brownian;
extern cs_lagr_boundary_interactions_t       *cs_glob_lagr_boundary_interactions;

extern cs_lagr_extra_module_t                *cs_glob_lagr_extra_module;
extern cs_lagr_coal_comb_t                   *cs_glob_lagr_coal_comb;

extern const cs_lagr_zone_data_t             *cs_glob_lagr_boundary_conditions;
extern const cs_lagr_zone_data_t             *cs_glob_lagr_volume_conditions;
extern cs_lagr_internal_condition_t          *cs_glob_lagr_internal_conditions;

/* Projection matrices for global to local coordinates on boundary faces */
extern cs_real_33_t  *cs_glob_lagr_b_face_proj;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to injection set structure.
 *
 * This access method ensures the strucure is initialized for the given
 * zone and injection set.
 *
 * \param[in]  zone_data  pointer to boundary or volume conditions structure
 * \param[in]  zone_id    zone id
 * \param[in]  set_id     injection set id
 *
 * \return pointer to injection set data structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_injection_set_t *
cs_lagr_get_injection_set(cs_lagr_zone_data_t  *zone_data,
                          int                   zone_id,
                          int                   set_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize injection set data structure fields to defaults.
 *
 * \param[in, out]   zis  pointer to structure to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_injection_set_default(cs_lagr_injection_set_t  *zis);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get read/write pointer to global particle counter
 *
 * \return  pointer to lagrangian particle counter structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_counter_t *
cs_lagr_get_particle_counter(void);

/*----------------------------------------------------------------------------*/
/*!
  \brief Update global particle counter
 *
 * All fields handled in the local particle set are updated relative
 * to that data (using global sums).
 *
 * \return  pointer to lagrangian particle counter structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_counter_t *
cs_lagr_update_particle_counter(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_particle_counter_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_specific_physics_t *
cs_get_lagr_specific_physics(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_reentrained_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_reentrained_model_t *
cs_get_lagr_reentrained_model(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_precipitation_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_precipitation_model_t *
cs_get_lagr_precipitation_model(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_clogging_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_clogging_model_t *
cs_get_lagr_clogging_model(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_consolidation_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_consolidation_model_t *
cs_get_lagr_consolidation_model(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_time_step_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_time_step_t *
cs_get_lagr_time_step(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_source_terms_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_source_terms_t *
cs_get_lagr_source_terms(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_encrustation_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_encrustation_t *
cs_get_lagr_encrustation(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_physico_chemical_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_physico_chemical_t *
cs_get_lagr_physico_chemical(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_brownian_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_brownian_t *
cs_get_lagr_brownian(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main boundary conditions structure.
 *
 * \return  pointer to current boundary zone data structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_zone_data_t  *
cs_lagr_get_boundary_conditions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main volume conditions structure.
 *
 * \return pointer to current volume zone data structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_zone_data_t  *
cs_lagr_get_volume_conditions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main internal conditions structure.
 *
 * \return pointer to current internal conditions structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_internal_condition_t  *
cs_lagr_get_internal_conditions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the global boundary and volume condition structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_finalize_zone_conditions(void);

/*----------------------------------------------------------------------------
 * Destroy finalize the global cs_lagr_internal_condition_t structure.
 *----------------------------------------------------------------------------*/

void
cs_lagr_finalize_internal_cond(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_boundary_interactions_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_boundary_interactions_t *
cs_get_lagr_boundary_interactions(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_extra_module_t
 *----------------------------------------------------------------------------*/

cs_lagr_extra_module_t *
cs_get_lagr_extra_module(void);

/*----------------------------------------------------------------------------
 * Prepare for execution of the Lagrangian model.
 *
 * This should be called before the fist call to cs_lagr_solve_time_step.
 *
 *  parameters:
 *    dt     <-- time step (per cell)
 *----------------------------------------------------------------------------*/

void
cs_lagr_solve_initialize(const cs_real_t  *dt);

/*--------------------------------------------------------------------
 * Execute one time step of the Lagrangian model.
 *
 * This is the main function for that model.
 *
 *  parameters:
 *    itypfb <-- boundary face types
 *    dt     <-- time step (per cell)
 *-------------------------------------------------------------------- */

void
cs_lagr_solve_time_step(const int         itypfb[],
                        const cs_real_t  *dt);

/*----------------------------------------------------------------------------
 * Return pointers to lagrangian arrays
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   dim_bound_stat   --> dimensions for bound_stat pointer
 *   p_bound_stat     --> bound_stat pointer
 *----------------------------------------------------------------------------*/

void
cs_lagr_init_c_arrays(int          dim_cs_glob_lagr_source_terms[2],
                      cs_real_t  **p_cs_glob_lagr_source_terms);

/*----------------------------------------------------------------------------
 * Free lagrangian arrays
 *
 * This function is intended for use by Fortran wrappers.
 *----------------------------------------------------------------------------*/

void
cs_lagr_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_H__ */
