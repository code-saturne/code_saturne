#ifndef __CS_LAGR_H__
#define __CS_LAGR_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "assert.h"
#include "cs_base.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Lagrangian boundary condition types */

typedef enum {

  CS_LAGR_INLET = 1,
  CS_LAGR_OUTLET = 2,
  CS_LAGR_REBOUND = 3,
  CS_LAGR_DEPO1 = 4,
  CS_LAGR_DEPO2 = 5,
  CS_LAGR_FOULING = 7,
  CS_LAGR_JBORD1 = 8,
  CS_LAGR_JBORD2 = 9,
  CS_LAGR_JBORD3 = 10,
  CS_LAGR_JBORD4 = 11,
  CS_LAGR_JBORD5 = 12,
  CS_LAGR_DEPO_DLVO = 13,
  CS_LAGR_SYM = 14

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

/*! Fxed maximum sizes */
/*---------------------*/

typedef struct {

  int nusbrd;  /*!< maximum number of additional user
                    particle/boundary interactions */

  int nflagm;  /*!< maximum number of boundary zones */
  int ndlaim;  /*!< maximum number of particle integer data */

  int ncharm2; /*!< maximum number of coal classes */
  int nlayer;  /*!< maximal number of coal layers */

} cs_lagr_const_dim_t;

/*! General dimensions */
/*---------------------*/

typedef struct {

  int   ntersl;  /*!< number of source terms for return coupling */
  int   nvisbr;  /*!< number of volume statistics */

} cs_lagr_dim_t;

/*! Time and coupling scheme for the Lagrangian module */
/*-----------------------------------------------------*/

typedef struct {

  /*! Lagrangian module status.
     the different values correspond to the following coupling:
     - 0: Lagrangian module off
     - 1: Lagrangian two-phase flow in one-way coupling (no influence of
          the particles on the continuous phase)
     - 2 Lagrangian two-phase flow with two-way coupling (influence of
         the particles on the dynamics of the continuous phase). Dynamics,
         temperature and mass may be coupled independently.
     - 3 Lagrangian two-phase flow on frozen continuous phase. This option
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
      Useful if \ref iilagr=1 or \ref iilagr=2 (if \ref iilagr=3,
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

  int  n_stat_classes;
  int  n_user_variables;

} cs_lagr_model_t;

/* ========================================================================== */

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

/* ========================================================================== */

typedef struct {

  /*  activation (=1) or not (=0) of an evolution equation on the particle
      temperature (in degrees Celsius).
      Useful if \ref physical_model=1 and if there is a thermal scalar associated with
      the continuous phase
  */
  int   itpvar;

  /*  activation (=1) or not (=0) of an evolution equation on the particle
      diameter. Useful if \ref physical_model = 1
  */
  int   idpvar;

  /*  activation (=1) or not (=0) of an evolution equation on the particle mass
      Useful if \ref physical_model = 1
  */
  int   impvar;

  /*  initialization temperature (in degree Celsius) for the particles already
      present in the calculation domain when an evolution equation on
      the particle temperature is activated during a calculation (\ref physical_model =
      1 and \ref itpvar = 1).
      Useful if \ref isuila = 1 and \ref itpvar = 0 in the previous calculation
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

/* ========================================================================== */

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

/* ========================================================================== */

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

/* ========================================================================== */

typedef struct {

  /* Parameter of the particle clogging model */
  cs_real_t          jamlim;
  cs_real_t          mporos;
  cs_real_t          csthpp;
  cs_real_t          diam_mean;

} cs_lagr_clogging_model_t;

/* ========================================================================== */

typedef struct {

  /* Parameter of the particle consolidation model */
  cs_lnum_t          iconsol;
  cs_real_t          rate_consol;
  cs_real_t          slope_consol;
  cs_real_t          force_consol;

} cs_lagr_consolidation_model_t;

/* ========================================================================== */

typedef struct {

  /* current step id (for 2nd order scheme) */
  int    nor;

  /* duration of a Lagrangian iteration */
  cs_real_t          dtp;

  /* physical time of the Lagrangian simulation */
  cs_real_t          ttclag;

} cs_lagr_time_step_t;

/* ========================================================================== */

typedef struct {

  /*! number of particles per class and per boundary zone */
  cs_lnum_t  nb_part;

  /*! injection frequency
    (if < 0, particles are introduced only at first iteration) */
  int        injection_frequency;

  /*! velocity condition type:
    - -1 imposed fluid velocity (from cell velocity)
    -  0 imposed fluid velocity along the normal of the boundary face, with \ref iuno norm.
    -  1 imposed velocity: \ref iupt \ref ivpt \ref iwpt must be given.
    -  2 velocity profile given by user.*/
  int        velocity_profile;

  /*! distribution profile:
    - 1 uniform distribution,
    - 2 presence rate profile given by user.*/
  int        distribution_profile;

  /*! temperature profile:
    - 1 constant temperature profile
    - 2 temperature profile given by the user */
  int        temperature_profile;

  /*! type of user profiles:
    - 1: flat diameter profile
    - 2: user profile to be defined in cs_user_lagr_boundary_conditions */
  int        diameter_profile;

  /*! type of coal initial composition (if \ref physical_model=2)
    - 1: coal initial composition is given by DP_FCP
    - 0: user profile to be defined cs_user_lagr_boundary_conditions */
  int        coal_profile;

  /*! coal number of the particle (if \ref physical_model=2)*/
  int        coal_number;

  /*! statistics group number */
  int        cluster;

  /*! particle velocity magnitude */
  cs_real_t  velocity_magnitude;

  /*! particle velocity components by class and zone */
  cs_real_t  velocity[3];

  /*! particle temperature (size: n_layer) */
  cs_real_t  *temperature;

  /*! particle diameter */
  cs_real_t   diameter;

  /*! particle diameter variance */
  cs_real_t   diameter_variance;

  /*! density */
  cs_real_t   density;

  /*! fouling index */
  cs_real_t   foul_index;

  /*! particle specific heat */
  cs_real_t   cp;

  /*! particle weight */
  cs_real_t   stat_weight;

  /*! flow rate */
  cs_real_t   flow_rate;

  /*! particle emissivity */
  cs_real_t   emissivity;

  /*! water mass fraction in coal particles */
  cs_real_t   water_mass_fraction;

  /*! active coal mass fraction in coal particles */
  cs_real_t  *coal_mass_fraction;

  /*! coke mass fraction in coal particles */
  cs_real_t  *coke_mass_fraction;

  /*! diameter of shrinking core */
  cs_real_t   shrinking_diameter;

  /*! initial particle diameter (for coal particles) */
  cs_real_t   initial_diameter;

  /*! coke density after pyrolysis (for coal particles) */
  cs_real_t  *coke_density;

} cs_lagr_zone_class_data_t;

/* ========================================================================== */

/* ========================================================================== */

typedef struct {

  /*! activation (=1) or not (=0) of the two-way coupling on the dynamics
    of the continuous phase.
    Useful if \ref iilagr = 2 and \ref iccvfg = 0 */
  int  ltsdyn;

  /*! activation (=1) or not (=0) of the two-way coupling on the mass.
    Useful if \ref iilagr = 2, \ref physical_model = 1 and \ref impvar = 1 */
  int  ltsmas;

  /*  if \ref physical_model = 1 and \ref itpvar = 1, \ref ltsthe
   activates (=1) or not (=0) the two-way coupling on temperature.
   if \ref physical_model = 2, \ref ltsthe activates (=1) or not (=0) the
   two-way coupling on the eulerian variables related to pulverised
   coal combustion.
   Useful if \ref iilagr = 2 */
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
//TODO
  int  *itsmv1;//ncharm2

  /*  source term for the heavy volatile matters */
//TODO
  int  *itsmv2;//ncharm2

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
    Useful if \ref iilagr = 2 and \ref isttio = 1 */
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

/* ========================================================================== */

/* Structures useful to deal with boundary conditions
   For USLABO => _boundary_track_treatment */

typedef struct {

  int         n_b_zones;       /* NFRLAG */
  int         n_b_max_zones;

  cs_lnum_t  *b_zone_id;       /* ILFLAG */
  int        *b_zone_classes;  /* IUSNCL */
  int        *b_zone_natures;  /* IUSCLB */

  int        *b_face_zone_id;  /* IFRLAG */

  bool        steady_bndy_conditions;

  cs_real_t  *particle_flow_rate; /* DEBLAG -> post-processing use */

} cs_lagr_bdy_condition_t;

/* Structures useful to deal with iternal conditions */

typedef struct {

  int  *i_face_zone_id;

} cs_lagr_internal_condition_t;

/* ========================================================================== */

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

/* ========================================================================== */

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

/* ========================================================================== */

typedef struct {

  /* brownnian motion activation */
  int  lamvbr;

} cs_lagr_brownian_t;

/* ========================================================================== */

typedef struct {

  /*! number of additional user data to record for the calculation
    of additional boundary statistics in \ref bound_stat */
  int  nusbor;

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
    \ref inbrbd = 1 is a compulsory condition to use the particulate average
    \ref imoybr = 2.  */
  int  inbrbd;

  /*!  activation (=1) or not (=0) of the recording of the particulate mass flow
    related to the particle/boundary interactions, and of the calculation of
    the associated boundary statistics.
    \ref inbrbd = 1 is a compulsory condition to use \ref iflmbd=1.
    Useful if \ref inbrbd=1 */
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

  /*!  supplementary user boundary statistics */
  int  *iusb;

  /*  the recordings in \ref bound_stat at every particle/boundary interaction are
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
             This average can only be calculated when \ref inbrbd=1.
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

  /*!  if the recording of the boundary statistics is steady, \ref tstatp
    contains the cumulated physical duration of the recording of the boundary
    statistics.
    if the recording of the boundary statisticss is unsteady, then
    \ref tstat=dtp (it is the Lagrangian time step, because the
    statistics are reset to zero at every time step). */
  cs_real_t        tstatp;

  /*!  name of the boundary statistics, displayed in the log
    and the post-processing files.
    Warning: this name is also used to reference information in the restart
    file (\ref isuist =1). If the name of a variable is changed between two
    calculations, it will not be possible to read its value from the restart
    file */
  char           **nombrd;

} cs_lagr_boundary_interactions_t;

/* ========================================================================== */

typedef struct {

  /* Turbulence model */
  int iturb;
  int itytur;

  /* cpincl */
  int ncharb;

  /* ppppar */
  int ncharm;

  /* radiation */
  int iirayo;

  /* icp */
  int icp;

  /* diftl0 */
  cs_real_t diftl0;

  /* cmu */
  cs_real_t cmu;

  /* visls0 */
  cs_real_t visls0;

  /*****************
   * Useful fields *
   *****************/

  /* wall ustar */
  cs_real_t *uetbor;

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

} cs_lagr_extra_module_t;

/* external data relative to coal combustion */

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

extern const cs_lagr_const_dim_t   *cs_glob_lagr_const_dim;

/*! General dimensions */

extern cs_lagr_dim_t           *cs_glob_lagr_dim;

/*! Time and Lagrangian-Eulerian coupling scheme */
extern cs_lagr_time_scheme_t   *cs_glob_lagr_time_scheme;

/*! Main Lagragian physical model parameters */
extern cs_lagr_model_t         *cs_glob_lagr_model;

/*! Read-only pointer to global particle counter */
extern const cs_lagr_particle_counter_t      *cs_glob_lagr_particle_counter;

/* Lagrangian log output every frequency_n time steps */

extern int cs_glob_lagr_log_frequency_n;

/* Statisics on borders*/
extern cs_real_t *bound_stat;

extern int cs_glob_lagr_nzone_max;
extern int cs_glob_lagr_nclass_max;

extern cs_lagr_specific_physics_t            *cs_glob_lagr_specific_physics;
extern cs_lagr_reentrained_model_t           *cs_glob_lagr_reentrained_model;
extern cs_lagr_precipitation_model_t         *cs_glob_lagr_precipitation_model;
extern cs_lagr_clogging_model_t              *cs_glob_lagr_clogging_model;
extern cs_lagr_consolidation_model_t         *cs_glob_lagr_consolidation_model;
extern cs_lagr_time_step_t                   *cs_glob_lagr_time_step;
extern cs_lagr_source_terms_t                *cs_glob_lagr_source_terms;
extern cs_lagr_encrustation_t                *cs_glob_lagr_encrustation;
extern cs_lagr_physico_chemical_t            *cs_glob_lagr_physico_chemical;
extern cs_lagr_brownian_t                    *cs_glob_lagr_brownian;
extern cs_lagr_boundary_interactions_t       *cs_glob_lagr_boundary_interactions;

extern cs_lagr_extra_module_t                *cs_glob_lagr_extra_module;
extern cs_lagr_coal_comb_t                   *cs_glob_lagr_coal_comb;

extern cs_lagr_bdy_condition_t               *cs_glob_lagr_bdy_conditions;
extern cs_lagr_internal_condition_t          *cs_glob_lagr_internal_conditions;
extern cs_lagr_zone_class_data_t             *lagr_zone_class_data;

/* Unit normals and offsets of boundary faces */
extern cs_real_4_t   *cs_glob_lagr_b_u_normal;

/* Projection matrices for global to local coordinates on boundary faces */
extern cs_real_33_t  *cs_glob_lagr_b_face_proj;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to class/boundary zone parameters structure
 *
 * \param[in]    iclass     particle class number
 * \param[in]    izone      boundary zone number
 *
 * \return
 *   pointer to particle class and boundary zone structure of parameters
 */
/* ----------------------------------------------------------------------------*/

cs_lagr_zone_class_data_t *
cs_lagr_get_zone_class_data(int  iclass,
                            int  izone);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set injection parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   number     pointer to number of particles to inject
 * \param[in]   freq       pointer to injection frequency
 * \param[in]   stat       pointer to statistical groups id
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_injection(int  iclass,
                                 int  izone,
                                 int  number,
                                 int  freq,
                                 int  stat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set temperature parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    temperature profile
 * \param[in]   temp       pointer to temperature values
 * \param[in]   emissivity emissivity value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_temperature(int         iclass,
                                   int         izone,
                                   int         profile,
                                   cs_real_t  *temp,
                                   cs_real_t   emissivity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set temperature parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   cp         pointer to specific heat value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_cp(int        iclass,
                          int        izone,
                          cs_real_t  cp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set coal parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    coal profile
 * \param[in]   number     coal number
 * \param[in]   temp       pointer to temperature array
 * \param[in]   coal_mf    pointer to coal mass fraction
 * \param[in]   coke_mf    pointer to coke mass fraction
 * \param[in]   coke_density  pointer to coke density after pyrolysis
 * \param[in]   water_mf   pointer to water mass fraction
 * \param[in]   shrink_diam  pointer to coke shrinking diameter
 * \param[in]   init_diam  pointer to initial particle diameter
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_coal(int         iclass,
                            int         izone,
                            int         profile,
                            int         number,
                            cs_real_t  *temp,
                            cs_real_t  *coal_mf,
                            cs_real_t  *coke_mf,
                            cs_real_t  *coke_density,
                            cs_real_t   water_mf,
                            cs_real_t   shrink_diam,
                            cs_real_t   init_diam);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set coal parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    pointer to flag for flow and stat weight profile
 * \param[in]   weight     pointer to stat weight value
 * \param[in]   flow       pointer to mass flow rate value
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_stat(int        iclass,
                            int        izone,
                            int        profile,
                            cs_real_t  weight,
                            cs_real_t  flow);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set diameter parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    pointer to flag for diameter profile
 * \param[in]   diam       pointer to diameter value
 * \param[in]   diam_dev   pointer to diameter standard deviation value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_diam(int        iclass,
                            int        izone,
                            int        profile,
                            cs_real_t  diam,
                            cs_real_t  diam_dev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set density for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   density    pointer to density value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_density(int        iclass,
                               int        izone,
                               cs_real_t  density);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set density for a given class and boundary zone
 *
 * \param[in]   iclass      class number
 * \param[in]   izone       boundary zone number
 * \param[in]   foul_index  pointer to fouling index value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_foul_index(int        iclass,
                                  int        izone,
                                  cs_real_t  foul_index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set velocity parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    pointer to velocity profile
 * \param[in]   velocity   pointer to velocity values array
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_velocity(int        iclass,
                                int        izone,
                                int        profile,
                                cs_real_t  velocity[]);

/*----------------------------------------------------------------------------
 * \brief Initialize Lagrangian module parameters for a given set of data
 *
 *
 *----------------------------------------------------------------------------*/

cs_lagr_zone_class_data_t *
cs_lagr_init_zone_class_new(int       iclass,
                            int       izone);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get read/write pointer to global particle counter
 *
 * \return
 *   pointer to lagrangian particle counter structure
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
 * \brief Return pointer to the main internal conditions structure.
 *
 * \return
 *   pointer to current internal_contditions or NULL
 */
/*----------------------------------------------------------------------------*/

cs_lagr_internal_condition_t  *
cs_lagr_get_internal_conditions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main boundary conditions structure.
 *
 * \return
 *   pointer to current bdy_conditions or NULL
 */
/*----------------------------------------------------------------------------*/

cs_lagr_bdy_condition_t *
cs_lagr_get_bdy_conditions(void);

/*----------------------------------------------------------------------------
 * Destroy finalize the global cs_lagr_bdy_condition_t structure.
 *----------------------------------------------------------------------------*/

void
cs_lagr_finalize_bdy_cond(void);

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
 *
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
