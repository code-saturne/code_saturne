/*============================================================================
 * Lagrangian model options.
 *============================================================================*/

/* VERS */

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

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_lagr.h"
#include "cs_lagr_post.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_prototypes.h"
#include "cs_notebook.h"
#include "cs_prototypes.h"

/*---------------------------------------------------------------------------*/
/*
 * \brief User function of the Lagrangian particle-tracking module
 *
 *  User input of physical, numerical and post-processing options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_model(void)
{
  /* Particle-tracking mode
   * ====================== */

  /* iilagr = CS_LAGR_OFF: no particle tracking (default)
   *        = CS_LAGR_ONEWAY_COUPLING: particle-tracking one-way coupling
   *        = CS_LAGR_TWOWAY_COUPLING: particle-tracking two-way coupling
   *        = CS_LAGR_FROZEN_CONTINUOUS_PHASE: particle tracking on frozen field
   *     (this option requires a calculation restart isuite=1,
   *     all Eulerian fields are frozen (pressure, velocities,
   *     scalars). This option is stronger than iccvfg)     */

  cs_glob_lagr_time_scheme->iilagr = CS_LAGR_ONEWAY_COUPLING;

  /* Particle-tracking calculation restart
   * ===================================== */

  /* isuila = 0 : no restart (default)
     = 1 : restart (this value requires a restart on the continuous
     phase too, i.e. isuite = 1)    */

  cs_glob_lagr_time_scheme->isuila = 0;

  /* Restart on volume and boundary statistics, and two-way coupling terms; */
  /* useful if isuila = 1 (defaul off: 0 ; on: 1)  */

  if (cs_glob_lagr_time_scheme->isuila == 1)
    cs_glob_lagr_stat_options->isuist = 0;

  /* Particle tracking: specific models
   * ================================== */

  /* iphyla = 0 : only transport modeling (default)
   *         = 1 : equation on temperature (in Celsius degrees), diameter or mass
   *         = 2 : pulverized coal combustion (only available if the continuous
   *               phase is a flame of pulverized coal)     */

  cs_glob_lagr_model->physical_model = 0;

  /* 3.1 equation on temperature, diameter or mass */
  if (cs_glob_lagr_model->physical_model == 1) {
    /* equation on diameter */
    /* (default off: 0 ; on: 1)  */
    cs_glob_lagr_specific_physics->idpvar   = 0;

    /* equation on temperature (in Celsius degrees)  */
    /* (default off: 0 ; on: 1)  */
    /* This option requires a thermal scalar for the continuous phase.   */
    cs_glob_lagr_specific_physics->itpvar   = 0;

    /* equation on mass     */
    /* (default off: 0 ; on: 1)  */
    cs_glob_lagr_specific_physics->impvar   = 0;
  }

  /* Coal fouling
   * ---------------------------------------------------------------------
   * Reference internal reports EDF/R&D: HI-81/00/030/A and HI-81/01/033/A
   *
   *  Evaluation of the probability for a particle to stick to a wall.
   *  This probability is the ratio of a critical viscosity on the
   *  viscosity of coal ashes
   *
   *           visref
   *  P(Tp) = --------      for viscen >= visref
   *           viscen
   *
   *        = 1             otherwise
   *
   *
   *  The expression of J.D. Watt and T.Fereday (J.Inst.Fuel-Vol42-p99)
   *  is used to evaluate the viscosity of the ashes
   *
   *                     Enc1 * 1.0d+7
   *  Log (10*viscen) = --------------- + Enc2
   *    10                           2
   *                    (Tp(C) - 150)
   *
   *  In literature, the range of the critical viscosity visref is between
   *  8 Pa.s and 1.D7 Pa.s For general purpose 1.0D+4 Pa.s is chosen
   *----------------------------------------------------------------------- */

  if (cs_glob_lagr_model->physical_model == 2) {
    /* iencra = 0 no fouling (default)
       = 1 fouling
       The boundary on which the fouling can occur must be specified with
       boundary condition definitions.

       * Post-processing:
       * iencnbbd = 1 / iencmabd = 1 / iencdibd = 1 /iencckbd = 1 (10.2) */

    cs_glob_lagr_model->fouling = 0;

    /* Example of definition of fouling criteria for each coal first
       (and single) coal icha = 1    */
    int icha = 0;

    /* tprenc : threshold temperature below which no fouling occurs
       (in degrees Celcius) */
    cs_glob_lagr_encrustation->tprenc[icha] = 600.0;

    /* visref : critical viscosity (Pa.s) */
    cs_glob_lagr_encrustation->visref[icha] = 10000.0;

    /* > coal composition in mineral matters:
       (with SiO2 + Al2O3 + Fe2O3 + CaO + MgO = 100% in mass)  */
    cs_real_t sio2  = 36.0;
    cs_real_t al2o3 = 20.8;
    cs_real_t fe2o3 = 4.9;
    cs_real_t cao   = 13.3;

    /* Enc1 and Enc2 : coefficients in Watt and Fereday expression  */
    cs_glob_lagr_encrustation->enc1[icha]
      = 0.00835 * sio2 + 0.00601 * al2o3 - 0.109;
    cs_glob_lagr_encrustation->enc2[icha]
      =  0.0415  * sio2 + 0.0192  * al2o3 + 0.0276 * fe2o3 + 0.016 * cao - 3.92;
  }

  /* Calculation features for the dispersed phases
   * ============================================= */

  /* Additional variables
   * --------------------
   *
   *   Additional variables may be accessed using the (CS_LAGR_USER + i)
   *   attribute, where 0 <= i < lagr_params->n_user_variables
   *   is the additional variable index.
   *
   *   The integration of the associated differential stochastic equation
   *   requires a user intervention in cs_user_lagr_sde() function */

  cs_lagr_set_n_user_variables(0);

  /* Steady or unsteady continuous phase
   * -----------------------------------
   *   if steady: isttio = 1
   *   if unsteady: isttio = 0
   *   if iilagr = CS_LAGR_FROZEN_CONTINUOUS_PHASE then isttio = 1

   Remark: if isttio = 0, then the statistical averages are reset
   at each time step   */

  if (cs_glob_lagr_time_scheme->iilagr != CS_LAGR_FROZEN_CONTINUOUS_PHASE)
    cs_glob_lagr_time_scheme->isttio   = 0;

  /* Two-way coupling: (iilagr = CS_LAGR_TWOWAY_COUPLING)
     ------------------------------ */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {
    /* * number of absolute time step (i.e. with restart)
       from which a time average for two-way coupling source terms is
       computed (steady source terms)
       * if the time step is lower than NSTITS, source terms are
       unsteady: they are reset at each time step
       * useful only if ISTTIO = 1.
       * the min value for NSTITS is 1 */

    cs_glob_lagr_source_terms->nstits   = 1;

    /* two-way coupling for dynamic (velocities and turbulent scalars)   */
    /* (default off: 0 ; on: 1)  */
    /* (useful if ICCVFG = 0)    */

    cs_glob_lagr_source_terms->ltsdyn   = 0;

    /* two-way coupling for mass (if IPHYLA = 1 and IMPVAR = 1)     */
    /* (default off: 0 ; on: 1)  */

    if (   cs_glob_lagr_model->physical_model == 1
        && (   cs_glob_lagr_specific_physics->impvar == 1
            || cs_glob_lagr_specific_physics->idpvar == 1))
      cs_glob_lagr_source_terms->ltsmas     = 0;

    /* two-way coupling for thermal scalar */
    /* (if iphyla = 1 and impvar = 1, or iphyla = 2) */
    /* or for coal variables (if IPHYLA = 2)    */
    /* (default off: 0 ; on: 1)  */

    if (   (   cs_glob_lagr_model->physical_model == 1
            && cs_glob_lagr_specific_physics->itpvar == 1)
        || cs_glob_lagr_model->physical_model == 2)
      cs_glob_lagr_source_terms->ltsthe     = 0;

  }

  /* Volume statistics
     ----------------- */

  /* Threshold for the management of volume statistics
     -------------------------------------------------
     * the value of the seuil variable is a statistical weight.
     * each cell of the mesh contains a statistical weight
     (sum of the statistical weights of all the particles
     located in the cell); seuil is the minimal value from
     which the contribution in statistical weight of a particle
     is not taken into account anymore in the full model
     of turbulent dispersion, in the resolution of the
     Poisson equation of correction of the mean velocities, and
     in the writing of the log and post-processing. */

  cs_glob_lagr_stat_options->threshold = 0.0;

  /* Calculation of the volume statistics from the absolute number
   * of time steps
   * * idstnt is a absolute number of time steps
   * (i.e. including calculation restarts) */

  cs_glob_lagr_stat_options->idstnt = 1;

  /* Steady calculation from the absolute time step nstist
   *   * nstist is a absolute number of time steps
   *     (i.e. including calculation restarts) from which the statistics
   *     are averaged in time.
   *   * useful if the calculation is steady (isttio=1)
   *   * if the number of time steps is lower than nstits,
   *     the transmitted source terms are unsteady (i.e. they are reset to
   *     zero ar each time step)
   *   * the minimal value acceptable for nstist is 1.    */

  cs_glob_lagr_stat_options->nstist = cs_glob_lagr_stat_options->idstnt;

  /* Volume statistical variables
     ---------------------------- */

  /* Activation of the calculation of the particle volume fraction */

  cs_lagr_stat_activate(CS_LAGR_STAT_VOLUME_FRACTION);

  /* Activation of the calculation of the particle velocity */

  cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY);

  /* Activation of the calculation of the particle residence time */

  cs_lagr_stat_activate_attr(CS_LAGR_RESIDENCE_TIME);

  /* Activation of the calculation of the weight */

  cs_lagr_stat_activate_attr(CS_LAGR_STAT_WEIGHT);

  /*Specific models (iphyla = 1) following the chosen options: */
  /* Mean and variance of the temperature     */
  /* Mean and variance of the diameter   */
  /* Mean and variance of the mass  */

  /* Statistics per class
   * -------------------- */

  cs_glob_lagr_model->n_stat_classes = 0;

  /* Options concerning the numerical treatment of the dispersed phase
   * ================================================================= */

  /* Integration order of the stochastic differential equations   */
  /* (default 2; acceptable values 1 or 2)    */

  cs_glob_lagr_time_scheme->t_order = 1;

  /* Options concerning the treatment of the dispersed phase
   * ======================================================= */

  /*  Caution: In this version, the turbulent dispersion works only if
               the continuous phase is calculated with a k-eps or
               a Rij-eps model

      Activation of the turbulent dispersion
      (default on: 1 ; off: 0) */

  cs_glob_lagr_time_scheme->idistu = 1;

  /* Turbulent dispersion imposed to the fluid one.

     If activated, then particle turbulent dispersion is
     equal to the fluid-particle one. The crossing-trajectory effects
     are suppressed ; it is then a case of turbulent diffusion. If the
     simulated particle density is equal to the fluid density, then
     we are simulating the displacement of fluid particles.
     (default off: 0 ; on: 1) */

  cs_glob_lagr_time_scheme->idiffl = 0;

  /* modcpl :
     = 0 for the incomplete model (default value)
     > 0 for the full model, is equal the absolute number
     of time steps from which the full model is activated
     modcpl must not be lower than idstnt */

  cs_glob_lagr_time_scheme->modcpl = 0;

  /* idirla (=1 or 2 or 3) : 1st, 2nd or 3rd direction
     of the full model. Corresponds to the main direction
     of the flow. Allow to calculate a non-isotropic Lagrangian timescale
     (default idirla=1) */

  if (cs_glob_lagr_time_scheme->modcpl > 0)
    cs_glob_lagr_time_scheme->idirla = 1;

  /* Options concerning the treatment of specific forces
   * =================================================== */

  /* idlvo = 0
     = 1 dlvo deposition conditions are activated for the
     wall with appropriate conditions idepfa (see cs_user_lag2.c) */

  cs_glob_lagr_model->dlvo  = 0;

  if (cs_glob_lagr_model->dlvo == 1) {
    /* Constants for the van der Waals forces
       --------------------------------------
       Hamaker constant for the particle/fluid/substrate system:*/
    cs_glob_lagr_physico_chemical->cstham = 6e-20;

    /* Retardation wavelength for the particle/fluid/substrate system:*/
    cs_glob_lagr_physico_chemical->lambda_vdw = 1000.0;

    /* Constants for the elecstrostatic forces
       ---------------------------------------
       Dielectric constant of the fluid (example: water at 293 K)*/
    cs_glob_lagr_physico_chemical->epseau = 80.1;

    /* Electrokinetic potential of the first solid - particle (Volt)*/
    cs_glob_lagr_physico_chemical->phi_p  = 0.05;

    /* Electrokinetic potential of the second solid - surface (Volt)*/
    cs_glob_lagr_physico_chemical->phi_s  =  -0.05;

    /* Valency of ions in the solution (used for EDL forces)*/
    cs_glob_lagr_physico_chemical->valen  = 1.0;

    /* Ionic force (mol/l)*/
    cs_glob_lagr_physico_chemical->fion   = 0.01;
  }

  /* Activation of Brownian motion
   * ============================= */

  /* Activation of Brownian motion:
     (default off: 0 ; on: 1)
     Caution: OPTION FOR DEVELOPERS ONLY
     ======== */
  cs_glob_lagr_brownian->lamvbr = 0;

  /* Activation of deposition model
   * ============================== */

  /* Activation of the deposition model (default off: 0 ; on: 1) */
  cs_glob_lagr_model->deposition = 0;

  /* Activation of roughness and resuspension model
   * ============================================== */

  /* Activation of the resuspension model (default off: 0 ; on: 1) */
  cs_glob_lagr_model->resuspension = 0;

  /* Caution: OPTION FOR DEVELOPERS ONLY
     ========
     dlvo deposition conditions for roughness surface */

  cs_glob_lagr_model->roughness = 0;

  /* Parameters of the particle resuspension model for the roughness */

  /* average distance between two large-scale asperities */
  cs_glob_lagr_reentrained_model->espasg = 2e-05;

  /* density of the small-scale asperities */
  cs_glob_lagr_reentrained_model->denasp = 63600000000000.0;

  /* radius of small asperities */
  cs_glob_lagr_reentrained_model->rayasp = 5e-09;

  /* radius of large asperities */
  cs_glob_lagr_reentrained_model->rayasg = 2e-06;

  /* Young's modulus (GPa) */
  cs_glob_lagr_reentrained_model->modyeq = 266000000000.0;

  /* Activation of the clogging model
   * ================================ */

  /* Activation of the clogging model
     (default off: 0 ; on: 1)
     Caution: OPTION FOR DEVELOPERS ONLY
     ======== */

  cs_glob_lagr_model->clogging = 0;

  /* Parameters for the particle clogging model */

  /* Mean diameter*/
  cs_glob_lagr_clogging_model->diam_mean = 1.0e-6;

  /* Jamming limit */
  cs_glob_lagr_clogging_model->jamlim      = 0.74;

  /* Minimal porosity
   * from 0.366 to 0.409 for random packings
   * equal to 0.26 for close packings */
  cs_glob_lagr_clogging_model->mporos      = 0.366;

  /* Hamaker constant for the particle/fluid/particle system */
  cs_glob_lagr_clogging_model->csthpp      = 5e-20;

  /* Influence of the deposit on the flow
   * ==================================== */

  /* Activation of the influence of the deposit on the flow
     by the head losses calculation (with clogging model only)
     (default off: 0 ; on: 1) */

  cs_glob_lagr_reentrained_model->iflow  = 0;

  if (cs_glob_lagr_reentrained_model->iflow == 1) {

    /* One-way coupling */
    cs_glob_lagr_time_scheme->iilagr  = CS_LAGR_ONEWAY_COUPLING;

    /* The statistical averages are not reset
       at each time step */
    cs_glob_lagr_time_scheme->isttio   = 1;

  }

  /* Activation of the consolidation model
   * ===================================== */

  /* Activation of the consolidation model
     (default off: 0 ; on: 1) */

  /* Caution: valid only for multilayer deposition: */
  if (cs_glob_lagr_model->clogging > 0)
    cs_glob_lagr_model->consolidation = 0;

  /* Parameters for the particle consolidation model */

  /* Consolidated height hconsol calculated using the deposit time
   * hconsol = t_depo * rconsol
   * Adhesion calculated using the following formula:
   * Fadh = F_consol + (F_DLVO - F_consol)
   *        * (0.5+0.5*tanh((h-hconsol)/kconsol/hconsol))
   */

  /* Consolidated force (N) */
  cs_glob_lagr_consolidation_model->force_consol = 3.0e-8;

  /* Slope of consolidation (->0 for a two-layer system) */
  cs_glob_lagr_consolidation_model->slope_consol = 0.1;

  /* Consolidation rate (m/s) */
  cs_glob_lagr_consolidation_model->rate_consol  = 4.0e-3;

  /* Activation of the precipitation/disolution model
   * ================================================ */

  /* Activation of the precipitation/dissolution model
     (default off: 0 ; on: 1)
     Caution: OPTION FOR DEVELOPERS ONLY */
  cs_glob_lagr_model->precipitation = 0;

  /* Diameter of particles formed by precipitation */
  cs_glob_lagr_precipitation_model->diameter  = 2e-06;

  /* Diameter of particles formed by precipitation */
  cs_glob_lagr_precipitation_model->rho       = 5200.0;

  /* Number of particle classes */
  cs_glob_lagr_precipitation_model->nbrclas   = 2;

  /* Boundary statistics
   * =================== */

  /* The boundary statistic 'number of particle/boundary interactions' must be
     selected to activate the particle average imoybr(...) = 2 */

  /* Number of particle/boundary interactions
     (default off: 0 ; on: 1) */
  cs_glob_lagr_boundary_interactions->has_part_impact_nbr      = 1;

  /* Particle mass flux associated to particle/boundary interactions
     (default off: 0 ; on: 1)*/
  cs_glob_lagr_boundary_interactions->iflmbd      = 1;

  /* Angle between particle velocity and the plan of the boundary face
     (default off: 0 ; on: 1) */
  cs_glob_lagr_boundary_interactions->iangbd      = 0;

  /* Norm of particle velocity during the interation with the boundary face
     (default off: 0 ; on: 1) */
  cs_glob_lagr_boundary_interactions->ivitbd      = 0;

  /* (default off: 0 ; on: 1) */
  if (   cs_glob_lagr_model->physical_model == 2
      && cs_glob_lagr_model->fouling == 1) {

    /* Number of particle/boundary interactions with fouling */
    cs_glob_lagr_boundary_interactions->iencnbbd      = 0;

    /* Mass of fouled coal particles */
    cs_glob_lagr_boundary_interactions->iencmabd      = 0;

    /* Diameter of fouled coal particles */
    cs_glob_lagr_boundary_interactions->iencdibd      = 0;

    /* Coke fraction of fouled coal particles */
    cs_glob_lagr_boundary_interactions->iencckbd      = 0;
  }

  /* Name of the recordings for display,
     Average in time of particle average
     of the boundary statistics
     -----------------------------------

     * A priori the user intervenes only in the additional user information
     to be recorded: he must prescribe the name of the recording as well as
     the type of average that he wishes to apply to it for the writing
     of the log and the post-processing.

     * The applied average is prescribed through the imoybr array:
     - if imoybr(iusb(ii)) = 0 -> no average applied
     - if imoybr(iusb(ii)) = 1 -> a time average is applied, i.e. the
         statistic is divided by the last time step in the case of an unsteady
         calculation with a number of iterations lower than nstist; or that
         the statistic is divided by the recording time in the case of a
         steady calculation.
     - if imoybr(iusb(ii)) = 2 -> a particle average is applied, i.e. the
         statistic is divided by the number of recorded particle/boundary
         interactions (in terms of statistical weight) in bound_stat(nfabor,inbr)
         To use this average, has_part_impact_nbr is set to 1.
     - if imoybr(iusb(ii)) = 3 -> (coal fouling only) a particle average
         is applied, i.e. the statistic is divided by the number of recorded
         particle/boundary interactions with fouling (in terms of statistical
         weight) in bound_stat(nfabor,inbr), To use this average, iencnbbd must be
         set to 1.
     * The back-ups in the restart file are performed without applying
       this average. */

  /* Frequency for the output of the Lagrangian log
   * ============================================== */

  cs_glob_lagr_log_frequency_n = 1;

  /* Post-process particle attributes
   * ================================ */

  cs_lagr_post_set_attr(CS_LAGR_STAT_CLASS, true);


}

/*----------------------------------------------------------------------------*/

END_C_DECLS
