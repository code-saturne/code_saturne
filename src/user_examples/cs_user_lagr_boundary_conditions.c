/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_math.h"
#include "cs_selector.h"
#include "cs_parameters.h"

#include "cs_mesh.h"

#include "cs_lagr.h"
#include "cs_lagr_new.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_prototypes.h"

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle boundary conditions.
 *
 * This is used definition of for inlet and of the other boundaries
 *
 * \param[in] bc_type    type of the boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_boundary_conditions(const int  bc_type[])
{
  cs_mesh_t *mesh = cs_glob_mesh;

  cs_lnum_t  nlelt = mesh->n_b_faces;

  cs_lagr_bdy_condition_t *lagr_bdy_cond = cs_lagr_get_bdy_conditions();
  cs_lagr_coal_comb_t *lag_cc = cs_glob_lagr_coal_comb;

  /* Allocate a temporary array for boundary faces selection */

  cs_lnum_t *lstelt;
  BFT_MALLOC(lstelt, nlelt, cs_lnum_t);

  /* Construction of the boundary zones */

  /* Definition of the boundary zones
     For the Lagrangian module, the user defines nfrlag boundary zones
     from the properties (groups, ...) of the boundary faces, from the
     boundary conditions, or even from their coordinates. To do this,
     we fill the ifrlag(mesh->n_b_faces) array which gives for every boundary
     face the number (id+1) of the zone to which it belongs ifrlag(face_id)

     Be careful, all the boundary faces must have been assigned.

     The number of the zones (thus the values of ifrlag(face_id)) is arbitrarily
     chosen by the user, but must be a positive integer and inferior or equal
     to nflagm.

     Afterwards, we assign to every zone a type named itylag that will be used
     to prescribe global boundary conditions. */

  /* First zone numbered izone = 0 */
  cs_selector_get_b_face_list("10",
                              &nlelt,
                              lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {
    cs_lnum_t face_id = lstelt[ilelt];
    lagr_bdy_cond->b_face_zone_id[face_id] = 0; /* zone_id */
  }

  /* Second zone numbered izone = 1 */
  cs_selector_get_b_face_list("4 and y < 1.0",
                              &nlelt,
                              lstelt);

  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {
    cs_lnum_t face_id = lstelt[ilelt];
    lagr_bdy_cond->b_face_zone_id[face_id] = 1; /* zone_id */
  }

  /* Third zone numbered izone = 3 (inlet) */

  for (cs_lnum_t face_id = 0; face_id < nlelt; face_id++) {
    if (bc_type[face_id] == CS_INLET) {
      lagr_bdy_cond->b_face_zone_id[face_id] = 3; /* zone_id */
    }
  }

  /* Nth zone numbered izone = 4 */
  cs_selector_get_b_face_list("4 and y < 1.0",
                              &nlelt,
                              lstelt);

  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {
    cs_lnum_t face_id = lstelt[ilelt];
    lagr_bdy_cond->b_face_zone_id[face_id] = 4; /* zone_id */
  }

  /* Injection per particle class into the calculation domain
     ======================================================== */

  /* To provide information about the particle classes,
     we follow a two-step procedure:
     1) first, the number of particle classes is prescribed
     for each boundary zone: iusncl (by default, this parameter is set to zero)

     2) afterwards, for each zone and for each class, we prescribe
     the physical properties of the particles
     Number of particle classes entering the domain
     We assign here the number of classes for each zone previously identified.

     This number is zero by default.
     The maximal number of classes is nclagm
     ---> First zone numbered izone = 0: 1 class injected  */

  int nbclas = 1;
  lagr_bdy_cond->b_zone_classes[0] = nbclas;

  /* ---> Second zone numbered izone = 1: 0 class injected */
  nbclas = 0;
  lagr_bdy_cond->b_zone_classes[1] = nbclas;

  /* ---> Third zone numbered izone = 3 : 0 class injected */
  nbclas = 0;
  lagr_bdy_cond->b_zone_classes[3] = nbclas;

  /* ---> Zone numbered izone = 4 : 0 class injected */
  nbclas = 0;
  lagr_bdy_cond->b_zone_classes[4] = nbclas;

  /* For every class associated with a zone,
     we provide the following information:

     b_zone_classes number of classes per zone
     b_zone_natures boundary conditions for the particles
       = CS_LAGR_INLET     -> zone of particle inlet
       = CS_LAGR_OUTLET    -> particle outlet
       = CS_LAGR_REBOUND   -> rebound of the particles
       = CS_LAGR_DEPO1     -> definitive deposition
       = CS_LAGR_DEPO2     -> definitive deposition, but the particle remains
                              in memory (useful only if iensi2 = 1)
       = CS_LAGR_DEPO_DLVO -> deposition of the particle with DLVO forces
       = CS_LAGR_FOULING   -> fouling (coal only physical_model = 2)
       = CS_LAGR_SYM       -> symmetry condition for the particles (zero flux)

     *   nb_part : number of particles per class and per zone
     *   injection_frequency : injection frequency. If injection_frequency = 0,
                               then the injection occurs only at the first
                               absolute iteration.

     *   cluster : number of the group to which the particle belongs
                   (only if one wishes to calculate statistics per group)

     *   velocity_profile : type of condition on the velocity
              = -1 imposed flow velocity
              =  0 imposed velocity along the normal direction of the
                    boundary face, with norm equal to velocity[0] (m/s)
              =  1 imposed velocity: we prescribe velocity[0] (m/s)
                                                  velocity[1] (m/s)
                                                  velocity[2] (m/s)
              =  2 user-defined profile

     *   distribution_profile : type of statistical weight
              = 1 automatic: we prescribe
                    flow_rate:   mass flow rate (kg/s)
                    stat_weight: statistical weight (number of samples) associated
                                 to the particle (automatically computed to respect
                                 a mass flow rate if it is defined)
              = 2 user-defined profile

     *   temperature_profile : type of temperature condition
              = 1 imposed temperature: we prescribe temperature[] in Kelvin:
                        - temperature[0]       if physical_model != 2
                        - temperature[nlayers] if physical_model = 2
              = 2 user-defined profile

     *   diameter_profile : type of diameter condition
              = 1 imposed diameter: we prescribe diameter (m) and diameter_variance
                                    (standard deviation, in m)
              = 2 user-defined profile
     *   coal_number : number of the coal of the particle (only if physical_model = 2)

     *   coal_profile : type of coal injection composition (only if physical_model = 2)
              = 0 coal injected with an user-defined composition
              = 1 raw coal injection

     density is 2500.0 */

  int  izone = 0;
  nbclas = lagr_bdy_cond->b_zone_classes[izone];
  lagr_bdy_cond->b_zone_natures[izone] = CS_LAGR_INLET;

  for (cs_lnum_t iclas = 0; iclas < nbclas; iclas++) {

    /* Ensure defaults are set   */
    cs_lagr_zone_class_data_t *zone_class_data = cs_lagr_init_zone_class_new(iclas, izone);

    /* Now define parameters for this class and zone */

    zone_class_data->nb_part = 10;
    zone_class_data->injection_frequency = 2;
    if (cs_glob_lagr_model->n_stat_classes > 0)
      zone_class_data->cluster = iclas + 1;

    int vel_profile = 0;
    cs_real_t vel[3];
    vel[0] = 1.1;
    vel[1] = 0.0;
    vel[2] = 0.0;

    cs_lagr_set_zone_class_velocity(iclas,
                                    izone,
                                    vel_profile,
                                    vel);

    int stat_profile = 1;
    cs_real_t stat_weight = 1.0;
    cs_real_t flow_rate   = 0.0;
    cs_lagr_set_zone_class_stat(iclas,
                                izone,
                                stat_profile,
                                stat_weight,
                                flow_rate);

    /* if the physics is " simple"    */

    if (   cs_glob_lagr_model->physical_model == 0
        || cs_glob_lagr_model->physical_model == 1) {

      /* Mean value and standard deviation of the diameter */
      int diam_profile = 1;
      cs_real_t diam = 5e-05;
      cs_real_t diam_dev  = 0.0;
      cs_lagr_set_zone_class_diam(iclas,
                                  izone,
                                  diam_profile,
                                  diam,
                                  diam_dev);

      /* Density    */
      cs_real_t density = 2500.0;
      cs_lagr_set_zone_class_density(iclas,
                                     izone,
                                     density);

      cs_real_t foul_index = 100.0;
      cs_lagr_set_zone_class_foul_index(iclas,
                                        izone,
                                        foul_index);

      if (cs_glob_lagr_model->physical_model == 1) {

        /* Temperature and Cp   */
        if (cs_glob_lagr_specific_physics->itpvar == 1) {

          int temp_profile  = 1;
          cs_real_t temp[1] = {20.0};
          cs_real_t cp    = 1400.0;
          cs_real_t emissivity   = 0.7;
          cs_lagr_set_zone_class_temperature(iclas,
                                             izone,
                                             temp_profile,
                                             temp,
                                             emissivity);

          cs_lagr_set_zone_class_cp(iclas,
                                    izone,
                                    cp);

        }

      }

    }

    /* Coal  */
    else if (cs_glob_lagr_model->physical_model == 2) {
      int nlayer = cs_glob_lagr_model->n_temperature_layers;
      /* CAUTION :
         1) To transport and burn coal particles with the Lagrangian
         module, a specific physics for the dispersed phase must
         be activated for the carrier phase.
         2) The physical properties of the coal particles are known
         from the thermo-chemical file: dp_FCP
         3) For the current phase ICLAS, and for the current boundary
         zone IZONE, we assign to the coal particles the properties of
         the coal ICOAL of the ICOAL class taken from the file dp_FCP.
         4) icoal : number of the coal between 1 and ncharb defined by
         the user in the file dp_FCP.
         Mean value and standard deviation of the diameter  */

      int diam_profile = 1;
      cs_real_t diam = 5e-05;
      cs_real_t diam_dev  = 0.0;
      cs_lagr_set_zone_class_diam(iclas,
                                  izone,
                                  diam_profile,
                                  diam,
                                  diam_dev);

      /* Temperature (in K)   */
      cs_real_t temp[nlayer] ;
      for (int ilayer = 0; ilayer < nlayer; ilayer++)
        temp[ilayer] = 800.0 ;
      cs_real_t cp;

      /* Number of the coal   */
      int icoal = 0;

      /* Raw coal or user defined coal injection condition  */
      int coal_profile = 1;
      if (coal_profile == 0) {

        /* Example of user-defined injection, coal after devolatilisation */
        /* Specific heat   */
        cp = lag_cc->cp2ch[icoal];

        /* water mass fraction in the particle */
        cs_real_t water_mass_fraction  = 0.0;

        /* Density    */
        cs_real_t density
          =   lag_cc->xashch[icoal] * lag_cc->rho0ch[icoal]
            + (1.0 - lag_cc->xwatch[icoal] - lag_cc->xashch[icoal])
            * lag_cc->rho0ch[icoal]
            * (1.0 - (lag_cc->y1ch[icoal] + lag_cc->y2ch[icoal]) / 2.0);

        cs_real_t coal_mass_fraction[nlayer];
        cs_real_t coke_mass_fraction[nlayer];
        cs_real_t coke_density[nlayer];

        for (int ilayer = 0; ilayer < nlayer; ilayer++) {
          /* reactive coal mass fraction in the particle   */
          coal_mass_fraction[ilayer] = 0.0;

          /* coke density after pyrolysis   */
          coke_density[ilayer]
            =   (1.0 - lag_cc->xwatch[icoal] - lag_cc->xashch[icoal])
              * lag_cc->rho0ch[icoal]
              * (1.0 - (lag_cc->y1ch[icoal] + lag_cc->y2ch[icoal]) / 2.0);

          /* coke mass fraction in the particle  */
          coke_mass_fraction[ilayer] = coke_density[ilayer] / density;

        }

        /* coke diameter   */
        cs_real_t shrinking_diameter = diam;

        /* initial particle diameter */
        cs_real_t initial_diameter = diam;

        cs_lagr_set_zone_class_coal(iclas,
                                    izone,
                                    coal_profile,
                                    icoal,
                                    temp,
                                    coal_mass_fraction,
                                    coke_mass_fraction,
                                    coke_density,
                                    water_mass_fraction,
                                    shrinking_diameter,
                                    initial_diameter);

      }

      cs_lagr_set_zone_class_cp(iclas,
                                izone,
                                cp);

    }

    /* Complete definition of parameters for this class and zone    */
    //    lagr_define_zone_class_param (&iclas, &izone, iczpar, rczpar);
    /* ===============================     */
  }

  /* ---> Second zone, numbered izone = 1 */
  /* IUSCLB : rebound of the particle    */
  izone  = 1;
  lagr_bdy_cond->b_zone_natures[izone] = CS_LAGR_REBOUND;

  /* same procedure for the other zones...    */

  /* Deallocate the temporary array */
  BFT_FREE(lstelt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
