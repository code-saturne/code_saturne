/*============================================================================
 * Methods for lagrangian module
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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_math.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_thermal_model.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_physical_constants.h"
#include "cs_time_step.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_random.h"
#include "cs_prototypes.h"

#include "cs_gui_particles.h"
#include "cs_gui_util.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_geom.h"
#include "cs_lagr_new.h"
#include "cs_lagr_precipitation_model.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_injection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inject particles in the computational domain.
 *
 * \param[in] time_id     time step indicator for fields
 *                         0: use fields at current time step
 *                         1: use fields at previous time step
 * \param[in] itypfb      boundary face types
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_injection(int        time_id,
                  const int  itypfb[],
                  cs_real_t  vislen[])
{
  cs_real_t dnbpnw_preci = 0.;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  cs_real_t  *xashch = cs_glob_lagr_coal_comb->xashch;
  cs_real_t  *cp2ch  = cs_glob_lagr_coal_comb->cp2ch;
  cs_real_t  *xwatch = cs_glob_lagr_coal_comb->xwatch;
  cs_real_t  *rho0ch = cs_glob_lagr_coal_comb->rho0ch;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  cs_real_t tkelvn = -cs_physical_constants_celsius_to_kelvin;

  /* Non-lagrangian fields */
  cs_real_t *vela = extra->vel->vals[time_id];
  cs_real_t *cscalt = NULL, *temp = NULL, *temp1 = NULL;

  cs_lagr_particle_counter_t *pc = cs_lagr_get_particle_counter();
  const cs_time_step_t *ts = cs_glob_time_step;

  if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
      || cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY)
    cscalt = extra->scal_t->vals[time_id];

  if (extra->temperature != NULL)
    temp = extra->temperature->val;

  if (extra->t_gaz != NULL)
    temp1 = extra->t_gaz->val;

  /* Memory management */

  int *ilftot;
  BFT_MALLOC(ilftot, cs_glob_lagr_const_dim->nflagm, int);

  /* Initialization */

  cs_real_t pis6 = cs_math_pi / 6.0;

  cs_lagr_bdy_condition_t  *bdy_cond = cs_lagr_get_bdy_conditions();
  cs_lagr_get_internal_conditions();

  /* User initialization by class and boundary */

  if (cs_gui_file_is_loaded() && ts->nt_cur == ts->nt_prev +1)
    cs_gui_particles_bcs();

  cs_user_lagr_boundary_conditions(itypfb);

  /* setup BCs */

  /* Faces must all belong to a valid boundary zone */

  int *ifvu;
  BFT_MALLOC(ifvu, cs_glob_lagr_const_dim->nflagm, int);
  for (int ii = 0; ii < cs_glob_lagr_const_dim->nflagm; ii++)
    ifvu[ii] = 0;

  for (cs_lnum_t ifac = 0; ifac < mesh->n_b_faces; ifac++) {
    if (   bdy_cond->b_face_zone_id[ifac] < 0
        || bdy_cond->b_face_zone_id[ifac] >= cs_glob_lagr_const_dim->nflagm)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian module: \n"
                  "  the zone id associated to face %d must be "
                  "an integer >= 0 and < nflagm = %d.\n"
                  "  This number is here %d."),
                (int)ifac,
                (int)cs_glob_lagr_const_dim->nflagm,
                (int)bdy_cond->b_face_zone_id[ifac]);
    else
      ifvu[bdy_cond->b_face_zone_id[ifac]] = 1;
  }

  cs_parall_max(cs_glob_lagr_const_dim->nflagm, CS_INT_TYPE, ifvu);

  /* Build a list of boundary zone ids. */

  bdy_cond->n_b_zones = 0;

  for (int ii = 0; ii < cs_glob_lagr_const_dim->nflagm; ii++) {
    if (ifvu[ii]) {
      bdy_cond->b_zone_id[bdy_cond->n_b_zones] = ii;
      bdy_cond->n_b_zones += 1;
    }
  }

  BFT_FREE(ifvu);

  /* --> Calculation of the surfaces of the Lagrangian boundary zones  */

  cs_real_t *surflag = NULL;
  cs_real_t *surlgrg = NULL;
  cs_lnum_t *ninjrg = NULL;
  cs_lnum_t  nfrtot;

  if (cs_glob_rank_id >= 0) {

    BFT_MALLOC(surflag, cs_glob_lagr_const_dim->nflagm, cs_real_t);
    BFT_MALLOC(surlgrg, cs_glob_lagr_const_dim->nflagm * cs_glob_n_ranks, cs_real_t);
    BFT_MALLOC(ninjrg, cs_glob_n_ranks, cs_lnum_t);

    for (cs_lnum_t kk = 0; kk < cs_glob_lagr_const_dim->nflagm; kk++) {

      surflag[kk] = 0.0;

      for (cs_lnum_t jj = 0; jj < cs_glob_n_ranks; jj++)
        surlgrg[cs_glob_n_ranks * kk + jj] = 0.0;

    }

    for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {

      int izone = bdy_cond->b_zone_id[ii];

      surflag[izone] = 0.0;

      for (cs_lnum_t ifac = 0; ifac < mesh->n_b_faces; ifac++) {

        if (izone == bdy_cond->b_face_zone_id[ifac]) {

          surflag[izone] += fvq->b_face_surf[ifac];
          surlgrg[(izone) * cs_glob_n_ranks + cs_glob_rank_id]
            += fvq->b_face_surf[ifac];

        }

      }

    }

    for (cs_lnum_t kk = 0; kk < cs_glob_lagr_const_dim->nflagm; kk++) {

      cs_parall_sum(1, CS_DOUBLE, &(surflag[kk]));

      for (cs_lnum_t jj = 0; jj < cs_glob_n_ranks; jj++)
        cs_parall_sum(1, CS_DOUBLE, &(surlgrg[kk * cs_glob_n_ranks + jj]));

    }

    if (cs_glob_rank_id == 0) {

      nfrtot = 0;
      cs_lnum_t jj = 0;

      for (cs_lnum_t kk = 0; kk < cs_glob_lagr_const_dim->nflagm; kk++) {

        if (surflag[kk] > 1e-15) {

          nfrtot++;
          ilftot[jj] = kk;
          jj++;

        }

      }

    }

    cs_parall_bcast(0, 1, CS_LNUM_TYPE, &nfrtot);
    cs_parall_bcast(0, nfrtot, CS_LNUM_TYPE, ilftot);

  }

  else {

    nfrtot = bdy_cond->n_b_zones;

    for (int ii = 0; ii < bdy_cond->n_b_zones; ii++)
      ilftot[ii] = bdy_cond->b_zone_id[ii];

  }

  /* --> Nombre de classes.    */

  for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {

    int izone = bdy_cond->b_zone_id[ii];

    if (bdy_cond->b_zone_classes[izone] < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian module: \n"
                  "  number of particle classes for zone %d "
                  "is not defined (=%d)\n"),
                (int)izone + 1,
                (int)bdy_cond->b_zone_classes[izone]);

  }

  /* Verification of particle classes de particules: just a warning */
  if (cs_glob_lagr_model->n_stat_classes > 0) {

    for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {

      int izone = bdy_cond->b_zone_id[ii];

      for (int iclas = 0;
           iclas < bdy_cond->b_zone_classes[izone];
           iclas++) {

        cs_lagr_zone_class_data_t *userdata
          = cs_lagr_get_zone_class_data(iclas, izone);

        if (   userdata->cluster <= 0
            || userdata->cluster > cs_glob_lagr_model->n_stat_classes)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian module: \n"
                      "  number of clusters = %d is either not defined (negative)\n"
                      "   or > to the number of statistical classes %d "
                      "for zone %d and class %d."),
                    (int)userdata->cluster,
                    (int)cs_glob_lagr_model->n_stat_classes,
                    (int)izone + 1,
                    (int)iclas);

      }

    }

  }

  /* --> Boundary conditions */
  for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {

    int izone = bdy_cond->b_zone_id[ii];

    if (   bdy_cond->b_zone_natures[izone] != CS_LAGR_INLET
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_OUTLET
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_REBOUND
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_DEPO1
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_DEPO2
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_SYM
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_FOULING
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_JBORD1
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_JBORD2
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_JBORD3
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_JBORD4
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_JBORD5
        && bdy_cond->b_zone_natures[izone] != CS_LAGR_DEPO_DLVO) {

      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian boundary zone %d nature %d is unknown."),
                izone + 1,
                (int)bdy_cond->b_zone_natures[izone]);
    }

  }

  for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {
    int izone = bdy_cond->b_zone_id[ii];
    if (   bdy_cond->b_zone_natures[izone] == CS_LAGR_FOULING
        && cs_glob_lagr_model->physical_model != 2)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian boundary zone %d nature is of type CS_LAGR_FOULING,\n"
           "but cs_glob_lagr_model->physical_model is not equal to 2."),
         (int)izone + 1);
  }

  for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {
    int izone = bdy_cond->b_zone_id[ii];
    if (   bdy_cond->b_zone_natures[izone] == CS_LAGR_FOULING
        && cs_glob_lagr_model->fouling != 1)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian boundary zone %d nature is of type CS_LAGR_FOULING,\n"
           "but fouling is not activated."),
         (int)izone + 1);
  }

  /* --> Type de condition pour le diametre.  */

  if (   cs_glob_lagr_model->physical_model == 1
      && (   cs_glob_lagr_specific_physics->itpvar == 1
          || cs_glob_lagr_specific_physics->idpvar == 1
          || cs_glob_lagr_specific_physics->impvar == 1)) {

    for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {
      int izone = bdy_cond->b_zone_id[ii];
      for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {
        cs_lagr_zone_class_data_t *userdata = cs_lagr_get_zone_class_data(iclas, izone);
        if (   userdata->temperature_profile < 1
            || userdata->temperature_profile > 2)
          bft_error
            (__FILE__, __LINE__, 0,
             _("Lagrangian boundary zone %d temperature profile "
               "value is invalid (=%d)\n"),
             (int)izone + 1,
             (int)userdata->temperature_profile);
      }
    }

  }

  /* local copy of parameters */

  cs_lagr_zone_class_data_t *local_zone_class_data = NULL;
  BFT_MALLOC(local_zone_class_data,
             cs_glob_lagr_nzone_max * cs_glob_lagr_nclass_max,
             cs_lagr_zone_class_data_t);

  for (int izone = 0; izone < cs_glob_lagr_nzone_max; izone++) {

    for (int iclas = 0; iclas < cs_glob_lagr_nclass_max; iclas++) {

      cs_lagr_zone_class_data_t *userdata
        = cs_lagr_get_zone_class_data(iclas, izone);

      local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone] = *userdata;

    }

  }

  for (int ii = 0; ii < bdy_cond->n_b_zones; ii++) {

    int izone = bdy_cond->b_zone_id[ii];

    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *userdata = cs_lagr_get_zone_class_data(iclas, izone);

      if (userdata->nb_part < 0)
        bft_error
          (__FILE__, __LINE__, 0,
           _("Number of particles for zone %d and class %d is not defined (=%d)\n"),
           (int)izone + 1,
           (int)iclas,
           (int)bdy_cond->b_zone_classes[izone]);

      /* --> Type de condition pour le taux de presence.    */
      if (   userdata->distribution_profile < 1
          || userdata->distribution_profile > 2)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d, class %d:\n"
                    "  distribution profile value is invalid (=%d)\n"),
                  (int)izone + 1,
                  (int)iclas,
                  (int)userdata->distribution_profile);

      /* --> Type de condition pour la vitesse.   */
      if (   userdata->velocity_profile <  -1
          || userdata->velocity_profile > 2)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d, class %d:\n"
                    "  velocity profile value is invalid (=%d)\n"),
                  (int)izone + 1,
                  (int)iclas,
                  (int)userdata->velocity_profile);

      /* --> Type de condition pour le diametre.  */
      if (   userdata->diameter_profile < 1
          || userdata->diameter_profile > 2)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d:\n"
                    "  class %d diameter profile value is invalid (=%d)\n"),
                  (int)izone + 1,
                  (int)iclas,
                  (int)userdata->diameter_profile);

      /* statistical weight */
      if (   userdata->stat_weight <= 0.0
          && userdata->flow_rate <= 0.0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d, class %d:\n"
                    "  statistical weight value is invalid (=%e10.3)\n"),
                  (int)izone + 1,
                  (int)iclas,
                  (double)userdata->stat_weight);

      /* particle mass flow rate */
      if (   userdata->flow_rate > 0.0
          && userdata->nb_part  == 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d, class %d:\n"
                    " flow rate is positive (%e10.3)\n"
                    " while number of particle is 0."),
                  (int)izone + 1,
                  (int)iclas,
                  (double)userdata->flow_rate);

      /* --> Proprietes des particules : le diametre, son ecart-type, et rho    */
      if (cs_glob_lagr_model->physical_model != 2) {
        if (   userdata->density  < 0.0
            || userdata->diameter < 0.0
            || userdata->diameter_variance < 0.0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian boundary zone %d, class %d:\n"
                      "  error on particles properties definition:\n"
                      "  rho = %e10.3, diameter = %e10.3,\n"
                      "  diameter standard deviation = %e10.3."),
                    (int)izone + 1,
                    (int)iclas,
                    (double)userdata->density,
                    (double)userdata->diameter,
                    (double)userdata->diameter_variance);
      }

      if (userdata->diameter < 3.0 * userdata->diameter_variance)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d, class %d:\n"
                    "  diameter (%e10.3) is smaller than 3 times\n"
                    "  its standard deviation (%e10.3)."),
                  (int)izone + 1,
                  (int)iclas,
                  (double)userdata->diameter,
                  (double)userdata->diameter_variance);

      /* --> Proprietes des particules : Temperature et CP  */
      if (   cs_glob_lagr_model->physical_model == 1
          && cs_glob_lagr_specific_physics->itpvar == 1) {

        if (   userdata->cp < 0.0
            || userdata->temperature[0] < tkelvn)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian boundary zone %d, class %d:\n"
                      "  specific heat capacity is negative (%e10.3)\n"
                      "  or temperature (%e10.3) is lower than %e10.3."),
                    (int)izone + 1,
                    (int)iclas,
                    (double)userdata->cp,
                    (double)userdata->temperature[0],
                    (double)tkelvn);

      }

      /* --> Proprietes des particules : Emissivite    */
      if (   cs_glob_lagr_model->physical_model == 1
          && cs_glob_lagr_specific_physics->itpvar == 1
          && extra->iirayo > 0) {

        if (   userdata->emissivity < 0.0
            || userdata->emissivity > 1.0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian boundary zone %d, class %d:\n"
                      "  particle emissivity is not properly set = %e10.3."),
                    (int)izone + 1,
                    (int)iclas,
                    (double)userdata->emissivity);

      }

      /* Charbon    */
      if (cs_glob_lagr_model->physical_model == 2) {

        if (   userdata->coal_number < 1
            && userdata->coal_number > extra->ncharb)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian boundary zone %d, class %d:\n"
                      "  the coal number %d for the injected particle is either negative\n"
                      "  or over the maximum number of coal given in dp_FCP (ncharb = %d)."),
                    (int)izone + 1,
                    (int)iclas,
                    (int)userdata->coal_number,
                    (int)extra->ncharb);

        int coal_id = userdata->coal_number - 1;

        /* Properties of coal particles */
        if (   userdata->coal_profile < 0
            || userdata->coal_profile > 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian boundary zone %d, class %d:\n"
                      "  the coal profile flag is not properly set.\n"
                      "  I must be equal to 0 (user definition) or 1 "
                      " (composition is set identical to fresh coal).\n"
                      "  Coal profile is %d."),
                    (int)izone + 1,
                    (int)iclas,
                    (int)userdata->coal_profile);

        else if (   userdata->coal_profile     == 0
                 && userdata->diameter_profile == 2) {
          bft_printf(_("\nLagrangian module warning:\n"));
          bft_printf(_("  Both coal and diameter profiles are user\n"
                       " defined for zone %d, class %d.\n"
                       "  This may cause unwanted problems:\n"
                       "  coal initial and shrinking diameters are set before\n"
                       "  the definition of the particle diameter."),
                     (int)izone + 1,
                     (int)iclas);
        }

        else if (   userdata->coal_profile     == 0
                 && userdata->diameter_profile == 1
                 && userdata->diameter_variance > 0.0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian boundary zone %d, class %d:\n"
                      "  both coal and diameter profiles are set to be user defined\n"
                      "  This may cause unwanted problems: coal initial and shrinking\n"
                      "  diameters are set before the definition of the particle diameter."),
                    (int)izone + 1,
                    (int)iclas);

        for (int ilayer = 0;
             ilayer < cs_glob_lagr_model->n_temperature_layers;
             ilayer++) {

          if (userdata->temperature[ilayer] < tkelvi)
            bft_error(__FILE__, __LINE__, 0,
                      _("Lagrangian boundary zone %d, class %d:\n"
                        "  temperature is not properly set for layer %d: %e10.3."),
                      (int)izone + 1,
                      (int)iclas,
                      (int)ilayer,
                      (double)userdata->temperature[ilayer]);

        }

        /* Properties of coal particles */

        if (userdata->coal_profile == 0) {

          if (   userdata->density < 0.0
              || userdata->cp < 0.0
              || userdata->water_mass_fraction < 0.0
              || userdata->water_mass_fraction > 1.0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Lagrangian boundary zone %d, class %d:\n"
                        "  wrong conditions, with\n"
                        "    density = %e10.3\n"
                        "    Cp = %e10.3\n"
                        "    steam mass fraction = %e10.3."),
                      (int)izone + 1,
                      (int)iclas,
                      (double)userdata->density,
                      (double)userdata->cp,
                      (double)userdata->water_mass_fraction);

          for (int ilayer = 0;
               ilayer < cs_glob_lagr_model->n_temperature_layers;
               ilayer++) {

            if (   userdata->coal_mass_fraction[ilayer] < 0.0
                || userdata->coal_mass_fraction[ilayer] > 1.0
                || userdata->coke_mass_fraction[ilayer] < 0.0
                || userdata->coke_mass_fraction[ilayer] > 1.0
                || userdata->coke_density[ilayer] < 0.0)
              bft_error(__FILE__, __LINE__, 0,
                        _("Lagrangian boundary zone %d, class %d:\n"
                          "  wrong conditions on layer %d.\n"
                          "    coal mass fraction = %e10.3\n"
                          "    coke mass fraction = %e10.3\n"
                          "    coke density after pyrolysis = %e10.3."),
                        (int)izone + 1,
                        (int)iclas,
                        (int)ilayer,
                        (double)userdata->coal_mass_fraction[ilayer],
                        (double)userdata->coke_mass_fraction[ilayer],
                        (double)userdata->coke_density[ilayer]);

          }

          if (   userdata->shrinking_diameter < 0.0
              || userdata->initial_diameter < 0.0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Lagrangian boundary zone %d, class %d:\n"
                        "  wrong conditions, with\n"
                        "    coke diameter = %e10.3\n"
                        "    initial diameter = %e10.3\n"),
                      (int)izone + 1,
                      (int)iclas,
                      (double)userdata->shrinking_diameter,
                      (double)userdata->initial_diameter);

        }

        /* irawcl = 1 --> Composition du charbon definie dans le fichier XML (DP_FCP)
         * on verifie les donnes contenues dans le XML   */
        else if (userdata->coal_profile == 1) {

          if (   rho0ch[coal_id] < 0.0
              || cp2ch[coal_id]  < 0.0
              || xwatch[coal_id] < 0.0
              || xwatch[coal_id] > 1.0
              || xashch[coal_id] < 0.0
              || xashch[coal_id] > 1.0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Lagrangian boundary zone %d, class %d:\n"
                        "  wrong conditions for coal number %d.\n"
                        "    density RHO0CH = %e10.3\n"
                        "    Cp CP2CH = %e10.3\n"
                        "    water mass fraction XWATCH = %e10.3\n"
                        "    ashes mass fraction XASHCH = %e10.3."),
                      (int)izone + 1,
                      (int)iclas,
                      (int)coal_id,
                      (double)rho0ch[coal_id],
                      (double)cp2ch[coal_id],
                      (double)xwatch[coal_id],
                      (double)xashch[coal_id]);

          if (xwatch[coal_id] + xashch[coal_id] > 1.0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Lagrangian boundary zone %d, class %d:\n"
                        "  wrong conditions for coal number %d.\n"
                        "    water mass fraction XWATCH = %e10.3\n"
                        "    ashes mass fraction XASHCH = %e10.3\n"
                        "    mass fraction is larger than 1: %e10.3."),
                      (int)izone + 1,
                      (int)iclas,
                      (int)userdata->coal_number,
                      (double)xwatch[coal_id],
                      (double)xashch[coal_id],
                      (double)(  xwatch[coal_id]
                               + xashch[coal_id]));

          if (   userdata->density >= 0.0
              || userdata->water_mass_fraction >= 0.0
              || userdata->cp >= 0.0
              || userdata->shrinking_diameter >= 0.0
              || userdata->initial_diameter >= 0.0)
            bft_error
              (__FILE__, __LINE__, 0,
               _("Lagrangian boundary zone %d, class %d:\n"
                 "  conditions set with the contents of a DP_FCP file,\n"
                 "  but one is initialized to a different value than %e10.3.\n"
                 "  density = %e10.3\n"
                 "  water mass fraction = %e10.3\n"
                 "  Cp = %e10.3\n"
                 "  coke diameter = %e10.3\n"
                 "  initial diameter = %e10.3."),
               (int)izone + 1,
               (int)iclas,
               (double)-cs_math_big_r,
               (double)userdata->density,
               (double)userdata->water_mass_fraction,
               (double)userdata->cp,
               (double)userdata->shrinking_diameter,
               (double)userdata->initial_diameter);

          for (int ilayer = 0;
               ilayer < cs_glob_lagr_model->n_temperature_layers;
               ilayer++) {

            if (   userdata->coal_mass_fraction[ilayer] >= 0.0
                || userdata->coke_mass_fraction[ilayer] >= 0.0
                || userdata->coke_density[ilayer] >= 0.0)
              bft_error
                (__FILE__, __LINE__, 0,
                 _("Lagrangian boundary zone %d, class %d:\n"
                   "  wrong conditions for layer %d.\n"
                   "    Conditions are set with the contents of a DP_FCP file,\n"
                   "    but one is initialized to a different value than %e10.3.\n"
                   "  coal mass fraction IFRMCH = %e10.3\n"
                   "  coke mass fraction IFRMCK = %e10.3\n"
                   "  initial coke mass fraction IRHOCK0 = %e10.3\n"
                   "  coke diameter = %e10.3\n"
                   "  initial diameter = %e10.3."),
                 (int)izone + 1,
                 (int)iclas,
                 (int)ilayer,
                 (double)-cs_math_big_r,
                 (double)userdata->coal_mass_fraction[ilayer],
                 (double)userdata->coke_mass_fraction[ilayer],
                 (double)userdata->coke_density[ilayer],
                 (double)userdata->shrinking_diameter,
                 (double)userdata->initial_diameter);

          }

        }

      }

    }

  }

  /* ==============================================================================
   * 4. Transformation des donnees utilisateur
   * ============================================================================== */

  /* Compute number of particles to inject for this iteration */

  p_set->n_part_new = 0;

  for (int ii = 0; ii < nfrtot; ii++) {

    int izone = ilftot[ii];

    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *userdata
        = cs_lagr_get_zone_class_data(iclas, izone);
      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);

      /* Inject only at first time step if injection frequency is zero */

      if (local_userdata->injection_frequency <= 0) {
        if (ts->nt_cur == ts->nt_prev+1 && pc->n_g_cumulative_total == 0)
          local_userdata->injection_frequency = ts->nt_cur;
        else
          local_userdata->injection_frequency = ts->nt_cur+1;
      }

      if (ts->nt_cur % local_userdata->injection_frequency == 0)
        p_set->n_part_new += userdata->nb_part;

    }

  }

  /* ==============================================================================
   * 5. Precipitation/Dissolution
   * ============================================================================== */

  if (cs_glob_lagr_model->precipitation == 1)
    precdi(vela, &dnbpnw_preci);

  /* --> Limite du nombre de particules  */
  cs_lnum_t tmp = cs_lagr_particle_set_resize(p_set->n_particles + p_set->n_part_new);

  if (tmp < 0) {

    bft_printf(_("\n Lagrangian module: \n"));
    bft_printf
      (_("  If particles are injected according to boundary conditions,\n"
         "  the total number of particle in the domain would exceed the total\n"
         "  number admissible set by cs_lagr_set_n_g_particles_max.\n"
         "  No particle are injected at iteration %d."),
       ts->nt_cur);
    p_set->n_part_new = 0;

  }

  /* In no new particles are injected, return */

  if (p_set->n_part_new == 0) {

    BFT_FREE(surflag);
    BFT_FREE(surlgrg);
    BFT_FREE(ninjrg);

    BFT_FREE(ilftot);
    BFT_FREE(local_zone_class_data);

    return;
  }

  /* ----------------------------------------------------------------------
   * --> Tirage aleatoire des positions des P_SET->N_PART_NEW nouvelles particules
   * au niveau des zones de bord et reperage des cellules correspondantes
   * ---------------------------------------------------------------------- */

  /* Initialize local number of particles injected by rank */
  cs_lnum_t nlocnew = 0;

  /* Distribute new particles  */
  /* For each boundary zone    */
  for (int ii = 0; ii < nfrtot; ii++) {

    int izone = ilftot[ii];

    /* for each class  */
    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *userdata = cs_lagr_get_zone_class_data(iclas, izone);
      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);

      /* if new particles must be added */
      if (ts->nt_cur % local_userdata->injection_frequency == 0) {

        /* Compute on rank 0 the number of particles to inject for each rank
         * based on the relative surface of each injection zone present on
         * each rank: fill array ninjrg(cs_glob_n_ranks) */
        if (cs_glob_rank_id == 0) {

          for (int irp = 0; irp < cs_glob_n_ranks; irp++)
            ninjrg[irp] = 0;

          for (cs_lnum_t ipart = 0; ipart < userdata->nb_part; ipart++) {

            cs_real_t random;
            cs_random_uniform(1, &random);

            /* blindage   */
            random = random + 1e-09;
            int irp = 0;
            cs_real_t offset = surlgrg[izone * cs_glob_n_ranks + irp] / surflag[izone];

            while (random > offset) {

              irp++;
              offset += surlgrg[izone * cs_glob_n_ranks + irp] / surflag[izone];

            }

            ninjrg[irp]++;

          }

        }

        /* Broadcast to all ranks */
        cs_parall_bcast(0, cs_glob_n_ranks, CS_LNUM_TYPE, ninjrg);

        /* End of the computation of the number of particles to inject */
        if (cs_glob_rank_id >= 0) {
          local_userdata->nb_part = ninjrg[cs_glob_rank_id];
          nlocnew += ninjrg[cs_glob_rank_id];
        }
        else {
          local_userdata->nb_part = userdata->nb_part;
          nlocnew += userdata->nb_part;
        }

      }

    }

  }

  tmp = cs_lagr_particle_set_resize(p_set->n_particles + nlocnew);

  if (tmp < 0) {

    bft_printf(_("\n Lagrangian module: \n"));
    bft_printf
      (_("If particles are injected according to boundary conditions,\n"
         " the total number of particle in the domain would exceed the\n"
         " total number admissible set by cs_lagr_set_n_g_particles_max.\n"
         " No particle are injected at iteration %d.\n"),
       (int)ts->nt_cur);

    BFT_FREE(surflag);
    BFT_FREE(surlgrg);
    BFT_FREE(ninjrg);

    BFT_FREE(ilftot);
    BFT_FREE(local_zone_class_data);

    return;
  }

  /* Allocate a work array     */
  cs_lnum_t *iwork;
  BFT_MALLOC(iwork, p_set->n_particles + nlocnew, cs_lnum_t);

  /* Now define particles */
  /* initialize new particles counter    */
  cs_lnum_t npt = p_set->n_particles;

  /* For each boundary zone    */
  for (int ii = 0; ii < nfrtot; ii++) {

    int izone = ilftot[ii];

    /* for each class  */
    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);

      /* if new particles must be added */
      if (ts->nt_cur % local_userdata->injection_frequency == 0) {

        if (local_userdata->nb_part > 0) {

          cs_lagr_new(&npt,
                      local_userdata->nb_part,
                      izone,
                      bdy_cond->b_face_zone_id,
                      iwork);

        }

      }

    }

  }

  /* Injected particles count check */

  if ((p_set->n_particles + nlocnew) != npt) {

    bft_printf(_("\n Lagrangian module: \n"));
    bft_printf(_("  Bad boundary conditions.\n The number of injected particles\n"
                 "  for this time step does not match the one specified\n"
                 "  in the boundary conditions.\n"
                 "  number set for injection NBPNEW = %d\n"
                 "  number affectively injected NPT-NBPART = %d\n"),
               (int)nlocnew,
               (int)npt-p_set->n_particles);
    cs_exit(1);

  }

  /* reinitialisation du compteur de nouvelles particules */
  npt = p_set->n_particles;

  /* pour chaque zone de bord: */
  for (int ii = 0; ii < nfrtot; ii++) {

    int izone = ilftot[ii];

    /* pour chaque classe : */
    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);
      cs_lagr_zone_class_data_t *userdata
        = cs_lagr_get_zone_class_data(iclas, izone);

      /* si de nouvelles particules doivent entrer */
      if (ts->nt_cur % local_userdata->injection_frequency == 0) {

        for (cs_lnum_t ip = npt; ip < npt + local_userdata->nb_part; ip++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

          cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);
          cs_lnum_t ifac = iwork[ip];

          /* Random value associated with each particle */

          cs_real_t part_random = -1;
          cs_random_uniform(1, &part_random);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RANDOM_VALUE,
                                    part_random);

          /* Particle velocity components */

          cs_real_t *part_vel = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);

          /* prescribed components */
          if (userdata->velocity_profile == 1) {
            for (cs_lnum_t i = 0; i < 3; i++)
              part_vel[i] = userdata->velocity[i];
          }

          /* prescribed norm */
          else if (userdata->velocity_profile == 0) {
            for (cs_lnum_t i = 0; i < 3; i++)
              part_vel[i] = -   fvq->b_face_normal[ifac * 3 + i]
                              / fvq->b_face_surf[ifac]
                              * userdata->velocity_magnitude;
          }

          /* velocity as seen from fluid */
          else if (userdata->velocity_profile ==  -1) {
            for (cs_lnum_t i = 0; i < 3; i++)
              part_vel[i] = vela[cell_id * 3  + i];
          }

          /* prescribed velocity profile */
          else if (userdata->velocity_profile == 2) {
            cs_user_lagr_new_p_attr(particle,
                                    p_am,
                                    ifac,
                                    CS_LAGR_VELOCITY);
          }

          /* fluid velocity seen */
          cs_real_t *part_seen_vel = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY_SEEN);
          for (int i = 0; i < 3; i++)
            part_seen_vel[i] = vela[cell_id * 3 + i];

          /* Residence time (may be negative to ensure continuous injection) */
          if (local_userdata->injection_frequency == 1) {
            cs_real_t res_time = - part_random *cs_glob_lagr_time_step->dtp;
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RESIDENCE_TIME,
                                      res_time);
          }
          else
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RESIDENCE_TIME,
                                      0.0);

          /* Diameter */
          if (userdata->diameter_profile == 1) {

            if (userdata->diameter_variance > 0.0) {

              double    random;
              cs_random_normal(1, &random);

              cs_real_t diam =   userdata->diameter
                               + random * userdata->diameter_variance;
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, diam);

              /* On verifie qu'on obtient un diametre dans la gamme des 99,7% */
              cs_real_t d3 = 3.0 * userdata->diameter_variance;

              if (  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER)
                  < userdata->diameter - d3)
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                                          userdata->diameter);

              if (  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER)
                  > userdata->diameter + d3)
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                                          userdata->diameter);

            }

            else
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                                        userdata->diameter);

          }

          /* in case of profile for diameter */
          else if (userdata->diameter_profile == 2) {
            cs_user_lagr_new_p_attr(particle,
                                    p_am,
                                    ifac,
                                    CS_LAGR_DIAMETER);
          }

          /* Other parameters */
          cs_real_t diam = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);
          cs_real_t mporos = cs_glob_lagr_clogging_model->mporos;
          if (cs_glob_lagr_model->clogging == 1) {
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                                      diam/(1.-mporos));
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, diam);
          }

          /* Other variables (mass, ...) depending on physical model  */
          cs_real_t d3 = pow(diam, 3.0);

          if (cs_glob_lagr_model->n_stat_classes > 0)
            cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_STAT_CLASS,
                                      userdata->cluster);

          /* used for 2nd order only */
          if (p_am->displ[0][CS_LAGR_TAUP_AUX] > 0)
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TAUP_AUX, 0.0);

          if (   cs_glob_lagr_model->physical_model == 0
              || cs_glob_lagr_model->physical_model == 1) {

            if (cs_glob_lagr_model->clogging == 0)
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                        userdata->density * pis6 * d3);
            else
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                        userdata->density * pis6 * d3
                                        * pow(1.0-mporos,3));

            if (   cs_glob_lagr_model->physical_model == 1
                && cs_glob_lagr_specific_physics->itpvar == 1) {

              /* si Temperature constante imposee :  */
              if (userdata->temperature_profile == 1)
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE,
                                          userdata->temperature[0]);

              /* si profil pour la temperature :     */
              else if (userdata->temperature_profile == 2) {
                cs_user_lagr_new_p_attr(particle,
                                        p_am,
                                        ifac,
                                        CS_LAGR_TEMPERATURE);
              }

              if (   cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0
                  || cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 0
                  || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
                cs_lagr_particle_set_real(particle, p_am,
                                          CS_LAGR_FLUID_TEMPERATURE,
                                          temp1[cell_id] - tkelvi);

              else if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0
                       || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0
                       || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 0
                       || cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 0)
                cs_lagr_particle_set_real(particle, p_am,
                                          CS_LAGR_FLUID_TEMPERATURE,
                                          temp[cell_id] - tkelvi);

              else if (cs_glob_thermal_model->itherm ==
                               CS_THERMAL_MODEL_TEMPERATURE) {

                /* Kelvin */
                if (cs_glob_thermal_model->itpscl == 1)
                  cs_lagr_particle_set_real(particle, p_am,
                                            CS_LAGR_FLUID_TEMPERATURE,
                                            cscalt[cell_id] - tkelvi);

                /* Celsius    */
                else if (cs_glob_thermal_model->itpscl == 2)
                  cs_lagr_particle_set_real(particle, p_am,
                                            CS_LAGR_FLUID_TEMPERATURE,
                                            cscalt[cell_id]);

              }

              else if (cs_glob_thermal_model->itherm ==
                               CS_THERMAL_MODEL_ENTHALPY) {

                int mode = 1;
                temp[0] = cs_lagr_particle_get_real(particle, p_am,
                                                    CS_LAGR_FLUID_TEMPERATURE);
                CS_PROCF(usthht, USTHHT)(&mode, &(cscalt[cell_id]), temp);

              }

              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP,
                                        userdata->cp);
              if (extra->iirayo > 0)
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_EMISSIVITY,
                                          userdata->emissivity);

            }

          }

          else if (cs_glob_lagr_model->physical_model == 2) {

            int coal_id = userdata->coal_number - 1;

            cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_COAL_NUM, coal_id);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FLUID_TEMPERATURE,
                                      temp1[cell_id] - tkelvi);

            cs_real_t *particle_temp
              = cs_lagr_particle_attr(particle, p_am, CS_LAGR_TEMPERATURE);
            for (int ilayer = 0;
                 ilayer < cs_glob_lagr_model->n_temperature_layers;
                 ilayer++)
              particle_temp[ilayer] = userdata->temperature[ilayer];

            /* user-defined composition (cs_user_lagr_boundary_conditions) */
            if (userdata->coal_profile == 0) {

              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP, userdata->cp);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                        userdata->density * pis6 * d3);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_WATER_MASS,
                                          userdata->water_mass_fraction
                                        * cs_lagr_particle_get_real(particle, p_am,
                                                                    CS_LAGR_MASS));

              cs_real_t *particle_coal_mass
                = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COAL_MASS);
              cs_real_t *particle_coke_mass
                = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COKE_MASS);
              for (int ilayer = 0;
                   ilayer < cs_glob_lagr_model->n_temperature_layers;
                   ilayer++) {

                particle_coal_mass[ilayer] = userdata->coal_mass_fraction[ilayer]
                  * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                  / cs_glob_lagr_model->n_temperature_layers;
                particle_coke_mass[ilayer] = userdata->coke_mass_fraction[ilayer]
                  * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                  / cs_glob_lagr_model->n_temperature_layers;

              }

              cs_lagr_particle_set_real(particle, p_am,
                                        CS_LAGR_SHRINKING_DIAMETER,
                                        userdata->shrinking_diameter);
              cs_lagr_particle_set_real(particle, p_am,
                                        CS_LAGR_INITIAL_DIAMETER,
                                        userdata->initial_diameter);

              cs_real_t *particle_coal_density
                = cs_lagr_particle_attr(particle, p_am,
                                        CS_LAGR_COAL_DENSITY);
              for (int ilayer = 0;
                   ilayer < cs_glob_lagr_model->n_temperature_layers;
                   ilayer++)
                particle_coal_density[ilayer] = userdata->coke_density[ilayer];

            }

            /* composition from DP_FCP */
            else if (userdata->coal_profile == 1) {

              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP,
                                        cp2ch[coal_id]);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                        rho0ch[coal_id] * pis6 * d3);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_WATER_MASS,
                                        xwatch[coal_id]
                                        * cs_lagr_particle_get_real(particle, p_am,
                                                                    CS_LAGR_MASS));

              cs_real_t *particle_coal_mass
                = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COAL_MASS);
              cs_real_t *particle_coke_mass
                = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COKE_MASS);
              for (int ilayer = 0;
                   ilayer < cs_glob_lagr_model->n_temperature_layers;
                   ilayer++) {

                particle_coal_mass[ilayer]
                  =    (1.0 - xwatch[coal_id]
                            - xashch[coal_id])
                    * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                   / cs_glob_lagr_model->n_temperature_layers;
                particle_coke_mass[ilayer] = 0.0;

              }

              cs_lagr_particle_set_real
                (particle, p_am,
                 CS_LAGR_SHRINKING_DIAMETER,
                 cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER));
              cs_lagr_particle_set_real
                (particle, p_am,
                 CS_LAGR_INITIAL_DIAMETER,
                 cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER));

              cs_real_t *particle_coal_density
                = cs_lagr_particle_attr(particle, p_am,
                                        CS_LAGR_COAL_DENSITY);
              for (int ilayer = 0;
                   ilayer < cs_glob_lagr_model->n_temperature_layers;
                   ilayer++)
                particle_coal_density[ilayer] = rho0ch[coal_id];

            }
          }

          /* statistical weight */
          if (userdata->distribution_profile == 1)
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_STAT_WEIGHT,
                                      userdata->stat_weight);

          else if (userdata->distribution_profile == 2) {
            cs_user_lagr_new_p_attr(particle,
                                    p_am,
                                    ifac,
                                    CS_LAGR_STAT_WEIGHT);
          }

          /* Fouling index */
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FOULING_INDEX,
                                    userdata->foul_index);

          /* Initialization of deposition model */

          if (cs_glob_lagr_model->deposition == 1) {

            cs_real_t random;
            cs_random_uniform(1, &random);
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_INTERF, 5.0 + 15.0 * random);
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_YPLUS, 1000.0);
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_MARKO_VALUE, -1);
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_NEIGHBOR_FACE_ID, -1);

          }

          /* Initialization of clogging model */

          if (cs_glob_lagr_model->clogging == 1) {

            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DEPO_TIME, 0.0);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT, 0.0);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART, 1.0);

          }

          /* Loop on the additional user variables */

          for (int i = 0;
              i < cs_glob_lagr_model->n_user_variables;
              i++) {

            cs_real_t *user_var = cs_lagr_particle_attr(particle, p_am, CS_LAGR_USER + i);
            *user_var = 0.;

            cs_user_lagr_new_p_attr(particle,
                                    p_am,
                                    ifac,
                                    CS_LAGR_USER + i);
          }

        }

        npt = npt + local_userdata->nb_part;

      }

    }

  }

  /* Control test */

  if ((p_set->n_particles + nlocnew) != npt)
    bft_error(__FILE__, __LINE__, 0,
              _("  Bad lagrangian boundary conditions.\n"
                "    the number of injected particles for this lagrangian\n"
                "    iteration does not match the one specified in the boundary\n"
                "    conditions.\n"
                "    number set for injection NBPNEW = %d\n"
                "    number effectively injected NPT-NBPART = %d"),
              (int)nlocnew,
              (int)npt-p_set->n_particles);

  /* Update weights to have the correct flow rate
     ============================================ */

  /* Reinitialize new particles counter */
  npt = p_set->n_particles;

  /* for each boundary zone */
  for (int z_id = 0; z_id < bdy_cond->n_b_zones; z_id++) {

    int izone = bdy_cond->b_zone_id[z_id];

    /* pour chaque classe : */
    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);
      cs_lagr_zone_class_data_t *userdata = cs_lagr_get_zone_class_data(iclas, izone);

      /* si de nouvelles particules sont entrees, */
      /* et si on a un debit non nul :  */
      if (   ts->nt_cur % local_userdata->injection_frequency == 0
          && userdata->flow_rate > 0.0
          && userdata->nb_part > 0) {

        cs_real_t dmass = 0.0;

        for (cs_lnum_t ip = npt; ip < npt + local_userdata->nb_part; ip++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
          dmass += cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

        }

        cs_parall_sum(1, CS_REAL_TYPE, &dmass);

        /* Compute weights */

        if (dmass > 0.0) {

          for (cs_lnum_t ip = npt; ip < npt + local_userdata->nb_part; ip++) {

            unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
            cs_lagr_particle_set_real
              (particle, p_am, CS_LAGR_STAT_WEIGHT,
               (userdata->flow_rate * cs_glob_lagr_time_step->dtp) / dmass);

          }

        }

        else {

          bft_printf(_("\n Lagrangian module: \n"));
          bft_printf
            (_(" In zone %d, class %d conditions are erroneous.\n"
               "   Imposed Flow rate value is =%e10.3 "
               "while number of particles is 0."),
             (int)izone + 1,
             (int)iclas,
             (double)userdata->flow_rate);
          cs_exit(1);

        }

        npt = npt + local_userdata->nb_part;

      }

    }

  }

  /* ==============================================================================
   * 6. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANEES VUES
   * PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
   * ============================================================================== */

  /* si de nouvelles particules doivent entrer :   */

  cs_lnum_t npar1 = p_set->n_particles;
  cs_lnum_t npar2 = p_set->n_particles + nlocnew;
  cs_lagr_new_particle_init(npar1, npar2, time_id, vislen);

  /* ==============================================================================
   * 7. MODIFICATION DES TABLEAUX DE DONNEES PARTICULAIRES
   * ============================================================================== */

  cs_user_lagr_in(time_id, iwork, local_zone_class_data, vislen);

  /* ==============================================================================
   * 7bis. Random id associated with particles (to be initialized later)
   * ============================================================================== */

  srand(cs_glob_rank_id + 1);

  /* for safety, build values at previous time step */
  for (cs_lnum_t ip = npar1; ip < npar2; ip++)
    cs_lagr_particles_current_to_previous(p_set, ip);

  /* reinitialisation du compteur de nouvelles particules    */
  npt = p_set->n_particles;

  /* pour chaque zone de bord: */
  for (int z_id = 0; z_id < bdy_cond->n_b_zones; z_id++) {

    int izone = bdy_cond->b_zone_id[z_id];

    /* loop on classes */
    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);
      cs_lagr_zone_class_data_t *userdata
        = cs_lagr_get_zone_class_data(iclas, izone);

      /* if new particles have been added */
      if (ts->nt_cur % local_userdata->injection_frequency == 0) {

        for (cs_lnum_t ip = npt; ip < npt + local_userdata->nb_part; ip++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

          if (   cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER) < 0.0
              && userdata->diameter_variance > 0.0){

            bft_printf(_("\nLagrangian module warning: \n"));
            bft_printf
              (_("  In zone %d, class %d conditions are erroneous.\n"
                 "    Computation of a particle diameter from mean diameter\n"
                 "    and standard deviation yields a negative value,\n"
                 "    due to a stochastic selection of 'gaussian border'\n"
                 "  mean diameter = %e10.3\n"
                 "  standard deviation = %e10.3\n"
                 "  computed diameter = %e10.3\n"),
               (int)izone + 1,
               (int)iclas,
               (double)userdata->diameter,
               (double)userdata->diameter_variance,
               (double)cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER));

          }

        }

        npt = npt + local_userdata->nb_part;

      }

    }

  }

  /* Control test */

  if ((p_set->n_particles + nlocnew) != npt)
    bft_error(__FILE__, __LINE__, 0,
              _("  Bad lagrangian boundary conditions.\n"
                "    the number of injected particles for this time step\n"
                "    does not match the one specified in the boundary conditions.\n"
                "  number set for injection NBPNEW = %d\n"
                "  number affectively injected NPT-NBPART = %d\n"),
              (int)nlocnew,
              (int)npt-p_set->n_particles);

  /* Free memory */
  BFT_FREE(iwork);

  /* ==============================================================================
   * 9. CALCUL DE LA MASSE TOTALE INJECTES EN CHAQUE ZONE
   * Attention cette valeur est modifie dans USLABO pour tenir compte
   * des particules qui sortent
   * + calcul du nombres physiques de particules qui rentrent (tenant
   * compte des poids)
   * ============================================================================== */

  /* reinitialisation du compteur de nouvelles particules    */
  npt = p_set->n_particles;
  p_set->weight_new = 0.0;

  /* pour chaque zone de bord :     */
  for (int z_id = 0; z_id < bdy_cond->n_b_zones; z_id++) {

    int izone = bdy_cond->b_zone_id[z_id];
    bdy_cond->particle_flow_rate[izone] = 0.0;

    /* pour chaque classe : */
    for (int iclas = 0; iclas < bdy_cond->b_zone_classes[izone]; iclas++) {

      cs_lagr_zone_class_data_t *local_userdata
        = &(local_zone_class_data[iclas * cs_glob_lagr_nzone_max + izone]);

      /* si de nouvelles particules sont entrees, */
      if (ts->nt_cur % local_userdata->injection_frequency == 0) {

        for (cs_lnum_t ip = npt; ip < npt + local_userdata->nb_part; ip++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

          bdy_cond->particle_flow_rate[izone]
            += (  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT)
                * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS));
          p_set->weight_new += cs_lagr_particle_get_real(particle, p_am,
                                                         CS_LAGR_STAT_WEIGHT);

        }

      }

      npt = npt + local_userdata->nb_part;

    }

  }

  if (cs_glob_lagr_model->precipitation == 1)
    p_set->weight_new += dnbpnw_preci;

  /* Update total number of particles */

  p_set->n_particles += nlocnew;
  p_set->n_part_new = nlocnew;

  pc = cs_lagr_update_particle_counter();
  pc->n_g_total += pc->n_g_new;

  /**************
   * Free memory
   **************/

  BFT_FREE(surflag);
  BFT_FREE(surlgrg);
  BFT_FREE(ninjrg);

  BFT_FREE(ilftot);
  BFT_FREE(local_zone_class_data);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
