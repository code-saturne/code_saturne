/*============================================================================
 * Lagrangian module logging
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

#include "cs_log.h"

#include "cs_math.h"

#include "cs_mesh.h"
#include "cs_parall.h"

#include "cs_field.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_post.h"
#include "cs_lagr_stat.h"

#include "cs_lagr_prototypes.h"

#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_log.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_log.c
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

const char *_astat[2] = {N_("off"), N_("on")};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return string indicating on/off depending on integer value
 *
 * \param  i  input integer
 *
 * \return  status string, possibly translated
 */
/*----------------------------------------------------------------------------*/

static const char *
_status(int i)
{
  return (i > 0) ? _(_astat[1]) : _(_astat[0]);
}

/*----------------------------------------------------------------*/
/*!
 *\brief Computes min/max for boundary statistics.
 *
 * Parameters:
 * \param[in]  s_id      stat id
 * \param[out] nbrfac    number of particles used for the statistics
 * \param[out] gmin      min value
 * \param[out] gmax      max value
 * \param[in]  unsnbr    inverse of the number of particles impacting
 *                       the boundary
 * \param[in]  unsnbrfou inverse of the number of particles impacting
 *                       the boundary taking fooling into account.
 */
 /*----------------------------------------------------------------*/

static void
_lagr_min_max_boundary_stats(int         s_id,
                             cs_lnum_t  *nbrfac,
                             cs_real_t  *gmin,
                             cs_real_t  *gmax,
                             cs_real_t   unsnbr[],
                             cs_real_t   unsnbrfou[])
{
  cs_lagr_boundary_interactions_t *lagr_bd_i = cs_glob_lagr_boundary_interactions;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* Initializations */
  *nbrfac = 0;
  *gmax = -cs_math_big_r;
  *gmin =  cs_math_big_r;

  cs_real_t threshold = cs_glob_lagr_stat_options->threshold;

  if (lagr_bd_i->imoybr[s_id] == 3) {

    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
      if (bound_stat[ifac + n_b_faces * lagr_bd_i->iencnb] > threshold) {
        *nbrfac++;
        *gmax = CS_MAX(*gmax, bound_stat[ifac + n_b_faces * s_id] * unsnbrfou[ifac]);
        *gmin = CS_MIN(*gmin, bound_stat[ifac + n_b_faces * s_id] * unsnbrfou[ifac]);
      }
    }

  }
  else if (lagr_bd_i->imoybr[s_id] == 2) {

    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
      if (bound_stat[ifac + n_b_faces * lagr_bd_i->inbr] > threshold) {
        *nbrfac++;
        *gmax = CS_MAX(*gmax, bound_stat[ifac + n_b_faces * s_id] * unsnbr[ifac]);
        *gmin = CS_MIN(*gmin, bound_stat[ifac + n_b_faces * s_id] * unsnbr[ifac]);
      }
    }

  }
  else if (lagr_bd_i->imoybr[s_id] == 1) {

    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
      if (bound_stat[ifac + n_b_faces * lagr_bd_i->inbr] > threshold) {
        *nbrfac++;
        *gmax = CS_MAX(*gmax, bound_stat[ifac + n_b_faces * s_id] / lagr_bd_i->tstatp);
        *gmin = CS_MIN(*gmin, bound_stat[ifac + n_b_faces * s_id] / lagr_bd_i->tstatp);
      }
    }

  }
  else if (lagr_bd_i->imoybr[s_id] == 0) {

    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
      if (bound_stat[ifac + n_b_faces * lagr_bd_i->inbr] > threshold) {
        *nbrfac++;
        *gmax = CS_MAX(*gmax, bound_stat[ifac + n_b_faces * s_id]);
        *gmin = CS_MIN(*gmin, bound_stat[ifac + n_b_faces * s_id]);
      }
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log Lagrangian module injection info.
 *
 * \param[in]  log  associated log file
 */
/*----------------------------------------------------------------------------*/

static void
_log_setup_injection(cs_log_t  log)
{
  /* Check if this call is needed */

  if (cs_glob_lagr_time_scheme == NULL)
    return;

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_OFF)
    return;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  cs_log_printf(log,
                _("\n"
                  "  Lagrangian particle injection\n"
                  "  -----------------------------\n\n"));

  for (int i_loc = 0; i_loc < 2; i_loc++) {

    cs_lagr_zone_data_t *zd = NULL;

    int n_zones = 0;

    if (i_loc == 0) {
      zd = cs_lagr_get_boundary_conditions();
      n_zones = cs_boundary_zone_n_zones();
    }
    else {
      zd = cs_lagr_get_volume_conditions();
      n_zones = cs_volume_zone_n_zones();
    }

    for (int z_id = 0; z_id < n_zones; z_id++) {

      const cs_zone_t  *z;
      if (i_loc == 0)
        z = cs_boundary_zone_by_id(z_id);
      else
        z = cs_volume_zone_by_id(z_id);

      for (int set_id = 0;
           set_id < zd->n_injection_sets[z_id];
           set_id++) {

        const cs_lagr_injection_set_t
          *zis = cs_lagr_get_injection_set(zd, z_id, set_id);

        cs_log_printf(log,
                      _("  zone: %d (%s), set:  %d\n"),
                      z->id, z->name, set_id);

        if (zis->n_inject > 0)
          cs_log_printf(log,
                        _("    n particles to inject: %llu\n"),
                        (unsigned long long)(zis->n_inject));

        if (zis->velocity_profile == -1)
          cs_log_printf(log,
                        _("    velocity from fluid\n"));
        else if (zis->velocity_profile == 0)
          cs_log_printf(log,
                        _("    velocity magnitude: %g (normal to boundary)\n"),
                        zis->velocity_magnitude);
        else if (zis->velocity_profile == 1)
          cs_log_printf(log,
                        _("    velocity: [%g, %g, %g]"),
                        zis->velocity[0], zis->velocity[1], zis->velocity[2]);

        cs_log_printf(log,
                      _("    diameter: %g; (variance: %g)\n"
                        "    density: %g\n"),
                      zis->diameter, zis->diameter_variance, zis->density);

        if (zis->flow_rate > 0)
          cs_log_printf(log,
                        _("    flow rate: %g\n"),
                        zis->flow_rate);

        cs_log_printf(log,
                      _("    statistical cluster id: %d\n"
                        "    statistical weight: %g\n"),
                      zis->cluster,
                      zis->stat_weight);

        if (cs_glob_lagr_model->deposition > 0)
          cs_log_printf(log,
                        _("    fouling index: %g\n"),
                        zis->fouling_index);

        if (   cs_glob_lagr_model->physical_model == 1
            && cs_glob_lagr_specific_physics->itpvar == 1) {
          if (zis->temperature_profile == 0)
            cs_log_printf(log,
                          _("    temperature from fluid\n"));
          else if (zis->temperature_profile == 1)
            cs_log_printf(log,
                          _("    temperature: %g\n"),
                          zis->temperature);
          cs_log_printf(log,
                        _("    Cp: %g\n"),
                        zis->cp);
          if (extra->radiative_model > 0)
            cs_log_printf(log,
                          _("    emissivity: %g\n"),
                          zis->emissivity);
        }
        else if (cs_glob_lagr_model->physical_model == 2) {
          cs_log_printf(log,
                        _("    coal number: %d\n"),
                        zis->coal_number);
        }

        cs_log_printf(log, "\n");
      }

    }

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log Lagrangian module output in the setup file.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_log_setup(void)
{
  /* Check if this call is needed */

  if (cs_glob_lagr_time_scheme == NULL)
    return;

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_OFF)
    return;

  /* Now add Lagrangian setup info */

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Lagrangian model options\n"
                  "------------------------\n"));

  cs_log_printf
    (CS_LOG_SETUP,
     _("  Continuous phase:\n"
       "    iilagr:                 %3d  (0: Lagrangian deactivated\n"
       "                                  1: one way coupling\n"
       "                                  2: two way coupling\n"
       "                                  3: on frozen fields)\n"
       "    restart: %s\n"
       "    statistics/return source terms restart: %s\n\n"
       "  Specific physics associated with particles\n"
       "    physical_model:         %3d  (0: no additional equations\n"
       "                                  1: equations on Dp Tp Mp\n"
       "                                  2: coal particles)\n"),
     cs_glob_lagr_time_scheme->iilagr,
     _status(cs_glob_lagr_time_scheme->isuila),
     _status(cs_glob_lagr_stat_options->isuist),
     cs_glob_lagr_model->physical_model);

  if (cs_glob_lagr_model->physical_model == 1)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    idpvar:                 %3d  (1: eqn diameter Dp,    or 0)\n"
         "    itpvar:                 %3d  (1: eqn temperature Tp, or 0)\n"
         "    impvar:                 %3d  (1: eqn mass Mp,        or 0)\n"),
       cs_glob_lagr_specific_physics->idpvar,
       cs_glob_lagr_specific_physics->itpvar,
       cs_glob_lagr_specific_physics->impvar);

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n  Global parameters:\n"
       "    user particle variables: %2d\n"
       "    isttio:                 %3d  (1: steady carrier phase)\n"),
     cs_glob_lagr_model->n_user_variables,
     cs_glob_lagr_time_scheme->isttio);

  if (cs_glob_lagr_model->physical_model == 2) {

    cs_log_printf
      (CS_LOG_SETUP,
       _("\n  Coal options:\n"
         "    fouling: %s\n"),
       _status(cs_glob_lagr_model->fouling));

    const cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

    for (int i = 0; i < extra->ncharb; i++)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    tprenc[%3d]:    %11.5e (threshold T for coal fouling %d)\n"),
         i, cs_glob_lagr_encrustation->tprenc[i], i);

    for (int i = 0; i < extra->ncharb; i++)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    visref[%3d]:    %11.5e (critical coal viscosity %d)\n"),
         i, cs_glob_lagr_encrustation->visref[i], i);

    for (int i = 0; i < extra->ncharb; i++)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    enc1[%3d]:      %11.5e (fouling coefficient 1 %d)\n"),
         i, cs_glob_lagr_encrustation->enc1[i], i);

    for (int i = 0; i < extra->ncharb; i++)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    enc2[%3d]:      %11.5e (fouling coefficient 2 %d)\n"),
         i, cs_glob_lagr_encrustation->enc2[i], i);

  }

  if (cs_glob_lagr_model->physical_model == 2) {

    cs_log_printf
      (CS_LOG_SETUP,
       _("\n  Return coupling options:\n"
         "    start iteration for time average:  %d\n"
         "    dynamic return coupling:           %s\n"
         "    mass return coupling:              %s\n"
         "    thermal return coupling:           %s\n"),
       cs_glob_lagr_source_terms->nstits,
       _status(cs_glob_lagr_source_terms->ltsdyn),
       _status(cs_glob_lagr_source_terms->ltsmas),
       _status(cs_glob_lagr_source_terms->ltsthe));
  }

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "  Statistics options:\n"
       "  starting iteration for statistics:        %d\n"
       "  starting iteration for steady statistics: %d\n"
       "  threshold for statistical meaning:        %11.3e\n"),
     cs_glob_lagr_stat_options->idstnt,
     cs_glob_lagr_stat_options->nstist,
     cs_glob_lagr_stat_options->threshold);

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n  Turbulent dispersion options:\n"
       "    lagrangian turbulent dispersion:              %s\n"
       "      identical to fluid turbulent diffusion:     %s\n"
       "    apply complete model from time step:          %d\n"),
     _status(cs_glob_lagr_time_scheme->idistu),
     _status(cs_glob_lagr_time_scheme->idiffl),
     cs_glob_lagr_time_scheme->modcpl);

  if (cs_glob_lagr_time_scheme->modcpl) {
    const char c_dir[] = "xyze";
    int _idirla = cs_glob_lagr_time_scheme->modcpl;
    cs_log_printf
      (CS_LOG_SETUP,
       _("    complete model main flow direction: %c\n"),
       c_dir[_idirla]);
  }

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n  Numerical options:\n"
       "    trajectory time scheme order:                 %d\n"
       "    Poisson correction for particle velocity:     %s\n"),
     cs_glob_lagr_time_scheme->t_order,
     _status(cs_glob_lagr_time_scheme->ilapoi));

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n  Trajectory/particle postprocessing options:\n"));
     for (int attr = 0; attr < CS_LAGR_N_ATTRIBUTES; attr++) {
       if (cs_lagr_post_get_attr(attr))
         cs_log_printf(CS_LOG_SETUP,
                       "    %s\n", cs_lagr_attribute_name[attr]);
     }

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n  Statistics for particles/boundary interaction:\n"));

  if (cs_glob_lagr_dim->n_boundary_stats == 0)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "none");
  if (cs_glob_lagr_boundary_interactions->has_part_impact_nbr)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "particle impact number");
  if (cs_glob_lagr_boundary_interactions->iflmbd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "particle mass flow");
  if (cs_glob_lagr_boundary_interactions->iangbd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "impact angle");
  if (cs_glob_lagr_boundary_interactions->ivitbd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "impact velocity");
  if (cs_glob_lagr_boundary_interactions->iencnbbd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "interactions with fouling");
  if (cs_glob_lagr_boundary_interactions->iencmabd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "fouling coal mass flux");
  if (cs_glob_lagr_boundary_interactions->iencdibd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "fouling coal diameter");
  if (cs_glob_lagr_boundary_interactions->iencckbd)
    cs_log_printf(CS_LOG_SETUP, "    %s\n", "fouling coal coke fraction");

  /* Volumic statistics   */

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Lagrangian statistics\n"
                  "---------------------\n\n"));

  cs_log_printf
    (CS_LOG_SETUP,
     _("  Start of calculation from absolute iteration number: %10d\n"),
     cs_glob_lagr_stat_options->idstnt);

  if (cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt) {

    if (cs_glob_lagr_time_scheme->isttio == 1) {
      cs_log_printf
        (CS_LOG_SETUP,
         _("  Start of steady-state statistics from Lagrangian "
           "iteration number: %10d\n"),
         cs_glob_lagr_stat_options->nstist);

    }
    cs_log_printf(CS_LOG_SETUP, "\n");

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log Lagrangian module output in the main log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_log_iteration(void)
{
  /* Check if this call is needed */

  if (cs_glob_lagr_time_scheme == NULL)
    return;

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_OFF)
    return;

  /* Return silently if called before Lagrangian start */
  if (cs_glob_lagr_particle_set == NULL)
    return;

  const cs_real_t  *b_stats = bound_stat;

  cs_log_printf(CS_LOG_DEFAULT,
                _("   ** INFORMATION ON THE LAGRANGIAN CALCULATION\n"));
  cs_log_separator(CS_LOG_DEFAULT);

  /* Log injection setup info on the first iteration.
     TODO this should be in the setup definition, but would require defining
     injections earlier in the setup phase, including with user functions */

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1)
    _log_setup_injection(CS_LOG_DEFAULT);

  /* Make sure counters are up-to-date */

  const cs_lagr_particle_counter_t *pc = cs_lagr_update_particle_counter();

  /* Number of particles  */
  cs_log_printf(CS_LOG_DEFAULT, "\n");
  cs_log_printf(CS_LOG_DEFAULT,
                _("   Current number of particles "
                  "(with and without statistical weight) :\n"));
  cs_log_printf(CS_LOG_DEFAULT, "\n");
  cs_log_printf
    (CS_LOG_DEFAULT,
     _("ln  newly injected                           %8llu   %14.5E\n"),
     (unsigned long long)(pc->n_g_new),
     pc->w_new);

  if (cs_glob_lagr_model->physical_model == 2 && cs_glob_lagr_model->fouling == 1)
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("ln  coal particles fouled                    %8llu   %14.5E\n"),
       (unsigned long long)(pc->n_g_fouling), pc->w_fouling);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("ln  out, or deposited and eliminated         %8llu   %14.5E\n"),
     (unsigned long long)(pc->n_g_exit), pc->w_exit);
  cs_log_printf
    (CS_LOG_DEFAULT,
     _("ln  deposited                                %8llu   %14.5E\n"),
     (unsigned long long)(pc->n_g_deposited), pc->w_deposited);

  if (cs_glob_lagr_model->resuspension > 0)
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("ln  resuspended                              %8llu   %14.5E\n"),
       (unsigned long long)(pc->n_g_resuspended),
       pc->w_resuspended);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("ln  lost in the location stage               %8llu\n"),
     (unsigned long long)(pc->n_g_failed));
  cs_log_printf
    (CS_LOG_DEFAULT,
     _("ln  total number at the end of the time step %8llu   %14.5E\n"),
     (unsigned long long)pc->n_g_total, pc->w_total);

  if (pc->n_g_cumulative_total > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  _("%% of lost particles (restart(s) included): %13.4E\n"),
                  pc->n_g_cumulative_failed * 100. / pc->n_g_cumulative_total);
      cs_log_separator(CS_LOG_DEFAULT);

  /* Flow rate for each zone   */
  cs_log_printf(CS_LOG_DEFAULT,
                _("   Zone  Class  Mass flow rate(kg/s)      Name (type)\n"));

  cs_lagr_zone_data_t *bdy_cond = cs_lagr_get_boundary_conditions();

  /* TODO log for volume conditions and internal zone */
#if 0
  cs_lagr_internal_condition_t *internal_cond = cs_lagr_get_internal_conditions();
#endif

  int n_stats = cs_glob_lagr_model->n_stat_classes + 1;

  cs_real_t *flow_rate;
  int flow_rate_size = bdy_cond->n_zones*n_stats;
  BFT_MALLOC(flow_rate, flow_rate_size, cs_real_t);

  for (int i = 0; i < flow_rate_size; i++)
    flow_rate[i] = bdy_cond->particle_flow_rate[i];

  cs_parall_sum(flow_rate_size, CS_REAL_TYPE, flow_rate);

  for (cs_lnum_t z_id = 0; z_id < bdy_cond->n_zones; z_id++) {

    if (CS_ABS(flow_rate[z_id*n_stats]) > 0.) {

      const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);
      const char *chcond;

      if (bdy_cond->zone_type[z_id] == CS_LAGR_INLET)
        chcond = _("inlet");

      else if (bdy_cond->zone_type[z_id] == CS_LAGR_REBOUND)
        chcond = _("rebound");

      else if (bdy_cond->zone_type[z_id] == CS_LAGR_OUTLET)
        chcond = _("outlet");

      else if (   bdy_cond->zone_type[z_id] == CS_LAGR_DEPO1
               || bdy_cond->zone_type[z_id] == CS_LAGR_DEPO2)
        chcond = _("deposition");

      else if (bdy_cond->zone_type[z_id] == CS_LAGR_FOULING)
        chcond = _("fouling");

      else if (bdy_cond->zone_type[z_id] == CS_LAGR_DEPO_DLVO)
        chcond = _("dlvo conditions");

      else if (bdy_cond->zone_type[z_id] == CS_LAGR_SYM)
        chcond = _("symmetry");

      else
        chcond = _("user");

      cs_log_printf(CS_LOG_DEFAULT,
                    "   %3d         %12.5e               %s (%s)\n",
                    z_id,
                    flow_rate[z_id*n_stats]/cs_glob_lagr_time_step->dtp,
                    z->name, chcond);

      for (int j = 1; j < n_stats; j++) {
        if (CS_ABS(flow_rate[z_id*n_stats + j]) > 0)
          cs_log_printf(CS_LOG_DEFAULT,
                        "         %3d   %12.5e\n",
                        j,
                        flow_rate[z_id*n_stats + j]/cs_glob_lagr_time_step->dtp);
      }

    }

  }

  cs_log_separator(CS_LOG_DEFAULT);

  BFT_FREE(flow_rate);

  /* Boundary statistics  */

  if (cs_glob_lagr_dim->n_boundary_stats > 0) {

    cs_log_printf(CS_LOG_DEFAULT,
                  _("   Boundary statistics :\n"));
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n");

    if (cs_glob_lagr_time_scheme->isttio == 1) {

      if (cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->nstist)
        cs_log_printf(CS_LOG_DEFAULT,
                      _("Number of iterations in steady-state statistics: %10d\n"),
                      cs_glob_lagr_boundary_interactions->npstf);

      else
        cs_log_printf(CS_LOG_DEFAULT,
                      _("Start of steady-state statistics from time step n°: %8d\n"),
                      cs_glob_lagr_stat_options->nstist);

    }

    cs_log_printf(CS_LOG_DEFAULT,
                  _("Total number of iterations in the statistics:%10d\n"),
                  cs_glob_lagr_boundary_interactions->npstft);
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n");

    cs_log_printf(CS_LOG_DEFAULT,
                  _("                           Min value    Max value    \n"));

    const cs_real_t _threshold = 1.e-30;

    cs_real_t *tabvr = NULL;
    cs_real_t *tabvrfou = NULL;

    if (cs_glob_lagr_boundary_interactions->has_part_impact_nbr == 1) {

      /* Allocate a work array */
      BFT_MALLOC(tabvr, cs_glob_mesh->n_b_faces, cs_real_t);

      int s_id = cs_glob_lagr_boundary_interactions->inbr;

      for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {

        if (b_stats[ifac + cs_glob_mesh->n_b_faces * s_id] > _threshold)
          tabvr[ifac] = 1.0 / b_stats[ifac + cs_glob_mesh->n_b_faces * s_id];
        else
          tabvr[ifac] = 0.0;

      }

    }

    if (cs_glob_lagr_boundary_interactions->iencnbbd == 1) {

      BFT_MALLOC(tabvrfou, cs_glob_mesh->n_b_faces, cs_real_t);

      int s_id = cs_glob_lagr_boundary_interactions->iencnb;

      for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {

        if (b_stats[ifac + cs_glob_mesh->n_b_faces * s_id] > _threshold)
          tabvrfou[ifac] = 1.0 / b_stats[ifac + cs_glob_mesh->n_b_faces * s_id];
        else
          tabvrfou[ifac] = 0.0;

      }
    }

    for (int s_id = 0; s_id < cs_glob_lagr_dim->n_boundary_stats; s_id++) {

      cs_real_t gmin;
      cs_real_t gmax;
      cs_lnum_t nbrfac;

      _lagr_min_max_boundary_stats(s_id, &nbrfac, &gmin, &gmax, tabvr, tabvrfou);

      cs_parall_max(1, CS_INT_TYPE, &nbrfac);
      cs_parall_min(1, CS_REAL_TYPE, &gmin);
      cs_parall_max(1, CS_REAL_TYPE, &gmax);

      /* If there is no particles, no statistics */
      if (nbrfac > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      "lp  %20s  %12.5E  %12.5E\n",
                      cs_glob_lagr_boundary_interactions->nombrd[s_id],
                      gmin,
                      gmax);
      else
        cs_log_printf(CS_LOG_DEFAULT,
                      "lp  %20s\n",
                      cs_glob_lagr_boundary_interactions->nombrd[s_id]);

    }

    /* Free memory */
    if (tabvr != NULL)
      BFT_FREE(tabvr);
    if (tabvrfou != NULL)
      BFT_FREE(tabvrfou);

    cs_log_separator(CS_LOG_DEFAULT);

  }

  /* Information about two-way coupling  */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

    if (cs_glob_lagr_time_scheme->isttio == 0) {

      cs_log_printf(CS_LOG_DEFAULT,
                    _("   Unsteady two-way coupling source terms:\n"));
      cs_log_separator(CS_LOG_DEFAULT);

    }
    else if (cs_glob_lagr_time_scheme->isttio == 1) {

      cs_log_printf(CS_LOG_DEFAULT,
                    _("   Two-way coupling source terms:\n"));
      cs_log_separator(CS_LOG_DEFAULT);

      if (cs_glob_time_step->nt_cur < cs_glob_lagr_source_terms->nstits)
        cs_log_printf(CS_LOG_DEFAULT,
                      _("Reset of the source terms (Start of steady-state at:): %10d\n"),
                      cs_glob_lagr_source_terms->nstits);

      else if (cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->nstist)
        cs_log_printf(CS_LOG_DEFAULT,
                      _("Number of iterations for the steady-state source terms:%10d\n"),
                      cs_glob_lagr_source_terms->npts);

    }

    cs_real_t vmax[2] = {cs_glob_lagr_source_terms->vmax,
                         cs_glob_lagr_source_terms->tmamax};
    cs_gnum_t cpt = cs_glob_lagr_source_terms->ntxerr;

    cs_parall_max(2, CS_REAL_TYPE, vmax);
    cs_parall_counter(&cpt, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("Maximum particle volume fraction: %14.5e\n"), vmax[0]);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("Maximum particle mass fraction:   %14.5e\n"), vmax[1]);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("Number of cells with a particle volume fraction "
                    "greater than 0.8: %10llu\n"), (unsigned long long)cpt);
    cs_log_separator(CS_LOG_DEFAULT);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
