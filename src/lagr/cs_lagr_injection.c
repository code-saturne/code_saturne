/*============================================================================
 * Particle injection for lagrangian module.
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

#include "cs_base.h"
#include "cs_math.h"

#include "cs_boundary_zone.h"
#include "cs_volume_zone.h"

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_thermal_model.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_random.h"

#include "cs_lagr.h"
#include "cs_lagr_gradients.h"
#include "cs_lagr_tracking.h"
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
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given local id in a given array of
 *        ordered values.
 *
 * We assume the id is present in the array.
 *
 * \param[in]  n   number of values
 * \param[in]  x   value to locate
 * \param[in]  a   array of ordered values (size n)
 *
 * \return  index of x in array (smallest i such that a[i] >= x)
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_segment_binary_search(cs_lnum_t     n,
                       double        x,
                       const double  a[])
{
  cs_lnum_t start_id = 0;
  cs_lnum_t end_id = n-1;
  cs_lnum_t mid_id = (end_id - start_id) / 2;

  x = CS_MIN(x, a[end_id]); /* precaution: force in range */

  while (start_id < end_id) {
    if (a[mid_id] < x)
      start_id = mid_id + 1;
    else /* if (a[mid_id] >= x) */
      end_id = mid_id;
    mid_id = start_id + ((end_id - start_id) / 2);
  }

  assert(mid_id >= 0 && mid_id < n);

  return mid_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute new particles in a given region.
 *
 * \param[in]   n_g_particles     global number of particles to inject
 * \param[in]   n_elts            number of elements in region
 * \param[in]   elt_id            element ids (or NULL)
 * \param[in]   elt_weight        parent element weights
 *                                (i.e. all local surfaces or volumes)
 * \param[in]   elt_profile       optional profile values for elements (or NULL)
 * \param[out]  elt_particle_idx  start index of added particles for each
 *                                element (size: n_elts + 1)
 *
 * \return  number of particles added on local rank
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_distribute_particles(cs_gnum_t         n_g_particles,
                      cs_lnum_t         n_elts,
                      const cs_lnum_t   elt_id[],
                      const cs_real_t   elt_weight[],
                      const cs_real_t  *elt_profile,
                      cs_lnum_t         elt_particle_idx[])
{
  cs_lnum_t n_particles = (cs_glob_n_ranks > 1) ? 0 : n_g_particles;

  /* Compute local element weight */

  cs_real_t *elt_cm_weight = NULL;

  BFT_MALLOC(elt_cm_weight, n_elts, cs_real_t);

  if (elt_id != NULL) {
    if (elt_profile != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++)
        elt_cm_weight[i] = elt_weight[elt_id[i]]*elt_profile[i];
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++)
        elt_cm_weight[i] = elt_weight[elt_id[i]];
    }
  }
  else {
    if (elt_profile != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++)
        elt_cm_weight[i] = elt_weight[i]*elt_profile[i];
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++)
        elt_cm_weight[i] = elt_weight[i];
    }
  }

  /* Transform to cumulative weight using Kahan summation */

  double l_weight = 0;
  {
    double d = 0., c = 0.;
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      double z = elt_cm_weight[i] - c;
      double t = d + z;
      c = (t - d) - z;
      d = t;
      elt_cm_weight[i] = d;
    }
    l_weight = d;
  }

#if defined(HAVE_MPI)

  /* Pre_distribution to various ranks; we assume that the number of
     injected particles at a given time is not huge, so it is cheaper
     to precompute the distribution on a single rank and broadcast it.
     For a higher number of particles, computing by blocs and then
     redistributing (with "all to all" operations) could be more efficient. */

  if (cs_glob_n_ranks > 1) {

    int n_ranks = cs_glob_n_ranks;
    int l_rank = cs_glob_rank_id;
    int r_rank = 0; /* Root rank for serialized operations */

    cs_lnum_t  *n_rank_particles = NULL;
    double     *cm_weight = NULL;

    if (l_rank == r_rank) {

      BFT_MALLOC(n_rank_particles, n_ranks, cs_lnum_t);
      BFT_MALLOC(cm_weight, n_ranks, double);

      for (int i = 0; i < n_ranks; i++)
        n_rank_particles[i] = 0;

    }

    MPI_Gather(&l_weight, 1, MPI_DOUBLE, cm_weight, 1, MPI_DOUBLE,
               r_rank, cs_glob_mpi_comm);

    if (l_rank == r_rank) {

      /* Scan (cumulative sum) operation */
      for (int i = 1; i < n_ranks; i++)
        cm_weight[i] += cm_weight[i-1];

      /* Scale to [0, 1] */
      double tot_weight = cm_weight[n_ranks-1];

      if (tot_weight > 0.) {

        for (int i = 0; i < n_ranks; i++)
          cm_weight[i] /= tot_weight;

        /* Compute distribution */

        for (cs_gnum_t i = 0; i < n_g_particles; i++) {
          cs_real_t r;
          cs_random_uniform(1, &r);
          int r_id = _segment_binary_search(n_ranks, r, cm_weight);
          n_rank_particles[r_id] += 1;
        }

      }

      BFT_FREE(cm_weight);
    }

    MPI_Scatter(n_rank_particles, 1, CS_MPI_LNUM,
                &n_particles, 1, CS_MPI_LNUM,
                r_rank, cs_glob_mpi_comm);

    BFT_FREE(n_rank_particles);
  }

#endif /* defined(HAVE_MPI) */

  /* Check for empty zones */

  if (n_particles > 0 && n_elts < 1)
    n_particles = 0;

  /* Now distribute locally */

  for (cs_lnum_t i = 0; i < n_elts; i++)
    elt_particle_idx[i] = 0;
  elt_particle_idx[n_elts] = 0;

  for (cs_lnum_t i = 0; i < n_elts; i++)
    elt_cm_weight[i] /= l_weight;

  /* Compute distribution */

  for (cs_lnum_t i = 0; i < n_particles; i++) {
    cs_real_t r;
    cs_random_uniform(1, &r);
    cs_lnum_t e_id = _segment_binary_search(n_elts, r, elt_cm_weight);
    elt_particle_idx[e_id+1] += 1;
  }

  BFT_FREE(elt_cm_weight);

  /* transform count to index */

  for (cs_lnum_t i = 0; i < n_elts; i++)
    elt_particle_idx[i+1] += elt_particle_idx[i];

  assert(elt_particle_idx[n_elts] == n_particles);

  return n_particles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check injection parameters are valid.
 *
 * \param[in]  zis  pointer to injection data for a given zone and set
 */
/*----------------------------------------------------------------------------*/

static void
_injection_check(const cs_lagr_injection_set_t  *zis)
{
  const char _profile_err_fmt_i[]
    = N_("Lagrangian %s zone %d, set %d\n"
         "  %s profile value (%d) is invalid.");
  const char _profile_err_fmt_d[]
    = N_("Lagrangian %s zone %d, set %d\n"
         "  %s profile value (%g) is invalid.");

  char z_type_name[32] = "unknown";
  if (zis->location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
    strncpy(z_type_name, _("boundary"), 31);
  else if (zis->location_id == CS_MESH_LOCATION_CELLS)
    strncpy(z_type_name, _("volume"), 31);
  z_type_name[31] = '\0';

  int z_id = zis->zone_id;
  int set_id = zis->set_id;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Verification of particle classes */

  if (cs_glob_lagr_model->n_stat_classes > 0) {
    if (   zis->cluster < 0
        || zis->cluster > cs_glob_lagr_model->n_stat_classes)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian module: \n"
                  "  number of clusters = %d is either not defined (negative)\n"
                  "  or > to the number of statistical classes %d\n"
                  "  for zone %d and set %d."),
                (int)zis->cluster,
                (int)cs_glob_lagr_model->n_stat_classes,
                z_id,
                set_id);
  }

  if (cs_glob_lagr_model->n_particle_aggregates > 0) {
    if (zis->particle_aggregate < 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian module: \n"
                  "  number of particles = %d is either not defined (negative)\n"
                  "  or smaller than 1 than \n"
                  "  for zone %d and set %d."),
                (int)zis->particle_aggregate,
                z_id,
                set_id);
  }

  /* temperature */
  if (   cs_glob_lagr_model->physical_model == 1
      && (   cs_glob_lagr_specific_physics->itpvar == 1
          || cs_glob_lagr_specific_physics->idpvar == 1
          || cs_glob_lagr_specific_physics->impvar == 1)) {
    if (zis->temperature_profile < 1 || zis->temperature_profile > 1)
      bft_error(__FILE__, __LINE__, 0, _profile_err_fmt_i,
                z_type_name, z_id, set_id,
                _("temperature"), (int)zis->temperature_profile);
  }

  /* velocity */
  if (   zis->location_id != CS_MESH_LOCATION_BOUNDARY_FACES
      && zis->velocity_profile == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Lagrangian %s zone %d, set %d:\n"
                " velocity profile type 0 may not be used\n"
                " for volume zones, as it requires surface normals."),
              z_type_name, z_id, set_id);
  else if (zis->velocity_profile <  -1 || zis->velocity_profile > 1)
    bft_error(__FILE__, __LINE__, 0, _profile_err_fmt_i,
              z_type_name, z_id, set_id,
              _("velocity"), (int)zis->velocity_profile);

  /* statistical weight */
  if (zis->stat_weight <= 0.0 && zis->flow_rate <= 0.0)
    bft_error(__FILE__, __LINE__, 0, _profile_err_fmt_d,
              z_type_name, z_id, set_id,
              _("statistical weight"), (double)zis->stat_weight);

  /* mass flow rate */
  if (zis->flow_rate > 0.0 && zis->n_inject  == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Lagrangian %s zone %d, set %d:\n"
                " flow rate is positive (%g)\n"
                " while number injected particles is 0."),
              z_type_name, z_id, set_id,
              (double)zis->flow_rate);

  /* particle properties: diameter, variance, and rho */
  if (cs_glob_lagr_model->physical_model != 2) {
    if (   zis->density  < 0.0
        || zis->diameter < 0.0
        || zis->diameter_variance < 0.0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian %s zone %d, set %d:\n"
           "  error on particle properties definition:\n"
           "  rho = %g, diameter = %g,\n"
           "  diameter standard deviation = %g\n"
           "This may lead to injection of  particles with negative diameters."),
         z_type_name, z_id, set_id,
         (double)zis->density,
         (double)zis->diameter,
         (double)zis->diameter_variance);
  }

  if (zis->diameter < 3.0 * zis->diameter_variance)
    bft_error(__FILE__, __LINE__, 0,
              _("Lagrangian %s zone %d, set %d:\n"
                "  diameter (%g) is smaller than 3 times\n"
                "  its standard deviation (%g)."),
              z_type_name, z_id, set_id,
              (double)zis->diameter,
              (double)zis->diameter_variance);

  /* Ellipsoidal particle properties: radii */
  if (cs_glob_lagr_model->shape == 2) {
    if (   zis->radii[0] < 0.0
        || zis->radii[1] < 0.0
        || zis->radii[2] < 0.0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian %s zone %d, set %d:\n"
           "  error on particle properties definition:\n"
           "  Ellispoid radii = %g, %g, %g\n"
           "This may lead to injection of  particles with negative radii."),
         z_type_name, z_id, set_id,
         (double)zis->radii[0],
         (double)zis->radii[1],
         (double)zis->radii[2]);
  }


  /* temperature and Cp */
  if (   cs_glob_lagr_model->physical_model == 1
      && cs_glob_lagr_specific_physics->itpvar == 1) {
    cs_real_t tkelvn = -cs_physical_constants_celsius_to_kelvin;
    if (zis->cp < 0.0 || zis->temperature < tkelvn)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  "  specific heat capacity (%g) is negative\n"
                  "  or temperature (%g) is lower than %g."),
                z_type_name, z_id, set_id,
                (double)zis->cp,
                (double)zis->temperature,
                (double)tkelvn);
  }

  /* emissivity */
  if (   cs_glob_lagr_model->physical_model == 1
      && cs_glob_lagr_specific_physics->itpvar == 1
      && extra->radiative_model > 0) {

    if (zis->emissivity < 0.0 || zis->emissivity > 1.0)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  "  particle emissivity (%g) is not properly set."),
                z_type_name, z_id, set_id,
                (double)zis->emissivity);

  }

  /* Coal */

  if (cs_glob_lagr_model->physical_model == 2) {

    cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

    if (zis->coal_number < 1 && zis->coal_number > extra->ncharb)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian %s zone %d, set %d:\n"
           "  the coal number %d for the injected particle is either negative\n"
           "  or greater than the maximum number of coals defined (%d)."),
         z_type_name, z_id, set_id,
         (int)zis->coal_number, (int)extra->ncharb);

    int coal_id = zis->coal_number - 1;

    /* properties of coal particles */
    if (zis->temperature < tkelvi)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  "  temperature is not properly set: %g."),
                z_type_name, z_id, set_id,
                (double)zis->temperature);

    /* Properties of coal particles */

    /* Composition of coal defined in XML file (DP_FCP) */

    cs_real_t *xashch = cs_glob_lagr_coal_comb->xashch;
    cs_real_t *cp2ch  = cs_glob_lagr_coal_comb->cp2ch;
    cs_real_t *xwatch = cs_glob_lagr_coal_comb->xwatch;
    cs_real_t *rho0ch = cs_glob_lagr_coal_comb->rho0ch;

    if (   rho0ch[coal_id] < 0.0
        || cp2ch[coal_id]  < 0.0
        || xwatch[coal_id] < 0.0
        || xwatch[coal_id] > 1.0
        || xashch[coal_id] < 0.0
        || xashch[coal_id] > 1.0)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  "  wrong conditions for coal number %d.\n"
                  "    coal density = %g\n"
                  "    Cp CP2CH = %g\n"
                  "    water mass fraction = %g\n"
                  "    ashes mass fraction = %g."),
                z_type_name, z_id, set_id,
                (int)coal_id,
                (double)rho0ch[coal_id],
                (double)cp2ch[coal_id],
                (double)xwatch[coal_id],
                (double)xashch[coal_id]);

    if (xwatch[coal_id] + xashch[coal_id] > 1.0)
      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  "  wrong conditions for coal number %d.\n"
                  "    water mass fraction = %g\n"
                  "    ashes mass fraction = %g\n"
                  "    mass fraction is larger than 1: %g."),
                z_type_name, z_id, set_id,
                (int)zis->coal_number,
                (double)xwatch[coal_id],
                (double)xashch[coal_id],
                (double)(xwatch[coal_id] + xashch[coal_id]));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build particle injection face ids array for a given boundary
 *        zone and set.
 *
 * The caller is responsible for freeing the returned array.
 *
 * \param[in]  n_faces            number of elements in zone
 * \param[in]  face_ids           matching face ids
 * \param[in]  face_particle_idx  starting id of new particles for a given
 *                                face (size: n_faces+1)
 *
 * \return array of ids of faces for injected particles
 *         (size: face_particle_idx[n_faces])
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t *
_get_particle_face_ids(cs_lnum_t         n_faces,
                       const cs_lnum_t   face_ids[],
                       const cs_lnum_t   face_particle_idx[])
{
  cs_lnum_t  *particle_face_id = NULL;

  cs_lnum_t n_p_new = face_particle_idx[n_faces];

  BFT_MALLOC(particle_face_id, n_p_new, cs_lnum_t);

  /* Loop on zone elements where particles are injected */

  n_p_new = 0;

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    /* Loop on particles added for this face */

    for (cs_lnum_t j = face_particle_idx[i]; j < face_particle_idx[i+1]; j++)
      particle_face_id[j] = face_ids[i];

  }

  return(particle_face_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize particle values
 *
 * \param[in,out]  p_set             particle set
 * \param[in]      zis               injection data this zone and set
 * \param[in]      time_id           time step indicator for fields
 *                                     0: use fields at current time step
 *                                     1: use fields at previous time step
 * \param[in]      n_elts            number of elements in zone
 * \param[in]      face_ids          matching face ids if zone is a boundary
 * \param[in]      elt_particle_idx  starting id of new particles for a given
 *                                   element (size: n_elts+1)
 */
/*----------------------------------------------------------------------------*/

static void
_init_particles(cs_lagr_particle_set_t         *p_set,
                const cs_lagr_injection_set_t  *zis,
                int                             time_id,
                cs_lnum_t                       n_elts,
                const cs_lnum_t                *face_ids,
                const cs_lnum_t                 elt_particle_idx[])
{
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Non-lagrangian fields */

  const cs_real_t  *xashch = cs_glob_lagr_coal_comb->xashch;
  const cs_real_t  *cp2ch  = cs_glob_lagr_coal_comb->cp2ch;
  const cs_real_t  *xwatch = cs_glob_lagr_coal_comb->xwatch;
  const cs_real_t  *rho0ch = cs_glob_lagr_coal_comb->rho0ch;

  const cs_real_t *vela = extra->vel->vals[time_id];
  const cs_real_t *cval_h = NULL, *cval_t = NULL;
  cs_real_t tscl_shift = 0;

  /* Initialize pointers (used to simplify future tests) */

  if (   (   cs_glob_lagr_model->physical_model == 1
          && cs_glob_lagr_specific_physics->itpvar == 1)
      || cs_glob_lagr_model->physical_model == 2) {

    if (   cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0
        || cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 0
        || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
      cval_t = cs_field_by_name("t_gas")->val;

    else {
      const cs_field_t *f = cs_field_by_name_try("temperature");
      if (f != NULL)
        cval_t = f->val;
      else if (   cs_glob_thermal_model->itherm
               == CS_THERMAL_MODEL_ENTHALPY)
        cval_h = cs_field_by_name("enthalpy")->val;
    }

    if (cs_glob_thermal_model->itpscl == 1) /* Kelvin */
      tscl_shift = - cs_physical_constants_celsius_to_kelvin;
  }

  const cs_real_t pis6 = cs_math_pi / 6.0;

  /* Loop on zone elements where particles are injected */

  for (cs_lnum_t li = 0; li < n_elts; li++) {

    cs_lnum_t n_e_p = elt_particle_idx[li+1] - elt_particle_idx[li];

    if (n_e_p < 1)
      continue;

    cs_lnum_t p_s_id = p_set->n_particles +  elt_particle_idx[li];
    cs_lnum_t p_e_id = p_s_id + n_e_p;

    const cs_lnum_t face_id = (face_ids != NULL) ? face_ids[li] : -1;

    /* Loop on particles added for this face */

    for (cs_lnum_t p_id = p_s_id; p_id < p_e_id; p_id++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

      cs_lnum_t cell_id = cs_lagr_particles_get_lnum(p_set, p_id,
                                                     CS_LAGR_CELL_ID);

      /* Random value associated with each particle */

      cs_real_t part_random = -1;
      cs_random_uniform(1, &part_random);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RANDOM_VALUE,
                                part_random);

      /* Particle velocity components */

      cs_real_t *part_vel = cs_lagr_particle_attr(particle, p_am,
                                                  CS_LAGR_VELOCITY);

      /* prescribed components */
      if (zis->velocity_profile == 1) {
        for (cs_lnum_t i = 0; i < 3; i++)
          part_vel[i] = zis->velocity[i];
      }

      /* prescribed norm */
      else if (zis->velocity_profile == 0) {
        assert(face_id >= 0);
        for (cs_lnum_t i = 0; i < 3; i++)
          part_vel[i] = -   fvq->b_face_normal[face_id * 3 + i]
                          / fvq->b_face_surf[face_id]
                          * zis->velocity_magnitude;
      }

      /* velocity as seen from fluid */
      else if (zis->velocity_profile ==  -1) {
        for (cs_lnum_t i = 0; i < 3; i++)
          part_vel[i] = vela[cell_id * 3  + i];
      }

      /* fluid velocity seen */
      cs_real_t *part_seen_vel = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_VELOCITY_SEEN);
      for (cs_lnum_t i = 0; i < 3; i++)
        part_seen_vel[i] = vela[cell_id * 3 + i];

      /* Residence time (may be negative to ensure continuous injection) */
      if (zis->injection_frequency == 1) {
        cs_real_t res_time = - part_random *cs_glob_lagr_time_step->dtp;
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RESIDENCE_TIME,
                                  res_time);
      }
      else
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RESIDENCE_TIME,
                                  0.0);

      /* Diameter (always set base) */

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                                zis->diameter);

      /* Shape for spheroids without inertia */
      if (cs_glob_lagr_model->shape == 1) {

        /* Spherical radii a b 0 */
        cs_real_t *radii = cs_lagr_particle_attr(particle, p_am,
                                                        CS_LAGR_RADII);
        /* Shape lambda, gamma , 0 , 0 */
        cs_real_t *shape_param = cs_lagr_particle_attr(particle, p_am,
                                                        CS_LAGR_SHAPE_PARAM);
        /* Orientation */
        cs_real_t *orientation = cs_lagr_particle_attr(particle, p_am,
                                                        CS_LAGR_ORIENTATION);
        for (cs_lnum_t i = 0; i < 3; i++) {
          radii[i] = zis->radii[i];
          orientation[i] = zis->orientation[i];
        }

        /* Compute shape parameters from radii */
        //FIXME do note divide by 0...
        //shape_param[0] = 1 - radii[2] / radii[1] ;  /* lambda */
        //shape_param[1] =  (cs_math_pow2(shape_param[0]) -1)/
        //                      (cs_math_pow2(shape_param[0]) +1)  ;  /* gamma */
        //Let start with the rod case
        shape_param[0] = 1000000. ; /* lambda is infinity */
        shape_param[1] = 1.; /* gamma is one */
        shape_param[2] = 0;
        shape_param[3] = 0;

        /* Compute orientation from uniform orientation on a unit-sphere */
        cs_real_t theta0 ;
        cs_real_t phi0 ;
        cs_random_uniform(1, &theta0) ;
        cs_random_uniform(1, &phi0) ;
        theta0   = acos(2.0*theta0-1.0) ;
        phi0     = phi0*2.0*cs_math_pi ;
        orientation[0] = sin(theta0)*cos(phi0) ;
        orientation[1] = sin(theta0)*sin(phi0) ;
        orientation[2] = cos(theta0) ;
        /* TODO initialize other things */

        }

      /* Shape general ellipsoid with inertia */
      if (cs_glob_lagr_model->shape == 2) {

        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_SHAPE,
                                  zis->shape);

        /* Ellipsoid radii a b c */
        cs_real_t *radii = cs_lagr_particle_attr(particle, p_am,
                                                     CS_LAGR_RADII);

        cs_real_t *shape_param = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_SHAPE_PARAM);
        /* Euler parameters */
        cs_real_t *euler = cs_lagr_particle_attr(particle, p_am,
                                                 CS_LAGR_EULER);

        for (cs_lnum_t i = 0; i < 3; i++)
          radii[i] = zis->radii[i];

        for (cs_lnum_t i = 0; i < 4; i++)
          euler[i] = zis->euler[i];


        /* TODO initialize other things */

        /* Compute shape parameters from radii */
        // FIXME is it valid for all ellispoids or for spheroids only?
        cs_real_t lamb = radii[2] / radii[1];//FIXME do note divide by 0...
        cs_real_t lamb_m1 = (radii[2] - radii[1]) / radii[1];
        cs_real_t lamb_p1 = (radii[2] + radii[1]) / radii[1];
        cs_real_t _a2 = radii[0] * radii[0];
        //TODO MF shape_param check development in series
        cs_real_t aux1 = lamb * lamb;
        cs_real_t aux2 = aux1 -1;
        if (lamb_m1 > 1e-10) {
          cs_real_t aux3 = sqrt(aux2 - 1);
          cs_real_t kappa = -log(lamb + aux3);
          shape_param[0] = aux1/aux2 + lamb*kappa/(aux2*aux3);
          shape_param[1] = shape_param[0];
          shape_param[2] = -2./aux2 - 2.*lamb*kappa/(aux2*aux3);
          shape_param[3] = -2. * _a2 *lamb*kappa/aux3;
        }
        else if (lamb_m1 < -1e-10) {
          cs_real_t aux3 = sqrt(1. - aux2);
          cs_real_t kappa = acos(lamb);
          shape_param[0] = aux1/aux2+lamb*kappa/(-aux2*aux3);
          shape_param[1] = shape_param[0];
          shape_param[2] = -2./aux2 - 2.*lamb*kappa/(-aux2*aux3);
          shape_param[3] = 2. * _a2 * lamb*kappa/aux3;
        }
        else {
          shape_param[0] = 2.0/3.0;
          shape_param[1] = 2.0/3.0;
          shape_param[2] = 2.0/3.0;
          shape_param[3] = 2. * _a2;
        }

        /* Compute Euler angles
           (random orientation with a uniform distribution in [-1;1]) */
        cs_real_33_t trans_m;
        // Generate the first two vectors
        for (cs_lnum_t id = 0; id < 3; id++) {
          cs_random_uniform(1, &trans_m[id][0]); /* (?,0) */
          cs_random_uniform(1, &trans_m[id][1]); /* (?,1) */
          cs_random_uniform(1, &trans_m[id][2]); /* (?,2) */
          cs_real_3_t loc_vector =  {-1.+2*trans_m[id][0],
            -1.+2*trans_m[id][1],
            -1.+2*trans_m[id][2]};
          cs_real_t norm_trans_m = cs_math_3_norm( loc_vector );
          while ( norm_trans_m > 1 )
          {
            cs_random_uniform(1, &trans_m[id][0]); /* (?,0) */
            cs_random_uniform(1, &trans_m[id][1]); /* (?,1) */
            cs_random_uniform(1, &trans_m[id][2]); /* (?,2) */
            loc_vector[0] = -1.+2*trans_m[id][0];
            loc_vector[1] = -1.+2*trans_m[id][1];
            loc_vector[2] = -1.+2*trans_m[id][2];
            norm_trans_m = cs_math_3_norm( loc_vector );
          }
          for (cs_lnum_t id1 = 0; id1 < 3; id1++)
            trans_m[id][id1] = (-1.+2*trans_m[id][id1]) / norm_trans_m;
        }
        // Correct 2nd vector (for perpendicularity to the 1st)
        cs_real_3_t loc_vector0 =  {trans_m[0][0],
          trans_m[0][1],
          trans_m[0][2]};
        cs_real_3_t loc_vector1 =  {trans_m[1][0],
          trans_m[1][1],
          trans_m[1][2]};
        cs_real_t scal_prod = cs_math_3_dot_product(loc_vector0, loc_vector1);
        for (cs_lnum_t id = 0; id < 3; id++)
          trans_m[1][id] -= scal_prod * trans_m[0][id];
        // Re-normalize
        loc_vector1[0] = trans_m[1][0];
        loc_vector1[1] = trans_m[1][1];
        loc_vector1[2] = trans_m[1][2];
        cs_real_t norm_trans_m = cs_math_3_norm( loc_vector1 );
        for (cs_lnum_t id = 0; id < 3; id++)
          trans_m[1][id] /= norm_trans_m;

        // Compute last vector (cross product of the two others)
        loc_vector1[0] = trans_m[1][0];
        loc_vector1[1] = trans_m[1][1];
        loc_vector1[2] = trans_m[1][2];
        cs_real_3_t loc_vector2 =  {trans_m[2][0],
          trans_m[2][1],
          trans_m[2][2]};
        cs_math_3_cross_product( loc_vector0, loc_vector1, loc_vector2);
        for (cs_lnum_t id = 0; id < 3; id++)
          trans_m[2][id] = loc_vector2[id];

        // Write Euler angles
        cs_real_t random;
        cs_random_uniform(1, &random);
        if (random >= 0.5)
          euler[0] = pow( 0.25*(trans_m[0][0]+trans_m[1][1]+trans_m[2][2]+1.) ,0.5);
        else
          euler[0] = -pow( 0.25*(trans_m[0][0]+trans_m[1][1]+trans_m[2][2]+1.) ,0.5);
        euler[1] = 0.25 * (trans_m[2][1] - trans_m[1][2]) / euler[0];
        euler[2] = 0.25 * (trans_m[0][2] - trans_m[2][0]) / euler[0];
        euler[3] = 0.25 * (trans_m[1][0] - trans_m[0][1]) / euler[0];

        /* Compute initial angular velocity */
        // Get velocity gradient
        cs_lagr_gradients(0, extra->grad_pr, extra->grad_vel);
        // Local reference frame
        cs_real_33_t grad_vf_r;
        cs_math_33_transform_a_to_r(extra->grad_vel[cell_id], trans_m, grad_vf_r);

        cs_real_t *ang_vel = cs_lagr_particle_attr(particle, p_am,
            CS_LAGR_ANGULAR_VEL);

        ang_vel[0] = 0.5*(grad_vf_r[2][1] - grad_vf_r[1][2]);
        ang_vel[1] = 0.5*(grad_vf_r[0][2] - grad_vf_r[2][0]);
        ang_vel[2] = 0.5*(grad_vf_r[0][1] - grad_vf_r[1][0]);

      }

      if (zis->diameter_variance > 0.0) {

        /* Randomize diameter, ensuring we obtain a
           positive diameter in the 99,7% range */

        cs_real_t d3   = 3.0 * zis->diameter_variance;

        int i_r = 0; /* avoid infinite loop in case of very improbable
                        random series... */

        for (i_r = 0; i_r < 20; i_r++) {
          double    random;
          cs_random_normal(1, &random);

          cs_real_t diam =   zis->diameter
                           + random * zis->diameter_variance;

          if (diam > 0 && (   diam >= zis->diameter - d3
                           && diam <= zis->diameter + d3)) {
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, diam);
            break;
          }
        }

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
                                  zis->cluster);

      if (cs_glob_lagr_model->n_particle_aggregates > 0)
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_PARTICLE_AGGREGATE,
                                  zis->particle_aggregate);

      /* used for 2nd order only */
      if (p_am->displ[0][CS_LAGR_TAUP_AUX] > 0)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TAUP_AUX, 0.0);

      if (   cs_glob_lagr_model->physical_model == 0
          || cs_glob_lagr_model->physical_model == 1) {

        if (cs_glob_lagr_model->clogging == 0)
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                    zis->density * pis6 * d3);
        else
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                    zis->density * pis6 * d3
                                    * pow(1.0-mporos, 3));

        if (   cs_glob_lagr_model->physical_model == 1
            && cs_glob_lagr_specific_physics->itpvar == 1) {

          if (cval_t != NULL)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_FLUID_TEMPERATURE,
                                      cval_t[cell_id] + tscl_shift);

          else if (cval_h != NULL) {

            int mode = 1;
            cs_real_t temp[1];
            CS_PROCF(usthht, USTHHT)(&mode, &(cval_h[cell_id]), temp);
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_FLUID_TEMPERATURE,
                                      temp[0]);

          }

          /* constant temperature set, may be modified later by user function */
          if (zis->temperature_profile == 1)
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE,
                                      zis->temperature);

          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP,
                                    zis->cp);
          if (extra->radiative_model > 0)
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_EMISSIVITY,
                                      zis->emissivity);

        }

      }

      else if (cs_glob_lagr_model->physical_model == 2) {

        int coal_id = zis->coal_number - 1;

        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_COAL_ID, coal_id);
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FLUID_TEMPERATURE,
                                  cval_t[cell_id] + tscl_shift);

        cs_real_t *particle_temp
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_TEMPERATURE);
        for (int ilayer = 0;
             ilayer < cs_glob_lagr_model->n_temperature_layers;
             ilayer++)
          particle_temp[ilayer] = zis->temperature;

        /* composition from DP_FCP */

        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP, cp2ch[coal_id]);

        cs_real_t mass = rho0ch[coal_id] * pis6 * d3;

        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mass);
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_WATER_MASS,
                                  xwatch[coal_id] * mass);

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
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COAL_DENSITY);
        for (int ilayer = 0;
             ilayer < cs_glob_lagr_model->n_temperature_layers;
             ilayer++)
          particle_coal_density[ilayer] = rho0ch[coal_id];

      }

      /* statistical weight */
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_STAT_WEIGHT,
                                zis->stat_weight);

      /* Fouling index */
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FOULING_INDEX,
                                zis->fouling_index);

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

      /* Initialize the additional user variables */

      for (int i = 0;
           i < cs_glob_lagr_model->n_user_variables;
           i++)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_USER + i, 0.0);

    }

  }

  /* Update weights to have the correct flow rate
     -------------------------------------------- */

  if (zis->flow_rate > 0.0 && zis->n_inject > 0) {

    cs_real_t dmass = 0.0;

    cs_lnum_t p_s_id = p_set->n_particles;
    cs_lnum_t p_e_id = p_s_id + elt_particle_idx[n_elts];

    for (cs_lnum_t p_id = p_s_id; p_id < p_e_id; p_id++)
      dmass += cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS);

    cs_parall_sum(1, CS_REAL_TYPE, &dmass);

    /* Compute weights */

    if (dmass > 0.0) {
      cs_real_t s_weight =   zis->flow_rate * cs_glob_lagr_time_step->dtp
                           / dmass;
      for (cs_lnum_t p_id = p_s_id; p_id < p_e_id; p_id++)
        cs_lagr_particles_set_real(p_set, p_id, CS_LAGR_STAT_WEIGHT, s_weight);
    }

    else {

      char z_type_name[32] = "unknown";
      if (zis->location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
        strncpy(z_type_name, _("boundary"), 31);
      else if (zis->location_id == CS_MESH_LOCATION_CELLS)
        strncpy(z_type_name, _("volume"), 31);
      z_type_name[31] = '\0';

      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  " imposed flow rate is %g\n"
                  " while mass of injected particles is 0."),
                z_type_name, zis->zone_id, zis->set_id,
                (double)zis->flow_rate);

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check injected particle values
 *
 * \param[in,out]  p_set             particle set
 * \param[in]      zis               injection data for this zone and set
 * \param[in]      n_elts            number of elements in zone
 * \param[in]      elt_particle_idx  starting id of new particles for a given
 *                                   element (size: n_elts+1)
 */
/*----------------------------------------------------------------------------*/

static void
_check_particles(cs_lagr_particle_set_t         *p_set,
                 const cs_lagr_injection_set_t  *zis,
                 cs_lnum_t                       n_elts,
                 const cs_lnum_t                 elt_particle_idx[])
{
  const cs_lnum_t s_id = p_set->n_particles;
  const cs_lnum_t e_id = s_id + elt_particle_idx[n_elts];

  char z_type_name[32] = "unknown";
  if (zis->location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
    strncpy(z_type_name, _("boundary"), 31);
  else if (zis->location_id == CS_MESH_LOCATION_CELLS)
    strncpy(z_type_name, _("volume"), 31);
  z_type_name[31] = '\0';

  int attrs[] = {CS_LAGR_DIAMETER, CS_LAGR_MASS, CS_LAGR_STAT_WEIGHT,
                 CS_LAGR_CP};

  for (cs_lnum_t p_id = s_id; p_id < e_id; p_id++) {

    for (int i_attr = 0; i_attr < 4; i_attr++) {

      int attr = attrs[i_attr];

      if (p_set->p_am->count[1][attr] > 0) {

        cs_real_t val  = cs_lagr_particles_get_real(p_set, p_id, attr);

        if (val <= 0.0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian %s zone %d, set %d:\n"
                      "  particle %d has a negative %s: %g"),
                    z_type_name, zis->zone_id, zis->set_id,
                    p_id, cs_lagr_attribute_name[attr], (double)val);

      }

    }

  }

  if (cs_glob_lagr_model->physical_model == 2) {

    int r01_attrs[] = {CS_LAGR_WATER_MASS, CS_LAGR_COAL_MASS, CS_LAGR_COKE_MASS,
                       CS_LAGR_COAL_DENSITY};
    int r00_attrs[] = {CS_LAGR_SHRINKING_DIAMETER, CS_LAGR_INITIAL_DIAMETER};

    for (cs_lnum_t p_id = s_id; p_id < e_id; p_id++) {

      for (int i_attr = 0; i_attr < 4; i_attr++) {

        int attr = r01_attrs[i_attr];
        int n_vals = p_set->p_am->count[1][attr];
        cs_real_t *vals = cs_lagr_particles_attr(p_set, p_id, attr);

        for (int l_id = 0; l_id < n_vals; l_id++) {
          if (vals[l_id] < 0.0) {
            if (n_vals == 1)
              bft_error(__FILE__, __LINE__, 0,
                        _("Lagrangian %s zone %d, set %d:\n"
                          "  particle %d has a negative %s: %g"),
                        z_type_name, zis->zone_id, zis->set_id,
                        p_id, cs_lagr_attribute_name[attr], (double)vals[0]);
            else
              bft_error(__FILE__, __LINE__, 0,
                        _("Lagrangian %s zone %d, set %d:\n"
                          "  particle %d has a negative %s\n"
                          "  in layer %d: %g"),
                        z_type_name, zis->zone_id, zis->set_id,
                        p_id, cs_lagr_attribute_name[attr], l_id, (double)vals[l_id]);
          }

        }

      }

      for (int i_attr = 0; i_attr < 2; i_attr++) {

        int attr = r00_attrs[i_attr];
        cs_real_t val = cs_lagr_particles_get_real(p_set, p_id, attr);

        if (val < 0) {
          bft_error(__FILE__, __LINE__, 0,
                    _("Lagrangian %s zone %d, set %d:\n"
                      "  particle %d has a negative %s: %g"),
                    z_type_name, zis->zone_id, zis->set_id,
                    p_id, cs_lagr_attribute_name[attr], (double)val);

        }

      }

    }

  }
}

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
  CS_UNUSED(itypfb);

  /* We may be mapped to an auxiliary field with no previous time id */

  cs_real_t dnbpnw_preci = 0.;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();
  time_id = CS_MIN(time_id, extra->vel->n_time_vals -1);

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;

  /* Non-lagrangian fields */
  cs_real_t *vela = extra->vel->vals[time_id];

  cs_lagr_particle_counter_t *pc = cs_lagr_get_particle_counter();
  const cs_time_step_t *ts = cs_glob_time_step;

  const int n_stats = cs_glob_lagr_model->n_stat_classes + 1;

  /* Initialization */

  cs_lagr_zone_data_t  *zda[2] = {cs_lagr_get_boundary_conditions(),
                                  cs_lagr_get_volume_conditions()};

  cs_lagr_get_internal_conditions();

  /* Boundary conditions */

  {
    cs_lagr_zone_data_t *zd = zda[0];

    for (int z_id = 0; z_id < zd->n_zones; z_id++) {

      if (zd->zone_type[z_id] > CS_LAGR_BC_USER)
        bft_error(__FILE__, __LINE__, 0,
                  _("Lagrangian boundary zone %d nature %d is unknown."),
                  z_id + 1,
                  (int)zd->zone_type[z_id]);

      if (   zd->zone_type[z_id] == CS_LAGR_FOULING
          && cs_glob_lagr_model->physical_model != 2)
        bft_error
          (__FILE__, __LINE__, 0,
           _("Lagrangian boundary zone %d nature is of type CS_LAGR_FOULING,\n"
             "but cs_glob_lagr_model->physical_model is not equal to 2."),
           z_id);
      if (   zd->zone_type[z_id] == CS_LAGR_FOULING
          && cs_glob_lagr_model->fouling != 1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("Lagrangian boundary zone %d nature is of type CS_LAGR_FOULING,\n"
             "but fouling is not activated."),
           z_id);

    }

  }

  /* Reset some particle counters */

  p_set->n_part_new = 0;
  p_set->weight_new = 0.0;

  for (int i_loc = 0; i_loc < 2; i_loc++) {
    cs_lagr_zone_data_t *zd = zda[i_loc];
    int fr_size = zd->n_zones * n_stats;
    for (int i = 0; i < fr_size; i++)
      zd->particle_flow_rate[i] = 0;
  }

  /* Injection due to precipitation/Dissolution
     ------------------------------------------ */

  if (cs_glob_lagr_model->precipitation == 1)
    cs_lagr_precipitation_injection(vela, &dnbpnw_preci);

  /* User-defined injection
     ---------------------- */

  /* Check various condition types and optional maximum particle limit */

  unsigned long long n_g_particles_next = pc->n_g_total;

  for (int i_loc = 0; i_loc < 2; i_loc++) {

    cs_lagr_zone_data_t *zd = zda[i_loc];

    /* compute global number of injected particles */

    for (int z_id = 0; z_id < zd->n_zones; z_id++) {
      for (int set_id = 0; set_id < zd->n_injection_sets[z_id]; set_id++) {
        cs_lagr_injection_set_t *zis
          = cs_lagr_get_injection_set(zd, z_id, set_id);
        _injection_check(zis);
        n_g_particles_next += (unsigned long long) (zis->n_inject);
      }
    }

  }

  /* Avoid injection if maximum defined number of particles reached */

  if (n_g_particles_next > cs_lagr_get_n_g_particles_max()) {

    bft_printf(_("\n Lagrangian module: \n"));
    bft_printf
      (_("  If particles are injected at time step %d,\n"
         "  the total number of particles in the domain would increase from\n"
         "  %llu to %llu, exceeding the maximums set by\n"
         "  cs_lagr_set_n_g_particles_max. (%llu).\n"
         "  No particles will be injected for this time step.\n"),
       ts->nt_cur,
       (unsigned long long)(pc->n_g_total),
       (unsigned long long)n_g_particles_next,
       (unsigned long long)(cs_lagr_get_n_g_particles_max()));

    return;

  }

  /* Now inject new particles
     ------------------------ */

  cs_lnum_t n_elts_m = CS_MAX(mesh->n_b_faces, mesh->n_cells);
  cs_lnum_t *elt_particle_idx = NULL;
  BFT_MALLOC(elt_particle_idx, n_elts_m+1, cs_lnum_t);

  /* Loop in injection type (boundary, volume) */

  for (int i_loc = 0; i_loc < 2; i_loc++) {

    cs_lagr_zone_data_t *zd = zda[i_loc];

    int n_zones = 0;

    const cs_real_t  *elt_weight = NULL;

    if (i_loc == 0) { /* boundary */
      elt_weight = fvq->b_face_surf;
      n_zones = cs_boundary_zone_n_zones();
    }
    else {            /* volume */
      elt_weight = fvq->cell_vol;
      n_zones = cs_volume_zone_n_zones();
    }

    /* Loop on injection zones */

    for (int z_id = 0; z_id < n_zones; z_id++) {

      /* Loop on injected sets */

      cs_lnum_t         n_z_elts = 0;
      const cs_lnum_t  *z_elt_ids = NULL;

      if (i_loc == 0) {
        const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);
        n_z_elts = z->n_elts;
        z_elt_ids = z->elt_ids;
      }
      else {
        const cs_zone_t  *z = cs_volume_zone_by_id(z_id);
        n_z_elts = z->n_elts;
        z_elt_ids = z->elt_ids;
      }

      for (int set_id = 0;
           set_id < zd->n_injection_sets[z_id];
           set_id++) {

        const cs_lagr_injection_set_t *zis = NULL;

        zis = cs_lagr_get_injection_set(zd, z_id, set_id);

        int injection_frequency = zis->injection_frequency;

        /* Inject only at first time step if injection frequency is zero */

        if (injection_frequency <= 0) {
          if (ts->nt_cur == ts->nt_prev+1 && pc->n_g_cumulative_total == 0)
            injection_frequency = ts->nt_cur;
          else
            injection_frequency = ts->nt_cur+1;
        }

        if (ts->nt_cur % injection_frequency != 0)
          continue;

        cs_real_t *elt_profile = NULL;
        if (zis->injection_profile_func != NULL) {
          BFT_MALLOC(elt_profile, n_z_elts, cs_real_t);
          zis->injection_profile_func(zis->zone_id,
                                      zis->location_id,
                                      zis->injection_profile_input,
                                      n_z_elts,
                                      z_elt_ids,
                                      elt_profile);
        }

        cs_lnum_t n_inject = _distribute_particles(zis->n_inject,
                                                   n_z_elts,
                                                   z_elt_ids,
                                                   elt_weight,
                                                   elt_profile,
                                                   elt_particle_idx);

        BFT_FREE(elt_profile);

        if (cs_lagr_particle_set_resize(p_set->n_particles + n_inject) < 0)
          bft_error(__FILE__, __LINE__, 0,
                    "Lagrangian module internal error: \n"
                    "  resizing of particle set impossible but previous\n"
                    "  size computation did not detect this issue.");

        /* Define particle coordinates and place on faces/cells */

        if (zis->location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
          cs_lagr_new(p_set,
                      n_z_elts,
                      z_elt_ids,
                      elt_particle_idx);
        else
          cs_lagr_new_v(p_set,
                        n_z_elts,
                        z_elt_ids,
                        elt_particle_idx);

        BFT_FREE(elt_profile);

        /* Initialize other particle attributes */

        _init_particles(p_set,
                        zis,
                        time_id,
                        n_z_elts,
                        z_elt_ids,
                        elt_particle_idx);

        assert(n_inject == elt_particle_idx[n_z_elts]);

        cs_lnum_t particle_range[2] = {p_set->n_particles,
                                       p_set->n_particles + n_inject};

        cs_lagr_new_particle_init(particle_range,
                                  time_id,
                                  vislen);

        /* Advanced user modification:

           WARNING: the user may change the particle coordinates but is
           prevented from changing the previous location (otherwise, if
           the particle is not in the same cell anymore, it would be lost).

           Moreover, a precaution has to be taken when calling
           "current to previous" in the tracking stage.
        */

        {
          cs_lnum_t *particle_face_ids = NULL;

          if (zis->location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
            particle_face_ids = _get_particle_face_ids(n_z_elts,
                                                       z_elt_ids,
                                                       elt_particle_idx);

          cs_lnum_t *saved_cell_id;
          cs_real_3_t *saved_coords;
          BFT_MALLOC(saved_cell_id, n_inject, cs_lnum_t);
          BFT_MALLOC(saved_coords, n_inject, cs_real_3_t);

          for (cs_lnum_t i = 0; i < n_inject; i++) {
            cs_lnum_t p_id = particle_range[0] + i;

            saved_cell_id[i] = cs_lagr_particles_get_lnum(p_set,
                                                           p_id,
                                                           CS_LAGR_CELL_ID);
            const cs_real_t *p_coords
              = cs_lagr_particles_attr_const(p_set,
                                             p_id,
                                             CS_LAGR_COORDS);
            for (cs_lnum_t j = 0; j < 3; j++)
              saved_coords[i][j] = p_coords[j];
          }

          cs_user_lagr_in(p_set,
                          zis,
                          particle_range,
                          particle_face_ids,
                          vislen);

          /* For safety, build values at previous time step, but reset saved values
             for previous cell number and particle coordinates */

          for (cs_lnum_t i = 0; i < n_inject; i++) {
            cs_lnum_t p_id = particle_range[0] + i;

            cs_lagr_particles_current_to_previous(p_set, p_id);

            cs_lagr_particles_set_lnum_n(p_set,
                                         p_id,
                                         1,
                                         CS_LAGR_CELL_ID,
                                         saved_cell_id[i]);
            cs_real_t *p_coords
              = cs_lagr_particles_attr_n(p_set,
                                         p_id,
                                         1,
                                         CS_LAGR_COORDS);
            for (cs_lnum_t j = 0; j < 3; j++)
              p_coords[j] = saved_coords[i][j];
          }

          BFT_FREE(saved_coords);
          BFT_FREE(saved_cell_id);

          /* Add particle tracking events for boundary injection */

          if (   particle_face_ids != NULL
              && cs_lagr_stat_is_active(CS_LAGR_STAT_GROUP_TRACKING_EVENT)) {

            cs_lagr_event_set_t  *events
              = cs_lagr_event_set_boundary_interaction();

            /* Event set "expected" size: n boundary faces*2 */
            cs_lnum_t events_min_size = mesh->n_b_faces * 2;
            if (events->n_events_max < events_min_size)
              cs_lagr_event_set_resize(events, events_min_size);

            for (cs_lnum_t i = 0; i < n_inject; i++) {
              cs_lnum_t p_id = particle_range[0] + i;

              cs_lnum_t event_id = events->n_events;
              events->n_events += 1;

              if (event_id >= events->n_events_max) {
                /* flush events */
                cs_lagr_stat_update_event(events,
                                          CS_LAGR_STAT_GROUP_TRACKING_EVENT);
                events->n_events = 0;
                event_id = 0;
              }

              cs_lagr_event_init_from_particle(events, p_set, event_id, p_id);

              cs_lnum_t face_id = particle_face_ids[i];
              cs_lagr_events_set_lnum(events,
                                      event_id,
                                      CS_LAGR_E_FACE_ID,
                                      face_id);

              cs_lnum_t *e_flag = cs_lagr_events_attr(events,
                                                      event_id,
                                                      CS_LAGR_E_FLAG);

              *e_flag = *e_flag | CS_EVENT_INFLOW;

            }

          }

          BFT_FREE(particle_face_ids);

        }

        /* check some particle attributes consistency */

        _check_particles(p_set, zis, n_z_elts, elt_particle_idx);

        /* update counters and balances */

        cs_real_t z_weight = 0.;

        for (cs_lnum_t p_id = particle_range[0];
             p_id < particle_range[1];
             p_id++) {
          cs_real_t s_weight = cs_lagr_particles_get_real(p_set, p_id,
                                                          CS_LAGR_STAT_WEIGHT);
          cs_real_t flow_rate = (  s_weight
                                 * cs_lagr_particles_get_real(p_set, p_id,
                                                              CS_LAGR_MASS));

          zd->particle_flow_rate[z_id*n_stats] += flow_rate;

          if (n_stats > 1) {
            int class_id = cs_lagr_particles_get_lnum(p_set, p_id,
                                                      CS_LAGR_STAT_CLASS);
            if (class_id > 0 && class_id < n_stats)
              zd->particle_flow_rate[z_id*n_stats + class_id] += flow_rate;
          }

          z_weight += s_weight;
        }

        p_set->n_particles += n_inject;
        p_set->n_part_new += n_inject;
        p_set->weight_new += z_weight;

      } /* end of loop on sets */

    } /* end of loop on zones */

  } /* end of loop on zone types (boundary/volume) */

  BFT_FREE(elt_particle_idx);

  /* Update global particle counters */

  pc = cs_lagr_update_particle_counter();
  pc->n_g_total += pc->n_g_new;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
