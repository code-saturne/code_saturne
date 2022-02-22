/*============================================================================
 * Methods for particle fragmentation modeling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * Functions dealing with the particle fragmentation modeling
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_prototypes.h"

#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_physical_constants.h"
#include "cs_parall.h"

#include "cs_lagr.h"

#include "cs_lagr_agglo.h"
#include "cs_lagr_clogging.h"
#include "cs_lagr_injection.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"

#include "cs_random.h"
#include "cs_lagr.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_fragmentation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert a new particle at the end of the particle set
 *
 * \param[in]  newpart                 number of particles previously inserted
 * \param[in]  vp                      statistical weight of inserted particle
 * \param[in]  corr                    particles indices of the particle set
 * \param[in]  newclass                class of the particle to be inserted
 * \param[in]  minimum_particle_diam   diameter of the monomere
 * \param[in]  idx                     index of the particle local to the
 *                                     current cell
 * \param[in]  mass                    mass of the particles
 *
 */
/*----------------------------------------------------------------------------*/

static void
_insert_particles(cs_lnum_t   newpart,
                  cs_lnum_t   vp,
                  cs_lnum_t  *corr,
                  cs_lnum_t   idx,
                  cs_lnum_t   newclass,
                  cs_real_t   minimum_particle_diam,
                  cs_real_t   mass)
{
  /* Create new particles (indices > n_particles) */
  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;
  cs_lnum_t inserted_parts = p_set->n_particles + newpart;

  cs_lagr_particle_set_resize(inserted_parts);
  cs_lagr_part_copy(inserted_parts-1, corr[idx]);

  /* Set properties of particles
   * (diameter, statistical weight, mass, velocity, velocity seen, cell id) */
  cs_lagr_particles_set_real(p_set, inserted_parts-1, CS_LAGR_STAT_WEIGHT, vp);

  cs_real_t fractal_dim = cs_lagr_particles_get_real(p_set, inserted_parts-1,
                                                     CS_LAGR_AGGLO_FRACTAL_DIM);
  cs_real_t diam = minimum_particle_diam * pow((cs_real_t)newclass,
                                               1./fractal_dim);

  cs_lagr_particles_set_real(p_set, inserted_parts-1, CS_LAGR_DIAMETER, diam);

  cs_lagr_particles_set_real(p_set, inserted_parts-1, CS_LAGR_MASS, mass);

  cs_real_t * inserted_vel = cs_lagr_particles_attr(p_set, inserted_parts-1,
                                                    CS_LAGR_VELOCITY);
  cs_real_t * idx_vel = cs_lagr_particles_attr(p_set, corr[idx],
                                               CS_LAGR_VELOCITY);
  inserted_vel[0] = idx_vel[0];
  inserted_vel[1] = idx_vel[1];
  inserted_vel[2] = idx_vel[2];

  cs_real_t * inserted_vel_f = cs_lagr_particles_attr(p_set, inserted_parts-1,
                                                      CS_LAGR_VELOCITY_SEEN);
  cs_real_t * idx_vel_f = cs_lagr_particles_attr(p_set, corr[idx],
                                                 CS_LAGR_VELOCITY_SEEN);
  inserted_vel_f[0] = idx_vel_f[0];
  inserted_vel_f[1] = idx_vel_f[1];
  inserted_vel_f[2] = idx_vel_f[2];

  cs_lnum_t iep = cs_lagr_particles_get_lnum(p_set, corr[idx], CS_LAGR_CELL_ID);
  cs_lagr_particles_set_lnum(p_set, inserted_parts-1, CS_LAGR_CELL_ID, iep);
  cs_lagr_particles_set_lnum(p_set, inserted_parts-1,
                             CS_LAGR_AGGLO_CLASS_ID,  newclass);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search in a sorted array
 *
 * \param[in]  array       sorted array
 * \param[in]  array_size  size of the sorted array
 * \param[in]  search      element that is searched in the array
 *
 * \return the index of the searched for element in the array (-1 if not found)
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_find_class(cs_lnum_t  array[][2],
            cs_lnum_t  array_size,
            cs_lnum_t  search)
{
  cs_lnum_t first = 0;
  cs_lnum_t last = array_size - 1;
  cs_lnum_t middle = (first+last)/2;

  while (first <= last) {
    if (array[middle][0] < search)
      first = middle + 1;
    else if (array[middle][0] == search) {
      return middle;
    }
    else
      last = middle - 1;

    middle = (first + last)/2;
  }

  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Adds a new physical particle to the simulation
 *        If a particle with the same size (class) already exists, merge them
 *        if not, create a new parcel at the end of particle set.
 *
 * \param[in]  lnum_particles          number of particles in cell
 *                                     before fragmentation
 * \param[in]  newpart                 pointer to number of particles
 *                                     created by fragmentation
 * \param[in]  vp                      statistical weight of the fragmented
 *                                     particle to be added.
 * \param[in]  corr                    array of global indices in particle set
 * \param[in]  frag_idx                local to the cell, index of particle
 *                                     to add to the particle set
 * \param[in]  newclass                class of the new fragment
 * \param[in]  minimum_particle_diam   minumum diameter (monomere diameter)
 * \param[in]  mass                    mass of the particles
 * \param[in]  agglo_max_weight                 maximum statistical weight that a
 *                                     particle can have
 */
/*----------------------------------------------------------------------------*/

static void
_add_particle(cs_lnum_t   lnum_particles,
              cs_lnum_t  *newpart,
              cs_lnum_t   vp,
              cs_lnum_t  *corr,
              cs_lnum_t   frag_idx,
              cs_lnum_t   newclass,
              cs_real_t   minimum_particle_diam,
              cs_real_t   mass,
              cs_real_t   agglo_max_weight,
              cs_lnum_t   interf[][2])
{
  /* Get information on the new fragment*/
  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;

  /* Find an existing particle in the same class*/
  cs_lnum_t position = _find_class(interf, lnum_particles, newclass);

  /* Merge with existing particle (if possible) */
  if (position >= 0) {
    cs_lnum_t found_idx = interf[position][1];
    cs_real_t stat_weight = cs_lagr_particles_get_real(p_set, found_idx,
                                                       CS_LAGR_STAT_WEIGHT);
    if (stat_weight + vp <= agglo_max_weight) {
      long long int auxx = round(stat_weight);
      cs_lagr_particles_set_real(p_set, found_idx, CS_LAGR_STAT_WEIGHT, auxx+vp);
      return;
    }
  }

  /* Add a new particle at the end of the set (otherwise)*/
  cs_lnum_t add_to_end = 1;
  for (cs_lnum_t indx = p_set->n_particles;
       indx < p_set->n_particles + *newpart;
       indx++) {
    int stat_class = cs_lagr_particles_get_lnum(p_set, indx,
                                                CS_LAGR_AGGLO_CLASS_ID);
    cs_real_t stat_weight = cs_lagr_particles_get_real(p_set, indx,
                                                       CS_LAGR_STAT_WEIGHT);
    if ((stat_class == newclass)
       && (stat_weight + vp <= agglo_max_weight)) {
      long long int auxx = round(stat_weight);
      cs_lagr_particles_set_real(p_set, indx, CS_LAGR_STAT_WEIGHT, auxx+vp);

      add_to_end = 0;
      return;
    }
  }

  if (add_to_end) {
    (*newpart)++;
    _insert_particles(*newpart, vp, corr, frag_idx, newclass,
                      minimum_particle_diam, mass);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Obtain the number of mesh cells occupied by at least one particle.
 *
 * \param[in]  array       TODO describe this
 * \param[in]  array_size  size of array
 *
 * \returns:
 *   integer that gives the number of cells occupied by at least
 *   one particle in a particle sub-set
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_get_nb_classes(cs_lnum_t  array[][2],
                cs_lnum_t  array_size)
{
  /* Count number of cells that contain particles */
  cs_lnum_t count_cls = 0;
  if (array_size > 1) {
    count_cls = 1;
    cs_lnum_t prev_cls = array[0][0];
    cs_lnum_t curr_cls = prev_cls;
    for (cs_lnum_t p = 1; p < array_size; ++p) {
      curr_cls = array[p][0];
      if (prev_cls != curr_cls) {
        count_cls++;
      }
      prev_cls = curr_cls;
    }
  }
  else if (array_size == 1) {
    count_cls = 1;
  }

  return count_cls;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Obtain the index of cells that are occupied by at least one
 *         particles and the list of indices of the particle set that contain
 *         the first element in a cell
 *
 * The two arrays that contain this information need to be pre-allocated to size
 * n_particle_cells and n_particle_cells+1 respectively.
 *
 * The value of particle_gaps at index i contains the first particle that
 * is found in cell (i+1). The last element of particle_gaps contains the size
 * of the particle set.
 *
 * \param[in]   array       TODO describe this
 * \param[in]   array_size  size of array
 * \param[in]   count_cls   number of cells that are are occupied by particles
 * \param[out]  cls_gaps    preallocated list of size count_csl+1
 *                          that will contain the starting indices of particles
 *                          that are in the current cell
 */
/*----------------------------------------------------------------------------*/

static void
_gaps_classes(cs_lnum_t  array[][2],
              cs_lnum_t  array_size,
              cs_lnum_t  count_cls,
              cs_lnum_t  cls_gaps[])
{
  if (array_size >= 1) {
    cs_lnum_t prev_cls = array[0][0];
    cs_lnum_t curr_cls = prev_cls;
    cs_lnum_t counter = 0;
    cls_gaps[0] = 0;

    counter = 1;

    for (cs_lnum_t part = 1; part < array_size; ++part) {
      curr_cls = array[part][0];
      if (prev_cls != curr_cls) {
        cls_gaps[counter] = part;
        counter++;
      }
      prev_cls = curr_cls;
    }
    cls_gaps[count_cls] = array_size;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compare the first elements of a cs_lnum_2_t. Used in sorting
 */
/*----------------------------------------------------------------------------*/

static int
_compare_interface(const void  *a,
                   const void  *b)
{
  /* Map as 1-d array for simplification and to avoid const warnings
     ((* (const cs_lnum_2_t *)a)[0] - (* (const cs_lnum_2_t *) b)[0]) */

  return ((*(const cs_lnum_t *)a) - (*(const cs_lnum_t *)b));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Particle fragmentation algorithm.
 *
 * Parcels represent physical particles with similar properties (size).
 * The number of physical particles in a parcel is represented by the
 * statistical weight.
 *
 * For each parcel, the number of fragmentation events is generated with
 * a Poisson distribution that depends on the fragmentation kernel.
 *
 * Working hypotheses:
 *  1) Discrete diameters
 *     - Minimal particle size , called a monomere (unbreakable particle)
 *  2) Fragmentation is binary
 *     - Equally sized fragments (N0/2, or (N0+-1)/2)
 *  3) Fragmentation is treated after agglomeration
 *
 * Warning:  particle indices are not necessarily contiguous
 *           (due to agglomeration occuring before).
 *
 * Two contiguous groups: particles present before agglomeration
 *                        newly created particles (due to agglomeration)
 *
 * \param[in]  dt                     time step
 * \param[in]  minimum_particle_diam  minumum diameter (monomere diameter)
 * \param[in]  main_start             index of the first particle in cell
 *                                    present before the agglomeration
 * \param[in]  main_end               index after the last particle in cell,
 *                                    present before the agglomeration
 * \param[in]  agglo_start            index of the first particle in cell,
 *                                    created by the agglomeration
 * \param[in]  agglo_end              index after the last particle in cell,
 *                                    created by the agglomeration
 *
 * \returns  modified list of particles containing newly created parcels
 *           at the end
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_lagr_fragmentation(cs_real_t  dt,
                      cs_real_t  minimum_particle_diam,
                      cs_lnum_t  main_start,
                      cs_lnum_t  main_end,
                      cs_lnum_t  agglo_start,
                      cs_lnum_t  agglo_end)
{
  /* Initialization */
  cs_lnum_t ret_val = 0;
  cs_real_t agglo_min_weight
    = cs_glob_lagr_agglomeration_model->min_stat_weight; // Minimum parcel size
  cs_real_t agglo_max_weight
    = cs_glob_lagr_agglomeration_model->max_stat_weight; // Maximum parcel size

  cs_lnum_t newpart = 0;

  /* Get mesh properties and particle set */
  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;

  /* Number of particles in the cell */
  cs_lnum_t lnum_particles = main_end - main_start;
  lnum_particles += agglo_end - agglo_start;

  /* Exit routine if no particles are in the cell */
  if (lnum_particles <=0) {
    return 0;
  }

  /* Create local array (containing the particle index) */
  cs_lnum_t *corr;
  BFT_MALLOC(corr, lnum_particles, cs_lnum_t);

  /* Create local list (containing the class and the global index) */
  cs_lnum_2_t * interf;
  BFT_MALLOC(interf, lnum_particles, cs_lnum_2_t);

  /* Browse the list of particles already existing */
  for (cs_lnum_t idx = main_start; idx < main_end; ++idx) {
    corr[idx - main_start] = idx;

    cs_lnum_t curr_cls = cs_lagr_particles_get_lnum(p_set, idx,
                                                    CS_LAGR_AGGLO_CLASS_ID);
    interf[idx-main_start][0] = curr_cls;
    interf[idx-main_start][1] = idx;
  }

  /* Browse the list of particle newly created by agglomeration */
  for (cs_lnum_t idx = agglo_start; idx < agglo_end; ++idx) {
    corr[idx - agglo_start + (main_end-main_start)] = idx;

    cs_lnum_t curr_cls = cs_lagr_particles_get_lnum(p_set, idx,
                                                    CS_LAGR_AGGLO_CLASS_ID);
    interf[idx - agglo_start + main_end - main_start][0] = curr_cls;
    interf[idx - agglo_start + main_end - main_start][1] = idx;
  }

  /* Sort particles by class */
  qsort(interf, lnum_particles, sizeof(cs_lnum_2_t), _compare_interface);

  /* Get fragmentation kernel */
  cs_real_t cker = 0.;
  cker = cs_glob_lagr_fragmentation_model->scalar_kernel;

  for (cs_lnum_t i = 0; i < lnum_particles; ++i) {

    if (cs_lagr_particles_get_flag(p_set, corr[i], CS_LAGR_PART_TO_DELETE))
      continue;

    /* Generate the number of fragmentation events */
    int vp = 0;
    cs_real_t rand;

    cs_lnum_t class_nb = cs_lagr_particles_get_lnum(p_set, corr[i],
                                                    CS_LAGR_AGGLO_CLASS_ID);
    if (class_nb > 1) {
      cs_real_t stat_weight = cs_lagr_particles_get_real(p_set, corr[i],
                                                         CS_LAGR_STAT_WEIGHT);

      cs_real_t lambda = cker * stat_weight * dt;
      // FIXME: change poisson random generator
      if (lambda > 700.) {
        cs_random_normal(1, &rand);
        vp = floor(lambda + sqrt(lambda) * rand);
      }
      else {
        cs_random_poisson(1, lambda, &vp);
      }

      /* Treat fragmentation events */
      if (vp > 0) {
        if (stat_weight < vp) {
          vp = stat_weight;
        }
        cs_real_t new_stat_weight = stat_weight - vp;
        cs_lagr_particles_set_real(p_set, corr[i],
                                   CS_LAGR_STAT_WEIGHT, new_stat_weight);
        cs_real_t mass = cs_lagr_particles_get_real(p_set, corr[i], CS_LAGR_MASS);
        /* Binary fragmentation (fragments with the same size) */
        if (class_nb%2 != 0) {
          cs_lnum_t class_nb_1 = class_nb / 2;
          cs_lnum_t class_nb_2 = class_nb - class_nb_1;

          _add_particle(lnum_particles, &newpart, vp, corr, i, class_nb_1,
                        minimum_particle_diam, mass*class_nb_1/class_nb,
                        agglo_max_weight, interf);
          _add_particle(lnum_particles, &newpart, vp, corr, i, class_nb_2,
                        minimum_particle_diam, mass*class_nb_2/class_nb,
                        agglo_max_weight, interf);
        }
        else {
          cs_lnum_t class_nb_even = class_nb / 2;
          _add_particle(lnum_particles, &newpart, 2*vp, corr, i, class_nb_even,
                        minimum_particle_diam, mass*0.5, agglo_max_weight,
                        interf);
        }
      }
    }
  }

  /* Local array to save new fragments (class, index) */
  cs_lnum_2_t *interf_frag;
  BFT_MALLOC(interf_frag, newpart, cs_lnum_2_t);

  for (cs_lnum_t i=0; i<newpart; ++i) {
    cs_lnum_t curr_class
      = cs_lagr_particles_get_lnum(p_set, p_set->n_particles+i,
                                   CS_LAGR_AGGLO_CLASS_ID);
    interf_frag[i][0] = curr_class;
    interf_frag[i][1] = p_set->n_particles+i;
  }

  /* Sort new fragments by class */
  qsort(interf_frag, newpart, sizeof(cs_lnum_2_t), _compare_interface);

  /* Merge new parcels with low statistical weight */
  cs_lnum_t delta_part = main_end - main_start;
  cs_lnum_t delta_agglo = agglo_end - agglo_start;

  cs_lnum_t tot_size = delta_part + delta_agglo + newpart;

  /* Sort all particles in a cell by class */
  cs_lnum_2_t *interf_tot;
  BFT_MALLOC(interf_tot, tot_size, cs_lnum_2_t);

  cs_lagr_agglo_merge_arrays(interf, interf_frag,
                             lnum_particles, newpart, interf_tot);

  BFT_FREE(interf);
  BFT_FREE(interf_frag);

  /* Merge new parcels with low statistical weight */
  cs_lnum_t nb_cls = _get_nb_classes(interf_tot, tot_size);
  /* Split the local particle set by class to facilitate merging */
  cs_lnum_t* cls_gaps;

  BFT_MALLOC(cls_gaps, tot_size+1, cs_lnum_t);

  _gaps_classes(interf_tot, tot_size,
                nb_cls, cls_gaps);

  /* Try to merge particles of the same class */
  for (cs_lnum_t i = 0; i < nb_cls; ++i) {
    cs_lnum_t start_gap = cls_gaps[i];
    cs_lnum_t end_gap = cls_gaps[i+1];

    if (end_gap - start_gap == 1) {
      cs_lnum_t part_idx = interf_tot[start_gap][1];
      cs_real_t weight = cs_lagr_particles_get_real(p_set, part_idx,
                                                    CS_LAGR_STAT_WEIGHT);
      /* Eliminate if statistical particle 0 */
      if (weight <= 0.) {
        cs_lagr_particles_set_lnum(p_set, part_idx,
                                   CS_LAGR_AGGLO_CLASS_ID, 0);
        cs_lagr_particles_set_real(p_set, part_idx,
                                   CS_LAGR_STAT_WEIGHT, 0);
        cs_lagr_particles_set_flag(p_set, part_idx, CS_LAGR_PART_TO_DELETE);
        // CELL NUM TODO ask value to delete
      }
      continue;
    }

    /* Add all small weights */
    cs_real_t sum = 0.;
    cs_lnum_t last_small = -1;
    cs_lnum_t found_small = 0;

    for (cs_lnum_t idx = start_gap; idx < end_gap; ++idx) {
      cs_lnum_t part_idx = interf_tot[idx][1];
      cs_real_t weight = cs_lagr_particles_get_real(p_set, part_idx,
                                                    CS_LAGR_STAT_WEIGHT);

      if (weight > 0. && weight < agglo_min_weight) {
        last_small = idx;
        sum += weight;
        found_small = 1;
      }
    }

    if (found_small) {
      /* Put small particles in a large one (if possible) */
      cs_lnum_t put_in_large = 0;
      for (cs_lnum_t idx = start_gap; idx < end_gap; ++idx) {
        cs_lnum_t part_idx = interf_tot[idx][1];
        cs_real_t weight = cs_lagr_particles_get_real(p_set, part_idx,
                                                      CS_LAGR_STAT_WEIGHT);

        if (weight >= agglo_min_weight && sum + weight < agglo_max_weight) {
          put_in_large = 1;
          cs_lagr_particles_set_real(p_set, part_idx, CS_LAGR_STAT_WEIGHT,
                                     sum+weight);
          break;
        }
      }

      /* if small particles put in a large one, then all small particles removed,
         else put all small particles in the last one */
      if (put_in_large) {
        last_small = -1;
      }
      else {
        cs_lagr_particles_set_real(p_set, interf_tot[last_small][1],
                                   CS_LAGR_STAT_WEIGHT, sum);
      }

      for (cs_lnum_t idx = start_gap; idx < end_gap; ++idx) {
        cs_lnum_t part_idx = interf_tot[idx][1];
        cs_real_t weight = cs_lagr_particles_get_real(p_set, part_idx,
                                                      CS_LAGR_STAT_WEIGHT);

        if (weight > 0. && weight < agglo_min_weight && idx != last_small) {
          cs_lagr_particles_set_lnum(p_set, part_idx,
                                     CS_LAGR_AGGLO_CLASS_ID, 0);
          cs_lagr_particles_set_real(p_set, part_idx,
                                     CS_LAGR_STAT_WEIGHT, 0);
          cs_lagr_particles_set_flag(p_set, part_idx, CS_LAGR_PART_TO_DELETE);
          // CELL NUM TODO ask value to delete
        }
      }
    }

    /* Eliminate particles (if statistical weight < 0) */
    for (cs_lnum_t idx = start_gap; idx < end_gap; ++idx) {
      cs_lnum_t part_idx = interf_tot[idx][1];
      cs_real_t weight = cs_lagr_particles_get_real(p_set, part_idx,
                                                    CS_LAGR_STAT_WEIGHT);
      if (weight <= 0.) {
        cs_lagr_particles_set_lnum(p_set, part_idx,
                                   CS_LAGR_AGGLO_CLASS_ID, 0);
        cs_lagr_particles_set_real(p_set, part_idx,
                                   CS_LAGR_STAT_WEIGHT, 0);
        cs_lagr_particles_set_flag(p_set, part_idx, CS_LAGR_PART_TO_DELETE);
        //CELL NUM TODO ask value to delete
      }
    }
  }

  BFT_FREE(cls_gaps);
  BFT_FREE(interf_tot);

  p_set->n_particles += newpart;

  BFT_FREE(corr);

  return ret_val;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
