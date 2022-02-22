/*============================================================================
 * Methods for particle agglomeration modeling
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
 * Functions dealing with the particle clogging modeling
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
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_agglo.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
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
#include "cs_prototypes.h"

#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"

#include "cs_random.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_agglo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Obtain the (i,j) index of an upper triangular matrix
 *        from the index of the flattened matrix
 *
 * \param[out]  res     result that contains the row column indices
 * \param[in]   idx     index of the flattened array
 * \param[in]   size    number of rows (same as the number of columns
 *                      of the upper triangular matrix
 */
/*----------------------------------------------------------------------------*/

static void
_get_row_col(cs_lnum_t  res[2],
             cs_gnum_t  idx,
             cs_lnum_t  size)
{
  unsigned long long aux =   (unsigned long long)size
                           * ((unsigned long long)size + 1) / 2 - idx;
  /* check if aux is triangular number */
  bool is_triag = false;
  unsigned long long square_root = (unsigned long long)round((sqrt(8*aux+1)));
  if (square_root*square_root == (8*aux+1) ) {
    is_triag = true;
  }

  cs_lnum_t k = -1;
  if (is_triag) {
    k = (-1 + square_root) / 2-1;
  }
  else {
    aux = sqrt(8*aux + 1);
    k = floor((-1 + aux) / 2);
  }

  cs_lnum_t i = size - k - 1;
  aux = idx + (unsigned long long)i    * ((unsigned long long)i+1) / 2
            - (unsigned long long)size *  (unsigned long long)i;
  cs_lnum_t j = (cs_lnum_t)aux;

  res[0] = i;
  res[1] = j;
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary serch in a sorted array
 *
 * \param[in]  array       sorted array
 * \param[in]  array_size  size of the sorted array
 * \param[in]  search      element that is searched in the array
 *
 * \return the index of the searched for element in the array.
 *         if not found, returns -1
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
    if (array[middle][0] < search) {
      first = middle + 1;
    }
    else if (array[middle][0] == search) {
      return middle;
    }
    else
      last = middle - 1;

    middle = (first + last) /2;
  }
  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Obtain the number of mesh cells occupied by at least one particle.
 *
 * \param[in]        p_set             pointer to particle data structure
 * \param[in]        start             start position in p_set
 * \param[in]        end               end position in p_se
 *
 * \returns:
 *   integer that gives the number of cells occupied by at least one
 *           particle in a particle sub-set
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
  }else if (array_size == 1) {
    count_cls = 1;
  }

  return count_cls;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Obtain the index of cells that are occupied by at least one
 *         particle and the list of indices of the particle set that contain
 *         the first element in a cell
 *
 * The two arrays that contain this information need to be pre-allocated to size
 * counter_particle_cells and counter_particle_cells+1 respectively.
 *
 * The value of cls_gaps at index i contains the first particle in
 * p_set that is found in cell[i+1]. The last element of particle_gaps contains
 * the size of p_set
 *
 * \param[in]   array       per-cell with particles info
 * \param[in]   array_size  array size
 * \param[in]   count_cls   number of cells occupied by particles
 * \param[out]  cls_gaps    preallocated list of size count_cls
 *                          that will contain the starting indices
 *                          of particles that are in the current cell
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_gaps_classes(cs_lnum_2_t  *array,
              cs_lnum_t     array_size,
              cs_lnum_t     count_cls,
              cs_lnum_t    *cls_gaps)
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

  return 0;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge two sorted arrays in a third sorted array
 *
 * \param[in]       arr1   first sorted array
 * \param[in]       arr2   second sorted array
 * \param[in]       n1     size of first sorted array
 * \param[in]       n2     size of second sorted array
 * \param[in, out]  arr3   preallocated array that will contain the sorted
 *                         merge of the two previous arrays
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_agglo_merge_arrays(cs_lnum_2_t  arr1[],
                           cs_lnum_2_t  arr2[],
                           cs_lnum_t    n1,
                           cs_lnum_t    n2,
                           cs_lnum_2_t  arr3[])
{
  cs_lnum_t i = 0, j = 0, k = 0;
  /* Browse both arrays  */
  while (i<n1 && j <n2) {
    if (arr1[i][0] < arr2[j][0]) {
      arr3[k][0] = arr1[i][0];
      arr3[k][1] = arr1[i][1];
      k++;
      i++;
    }
    else {
      arr3[k][0] = arr2[j][0];
      arr3[k][1] = arr2[j][1];
      k++;
      j++;
    }
  }

  /* Store remaining elements of first array  */
  while (i < n1) {
    arr3[k][0] = arr1[i][0];
    arr3[k][1] = arr1[i][1];
    k++;
    i++;
  }

  /* Store remaining elements of second array */
  while (j < n2) {
    arr3[k][0] = arr2[j][0];
    arr3[k][1] = arr2[j][1];
    k++;
    j++;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Agglomeration algorithm based on algorithms used in
 *        rare gas modelling.
 *
 * Parcels represent physical particles with similar properties (size).
 * The number of physical particles in a parcel is represented by
 * the statistical weight.
 *
 * - We select randomly a number of parcels within which agglomeration
 *   is treated.
 * - We selected randomly pairs of particles and, for each pair,
 *   the number of agglomeration events is generated with
 *   a Poisson distribution that depends on the agglomeration kernel and
 *   number of each particle (i.e. statistical weight).
 *
 * Working hypotheses:
 *
 * 1) Discrete diameters
 *    - Minimal particle size , called a monomere (unbreakable particle)
 *    - Aggregates correponds to group of monomeres sticking together.
 *    - Each class of the parcel represent one size (stored in CS_LAGR_AGGREGATE)
 * 2) Agglomeration happens between two parcels
 * 3) Particles in the same cell are contiguous in the particle list
 *
 * \param[in]  cell_id                current cell id
 * \param[in]  dt                     time step
 * \param[in]  minimum_particle_diam  minumum diameter (monomere diameter)
 * \param[in]  start_particle         index of the first particle
 * \param[in]  end_particle           index after the last particle
 *
 * \returns a modified list of particles, containing newly
 *          created parcels at the end
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_lagr_agglomeration(cs_lnum_t  cell_id,
                      cs_real_t  dt,
                      cs_real_t  minimum_particle_diam,
                      cs_lnum_t  start_particle,
                      cs_lnum_t  end_particle)
{
  /* Initialisation */
  cs_lnum_t ret_val = 0;

  // FIXME: call DLVO routine
  cs_real_t alp = 1.;                    /* Efficiency of agglomeration */

  cs_real_t agglo_min_weight
    = cs_glob_lagr_agglomeration_model->min_stat_weight; // Minimum parcel size
  cs_real_t agglo_max_weight
    = cs_glob_lagr_agglomeration_model->max_stat_weight; // Maximum parcel size
  cs_lnum_t newpart = 0;

  /* Get fluid and particle properties */
  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  /* Get cell information */
  cs_real_t vol_cell = fvq->cell_vol[cell_id];
  /* local number of particles  */
  cs_lnum_t lnum_particles = end_particle - start_particle;

  /* Exit routine if no particles are in the cell */
  if (lnum_particles <=0) {
    return 0;
  }

  /* Create local array (containing the class and particle index) */
  cs_lnum_2_t *interf;

  BFT_MALLOC(interf, lnum_particles, cs_lnum_2_t);

  for (cs_lnum_t i = start_particle; i < end_particle; ++i) {
    interf[i-start_particle][0]
      = cs_lagr_particles_get_lnum(p_set, i, CS_LAGR_AGGLO_CLASS_ID);
    interf[i-start_particle][1] = i;
  }

  /* Sort by class */
  qsort(interf, lnum_particles, sizeof(cs_lnum_2_t), _compare_interface);

  /* Select pairs of particles (for agglomeration) */
  cs_gnum_t _gn_particles = lnum_particles;
  cs_gnum_t lnum_maxpairs = _gn_particles*(_gn_particles+1) / 2;

  cs_gnum_t lnum_colpairs = 1;
  if (floor(0.5 * lnum_particles * log(lnum_particles+0.1))  > 0) {
    // magic constant 0.01 to avoid log(0)
    lnum_colpairs = floor(0.5 * lnum_particles * log(lnum_particles + 0.1));
  }

  signed long long kk = 0;  /* Number of particle pairs */
  cs_gnum_t selected_val = lnum_colpairs;
  if (selected_val <= 0) {
    selected_val = 1;
  }

  kk = selected_val - 1;
  long long int vp = 0;
  cs_lnum_t n_classes_new = 0;

  /* Treat agglomeration between pairs*/
  while (kk >= 0) {
    cs_real_t rand;

    /* Evaluate number of agglomeration events */
    vp = 0;

    cs_random_uniform (1, &rand);
    cs_gnum_t pp = floor(lnum_maxpairs * rand);

    // TODO check if cs_random_uniform < 1 and not <= 1
    if (pp == lnum_maxpairs) {
      pp=0;
    }

    cs_lnum_t p1, p2;

    cs_lnum_2_t pps;
    _get_row_col(pps, pp, lnum_particles);

    p1 = pps[0];
    p2 = pps[1];

    cs_real_t cker = 0, lambda = 0;
    cs_real_2_t stat_weights;  /* Weights of the two particles */

    stat_weights[0] = cs_lagr_particles_get_real
                        (p_set, p1+start_particle, CS_LAGR_STAT_WEIGHT);
    stat_weights[1] = cs_lagr_particles_get_real
                        (p_set, p2+start_particle, CS_LAGR_STAT_WEIGHT);

    cker = cs_glob_lagr_agglomeration_model->scalar_kernel;

    lambda = cker * (lnum_maxpairs / selected_val)
                  * stat_weights[0] * stat_weights[1] * dt / vol_cell;

    /* If auto agglomeration, divide by half */
    lambda *= 0.5 * (p1==p2) + (p1 != p2);

    /* Efficiency of agglomeration */
    lambda *= alp;

    // FIXME: change Poisson random generator
    if (lambda > 700.) {
      cs_random_normal(1, &rand);
      vp = floor(lambda + sqrt(lambda) * rand);
    }
    else {
      int _vp;
      cs_random_poisson(1, lambda, &_vp);
      vp = _vp;
    }

    /* Clip if vp > stat_weight */
    if ( vp > 0 ) {
      cs_lnum_t max_vp = round( cs_math_fmin(stat_weights[0], stat_weights[1]) );
      /* Auto-agglomeration */
      if (p1 == p2 && 2*vp > max_vp )
        vp = floor( 0.5*max_vp );
      /* Agglomeration between different clusters */
      if ( p1 != p2 && vp > max_vp )
        vp = max_vp ;
    }
    ret_val = 1;

    /* Treat agglomeration events */
    if (vp > 0) {

      /* Remove elements from parcels p1 and p2 */
      cs_lagr_particles_set_real(p_set, p1+start_particle,
                                 CS_LAGR_STAT_WEIGHT, round(stat_weights[0])- vp );
      stat_weights[1] = cs_lagr_particles_get_real(p_set, p2+start_particle,
                                                   CS_LAGR_STAT_WEIGHT);
      cs_lagr_particles_set_real(p_set, p2+start_particle,
                                 CS_LAGR_STAT_WEIGHT, round(stat_weights[1]) - vp);

      /* Add elements for parcel (p1+p2)
       * --> either add to an existing parcel (if it exists and is not too big)
       * --> or create a new parcel (otherwise) */
      n_classes_new =   cs_lagr_particles_get_lnum
                          (p_set, p1+start_particle, CS_LAGR_AGGLO_CLASS_ID)
                     + cs_lagr_particles_get_lnum(p_set, p2+start_particle,
                                                  CS_LAGR_AGGLO_CLASS_ID);

      if (n_classes_new > cs_glob_lagr_agglomeration_model->n_max_classes) {
        bft_error(__FILE__, __LINE__, 0,
                  _(" ** Lagrangian module:\n"
                    "    number of classes > n_max_classes (%d)\n"
                    "    --------------------------------\n\n"
                    " To fix this, increase "
                    " cs_glob_lagr_agglomeration_model->n_max_classes"),
                  cs_glob_lagr_agglomeration_model->n_max_classes);
      }

      /* Find one existing parcel with the same class (to merge with it) */
      cs_lnum_t position = _find_class(interf, lnum_particles, n_classes_new);

      if (position >= 0) {
        cs_lnum_t found_idx = interf[position][1];
        cs_real_t stat_weight
          = cs_lagr_particles_get_real(p_set, found_idx, CS_LAGR_STAT_WEIGHT);

        if (stat_weight + vp <= agglo_max_weight) {
          cs_lagr_particles_set_real(p_set, found_idx,
                                     CS_LAGR_STAT_WEIGHT, round(stat_weight)+vp);

          kk--;
          continue;
        }
      }

      /* Else, find new aggregated particles with the same class
         (to merge with it) */

      cs_lnum_t add_to_end = 1;

      for (cs_lnum_t indx = p_set->n_particles;
          indx < p_set->n_particles + newpart;
          indx++) {
        cs_lnum_t stat_class
          = cs_lagr_particles_get_lnum(p_set, indx, CS_LAGR_AGGLO_CLASS_ID);
        cs_real_t stat_weight
          = cs_lagr_particles_get_real(p_set, indx, CS_LAGR_STAT_WEIGHT);
        if (   (stat_class == n_classes_new)
            && (stat_weight + vp <= agglo_max_weight)) {
          cs_lagr_particles_set_real(p_set, indx, CS_LAGR_STAT_WEIGHT,
                                     round(stat_weight)+vp);

          add_to_end = 0;
          break;
        }
      }

      /* Else, create a new parcel at the end
         Principle: copy parcel p1 and modify its properties */
      if ( add_to_end == 1 ) {
        newpart++;

        /* Copy parcel p1 into a new parcel */
        cs_lnum_t inserted_parts = p_set->n_particles + newpart;

        cs_lagr_particle_set_resize(inserted_parts);

        cs_lagr_part_copy(inserted_parts-1, p1+start_particle);

        /* Set statistical weight*/
        cs_lagr_particles_set_real(p_set, inserted_parts-1,
                                   CS_LAGR_STAT_WEIGHT, vp);

        /* Set diameter (using a law based on fractal dimension of aggregates) */
        cs_real_t fractal_dim
          = cs_lagr_particles_get_real(p_set, inserted_parts-1,
                                       CS_LAGR_AGGLO_FRACTAL_DIM);
        cs_real_t diam = minimum_particle_diam * pow((cs_real_t)n_classes_new,
                                                     1./fractal_dim);
        cs_lagr_particles_set_real(p_set, inserted_parts-1,
                                   CS_LAGR_DIAMETER, diam);

        /* Set mass (equal to the sum of the two mass) */
        cs_real_t mass1
          = cs_lagr_particles_get_real(p_set, p1+start_particle, CS_LAGR_MASS);
        cs_real_t mass2
          = cs_lagr_particles_get_real(p_set, p2+start_particle, CS_LAGR_MASS);

        cs_lagr_particles_set_real(p_set, inserted_parts-1, CS_LAGR_MASS,
                                   mass1 +mass2);

        /* Set cell_id and Class_id */
        cs_lagr_particles_set_lnum(p_set, inserted_parts-1,
                                   CS_LAGR_CELL_ID, cell_id);

        cs_lagr_particles_set_lnum(p_set, inserted_parts-1,
                                   CS_LAGR_AGGLO_CLASS_ID, n_classes_new);

        /* Set particle velocity */

        cs_real_t * inserted_vel
          = cs_lagr_particles_attr(p_set, inserted_parts-1,
                                   CS_LAGR_VELOCITY);
        cs_real_t * p1_vel = cs_lagr_particles_attr(p_set, p1+start_particle,
                                                    CS_LAGR_VELOCITY);
        inserted_vel[0] = p1_vel[0];
        inserted_vel[1] = p1_vel[1];
        inserted_vel[2] = p1_vel[2];

        /* Set particle velocity seen*/
        cs_real_t * inserted_vel_seen
          = cs_lagr_particles_attr(p_set, inserted_parts-1,
                                   CS_LAGR_VELOCITY_SEEN);
        cs_real_t * p1_vel_seen
          = cs_lagr_particles_attr(p_set, p1+start_particle,
                                   CS_LAGR_VELOCITY_SEEN);
        inserted_vel_seen[0] = p1_vel_seen[0];
        inserted_vel_seen[1] = p1_vel_seen[1];
        inserted_vel_seen[2] = p1_vel_seen[2];

      }
    }
    kk--;
  }

  /* Store class and index of newly created particles */
  cs_lnum_2_t *interf_agglo;
  BFT_MALLOC(interf_agglo, newpart, cs_lnum_2_t);

  for (cs_lnum_t i = 0; i < newpart; i++) {
    cs_lnum_t curr_class = cs_lagr_particles_get_lnum
                             (p_set, p_set->n_particles+i,
                              CS_LAGR_AGGLO_CLASS_ID);
    interf_agglo[i][0] = curr_class;
    interf_agglo[i][1] = p_set->n_particles+i;
  }

  /* Sort by the class */
  qsort(interf_agglo, newpart, sizeof(cs_lnum_2_t), _compare_interface);

  /* Local array, containing all particles in the current cell,
     sorted by their class */
  cs_lnum_2_t  *interf_tot;
  cs_lnum_t tot_size = lnum_particles+newpart;

  BFT_MALLOC(interf_tot, tot_size, cs_lnum_2_t);

  /* Merge arrays of existing particles and particles
     created by agglomeration */
  cs_lagr_agglo_merge_arrays(interf, interf_agglo,
                             lnum_particles, newpart,
                             interf_tot);
  BFT_FREE(interf_agglo);
  BFT_FREE(interf);

  cs_lnum_t nb_cls = _get_nb_classes(interf_tot, tot_size);

  cs_lnum_t* cls_gaps;

  BFT_MALLOC(cls_gaps, tot_size+1, cs_lnum_t);

  _gaps_classes(interf_tot, tot_size,
                nb_cls, cls_gaps);

  /* Loop to merge particles of same class
     --> Keep only particles with large weight */
  for (cs_lnum_t i = 0; i < nb_cls; ++i) {
    cs_lnum_t start_gap = cls_gaps[i];
    cs_lnum_t end_gap = cls_gaps[i+1];

    /* Only one particle of current class */
    if (end_gap - start_gap == 1) {
      cs_lnum_t part_idx = interf_tot[start_gap][1];
      cs_real_t weight = cs_lagr_particles_get_real
                          (p_set, part_idx, CS_LAGR_STAT_WEIGHT);
      /* Delete particle (if weight < 0) */
      if (weight <= 0.) {
        cs_lagr_particles_set_lnum(p_set, part_idx,
                                   CS_LAGR_AGGLO_CLASS_ID, 0);
        cs_lagr_particles_set_real(p_set, part_idx, CS_LAGR_STAT_WEIGHT, 0);
        cs_lagr_particles_set_lnum(p_set, part_idx,
                                   CS_LAGR_P_FLAG, CS_LAGR_PART_TO_DELETE);


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
      /* Put the small particles in a large one (if possible) */
      cs_lnum_t put_in_large = 0;
      for (cs_lnum_t idx = start_gap; idx < end_gap; ++idx) {
        cs_lnum_t part_idx = interf_tot[idx][1];
        cs_real_t weight = cs_lagr_particles_get_real
                             (p_set, part_idx, CS_LAGR_STAT_WEIGHT);

        if (weight >= agglo_min_weight && sum + weight < agglo_max_weight) {
          put_in_large = 1;
          cs_lagr_particles_set_real(p_set, part_idx,
                                     CS_LAGR_STAT_WEIGHT, sum+weight);
          break;
        }
      }

      /* if small particles put in a large one, then all small particles
         removed, else put all small particles in the last one */
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
          cs_lagr_particles_set_real(p_set, part_idx, CS_LAGR_STAT_WEIGHT, 0);
          cs_lagr_particles_set_lnum(p_set, part_idx,
                                     CS_LAGR_P_FLAG, CS_LAGR_PART_TO_DELETE);

        }
      }
    }

    /* Eliminate particles (if statistical weight < 0) */
    for (cs_lnum_t idx = start_gap; idx < end_gap; ++idx) {
      cs_lnum_t part_idx = interf_tot[idx][1];
      cs_real_t weight = cs_lagr_particles_get_real
                           (p_set, part_idx, CS_LAGR_STAT_WEIGHT);

      if (weight <= 0.) {
        cs_lagr_particles_set_lnum(p_set, part_idx,
                                   CS_LAGR_AGGLO_CLASS_ID, 0);
        cs_lagr_particles_set_real(p_set, part_idx, CS_LAGR_STAT_WEIGHT, 0);
        cs_lagr_particles_set_lnum(p_set, part_idx,
                                   CS_LAGR_P_FLAG, CS_LAGR_PART_TO_DELETE);

      }
    }
  }

  BFT_FREE(cls_gaps);
  BFT_FREE(interf_tot);

  p_set->n_particles += newpart;

  return ret_val;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
