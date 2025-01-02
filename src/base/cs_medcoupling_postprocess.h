#ifndef __CS_MEDCOUPLING_POSTPROCESS_HXX__
#define __CS_MEDCOUPLING_POSTPROCESS_HXX__

/*============================================================================
 * Postprocessing utilities based on MEDCoupling functions
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _medcoupling_slice_t cs_medcoupling_slice_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Get pointer to a slice based on id
 *
 * return pointer to slice. Raises an error if index is out of
 * bounds.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*
 * Get pointer to slice based on name, raises an error
 * if not found.
 *
 * return pointer to slice, raises error if not found.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Get pointer to slice based on name. Returns NULL if
 * not found.
 *
 * return pointer to slice, NULL if not found.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_name_try(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Add a slice based on a plane.
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_postprocess_add_plane_slice(const char  *name,
                                           const char  *selection_criteria,
                                           const cs_real_t  origin[],
                                           const cs_real_t  normal[],
                                           const cs_real_t  length1,
                                           const cs_real_t  length2);

/*----------------------------------------------------------------------------*/
/*
 * Add a slice based on a disc
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_postprocess_add_disc_slice(const char  *name,
                                          const char  *selection_criteria,
                                          const cs_real_t  origin[],
                                          const cs_real_t  normal[],
                                          const cs_real_t  radius,
                                          const int        n_sectors);

/*----------------------------------------------------------------------------*/
/*
 * Add a slice based on an annulus
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_postprocess_add_annulus_slice(const char  *name,
                                             const char  *selection_criteria,
                                             const cs_real_t  origin[],
                                             const cs_real_t  normal[],
                                             const cs_real_t  radius1,
                                             const cs_real_t  radius2,
                                             const int        n_sectors);

/*----------------------------------------------------------------------------*/
/*
 * Get number cells that may be intersected by the slice.
 *
 * return Number of elements
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_medcoupling_slice_get_n_elts(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Get list of ids of the elements which may be intersected.
 *
 * return Pointer to list of ids (cs_lnum_t *). Do not deallocate!
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_medcoupling_slice_get_elt_ids(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Get list of intersection surfaces for each cell intersected.
 *
 * return Pointer to list of intersection surfaces (cs_real_t *)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_slice_get_surfaces(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Get total intersection surface between a slice and volume mesh
 *
 * return Value of total intersection surface
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_get_total_surface(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Activate postprocessing of intersected cells
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_slice_activate_postprocess
(
  const char *name
);

/*----------------------------------------------------------------------------*/
/*
 * Compute integral of a scalar over a slice.
 *
 * return Global integrated value. A cs_parall_sum is used.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_integral(const char       *name,
                                     const cs_real_t  *scalar);

/*----------------------------------------------------------------------------*/
/*
 * Compute mean value of a scalar over a slice.
 *
 * return Global integrated value. A cs_parall_sum is used.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_mean(const char       *name,
                                 const cs_real_t  *scalar);

/*----------------------------------------------------------------------------*/
/*
 * Compute integral of a scalar over a slice using a scalar and/or
 * vectorial weights. If NULL is provided for both weights,
 * the non-weighted function is called.
 *
 * return Computed integral value over entire slice (parallel)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_integral_weighted(const char        *name,
                                              const cs_real_t   *scalar,
                                              const cs_real_t   *weight_s,
                                              const cs_real_3_t *weight_v);

/*----------------------------------------------------------------------------*/
/*
 * Compute mean of a scalar over a slice using a scalar and/or vectorial
 * weights. If NULL is provided for both weights, the non-weighted
 * function is called.
 *
 * return Computed mean value over entire slice (parallel)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_mean_weighted(const char        *name,
                                          const cs_real_t   *scalar,
                                          const cs_real_t   *weight_s,
                                          const cs_real_3_t *weight_v);

/*----------------------------------------------------------------------------*/
/*
 * Destroy all slices
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_slice_destroy_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MEDCOUPLING_POSTPROCESS_HXX__ */
