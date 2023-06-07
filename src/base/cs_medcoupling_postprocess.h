#ifndef __CS_MEDCOUPLING_POSTPROCESS_HXX__
#define __CS_MEDCOUPLING_POSTPROCESS_HXX__

/*============================================================================
 * Postprocessing utilities based on MEDCoupling functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_defs.h"

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
/*!
 * \brief Get pointer to a slice based on id
 *
 * \param[in] id index of slice
 *
 * \return pointer to slice. Raises an error if index is out of
 * bounds.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to slice based on name, raises an error
 * if not found.
 *
 * \param[in] name  Name of the slice structure
 *
 * \return pointer to slice, raises error if not found.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to slice based on name. Returns NULL if
 * not found.
 *
 * \param[in] name  Name of the slice structure
 *
 * \return pointer to slice, NULL if not found.
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_name_try(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a slice based on a plane.
 *
 * \param[in] name                Name of the slice
 * \param[in] selection_criteria  Selection criteria for cells to intersect
 * \param[in] origin              Coordinates of origin point of slice
 * \param[in] normal              Normal vector of the slice
 * \param[in] length1             Length along the first axis of the plane
 * \param[in] length2             Length along the second axis of the plane
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
/*!
 * \brief Add a slice based on a disc
 *
 * \param[in] name                Name of the slice
 * \param[in] selection_criteria  Selection criteria for cells to intersect
 * \param[in] origin              Coordinates of origin point of slice
 * \param[in] normal              Normal vector of the slice
 * \param[in] radius              Radius of the disc
 * \param[in] n_sectors           Number of sectors for discretization.
 *                                If negative, default value (36) is used.
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
/*!
 * \brief Add a slice based on an annulus
 *
 * \param[in] name                Name of the slice
 * \param[in] selection_criteria  Selection criteria for cells to intersect
 * \param[in] origin              Coordinates of origin point of slice
 * \param[in] normal              Normal vector of the slice
 * \param[in] radius1             Inner radius of the annulus (hole)
 * \param[in] radius2             Outer radius of the annulus
 * \param[in] n_sectors           Number of sectors for discretization.
 *                                If negative, default value (36) is used.
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
/*!
 * \brief Get number cells that may be intersected by the slice.
 *
 * \param[in] name  Name of the slice
 *
 * \return Number of elements
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_medcoupling_slice_get_n_elts(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get list of ids of the elements which may be intersected.
 *
 * \param[in] name  Name of the slice
 *
 * \return Pointer to list of ids (cs_lnum_t *). Do not deallocate!
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_medcoupling_slice_get_elt_ids(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get list of intersection surfaces for each cell intersected.
 *
 * \param[in] name  Name of the slice
 *
 * \return Pointer to list of intersection surfaces (cs_real_t *)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_slice_get_surfaces(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get total intersection surface between a slice and volume mesh
 *
 * \param[in] name  Name of the slice
 *
 * \return Value of total intersection surface
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_get_total_surface(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute integral of a scalar over a slice.
 *
 * \param[in] name    Name of the slice
 * \param[in] scalar  Array of scalar values (size n_cells)
 *
 * \return Global integrated value. A cs_parall_sum is used.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_integral(const char       *name,
                                     const cs_real_t  *scalar);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute mean value of a scalar over a slice.
 *
 * \param[in] name    Name of the slice
 * \param[in] scalar  Array of scalar values (size n_cells)
 *
 * \return Global integrated value. A cs_parall_sum is used.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_mean(const char       *name,
                                 const cs_real_t  *scalar);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute integral of a scalar over a slice using a scalar and/or
 *        vectorial weights. If NULL is provided for both weights,
 *        the non-weighted function is called.
 *
 * \param[in] name      Name of the slice
 * \param[in] scalar    Array of scalar values (size n_cells)
 * \param[in] weight_s  Scalar weight array (size n_cells)
 * \param[in] weight_v  Vectorial weight array (size n_cells)
 *
 * \return Computed integral value over entire slice (parallel)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_integral_weighted(const char        *name,
                                              const cs_real_t   *scalar,
                                              const cs_real_t   *weight_s,
                                              const cs_real_3_t *weight_v);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute mean of a scalar over a slice using a scalar and/or vectorial
 *        weights. If NULL is provided for both weights, the non-weighted
 *        function is called.
 *
 * \param[in] name      Name of the slice
 * \param[in] scalar    Array of scalar values (size n_cells)
 * \param[in] weight_s  Scalar weight array (size n_cells)
 * \param[in] weight_v  Vectorial weight array (size n_cells)
 *
 * \return Computed mean value over entire slice (parallel)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_medcoupling_slice_scalar_mean_weighted(const char        *name,
                                          const cs_real_t   *scalar,
                                          const cs_real_t   *weight_s,
                                          const cs_real_3_t *weight_v);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all slices
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_slice_destroy_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MEDCOUPLING_POSTPROCESS_HXX__ */
