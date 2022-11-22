#ifndef __CS_FUNCTION_DEFAULT_H__
#define __CS_FUNCTION_DEFAULT_H__

/*============================================================================
 * Base predefined function objects.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_function.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Predefined function objects */

typedef enum {

  CS_FUNCTION_CELL_RANK_ID,      /*!< cell MPI rank id (integer) */
  CS_FUNCTION_B_FACE_RANK_ID     /*!< boundary face MPI rank id (integer) */

} cs_function_predefined_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define functions based on code_saturne case setup.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_default_define(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access a function for evaluation of element's MPI rank id.
 *
 * \param[in]   location_id  base associated mesh location id
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_mpi_rank_id(cs_mesh_location_type_t  location_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access a function for evaluation of mesh element's
 *        refinement generation (i.e. level).
 *
 * \param[in]   location_id  base associated mesh location id
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_refinement_generation(cs_mesh_location_type_t  location_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of normal boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress_normal(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of tangential boundary stress.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_stress_tangential(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function object for computation of boundary thermal flux.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_thermal_flux(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function for computation of cell Q criterion.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_q_criterion(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract optional boundary face class of element zone id.
 *
 * For boundary faces, if no face classes have been defined by
 * \ref cs_boundary_zone_face_class_id the highest boundary face zone id is
 *
 * For cells, the highest cell volume zone id is used.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to field
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_class_or_zone_id(int               location_id,
                             cs_lnum_t         n_elts,
                             const cs_lnum_t  *elt_ids,
                             void             *input,
                             void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute non-reconstructed cell-based field values at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to field
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_field_boundary_nr(int               location_id,
                              cs_lnum_t         n_elts,
                              const cs_lnum_t  *elt_ids,
                              void             *input,
                              void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute stress at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_stress(int               location_id,
                            cs_lnum_t         n_elts,
                            const cs_lnum_t  *elt_ids,
                            void             *input,
                            void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute normal stress at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_stress_normal(int               location_id,
                                   cs_lnum_t         n_elts,
                                   const cs_lnum_t  *elt_ids,
                                   void             *input,
                                   void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute tangential stress at boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_stress_tangential(int               location_id,
                                       cs_lnum_t         n_elts,
                                       const cs_lnum_t  *elt_ids,
                                       void             *input,
                                       void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define function for computation of boundary layer Nusselt.
 *
 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_boundary_nusselt(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute thermal flux at boundary (in \f$ W\,m^{-2} \f$),
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_thermal_flux(int               location_id,
                                  cs_lnum_t         n_elts,
                                  const cs_lnum_t  *elt_ids,
                                  void             *input,
                                  void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute local Nusselt number near boundary.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_boundary_nusselt(int               location_id,
                             cs_lnum_t         n_elts,
                             const cs_lnum_t  *elt_ids,
                             void             *input,
                             void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Q-criterion from Hunt et. al over each cell of a specified
 *        volume region.
 *
 * \f[
 *    Q = \tens{\Omega}:\tens{\Omega} -
 *    \deviator{ \left(\tens{S} \right)}:\deviator{ \left(\tens{S} \right)}
 * \f]
 * where \f$\tens{\Omega}\f$ is the vorticity tensor and
 * \f$\deviator{ \left(\tens{S} \right)}\f$ the deviatoric of the rate of strain
 * tensor.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_q_criterion(int               location_id,
                        cs_lnum_t         n_elts,
                        const cs_lnum_t  *elt_ids,
                        void             *input,
                        void             *vals);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FUNCTION_DEFAULT_H__ */
