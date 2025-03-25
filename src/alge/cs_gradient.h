#ifndef __CS_GRADIENT_H__
#define __CS_GRADIENT_H__

/*============================================================================
 * Gradient reconstruction.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_base.h"
#include "base/cs_halo.h"
#include "base/cs_internal_coupling.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Gradient reconstruction method
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_GRADIENT_GREEN_ITER,        /*!< Iterative */
  CS_GRADIENT_LSQ,               /*!< Least-squares */
  CS_GRADIENT_GREEN_LSQ,         /*!< Green-Gauss reconstruction with
                                   least squares gradient face values */
  CS_GRADIENT_GREEN_VTX,         /*!< Green-Gauss with vertex interpolation */

  CS_GRADIENT_GREEN_R            /*!< Green-Gauss times a renormalization
                                   matrix (to be exact on affine functions) */
} cs_gradient_type_t;

/*----------------------------------------------------------------------------
 * Gradient limiter mode
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_GRADIENT_LIMIT_NONE = -1,   /*!< no limiter */
  CS_GRADIENT_LIMIT_CELL = 0,    /*!< cell values extrapolated by cell gradient
                                   to neighbor cells can not exceed the
                                   limiting factor  times the actual value */
  CS_GRADIENT_LIMIT_FACE = 1,    /*!< cell values extrapolated by face gradient
                                   (mean of cell gradients) extrapolated to
                                   neighbor cells can not exceed the limiting
                                   factor times the actual value */

} cs_gradient_limit_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for gradient types */

extern const char *cs_gradient_type_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free saved gradient quantities.
 *
 * This is required when the mesh changes, so that the on-demand computation
 * will be updated.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_free_quantities(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * \param[in]       var_name       variable name
 * \param[in]       gradient_type  gradient type
 * \param[in]       halo_type      halo type
 * \param[in]       inc            if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps     if > 1, number of reconstruction sweeps
 *                                 (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       hyd_p_flag     flag for hydrostatic pressure
 * \param[in]       w_stride       stride for weighting coefficient
 * \param[in]       verbosity      verbosity level
 * \param[in]       clip_mode      clipping mode
 * \param[in]       epsilon        precision for iterative gradient calculation
 * \param[in]       clip_coeff     clipping coefficient
 * \param[in]       f_ext          exterior force generating the
 *                                 hydrostatic pressure
 * \param[in]       bc_coeffs      boundary condition structure
 * \param[in, out]  var            gradient's base variable
 * \param[in, out]  c_weight       cell variable weight, or nullptr
 * \param[in]       cpl            associated internal coupling, or nullptr
 * \param[out]      grad           gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar(const char                    *var_name,
                   cs_gradient_type_t             gradient_type,
                   cs_halo_type_t                 halo_type,
                   int                            inc,
                   int                            n_r_sweeps,
                   int                            hyd_p_flag,
                   int                            w_stride,
                   int                            verbosity,
                   cs_gradient_limit_t            clip_mode,
                   double                         epsilon,
                   double                         clip_coeff,
                   cs_real_3_t                    f_ext[],
                   const cs_field_bc_coeffs_t    *bc_coeffs,
                   cs_real_t                      var[],
                   cs_real_t                     *c_weight,
                   const cs_internal_coupling_t  *cpl,
                   cs_real_t                      grad[][3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute cell gradient of vector field.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_v     boundary condition structure
 * \param[in, out]  var             gradient's base variable
 * \param[in, out]  c_weight        cell variable weight, or nullptr
 * \param[in]       cpl             associated internal coupling, or nullptr
 * \param[out]      gradv           gradient
                                    (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector(const char                    *var_name,
                   cs_gradient_type_t             gradient_type,
                   cs_halo_type_t                 halo_type,
                   int                            inc,
                   int                            n_r_sweeps,
                   int                            verbosity,
                   cs_gradient_limit_t            clip_mode,
                   double                         epsilon,
                   double                         clip_coeff,
                   const cs_field_bc_coeffs_t    *bc_coeffs_v,
                   cs_real_t                      var[][3],
                   cs_real_t                     *c_weight,
                   const cs_internal_coupling_t  *cpl,
                   cs_real_t                      gradv[][3][3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_ts    boundary condition structure
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor(const char                  *var_name,
                   cs_gradient_type_t           gradient_type,
                   cs_halo_type_t               halo_type,
                   int                          inc,
                   int                          n_r_sweeps,
                   int                          verbosity,
                   cs_gradient_limit_t          clip_mode,
                   double                       epsilon,
                   double                       clip_coeff,
                   const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                   cs_real_6_t                 *var,
                   cs_real_63_t                *grad);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * This variant of the \ref cs_gradient_scalar function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]   var_name        variable name
 * \param[in]   gradient_type   gradient type
 * \param[in]   halo_type       halo type
 * \param[in]   inc             if 0, solve on increment; 1 otherwise
 * \param[in]   n_r_sweeps      if > 1, number of reconstruction sweeps
 *                              (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]   hyd_p_flag      flag for hydrostatic pressure
 * \param[in]   w_stride        stride for weighting coefficient
 * \param[in]   verbosity       verbosity level
 * \param[in]   clip_mode       clipping mode
 * \param[in]   epsilon         precision for iterative gradient calculation
 * \param[in]   clip_coeff      clipping coefficient
 * \param[in]   f_ext           exterior force generating the
 *                              hydrostatic pressure
 * \param[in]   bc_coeffs       boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[in]   cpl             associated internal coupling, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_synced_input(const char                    *var_name,
                                cs_gradient_type_t             gradient_type,
                                cs_halo_type_t                 halo_type,
                                int                            inc,
                                int                            n_r_sweeps,
                                int                            hyd_p_flag,
                                int                            w_stride,
                                int                            verbosity,
                                cs_gradient_limit_t            clip_mode,
                                double                         epsilon,
                                double                         clip_coeff,
                                cs_real_t                      f_ext[][3],
                                const cs_field_bc_coeffs_t    *bc_coeffs,
                                const cs_real_t                var[],
                                const cs_real_t               *c_weight,
                                const cs_internal_coupling_t  *cpl,
                                cs_real_t                      grad[][3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute cell gradient of vector field.
 *
 * This variant of the \ref cs_gradient_vector function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]   var_name        variable name
 * \param[in]   gradient_type   gradient type
 * \param[in]   halo_type       halo type
 * \param[in]   inc             if 0, solve on increment; 1 otherwise
 * \param[in]   n_r_sweeps      if > 1, number of reconstruction sweeps
 *                              (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]   verbosity       verbosity level
 * \param[in]   clip_mode       clipping mode
 * \param[in]   epsilon         precision for iterative gradient calculation
 * \param[in]   clip_coeff      clipping coefficient
 * \param[in]   bc_coeffs_v     boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
                                (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_synced_input(const char                    *var_name,
                                cs_gradient_type_t             gradient_type,
                                cs_halo_type_t                 halo_type,
                                int                            inc,
                                int                            n_r_sweeps,
                                int                            verbosity,
                                cs_gradient_limit_t            clip_mode,
                                double                         epsilon,
                                double                         clip_coeff,
                                const cs_field_bc_coeffs_t    *bc_coeffs_v,
                                const cs_real_t                var[][3],
                                const cs_real_t                val_f[][3],
                                const cs_real_t                c_weight[],
                                cs_real_t                      grad[][3][3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute cell gradient of tensor.
 *
 * This variant of the \ref cs_gradient_tensor function assumes ghost cell
 * values for input arrays (var and optionally c_weight)
 * have already been synchronized.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_GREEN_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeffs_ts    boundary condition structure
 * \param[in, out]  var             gradient's base variable:
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_synced_input(const char                  *var_name,
                                cs_gradient_type_t           gradient_type,
                                cs_halo_type_t               halo_type,
                                int                          inc,
                                int                          n_r_sweeps,
                                int                          verbosity,
                                cs_gradient_limit_t          clip_mode,
                                double                       epsilon,
                                double                       clip_coeff,
                                const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                                const cs_real_t              var[][6],
                                const cs_real_t              val_f[][6],
                                cs_real_63_t                *grad);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute the gradient of a scalar field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs       boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_field_bc_coeffs_t  *bc_coeffs,
                        const cs_real_t              var[],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute the gradient of a vector field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs       boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_field_bc_coeffs_t  *bc_coeffs_v,
                        const cs_real_t              var[][3],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[3][3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute the gradient of a tensor field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * When boundary conditions are provided, both the bc_coeff_a and bc_coeff_b
 * arrays must be given. If boundary values are known, bc_coeff_a
 * can point to the boundary values array, and bc_coeff_b set to nullptr.
 * If bc_coeff_a is nullptr, bc_coeff_b is ignored.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeffs_ts    boundary condition structure
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or nullptr
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                        const cs_real_t              var[][6],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[6][3]);

/*----------------------------------------------------------------------------
 * Determine gradient type by Fortran "imrgra" value
 *
 * parameters:
 *   imrgra         <-- Fortran gradient option
 *   gradient_type  --> gradient type
 *   halo_type      --> halo type
 *----------------------------------------------------------------------------*/

void
cs_gradient_type_by_imrgra(int                  imrgra,
                           cs_gradient_type_t  *gradient_type,
                           cs_halo_type_t      *halo_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief compute the steady balance due to porous modelling for the pressure
 *        gradient.
 *
 * \param[in]  inc  if 0, solve on increment; 1 otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_porosity_balance(int inc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GRADIENT__ */
