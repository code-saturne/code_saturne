#ifndef __CS_GRADIENT_H__
#define __CS_GRADIENT_H__

/*============================================================================
 * Gradient reconstruction.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_internal_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

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
  CS_GRADIENT_GREEN_LSQ          /*!< Green-Gauss reconstruction with
                                   least squares gradient face values */

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

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the steady balance due to porous modelling for the pressure
 * gradient
 *----------------------------------------------------------------------------*/

void CS_PROCF (grdpor, GRDPOR)
(
 const cs_int_t   *const inc          /* <-- 0 or 1: increment or not         */
);

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdcel, CGDCEL)
(
 const cs_int_t   *const f_id,        /* <-- field id, or -1                  */
 const cs_int_t   *const imrgra,      /* <-- gradient computation mode        */
 const cs_int_t   *const inc,         /* <-- 0 or 1: increment or not         */
 const cs_int_t   *const iccocg,      /* <-- 1 or 0: recompute COCG or not    */
 const cs_int_t   *const n_r_sweeps,  /* <-- >1: with reconstruction          */
 const cs_int_t   *const idimtr,      /* <-- 0, 1, 2: scalar, vector, tensor
                                             in case of rotation              */
 const cs_int_t   *const iphydp,      /* <-- use hydrosatatic pressure        */
 const cs_int_t   *const ipond,       /* <-- >0: weighted gradient computation*/
 const cs_int_t   *const iwarnp,      /* <-- verbosity level                  */
 const cs_int_t   *const imligp,      /* <-- type of clipping                 */
 const cs_real_t  *const epsrgp,      /* <-- precision for iterative gradient
                                             calculation                      */
 const cs_real_t  *const extrap,      /* <-- extrapolate gradient at boundary */
 const cs_real_t  *const climgp,      /* <-- clipping coefficient             */
       cs_real_3_t       f_ext[],     /* <-- exterior force generating the
                                             hydrostatic pressure             */
 const cs_real_t         coefap[],    /* <-- boundary condition term          */
 const cs_real_t         coefbp[],    /* <-- boundary condition term          */
       cs_real_t         pvar[],      /* <-- gradient's base variable         */
       cs_real_t         ktvar[],     /* <-- gradient coefficient variable    */
       cs_real_3_t       grad[]       /* <-> gradient                         */
);

/*----------------------------------------------------------------------------
 * Compute cell gradient of vector field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdvec, CGDVEC)
(
 const cs_int_t         *const f_id,      /* <-- field id, or -1              */
 const cs_int_t         *const imrgra,    /* <-- gradient computation mode    */
 const cs_int_t         *const inc,       /* <-- 0 or 1: increment or not     */
 const cs_int_t         *const n_r_sweeps,/* <-- >1: with reconstruction      */
 const cs_int_t         *const iwarnp,    /* <-- verbosity level              */
 const cs_int_t         *const imligp,    /* <-- type of clipping             */
 const cs_real_t        *const epsrgp,    /* <-- precision for iterative
                                                 gradient calculation         */
 const cs_real_t        *const climgp,    /* <-- clipping coefficient         */
 const cs_real_3_t             coefav[],  /* <-- boundary condition term      */
 const cs_real_33_t            coefbv[],  /* <-- boundary condition term      */
       cs_real_3_t             pvar[],    /* <-- gradient's base variable     */
       cs_real_33_t            gradv[]    /* <-> gradient of the variable
                                                 (du_i/dx_j : gradv[][i][j])  */
);

/*----------------------------------------------------------------------------
 * Compute cell gradient of tensor field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdts, CGDTS)
(
 const cs_int_t         *const f_id,
 const cs_int_t         *const imrgra,    /* <-- gradient computation mode    */
 const cs_int_t         *const inc,       /* <-- 0 or 1: increment or not     */
 const cs_int_t         *const n_r_sweeps,    /* <-- >1: with reconstruction  */
 const cs_int_t         *const iwarnp,    /* <-- verbosity level              */
 const cs_int_t         *const imligp,    /* <-- type of clipping             */
 const cs_real_t        *const epsrgp,    /* <-- precision for iterative
                                                 gradient calculation         */
 const cs_real_t        *const climgp,    /* <-- clipping coefficient         */
 const cs_real_6_t             coefav[],  /* <-- boundary condition term      */
 const cs_real_66_t            coefbv[],  /* <-- boundary condition term      */

       cs_real_6_t             pvar[],    /* <-- gradient's base variable     */
       cs_real_63_t            grad[]    /* <-> gradient of the variable
                                                 (du_i/dx_j : gradv[][i][j])  */
);

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
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 *         tensor field.
 *
 * \param[in]       var_name       variable name
 * \param[in]       gradient_type  gradient type
 * \param[in]       halo_type      halo type
 * \param[in]       inc            if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps     if > 1, number of reconstruction sweeps
 *                                 (only used by CS_GRADIENT_ITER)
 * \param[in]       tr_dim         2 for tensor with periodicity of rotation,
 *                                 0 otherwise
 * \param[in]       hyd_p_flag     flag for hydrostatic pressure
 * \param[in]       w_stride       stride for weighting coefficient
 * \param[in]       verbosity      verbosity level
 * \param[in]       clip_mode      clipping mode
 * \param[in]       epsilon        precision for iterative gradient calculation
 * \param[in]       extrap         boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff     clipping coefficient
 * \param[in]       f_ext          exterior force generating the
 *                                 hydrostatic pressure
 * \param[in]       bc_coeff_a     boundary condition term a
 * \param[in]       bc_coeff_b     boundary condition term b
 * \param[in, out]  var            gradient's base variable
 * \param[in, out]  c_weight       cell variable weight, or NULL
 * \param[in]       cpl            associated internal coupling, or NULL
 * \param[out]      grad           gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar(const char                    *var_name,
                   cs_gradient_type_t             gradient_type,
                   cs_halo_type_t                 halo_type,
                   int                            inc,
                   bool                           recompute_cocg,
                   int                            n_r_sweeps,
                   int                            tr_dim,
                   int                            hyd_p_flag,
                   int                            w_stride,
                   int                            verbosity,
                   cs_gradient_limit_t            clip_mode,
                   double                         epsilon,
                   double                         extrap,
                   double                         clip_coeff,
                   cs_real_3_t                    f_ext[],
                   const cs_real_t                bc_coeff_a[],
                   const cs_real_t                bc_coeff_b[],
                   cs_real_t                      var[restrict],
                   cs_real_t            *restrict c_weight,
                   const cs_internal_coupling_t  *cpl,
                   cs_real_t                      grad[restrict][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of vector field.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[in, out]  c_weight        cell variable weight, or NULL
 * \param[in]       cpl             associated internal coupling, or NULL
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
                   const cs_real_t                bc_coeff_a[][3],
                   const cs_real_t                bc_coeff_b[][3][3],
                   cs_real_t                      var[restrict][3],
                   cs_real_t        *restrict     c_weight,
                   const cs_internal_coupling_t  *cpl,
                   cs_real_t                      gradv[restrict][3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor(const char                *var_name,
                   cs_gradient_type_t         gradient_type,
                   cs_halo_type_t             halo_type,
                   int                        inc,
                   int                        n_r_sweeps,
                   int                        verbosity,
                   cs_gradient_limit_t        clip_mode,
                   double                     epsilon,
                   double                     clip_coeff,
                   const cs_real_6_t          bc_coeff_a[],
                   const cs_real_66_t         bc_coeff_b[],
                   cs_real_6_t      *restrict var,
                   cs_real_63_t     *restrict grad);

/*----------------------------------------------------------------------------*/
/*!
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
 * \param[in]   recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]   n_r_sweeps      if > 1, number of reconstruction sweeps
 *                              (only used by CS_GRADIENT_ITER)
 * \param[in]   tr_dim          2 for tensor with periodicity of rotation,
 *                              0 otherwise
 * \param[in]   hyd_p_flag      flag for hydrostatic pressure
 * \param[in]   w_stride        stride for weighting coefficient
 * \param[in]   verbosity       verbosity level
 * \param[in]   clip_mode       clipping mode
 * \param[in]   epsilon         precision for iterative gradient calculation
 * \param[in]   extrap          boundary gradient extrapolation coefficient
 * \param[in]   clip_coeff      clipping coefficient
 * \param[in]   f_ext           exterior force generating the
 *                              hydrostatic pressure
 * \param[in]   bc_coeff_a      boundary condition term a
 * \param[in]   bc_coeff_b      boundary condition term b
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   cpl             associated internal coupling, or NULL
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_synced_input(const char                 *var_name,
                                cs_gradient_type_t          gradient_type,
                                cs_halo_type_t              halo_type,
                                int                         inc,
                                bool                        recompute_cocg,
                                int                         n_r_sweeps,
                                int                         tr_dim,
                                int                         hyd_p_flag,
                                int                         w_stride,
                                int                         verbosity,
                                cs_gradient_limit_t         clip_mode,
                                double                      epsilon,
                                double                      extrap,
                                double                      clip_coeff,
                                cs_real_t                   f_ext[][3],
                                const cs_real_t             bc_coeff_a[],
                                const cs_real_t             bc_coeff_b[],
                                const cs_real_t             var[restrict],
                                const cs_real_t             c_weight[restrict],
                                const cs_internal_coupling_t  *cpl,
                                cs_real_t                   grad[restrict][3]);

/*----------------------------------------------------------------------------*/
/*!
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
 *                              (only used by CS_GRADIENT_ITER)
 * \param[in]   verbosity       verbosity level
 * \param[in]   clip_mode       clipping mode
 * \param[in]   epsilon         precision for iterative gradient calculation
 * \param[in]   clip_coeff      clipping coefficient
 * \param[in]   bc_coeff_a      boundary condition term a
 * \param[in]   bc_coeff_b      boundary condition term b
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   cpl             associated internal coupling, or NULL
 * \param[out]  gradv           gradient
                                (\f$ \der{u_i}{x_j} \f$ is gradv[][i][j])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_synced_input(const char                *var_name,
                                cs_gradient_type_t         gradient_type,
                                cs_halo_type_t             halo_type,
                                int                        inc,
                                int                        n_r_sweeps,
                                int                        verbosity,
                                cs_gradient_limit_t        clip_mode,
                                double                     epsilon,
                                double                     clip_coeff,
                                const cs_real_t            bc_coeff_a[][3],
                                const cs_real_t            bc_coeff_b[][3][3],
                                const cs_real_t            var[restrict][3],
                                const cs_real_t            c_weight[restrict],
                                const cs_internal_coupling_t  *cpl,
                                cs_real_33_t     *restrict gradv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor.
 *
 * \param[in]       var_name        variable name
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 *                                  (only used by CS_GRADIENT_ITER)
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       bc_coeff_a      boundary condition term a
 * \param[in]       bc_coeff_b      boundary condition term b
 * \param[in, out]  var             gradient's base variable
 * \param[out]      grad            gradient
                                    (\f$ \der{t_ij}{x_k} \f$ is grad[][ij][k])
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_synced_input(const char                *var_name,
                                cs_gradient_type_t         gradient_type,
                                cs_halo_type_t             halo_type,
                                int                        inc,
                                int                        n_r_sweeps,
                                int                        verbosity,
                                cs_gradient_limit_t        clip_mode,
                                double                     epsilon,
                                double                     clip_coeff,
                                const cs_real_t            bc_coeff_a[][6],
                                const cs_real_t            bc_coeff_b[][6][6],
                                const cs_real_t            var[restrict][6],
                                cs_real_63_t     *restrict grad);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a scalar field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_scalar_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_real_t              bc_coeff_a[],
                        const cs_real_t              bc_coeff_b[],
                        const cs_real_t              var[],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a vector field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_vector_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_real_t              bc_coeff_a[][3],
                        const cs_real_t              bc_coeff_b[][3][3],
                        const cs_real_t              var[restrict][3],
                        const cs_real_t              c_weight[],
                        cs_real_t                    grad[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient of a tensor field at a given cell
 *         using least-squares reconstruction.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   c_id            cell id
 * \param[in]   halo_type       halo type
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   var             gradient's base variable
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[out]  grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_tensor_cell(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *fvq,
                        cs_lnum_t                    c_id,
                        cs_halo_type_t               halo_type,
                        const cs_real_t              bc_coeff_a[][6],
                        const cs_real_t              bc_coeff_b[][6][6],
                        const cs_real_t              var[restrict][6],
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

END_C_DECLS

#endif /* __CS_GRADIENT__ */
