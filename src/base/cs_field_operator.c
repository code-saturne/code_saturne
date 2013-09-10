/*============================================================================
 * Field based algebraic operators.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_gradient_perio.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field_operator.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_field_operator.c
        Field based algebraic operators.
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void cs_f_field_gradient_scalar(int                    f_id,
                                int                    use_previous_t,
                                int                    imrgra,
                                int                    inc,
                                int                    recompute_cocg,
                                int                    n_r_sweeps,
                                int                    verbosity,
                                int                    clip_mode,
                                double                 epsilon,
                                double                 extrap,
                                double                 clip_coeff,
                                cs_real_3_t  *restrict grad);

void cs_f_field_gradient_potential(int                    f_id,
                                   int                    use_previous_t,
                                   int                    imrgra,
                                   int                    inc,
                                   int                    recompute_cocg,
                                   int                    n_r_sweeps,
                                   int                    hyd_p_flag,
                                   int                    verbosity,
                                   int                    clip_mode,
                                   double                 epsilon,
                                   double                 extrap,
                                   double                 clip_coeff,
                                   cs_real_3_t            f_ext[],
                                   cs_real_3_t  *restrict grad);

void cs_f_field_gradient_vector(int                     f_id,
                                int                     use_previous_t,
                                int                     imrgra,
                                int                     inc,
                                int                     n_r_sweeps,
                                int                     verbosity,
                                int                     clip_mode,
                                double                  epsilon,
                                double                  clip_coeff,
                                cs_real_33_t  *restrict grad);

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f_id           <-- field id
 *   use_previous_t <-- should we use values from the previous time step ?
 *   imrgra         <-- gradient reconstruction mode
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   n_r_sweeps     <-- if > 1, number of reconstruction sweeps
 *   verbosity      <-- verbosity level
 *   clip_mode      <-- clipping mode
 *   epsilon        <-- precision for iterative gradient calculation
 *   extrap         <-- boundary gradient extrapolation coefficient
 *   clip_coeff     <-- clipping coefficient
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_f_field_gradient_scalar(int                    f_id,
                                int                    use_previous_t,
                                int                    imrgra,
                                int                    inc,
                                int                    recompute_cocg,
                                int                    n_r_sweeps,
                                int                    verbosity,
                                int                    clip_mode,
                                double                 epsilon,
                                double                 extrap,
                                double                 clip_coeff,
                                cs_real_3_t  *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  bool _use_previous_t = use_previous_t ? true : false;
  bool _recompute_cocg = recompute_cocg ? true : false;

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  cs_field_gradient_scalar(f,
                           _use_previous_t,
                           gradient_type,
                           halo_type,
                           inc,
                           _recompute_cocg,
                           n_r_sweeps,
                           verbosity,
                           clip_mode,
                           epsilon,
                           extrap,
                           clip_coeff,
                           grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f_id           <-- field id
 *   use_previous_t <-- should we use values from the previous time step ?
 *   imrgra         <-- gradient reconstruction mode
 *   halo_type      <-- halo type
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   n_r_sweeps     <-- if > 1, number of reconstruction sweeps
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   verbosity      <-- verbosity level
 *   clip_mode      <-- clipping mode
 *   epsilon        <-- precision for iterative gradient calculation
 *   extrap         <-- boundary gradient extrapolation coefficient
 *   clip_coeff     <-- clipping coefficient
 *   f_ext          <-- exterior force generating the hydrostatic pressure
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_f_field_gradient_potential(int                    f_id,
                                   int                    use_previous_t,
                                   int                    imrgra,
                                   int                    inc,
                                   int                    recompute_cocg,
                                   int                    n_r_sweeps,
                                   int                    hyd_p_flag,
                                   int                    verbosity,
                                   int                    clip_mode,
                                   double                 epsilon,
                                   double                 extrap,
                                   double                 clip_coeff,
                                   cs_real_3_t            f_ext[],
                                   cs_real_3_t  *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  bool _use_previous_t = use_previous_t ? true : false;
  bool _recompute_cocg = recompute_cocg ? true : false;

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  cs_field_gradient_potential(f,
                              _use_previous_t,
                              gradient_type,
                              halo_type,
                              inc,
                              _recompute_cocg,
                              n_r_sweeps,
                              hyd_p_flag,
                              verbosity,
                              clip_mode,
                              epsilon,
                              extrap,
                              clip_coeff,
                              f_ext,
                              grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f_id           <-- field id
 *   use_previous_t <-- should we use values from the previous time step ?
 *   imrgra         <-- gradient reconstruction mode
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   n_r_sweeps     <-- if > 1, number of reconstruction sweeps
 *   verbosity      <-- verbosity level
 *   clip_mode      <-- clipping mode
 *   epsilon        <-- precision for iterative gradient calculation
 *   extrap         <-- boundary gradient extrapolation coefficient
 *   clip_coeff     <-- clipping coefficient
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_f_field_gradient_vector(int                     f_id,
                                int                     use_previous_t,
                                int                     imrgra,
                                int                     inc,
                                int                     n_r_sweeps,
                                int                     verbosity,
                                int                     clip_mode,
                                double                  epsilon,
                                double                  clip_coeff,
                                cs_real_33_t  *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  bool _use_previous_t = use_previous_t ? true : false;

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  cs_field_gradient_vector(f,
                           _use_previous_t,
                           gradient_type,
                           halo_type,
                           inc,
                           n_r_sweeps,
                           verbosity,
                           clip_mode,
                           epsilon,
                           clip_coeff,
                           grad);
}

/*! \endcond (end ignore by Doxygen) */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       extrap          boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff      clipping coefficient
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void cs_field_gradient_scalar(const cs_field_t          *f,
                              bool                       use_previous_t,
                              cs_gradient_type_t         gradient_type,
                              cs_halo_type_t             halo_type,
                              int                        inc,
                              bool                       recompute_cocg,
                              int                        n_r_sweeps,
                              int                        verbosity,
                              int                        clip_mode,
                              double                     epsilon,
                              double                     extrap,
                              double                     clip_coeff,
                              cs_real_3_t      *restrict grad)
{
  int tr_dim = 0;

  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  cs_gradient_perio_init_rij(f, &tr_dim, grad);

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     recompute_cocg,
                     n_r_sweeps,
                     tr_dim,
                     0, /* hyd_p_flag */
                     verbosity,
                     clip_mode,
                     epsilon,
                     extrap,
                     clip_coeff,
                     NULL, /* f_ext */
                     f->bc_coeffs->a,
                     f->bc_coeffs->b,
                     var,
                     NULL,
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       hyd_p_flag      flag for hydrostatic pressure
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       extrap          boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       f_ext           exterior force generating
 *                                  the hydrostatic pressure
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void cs_field_gradient_potential(const cs_field_t          *f,
                                 bool                       use_previous_t,
                                 cs_gradient_type_t         gradient_type,
                                 cs_halo_type_t             halo_type,
                                 int                        inc,
                                 bool                       recompute_cocg,
                                 int                        n_r_sweeps,
                                 int                        hyd_p_flag,
                                 int                        verbosity,
                                 int                        clip_mode,
                                 double                     epsilon,
                                 double                     extrap,
                                 double                     clip_coeff,
                                 cs_real_3_t                f_ext[],
                                 cs_real_3_t      *restrict grad)
{
  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     recompute_cocg,
                     n_r_sweeps,
                     0, /* tr_dim */
                     hyd_p_flag,
                     verbosity,
                     clip_mode,
                     epsilon,
                     extrap,
                     clip_coeff,
                     f_ext,
                     f->bc_coeffs->a,
                     f->bc_coeffs->b,
                     var,
                     NULL,
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       gradient_type   gradient type
 * \param[in]       halo_type       halo type
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       recompute_cocg  should COCG FV quantities be recomputed ?
 * \param[in]       n_r_sweeps      if > 1, number of reconstruction sweeps
 * \param[in]       verbosity       verbosity level
 * \param[in]       clip_mode       clipping mode
 * \param[in]       epsilon         precision for iterative gradient calculation
 * \param[in]       extrap          boundary gradient extrapolation coefficient
 * \param[in]       clip_coeff      clipping coefficient
 * \param[in]       f_ext           exterior force generating
 *                                  the hydrostatic pressure
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void cs_field_gradient_vector(const cs_field_t          *f,
                              bool                       use_previous_t,
                              cs_gradient_type_t         gradient_type,
                              cs_halo_type_t             halo_type,
                              int                        inc,
                              int                        n_r_sweeps,
                              int                        verbosity,
                              int                        clip_mode,
                              double                     epsilon,
                              double                     clip_coeff,
                              cs_real_33_t     *restrict grad)
{
  cs_real_3_t *var = (use_previous_t) ? (cs_real_3_t *)(f->val_pre)
                                      : (cs_real_3_t *)(f->val);

  cs_gradient_vector(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     n_r_sweeps,
                     verbosity,
                     clip_mode,
                     epsilon,
                     clip_coeff,
                     (const cs_real_3_t *)(f->bc_coeffs->a),
                     (const cs_real_33_t *)(f->bc_coeffs->b),
                     var,
                     grad);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
