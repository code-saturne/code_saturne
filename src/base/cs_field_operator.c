/*============================================================================
 * Field based algebraic operators.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include "cs_parameters.h"
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
                                cs_real_3_t  *restrict grad);

void cs_f_field_gradient_potential(int                    f_id,
                                   int                    use_previous_t,
                                   int                    imrgra,
                                   int                    inc,
                                   int                    recompute_cocg,
                                   int                    hyd_p_flag,
                                   cs_real_3_t            f_ext[],
                                   cs_real_3_t  *restrict grad);

void cs_f_field_gradient_vector(int                     f_id,
                                int                     use_previous_t,
                                int                     imrgra,
                                int                     inc,
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
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_f_field_gradient_scalar(int                    f_id,
                                int                    use_previous_t,
                                int                    imrgra,
                                int                    inc,
                                int                    recompute_cocg,
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
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   f_ext          <-- exterior force generating the hydrostatic pressure
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_f_field_gradient_potential(int                    f_id,
                                   int                    use_previous_t,
                                   int                    imrgra,
                                   int                    inc,
                                   int                    recompute_cocg,
                                   int                    hyd_p_flag,
                                   cs_real_3_t            f_ext[],
                                   cs_real_3_t  *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  bool _use_previous_t = use_previous_t ? true : false;
  bool _recompute_cocg = recompute_cocg ? true : false;

  if (imrgra >= 10)
    imrgra = 10;
  else if (imrgra > 0)
    imrgra = 0;
  if (imrgra < 0)
    imrgra = -imrgra;

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
                              hyd_p_flag,
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
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_f_field_gradient_vector(int                     f_id,
                                int                     use_previous_t,
                                int                     imrgra,
                                int                     inc,
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
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void cs_field_gradient_scalar(const cs_field_t          *f,
                              bool                       use_previous_t,
                              cs_gradient_type_t         gradient_type,
                              cs_halo_type_t             halo_type,
                              int                        inc,
                              bool                       recompute_cocg,
                              cs_real_3_t      *restrict grad)
{
  int tr_dim = 0;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  // Get the calculation option from the field
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  cs_gradient_perio_init_rij(f, &tr_dim, grad);

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     recompute_cocg,
                     var_cal_opt.nswrgr,
                     tr_dim,
                     0, /* hyd_p_flag */
                     var_cal_opt.iwarni,
                     var_cal_opt.imligr,
                     var_cal_opt.epsrgr,
                     var_cal_opt.extrag,
                     var_cal_opt.climgr,
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
 * \param[in]       hyd_p_flag      flag for hydrostatic pressure
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
                                 int                        hyd_p_flag,
                                 cs_real_3_t                f_ext[],
                                 cs_real_3_t      *restrict grad)
{
  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  // Get the calculation option from the field
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     recompute_cocg,
                     var_cal_opt.nswrgr,
                     0, /* tr_dim */
                     hyd_p_flag,
                     var_cal_opt.iwarni,
                     var_cal_opt.imligr,
                     var_cal_opt.epsrgr,
                     var_cal_opt.extrag,
                     var_cal_opt.climgr,
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
                              cs_real_33_t     *restrict grad)
{
  cs_real_3_t *var;

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  // Get the calculation option from the field
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  if (f->interleaved)
    var = (use_previous_t) ? (cs_real_3_t *)(f->val_pre)
                           : (cs_real_3_t *)(f->val);
  else {
    const int dim = f->dim;
    const cs_real_t *s =  (use_previous_t) ? f->val_pre : f->val;
    const cs_lnum_t *n_loc_elts
      = cs_mesh_location_get_n_elts(f->location_id);
    const cs_lnum_t _n_loc_elts = n_loc_elts[2];
    BFT_MALLOC(var, _n_loc_elts, cs_real_3_t);
    for (cs_lnum_t i = 0; i < _n_loc_elts; i++) {
      for (int j = 0; j < dim; j++)
        var[i][j] = s[j*_n_loc_elts + i];
    }
  }

  cs_gradient_vector(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     var_cal_opt.nswrgr,
                     var_cal_opt.iwarni,
                     var_cal_opt.imligr,
                     var_cal_opt.epsrgr,
                     var_cal_opt.climgr,
                     (const cs_real_3_t *)(f->bc_coeffs->a),
                     (const cs_real_33_t *)(f->bc_coeffs->b),
                     var,
                     grad);

  if (! f->interleaved)
    BFT_FREE(var);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
