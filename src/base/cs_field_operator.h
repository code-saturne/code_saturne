#ifndef __CS_FIELD_OPERATOR_H__
#define __CS_FIELD_OPERATOR_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"
#include "cs_gradient.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   gradient_type  <-- gradient type
 *   halo_type      <-- halo type
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   use_previous_t <-- should we use values from the previous time step ?
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   n_r_sweeps     <-- if > 1, number of reconstruction sweeps
 *   verbosity      <-- verbosity level
 *   clip_mode      <-- clipping mode
 *   epsilon        <-- precision for iterative gradient calculation
 *   extrap         <-- boundary gradient extrapolation coefficient
 *   clip_coeff     <-- clipping coefficient
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_field_gradient_scalar(const cs_field_t          *f,
                              cs_gradient_type_t         gradient_type,
                              cs_halo_type_t             halo_type,
                              int                        inc,
                              bool                       use_previous_t,
                              bool                       recompute_cocg,
                              int                        n_r_sweeps,
                              int                        verbosity,
                              int                        clip_mode,
                              double                     epsilon,
                              double                     extrap,
                              double                     clip_coeff,
                              cs_real_3_t      *restrict grad);

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   gradient_type  <-- gradient type
 *   halo_type      <-- halo type
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   use_previous_t <-- should we use values from the previous time step ?
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

void cs_field_gradient_potential(const cs_field_t          *f,
                                 cs_gradient_type_t         gradient_type,
                                 cs_halo_type_t             halo_type,
                                 int                        inc,
                                 bool                       use_previous_t,
                                 bool                       recompute_cocg,
                                 int                        n_r_sweeps,
                                 int                        hyd_p_flag,
                                 int                        verbosity,
                                 int                        clip_mode,
                                 double                     epsilon,
                                 double                     extrap,
                                 double                     clip_coeff,
                                 cs_real_3_t                f_ext[],
                                 cs_real_3_t      *restrict grad);

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f              <-- pointer to field
 *   gradient_type  <-- gradient type
 *   halo_type      <-- halo type
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   use_previous_t <-- should we use values from the previous time step ?
 *   recompute_cocg <-- should COCG FV quantities be recomputed ?
 *   n_r_sweeps     <-- if > 1, number of reconstruction sweeps
 *   verbosity      <-- verbosity level
 *   clip_mode      <-- clipping mode
 *   epsilon        <-- precision for iterative gradient calculation
 *   extrap         <-- boundary gradient extrapolation coefficient
 *   clip_coeff     <-- clipping coefficient
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void cs_field_gradient_vector(const cs_field_t          *f,
                              cs_gradient_type_t         gradient_type,
                              cs_halo_type_t             halo_type,
                              int                        inc,
                              bool                       use_previous_t,
                              int                        n_r_sweeps,
                              int                        verbosity,
                              int                        clip_mode,
                              double                     epsilon,
                              double                     clip_coeff,
                              cs_real_33_t     *restrict grad);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_OPERATOR_H__ */
