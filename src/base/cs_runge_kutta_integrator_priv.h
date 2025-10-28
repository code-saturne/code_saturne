#ifndef __RK_INTEGRATOR_PRIV_H__
#define __RK_INTEGRATOR_PRIV_H__

/*============================================================================
 * Explicit Runge-Kutta integrator utilities.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "alge/cs_blas.h"

#include "base/cs_base_accel.h"
#include "base/cs_dispatch.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define RK_HIGHEST_ORDER 4

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Available Runge-Kutta schemes */
/* ----------------------------- */

typedef enum {
  CS_RK_NONE,
  CS_RK1,
  CS_RK2,
  CS_RK3,
  CS_RK4
} cs_runge_kutta_scheme_t;

/* Coeffients adapted from the Butcher tableau */
/* ------------------------------------------- */

typedef struct {

  double *a;
  double *c;

} cs_runge_kutta_coeff_t;

/*  Descriptor of a Runge-Kutta integrator */
/* --------------------------------------- */

typedef struct {

  int rk_id;
  cs_runge_kutta_scheme_t scheme;

} cs_runge_kutta_def_t;

/* Generic Runge-Kutta integrator */
/* -------------------------------*/

typedef struct {

  cs_runge_kutta_scheme_t scheme;     // Selected RK scheme
  char                   *name;       // associated equation's or field's name
  int                     n_stages;   // number of stages
  int                     i_stage;    // Current stage index

  const cs_real_t        *dt;         // time step array
  cs_real_t              *scaled_dt;  // scaled time step used in projection
                                      // step per stage
  cs_lnum_t               n_elts;     // number of computational elements

  /* variable storage */
  cs_real_t *u_old;       // variable at beginning of the time stepping
  cs_real_t *u_new;       // updated variable

  cs_real_t *rhs_stages;  // RHS temporary storage


  cs_real_t *mass;         // mass array for the variables, equivalent to
                           // diag elements of mass matrix in an explicit
                           // time scheme

  cs_runge_kutta_coeff_t rk_coeff;

} cs_runge_kutta_integrator_t;

#endif /* __RK_INTEGRATOR_PRIV_H__ */
