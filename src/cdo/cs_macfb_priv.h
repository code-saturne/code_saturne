#ifndef __CS_MACFB_PRIV_H__
#define __CS_MACFB_PRIV_H__

/*============================================================================
 * Definition of cs_macfb_scaleq_t and cs_macfb_vecteq structures
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_hodge.h"
#include "cs_macfb_advection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for MAC face-based discretization */

struct _cs_macfb_t {

  /* Ids related to the variable field and to the boundary flux field */

  int var_field_id;
  int bflux_field_id;

  /* System size (n_faces) */

  cs_lnum_t n_faces;
  cs_lnum_t n_dofs;

  /* Solution of the algebraic system DoF unknowns (x) + BCs */

  cs_real_t *face_values;     /* At the last iteration */
  cs_real_t *face_values_pre; /* At the previous iteration */

  /* Array storing the value arising from the contribution of all source
     terms (only allocated to n_cells) */

  cs_real_t *source_terms;

  /* Pointer of function to build the diffusion term */

  cs_hodge_t         **diffusion_hodge;
  cs_cdo_enforce_bc_t *enforce_dirichlet;
  cs_cdo_enforce_bc_t *enforce_robin_bc;
  cs_cdo_enforce_bc_t *enforce_sliding;

  /* Pointer of functions to define the advection term:
   * advection_open is called first, then advection_main which calls
   * advection_scheme and after the build step, advection_close is called last
   */

  cs_macfb_adv_open_hook_t  *advection_open;
  cs_macfb_adv_build_t      *advection_main;
  cs_macfb_adv_close_hook_t *advection_close;
  cs_macfb_adv_scheme_t     *advection_scheme;

  void *advection_input;

  /* If one needs to build a local hodge op. for time and reaction */

  cs_hodge_param_t    mass_hodgep;
  cs_hodge_t        **mass_hodge;
  cs_hodge_compute_t *get_mass_matrix;
};

typedef struct _cs_macfb_t cs_macfb_priv_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the advection-related parameters in the context structure of
 *         MAC face-based schemes
 *
 * \param[in]      eqp    pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb    pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] eqc    pointer to a cs_macfb_priv_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_set_advection_function(const cs_equation_param_t *eqp,
                                     cs_equation_builder_t     *eqb,
                                     cs_macfb_priv_t           *eqc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MACFB_PRIV_H__ */
