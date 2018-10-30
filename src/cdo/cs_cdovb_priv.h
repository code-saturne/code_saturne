#ifndef __CS_CDOVB_PRIV_H__
#define __CS_CDOVB_PRIV_H__

/*============================================================================
 * Definition of cs_cdovb_scaleq_t and cs_cdovb_vecteq structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_defs.h"
#include "cs_hodge.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_advection.h"
#include "cs_cdo_time.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdovb_priv.h

  \brief Structures for building an algebraic CDO vertex-based system for
         unsteady convection-diffusion-reaction equations with source terms
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Structure definitions
 *============================================================================*/

/* Algebraic system for CDO vertex-based discretization */

struct _cs_cdovb_t {

  /* Ids related to the variable field and to the boundary flux field */
  int          var_field_id;
  int          bflux_field_id;

  /* System size */
  cs_lnum_t    n_dofs;

  /* Array storing the value arising from the contribution of all source
     terms */
  cs_real_t   *source_terms;

  /* Array for extra-operations */
  cs_real_t   *cell_values;     /* NULL if not requested */

  /* Pointer of function to build the diffusion term */
  cs_hodge_t                      *get_stiffness_matrix;
  cs_cdo_diffusion_enforce_bc_t   *enforce_dirichlet;
  cs_cdo_diffusion_enforce_bc_t   *apply_robin_bc;
  cs_flag_t                       *vtx_bc_flag;

  /* Pointer of function to build the advection term */
  cs_cdo_advection_t              *get_advection_matrix;
  cs_cdo_advection_bc_t           *add_advection_bc;

  /* Pointer of function to apply the time scheme */
  cs_cdo_time_scheme_t            *apply_time_scheme;

  /* If one needs to build a local hodge op. for time and reaction */
  cs_param_hodge_t                 hdg_mass;
  cs_hodge_t                      *get_mass_matrix;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_PRIV_H__ */
