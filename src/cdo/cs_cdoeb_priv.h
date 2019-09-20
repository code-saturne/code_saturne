#ifndef __CS_CDOEB_PRIV_H__
#define __CS_CDOEB_PRIV_H__

/*============================================================================
 * Definition of cs_cdovb_scaleq_t and cs_cdovb_vecteq structures
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_hodge.h"
#include "cs_equation_assemble.h"
#include "cs_equation_bc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdoeb_priv.h

  \brief Structures for building an algebraic CDO edge-based system for
         unsteady diffusion-reaction equations with source terms
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

/* Algebraic system for CDO edge-based discretization */

struct _cs_cdoeb_t {

  /* Ids related to the variable field and to the boundary flux field */
  int          var_field_id;
  int          bflux_field_id;

  /* System size */
  cs_lnum_t    n_dofs;

  /* Array storing the computed values */
  cs_real_t   *edge_values;
  cs_real_t   *edge_values_pre;

  /* Array storing the value arising from the contribution of all source
     terms */
  cs_real_t   *source_terms;

  /* Assembly process */
  cs_equation_assembly_t   *assemble;

  /* Boundary conditions */
  cs_flag_t                *edge_bc_flag;
  cs_cdo_enforce_bc_t      *enforce_essential_bc;

  /* Pointer of function to build the diffusion term */
  cs_hodge_t               *get_curlcurl_hodge;

  /* Mass matrix settings (useful for the unsteady and reaction terms) */
  cs_param_hodge_t          hdg_mass;
  cs_hodge_t               *get_mass_matrix;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOEB_PRIV_H__ */
