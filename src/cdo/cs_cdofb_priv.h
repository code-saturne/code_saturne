#ifndef __CS_CDOFB_PRIV_H__
#define __CS_CDOFB_PRIV_H__

/*============================================================================
 * Definition of cs_cdofb_scaleq_t and cs_cdofb_vecteq structures
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

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO face-based discretization */
struct  _cs_cdofb_t {

  /* System size (n_faces + n_cells) */
  cs_lnum_t                        n_dofs;

  /* Solution of the algebraic system at the last iteration
     DoF unknowns (x) + BCs */
  cs_real_t                       *face_values;

  /* Right-hand side related to cell dofs */
  cs_real_t                       *cell_rhs;

  /* Inverse of a diagonal matrix (block cell-cell) */
  cs_real_t                       *acc_inv;

  /* Lower-Left block of the full matrix (block cell-vertices).
     Access to the values thanks to the c2f connectivity */
  cs_real_t                       *acf;

  /* Array storing the value arising from the contribution of all source
     terms */
  cs_real_t                       *source_terms;

  /* Pointer of function to build the diffusion term */
  cs_hodge_t                      *get_stiffness_matrix;
  cs_hodge_t                      *get_diffusion_hodge;
  cs_cdo_diffusion_enforce_dir_t  *enforce_dirichlet;
  cs_cdo_diffusion_flux_trace_t   *boundary_flux_op;

  /* Pointer of function to build the advection term */
  cs_cdo_advection_t              *get_advection_matrix;
  cs_cdo_advection_bc_t           *add_advection_bc;

  /* Pointer of function to apply the time scheme */
  cs_cdo_time_scheme_t            *apply_time_scheme;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_PRIV_H__ */

