#ifndef __CS_CDOCB_PRIV_H__
#define __CS_CDOCB_PRIV_H__

/*============================================================================
 * Definition of cs_cdocb_scaleq_t and cs_cdocb_vecteq structures
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Context related to the resolution of a saddle point problem */

typedef struct {

  cs_real_t     *div_op;    /* Block related to the -divergence (block
                               A_{10}) */

  /* Arrays split according to the block shape. U is interlaced or not
   * according to the SLES strategy */

  cs_lnum_t      n_faces;       /* local number of DoFs for each component
                                 * of the velocity */
  cs_lnum_t      n_cells;       /* local number of DoFs for the pressure */

  cs_real_t     *flux;          /* flux values at faces */
  cs_real_t     *potential;     /* pressure values at cells */

  cs_sles_t     *sles;          /* main SLES structure */
  cs_sles_t     *schur_sles;    /* auxiliary SLES for the Schur complement
                                 * May be NULL */

  cs_real_t      graddiv_coef;  /* value of the grad-div coefficient in case
                                 * of augmented system */

} cs_cdocb_monolithic_sles_t;

/* Algebraic system for CDO cell-based discretization */

struct _cs_cdocb_t {

  /* Ids related to the variable field and to the boundary flux field */

  int          var_field_id;
  int          bflux_field_id;

  /* System size (n_faces + n_cells) */

  cs_lnum_t    n_faces;
  cs_lnum_t    n_cells;
  cs_lnum_t    n_dofs;

  /* Solution of the algebraic system DoF unknowns (x) + BCs */

  cs_real_t   *flux;       /* At the last iteration */
  cs_real_t   *flux_pre;   /* At the previous iteration */

  /* Array storing the value arising from the contribution of all source
     terms (only allocated to n_cells) */

  cs_real_t                   *source_terms;

  /* Members used for the cell-wise building of the linear system */
  /* ------------------------------------------------------------ */

  /* Array storing the local cell-wise divergence operator (one by thread) */

  cs_real_t                  **div_op_cw;

  /* Pointer of function to build the diffusion term */

  cs_hodge_t                 **diff_hodge;
  cs_hodge_compute_t          *compute_diff_hodge;

  /* Boundary conditions */
  /* ------------------- */

  cs_cdo_enforce_bc_t         *enforce_dirichlet;
  cs_cdo_enforce_bc_t         *enforce_neumann;
  cs_cdo_enforce_bc_t         *enforce_robin_bc;

  /* Linear system */
  /* ------------- */

  cs_cdocb_monolithic_sles_t  *msles;
  cs_sles_t                   *schur_sles;

  cs_real_t                    graddiv_coef;

};

typedef struct _cs_cdocb_t  cs_cdocb_priv_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOCB_PRIV_H__ */
