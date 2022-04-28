#ifndef __CS_CDOVB_PRIV_H__
#define __CS_CDOVB_PRIV_H__

/*============================================================================
 * Definition of cs_cdovb_scaleq_t and cs_cdovb_vecteq structures
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_cdo_advection.h"
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the reaction term for a CDO vertex-based scheme
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      mass_hodge  pointer to a Hodge structure or NULL if useless
 * \param[in, out] eqb         pointer to the equation builder structure
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_reaction_t)(const cs_equation_param_t    *eqp,
                      const cs_cell_mesh_t         *cm,
                      const cs_hodge_t             *mass_hodge,
                      const cs_equation_builder_t  *eqb,
                      cs_cell_builder_t            *cb,
                      cs_cell_sys_t                *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the time scheme for a CDO vertex-based scheme
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      mass_hodge  pointer to a Hodge structure or NULL if useless
 * \param[in]      inv_dtcur   value of 1./dt for the current time step
 * \param[in, out] eqb         pointer to the equation builder structure
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_time_t)(const cs_equation_param_t   *eqp,
                  const cs_cell_mesh_t        *cm,
                  const cs_hodge_t            *mass_hodge,
                  const double                 inv_dtcur,
                  cs_equation_builder_t       *eqb,
                  cs_cell_builder_t           *cb,
                  cs_cell_sys_t               *csys);

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

  /* Boundary conditions */

  cs_flag_t                *vtx_bc_flag;
  cs_cdo_enforce_bc_t      *enforce_dirichlet;
  cs_cdo_enforce_bc_t      *enforce_robin_bc;

  /* Only for vector-valued variables */

  cs_cdo_enforce_bc_t      *enforce_sliding;

  /* Pointer of function to build the diffusion term */

  cs_hodge_t              **diffusion_hodge;
  cs_hodge_compute_t       *get_stiffness_matrix;

  /* Pointer of function to build the advection term */

  cs_cdovb_advection_t     *get_advection_matrix;
  cs_cdovb_advection_bc_t  *add_advection_bc;

  /* Pointer of function to build the unsteady term */

  cs_cdovb_reaction_t      *add_reaction_term;

  /* Pointer of function to build the unsteady term */

  cs_cdovb_time_t          *add_unsteady_term;

  /* If one needs to build a local Hodge operator for the unsteady and/or the
     reaction term(s) */

  cs_hodge_param_t          mass_hodgep;
  cs_hodge_t              **mass_hodge;
  cs_hodge_compute_t       *get_mass_matrix;

};

/* Scalar-valued and vector-valued equations rely on the same pattern */

typedef struct _cs_cdovb_t cs_cdovb_scaleq_t;
typedef struct _cs_cdovb_t cs_cdovb_vecteq_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_PRIV_H__ */
