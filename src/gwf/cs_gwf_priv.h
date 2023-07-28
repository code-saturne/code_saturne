#ifndef __CS_GWF_PRIV_H__
#define __CS_GWF_PRIV_H__

/*============================================================================
 * Structures/types related to the groundwater flow module
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

#include "cs_advection_field.h"
#include "cs_equation.h"
#include "cs_gwf_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _gwf_darcy_flux_t  cs_gwf_darcy_flux_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux.
 *        cell_values could be set to NULL when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or NULL
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      cur2prev      true or false
 * \param[in, out] darcy         pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_darcy_update_t)(const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *cdoq,
                        const cs_real_t             *dof_values,
                        const cs_real_t             *cell_values,
                        cs_real_t                    t_eval,
                        bool                         cur2prev,
                        cs_gwf_darcy_flux_t         *darcy);

/*! \struct cs_gwf_darcy_flux_t
 *
 * \brief Structure to handle the Darcy flux
 */

struct _gwf_darcy_flux_t {

  /*!
   * \var adv_field
   * Pointer to a \ref cs_adv_field_t structure. Darcy advective flux in the
   * liquid phase. This structure is used to define the advective term in
   * tracer equations for instance.
   *
   * \var flux_location
   * Indicate where the arrays defining the Darcy fluxes are located
   *
   * \var flux_val
   * Array storing the liquid Darcian flux in each location (for instance the
   * dual faces associated to each cell)
   *
   * \var boundary_flux_val
   * Array storing the normal Darcian flux across the boundary of the
   * computational domain for the liquid phase. This is an optional array.
   *
   * \var simplified_boundary_update
   * If true, an constant approximation of the gradient is computed to define
   * the boundary flux. Otherwise, one tries to reproduce exatcly the
   * divergence of the advection field, the difference is distributed according
   * the area of each boundary face attached to a boundary vertex.
   * This option is only useful for vertex-based scheme.
   *
   * \var update_input
   * Context pointer for the update function or NULL if useless
   *
   * \var update_func
   * Pointer to the function which performs the update of the advection field
   */

  cs_adv_field_t               *adv_field;

  cs_flag_t                     flux_location;
  cs_real_t                    *flux_val;
  cs_real_t                    *boundary_flux_val;
  bool                          simplified_boundary_update;

  void                         *update_input;
  cs_gwf_darcy_update_t        *update_func;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the values of (potential) fields needed for the update of
 *        the Darcy velocity/fluxes.
 *
 * \param[in]  eq         pointer to an equation structure
 * \param[out] dof_vals   double pointer to the values (degrees of freedom)
 * \param[out] cell_vals  double pointer to the values (cell values)
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_get_value_pointers(const cs_equation_t       *eq,
                          cs_real_t                **p_dof_vals,
                          cs_real_t                **p_cell_vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in]      loc_flag   flag to define where the flux is defined
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_darcy_flux_t *
cs_gwf_darcy_flux_create(cs_flag_t         loc_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] p_darcy   pointer of pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_free(cs_gwf_darcy_flux_t  **p_darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] darcy   pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_log(cs_gwf_darcy_flux_t    *darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the definition of the advection field attached to a
 *        \ref cs_gwf_darcy_flux_t structure
 *        If the function pointer is set to NULL, then an automatic settings
 *        is done.
 *
 * \param[in]      connect         pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq            pointer to a cs_cdo_quantities_t structure
 * \param[in]      space_scheme    space discretization using this structure
 * \param[in]      update_context  pointer to the context for the update step
 * \param[in]      update_func     pointer to an update function or NULL
 * \param[in, out] darcy           pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_define(const cs_cdo_connect_t       *connect,
                         const cs_cdo_quantities_t    *quant,
                         cs_param_space_scheme_t       space_scheme,
                         void                         *update_context,
                         cs_gwf_darcy_update_t        *update_func,
                         cs_gwf_darcy_flux_t          *darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate the balance by zone (relying on the splitting arising from
 *         the boundary settings) for the advection field attached to a \ref
 *         cs_gwf_darcy_flux_t structure
 *
 * \param[in]       connect       pointer to a cs_cdo_connect_t structure
 * \param[in]       quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]       eqp           pointer to the set of equation parameters
 * \param[in, out]  darcy         pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_balance(const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant,
                          const cs_equation_param_t    *eqp,
                          cs_gwf_darcy_flux_t          *darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the associated Darcy flux over the boundary of the domain for
 *        each vertex of a boundary face.  Case of a vertex-based
 *        discretization and single-phase flows in porous media (saturated or
 *        not).
 *
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      eq       pointer to the equation related to this Darcy flux
 * \param[in, out] adv      pointer to the Darcy advection field
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_update_on_boundary(cs_real_t                t_eval,
                                     const cs_equation_t     *eq,
                                     cs_adv_field_t          *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the associated Darcy flux over the boundary of the domain for
 *        each vertex of a boundary face without using an equation (i.e. there
 *        is no associated boundary condition).
 *        Case of a vertex-based discretization.
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      cell_vel   Darcy velocity in each cell
 * \param[in, out] adv        pointer to the Darcy advection field to update
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_update_on_boundary_wo_eq(const cs_cdo_connect_t     *connect,
                                           const cs_cdo_quantities_t  *cdoq,
                                           cs_real_t                  *cell_vel,
                                           cs_adv_field_t             *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update head values (pressure head or head values for laws)
 *        Up to now, this is only used for single-phase flows in porous media
 *        (saturated or not case).
 *
 * \param[in]      connect         pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq            pointer to a cs_cdo_quantities_t structure
 * \param[in]      richards        pointer to the Richards equation
 * \param[in]      option_flag     calculation option related to the GWF module
 * \param[in, out] pressure_head   pressure head field
 * \param[in, out] head_in_law     values of the head used in law
 * \param[in]      cur2prev        true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_update_head(const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *cdoq,
                   const cs_equation_t         *richards,
                   cs_flag_t                    option_flag,
                   cs_field_t                  *pressure_head,
                   cs_real_t                    head_in_law[],
                   bool                         cur2prev);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PRIV_H__ */
