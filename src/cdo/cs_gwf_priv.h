#ifndef __CS_GWF_PRIV_H__
#define __CS_GWF_PRIV_H__

/*============================================================================
 * Structures/types related to the groundwater flow module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#include "cs_gwf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_gwf_single_phase_t
 *
 * \brief Structure to handle the modelling of a single-phase flows in a porous
 *        media considered as saturated or not. Several simplifications can be
 *        be operated in this context. Only the liquid phase is taken into
 *        account.
 */

typedef struct {

    /*!
   * @name Equation
   * @{
   *
   * \var richards

   * The Richards equation is the governing equation which corresponds to the
   * mass conservation of water. "hydraulic_head" is the associated variable
   * "permeability" is the diffusion property related to this equation.
   */

  cs_equation_t                *richards;

  /*!
   * @}
   * @name Darcy advection field
   * @{
   *
   * \var adv_field
   * Pointer to a \ref cs_adv_field_t structure. Darcy advective flux in the
   * liquid phase. This structure is used to define the advective term in
   * tracer equations.
   *
   * \var flux location
   * Indicate where the arrays defining the Darcy fluxes are located
   *
   * \var darcian_flux
   * Array storing the liquid Darcian flux in each location (for instance the
   * dual faces associated to each cell)
   *
   * \var darcian_boundary_flux
   * Array storing the normal Darcian flux across the boundary of the
   * computational domain for the liquid phase. This is an optional array.
   */

  cs_adv_field_t               *adv_field;
  cs_flag_t                     flux_location;
  cs_real_t                    *darcian_flux;
  cs_real_t                    *darcian_boundary_flux;

  /*!
   * @}
   * @name Properties related to the model
   * @{
   *
   * \var moisture_content
   * This quantity describes the level of saturation in a soil. This is a
   * constant value in case of a satured soil and variable one in case of a
   * unsaturated soil.
   *
   * \var soil_capacity
   * property attached to the unsteady term in the Richards equation
   */

  cs_property_t                *moisture_content;
  cs_property_t                *soil_capacity;

  /*!
   * @}
   * @name Additional fields/arrays
   * @{
   *
   * \var moisture_field
   * Pointer to a \ref cs_field_t structure. Structure storing the value of the
   * moisture content in each cell. This is an optional structure i.e. it may
   * be set to NULL (for instance in case of a satured soil on the full
   * domain).
   *
   * \var capacity_field
   * Pointer to a \ref cs_field_t structure. Structure storing the value of the
   * soil capacity in each cell. This is an optional structure i.e. it is
   * set to NULL in case of a satured soil on the full domain.
   *
   * \var pressure_head
   * Pointer to a \ref cs_field_t structure. Allocated only if the gravitation
   * effect is active. Location of this field depends on the discretization
   * scheme used to solve the Richards equation.
   * The pressure head is denoted by h, hydraulic head (unknowns solved in the
   * Richards eq.) is denoted by H and there are linked by:
   * h = H - gravity_potential
   *
   * \var head_in_law
   * Array of values located at the same location as the hydraulic head solved
   * in the Richards equation. The values stored in this array are used to
   * update related quantities with law such as Van Genuchten/Mualen.
   */

  cs_field_t                   *moisture_field;
  cs_field_t                   *capacity_field;
  cs_field_t                   *pressure_head;
  cs_real_t                    *head_in_law;

  /*!
   * @}
   */

} cs_gwf_single_phase_t;

/* --------------------------------------------------------------------------
 * Modelling context for two-phase flows in porous with two components (water
 * and another one like the hydrogen)
 * -------------------------------------------------------------------------- */

typedef struct {

  /* Set of equations associated to this module */
  cs_equation_t     *water_eq;  /* equation of conservation for the water
                                   specie */

  cs_equation_t     *gcomp_eq;   /* equation of conservation for the component
                                    denoted by "comp". For instance hydrogen */

  /* Additional fields/array related to heads */

  /* TODO */

  /* Properties */
  /* ---------- */

  /* TODO */

  /* Advection field = Darcy flux/velocity */
  /* ------------------------------------- */

  cs_adv_field_t  *darcy_l_field; /* Darcy advective flux in the liquid phase */
  cs_adv_field_t  *darcy_g_field; /* Darcy advective flux in the gaz phase */

  /* Additional parameters related to the Darcy advection fields */
  cs_flag_t        flux_location;   /* Indicate where the arrays are defined */

  /* Array defining the advection field (optional) */
  cs_real_t       *darcian_l_flux;
  cs_real_t       *darcian_g_flux;

  /* Array defining the normal flux of the advection field across the domain
     boundary (optional) */
  cs_real_t       *darcian_l_boundary_flux;
  cs_real_t       *darcian_g_boundary_flux;

} cs_gwf_two_phase_t;

/* --------------------------------------------------------------------------
 * Main set of parameters/structures related to the groundwater flow (GWF)
 * module
 * -------------------------------------------------------------------------- */

struct _gwf_t {

  cs_gwf_model_type_t   model;     /* Model (system of equations related to the
                                      chosen physical modelling) */
  cs_gwf_option_flag_t  flag;      /* Flag dedicated to general options to
                                    * handle the GWF module*/
  cs_flag_t             post_flag; /* Flag dedicated to the post-processing of
                                    * the GWF module */

  /* Properties */
  /* ---------- */

  /* Permeability characterizes the behavior of a soil w.r.t. to a component */
  cs_property_t   *permeability;
  cs_field_t      *permea_field;     /* Related cs_field_t structure at cells */

  /* Pointer to a structure cast on-the-fly which depends on the choice of
     model (for instance single-phase or two-phase flows in porous media) */
  void            *model_context;

  /* Associated tracers */
  /* ------------------ */

  /* Members related to the associated tracer equations */
  int                          n_tracers;
  cs_gwf_tracer_t            **tracers;
  cs_gwf_tracer_setup_t      **finalize_tracer_setup;  /* Function pointers */
  cs_gwf_tracer_add_terms_t  **add_tracer_terms;       /* Function pointers */

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PRIV_H__ */
