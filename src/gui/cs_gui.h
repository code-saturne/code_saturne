#ifndef __CS_GUI_H__
#define __CS_GUI_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Arguments passed by context pointer to cs_meg_* functions */

typedef struct {

  const cs_zone_t    *zone;    /*<! Pointer to zone */
  const cs_field_t  **fields;  /*<! Array of field pointers (NULL-terminated) */

} cs_gui_volume_meg_context_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Specific heat variable or constant indicator.
 *----------------------------------------------------------------------------*/

void
cs_gui_cp_params (void);

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *----------------------------------------------------------------------------*/

void
cs_gui_laminar_viscosity(void);

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *----------------------------------------------------------------------------*/

void
cs_gui_dt(void);

/*----------------------------------------------------------------------------
 * Hydrostatic pressure parameter.
 *----------------------------------------------------------------------------*/

void
cs_gui_hydrostatic_pressure(void);

/*----------------------------------------------------------------------------
 * Hydrostatic equilibrium parameter.
 *----------------------------------------------------------------------------*/

void
cs_gui_hydrostatic_equ_param(void);

/*----------------------------------------------------------------------------
 * Time passing parameters.
 *----------------------------------------------------------------------------*/

void
cs_gui_dt_param(void);

/*----------------------------------------------------------------------------
 * Define porosity.
 *----------------------------------------------------------------------------*/

void
cs_gui_porosity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read GUi-defined Checkpoint parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_checkpoint_parameters(void);

/*----------------------------------------------------------------------------
 * Space scheme options, linear solver precision and time step factor
 *----------------------------------------------------------------------------*/

void
cs_gui_equation_parameters(void);

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *----------------------------------------------------------------------------*/

void
cs_gui_finalize(void);

/*----------------------------------------------------------------------------
 * Return a pointer to equation parameters based on a field or equation name.
 *
 * parameters:
 *   name <-- field or equation name
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_gui_get_equation_param(const char  *name);

/*-----------------------------------------------------------------------------
 * Get value of reference fluid properties parameter.
 *
 * parameters:
 *   name            <--   parameter name
 *   value           -->   parameter value
 *----------------------------------------------------------------------------*/

void
cs_gui_fluid_properties_value(const char  *param,
                              double      *value);

/*----------------------------------------------------------------------------
 * Groundwater model : read laws for capacity, saturation and permeability
 *
 * parameters:
 *   permeability  <--  permeability type
 *   diffusion     <--  diffusion type
 *   unsaturated   <--  unsaturated zone taken into account
 *----------------------------------------------------------------------------*/

void
cs_gui_groundwater_property_laws(int  permeability,
                                 int  diffusion,
                                 int  unsaturated);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute GUI-defined head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * cku11, cku22, cku33, cku12, cku13, cku23.
 *
 * \param[in]       zone       pointer to zone structure
 * \param[in]       cvara_vel  velocity values at the  previous time step
 * \param[in, out]  cku        head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_head_losses(const cs_zone_t   *zone,
                   const cs_real_3_t *cvara_vel,
                   cs_real_t          cku[][6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply initial conditions based on GUI-defined settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_initial_conditions(void);

/*-----------------------------------------------------------------------------
 * Selection of linear solvers.
 *----------------------------------------------------------------------------*/

void
cs_gui_linear_solvers(void);

/*----------------------------------------------------------------------------
 * User momentum source terms.
 *
 * parameters:
 *   vel      <--  fluid velocity
 *   tsexp    -->  explicit source terms
 *   tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void
cs_gui_momentum_source_terms(const cs_real_3_t  *vel,
                             cs_real_3_t        *tsexp,
                             cs_real_33_t       *tsimp);

/*-----------------------------------------------------------------------------
 * Define global numerical options.
 *----------------------------------------------------------------------------*/

void
cs_gui_numerical_options(void);

/*-----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_parallel_io(void);

/*-----------------------------------------------------------------------------
 * Set partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_gui_partition(void);

/*-----------------------------------------------------------------------------
 * Set MPI related algorithm options
 *----------------------------------------------------------------------------*/

void
cs_gui_mpi_algorithms(void);

/*----------------------------------------------------------------------------
 * Treatment of physical constants (gravity and Coriolis).
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_constants(void);

/*----------------------------------------------------------------------------
 * Treatment of gravity and fluid physical properties
 * Initialize reference pressure and temperature if present
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_properties(void);

/*----------------------------------------------------------------------------
 * User law for material properties
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_variable(void);

/*----------------------------------------------------------------------------
 * Determine porosity model type
 *----------------------------------------------------------------------------*/

void
cs_gui_porous_model(void);

/*-----------------------------------------------------------------------------
 * Get initial value from property markup.
 *
 * parameters:
 *   property_name      <--  name of the property
 *   value              -->  new initial value of the property
 *----------------------------------------------------------------------------*/

void
cs_gui_properties_value(const char  *property_name,
                        double      *value);

/*-----------------------------------------------------------------------------
 * Get value of property markup for fluid of given id
 *
 * parameters:
 *   fluid_id           <--   fluid index
 *   property_name      <--   name of the property
 *   value              -->   new initial value of the property
 *----------------------------------------------------------------------------*/

void
cs_gui_properties_value_by_fluid_id(const int    fluid_id,
                                    const char  *property_name,
                                    double      *value);

/*----------------------------------------------------------------------------
 * Read minimum / maximum values (used in clipping) and turbulent flux model
 * for additional user or model variables.
 *
 * Also read reference dynamic and user scalar viscosity
 *----------------------------------------------------------------------------*/

void
cs_gui_scalar_model_settings(void);

/*----------------------------------------------------------------------------
 * Thermal model.
 *----------------------------------------------------------------------------*/

void
cs_gui_thermal_model(void);

/*----------------------------------------------------------------------------
 * Get thermal scalar model.
 *
 * return:
 *   value of itherm*10 + (temperature variant flag), or -1 if not defined
 *----------------------------------------------------------------------------*/

int
cs_gui_thermal_model_code(void);

/*----------------------------------------------------------------------------
 * Time moments definition
 *----------------------------------------------------------------------------*/

void
cs_gui_time_moments(void);

/*-----------------------------------------------------------------------------
 * Set turbomachinery model
 *----------------------------------------------------------------------------*/

void
cs_gui_turbomachinery(void);

/*-----------------------------------------------------------------------------
 * Set turbomachinery options.
 *----------------------------------------------------------------------------*/

void
cs_gui_turbomachinery_rotor(void);

/*----------------------------------------------------------------------------
 * Turbulence model
 *----------------------------------------------------------------------------*/

void
cs_gui_turb_model(void);

/*----------------------------------------------------------------------------
 * Define reference length and reference velocity for the initialization of
 * the turbulence variables
 *----------------------------------------------------------------------------*/

void
cs_gui_turb_ref_values(void);

/*----------------------------------------------------------------------------
 * Logging output for MEI usage.
 *----------------------------------------------------------------------------*/

void
cs_gui_usage_log(void);

/*----------------------------------------------------------------------------
 * Define user variables through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_user_variables(void);

/*----------------------------------------------------------------------------
 * Define user arrays through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_user_arrays(void);

/*----------------------------------------------------------------------------
 * Define user calculator functions through the GUI
 *----------------------------------------------------------------------------*/

void
cs_gui_calculator_functions(void);

/*----------------------------------------------------------------------------
 * Define balance by zone through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_balance_by_zone(void);

/*----------------------------------------------------------------------------
 * Define pressure drop through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_pressure_drop_by_zone(void);

/*----------------------------------------------------------------------------
 * Define fans through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_define_fans(void);

/*----------------------------------------------------------------------------
 * Define error estimator through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_error_estimator(void);

/*----------------------------------------------------------------------------
 * Define volume and boundary zones through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_zones(void);

/*----------------------------------------------------------------------------
 * Define internal coupling through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_internal_coupling(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add new volume MEG function context info.
 *
 * \param[in]  zone      pointer to associated zone
 * \param[in]  fields    array of field pointers
 * \param[in]  n_fields  number gof field pointers
 *
 * \return: pointer to MEG context info
 */
/*----------------------------------------------------------------------------*/

const cs_gui_volume_meg_context_t *
cs_gui_add_volume_meg_context(const  cs_zone_t   *zone,
                              const  cs_field_t  *fields[],
                              const  int          n_fields);

/*----------------------------------------------------------------------------
 * Define user scalar source terms.
 *----------------------------------------------------------------------------*/

void
cs_gui_scalar_source_terms(cs_field_t        *f,
                           const cs_real_t   *pvar,
                           cs_real_t         *tsexp,
                           cs_real_t         *tsimp);

/*----------------------------------------------------------------------------
 * Define user thermal scalar source terms
 *----------------------------------------------------------------------------*/

void
cs_gui_thermal_source_terms(cs_field_t        *f,
                            const cs_real_t   *pvar,
                            cs_real_t         *tsexp,
                            cs_real_t         *tsimp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read GUI defined time tables
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_time_tables(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_H__ */
