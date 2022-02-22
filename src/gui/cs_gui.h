#ifndef __CS_GUI_H__
#define __CS_GUI_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
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
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Specific heat variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCPVA
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscpva, CSCPVA) (void);

/*----------------------------------------------------------------------------
 * Volumic viscosity variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCVVVA (ICP)
 * *****************
 *
 * INTEGER          IVISCV     -->   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvvva, CSVVVA) (int *iviscv);

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * subroutine csivis
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (void);

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIDTV ()
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(csidtv, CSIDTV) (void);

/*----------------------------------------------------------------------------
 * Hydrostatic pressure parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIPHY ()
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (csiphy, CSIPHY) (void);

/*----------------------------------------------------------------------------
 * Hydrostatic equilibrium parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCFGP (ICFGRP)
 * *****************
 *
 * INTEGER          ICFGRP  -->   hydrostatic equilibrium
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscfgp, CSCFGP) (int *icfgrp);

/*----------------------------------------------------------------------------
 * Restart parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISUI
 * *****************
 *
 * INTEGER          NTSUIT  -->   checkpoint frequency
 * INTEGER          ICCFVG  -->   restart with frozen field
 *----------------------------------------------------------------------------*/


void CS_PROCF (csisui, CSISUI) (int *ntsuit,
                                int *iccvfg);

/*----------------------------------------------------------------------------
 * Time passing parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTIME
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstime, CSTIME) (void);

/*----------------------------------------------------------------------------
 *
 * Fortran Interface:
 *
 * SUBROUTINE UINUM1
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (double *cdtvar);

/*----------------------------------------------------------------------------
 * Global numerical parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNUM2
 * *****************
 *
 * INTEGER          RELAXP  -->   pressure relaxation
 * INTEGER          IMRGRA  -->   gradient reconstruction
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2) (double  *relaxp,
                                int     *imrgra);

/*----------------------------------------------------------------------------
 * Turbulence initialization parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTINI
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstini, CSTINI) (void);

/*----------------------------------------------------------------------------
 * Define porosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPORO
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiporo, UIPORO) (void);

/*----------------------------------------------------------------------------
 * User momentum source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsnv (ncelet, vel, tsexp, tsimp)
 * *****************
 *
 * integer          ncelet   <--  number of cells with halo
 * double precision vel      <--  fluid velocity
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsnv, UITSNV)(const cs_real_3_t *restrict vel,
                              cs_real_3_t       *restrict tsexp,
                              cs_real_33_t      *restrict tsimp);

/*----------------------------------------------------------------------------
 * User scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitssc (f_id, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          idarcy   <--  groundwater module activation
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitssc, UITSSC)(const int                  *idarcy,
                              const int                  *f_id,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp);

/*----------------------------------------------------------------------------
 * Thermal scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsth (f_id, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsth, UITSTH)(const int                  *f_id,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp);

/*----------------------------------------------------------------------------
 * User law for material Properties
 *
 * Fortran Interface:
 *
 * subroutine uiphyv
 * *****************
 *
 * integer          iviscv   <--  pointer for volumic viscosity viscv
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const int       *iviscv);

/*----------------------------------------------------------------------------
 * extra operations
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIEXOP
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiexop, UIEXOP)(void);

/*----------------------------------------------------------------------------
 * groundwater model : read laws for capacity, saturation and permeability
 *
 * Fortran Interface:
 *
 * subroutine uidapp
 * *****************
 * integer         permeability    <--  permeability type
 * integer         diffusion       <--  diffusion type
 * integer         gravity         <--  check if gravity is taken into account
 * double          gravity_x       <--  gravity direction
 * double          gravity_y       <--  gravity direction
 * double          gravity_z       <--  gravity direction
 * integer         unsaturated     <--  unsaturated zone taken into account
 *----------------------------------------------------------------------------*/

void CS_PROCF (uidapp, UIDAPP) (const int       *permeability,
                                const int       *diffusion,
                                const int       *gravity,
                                const cs_real_t *gravity_x,
                                const cs_real_t *gravity_y,
                                const cs_real_t *gravity_z,
                                const int       *unsaturated);

/*----------------------------------------------------------------------------
 * Define error estimators
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIERES
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uieres, UIERES) (int *iescal,
                                int *iespre,
                                int *iesder,
                                int *iescor,
                                int *iestot);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *----------------------------------------------------------------------------*/

void
cs_gui_finalize(void);

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
 * Treatment of gravity and fluid physical properties
 * Initialize reference pressure and temperature if present
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_properties(void);

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
cs_gui_error_estimator(int *iescal,
                       int *iespre,
                       int *iesder,
                       int *iescor,
                       int *iestot);

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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_H__ */
