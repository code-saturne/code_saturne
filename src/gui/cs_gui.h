#ifndef __CS_GUI_H__
#define __CS_GUI_H__

/*============================================================================
 * Management of the GUI parameters file: main parameters
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_volume_zone.h"

#include "mei_evaluate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Thermal model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTHER ()
 * *****************
 *
 *----------------------------------------------------------------------------*/


void CS_PROCF (csther, CSTHER) (void);

/*----------------------------------------------------------------------------
 * Turbulence model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTURB
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (csturb, CSTURB) (void);

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
 * User thermal scalar.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UITHSC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uithsc, UITHSC) (void);

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
 * INTEGER          ILEAUX  -->   restart with auxiliary
 * INTEGER          ICCFVG  -->   restart with frozen field
 *----------------------------------------------------------------------------*/


void CS_PROCF (csisui, CSISUI) (int *ntsuit,
                                int *ileaux,
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
 * INTEGER          EXTRAG  -->   wall pressure extrapolation
 * INTEGER          IMRGRA  -->   gradient reconstruction
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2) (double *relaxp,
                                double *extrag,
                                   int *imrgra);

void CS_PROCF (csphys, CSPHYS) (
                                double       *viscv0,
                                double       *visls0,
                                const    int *itempk);

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * subroutine cssca2
 * *****************
 *
 * integer          iturt    -->  turbulent flux model
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (int        *iturt);

void CS_PROCF (cssca3, CSSCA3) (double     *visls0);

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
 * Solver taking a scalar porosity into account
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIIPSU
 * *****************
 *
 * INTEGER          IPOROS     -->   porosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiipsu, UIIPSU) (int *iporos);

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
 * Variables and user scalars initialization.
 *
 * Fortran Interface:
 *
 * subroutine uiiniv
 * *****************
 *
 * integer          isuite   <--  restart indicator
 * integer          idarcy   <--  groundwater module activation
 * integer          iccfth   <--  type of initialization (compressible model)
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int          *isuite,
                              const int          *idarcy,
                                    int          *iccfth);

/*----------------------------------------------------------------------------
 * User law for material Properties
 *
 * Fortran Interface:
 *
 * subroutine uiphyv
 * *****************
 *
 * integer          iviscv   <--  pointer for volumic viscosity viscv
 * integer          itempk   <--  pointer for temperature (in K)
 * double precision visls0   <--  diffusion coefficient of the scalars
 * double precision viscv0   <--  volumic viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const cs_int_t  *iviscv,
                              const cs_int_t  *itempk,
                              const cs_real_t *visls0,
                              const cs_real_t *viscv0);

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
 * Define fans with GUI
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIFANS
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uifans, UIFANS) (void);

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

/*----------------------------------------------------------------------------
 * Initialize GUI reader structures.
 *----------------------------------------------------------------------------*/

void
cs_gui_init(void);

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *----------------------------------------------------------------------------*/

void
cs_gui_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add notebook variable to a formula.
 *
 * \param[in, out]  ev_law  pointer to MEI formula structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_add_notebook_variables(mei_tree_t  *ev_law);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute GUI-defined head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * cku11, cku22, cku33, cku12, cku13, cku23.
 *
 * \param[in]       zone  pointer to zone structure
 * \param[in, out]  cku   head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_head_losses(const cs_zone_t  *zone,
                   cs_real_t         cku[][6]);

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

/*----------------------------------------------------------------------------
 * 1D profile postprocessing
 *----------------------------------------------------------------------------*/

void
cs_gui_profile_output(void);

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
 * Initialization choice of the reference variables parameters.
 *
 * parameters:
 *   name            <--   parameter name
 *   value           -->   parameter value
 *----------------------------------------------------------------------------*/

void
cs_gui_reference_initialization(const char  *param,
                                double      *value);

/*----------------------------------------------------------------------------
 * Get thermal scalar model.
 *
 * return:
 *   value of itherm
 *----------------------------------------------------------------------------*/

int
cs_gui_thermal_model(void);

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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_H__ */
