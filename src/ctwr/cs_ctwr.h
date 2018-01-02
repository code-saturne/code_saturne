#ifndef __CS_CTWR_H__
#define __CS_CTWR_H__

/*============================================================================
 * Main for cooling towers related functions
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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Type of cooling tower model */

typedef enum {

  CS_CTWR_NONE,              /*!< no cooling tower model */
  CS_CTWR_POPPE,             /*!< Poppe's model */
  CS_CTWR_MERKEL             /*!< Merkel's model */

} cs_ctwr_model_t;

/*! Type of cooling tower exchange zone */

typedef enum {

  CS_CTWR_COUNTER_CURRENT = 1,   /*!< counter-current zone */
  CS_CTWR_CROSS_CURRENT,         /*!< cross-current zone */
  CS_CTWR_RAIN                   /*!< rain zone */

} cs_ctwr_zone_type_t;

typedef struct _cs_ctwr_zone_t cs_ctwr_zone_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cooling tower exchange zone
 *
 * \param[in]   zone_criteria   Zone selection criteria
 * \param[in]   model           model type
 * \param[in]   zone_type       exchange zone type
 * \param[in]   delta_t         Imposed delta temperature delta between inlet
 *                              and oulet of the zone
 * \param[in]   relax           Relaxation of the imposed delta temperature
 * \param[in]   t_l_bc          Liquid water temperature at the inlet
 * \param[in]   q_l_bc          Mass flow rate at the inlet
 * \param[in]   xap             Beta_x_0 of the exchange law
 * \param[in]   xnp             Exponent n of the exchange law
 * \param[in]   surface         Total Surface of ingoing water
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define(const char           zone_criteria[],
               cs_ctwr_model_t      model,
               cs_ctwr_zone_type_t  zone_type,
               cs_real_t            delta_t,
               cs_real_t            relax,
               cs_real_t            t_l_bc,
               cs_real_t            q_l_bc,
               cs_real_t            xap,
               cs_real_t            xnp,
               cs_real_t            surface);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map fields used by the cooling tower module to pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_field_pointer_map(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change mass source term from the evaporating liquid to the
 *        bulk, humid air.
 *
 * Careful, this is different from an injection source term, which would
 * normally be handled with 'cs_user_mass_source_term'
 *
 * \param[in]   p0              Reference pressure
 * \param[in]   molmassrat      Dry air to water vapor molecular mass ratio
 * \param[in]   mass_source     Mass source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_bulk_mass_source_term(const cs_real_t   p0,
                              const cs_real_t   molmassrat,
                              cs_real_t         mass_source[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define zones.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different packing zones.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy cs_ctwr_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_all_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log Packing zone definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform balances in packing zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_balance(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the field variables
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_init_field_vars(cs_real_t  rho0,
                        cs_real_t  t0,
                        cs_real_t  p0,
                        cs_real_t  molmassrat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the flow variables relevant to the cooling tower scalars
 * inside the packing zones
 *
 * \param[in,out] liq_mass_flow Liquid mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_init_flow_vars(cs_real_t  liq_mass_flow[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the thermo physical properties fields for the humid air and
 *        the liquid
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_phyvar_update(cs_real_t  rho0,
                      cs_real_t  t0,
                      cs_real_t  p0,
                      cs_real_t  molmassrat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in]     p0            Reference pressure
 * \param[in]     molmassrat    dry air to water vapor molecular mass ratio
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_source_term(int              f_id,
                    const cs_real_t  p0,
                    const cs_real_t  molmassrat,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Convert injected liquid scalars from and to their transported form.
 *
 * \param[in]   iflag     1: Convert transported variables to physical variables
 *                        2: Convert physical variables to
 *                           transported variables
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_transport_vars(int  iflag);

/*----------------------------------------------------------------------------
 * Get pointer to exchange area.
 *
 * parameters:
 *   ct_id  <--  Id (0 to n-1) of exchange area
 *
 * returns:
 *   pointer to exchange area structure
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t *
cs_ctwr_by_id(int ct_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_H__ */
