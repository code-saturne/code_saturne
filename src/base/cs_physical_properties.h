#ifndef __CS_PHYSICAL_PROPERTIES_H__
#define __CS_PHYSICAL_PROPERTIES_H__

/*============================================================================
 * Physical properties computation and management.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

typedef enum {

  CS_PHYS_PROP_PLANE_PH,
  CS_PHYS_PROP_PLANE_PT,
  CS_PHYS_PROP_PLANE_PS,
  CS_PHYS_PROP_PLANE_PU,
  CS_PHYS_PROP_PLANE_PV,
  CS_PHYS_PROP_PLANE_TS,
  CS_PHYS_PROP_PLANE_TX,

} cs_phys_prop_thermo_plane_type_t;

typedef enum {

  CS_PHYS_PROP_PRESSURE,
  CS_PHYS_PROP_TEMPERATURE,
  CS_PHYS_PROP_ENTHALPY,
  CS_PHYS_PROP_ENTROPY,
  CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY,
  CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY,
  CS_PHYS_PROP_SPECIFIC_VOLUME,
  CS_PHYS_PROP_DENSITY,
  CS_PHYS_PROP_INTERNAL_ENERGY,
  CS_PHYS_PROP_QUALITY,
  CS_PHYS_PROP_THERMAL_CONDUCTIVITY,
  CS_PHYS_PROP_DYNAMIC_VISCOSITY,
  CS_PHYS_PROP_SPEED_OF_SOUND

} cs_phys_prop_type_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define thermal table.
 *----------------------------------------------------------------------------*/

void
cs_thermal_table_set(const char                        *material,
                     const char                        *method,
                     const char                        *reference,
                     cs_phys_prop_thermo_plane_type_t   thermo_plane,
                     cs_temperature_scale_t             temp_scale);

/*----------------------------------------------------------------------------
 * Finalize thermal table
 *----------------------------------------------------------------------------*/

void
cs_thermal_table_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get backend set for CoolProp.
 *
 * Returns NULL if CoolProp is not used or backend not set yet.
 *
 * \return  pointer to CoolProp backend.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_physical_properties_get_coolprop_backend(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set backend for CoolProp.
 *
 * Ignored if CoolProp is not used.
 *
 * When called from user-defined functions, this should be set from
 * cs_user_model rather than cs_user_parameters, as some reference property
 * values may be computed before calling cs_user_parameters.
 *
 * A few primary backends in CoolProp are:
 *
 * - "HEOS": The Helmholtz Equation of State backend for use with pure and
 *   pseudo-pure fluids, and mixtures, all of which are based on multi-parameter
 *   Helmholtz Energy equations of state.
 * - "REFPROP": only if REFPROP library is available
 *   (set ALTERNATIVE_REFPROP_PATH environment variable if needed)
 * - "INCOMP": Incompressible backend (for pure fluids)
 * - "TTSE&XXXX": TTSE backend, with tables generated using the XXXX backend
 *   where XXXX is one of the base backends("HEOS", "REFPROP", etc.)
 * - "BICUBIC&XXXX": Bicubic backend, with tables generated using the XXXX
 *   backend where XXXX is one of the base backends("HEOS", "REFPROP", etc.)
 *
 * \param[in]  backend  backend name
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_properties_set_coolprop_backend(const char  *backend);

/*----------------------------------------------------------------------------
 * Compute a physical property.
 *
 * For values var1 and var2, we can use a stride so that accesses for a given
 * element with id i will be of the form: var[i*stride]; this allows regular
 * access with stride 1, and access to constant variables stored as a
 * single-valued array with a stride of 0.
 *
 * parameters:
 *   property     <-- property queried
 *   n_vals       <-- number of values
 *   var1_stride  <-- stride between successive values of var1
 *   var2_stride  <-- stride between successive values of var2
 *   var1         <-- values on first plane axis
 *   var2         <-- values on second plane axis
 *   val          --> resulting property values
 *----------------------------------------------------------------------------*/

void
cs_phys_prop_compute(cs_phys_prop_type_t          property,
                     cs_lnum_t                    n_vals,
                     cs_lnum_t                    var1_stride,
                     cs_lnum_t                    var2_stride,
                     const cs_real_t              var1[],
                     const cs_real_t              var2[],
                     cs_real_t                    val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get reference value of a physical property
 *
 * \param[in] name  property name
 *
 * \return reference value (cs_real_t)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_physical_property_get_ref_value(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set reference value for a physical property
 *
 * \param[in] name  property name
 * \param[in] val   new value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_set_ref_value(const char       *name,
                                   const cs_real_t   val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get property reference values for a given zone
 *
 * \param[in] name    property name
 * \param[in] zname   zone name
 * \param[in] retval  array of values to return
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_get_zone_values(const char  *name,
                                     const char  *zname,
                                     cs_real_t    retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a physical property
 *
 * \param[in] name    property name
 * \param[in] dim     property dimension
 * \param[in] refval  reference value
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_create(const char      *name,
                            const int        dim,
                            const cs_real_t  refval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property definition on a given zone using a single value
 *
 * \param[in] name   property name
 * \param[in] zname  zone name
 * \param[in] dim    property dimension
 * \param[in] val    reference value for the zone
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_define_from_value(const char       *name,
                                       const char       *zname,
                                       const int         dim,
                                       const cs_real_t   val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property multi-diemnsional definition on a given zone
 *
 * \param[in] name   property name
 * \param[in] zname  zone name
 * \param[in] dim    property dimension (>1)
 * \param[in] vals   array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_define_from_values(const char  *name,
                                        const char  *zname,
                                        const int    dim,
                                        cs_real_t    vals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property definition based on a cs_field_t.
 *
 * The field is created if needed
 *
 * \param[in] name          property name
 * \param[in] type_flag     field type flag
 * \param[in] location_id   location id flag
 * \param[in] dim           field dimension
 * \param[in] has_previous  does the field has val_pre
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_define_from_field(const char  *name,
                                       int          type_flag,
                                       int          location_id,
                                       int          dim,
                                       bool         has_previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return id of field associated to property
 *
 * \param[in] name  property name
 *
 * \return field id (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_physical_property_field_id_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update reference values for a property on a given zone
 *
 * \param[in] name   property name
 * \param[in] zname  zone name
 * \param[in] vals   array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_update_zone_values(const char       *name,
                                        const char       *zname,
                                        const cs_real_t   vals[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PHYSICAL_PROPERTIES_H__ */
