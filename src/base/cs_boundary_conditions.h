#ifndef __CS_BOUNDARY_CONDITIONS_H__
#define __CS_BOUNDARY_CONDITIONS_H__

/*============================================================================
 * Boundary condition handling.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

#include "cs_base.h"
#include "cs_field.h"
#include "cs_function.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_time_control.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*! Maximum number of physical model zones
 *
 * \deprecated This is used for Fortran compatibilty (and maps to the Fortran
 * \ref ppppar:nozppm  "nozppm" parameter). In C, we should move to high
 * level boundary condition definitions not requiring indexing by legacy
 * zone numbers.
*/
#define  CS_MAX_BC_PM_ZONE_NUM  2000

/*! Maximum number of boundary condition code types */
#define  CS_MAX_BC_TYPE  200

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*! Boundary condition type (code) associated with each boundary face */

extern const int  *cs_glob_bc_type;

/*! boundary zone number associated with each boundary face
 *  (specific physical models)
 *
 * \deprecated This is used for \ref cs_boundary_condition_pm_info_t only.
*/

/*----------------------------------------------------------------------------*/

/*! Legacy physical model boundary conditions.
 *
 * \remark The amppings of member arrays of this structure are shifted
 * by 1 when mapped to Fortran, so that both in Fortran and C, we can use
 * the natural indexing (1-based and 0-based respectively) without needing
 * to shift the zone id/number.
 *
 * \deprecated This should be used for Fortran compatibilty and migration
 * to C only. In C, we should then move to high level boundary condition
 * definitions not requiring indexing by legacy zone numbers indexed by
 * \ref cs_glob_bc_face_zone. */

typedef struct {

  /*! Legacy physical model zone id per boundary face */

  int  *izfppp;

  /*! indirection array allowing to sort the boundary faces
   *  according to their boundary condition type \c bc_type */
  int  *itrifb;

  /*! Imposed flow zone indicator (for inlet zones).
   * If the mass flow is imposed (\c iqimp(z_id) = 1), the matching
   * \c qimp value must be set, and the defined velocity boundary condition
   * will be rescaled so as to match the given mass flow (i.e. only its original
   * direction is used. Otherwise, the given velocity boundary condition
   * given by \c rcodcl1 is unchanged. */
  int   iqimp[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Turbulence inlet type:
    * - 0: given by the user
    * - 1: automatic, from hydraulic diameter and input velocity performed.
    * - 2: automatic, from turbulent intensity and input velocity performed.
    */
  int   icalke[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Imposed flow value (for inlet zones).
   * If the mass flow is imposed (\c iqimp(z_num - 1) = 1), the matching \c qimp
   * value must be set, and the defined velocity boundary condition will be
   * rescaled so as to match the given mass flow (i.e. only its original
   * direction is used. Otherwise, the given velocity boundary condition
   * given by \c rcodcl1 is unchanged. */
  cs_real_t  qimp[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! hydraulic diameter */
  cs_real_t  dh[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! turbulent intensity */
  cs_real_t  xintur[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! gas combustion (cogz) */
  int ientfu[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for fuel flow inlet

  int ientox[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for an air fow inlet

  int ientgb[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for burned gas inlet

  int ientgf[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for unburned gas inlet

  /*! inlet temperature (gas combustion) */
  double tkent[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Mean Mixture Fraction at Inlet (gas combustion) */
  double fment[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! atmo */
  /* atmospheric flows: auto inlet/outlet flag */
  int *iautom;

  /* atmospheric flows: on/off for profile from data */
  int iprofm[CS_MAX_BC_PM_ZONE_NUM+1];

} cs_boundary_condition_pm_info_t;

extern cs_boundary_condition_pm_info_t  *cs_glob_bc_pm_info;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Handling of boundary condition definition errors and associated output.
 *
 * This function checks for errors, and simply returns if no error is
 * encountered. In case of error, it outputs helpful information so as to
 * make it easier to locate the matching faces.
 *
 * For each boundary face, bc_type defines the boundary condition type.
 * As a convention here, zero values correspond to undefined types,
 * positive values to defined types (with no error), and negative values
 * to defined types with inconsistent or incompatible values, the
 * absolute value indicating the original boundary condition type.
 *
 * An optional label may be used if the error is related to another
 * attribute than the boundary type, for appropriate error reporting.
 *
 * parameters:
 *   bc_flag   <-- array of BC type ids
 *   type_name <-- name of attribute in error, or NULL
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const int   *bc_type,
                             const char  *type_name);

/*----------------------------------------------------------------------------
 * Locate shifted boundary face coordinates on possibly filtered
 * cells or boundary faces for later interpolation.
 *
 * parameters:
 *   location_type   <-- matching values location (CS_MESH_LOCATION_CELLS
 *                        or CS_MESH_LOCATION_BOUNDARY_FACES)
 *   n_location_elts <-- number of selected location elements
 *   n_faces         <-- number of selected boundary faces
 *   location_elts   <-- list of selected location elements (0 to n-1),
 *                       or NULL if no indirection is needed
 *   faces           <-- list of selected boundary faces (0 to n-1),
 *                       or NULL if no indirection is needed
 *   coord_shift     <-- array of coordinates shift relative to selected
 *                       boundary faces
 *   coord_stride    <-- access stride in coord_shift: 0 for uniform
 *                       shift, 1 for "per face" shift.
 *   tolerance       <-- relative tolerance for point location.
 *
 * returns:
 *   associated locator structure
 *----------------------------------------------------------------------------*/

ple_locator_t *
cs_boundary_conditions_map(cs_mesh_location_type_t    location_type,
                           cs_lnum_t                  n_location_elts,
                           cs_lnum_t                  n_faces,
                           const cs_lnum_t           *location_elts,
                           const cs_lnum_t           *faces,
                           cs_real_3_t               *coord_shift,
                           int                        coord_stride,
                           double                     tolerance);

/*----------------------------------------------------------------------------
 * Set mapped boundary conditions for a given field and mapping locator.
 *
 * parameters:
 *   field           <-- field whose boundary conditions are set
 *   locator         <-- associated mapping locator, as returned
 *                       by cs_boundary_conditions_map().
 *   location_type   <-- matching values location (CS_MESH_LOCATION_CELLS or
 *                       CS_MESH_LOCATION_BOUNDARY_FACES)
 *   normalize       <-- normalization option:
 *                         0: values are simply mapped
 *                         1: values are mapped, then multiplied
 *                            by a constant factor so that their
 *                            surface integral on selected faces
 *                            is preserved (relative to the
 *                            input values)
 *                         2: as 1, but with a boundary-defined
 *                            weight, defined by balance_w
 *                         3: as 1, but with a cell-defined
 *                            weight, defined by balance_w
 *   interpolate     <-- interpolation option:
 *                         0: values are simply based on matching
 *                            cell or face center values
 *                         1: values are based on matching cell
 *                            or face center values, corrected
 *                            by gradient interpolation
 *   n_faces         <-- number of selected boundary faces
 *   faces           <-- list of selected boundary faces (0 to n-1),
 *                       or NULL if no indirection is needed
 *   balance_w       <-- optional balance weight, or NULL
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_mapped_set(const cs_field_t          *f,
                                  ple_locator_t             *locator,
                                  cs_mesh_location_type_t    location_type,
                                  int                        normalize,
                                  int                        interpolate,
                                  cs_lnum_t                  n_faces,
                                  const cs_lnum_t           *faces,
                                  cs_real_t                 *balance_w);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add location of locate shifted boundary face coordinates on
 *        cells or boundary faces for automatic interpolation.
 *
 * \note
 * This function is currently restricted to mapping of boundary face
 * locations (usually from boundary zones) to cell of boundary face
 * locations, but could be extended to other location types in the future.
 *
 * \param[in]  bc_location_id      id of selected boundary mesh location;
 *                                 currently restricted to subsets of
 *                                 boundary faces (i.e. boundary zone
 *                                 location ids).
 * \param[in]  source_location_id  id of selected location  mesh location
 *                                 (usually CS_MESH_LOCATION_CELLS but can be
 *                                 a more restricted cell or boundary face zone
 *                                 location location id).
 * \param[in]  coord_shift      coordinates shift relative to selected
 *                              boundary faces
 * \param[in]  tolerance        relative tolerance for point location.
 *
 * \return  id of added map
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_conditions_add_map(int         bc_location_id,
                               int         source_location_id,
                               cs_real_t   coord_shift[3],
                               double      tolerance);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the legacy boundary conditions zone data arrays
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create_legacy_zone_data(void);

/*----------------------------------------------------------------------------
 * Create the boundary conditions face type and face zone arrays
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create(void);

/*----------------------------------------------------------------------------
 * Free the boundary conditions face type and face zone arrays.
 *
 * This also frees boundary condition mappings which may have been defined.
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare (reset) condition coefficients for all variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_reset(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update per variable boundary condition codes.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_compute(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define automatic turbulence values for specific physical modules.
 *
 * The definitions are similar to those of the standard case, though wall
 * shear direction is not computed for second-order models, and determination
 * of face BC types is done using the legacy physical model zone info
 * (cs_glob_bc_pm_info->izfpp, ...).
 *
 * \deprecated  Code should migrate to the "per zone" open boundary condition
 * definitions.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_legacy_turbulence(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic adjustments for boundary condition codes.
 *
 * Currently handles mapped inlets, after the call to \ref stdtcl.
 * As portions of stdtcl are migrated to C, they should be called here,
 * before mapped inlets.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_complete(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to model inlet structure associated with a given
 *        open (inlet/outlet) boundary.
 *
 * The returned pointer is of type void * as it should be cast to the
 * appropriate (model-dependent) type.

 * If no matching parent open boundary has been created yet, it is created.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

void *
cs_boundary_conditions_get_model_inlet(const cs_zone_t  *zone);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return legacy zone number related to a given zone, if available.
 *
 * \param[in]  z  pointer to associated zone
 *
 * \return  number associated with legacy zone, or 0 if unavailable.
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_conditions_get_legacy_zone_num(const  cs_zone_t  *z);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to boundary head losses array.
 *
 * The array is allocated if not previously available.
 *
 * \param[in]  alloc_if_null   do we need to allocate this if not present ?

 * \return  b_head_loss  pointer to boundary head losses array, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_boundary_conditions_get_b_head_loss(bool  alloc_if_null);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign pointer to model inlet structure associated with a given
 *        open (inlet/outlet) boundary.
 *
 * The returned pointer is of type void * as it should be cast to the
 * appropriate (model-dependent) type.

 * If no matching parent open boundary has been created yet, it is created.
 *
 * \param[in]  zone   pointer to associated zone
 * \param[in]  s_ptr  pointer to associated structure
 * \param[in]  s_del  destructor for associated structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_assign_model_inlet(const cs_zone_t  *zone,
                                          void             *s_ptr,
                                          void             *s_del);

/*----------------------------------------------------------------------------*/
/*
 * \brief Acess the time control structure of an inlet.
 *
 * This allows modifying that structure, for example updating the inlet
 * velocity values only in a certain time range, and avoiding
 * uneeded recomputations outside that range.
 *
 * \param[in]  zone  pointer to associated zone
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t *
cs_boundary_conditions_open_get_time_control(const  cs_zone_t  *zone);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant velocity to an open (inlet/outlet) boundary.
 *
 * This function may also be used to define the flow direction if called
 * before one of the \c cs_boundary_conditions_open_set_mass_flow_rate
 * or \c cs_boundary_conditions_open_set_volume_flow_rate functions.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  u  associated velocity value
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_value(const cs_zone_t  *z,
                                                  const cs_real_t   u[3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant velocity normal to an inlet.
 *
 * \param[in]  z       pointer to associated zone
 * \param[in]  u_norm  associated constant normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_normal_value(const  cs_zone_t  *z,
                                                         cs_real_t     u_norm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a normal velocity to an inlet using a provided function.
 *
 * Reminder: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated velocity vector evaluation function
 *                    at zone faces
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_func(const  cs_zone_t       *z,
                                                 cs_eval_at_location_t  *func,
                                                 void                   *input);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return the volume flow rate to an inlet or outlet.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * \param[in]  z  pointer to associated zone
 *
 * \return  volume flow rate associated with open boundary
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_boundary_conditions_open_get_mass_flow_rate(const  cs_zone_t  *z);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant mass flow rate to an inlet or outlet.
 *
 * By default, the flow direction is considered normal to the boundary.
 * The flow direction may be specified by calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func
 * for the appropriate zone before calling this function.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * mass flow rate.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_mass_flow_rate_by_value(const  cs_zone_t  *z,
                                                        cs_real_t          q);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a mass flow rate to an inlet or outlet  based on
 *        provided function.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * \ref cs_boundary_conditions_open_set_velocity_by_normal_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * mass flow rate.
 *
 * Since the flow rate is a global value, the provided function should
 * be associated with the CS_MESH_LOCATION_NONE location.
 *
 * Note also that during updates, this function will be called before
 * the velocity vector update, so in complex cases where flow rate computation
 * would require feedback from the velocity at this boundary, the user
 * must be aware that values from the previous time step or update will
 * be used, handle this in another manner.
 *
 * Reminder: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (mass flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_mass_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant volume flow rate to an inlet or outlet.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant volume flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_volume_flow_rate_by_value(const  cs_zone_t  *z,
                                                          cs_real_t          q);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a volume flow rate to an inlet or outlet based on
 *        provided function.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * \ref cs_boundary_conditions_open_set_velocity_by_normal_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * Since the flow rate is a global value, the provided function should
 * be associated with the CS_MESH_LOCATION_NONE location.
 *
 * Note also that during updates, this function will be called before
 * the velocity vector update, so in complex cases where flow rate computation
 * would require feedback from the velocity at this boundary, the user
 * must be aware that values from the previous time step or update will
 * be used, handle this in another manner.
 *
 * Reminder: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (volume flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_volume_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input);

/*----------------------------------------------------------------------------*/
/*
 * \brief Base the inlet turbulence values on a a circular duct with smooth
 *        wall (see ref cs_turbulence_bc_ke_hyd_diam).
 *
 * \param[in]  zone  pointer to associated zone
 * \param[in]  hd    associated hydraulic diameter
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_inlet_set_turbulence_hyd_diam(const  cs_zone_t  *zone,
                                                     cs_real_t          hd);

/*----------------------------------------------------------------------------*/
/*
 * \brief Base the inlet turbulence values on a a circular duct with smooth
 *        wall (see ref cs_turbulence_bc_ke_hyd_diam).
 *
 * \param[in]  zone  pointer to associated zone
 * \param[in]  ti    associated turbulence intensity
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_inlet_set_turbulence_intensity(const  cs_zone_t  *zone,
                                                      cs_real_t          ti);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_H__ */
