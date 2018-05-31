#ifndef __CS_XDEF_H__
#define __CS_XDEF_H__

/*============================================================================
 * Routines to handle extended definitions of quantities (cs_xdef_t structures)
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <string.h>

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_param.h"
#include "cs_quadrature.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \enum cs_xdef_type_t
 *
 * \var CS_XDEF_BY_ANALYTIC_FUNCTION
 * Definition relying on an analytic function (see \ref cs_analytic_func_t)
 *
 * \var CS_XDEF_BY_ARRAY
 * Definition based on an array
 *
 * \var CS_XDEF_BY_FIELD
 * Definition based on a field (see \ref cs_field_t)
 *
 * \var CS_XDEF_BY_FUNCTION
 * Definition relying on a generic user-defined function. TODO
 *
 * \var CS_XDEF_BY_QOV
 * QOV = Quantity Over a Volume
 * Definition which enables to spread a given quantity inside a volume.
 * Useful to initialized a tracer in a subdomain for instance.
 *
 * \var CS_XDEF_BY_TIME_FUNCTION
 * Definition relying on a function for setting the time step (see
 * \ref cs_timestep_func_t)
 *
 * \var CS_XDEF_BY_VALUE
 * Simple definition by a constant value
 */

typedef enum {

  CS_XDEF_BY_ANALYTIC_FUNCTION,
  CS_XDEF_BY_ARRAY,
  CS_XDEF_BY_FIELD,
  CS_XDEF_BY_FUNCTION,
  CS_XDEF_BY_QOV,
  CS_XDEF_BY_TIME_FUNCTION,
  CS_XDEF_BY_VALUE,

  CS_N_XDEF_TYPES

} cs_xdef_type_t;

/*!
 * \enum cs_xdef_support_t
 *
 * \var CS_XDEF_SUPPORT_TIME
 * Definition for the time step. No zone is attached.
 *
 * \var CS_XDEF_SUPPORT_BOUNDARY
 * Definition for a boundary zone. Zones are attached to a list of boundary
 * faces.
 *
 * \var CS_XDEF_SUPPORT_VOLUME
 * Definition for a volumic zone. Zones are attached to a list of cells.
 */

typedef enum {

  CS_XDEF_SUPPORT_TIME,      /* support for time step description */
  CS_XDEF_SUPPORT_BOUNDARY,  /* zones attached to boundary faces */
  CS_XDEF_SUPPORT_VOLUME,

  CS_N_XDEF_SUPPORTS

} cs_xdef_support_t;

/*!
 * \struct cs_xdef_t
 *  \brief Structure storing medata for defining a quantity in a very flexible
 *         way
 */

typedef struct {

  /*! \var dim
   * dimension of the values attached to this description
   *
   * \var type
   * type of definition (see \ref cs_xdef_type_t)
   *
   * \var z_id
   * id related to a zone (volume or boundary) for this definition
   *
   * \var support
   * support for this definition (see \ref cs_xdef_support_t)
   *
   * \var state
   * Flag storing state of the values related to this definition
   * Example: steady, uniform, cellwise...
   *
   * \var meta
   * Flag storing in a condensed way metadata about the description.
   * These metadata may vary according to the object on which the description
   * applies.
   *
   * \var qtype
   * type of quadrature to use for evaluating the description (see
   * \ref cs_quadrature_type_t)
   *
   * \var input
   * Pointer to a structure cast on-the-fly according to the type of description
   * May be set to NULL or \ref cs_xdef_array_input_t or
   * \ref cs_xdef_analytic_input_t or \ref cs_xdef_timestep_input_t
   */

  int                    dim;
  cs_xdef_type_t         type;
  int                    z_id;
  cs_xdef_support_t      support;

  cs_flag_t              state;
  cs_flag_t              meta;

  cs_quadrature_type_t   qtype;

  void                  *input;

} cs_xdef_t;

/*!
 * \struct cs_xdef_array_input_t
 * \brief Input structure when an array is used for the definition
 */

typedef struct {

  /*! * \var stride
   * stride to access the array values
   *
   * \var loc
   * flag to know where are defined array values
   *
   * \var values
   * array values
   *
   * \var index
   * optional index for accessing to the values
   */

  int           stride;
  cs_flag_t     loc;
  cs_real_t    *values;
  cs_lnum_t    *index;

} cs_xdef_array_input_t;

/*!
 * \struct cs_xdef_analytic_input_t
 * \brief Input structure when an analytic function is used for the definition
 */

typedef struct {

  /*! \var input
   * NULL or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */
  void                *input;

  /*! \var func
   * \ref cs_analytic_func_t to call
   */
  cs_analytic_func_t  *func;

} cs_xdef_analytic_input_t;

/*!
 * \struct cs_xdef_timestep_input_t
 * \brief Input structure when a time step function is used for the definition
 */

typedef struct {

  /*! \var input
   * NULL or pointer to a structure cast on-the-fly for additional information
   * used in the function
   *
   * \var func
   * \ref cs_timestep_func_t to call
   */
  void                *input;
  cs_timestep_func_t  *func;

} cs_xdef_timestep_input_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the volume zone if from the zone name (If name = NULL or
 *         has an empty length, all entities are selected)
 *
 * \param[in] z_name            name of the zone
 *
 * \return the id of the related zone
 */
/*----------------------------------------------------------------------------*/

static inline int
cs_get_vol_zone_id(const char   *z_name)
{
  int z_id = 0;
  if (z_name != NULL) {
    if (strlen(z_name) > 0) {
      const cs_zone_t  *z = cs_volume_zone_by_name(z_name);
      z_id = z->id;
    }
  }
  return z_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the boundary zone if from the zone name (If name = NULL or
 *         has an empty length, all entities are selected)
 *
 * \param[in] z_name            name of the zone
 *
 * \return the id of the related zone
 */
/*----------------------------------------------------------------------------*/

static inline int
cs_get_bdy_zone_id(const char   *z_name)
{
  int z_id = 0;
  if (z_name != NULL) {
    if (strlen(z_name) > 0) {
      const cs_zone_t  *z = cs_boundary_zone_by_name(z_name);
      z_id = z->id;
    }
  }
  return z_id;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on volumic
 *         elements
 *
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_volume_create(cs_xdef_type_t    type,
                      int               dim,
                      int               z_id,
                      cs_flag_t         state,
                      cs_flag_t         meta,
                      void             *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on boundary
 *         elements
 *
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_boundary_create(cs_xdef_type_t    type,
                        int               dim,
                        int               z_id,
                        cs_flag_t         state,
                        cs_flag_t         meta,
                        void             *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure for setting the
 *         time step
 *
 * \param[in]  type       type of definition
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure storing the parameters (cast
 *                        on-the-fly according to the type of definition)
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_timestep_create(cs_xdef_type_t             type,
                        cs_flag_t                  state,
                        cs_flag_t                  meta,
                        void                      *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_xdef_t structure
 *
 * \param[in, out] d    pointer to a cs_xdef_t structure
 *
 * \return NULL
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_free(cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  copy a cs_xdef_t structure
 *
 * \param[in]  src    pointer to a cs_xdef_t structure to copy
 *
 * \return a pointer to a new allocated cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_copy(cs_xdef_t     *src);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the array after having added
 *         this definition
 *
 * \param[in, out]  d          pointer to a cs_xdef_t structure
 * \param[in]       is_owner   manage or not the lifecycle of the array values
 * \param[in]       array      values
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array(cs_xdef_t     *d,
                  bool           is_owner,
                  cs_real_t     *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the index to get access to the
 *         array values.
 *
 * \param[in, out]  d             pointer to a cs_xdef_t structure
 * \param[in]       array_index   index on array values
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array_index(cs_xdef_t     *d,
                        cs_lnum_t     *array_index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for evaluating the given
 *         description
 *
 * \param[in, out]  d       pointer to a cs_xdef_t structure
 * \param[in]       qtype   type of quadrature
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_quadrature(cs_xdef_t              *d,
                       cs_quadrature_type_t    qtype);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the type of quadrature to use for evaluating the given
 *         description
 *
 * \param[in]  d       pointer to a cs_xdef_t structure
 *
 * \return the type of quadrature
 */
/*----------------------------------------------------------------------------*/

cs_quadrature_type_t
cs_xdef_get_quadrature(cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the flag dedicated to the state
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 *
 * \return the value of the flag
 */
/*----------------------------------------------------------------------------*/

cs_xdef_type_t
cs_xdef_get_type(const cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the flag dedicated to the state
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 *
 * \return the value of the flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_xdef_get_state_flag(const cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a a cs_xdef_t structure
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log(cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_H__ */
