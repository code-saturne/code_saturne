#ifndef __CS_XDEF_H__
#define __CS_XDEF_H__

/*============================================================================
 * Functions to handle extended definitions of quantities
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <string.h>

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_param_types.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy an input data structure.
 *         Complex data structure can be used when a \ref cs_xdef_t structure
 *         is defined by an analytic function, a DoF function or a time
 *         function. Please refer to \ref cs_xdef_analytic_context_t,
 *         \ref cs_xdef_time_func_context_t or \ref cs_xdef_dof_context_t
 *
 * \param[in, out]  input   pointer to an input structure associated to a
 *                          context structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_xdef_free_input_t)(void   *input);

/*!
 * \enum cs_xdef_type_t
 *
 * \var CS_XDEF_BY_ANALYTIC_FUNCTION
 * Definition relying on a \ref cs_analytic_func_t function pointer
 *
 * \var CS_XDEF_BY_ARRAY
 * Definition based on an array
 *
 * \var CS_XDEF_BY_DOF_FUNCTION
 * Definition relying on a \ref cs_dof_func_t function pointer
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
 * \var CS_XDEF_BY_SUB_DEFINITIONS
 * Definition relying on a combination of other CS_XDEF_***
 * This kind of definition is useful for the definition of a property
 * as the product of two existing ones.
 *
 * \var CS_XDEF_BY_TIME_FUNCTION
 * Definition relying on a function for setting the time step (see
 * \ref cs_time_func_t)
 *
 * \var CS_XDEF_BY_VALUE
 * Simple definition by a constant value
 */

typedef enum {

  CS_XDEF_BY_ANALYTIC_FUNCTION,
  CS_XDEF_BY_ARRAY,
  CS_XDEF_BY_DOF_FUNCTION,
  CS_XDEF_BY_FIELD,
  CS_XDEF_BY_FUNCTION,
  CS_XDEF_BY_QOV,
  CS_XDEF_BY_SUB_DEFINITIONS,
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
   * \var context
   * Pointer to a structure cast on-the-fly according to the type of description
   * May be set to NULL or \ref cs_xdef_array_context_t or
   * \ref cs_xdef_analytic_context_t or \ref cs_xdef_time_func_context_t or
   * \ref cs_xdef_dof_context_t
   */

  int                    dim;
  cs_xdef_type_t         type;
  int                    z_id;
  cs_xdef_support_t      support;

  cs_flag_t              state;
  cs_flag_t              meta;

  cs_quadrature_type_t   qtype;

  void                  *context;

} cs_xdef_t;

/*!
 * \struct cs_xdef_array_context_t
 * \brief Context structure when an array is used for the definition
 */

typedef struct {

  /*!
   * \var z_id
   * id related to a zone (volume or boundary) used for the size of the array.
   *
   * \var stride
   * Stride to access the array values
   *
   * \var loc
   * Flag to know where are defined array values
   *
   * \var values
   * Array values
   *
   * \var is_owner
   * If true the lifecycle of the values is managed by the cs_xdef_t structure.
   * Otherwise, the lifecycle is managed by the calling code.
   *
   * \var index
   * Optional index for accessing to the values. (shared pointer => One assumes
   * that the lifecycle of this buffer is managed outside (pointer to a
   * cs_adjacency_t stored either in the \ref cs_cdo_connect_t struct. or the
   * \ref cs_mesh_t struct.
   *
   * \var ids
   * Optional list of entity ids (shared pointer)
   * This can be either the list of ids associated to the given index (case of
   * a \ref cs_adjacency_t struct.) or simply an indirection list if index is
   * set to NULL
   */

  int                 z_id;
  int                 stride;
  cs_flag_t           loc;
  cs_real_t          *values;
  bool                is_owner;

  const cs_lnum_t    *index;
  const cs_lnum_t    *ids;

} cs_xdef_array_context_t;

/*!
 * \struct cs_xdef_analytic_context_t
 * \brief Context structure when a definition by analytic function is used
 */

typedef struct {

  /*! \var z_id
   * id related to a zone (volume or boundary) for this definition
   */

  int                    z_id;

  /*! \var func
   * pointer to a \ref cs_analytic_func_t to call
   */

  cs_analytic_func_t    *func;

  /*! \var input
   * NULL or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */

  void                  *input;

  /*! \var free_input
   * NULL or pointer to a function to free a given input structure
   */

  cs_xdef_free_input_t  *free_input;

} cs_xdef_analytic_context_t;

/*!
 * \struct cs_xdef_dof_context_t
 * \brief Context structure when a definition by DoF function is used
 */

typedef struct {

  /*! \var z_id
   * id related to a zone (volume or boundary) for this definition
   */

  int                    z_id;

  /*! \var loc
   *  Flag to know which type of entities are given as parameter in the
   *  \ref cs_dof_func_t
   */

  cs_flag_t              loc;

  /*! \var func
   * pointer to a \ref cs_dof_func_t to call
   */

  cs_dof_func_t         *func;

  /*! \var input
   * NULL or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */
  void                  *input;

  /*! \var free_input
   * NULL or pointer to a function to free a given input structure
   */

  cs_xdef_free_input_t  *free_input;

} cs_xdef_dof_context_t;

/*!
 * \struct cs_xdef_time_func_context_t
 * \brief Context structure when a time step function is used for the definition
 */

typedef struct {

  /*! \var func
   * pointer to a \ref cs_time_func_t to call
   */

  cs_time_func_t        *func;

  /*! \var input
   * NULL or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */

  void                  *input;

  /*! \var free_input
   * NULL or pointer to a function to free a given input structure
   */

  cs_xdef_free_input_t  *free_input;

} cs_xdef_time_func_context_t;

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the value associated to the given definition.
 *         This should be a definition by value and the dimension should be
 *         equal to one.
 *
 * \param[in]  def    pointer to a cs_xdef_t structure
 *
 * \return the value of the definition
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
cs_xdef_get_scalar_value(cs_xdef_t     *def)
{
  assert(def != NULL);
  assert(def->dim == 1);
  assert(def->type == CS_XDEF_BY_VALUE);

  cs_real_t  *value = (cs_real_t *)def->context;

  return value[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the values associated to the given definition.
 *         This should be a definition by array
 *
 * \param[in]  def    pointer to a cs_xdef_t structure
 *
 * \return the pointer to the array of values
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t *
cs_xdef_get_array(cs_xdef_t     *def)
{
  assert(def != NULL);
  assert(def->type == CS_XDEF_BY_ARRAY);

  cs_xdef_array_context_t  *ai = (cs_xdef_array_context_t *)def->context;

  return ai->values;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on volumic
 *         elements
 *
 * \param[in]  type        type of definition
 * \param[in]  dim         dimension of the values to define
 * \param[in]  z_id        volume zone id
 * \param[in]  state       flag to know if this uniform, cellwise, steady...
 * \param[in]  meta        metadata associated to this description
 * \param[in]  context     pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_volume_create(cs_xdef_type_t           type,
                      int                      dim,
                      int                      z_id,
                      cs_flag_t                state,
                      cs_flag_t                meta,
                      void                    *context);

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
 * \param[in]  context    pointer to a structure
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
                        void             *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure for setting the
 *         time step
 *
 * \param[in]  type       type of definition
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  context    pointer to a structure storing the parameters (cast
 *                        on-the-fly according to the type of definition)
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_timestep_create(cs_xdef_type_t       type,
                        cs_flag_t            state,
                        cs_flag_t            meta,
                        void                *context);

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
 * \brief In the case of a definition by an analytic function, a time function
 *        or a function relying on degrees of freedom (DoFs), this function
 *        allows one to set a more or less complex input data structure.  This
 *        call should be done before the first evaluation call of the
 *        associated cs_xdef_t structure.
 *
 * \param[in, out]  d         pointer to a cs_xdef_t structure
 * \param[in]       input     pointer to an input structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_input_context(cs_xdef_t       *d,
                          void            *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief In case of a definition by an analytic function, a time function or a
 *        function relying on degrees of freedom (DoFs). One can set a function
 *        to free a complex input data structure (please refer to \ref
 *        cs_xdef_free_input_t) for more details.
 *
 * \param[in, out]  d             pointer to a cs_xdef_t structure
 * \param[in]       free_input    pointer to a function which free the input
 *                                structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_free_input_function(cs_xdef_t               *d,
                                cs_xdef_free_input_t    *free_input);

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
 * \brief  In case of definition by array, set the zone id related to the size
 *         of the array. By default, the zone id is the same as the zone id
 *         related to the definition so that there is no need to call this
 *         function.
 *
 * \param[in, out]  d       pointer to a cs_xdef_t structure
 * \param[in]       z_id    zone id associated to the array size
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array_zone_id(cs_xdef_t     *d,
                          int            z_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the optional index and ids
 *         arrays that may be useful when operating on definitions by array
 *
 * \param[in, out]  d         pointer to a cs_xdef_t structure
 * \param[in]       index     optional pointer to an array of index values
 * \param[in]       ids       optional pointer to a list of entity ids
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array_pointers(cs_xdef_t            *d,
                           const cs_lnum_t      *index,
                           const cs_lnum_t      *ids);

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
 * \brief  Output the settings related to a cs_xdef_t structure in the setup
 *         logging file
 *
 * \param[in] prefix    optional string
 * \param[in] d         pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log_setup(const char          *prefix,
                  const cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a cs_xdef_t structure
 *
 * \param[in] log_type  related log file to consider
 * \param[in] prefix    optional string
 * \param[in] d         pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log(cs_log_t             log_type,
            const char          *prefix,
            const cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the cs_xdef_type's name string
 *
 * \param[in] xdef_type  type to query
 *
 * \return a pointer to mathing name string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_xdef_type_get_name(cs_xdef_type_t  xdef_type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_H__ */
