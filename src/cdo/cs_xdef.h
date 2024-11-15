#ifndef __CS_XDEF_H__
#define __CS_XDEF_H__

/*============================================================================
 * Functions to handle extended definitions of quantities
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_mesh_adjacencies.h"
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
 * \return a null pointer
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
   * May be set to null or \ref cs_xdef_array_context_t or
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
   * id related to a zone (volume or boundary).
   * If id = 0, then all cells (in case of volume zone) or all boundary faces
   * (in case of boundary zone) are selected. A full length array is thus
   * considered and the way to apply a definition by array is simpler.
   *
   * \var stride
   * Stride to access the array values
   *
   * \var value_location
   * Flag to know where are defined array values
   *
   * \var values
   * Array values
   *
   * \var is_owner
   * If true the lifecycle of the values is managed by the cs_xdef_t structure.
   * Otherwise, the lifecycle is managed by the calling code.
   *
   * \var full_length
   * The array describes only a part of the support. To know which part is
   * defined, one can relies on the elements of the zone when the support flag
   * corresponds to the (primal) cells or to the boundary faces. In other
   * cases, one needs the number of elements and its associated list.
   *
   * \var full2subset
   * (Optional) Array of size equal to the full support to get the position in
   * the subset list of an element. Allocated only if full_length is set to
   * false and only if needed.
   *
   * \var n_list_elts
   * (Optional) Number of element in the (sub)list of elements when the array
   * describes only a part of the full-length array (Case of value_location
   * which is neither the cells nor the boundary faces).
   *
   * \var elt_ids
   * (Optional) List of element ids. Useful when the array describes only a
   * part of the full-length array and the value location is neither the cells
   * for a volume definition nor the boundary faces for a boundary
   * definition. One assumes that the lifecycle of this array is managed
   * outside.
   *
   * \var adjacency
   * (Optional) Pointer to a shared adjacency structure (an indexed list). This
   * structure can be useful to manipulate arrays with advanced value location
   * (i.e. not the classical ones as vertices, cells or boundary faces). One
   * assumes that the lifecycle of this buffer is managed outside (pointer to a
   * cs_adjacency_t stored either in the \ref cs_cdo_connect_t struct. or the
   * \ref cs_mesh_t struct. for instance)
   */

  int                     z_id;
  int                     stride;
  cs_flag_t               value_location;
  bool                    is_owner;
  bool                    full_length;

  cs_real_t              *values;

  /* Automatic parameter */

  cs_lnum_t              *full2subset;

  /* Optional parameters */

  cs_lnum_t               n_list_elts;
  const cs_lnum_t        *elt_ids;
  const cs_adjacency_t   *adjacency;

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
   * null or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */

  void                  *input;

  /*! \var free_input
   * null or pointer to a function to free a given input structure
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

  /*! \var dof_location
   *  Flag to know which type of entities are given as parameter in the
   *  \ref cs_dof_func_t
   */

  cs_flag_t              dof_location;

  /*! \var func
   * pointer to a \ref cs_dof_func_t to call
   */

  cs_dof_func_t         *func;

  /*! \var input
   * null or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */

  void                  *input;

  /*! \var free_input
   * null or pointer to a function to free a given input structure
   */

  cs_xdef_free_input_t  *free_input;

} cs_xdef_dof_context_t;

/*!
 * \struct cs_xdef_time_func_context_t
 * \brief Context structure when a time step function is used for the definition
 */

typedef struct {

  /*! \var z_id
   * id related to a zone (volume or boundary) for this definition
   */

  int                    z_id;

  /*! \var func
   * pointer to a \ref cs_time_func_t to call
   */

  cs_time_func_t        *func;

  /*! \var input
   * null or pointer to a structure cast on-the-fly for additional information
   * used in the function
   */

  void                  *input;

  /*! \var free_input
   * null or pointer to a function to free a given input structure
   */

  cs_xdef_free_input_t  *free_input;

} cs_xdef_time_func_context_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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

static inline cs_real_t
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

static inline cs_real_t *
cs_xdef_array_get_values(const cs_xdef_t     *def)
{
  if (def == NULL)
    return NULL;

  assert(def->type == CS_XDEF_BY_ARRAY);

  cs_xdef_array_context_t  *ai = (cs_xdef_array_context_t *)def->context;

  return ai->values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief When the definition relies on a cs_field_t structure, return the
 *        pointer to the field structure
 *
 * \param[in]  def    pointer to a cs_xdef_t structure
 *
 * \return the pointer to the field structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_field_t *
cs_xdef_field_get(cs_xdef_t     *def)
{
  if (def == NULL)
    return NULL;

  assert(def->type == CS_XDEF_BY_FIELD);

  return (cs_field_t *)def->context;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize a new cs_xdef_t structure based on volumic
 *        elements
 *
 * \param[in] type        type of definition
 * \param[in] dim         dimension of the values to define
 * \param[in] z_id        volume zone id
 * \param[in] state       flag to know if this uniform, cellwise, steady...
 * \param[in] meta        metadata associated to this description
 * \param[in] context     pointer to a structure
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
 * \brief Allocate and initialize a new cs_xdef_t structure based on boundary
 *        elements
 *
 * \param[in] type       type of definition
 * \param[in] dim        dimension of the values to define
 * \param[in] z_id       volume zone id
 * \param[in] state      flag to know if this uniform, cellwise, steady...
 * \param[in] meta       metadata associated to this description
 * \param[in] context    pointer to a structure
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
 * \brief Allocate and initialize a new cs_xdef_t structure for setting the
 *        time step
 *
 * \param[in] type       type of definition
 * \param[in] state      flag to know if this uniform, cellwise, steady...
 * \param[in] meta       metadata associated to this description
 * \param[in] context    pointer to a structure storing the parameters (cast
 *                       on-the-fly according to the type of definition)
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
 * \return null pointer
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
 *        allows one to set a more or less complex input data structure. This
 *        call should be done before the first evaluation of the associated
 *        cs_xdef_t structure.
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
/*!
 * \brief  In case of definition by array, set the array values after having
 *         added this definition
 *
 * \param[in, out]  d          pointer to a cs_xdef_t structure
 * \param[in]       is_owner   manage or not the lifecycle of the array values
 * \param[in]       values     array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_array_set_values(cs_xdef_t     *d,
                         bool           is_owner,
                         cs_real_t     *values);

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
cs_xdef_array_set_zone_id(cs_xdef_t     *d,
                          int            z_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, build the full2subset array.
 *         The direct members of the cs_xdef_t structure are not modified but
 *         the context dedicated to definition by array is updated.
 *         d is declared as const to avoid a compiler warning
 *
 * \param[in, out]  d      pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_array_build_full2subset(const cs_xdef_t         *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the optional adjacency structure
 *
 * \param[in, out]  d      pointer to a cs_xdef_t structure
 * \param[in]       adj    pointer to the adjacency structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_array_set_adjacency(cs_xdef_t             *d,
                            const cs_adjacency_t  *adj);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the optional sub-list of
 *         elements used to link elements in the partial view and in the
 *         full-length view
 *
 * \param[in, out]  d        pointer to a cs_xdef_t structure
 * \param[in]       n_elts   number of elements in the sub-list
 * \param[in]       elt_ids  list of element ids constituting the sub-list
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_array_set_sublist(cs_xdef_t         *d,
                          cs_lnum_t          n_elts,
                          const cs_lnum_t    elt_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the current field values in case of definition by field
 *
 * \param[in]  def    pointer to a cs_xdef_t structure
 *
 * \return the pointer to the current field values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_xdef_field_get_values(cs_xdef_t     *def);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_H__ */
