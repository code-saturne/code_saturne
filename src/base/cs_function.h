#ifndef __CS_FUNCTION_H__
#define __CS_FUNCTION_H__

/*============================================================================
 * Function objects management.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_mesh_location.h"
#include "cs_param_types.h"
#include "cs_restart_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @defgroup function_flags Flags specifying general function attributes
 *
 * @{
 */

/*
 * Function property type
 */

/*! represents an intensive quantity */
#define CS_FUNCTION_INTENSIVE           (1 << 0)

/*! represents an extensive quantity */
#define CS_FUNCTION_EXTENSIVE           (1 << 1)

/*! represents a quantity which does not change over time */
#define CS_FUNCTION_TIME_INDEPENDENT    (1 << 2)

/*! user-defined */
#define CS_FUNCTION_USER                (1 << 3)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for evaluation of local function values.
 *
 * If the matching values are multidimensional, they must be interleaved.
 * The output values are assumed to use a dense storage (i.e. of size
 * \c n_elts * <value dimension> for the associated data type, in the same
 * order as \c elt_ids if present.)
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to optional (untyped) value
 *                               or structure.
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_eval_at_location_t) (int               location_id,
                         cs_lnum_t         n_elts,
                         const cs_lnum_t  *elt_ids,
                         void             *input,
                         void             *vals);

/* Function descriptor */
/*---------------------*/

typedef struct {

  const char             *name;         /*!< Canonical name */
  char                   *label;        /*!< Optional label (if NULL,
                                          name is used instead) */

  const int               id;           /*!< Function id */
  int                     type;         /*!< Function type flag */

  const int               dim;          /*!< Number of values per location
                                          element (i.e. field dimension) */

  const int               location_id;  /*!< Id of matching mesh location */

  const cs_datatype_t     datatype;     /*!< Associated data type */

  int                     post_vis;     /*!< postprocessing/visualization
                                          flag (same usage as cs_field_t
                                          keyword of the same name) */

  int                     log;          /*!< postprocessing/visualization
                                          flag (same usage as cs_field_t
                                          keyword of the same name) */

  int                     time_stamp;   /*!< Time stamp at last evaluation,
                                          or -1 if unused */

  cs_restart_file_t       restart_file;   /*!< type of associated checkpoint
                                            file if evaluation should be
                                            saved to checkpoint (for
                                            later postprocessing) */

  cs_eval_at_location_t  *eval_func;      /*!< Associated data evaluation
                                            function, or NULL */
  cs_analytic_func_t     *analytic_func;  /*!< Associated data evaluation
                                            function, or NULL */
  cs_dof_func_t          *dof_func;       /*!< Associated data evaluation
                                            function, or NULL */

  void                   *func_input;   /* Pointer to optional (untyped)
                                           value or structure, for use by
                                           evaluation function */

} cs_function_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a function whose data values will be computed using the
 *        provided evaluation function.
 *
 * If of dimension > 1, the evaluated values are always interleaved.
 *
 * \param[in]  name          name of associated function
 * \param[in]  location_id   id of associated mesh location
 * \param[in]  dim           dimension associated with element data
 * \param[in]  is_intensive  is the function intensive?
 * \param[in]  datatype      associated data values type
 * \param[in]  data_func     function used to define data values
 * \param[in]  data_input    pointer to optional (untyped) value or structure
 *                           to be used by data_func

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_by_func(const char             *name,
                           int                     location_id,
                           int                     dim,
                           bool                    is_intensive,
                           cs_datatype_t           datatype,
                           cs_eval_at_location_t  *data_func,
                           void                   *data_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a function whose data values will be computed using the
 *        provided "degree of freedom" type evaluation function.
 *
 * The provided function and optional associated input is of the same
 * form as and may be shared with some boundary condition or property
 * definitions.
 *
 * \param[in]  name          name of associated function
 * \param[in]  location_id   id of associated mesh location
 * \param[in]  dim           dimension associated with element data
 * \param[in]  is_intensive  is the function intensive?
 * \param[in]  data_func     function used to define data values
 * \param[in]  data_input    pointer to optional (untyped) value or structure
 *                           to be used by data_func

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_by_analytic_func(const char          *name,
                                    int                  location_id,
                                    int                  dim,
                                    bool                 is_intensive,
                                    cs_analytic_func_t  *data_func,
                                    void                *data_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a function whose data values will be computed using the
 *        provided "degree of freedom" type evaluation function.
 *
 * The provided function and optional associated input is of the same
 * form as and may be shared with some boundary condition or property
 * definitions.
 *
 * \param[in]  name          name of associated function
 * \param[in]  location_id   id of associated mesh location
 * \param[in]  dim           dimension associated with element data
 * \param[in]  is_intensive  is the function intensive?
 * \param[in]  data_func     function used to define data values
 * \param[in]  data_input    pointer to optional (untyped) value or structure
 *                           to be used by data_func

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_by_dof_func(const char      *name,
                               int              location_id,
                               int              dim,
                               bool             is_intensive,
                               cs_dof_func_t   *data_func,
                               void            *data_input);

/*----------------------------------------------------------------------------
 * Destroy all function management metadata.
 *----------------------------------------------------------------------------*/

void
cs_function_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined functions.
 *
 * \return  number of defined functions
 */
/*----------------------------------------------------------------------------*/

int
cs_function_n_functions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a function object based on its id.
 *
 * This function requires that a function of the given id is defined.
 *
 * \param[in]  id   function id
 *
 * \return  pointer to the function structure
 */
/*----------------------------------------------------------------------------*/

cs_function_t  *
cs_function_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a function object based on its name.
 *
 * This function requires that a function of the given name is defined.
 *
 * \param[in]  name  function name
 *
 * \return  pointer to the function structure
 */
/*----------------------------------------------------------------------------*/

cs_function_t  *
cs_function_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a function object based on its name if present.
 *
 * If no function of the given name is defined, NULL is returned.
 *
 * \param[in]  name  function name
 *
 * \return  pointer to the function structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_function_t  *
cs_function_by_name_try(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assig a label to a function object.
 *
 * \param[in, out]  f      pointer to associated function handle
 * \param[in]       label  associated label
 */
/*----------------------------------------------------------------------------*/

void
cs_function_set_label(cs_function_t   *f,
                      const char      *label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log function definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_log_defs(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all given function object settings
 *        to log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_log_all_settings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate function values.
 *
 * If the matching values are multidimensional, they must be interleaved.
 * The output values are assumed to use a dense storage (i.e. of size
 * \c n_elts * <value dimension> for the associated data type, in the same
 * order as \c elt_ids if present.)
 *
 * \param[in]       f            pointer to associated function handle
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_evaluate(const cs_function_t   *f,
                     const cs_time_step_t  *ts,
                     int                    location_id,
                     cs_lnum_t              n_elts,
                     const cs_lnum_t       *elt_ids,
                     void                  *vals);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FUNCTION_H__ */
