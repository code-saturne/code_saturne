#ifndef __CS_PROPERTY_H__
#define __CS_PROPERTY_H__

/*============================================================================
 * Manage the definition/setting of properties
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_cdo.h"
#include "cs_param.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of property considered */
typedef enum {

  CS_PROPERTY_ISO,     // isotropic
  CS_PROPERTY_ORTHO,   // orthotropic
  CS_PROPERTY_ANISO,   // anisotropic
  CS_PROPERTY_N_TYPES

} cs_property_type_t;

typedef struct _cs_property_t cs_property_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect,
                                const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize cs_timer_stats_t structure for monitoring purpose
 *
 * \param[in]  level      level of details requested
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_timer_stats(int   level);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name          name of the property
 * \param[in]  key_type      keyname of the type of property
 * \param[in]  n_subdomains  piecewise definition on n_subdomains
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_create(const char    *name,
                   const char    *key_type,
                   int            n_subdomains);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_property_t structure
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_free(cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given property has the name ref_name
 *
 * \param[in]  pty         pointer to a cs_property_t structure to test
 * \param[in]  ref_name    name of the property to find
 *
 * \return true if the name of the property is ref_name otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_property_check_name(const cs_property_t   *pty,
                       const char            *ref_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the property is uniform, otherwise false
 *
 * \param[in]    pty    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_property_is_uniform(const cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of a property
 *
 * \param[in]    pty    pointer to a property
 *
 * \return  the name of the related property
 */
/*----------------------------------------------------------------------------*/

const char *
cs_property_get_name(const cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of a property
 *
 * \param[in]    pty    pointer to a property
 *
 * \return  the type of the related property
 */
/*----------------------------------------------------------------------------*/

cs_property_type_t
cs_property_get_type(const cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a summary of a cs_property_t structure
 *
 * \param[in]  pty      pointer to a cs_property_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_property_summary(const cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure by value for entities attached to
 *         the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_value(cs_property_t    *pty,
                         const char       *ml_name,
                         const char       *key_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an isotropic cs_property_t structure by value for entities
 *         attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_iso_def_by_value(cs_property_t    *pty,
                             const char       *ml_name,
                             double            val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define orthotropic cs_property_t structure by value for entities
 *         attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       values to set (vector of size 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_ortho_def_by_value(cs_property_t    *pty,
                               const char       *ml_name,
                               const double      val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       tens      values to set (3x3 tensor)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_aniso_def_by_value(cs_property_t    *pty,
                               const char       *ml_name,
                               const double      tens[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function in
 *         a subdomain attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       func      pointer to a cs_analytic_func_t function
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_analytic(cs_property_t        *pty,
                            const char           *ml_name,
                            cs_analytic_func_t   *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to law depending on one
 *         scalar variable in a subdomain attached to the mesh location named
 *         ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       context   pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a law function defined by subdomain
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_law(cs_property_t             *pty,
                       const char                *ml_name,
                       const void                *context,
                       cs_onevar_law_func_t      *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a law depending on
 *         two scalars variables in a subdomain attached to the mesh location
 *         named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       context     pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_twovar_law(cs_property_t          *pty,
                              const char             *ml_name,
                              const void             *context,
                              cs_twovar_law_func_t   *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure by a law (using as arguments a
 *         scalar and a vector) for entities attached to the mesh location
 *         named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       context   pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a law function defined by subdomain
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_scavec_law(cs_property_t             *pty,
                              const char                *ml_name,
                              const void                *context,
                              cs_scavec_law_func_t      *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       desc      information about this array
 * \param[in]       array     pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_array(cs_property_t    *pty,
                         cs_desc_t         desc,
                         cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the "array" member of a cs_property_t structure
 *
 * \param[in, out]  pty          pointer to a cs_property_t structure
 * \param[in]       desc         information about this array
 * \param[in]       array        pointer to an array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_array(cs_property_t    *pty,
                      cs_desc_t         desc,
                      cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the second "array" member of a cs_property_t structure
 *
 * \param[in, out]  pty        pointer to a cs_property_t structure
 * \param[in]       desc       information about this array
 * \param[in]       array      pointer to an array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_second_array(cs_property_t    *pty,
                             cs_desc_t         desc,
                             cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *
 * \param[in]      c_id           id of the current cell
 * \param[in]      pty            pointer to a cs_property_t structure
 * \param[in]      do_inversion   true or false
 * \param[in, out] tensor         3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_cell_tensor(cs_lnum_t             c_id,
                            const cs_property_t  *pty,
                            bool                  do_inversion,
                            cs_real_3_t          *tensor);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *
 * \param[in]   c_id           id of the current cell
 * \param[in]   pty            pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_get_cell_value(cs_lnum_t              c_id,
                           const cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Fourier number in each cell
 *
 * \param[in]      pty        pointer to the diffusive property struct.
 * \param[in]      dt         value of the current time step
 * \param[in, out] fourier    pointer to an array storing Fourier numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_fourier(const cs_property_t     *pty,
                        double                   dt,
                        cs_real_t                fourier[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROPERTY_H__ */
