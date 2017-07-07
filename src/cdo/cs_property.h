#ifndef __CS_PROPERTY_H__
#define __CS_PROPERTY_H__

/*============================================================================
 * Manage the definition/setting of properties
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_xdef.h"
#include "cs_xdef_eval.h"

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

/* Set of parameters attached to a property */
typedef struct {

  char  *restrict  name;
  int              id;
  cs_flag_t        state_flag;

  /* The number of values to set depends on the type of property
      isotropic   = 1, orthotropic = 3, anisotropic = 9  */
  cs_property_type_t   type;

  /* Property is up to now only defined on the whole domain (volume) */
  int                  n_definitions;  /* Current number of definions used */
  cs_xdef_t          **defs;           /* List of definitions */
  short int           *def_ids;        /* Store the definition id for each cell,
                                          NULL is only one definition is set */

  /* Function pointers to handle generic tasks related to a property. There
     is one function related to each definition. Some functions may not be
     allocated according to the kind of property */

  /* Retrieve the evaluation of the property at the cell center for each
     definition */
  cs_xdef_eval_t     **get_eval_at_cell;

  /* Same thing as the previous one but now with the usage of cellwise algo.
     relying on a cs_cell_mesh_t structure */
  cs_xdef_eval_cw_t  **get_eval_at_cell_cw;

} cs_property_t;

/* List of available keys for setting an advection field */
typedef enum {

  CS_PTYKEY_POST,
  CS_PTYKEY_N_KEYS

} cs_property_key_t;

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
 * \brief  Retrieve the number of properties
 *
 * \return the number of properties
 */
/*----------------------------------------------------------------------------*/

int
cs_property_get_n_properties(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name          name of the property
 * \param[in]  type          type of property
 * \param[in]  n_subdomains  piecewise definition on n_subdomains
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add(const char            *name,
                cs_property_type_t     type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the related property definition from its name
 *
 * \param[in]  name    name of the property to find
 *
 * \return NULL if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_name(const char   *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the related property definition from its id
 *
 * \param[in]  id      id of the property to find
 *
 * \return NULL if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_id(int         id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all cs_property_t structures and the array storing all the
 *         structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage of the definition of a property based on several
 *         definitions (i.e. definition by subdomains)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_finalize_setup(void);

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
 * \brief  Define an isotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty          pointer to a cs_property_t structure
 * \param[in]       zone_name    name of an already defined zone
 * \param[in]       val          value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_iso_by_value(cs_property_t    *pty,
                             const char       *zone_name,
                             double            val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an orthotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty         pointer to a cs_property_t structure
 * \param[in]       zone_name   name of an already defined zone
 * \param[in]       val         values to set (vector of size 3)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_ortho_by_value(cs_property_t    *pty,
                               const char       *zone_name,
                               double            val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty         pointer to a cs_property_t structure
 * \param[in]       zone_name   name of an already defined zone
 * \param[in]       tens        values to set (3x3 tensor)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_aniso_by_value(cs_property_t    *pty,
                               const char       *zone_name,
                               cs_real_t         tens[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function in
 *         a subdomain attached to the mesh location named ml_name
 *
 * \param[in, out]  pty         pointer to a cs_property_t structure
 * \param[in]       zone_name   name of an already defined zone
 * \param[in]       func        pointer to a cs_analytic_func_t function
 * \param[in]       input       NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_analytic(cs_property_t        *pty,
                            const char           *zone_name,
                            cs_analytic_func_t   *func,
                            void                 *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to law depending on one
 *         scalar variable in a subdomain attached to the mesh location named
 *         ml_name
 *
 * \param[in, out] pty                  pointer to a cs_property_t structure
 * \param[in]      zone_name            name of an already defined zone
 * \param[in]      context              pointer to a structure (may be NULL)
 * \param[in]      get_eval_at_cell     pointer to a function (may be NULL)
 * \param[in]      get_eval_at_cell_cw  pointer to a function
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_func(cs_property_t         *pty,
                        const char            *zone_name,
                        void                  *context,
                        cs_xdef_eval_t        *get_eval_at_cell,
                        cs_xdef_eval_cw_t     *get_eval_at_cell_cw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_array(cs_property_t    *pty,
                         cs_flag_t         loc,
                         cs_real_t        *array,
                         cs_lnum_t        *index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_field(cs_property_t    *pty,
                         cs_field_t       *field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the value of the property at each cell. Store the
 *         evaluation in the given array.
 *
 * \param[in]       pty      pointer to a cs_property_t structure
 * \param[in, out]  array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_cells(const cs_property_t    *pty,
                          cs_real_t              *array);

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
cs_property_get_cell_tensor(cs_lnum_t               c_id,
                            const cs_property_t    *pty,
                            bool                    do_inversion,
                            cs_real_3_t            *tensor);

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
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      pty            pointer to a cs_property_t structure
 * \param[in]      do_inversion   true or false
 * \param[in, out] tensor         3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_tensor_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_property_t    *pty,
                           bool                    do_inversion,
                           cs_real_3_t            *tensor);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]   cm        pointer to a cs_cell_mesh_t structure
 * \param[in]   pty       pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_value_in_cell(const cs_cell_mesh_t   *cm,
                          const cs_property_t    *pty);

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
/*!
 * \brief  Print a summary of the settings for all defined cs_property_t
 *         structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROPERTY_H__ */
