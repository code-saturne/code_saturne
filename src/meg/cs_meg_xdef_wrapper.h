#ifndef __CS_MEG_XDEF_WRAPPER_H__
#define __CS_MEG_XDEF_WRAPPER_H__

/*============================================================================
 * MEG (Mathematical Expression Generator) functions xdef wrapper
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Enum for MEG function type */

typedef enum {
  CS_MEG_BOUNDARY_FUNC,
  CS_MEG_VOLUME_FUNC,
  CS_MEG_INITIALIZATION_FUNC,
  CS_MEG_SOURCE_TERM_FUNC,
  CS_MEG_IBM_FUNC,
  CS_MEG_FSI_STRUCT_FUNC,
  CS_MEG_POST_ACTIVATE_FUNC,
  CS_MEG_POST_PROFILES_FUNC,
  CS_MEG_CALCULATOR_FUNC,

  CS_MEG_N_FUNC_TYPES
} cs_meg_function_type_t;

typedef struct {
  /* type of meg function */
  cs_meg_function_type_t type;

  /* cs_zone_t id */
  int z_id;

  /* values stride size */
  int stride;

  /* Input name used for function */
  char name[512];

  /* boundary or source term data */
  char additional_data[512];

} cs_meg_xdef_input_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a MEG function xdef wrapper input data. Allocated memory is
 *  deleted by cs_meg_xdef_wrapper_finalize
 *
 * \param[in] type              type of meg function linked to this input
 * \param[in] z_id              id of zone on which this function is defined
 * \param[in] stride            stride of data
 * \param[in] name              name related to function
 * \param[in] additional_data   additional data (char *) provided to function,
 *                              such as condition or source type
 *
 * \returns pointer to newly allocated input data structure
 */
/*----------------------------------------------------------------------------*/

cs_meg_xdef_input_t *
cs_meg_xdef_wrapper_add_input(const cs_meg_function_type_t type,
                              const int                    z_id,
                              const int                    stride,
                              const char                  *name,
                              const char                  *additional_data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper function allowing to call MEG functions by xdef structres.
 * This is done by using the cs_xdef_analytic_function type.
 *
 * \param[in] time          when ?
 * \param[in] n_elts        number of elements to consider
 * \param[in] elt_ids       list of elements ids (in coords and retval)
 * \param[in] coords        where ?
 * \param[in] dense_output  perform an indirection in retval or not
 * \param[in] input         pointer to cs_meg_xdef_input_t
 * \param[in] retval        resultint value(s). Must be allocated
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_xdef_wrapper(cs_real_t         time,
                    cs_lnum_t         n_elts,
                    const cs_lnum_t  *elt_ids,
                    const cs_real_t  *coords,
                    bool              dense_output,
                    void             *input,
                    cs_real_t        *retval);
/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif
