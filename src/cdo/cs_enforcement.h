#ifndef __CS_ENFORCEMENT_H__
#define __CS_ENFORCEMENT_H__

/*============================================================================
 * Manage the list of solid cells and associated helper functions
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Structure and type definitions
 *============================================================================*/

/*! \enum cs_enforcement_selection_t
 *  \brief Type of entities used to define the selection where the enforcement
 *         takes place
 *
 * Two mechanisms are possible.
 *
 * 1) Cell selection: defined a selection of cells and then automatically built
 *    the related selection of degrees of freedom and assigned a value to each
 *    selected degrees of freedom
 *
 * 2) DoF selection (faces or vertices up to now): defined a selection of
 *    degrees of freedom (DoFs) and assign a values to a selection of degrees
 *    of freedom inside
 */

typedef enum {

  CS_ENFORCEMENT_SELECTION_CELLS,     /*!< List of cell ids */
  CS_ENFORCEMENT_SELECTION_FACES,     /*!< List of face ids */
  CS_ENFORCEMENT_SELECTION_EDGES,     /*!< List of edge ids */
  CS_ENFORCEMENT_SELECTION_VERTICES   /*!< List of vertex ids */

} cs_enforcement_selection_t;


/*! \enum cs_enforcement_type_t
 *  \brief Describe the way the values to enforce are defined
 */

typedef enum {

  CS_ENFORCEMENT_BY_CONSTANT,       /*!< The same constant value for each DoF */
  CS_ENFORCEMENT_BY_DOF_VALUES      /*!< A prescribed value for each DoF */

} cs_enforcement_type_t;


/*! \struct cs_enforcement_param_t
 *  \brief Set of data defining an enforcement
 */

typedef struct {

  cs_enforcement_selection_t    selection_type;
  cs_enforcement_type_t         type;

  cs_lnum_t     n_elts;         /*!< local number of selected entities */
  cs_lnum_t    *elt_ids;        /*!< local list of selected entities */
  int           stride;         /*!< number of values by entity */
  cs_real_t    *values;         /*!< associated values */

} cs_enforcement_param_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and define a cs_enforcement_param_t structure
 *
 * \param[in] sel_type   type of elements which have been selected
 * \param[in] type       way to set values for the selected elements
 * \param[in] stride     number of values to enforce by element
 * \param[in] n_elts     number of selected elements locally
 * \param[in] elt_ids    list of element ids
 * \param[in] values     array of values to enforce
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_enforcement_param_create(cs_enforcement_selection_t    sel_type,
                            cs_enforcement_type_t         type,
                            int                           stride,
                            cs_lnum_t                     n_elts,
                            const cs_lnum_t              *elt_ids,
                            const cs_real_t              *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset an existing cs_enforcement_param_t structure
 *
 * \param[in, out] efp        pointer to a cs_enforcement_param_t structure
 * \param[in]      sel_type   type of elements which have been selected
 * \param[in]      type       way to set values for the selected elements
 * \param[in]      stride     number of values to enforce by element
 * \param[in]      n_elts     number of selected elements locally
 * \param[in]      elt_ids    list of element ids
 * \param[in]      values     array of values to enforce
 */
/*----------------------------------------------------------------------------*/

void
cs_enforcement_param_reset(cs_enforcement_param_t       *efp,
                           cs_enforcement_selection_t    sel_type,
                           cs_enforcement_type_t         type,
                           int                           stride,
                           cs_lnum_t                     n_elts,
                           const cs_lnum_t              *elt_ids,
                           const cs_real_t              *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy a cs_enforcement_param_t structure
 *
 * \param[in] ref    reference structure to copy
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_enforcement_param_copy(const cs_enforcement_param_t   *ref);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_enforcement_param_t structure
 *
 * \param[in, out] p_efp    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_enforcement_param_free(cs_enforcement_param_t   **p_efp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a cs_enforcement_param_t structure
 *
 * \param[in] eqname  name of the related equation
 * \param[in] efp     pointer to a  cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_enforcement_param_log(const char                     *eqname,
                         const cs_enforcement_param_t   *efp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_enforcement_t structure for vertex-based scheme
 *
 * \param[in] connect    pointer to a cs_cdo_connect_t
 * \param[in] n_params   number of enforcement parameters
 * \param[in] efp_array  array of parameter structures defining the enforcement
 *
 * \return an array with the values to enforce
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_enforcement_define_at_vertices(const cs_cdo_connect_t     *connect,
                                  int                         n_params,
                                  cs_enforcement_param_t    **efp_array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_enforcement_t structure for face-based scheme
 *
 * \param[in] connect    pointer to a cs_cdo_connect_t
 * \param[in] n_params   number of enforcement parameters
 * \param[in] efp_array  array of parameter structures defining the enforcement
 *
 * \return an array with the values to enforce
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_enforcement_define_at_faces(const cs_cdo_connect_t     *connect,
                               int                         n_params,
                               cs_enforcement_param_t    **efp_array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_enforcement_t structure for edge-based scheme
 *
 * \param[in] connect    pointer to a cs_cdo_connect_t
 * \param[in] n_params   number of enforcement parameters
 * \param[in] efp_array  array of parameter structures defining the enforcement
 *
 * \return an array with the values to enforce
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_enforcement_define_at_edges(const cs_cdo_connect_t     *connect,
                               int                         n_params,
                               cs_enforcement_param_t    **efp_array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the cell-wise value to enforce
 *
 * \param[in]      forced_values     values to enforce or FLT_MAX
 * \param[in, out] csys              pointer to a cs_cell_sys_t structure
 * \param[in, out] cw_forced_values  local values to enforce
 *
 * \return true if at least one DoF has to be enforced
 */
/*----------------------------------------------------------------------------*/

bool
cs_enforcement_dofs_cw(const cs_real_t      *forced_values,
                       cs_cell_sys_t        *csys,
                       cs_real_t            *cw_forced_values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ENFORCEMENT_H__ */
