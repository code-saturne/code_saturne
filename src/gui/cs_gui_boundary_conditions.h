#ifndef __CS_GUI_BOUNDARY_CONDITION_H__
#define __CS_GUI_BOUNDARY_CONDITION_H__

/*============================================================================
 * Management of the GUI parameters file: boundary conditions
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Arguments passed by context pointer to cs_meg_* functions */

typedef struct {

  const  cs_zone_t    *zone;        /*<! Pointer to zone */

  const  char         *name;        /*<! Pointer to field or array name */
  const  char         *condition;   /*<! Pointer to condition name type */

  int                  dim;         /*<! Values dimension */

} cs_gui_boundary_meg_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * \param[in, out] itypfb   <-- type of boundary for each face
 *----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_processing(int  *itypfb);

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_verify(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define boundary conditions based on setup file.
 *
 * \param[in, out]  bdy   boundaries structure to update
 *                        (if NULL, default to cs_glob_domain->boundaries)
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_define(cs_boundary_t  *bdy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free GUI boundary condition structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_free_memory(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add new MEG-based cs_dof_func_t context info.
 *
 * \param[in]  zone       pointer to associated zone
 * \param[in]  name       name of associated field or array
 * \param[in]  condition  associated condition type
 * \param[in]  dim        associated dimension
 *
 * \return: pointer to cs_dof_func_t context info
 */
/*----------------------------------------------------------------------------*/

cs_gui_boundary_meg_context_t *
cs_gui_boundary_add_meg_context(const  cs_zone_t   *zone,
                                const  char        *name,
                                const  char        *condition,
                                int                 dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a boundary profiles
 *        using a MEG generated function.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_dof_func_meg(cs_lnum_t         n_elts,
                                        const cs_lnum_t  *elt_ids,
                                        bool              dense_output,
                                        void             *input,
                                        cs_real_t        *retval);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_BOUNDARY_CONDITION_H__ */
