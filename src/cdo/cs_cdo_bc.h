#ifndef __CS_CDO_BC_H__
#define __CS_CDO_BC_H__

/*============================================================================
 * Manage boundary conditions. Produce ready-to-use BC structure from the
 * user definition.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_time_step.h"

#include "cs_param.h"
#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* List of entities attached to a type of BC (for instance Dirichlet or Neumann)
   Two categories are considered
   First, the entities attached to a non-homogeneous BC, then those attached
   to a homogeneous BC.
   Only entities with a non-homogeneous used stride and values items
*/

typedef struct {

  cs_lnum_t   n_elts;
  cs_lnum_t   n_nhmg_elts; /* number of non-homogeneous elements */

  cs_lnum_t  *elt_ids;    /* size = n_elts */
  short int  *def_ids;    /* size = n_nhmg_elts */

} cs_cdo_bc_list_t;

/* Translation of the user-defined BCs setup into a computable-oriented
   structure */

typedef struct {

  cs_lnum_t   n_b_faces;

  // Selection of faces
  cs_cdo_bc_list_t  *dir;
  cs_cdo_bc_list_t  *neu;
  cs_cdo_bc_list_t  *rob;

} cs_cdo_bc_t;


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_cdo_bc_list_t structure
 *
 * \param[in]  n_elts      number of entries of the list
 * \param[in]  n_nhmg_elts  number of elements attached to a homogeneous BC
 *
 * \return  a new allocated pointer to a cs_cdo_bc_list_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_list_create(cs_lnum_t   n_elts,
                      cs_lnum_t   n_nhmg_elts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_cdo_bc_list_t structure
 *
 * \param[in]  bcl     pointer to the cs_cdo_bc_list_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_list_free(cs_cdo_bc_list_t   *bcl);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Prepare the treatment of the boundary conditions.
 *         Compile the information detailed in a cs_param_bc_t structure
 *         into the structure cs_cdo_bc_t (based on border faces).
 *         This is a primilary step to be ready to set the values of the BC
 *
 * \param[in] param_bc    pointer to the parameters related to BCs of an eq.
 * \param[in] n_b_faces   number of border faces
 *
 * \return a pointer to a new allocated cs_cdo_bc_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_t *
cs_cdo_bc_init(const cs_param_bc_t  *param_bc,
               cs_lnum_t             n_b_faces);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_bc_t structure
 *
 * \param[in, out]  face_bc   pointer to a cs_cdo_bc_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_t *
cs_cdo_bc_free(cs_cdo_bc_t   *face_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build cs_cdo_bc_list_t structures for Dirichlet BC on primal
 *         vertices.
 *         When there is a choice between homogeneous or non homogeneous BCs,
 *         we always set the homogeneous condition for Dirichlet BCs.
 *
 * \param[in] m         pointer to a cs_mesh_t structure
 * \param[in] face_bc   pointer to a cs_cdo_bc_t structure
 *
 * \return a pointer to a new allocated cs_cdo_bc_list_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_vtx_dir_create(const cs_mesh_t    *m,
                         const cs_cdo_bc_t  *face_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the Dirichlet values to enforce on the corresponding entities
 *
 * \param[in]      dof_flag  information about the corresponding DoF to treat
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      geom      structure storing geometric information
 * \param[in]      bc        pointer to a cs_param_bc_t structure
 * \param[in]      ent_dir   pointer to a cs_cdo_bc_list_t
 * \param[in, out] dir_val   array used to store Dirichlet values
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_bc_dirichlet_set(cs_flag_t                dof_flag,
                        const cs_time_step_t    *time_step,
                        const void              *geom,
                        const cs_param_bc_t     *bc,
                        const cs_cdo_bc_list_t  *ent_dir,
                        double                  *dir_val);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_BC_H__ */
