#ifndef __CS_CDO_BC_H__
#define __CS_CDO_BC_H__

/*============================================================================
 * Manage the low-level structure dedicated to boundary conditions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_time_step.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_CDO_BC_DEFAULT_DEF    -1

/*!
 * @defgroup cdo_bc_flags Flags specifying the type of boundary conditions
 *  associated to an element
 *
 * Homogeneous conditions are defined separately since a flag
 * CS_CDO_BC_HOMOGENEOUS would not enable to identify if it is associated to a
 * Dirichlet or a Neumann boundary condition
 *
 * @{
 */

/*!  1: Neumann boundary conditions */
#define CS_CDO_BC_NEUMANN               (1 << 0)
/*!  2: Homogeneous Neumann boundary conditions */
#define CS_CDO_BC_HMG_NEUMANN           (1 << 1)
/*!  4: Dirichlet boundary conditions */
#define CS_CDO_BC_DIRICHLET             (1 << 2)
/*!  8: Homogeneous Dirichlet boundary conditions */
#define CS_CDO_BC_HMG_DIRICHLET         (1 << 3)
/*! 16: Robin boundary conditions */
#define CS_CDO_BC_ROBIN                 (1 << 4)
/*! 32: Apply a sliding condition (for vector-valued equations) */
#define CS_CDO_BC_SLIDING               (1 << 5)
/*! 64: Apply a Dirichlet on the tangential part of a vector-valued quantity */
#define CS_CDO_BC_TANGENTIAL_DIRICHLET  (1 << 6)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure specific to store data related to the definition of boundary
 * conditions on boundary faces.
 *
 * For of scalar-valued equations, only some the classical (Dirichlet, Neumann
 * and Robin types are available. Other types of boundary conditions are
 * possible for vector-valued equations
 */

typedef struct {

  bool         is_steady;    /* Do we need to update BC faces during the
                                computation */
  cs_lnum_t    n_b_faces;    /* Number of boundary faces */

  /* Type of boundary conditions associated to a face. Size: n_b_faces */
  cs_flag_t   *flag;

  /* Id of the boundary condition definition or CS_BC_DEFAULT (=-1) if this face
     is related to the default boundary condition. Size = n_b_faces */
  short int   *def_ids;

  /* List of face ids by type of boundary conditions. Homogeneous types don't
   * need to rely on a definition since it can be the default bc. Moreover, some
   * optimizations can be performed that's why they are stored separately
   */

  /* Dirichlet */
  cs_lnum_t    n_hmg_dir_faces;
  cs_lnum_t   *hmg_dir_ids;
  cs_lnum_t    n_nhmg_dir_faces;
  cs_lnum_t   *nhmg_dir_ids;

  /* Neumann */
  cs_lnum_t    n_hmg_neu_faces;
  cs_lnum_t   *hmg_neu_ids;
  cs_lnum_t    n_nhmg_neu_faces;
  cs_lnum_t   *nhmg_neu_ids;

  /* Robin */
  cs_lnum_t    n_robin_faces;
  cs_lnum_t   *robin_ids;

  /* Sliding wall */
  cs_lnum_t    n_sliding_faces;
  cs_lnum_t   *sliding_ids;

} cs_cdo_bc_face_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Convert a flag into a description
 *
 * \param[in]       bc_flag  flag of boundary condition
 * \param[in, out]  desc     string storing the description of the BC
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_cdo_bc_get_desc(cs_flag_t   bc_flag,
                   char       *desc)
{
  if (desc == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Empty desciption buffer.", __func__);

  switch (bc_flag) {

  case CS_CDO_BC_HMG_DIRICHLET:
    sprintf(desc, "%s", "Homogenous Dirichlet");
    break;
  case CS_CDO_BC_DIRICHLET:
    sprintf(desc, "%s", "Dirichlet");
    break;
  case CS_CDO_BC_HMG_NEUMANN:
    sprintf(desc, "%s", "Homogeneous Neumann");
    break;
  case CS_CDO_BC_NEUMANN:
    sprintf(desc, "%s", "Neumann");
    break;
  case CS_CDO_BC_ROBIN:
    sprintf(desc, "%s", "Robin");
    break;
  case CS_CDO_BC_SLIDING:
    sprintf(desc, "%s", "Sliding");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid case. Please contact the support.\n", __func__);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Convert a cs_param_bc_type_t into a flag (enable multiple type for
 *          a same entity as required for vertices and edges)
 *
 * \param[in] bc_type   predefined type of boundary condition
 *
 * \return  a flag corresponding to the given type of boundary condition
 */
/*----------------------------------------------------------------------------*/

static inline cs_flag_t
cs_cdo_bc_get_flag(cs_param_bc_type_t   bc_type)
{
  cs_flag_t  ret_flag;

  switch (bc_type) {
  case CS_PARAM_BC_HMG_DIRICHLET:
    ret_flag = CS_CDO_BC_HMG_DIRICHLET;
    break;
  case CS_PARAM_BC_DIRICHLET:
    ret_flag = CS_CDO_BC_DIRICHLET;
    break;
  case CS_PARAM_BC_HMG_NEUMANN:
    ret_flag = CS_CDO_BC_HMG_NEUMANN;
    break;
  case CS_PARAM_BC_NEUMANN:
    ret_flag = CS_CDO_BC_NEUMANN;
    break;
  case CS_PARAM_BC_ROBIN:
    ret_flag = CS_CDO_BC_ROBIN;
    break;
  case CS_PARAM_BC_SLIDING:
    ret_flag = CS_CDO_BC_SLIDING;
    break;
  case CS_PARAM_BC_CIRCULATION:
    ret_flag = CS_CDO_BC_TANGENTIAL_DIRICHLET;
    break;

  default:
    ret_flag = 0;               /* Not handle automatically */
    break;
  }
  return ret_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check if a flag is associated to a Dirichlet BC (homogeneous or
 *          not)
 *
 * \param[in] flag     flag to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_cdo_bc_is_dirichlet(cs_flag_t    flag)
{
  if (flag & CS_CDO_BC_DIRICHLET)
    return true;
  else if (flag & CS_CDO_BC_HMG_DIRICHLET)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check if a flag is associated to a Neumann BC (homogeneous or not)
 *
 * \param[in] flag     flag to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_cdo_bc_is_neumann(cs_flag_t    flag)
{
  if (flag & CS_CDO_BC_NEUMANN)
    return true;
  else if (flag & CS_CDO_BC_HMG_NEUMANN)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check if a flag is associated to a sliding boundary
 *
 * \param[in] flag     flag to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_cdo_bc_is_sliding(cs_flag_t    flag)
{
  if (flag & CS_CDO_BC_SLIDING)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the structure which translates the BC definitions from the
 *         user viewpoint into a ready-to-use structure for setting the arrays
 *         keeping the values of the boundary condition to set.
 *
 * \param[in] default_bc   type of boundary condition to set by default
 * \param[in] is_steady    modification or not of the BC selection in time
 * \param[in] dim          dimension of the related equation
 * \param[in] n_defs       number of boundary definitions
 * \param[in] defs         list of boundary condition definition
 * \param[in] n_b_faces    number of border faces
 *
 * \return a pointer to a new allocated cs_cdo_bc_face_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_face_t *
cs_cdo_bc_face_define(cs_param_bc_type_t    default_bc,
                      bool                  is_steady,
                      int                   dim,
                      int                   n_defs,
                      cs_xdef_t           **defs,
                      cs_lnum_t             n_b_faces);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_bc_face_t structure
 *
 * \param[in, out]  face_bc   pointer to a cs_cdo_bc_face_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_face_t *
cs_cdo_bc_free(cs_cdo_bc_face_t   *face_bc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_BC_H__ */
