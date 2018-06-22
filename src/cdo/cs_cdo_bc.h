#ifndef __CS_CDO_BC_H__
#define __CS_CDO_BC_H__

/*============================================================================
 * Manage boundary conditions. Produce ready-to-use BC structure from the
 * user definition.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @defgroup cdo_bc_flags Flags specifying the type of boundary conditions
 *  associated to an element
 *
 * @{
 */

/*!  1: (non-homogeneous) Dirichlet boundary conditions */
#define CS_CDO_BC_DIRICHLET      (1 << 0)
/*!  2: Homogeneous Dirichlet boundary conditions */
#define CS_CDO_BC_HMG_DIRICHLET  (1 << 1)
/*!  4: (non-homogeneous) Neumann boundary conditions */
#define CS_CDO_BC_NEUMANN        (1 << 2)
/*!  8: Homogeneous Neumann boundary conditions */
#define CS_CDO_BC_HMG_NEUMANN    (1 << 3)
/*! 16: Robin boundary conditions */
#define CS_CDO_BC_ROBIN          (1 << 4)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* List of entities attached to a type of BC (for instance Dirichlet or Neumann)
   Two categories are considered
   First, the entities attached to a non-homogeneous BC, then those attached
   to a homogeneous BC.
*/

typedef struct {

  cs_lnum_t   n_elts;
  cs_lnum_t   n_nhmg_elts; /* number of non-homogeneous elements */

  cs_lnum_t  *elt_ids;     /* size = n_elts (first elements are those associated
                              to non-homogeneous BC then the homogeneous one */
  short int  *def_ids;     /* id related to the associated BC definition
                              Only for non homogeneous BCs (i.e. size is equal
                              to the number elements associated to a
                              non-homogeneous BC */

} cs_cdo_bc_list_t;

/* Translation of the user-defined BCs setup into a computable-oriented
   structure */

typedef struct {

  cs_lnum_t          n_elts;    /* Number of elements */
  cs_flag_t         *flag;      /* Type of boundary conditions associated to
                                   an element. For a face, one (and only one)
                                   type is set. For an edge and a vertex,
                                   several BCs can be set.
                                   size = n_elts */

  /* List of faces by type of boundary conditions */
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
  default:
    ret_flag = 0;
    break;
  }
  return ret_flag;
}

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
 * \brief  Define the structure which translates the BC definition from the
 *         user viewpoint into a ready-to-use structure for setting the arrays
 *         keeping the values of the boundary condition to set.
 *
 * \param[in] default_bc   type of boundary condition to set by default
 * \param[in] n_desc       number of boundary definitions
 * \param[in] desc         list of boundary condition definition
 * \param[in] n_b_faces    number of border faces
 *
 * \return a pointer to a new allocated cs_cdo_bc_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_t *
cs_cdo_bc_define(cs_param_bc_type_t    default_bc,
                 int                   n_desc,
                 cs_xdef_t           **desc,
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

END_C_DECLS

#endif /* __CS_CDO_BC_H__ */
