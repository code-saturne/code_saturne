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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_boundary_zone.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_bc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_cdo_bc_face_t structure
 *
 * \param[in] is_steady   true or false
 * \param[in] n_b_faces   number of boundary faces
 *
 * \return  a new allocated pointer to a cs_cdo_bc_face_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cdo_bc_face_t *
_cdo_bc_face_create(bool       is_steady,
                    cs_lnum_t  n_b_faces)
{
  cs_cdo_bc_face_t  *bc = NULL;

  BFT_MALLOC(bc, 1, cs_cdo_bc_face_t);

  bc->is_steady = is_steady;
  bc->n_b_faces = n_b_faces;

  /* Default initialization */
  bc->flag = NULL;
  BFT_MALLOC(bc->flag, n_b_faces, cs_flag_t);
  memset(bc->flag, 0, n_b_faces*sizeof(cs_flag_t));

  /* Default initialization */
  bc->def_ids = NULL;
  BFT_MALLOC(bc->def_ids, n_b_faces, short int);
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    bc->def_ids[i] = CS_CDO_BC_DEFAULT_DEF;

  /* Other lists of faces */
  bc->n_hmg_dir_faces = 0;
  bc->hmg_dir_ids = NULL;
  bc->n_nhmg_dir_faces = 0;
  bc->nhmg_dir_ids = NULL;

  bc->n_hmg_neu_faces = 0;
  bc->hmg_neu_ids = NULL;
  bc->n_nhmg_neu_faces = 0;
  bc->nhmg_neu_ids = NULL;

  bc->n_robin_faces = 0;
  bc->robin_ids = NULL;

  bc->n_sliding_faces = 0;
  bc->sliding_ids = NULL;

  bc->n_circulation_faces = 0;
  bc->circulation_ids = NULL;

  return bc;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
                      cs_lnum_t             n_b_faces)
{
  CS_UNUSED(dim); /* Only in debug */

  /* Set the default flag */
  cs_flag_t  default_flag = cs_cdo_bc_get_flag(default_bc);

  if (!(default_flag & CS_CDO_BC_HMG_DIRICHLET) &&
      !(default_flag & CS_CDO_BC_HMG_NEUMANN) &&
      !(default_flag & CS_CDO_BC_SLIDING))
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Incompatible type of boundary condition by default.\n"
                " Please modify your settings.\n"), __func__);

  cs_cdo_bc_face_t  *bc = _cdo_bc_face_create(is_steady, n_b_faces);

  if (n_b_faces == 0) /* In parallel run this situation may occur */
    return  bc;

  /* Loop on the definition of each boundary condition */
  for (int ii = 0; ii < n_defs; ii++) {

    const cs_xdef_t  *d = defs[ii];
    const cs_zone_t  *z = cs_boundary_zone_by_id(d->z_id);

    switch (d->meta) {

    case CS_CDO_BC_DIRICHLET:
      bc->n_nhmg_dir_faces += z->n_elts;
      break;
    case CS_CDO_BC_HMG_DIRICHLET:
      bc->n_hmg_dir_faces += z->n_elts;
      break;
    case CS_CDO_BC_NEUMANN:
      bc->n_nhmg_neu_faces += z->n_elts;
      break;
    case CS_CDO_BC_HMG_NEUMANN:
      bc->n_hmg_neu_faces += z->n_elts;
      break;
    case CS_CDO_BC_ROBIN:
      bc->n_robin_faces += z->n_elts;
      break;

      /* For vector-valued equations only */
    case CS_CDO_BC_SLIDING:
      assert(dim > 1);
      bc->n_sliding_faces += z->n_elts;
      break;

      /* For vector-valued equations only */
    case CS_CDO_BC_TANGENTIAL_DIRICHLET:
      assert(dim > 1);
      bc->n_circulation_faces += z->n_elts;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: This type of boundary condition is not handled.",
                __func__);
    } /* End of switch */

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t bf_id = z->elt_ids[i];

      assert(bc->flag[bf_id] < 1);                         /* Not already set */
      assert(bc->def_ids[bf_id] == CS_CDO_BC_DEFAULT_DEF); /* Not already set */
      bc->flag[bf_id] = d->meta;
      bc->def_ids[bf_id] = ii;  /* definition id */

    } /* Loop on faces associated to this definition */

  } /* Loop on definitions of boundary conditions */

  /* Set the default flag for all remaining unset faces */
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    if (bc->flag[i] == 0) { /* Not set yet --> apply the default settings */

      assert(bc->def_ids[i] == CS_CDO_BC_DEFAULT_DEF);
      bc->flag[i] = default_flag;
      if (default_flag & CS_CDO_BC_HMG_DIRICHLET)
        bc->n_hmg_dir_faces++;
      else if (default_flag & CS_CDO_BC_HMG_NEUMANN)
        bc->n_hmg_neu_faces++;
      else if (default_flag & CS_CDO_BC_SLIDING)
        bc->n_sliding_faces++;
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of default boundary condition", __func__);

    } /* Unset face */
  } /* Loop on border faces */

#if defined(DEBUG) && !defined(NDEBUG)
  /* Sanity check (there is no multiple definition or faces whitout settings */
  cs_lnum_t n_set_faces =
    bc->n_hmg_neu_faces + bc->n_nhmg_neu_faces +
    bc->n_hmg_dir_faces + bc->n_nhmg_dir_faces +
    bc->n_robin_faces + bc->n_sliding_faces + bc->n_circulation_faces;
  if (n_set_faces != bc->n_b_faces)
    bft_error(__FILE__, __LINE__, 0,
              " %s: There are %d faces without boundary conditions.\n"
              " Please check your settings.",
              __func__, bc->n_b_faces - n_set_faces);
#endif

  /* Allocate list of border faces by type of boundary conditions */
  BFT_MALLOC(bc->hmg_dir_ids, bc->n_hmg_dir_faces, cs_lnum_t);
  BFT_MALLOC(bc->nhmg_dir_ids, bc->n_nhmg_dir_faces, cs_lnum_t);
  BFT_MALLOC(bc->hmg_neu_ids, bc->n_hmg_neu_faces, cs_lnum_t);
  BFT_MALLOC(bc->nhmg_neu_ids, bc->n_nhmg_neu_faces, cs_lnum_t);
  BFT_MALLOC(bc->robin_ids, bc->n_robin_faces, cs_lnum_t);
  BFT_MALLOC(bc->sliding_ids, bc->n_sliding_faces, cs_lnum_t);
  BFT_MALLOC(bc->circulation_ids, bc->n_circulation_faces, cs_lnum_t);

  /* Fill the allocated lists */
  cs_lnum_t  shift[CS_PARAM_N_BC_TYPES];
  for (int ii = 0; ii < CS_PARAM_N_BC_TYPES; ii++)
    shift[ii] = 0;

  /* Loop on the border faces and append lists */
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {

    switch (bc->flag[i]) {

    case CS_CDO_BC_DIRICHLET:
      bc->nhmg_dir_ids[shift[CS_PARAM_BC_DIRICHLET]] = i;
      shift[CS_PARAM_BC_DIRICHLET] += 1;
      break;
    case CS_CDO_BC_HMG_DIRICHLET:
      bc->hmg_dir_ids[shift[CS_PARAM_BC_HMG_DIRICHLET]] = i;
      shift[CS_PARAM_BC_HMG_DIRICHLET] += 1;
      break;
    case CS_CDO_BC_NEUMANN:
      bc->nhmg_neu_ids[shift[CS_PARAM_BC_NEUMANN]] = i;
      shift[CS_PARAM_BC_NEUMANN] += 1;
    break;
    case CS_CDO_BC_HMG_NEUMANN:
      bc->hmg_neu_ids[shift[CS_PARAM_BC_HMG_NEUMANN]] = i;
      shift[CS_PARAM_BC_HMG_NEUMANN] += 1;
      break;
    case CS_CDO_BC_ROBIN:
      bc->robin_ids[shift[CS_PARAM_BC_ROBIN]] = i;
      shift[CS_PARAM_BC_ROBIN] += 1;
      break;

      /* For vector-valued equations only */
    case CS_CDO_BC_SLIDING:
      bc->sliding_ids[shift[CS_PARAM_BC_SLIDING]] = i;
      shift[CS_PARAM_BC_SLIDING] += 1;
      break;

      /* For vector-valued equations only */
    case CS_CDO_BC_TANGENTIAL_DIRICHLET:
      bc->circulation_ids[shift[CS_PARAM_BC_CIRCULATION]] = i;
      shift[CS_PARAM_BC_CIRCULATION] += 1;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: This type of boundary condition is not handled.",
                __func__);

    } /* End of switch */

  } /* Loop on boundary faces */

  return bc;
}

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
cs_cdo_bc_free(cs_cdo_bc_face_t   *face_bc)
{
  if (face_bc == NULL)
    return face_bc;

  BFT_FREE(face_bc->flag);
  BFT_FREE(face_bc->def_ids);

  BFT_FREE(face_bc->hmg_dir_ids);
  BFT_FREE(face_bc->nhmg_dir_ids);
  BFT_FREE(face_bc->hmg_neu_ids);
  BFT_FREE(face_bc->nhmg_neu_ids);
  BFT_FREE(face_bc->robin_ids);
  BFT_FREE(face_bc->sliding_ids);
  BFT_FREE(face_bc->circulation_ids);

  BFT_FREE(face_bc);

  return NULL;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
