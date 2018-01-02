/*============================================================================
 * Manage boundary conditions
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

#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_bc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
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
_get_bc_flag(cs_param_bc_type_t   bc_type)
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
  default:
    ret_flag = 0;
    break;
  }
  return ret_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_cdo_bc_t structure
 *
 * \param[in] n_elts   number of elements
 *
 * \return  a new allocated pointer to a cs_cdo_bc_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cdo_bc_t *
_cdo_bc_create(cs_lnum_t  n_elts)
{
  cs_cdo_bc_t  *bc = NULL;

  BFT_MALLOC(bc, 1, cs_cdo_bc_t);

  bc->n_elts = n_elts;

  BFT_MALLOC(bc->flag, n_elts, cs_flag_t);

  /* Default initialization */
  for (cs_lnum_t i = 0; i < n_elts; i++)
    bc->flag[i] = 0;     /* Nothing set. Process all the definition before
                            setting something. */

  bc->dir = NULL;
  bc->neu = NULL;
  bc->rob = NULL;

  return bc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the cs_cdo_bc_t structure with elements associated to the
 *          definition of id def_id
 *
 * \param[in]      def_id    id of the definition to add
 * \param[in]      n_faces   number of border faces to specify
 * \param[in]      elt_ids   list of border faces related to this def. (or NULL)
 * \param[in]      shift     shift to apply before adding new entries
 * \param[in, out] bc_defs   pointer to the list of definitions
 * \param[in, out] bc_list   pointer to the list dedicated to a specific type
 *                           of boundary condition
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_def_to_bc(const short int    def_id,
               const cs_lnum_t    n_faces,
               const cs_lnum_t   *elt_ids,
               cs_lnum_t          shift,
               short int         *bc_defs,
               cs_lnum_t         *bc_list)
{
  short int  *defs = bc_defs + shift;
  cs_lnum_t  *ids = bc_list + shift;

  if (elt_ids == NULL) { // This definition encompasses all the border faces

    assert(shift == 0); // Sanity check
    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {
      ids[f_id] = f_id;
      defs[f_id] = def_id;
    }

  }
  else { // List of faces related to a mesh location

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      const cs_lnum_t  f_id = elt_ids[i];
      ids[i] = f_id;
      defs[i] = def_id;
    }

  } /* elt_ids != NULL */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_cdo_bc_list_t structure
 *
 * \param[in]  n_elts      number of entries of the list
 * \param[in]  n_nhmg_elts number of elements attached to a homogeneous BC
 *
 * \return  a new allocated pointer to a cs_cdo_bc_list_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_list_create(cs_lnum_t   n_elts,
                      cs_lnum_t   n_nhmg_elts)
{
  cs_cdo_bc_list_t  *bcl = NULL;

  /* Sanity check */
  assert(n_elts >= n_nhmg_elts && n_elts > -1);

  BFT_MALLOC(bcl, 1, cs_cdo_bc_list_t);

  bcl->n_elts = n_elts;
  bcl->n_nhmg_elts = n_nhmg_elts;

  /* Allocate and initialize by default */
  bcl->elt_ids = NULL;
  if (n_elts > 0) {
    BFT_MALLOC(bcl->elt_ids, n_elts, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_elts; i++)
      bcl->elt_ids[i] = -1;
  }

  bcl->def_ids = NULL;
  if (n_nhmg_elts > 0) {
    BFT_MALLOC(bcl->def_ids, n_nhmg_elts, short int);
    for (cs_lnum_t i = 0; i < n_nhmg_elts; i++)
      bcl->def_ids[i] = -1;
  }

  return bcl;
}

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
cs_cdo_bc_list_free(cs_cdo_bc_list_t   *bcl)
{
  if (bcl == NULL)
    return bcl;

  if (bcl->n_elts > 0)
    BFT_FREE(bcl->elt_ids);

  if (bcl->n_nhmg_elts > 0)
    BFT_FREE(bcl->def_ids);

  BFT_FREE(bcl);
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the structure which translates the BC definition from the
 *         user viewpoint into a ready-to-use structure
 *         - Prepare the treatment of the boundary conditions.
 *         - Compile the information detailed in a cs_param_bc_t structure
 *
 * \param[in] param_bc    pointer to the parameters related to BCs of an eq.
 * \param[in] n_b_faces   number of border faces
 *
 * \return a pointer to a new allocated cs_cdo_bc_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_t *
cs_cdo_bc_define(const cs_param_bc_t  *param_bc,
                 cs_lnum_t             n_b_faces)
{
  /* Sanity checks */
  assert(param_bc != NULL);

  /* Set the default flag */
  cs_flag_t  default_flag = 0;
  if (param_bc->default_bc == CS_PARAM_BC_HMG_DIRICHLET)
    default_flag = CS_CDO_BC_HMG_DIRICHLET;
  else if (param_bc->default_bc == CS_PARAM_BC_HMG_NEUMANN)
    default_flag = CS_CDO_BC_HMG_NEUMANN;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Incompatible type of boundary condition by default.\n"
                " Please modify your settings.\n"));

  cs_cdo_bc_t  *bc = _cdo_bc_create(n_b_faces);

  if (n_b_faces == 0) {

    bc->dir = cs_cdo_bc_list_create(0, 0);
    bc->neu = cs_cdo_bc_list_create(0, 0);
    bc->rob = cs_cdo_bc_list_create(0, 0);

    return  bc;
  }

  /* Loop on the definition of each boundary condition */
  for (int ii = 0; ii < param_bc->n_defs; ii++) {

    const int  ml_id = param_bc->ml_ids[ii];
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
    const cs_flag_t  def_flag = _get_bc_flag(param_bc->types[ii]);

    if (elt_ids == NULL) {
      assert(n_elts[0] == n_b_faces);
      for (cs_lnum_t i = 0; i < n_elts[0]; i++)
        bc->flag[i] |= def_flag;
    }
    else
      for (cs_lnum_t i = 0; i < n_elts[0]; i++)
        bc->flag[elt_ids[i]] |= def_flag;

  } // Loop on boundary conditions

  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    if (bc->flag[i] == 0)
      bc->flag[i] = default_flag;

  /* Allocate and then define the cs_cdo_bc_list_t structures */
  cs_lnum_t  n_dir_faces = 0, n_hmg_dir_faces = 0;
  cs_lnum_t  n_neu_faces = 0, n_hmg_neu_faces = 0;
  cs_lnum_t  n_rob_faces = 0;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    if (bc->flag[i] & CS_CDO_BC_DIRICHLET)
      n_dir_faces++;
    if (bc->flag[i] & CS_CDO_BC_HMG_DIRICHLET)
      n_hmg_dir_faces++;
    if (bc->flag[i] & CS_CDO_BC_NEUMANN)
      n_neu_faces++;
    if (bc->flag[i] & CS_CDO_BC_HMG_NEUMANN)
      n_hmg_neu_faces++;
    if (bc->flag[i] & CS_CDO_BC_ROBIN)
      n_rob_faces++;
  }

  const cs_lnum_t  n_all_dir_faces = n_dir_faces + n_hmg_dir_faces;
  const cs_lnum_t  n_all_neu_faces = n_neu_faces + n_hmg_neu_faces;

  bc->dir = cs_cdo_bc_list_create(n_all_dir_faces,
                                  n_all_dir_faces - n_hmg_dir_faces);
  bc->neu = cs_cdo_bc_list_create(n_all_neu_faces,
                                  n_all_neu_faces - n_hmg_neu_faces);
  bc->rob = cs_cdo_bc_list_create(n_rob_faces, n_rob_faces);

  /* Sanity check (there is no multiple definition) */
  assert(n_b_faces == n_all_dir_faces + n_all_neu_faces + n_rob_faces);

  cs_lnum_t  shift[CS_PARAM_N_BC_TYPES];
  for (int ii = 0; ii < CS_PARAM_N_BC_TYPES; ii++)
    shift[ii] = 0;

  /* Manage the case of homogeneous BCs */
  cs_lnum_t  *dir_faces = bc->dir->elt_ids + bc->dir->n_nhmg_elts;
  cs_lnum_t  *neu_faces = bc->neu->elt_ids + bc->neu->n_nhmg_elts;

  /* Not useful to update the list of definitions for homogeneous BCs */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    if (bc->flag[f_id] & CS_CDO_BC_HMG_DIRICHLET) {
      dir_faces[shift[CS_PARAM_BC_HMG_DIRICHLET]] = f_id;
      shift[CS_PARAM_BC_HMG_DIRICHLET] += 1;
    }

    if (bc->flag[f_id] & CS_CDO_BC_HMG_NEUMANN) {
      neu_faces[shift[CS_PARAM_BC_HMG_NEUMANN]] = f_id;
      shift[CS_PARAM_BC_HMG_NEUMANN] += 1;
    }

  } // Loop on border faces

  /* Loop on the definition of each boundary condition */
  for (short int def_id = 0; def_id < param_bc->n_defs; def_id++) {

    const int  ml_id = param_bc->ml_ids[def_id];
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
    const cs_param_bc_type_t  type = param_bc->types[def_id];

    switch (type) {

    case CS_PARAM_BC_DIRICHLET:
      _add_def_to_bc(def_id, n_elts[0], elt_ids, shift[type],
                     bc->dir->def_ids, bc->dir->elt_ids);
      shift[type] += n_elts[0];
      break;

    case CS_PARAM_BC_NEUMANN:
      _add_def_to_bc(def_id, n_elts[0], elt_ids, shift[type],
                     bc->neu->def_ids, bc->neu->elt_ids);
      shift[type] += n_elts[0];
      break;

    case CS_PARAM_BC_ROBIN:
      _add_def_to_bc(def_id, n_elts[0], elt_ids, shift[type],
                     bc->rob->def_ids, bc->rob->elt_ids);
      shift[type] += n_elts[0];
      break;

    case CS_PARAM_BC_HMG_NEUMANN:
    case CS_PARAM_BC_HMG_DIRICHLET:
      break; // Nothing to do. Already treated.

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of boundary condition.\n"
                  " Stop generating the boundary condition structure."));
    }

  } // Loop on definitions of boundary conditions

  return bc;
}

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
cs_cdo_bc_free(cs_cdo_bc_t   *face_bc)
{
  if (face_bc == NULL)
    return face_bc;

  face_bc->dir = cs_cdo_bc_list_free(face_bc->dir);
  face_bc->neu = cs_cdo_bc_list_free(face_bc->neu);
  face_bc->rob = cs_cdo_bc_list_free(face_bc->rob);

  BFT_FREE(face_bc->flag);
  BFT_FREE(face_bc);

  return NULL;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
