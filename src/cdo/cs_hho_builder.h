#ifndef __CS_HHO_BUILDER_H__
#define __CS_HHO_BUILDER_H__

/*============================================================================
 * Low-level functions and structures used to build the algebraic system with
 * a cellwise process when Hybrid High Order schemes are set for the space
 * discretization
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
#include "cs_basis_func.h"
#include "cs_cdo_connect.h"
#include "cs_sdm.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Cellwise builder for HHO discretization */
typedef struct {

  /* Current and maximal number of face basis allocated. This number is equal
     to the max. number of faces for a cell */
  short int          n_face_basis;
  short int          n_max_face_basis;

  cs_basis_func_t  **face_basis;   /* P_(d-1)^k          polynomial basis */
  cs_basis_func_t   *cell_basis;   /* P_d^k              polynomial basis */
  cs_basis_func_t   *grad_basis;   /* P_d^(k+1) \ P_d^0  polynomial basis */

  cs_sdm_t   *grad_reco_op;  /* Gradient operator; Rectangular matrix */

  /* Temporary matrices defined by blocks */
  cs_sdm_t   *tmp;           /* Temporary block matrix (fs x ts) */
  cs_sdm_t   *bf_t;          /* Transposed  of Bf (used in stabilization) */
  cs_sdm_t   *jstab;         /* Stabilization part related to a face */

} cs_hho_builder_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate  a cs_hho_builder_t structure
 *
 * \param[in] order             order of the polynomial basis function
 * \param[in] n_fc              max. number of faces in a cell
 *
 * \return a pointer to a new allocated cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hho_builder_t *
cs_hho_builder_create(int     order,
                      int     n_fc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_hho_builder_t structure
 *
 * \param[in, out] p_builder  pointer of pointer on a cs_hho_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_free(cs_hho_builder_t  **p_builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set-up the basis functions related to a cell, its gradient and to
 *         the faces of this cell.
 *         Compute cell and face projection and its related modified Cholesky
 *         factorization.
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_cellwise_setup(const cs_cell_mesh_t    *cm,
                              cs_cell_builder_t       *cb,
                              cs_hho_builder_t        *hhob);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set-up the basis functions related to a cell only.
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_hho_builder_cellbasis_setup(const cs_cell_mesh_t    *cm,
                               cs_cell_builder_t       *cb,
                               cs_hho_builder_t        *hhob)
{
  if (hhob == NULL)
    return;

  /* Sanity checks */
  assert(cm != NULL);

  /* Setup cell basis functions */
  hhob->cell_basis->setup(hhob->cell_basis, cm, 0, cm->xc, cb);
  hhob->n_face_basis = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient operator stemming from the relation
 *            stiffness * grad_op = rhs
 *         where stiffness is a square matrix of size grd_size
 *               rhs is matrix of size (n_fc*f_size + c_size) * grd_size
 *         Hence, grad_op a matrix grd_size * (n_fc*f_size + c_size)
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_compute_grad_reco(const cs_cell_mesh_t    *cm,
                                 cs_cell_builder_t       *cb,
                                 cs_hho_builder_t        *hhob);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusion operator. The gradient reconstruction operator
 *         has to be built just before this call (cb->aux stores the rhs)
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_diffusion(const cs_cell_mesh_t    *cm,
                         cs_cell_builder_t       *cb,
                         cs_hho_builder_t        *hhob);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the polynomial spaces (cell and faces)
 *         of a function defined by an analytical expression depending on the
 *         location and the current time
 *         red array has to be allocated before calling this function
 *
 * \param[in]       def      pointer to a cs_xdef_t structure
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       t_eval   time at which one performs the evaluation
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 * \param[in, out]  red      vector containing the reduction
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_reduction_from_analytic(const cs_xdef_t         *def,
                                       const cs_cell_mesh_t    *cm,
                                       cs_real_t                t_eval,
                                       cs_cell_builder_t       *cb,
                                       cs_hho_builder_t        *hhob,
                                       cs_real_t                red[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the polynomial spaces (cell and faces)
 *         of a function defined by an analytical expression depending on the
 *         location and the current time
 *         This function handles the vector case.
 *
 *         red array has to be allocated before calling this function.
 *
 * \param[in]       def      pointer to a cs_xdef_t structure
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       t_eval   time at which one performs the evaluation
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 * \param[in, out]  red      vector containing the reduction
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_reduction_from_analytic_v(const cs_xdef_t         *def,
                                         const cs_cell_mesh_t    *cm,
                                         cs_real_t                t_eval,
                                         cs_cell_builder_t       *cb,
                                         cs_hho_builder_t        *hhob,
                                         cs_real_t                red[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projection of the Dirichlet boundary conditions onto
 *         the polynomial spaces on faces
 *
 * \param[in]       def      pointer to a cs_xdef_t structure
 * \param[in]       f        local face id in the cellwise view of the mesh
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       t_eval   time at which one performs the evaluation
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 * \param[in, out]  res      vector containing the result
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_compute_dirichlet(const cs_xdef_t         *def,
                                 short int                f,
                                 const cs_cell_mesh_t    *cm,
                                 cs_real_t                t_eval,
                                 cs_cell_builder_t       *cb,
                                 cs_hho_builder_t        *hhob,
                                 cs_real_t                res[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projection of the Dirichlet boundary conditions onto
 *         the polynomial spaces on faces. Vector case.
 *
 * \param[in]       def      pointer to a cs_xdef_t structure
 * \param[in]       f        local face id in the cellwise view of the mesh
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       t_eval   time at which one performs the evaluation
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 * \param[in, out]  res      vector containing the result
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_compute_dirichlet_v(const cs_xdef_t         *def,
                                   short int                f,
                                   const cs_cell_mesh_t    *cm,
                                   cs_real_t                t_eval,
                                   cs_cell_builder_t       *cb,
                                   cs_hho_builder_t        *hhob,
                                   cs_real_t                res[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HHO_BUILDER_H__ */
