#ifndef __CS_BASIS_FUNC_H__
#define __CS_BASIS_FUNC_H__

/*============================================================================
 * Build a set of basis functions for cells and faces and cell gradients
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"
#include "cs_quadrature.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Default basis function are rescale and modified according to the
   geometrical inertia tensor of the related entity unless a flag "monomial"
   is set */

#define CS_BASIS_FUNC_MONOMIAL  (1 << 0) /* 1: Use mononial basis functions */
#define CS_BASIS_FUNC_GRADIENT  (1 << 1) /* 2: Basis functions for gradient */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for evaluating all basis functions at a given
 *         point
 *
 * \param[in]      bf      pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_eval_all_at_point_t) (const void           *bf,
                                     const cs_real_t       coords[3],
                                     cs_real_t            *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for evaluating a set of basis functions at a
 *         given point.
 *
 * \param[in]      bf      pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id-1
 * \param[in]      end     ends evaluating basis function at this id
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_eval_at_point_t) (const void           *bf,
                                 const cs_real_t       coords[3],
                                 short int             start,
                                 short int             end,
                                 cs_real_t            *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for setting up a set of basis functions
 *
 * \param[in, out]  pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]       cm      pointer to a cs_cell_mesh_t
 * \param[in]       id      id of the element to consider
 * \param[in]       center  point used for centering the set of basis functions
 * \param[in, out]  cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_setup_t) (void                    *pbf,
                         const cs_cell_mesh_t    *cm,
                         const short int          id,
                         const cs_real_t          center[3],
                         cs_cell_builder_t       *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for computing the projector to the space spanned
 *         by the basis functions
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      id      id of the element to consider
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_compute_proj_t) (void                    *pbf,
                                const cs_cell_mesh_t    *cm,
                                const short int          id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for computing the Modified Choloesky factorization
 *         of the projection matrix (mass matrix) related to the basis function
 *
 * \param[in, out]  pbf    pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_compute_facto_t) (void              *pbf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for projecting an array on the polynomial basis
 *         function. This results from the application of a Modified Choloesky
 *         factorization which should be performed before calling this function.
 *         The input array is defined as follows:
 *         int_elem v.phi_i for all i in the basis. v is a function to estimate
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      array   array to project
 * \param[in, out] dof     projection of the array (= DoF values in this basis)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_project_t) (const void              *pbf,
                           const cs_real_t         *array,
                           cs_real_t               *dof);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for printing the projector matrix related to the
 *         given basis function
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_basis_func_dump_proj_t) (const void              *pbf);


/* Structure storing information to build and use a set of basis functions
 *
 * - When dealing with volumetric basis fucntions (eg cell- or grad-basis),
 *   the functions are ordered as follows:
 *   1, x, y, z, xx, xy, xz, yy, yz, zz
 * - When dealing with face-basis fucntions the functions are ordered as
 *   follows:
 *   1, x, y, xx, xy, yy
 */

typedef struct {

  cs_flag_t     flag;        // metadata
  short int     poly_order;  // max. polynomial order of the basis
  short int     dim;         // 2D or 3D associated geometrical entities
  int           size;        // number of elementary basis functions

  cs_real_t     phi0;        // constant basis function
  cs_nvec3_t   *axis;        // over_diam * main axis of the entity
  cs_real_3_t   center;      // center used for rescaling

  /* For poly_order > 1 we store the degree related to each monone for
     each basis function */
  int           n_deg_elts;
  short int    *deg;            /* size = n_deg_elts * dim */

  /* Function used for setting up the quantities used for defining
     the set of basis functions */
  cs_basis_func_setup_t               *setup;

  /* Function used for evaluating basis functions on a set of points */
  cs_basis_func_eval_all_at_point_t   *eval_all_at_point;
  cs_basis_func_eval_at_point_t       *eval_at_point;

  /* Projector is a mass matrix related to the basis function */
  cs_sdm_t                            *projector;
  cs_basis_func_compute_proj_t        *compute_projector;
  cs_basis_func_dump_proj_t           *dump_projector;

  /* Modified Cholesky factorization. To perform the inversion agains a set
     of right-hand sides */
  cs_basis_func_compute_facto_t       *compute_factorization;
  cs_basis_func_project_t             *project;
  cs_real_t                           *facto;
  int                                  facto_max_size;

  /* Quadrature function to use for computing integral over o volume (tet) or
     on a surface (tria) */
  int                                  n_gpts_tria;
  cs_quadrature_tria_t                *quadrature_tria;
  int                                  n_gpts_tetra;
  cs_quadrature_tet_t                 *quadrature_tetra;

} cs_basis_func_t;

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the polynomial order of the given basis
 *
 * \param[in]  bf  set of basis functions
 *
 * \return the polynomial order
 */
/*----------------------------------------------------------------------------*/

static inline short int
cs_basis_func_get_poly_order(const cs_basis_func_t   *bf)
{
  if (bf == NULL)
    return -1;
  else
    return bf->poly_order;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_basis_func_t structure
 *
 * \param[in]  flag     metadata related to the way of building basis functions
 * \param[in]  k        order
 * \param[in]  dim      2 or 3 geometrical dimension
 *
 * \return a pointer to the new cs_basis_func_t
 */
/*----------------------------------------------------------------------------*/

cs_basis_func_t *
cs_basis_func_create(cs_flag_t      flag,
                     short int      k,
                     short int      dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_basis_func_t structure which is associated to an
 *         existing set of basis functions.
 *         Up to now, only cell basis functions are handled.
 *         Building a projection matrix is not possible in this case.
 *
 * \param[in]  ref  set of basis function used as a reference
 *
 * \return a pointer to the new cs_basis_func_t for gradient of the current
 *         basis functions
 */
/*----------------------------------------------------------------------------*/

cs_basis_func_t *
cs_basis_func_grad_create(const cs_basis_func_t   *ref);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the center and the different axis from the reference basis
 *         Up to now, only cell basis functions are handled.
 *
 * \param[in]      ref   set of basis function used as a reference
 * \param[in, out] rcv   set of basis function where members are set
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_copy_setup(const cs_basis_func_t   *ref,
                         cs_basis_func_t         *rcv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_basis_func_t structure
 *
 * \param[in, out]  pbf   pointer to the cs_basis_func_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_basis_func_t *
cs_basis_func_free(cs_basis_func_t  *pbf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set options for basis functions when using HHO schemes
 *
 * \param[in]  face_flag    options related to face basis functinos
 * \param[in]  cell_flag    options related to cell basis functinos
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_set_hho_flag(cs_flag_t   face_flag,
                           cs_flag_t   cell_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get options for basis functions when using HHO schemes
 *
 * \param[out] face_flag   pointer to options related to face basis functinos
 * \param[out] cell_flag   pointer to options related to cell basis functinos
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_get_hho_flag(cs_flag_t   *face_flag,
                           cs_flag_t   *cell_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_basis_func_t structure
 *
 * \param[in]  pbf   pointer to the cs_basis_func_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_dump(const cs_basis_func_t  *pbf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a cs_basis_func_t structure
 *         Print into the file f if given otherwise open a new file named
 *         fname if given otherwise print into the standard output
 *
 * \param[in]  fp      pointer to a file structure or NULL
 * \param[in]  fname   filename or NULL
 * \param[in]  pbf     pointer to the cs_basis_func_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_fprintf(FILE                   *fp,
                      const char             *fname,
                      const cs_basis_func_t  *pbf);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASIS_FUNC_H__ */
