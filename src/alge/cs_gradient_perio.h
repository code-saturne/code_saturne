#ifndef __CS_GRADIENT_PERIO_H__
#define __CS_GRADIENT_PERIO_H__

/*============================================================================
 * Gradient reconstruction.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Periodicity management for INIMAS
 *
 * If INIMAS is called by NAVSTO :
 *    We assume that gradient on ghost cells given by a rotation is known
 *    and is equal to the velocity one for the previous time step.
 * If INIMAS is called by DIVRIJ
 *    We assume that (more justifiable than in previous case) gradient on
 *    ghost cells given by rotation is equal to Rij gradient for the previous
 *    time step.
 *
 * Fortran Interface:
 *
 * subroutine permas
 * *****************
 *
 * integer          iappel      :  -> : indicateur d'appel dans inimas
 *                                          = 1 si appel au debut
 *                                          = 2 si appel a la fin
 * double precision rom(ncelet) :  -> : masse volumique aux cellules
 *
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (permas, PERMAS)(const cs_int_t  *iappel,
                          cs_real_t        rom[]);

/*----------------------------------------------------------------------------
 * Preparation rotation periodicity for Reynolds stresses.
 *
 * Compute an estimation of the velocity gradient.
 * The gradient is not reconstructed (as we would need a gradient,
 * thus periodicity). I seems possible to use a least-squares gradient.
 *
 * The gradient is then saved in an array representing the ghost cells, then
 * rotated where necessary to be ready for use (cf. cs_gradient_perio_init_rij).
 *
 * Compute cell gradient of vector field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (perinr, PERINR)
(
 const cs_int_t   *const imrgra,  /* <-- gradient computation mode            */
 const cs_int_t   *const iwarnp,  /* <-- verbosity level                      */
 const cs_real_t  *const epsrgp,  /* <-- precision for iterative gradient
                                         calculation                          */
 const cs_real_t  *const extrap   /* <-- extrapolate gradient at boundary     */
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_perio_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_perio_finalize(void);

/*----------------------------------------------------------------------------
 * Update gradient rotational periodicity computation API in case of
 * mesh modification.
 *----------------------------------------------------------------------------*/

void
cs_gradient_perio_update_mesh(void);

/*----------------------------------------------------------------------------
 * Initialize ghost cell values for Reynolds stress tensor gradient.
 *
 * We retrieve the gradient given by perinr (phyvar) for the Reynolds
 * stress tensor in a buffer on ghost cells. then we define
 * dpdx, dpdy and dpdz gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a tensor
 * variable because we have to know all components.
 *
 * We set idimtr to 2 for the Reynolds stress tensor.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables when reconstructing gradients.
 *
 * parameters:
 *   f      <-- pointer to field
 *   tr_dim --> 2 for tensor (Rij) in case of rotation, 0 otherwise
 *   grad   <-> gradient of field
 *----------------------------------------------------------------------------*/

void
cs_gradient_perio_init_rij(const cs_field_t  *f,
                           int               *tr_dim,
                           cs_real_3_t        grad[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize ghost cell values for Reynolds stress tensor gradient.
 *
 * We retrieve the gradient given by perinr (phyvar) for the Reynolds
 * stress tensor in a buffer on ghost cells. then we define
 * dpdx, dpdy and dpdz gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a tensor
 * variable because we have to know all components.
 *
 * We set idimtr to 2 for the Reynolds stress tensor.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables when reconstructing gradients.
 *
 * \param[out]      tr_dim   2 for tensor (Rij) in case of rotation, 0 otherwise
 * \param[in, out]  grad     gradient of field
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_init_rij_tensor(int           *tr_dim,
                                  cs_real_63_t   grad[]);

/*----------------------------------------------------------------------------*/
/*!
 * Process grad buffers in case of rotation on Reynolds stress tensor.
 *
 * We retrieve the gradient given by cs_gradient_perio_init_rij (phyvar)
 * for the Reynolds stress tensor in a buffer on ghost cells. Then we define
 * grad gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a tensor
 * variable because we have to know all components.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables when reconstructing gradients.
 *
 * parameters:
 *   f_id      <--   field index
 *   grad      -->   gradient of field
 *
 * size of _drdxyz and _wdrdxy = n_ghost_cells*6*3
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_process_rij(const cs_int_t    *f_id,
                              cs_real_3_t        grad[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GRADIENT_PERIO__ */
