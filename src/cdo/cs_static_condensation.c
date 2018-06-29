/*============================================================================
 * Routines to handle the static condensation
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

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_static_condensation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_STATIC_CONDENSATION_DBG  0

/*============================================================================
 * Local private variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and store
 *          information inside the rc_tilda and acx_tilda to be able to compute
 *          the values at cell centers
 *          rc_tilda = Acc^-1 * cell_rhs
 *          Case of scalar-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in, out] rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in, out] acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_scalar_eq(const cs_adjacency_t    *c2x,
                                 cs_real_t               *rc_tilda,
                                 cs_real_t               *acx_tilda,
                                 cs_cell_builder_t       *cb,
                                 cs_cell_sys_t           *csys)
{
  const int  n_dofs = csys->n_dofs;
  const int  n_xc = n_dofs - 1;
  const double  cell_rhs = csys->rhs[n_xc];
  const double  *row_c = csys->mat->val + n_dofs*n_xc;
  assert(fabs(row_c[n_xc]) > cs_math_zero_threshold);
  const double  inv_acc = 1./row_c[n_xc];

  /* Compute and store rc_tilda */
  rc_tilda[csys->c_id] = inv_acc * cell_rhs;

  /* Compute and store acx_tilda */
  double  *acx = acx_tilda + c2x->idx[csys->c_id];
  for (int i = 0; i < n_xc; i++) acx[i] = inv_acc * row_c[i];

  /* Temporary storage of axc (part of the last column) */
  double  *axc = cb->values;
  for (int i = 0; i < n_xc; i++) axc[i] = csys->mat->val[n_dofs*i + n_xc];

  /* Update matrix and rhs */
  csys->n_dofs = n_xc;
  csys->mat->n_rows = csys->mat->n_cols = n_xc;

  for (short int i = 0; i < n_xc; i++) {

    double  *old_i = csys->mat->val + n_dofs*i; /* Old "i" row  */
    double  *new_i = csys->mat->val + n_xc*i;   /* New "i" row */

    /* Condensate the local matrix Axx:
       Axx --> Axx - Axc.Acc^-1.Acx */
    for (short int j = 0; j < n_xc; j++)
      new_i[j] = old_i[j] - axc[i]*acx[j];

    /* Update RHS_x: RHS_x = RHS_x - Axc*Acc^-1*RHS_c */
    csys->rhs[i] -= rc_tilda[csys->c_id] * axc[i];

  } // Loop on vi cell vertices

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Opposite process of the static condensation.
 *          Define the field at cells given the field at x locations and arrays
 *          storing the static condensation.
 *          Case of scalar-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in]      rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in]      acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in]      px          values of the fields at x locations
 * \param[in, out] pc          values of the field at cells
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_recover_scalar(const cs_adjacency_t    *c2x,
                                      const cs_real_t         *rc_tilda,
                                      const cs_real_t         *acx_tilda,
                                      const cs_real_t         *px,
                                      cs_real_t               *pc)
{
  assert(pc != NULL && px != NULL);

  const cs_lnum_t  n_cells = c2x->n_elts;;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_lnum_t  shift_c = c2x->idx[c_id];
    const int  n_xc = c2x->idx[c_id+1] - shift_c;
    const cs_lnum_t  *x_ids = c2x->ids + shift_c;
    const cs_real_t  *acx = acx_tilda + shift_c;

    double  acx_px = 0.;
    for (int i = 0; i < n_xc; i++) acx_px += acx[i]*px[x_ids[i]];

    pc[c_id] = rc_tilda[c_id] - acx_px;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and store
 *          information inside the rc_tilda and acx_tilda to be able to compute
 *          the values at cell centers
 *          rc_tilda = Acc^-1 * cell_rhs
 *          Case of vector-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in, out] rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in, out] acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_vector_eq(const cs_adjacency_t    *c2x,
                                 cs_real_t               *rc_tilda,
                                 cs_real_t               *acx_tilda,
                                 cs_cell_builder_t       *cb,
                                 cs_cell_sys_t           *csys)
{
  cs_sdm_t  *m = csys->mat;
  cs_sdm_block_t  *bd = m->block_desc;

  const int  stride = 3;
  const int  diag = stride + 1;
  const int  n_dofs = bd->n_row_blocks;
  const int  n_xc = n_dofs - 1;
  const double  *cell_rhs = csys->rhs + stride*n_xc;

  /* mCC is a small square matrix of size stride x stride
     (should be diagonal) */
  const cs_sdm_t  *mcc = cs_sdm_get_block(m, n_xc, n_xc);

  /* Compute and store rc_tilda (allocated to stride*n_cells) */
  cs_real_t  *_rc_tilda = rc_tilda + stride*csys->c_id;
  for (int i = 0; i < stride; i++)
    _rc_tilda[i] = cell_rhs[i]/mcc->val[diag*i];

  /* Compute and store acx_tilda */
  cs_real_t  *acx = acx_tilda + stride*c2x->idx[csys->c_id]; /* Only diagonal */
  for (int ix = 0; ix < n_xc; ix++) {
    const cs_sdm_t  *mcx = cs_sdm_get_block(m, n_xc, ix);
    for (int k = 0; k < stride; k++)
      acx[stride*ix + k] = mcx->val[diag*k]/mcc->val[diag*k];
  }

  /* Temporary storage of axc (part of the last column) */
  cs_real_t  *axc = cb->values;
  for (int ix = 0; ix < n_xc; ix++) {
    const cs_sdm_t  *mxc = cs_sdm_get_block(m, ix, n_xc);
    for (int k = 0; k < stride; k++)
      axc[stride*ix + k] = mxc->val[diag*k];
  }

  /* Update matrix and rhs */
  csys->n_dofs = stride*n_xc;
  for (short int bfi = 0; bfi < n_xc; bfi++) {

    const cs_real_t  *axc_i = axc + stride*bfi;
    for (short int bfj = 0; bfj < n_xc; bfj++) {

      cs_sdm_t  *mxx = cs_sdm_get_block(m, bfi, bfj);

      /* Condensate the local block mxx:
         mxx --> mxx - mxc.mcc^-1.mcx */
      for (int k = 0; k < stride; k++)
        mxx->val[diag*k] -= axc_i[k] * acx[stride*bfj+k];

    } /* Loop on blocks for face fj */

    /* Update RHS: RHS_x = RHS_x - mxc*mcc^-1*RHS_c */
    for (int k = 0; k < stride; k++)
      csys->rhs[3*bfi+k] -= _rc_tilda[k] * axc_i[k];

  } /* Loop on blocks for face fi */

  /* Reshape matrix */
  int  shift = n_xc;
  for (short int bfi = 1; bfi < n_xc; bfi++) {
    for (short int bfj = 0; bfj < n_xc; bfj++) {

      cs_sdm_t  *mxx_old = cs_sdm_get_block(m, bfi, bfj);

      /* Set the block (i,j) */
      cs_sdm_t  *mxx = bd->blocks + shift;

      cs_sdm_copy(mxx, mxx_old);
      shift++;

    }
  }

  m->n_rows = m->n_cols = stride * n_xc;
  bd->n_row_blocks = n_xc;      /* instead of n_xc + 1 */
  bd->n_col_blocks = n_xc;      /* instead of n_xc + 1 */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Opposite process of the static condensation.
 *          Define the field at cells given the field at x locations and arrays
 *          storing the static condensation.
 *          Case of vector-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in]      rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in]      acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in]      px          values of the fields at x locations
 * \param[in, out] pc          values of the field at cells
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_recover_vector(const cs_adjacency_t    *c2x,
                                      const cs_real_t         *rc_tilda,
                                      const cs_real_t         *acx_tilda,
                                      const cs_real_t         *px,
                                      cs_real_t               *pc)
{
  assert(pc != NULL && px != NULL);
  const cs_lnum_t  n_cells = c2x->n_elts;;
  const int  stride = 3;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Compute acx_px */
    cs_real_t  acx_px[3] = {0., 0., 0.}; /* Allocated to stride */
    for (cs_lnum_t i = c2x->idx[c_id]; i < c2x->idx[c_id+1]; i++) {
      for (int k = 0; k < stride; k++)
        acx_px[k] += px[stride*c2x->ids[i]+k] * acx_tilda[stride*i+k];
    }

    for (int k = 0; k < stride; k++)
      pc[stride*c_id+k] = rc_tilda[stride*c_id+k] - acx_px[k];

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
