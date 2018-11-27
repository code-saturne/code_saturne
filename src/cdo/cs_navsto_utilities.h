#ifndef __NAVSTO_UTILITIES_H__
#define __NAVSTO_UTILITIES_H__

/*============================================================================
 * Routines common to all the velocity-pressure couplings
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_math.h"
#include "cs_sdm.h"

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence of a cell using the \ref cs_cdo_quantities_t
 *         structure
 *
 * \param[in]     c_id         cell ID
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_dof        values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_navsto_get_cell_divergence(const cs_lnum_t               c_id,
                              const cs_cdo_quantities_t    *quant,
                              const cs_adjacency_t         *c2f,
                              const cs_real_t              *f_dof)
{
  cs_real_t  div = 0.0;

  for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

    const cs_lnum_t  f_id = c2f->ids[f];
    const cs_real_t  *_val = f_dof + 3*f_id;

    if (f_id < quant->n_i_faces) {
      const cs_real_t *_nuf = quant->i_face_normal + 3*f_id;

      div += c2f->sgn[f]*quant->i_face_surf[f_id]*
        cs_math_3_dot_product(_val, _nuf) / cs_math_3_norm(_nuf);

    }
    else {

      const cs_lnum_t  bf_id = f_id - quant->n_i_faces;
      const cs_real_t  *_nuf = quant->b_face_normal + 3*bf_id;

      div += c2f->sgn[f]*quant->b_face_surf[bf_id]*
        cs_math_3_dot_product(_val, _nuf) / cs_math_3_norm(_nuf);

    } /* Boundary face */

  } /* Loop on cell faces */

  div /= quant->cell_vol[c_id];

  return div;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence vector associated to the current cell.
 *         WARNING: mind that, differently form the original definition, the
 *         result here is not divided by the cell volume
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out] div        array related to the divergence operator
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_navsto_get_divergence_vect(const cs_cell_mesh_t  *cm,
                              cs_real_t              div[])
{
  /* D(\hat{u}) = \frac{1}{|c|} \sum_{f_c} \iota_{fc} u_f.f
   * But, when integrating [[ p, q ]]_{P_c} = |c| p_c q_c
   * Thus, the volume in the divergence drops
   */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_real_t  i_f = cm->f_sgn[f] * pfq.meas;

    cs_real_t  *_div_f = div + 3*f;
    _div_f[0] = i_f * pfq.unitv[0];
    _div_f[1] = i_f * pfq.unitv[1];
    _div_f[2] = i_f * pfq.unitv[2];

  } /* Loop on cell faces */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the grad-div part to the local matrix (i.e. for the current
 *         cell)
 *
 * \param[in]      n_fc       local number of faces for the current cell
 * \param[in]      zeta       scalar coefficient for the grad-div operator
 * \param[in]      div        divergence
 * \param[in, out] mat        local system matrix to update
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_grad_div(short int          n_fc,
                       const cs_real_t    zeta,
                       const cs_real_t    div[],
                       cs_sdm_t          *mat);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __NAVSTO_UTILITIES_H__ */
