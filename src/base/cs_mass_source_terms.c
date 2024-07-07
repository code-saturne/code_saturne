/*============================================================================
 * Mass source terms computation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mass_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_mass_source_terms.c
        Mass source terms computation.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Implicit and explicit mass source terms computation.
 *
 * Arrays st_exp, st_imp, and gapinj are incremented, so should be initialized.
 *
 * \param[in]     iterns        iteration number on Navier-Stoke
 * \param[in]     dim           associated field dimension
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     icetsm        source mass cells pointer (1-based numbering)
 * \param[in]     itpsmp        mass source type for the working variable
 * \param[in]     volume        cells volume
 * \param[in]     pvara         variable value at time step beginning
 * \param[in]     smcelp        value of the variable associated with mass source
 * \param[in]     gamma         flow mass value
 * \param[in,out] st_exp        explicit source term part linear in the variable
 * \param[in,out] st_imp        associated value with \c tsexp
 *                              to be stored in the matrix
 * \param[in,out] gapinj        explicit source term part independent
 *                              of the variable
 */
/*----------------------------------------------------------------------------*/

void
cs_mass_source_terms(int                   iterns,
                     int                   dim,
                     cs_lnum_t             ncesmp,
                     const cs_lnum_t       icetsm[],
                     int                   itpsmp[],
                     const cs_real_t       volume[],
                     const cs_real_t       pvara[],
                     const cs_real_t       smcelp[],
                     const cs_real_t       gamma[],
                     cs_real_t             st_exp[],
                     cs_real_t             st_imp[],
                     cs_real_t             gapinj[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;

  if (gamma == NULL) /* No mass source term here; should not occur */
    return;

  /* Remark for tests on gamma[i] > O && itpsmp[i] == 1 :
     *
     * If we remove matter or enter with the cell value
     * then the equation on the variable has not been modified.
     *
     * Otherwise, we add the term gamma*(f_i-f^(n+1))
     *
     * In st_imp, we add the term that will go on the diagonal, which is Gamma.
     * In st_exp, we add the term for the right-hand side, which is
     *   Gamma * Pvar
     *
     * In gapinj, we place the term Gamma Pinj which will go to the right-hand side.
     *
     * The distinction between st_exp and w1 (which both go finally to the
     * right-hand side) is used for the 2nd-order time scheme. */

  if (iterns == 1) {

    if (smcelp != NULL && gapinj != NULL) {
      if (dim == 1) {
        for (cs_lnum_t i = 0; i < ncesmp; i++) {
          cs_lnum_t c_id = icetsm[i];
          if (gamma[i] > 0. && itpsmp[i] == 1) {
            st_exp[c_id] -= volume[c_id]*gamma[i] * pvara[c_id];
            gapinj[c_id] += volume[c_id]*gamma[i] * smcelp[i];
          }
        }
      }
      else {
        cs_lnum_t _dim = dim;
        for (cs_lnum_t i = 0; i < ncesmp; i++) {
          cs_lnum_t c_id = icetsm[i];
          if (gamma[i] > 0. && itpsmp[i] == 1) {
            for (cs_lnum_t j = 0; j < _dim; j++) {
              cs_lnum_t k = c_id*_dim + j;
              st_exp[k] -= volume[c_id]*gamma[i] * pvara[k];
              gapinj[k] += volume[c_id]*gamma[i] * smcelp[i*_dim + j];
            }
          }
        }
      }
    }

    else { /* smcelp == NULL && gapinj == NULL */
      if (dim == 1) {
        for (cs_lnum_t i = 0; i < ncesmp; i++) {
          cs_lnum_t c_id = icetsm[i];
          if (gamma[i] > 0. && itpsmp[i] == 1)
            st_exp[c_id] -= volume[c_id]*gamma[i] * pvara[c_id];
        }
      }
      else {
        cs_lnum_t _dim = dim;
        for (cs_lnum_t i = 0; i < ncesmp; i++) {
          cs_lnum_t c_id = icetsm[i];
          if (gamma[i] > 0. && itpsmp[i] == 1) {
            for (cs_lnum_t j = 0; j < _dim; j++) {
              cs_lnum_t k = c_id*_dim + j;
              st_exp[k] -= volume[c_id]*gamma[i] * pvara[k];
            }
          }
        }
      }
    }

  }

  /* On the diagonal */

  if (dim == 1) {
    for (cs_lnum_t i = 0; i < ncesmp; i++) {
      cs_lnum_t c_id = icetsm[i];
      if (gamma[i] > 0. && itpsmp[i] == 1) {
        st_imp[c_id] += volume[c_id]*gamma[i];
      }
    }
  }
  else {
    cs_lnum_t _dim = dim, _dim2 = dim*dim;
    for (cs_lnum_t i = 0; i < ncesmp; i++) {
      cs_lnum_t c_id = icetsm[i];
      if (gamma[i] > 0. && itpsmp[i] == 1) {
        for (cs_lnum_t j = 0; j < _dim; j++) {
          cs_lnum_t k = c_id*_dim2 + j*_dim + j;
          st_imp[k] += volume[c_id]*gamma[i];
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
