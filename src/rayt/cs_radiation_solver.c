/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_radiation_solver.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_convection_diffusion.c

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Descend binary tree for the lexicographical ordering of axis coordinates.
 *
 * parameters:
 *   s     <-- axial coordinate
 *   level <-- level of the binary tree to descend
 *   n     <-- number of entities in the binary tree to descend
 *   order <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_axis_descend_tree(const cs_real_t   s[],
                         size_t            level,
                         const size_t      n,
                         cs_lnum_t         order[])
{
  const cs_real_t eps = 1e-24;

  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (s[i1] > s[i2])
        lv_cur++;
      else if (s[i1] + eps > s[i2] && i1 > i2)
        lv_cur++;

    }

    if (lv_cur >= n) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (s[i1] > s[i2])
      break;
    if (s[i1] + eps >= s[i2] && i1 >= i2)
      break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order cells based on an axis coordinate and original numbering.
 *
 * parameters:
 *   s     <-- axial coordinate
 *   order --> pointer to pre-allocated ordering table
 *   n     <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_axis(const cs_real_t  s[],
            cs_lnum_t        order[],
            const size_t     n)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0; i < n; i++)
    order[i] = i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2);
  do {
    i--;
    _order_axis_descend_tree(s, i, n, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_axis_descend_tree(s, 0, i, order);
  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Order linear solvers for DOM radiative model.
 *----------------------------------------------------------------------------*/

void CS_PROCF (rayord, RAYORD)
(
 const cs_int_t          *ndirs,
 const cs_real_t          sx[],
 const cs_real_t          sy[],
 const cs_real_t          sz[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  int kdir = 0;

  cs_real_t *s;
  BFT_MALLOC(s, n_cells, cs_real_t);

  for (int ii = -1; ii < 2; ii+=2) {
    for (int jj = -1; jj < 2; jj+=2) {
      for (int kk = -1; kk < 2; kk+=2) {
        for (int idir = 0; idir < *ndirs; idir++) {

          cs_real_t v[3] = {ii*sx[idir], jj*sy[idir], kk*sz[idir]};

          kdir = kdir + 1;

          char name[32];
          sprintf(name, "radiation_%03d", kdir);

          cs_sles_t *sles = cs_sles_find(-1, name);

          if (sles == NULL) {
            (void)cs_sles_it_define(-1,
                                    name,
                                    CS_SLES_B_GAUSS_SEIDEL,
                                    0,      /* poly_degree */
                                    1000);  /* n_max_iter */
            sles = cs_sles_find(-1, name);
          }

          if (strcmp(cs_sles_get_type(sles), "cs_sles_it_t") != 0)
            continue;

          cs_sles_it_t *sc = cs_sles_get_context(sles);

          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            s[c_id] =   v[0]*cell_cen[c_id][0]
                      + v[1]*cell_cen[c_id][1]
                      + v[2]*cell_cen[c_id][2];

          cs_lnum_t *order;
          BFT_MALLOC(order, n_cells, cs_lnum_t);

          _order_axis(s, order, n_cells);

          cs_sles_it_assign_order(sc, &order); /* becomes owner of order */

        }
      }
    }
  }

  BFT_FREE(s);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
