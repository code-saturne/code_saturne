/*============================================================================
 * Explicit matrix building
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_building.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_matrix_building.c

*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrix, MATRIX)
(
 const cs_int_t  *const   n_cells_ext,
 const cs_int_t  *const   n_cells,
 const cs_int_t  *const   n_i_faces,
 const cs_int_t  *const   n_b_faces,
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   ndircp,
 const cs_int_t  *const   isym,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   thetap,
 const cs_int_t  *const   imucpp,
 const cs_lnum_2_t        i_face_cells[],
 const cs_int_t           b_face_cells[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          rovsdt[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          xcpp[],
 cs_real_t                da[],
 cs_real_t                xa[][*n_i_faces])
{
  cs_matrix_scalar(*n_cells_ext,
                   *n_cells,
                   *n_i_faces,
                   *n_b_faces,
                   *iconvp,
                   *idiffp,
                   *ndircp,
                   *isym,
                   *thetap,
                   *imucpp,
                   i_face_cells,
                   b_face_cells,
                   coefbp,
                   cofbfp,
                   rovsdt,
                   i_massflux,
                   b_massflux,
                   i_visc,
                   b_visc,
                   xcpp,
                   da,
                   xa);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrxv
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrxv, MATRXV)
(
 const cs_int_t  *const   n_cells_ext,
 const cs_int_t  *const   n_cells,
 const cs_int_t  *const   n_i_faces,
 const cs_int_t  *const   n_b_faces,
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   ndircp,
 const cs_int_t  *const   isym,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   thetap,
 const cs_lnum_2_t        i_face_cells[],
 const cs_int_t           b_face_cells[],
 const cs_real_33_t       coefbu[],
 const cs_real_33_t       cofbfu[],
 const cs_real_33_t       fimp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_33_t             da[],
 cs_real_2_t              xa[])
{
  cs_matrix_vector(*n_cells_ext,
                   *n_cells,
                   *n_i_faces,
                   *n_b_faces,
                   *iconvp,
                   *idiffp,
                   *ndircp,
                   *isym,
                   *thetap,
                   i_face_cells,
                   b_face_cells,
                   coefbu,
                   cofbfu,
                   fimp,
                   i_massflux,
                   b_massflux,
                   i_visc,
                   b_visc,
                   da,
                   xa);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrdr
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrdt, MATRDT)
(
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   isym,
 const cs_real_t          coefbp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_t                da[])
{
  cs_matrix_time_step(*iconvp,
                      *idiffp,
                      *isym,
                      coefbp,
                      cofbfp,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      da);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrvv
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrvv, MATRVV)
(
 const cs_int_t  *const   n_cells_ext,
 const cs_int_t  *const   n_cells,
 const cs_int_t  *const   n_i_faces,
 const cs_int_t  *const   n_b_faces,
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   ndircp,
 const cs_int_t  *const   isym,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   thetap,
 const cs_lnum_2_t        i_face_cells[],
 const cs_int_t           b_face_cells[],
 const cs_real_33_t       coefbu[],
 const cs_real_33_t       cofbfu[],
 const cs_real_33_t       fimp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_33_t       i_visc[],
 const cs_real_t          b_visc[],
 cs_real_33_t             da[],
 cs_real_332_t            xa[])
{
  cs_matrix_tensorial_diffusion(*n_cells_ext,
                                *n_cells,
                                *n_i_faces,
                                *n_b_faces,
                                *iconvp,
                                *idiffp,
                                *ndircp,
                                *isym,
                                *thetap,
                                i_face_cells,
                                b_face_cells,
                                coefbu,
                                cofbfu,
                                fimp,
                                i_massflux,
                                b_massflux,
                                i_visc,
                                b_visc,
                                da,
                                xa);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*! \brief This function builds the matrix of advection/diffusion for a scalar
  field.

  The advection is upwind, the diffusion is not reconstructed.
  The matrix is splitted into a diagonal block (number of cells)
  and an extra diagonal part (of dimension 2 time the number of internal
  faces).

*/
/*-------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 ______________________________________________________________________________*/
/*!
 * \param[in]     n_cells_ext   number of extended (real + ghost) cells
 * \param[in]     n_cells       number of cells
 * \param[in]     n_i_faces     number of interior faces
 * \param[in]     n_b_faces     number of boundary faces
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     isym          indicator
 *                               - 1 symmetric matrix
 *                               - 2 non symmmetric matrix
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     i_face_cells  cell indices of interior faces
 * \param[in]     b_face_cells  cell indices of boundary faces
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (Impplicit part)
 * \param[in]     rovsdt        working array
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*-------------------------------------------------------------------------------*/

void
cs_matrix_scalar(
                 int                       n_cells_ext,
                 int                       n_cells,
                 int                       n_i_faces,
                 int                       n_b_faces,
                 int                       iconvp,
                 int                       idiffp,
                 int                       ndircp,
                 int                       isym,
                 double                    thetap,
                 int                       imucpp,
                 const cs_lnum_2_t         i_face_cells[],
                 const cs_int_t            b_face_cells[],
                 const cs_real_t           coefbp[],
                 const cs_real_t           cofbfp[],
                 const cs_real_t           rovsdt[],
                 const cs_real_t           i_massflux[],
                 const cs_real_t           b_massflux[],
                 const cs_real_t           i_visc[],
                 const cs_real_t           b_visc[],
                 const cs_real_t           xcpp[],
                 cs_real_t       *restrict da,
                 cs_real_t                 xa[][n_i_faces])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  /* Local variables */

  int face_id, ii, jj, cell_id, g_id, t_id;
  double flui, fluj, epsi;

  /*===============================================================================*/

  /*===============================================================================
    1. Initialization
    ===============================================================================*/

  if (isym != 1 && isym != 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of isym"));
  }

  epsi = 1.e-7;

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = rovsdt[cell_id];
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > THR_MIN)
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }

  if (isym == 2) {
#   pragma omp parallel for
    for (face_id = 0; face_id < n_i_faces; face_id++) {
      xa[0][face_id] = 0.;
      xa[1][face_id] = 0.;
    }
  } else {
#   pragma omp parallel for
    for (face_id = 0; face_id < n_i_faces; face_id++) {
      xa[0][face_id] = 0.;
    }
  }

  /* When solving the temperature, the convective part is multiplied by Cp */
  if (imucpp == 0) {

    /* 2. Computation of extradiagonal terms */

    if (isym == 2) {

#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                              private(flui, fluj)
      for (face_id = 0; face_id < n_i_faces; face_id++) {

        flui = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));
        fluj =-0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));

        xa[0][face_id] = thetap*(iconvp*flui -idiffp*i_visc[face_id]);
        xa[1][face_id] = thetap*(iconvp*fluj -idiffp*i_visc[face_id]);

      }

    } else {

#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                              private(flui)
      for (face_id = 0; face_id < n_i_faces; face_id++) {

        flui = 0.5*( i_massflux[face_id] -fabs(i_massflux[face_id]) );

        xa[0][face_id] = thetap*(iconvp*flui -idiffp*i_visc[face_id]);

      }

    }

    /* 3. Contribution of the extra-diagonal terms to the diagonal */

    if (isym == 2) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#      pragma omp parallel for private(face_id, ii, jj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0] - 1;
            jj = i_face_cells[face_id][1] - 1;

            da[ii] -= xa[0][face_id];
            da[jj] -= xa[1][face_id];

          }
        }
      }

    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0] - 1;
            jj = i_face_cells[face_id][1] - 1;

            da[ii] -= xa[0][face_id];
            da[jj] -= xa[0][face_id];

          }
        }
      }

    }

    /* 4. Contribution of border faces to the diagonal */

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                 private(face_id, ii, flui) if(n_b_faces > THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id] - 1;

          flui = 0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));

          da[ii] += thetap*( iconvp*flui*(coefbp[face_id]-1.)
                            + idiffp*b_visc[face_id]*cofbfp[face_id]);

        }
      }
    }

    /* When solving the temperature, the convective part is multiplied by Cp */
  } else {

    /* 2. Computation of extradiagonal terms */

    if (isym == 2) {

#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                              private(flui, fluj)
      for (face_id = 0; face_id < n_i_faces; face_id++) {

        flui = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));
        fluj =-0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));

        xa[0][face_id] = thetap*( iconvp*xcpp[i_face_cells[face_id][0] - 1]*flui
                                 -idiffp*i_visc[face_id]);
        xa[1][face_id] = thetap*( iconvp*xcpp[i_face_cells[face_id][1] - 1]*fluj
                                 -idiffp*i_visc[face_id]);

      }

    } else {

#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) private(flui)
      for (face_id = 0; face_id < n_i_faces; face_id++) {

        flui = 0.5*( i_massflux[face_id] -fabs(i_massflux[face_id]) );

        xa[0][face_id] = thetap*( iconvp*xcpp[i_face_cells[face_id][0] - 1]*flui
                                 -idiffp*i_visc[face_id]);

      }

    }

    /* 3. Contribution of the extra-diagonal terms to the diagonal */

    if (isym == 2) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0] - 1;
            jj = i_face_cells[face_id][1] - 1;

            da[ii] -= xa[0][face_id];
            da[jj] -= xa[1][face_id];

          }
        }
      }

    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0] - 1;
            jj = i_face_cells[face_id][1] - 1;

            da[ii] -= xa[0][face_id];
            da[jj] -= xa[0][face_id];
          }
        }
      }

    }

    /* 4. Contribution of boundary faces to the diagonal */

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                              private(face_id, ii, flui) \
                 if(n_b_faces > THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id] - 1;
          flui = 0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));
          da[ii] += thetap*(  iconvp*flui*xcpp[ii]*(coefbp[face_id]-1.)
                            + idiffp*b_visc[face_id]*cofbfp[face_id]);
        }
      }
    }

  }

  /* 5. If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
#   pragma omp parallel for firstprivate(epsi)
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      da[cell_id] = (1.+epsi)*da[cell_id];
    }
  }

}

/*----------------------------------------------------------------------------*/

/*! \brief This function builds the matrix of advection/diffusion for a vector
  field.

  The advection is upwind, the diffusion is not reconstructed.
  The matrix is splitted into a diagonal block (3x3 times number of cells)
  and an extra diagonal part (of dimension 2 time the number of internal
  faces).

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     n_cells_ext   number of extended (real + ghost) cells
 * \param[in]     n_cells       number of cells
 * \param[in]     n_i_faces     number of interior faces
 * \param[in]     n_b_faces     number of boundary faces
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     isym          indicator
 *                               - 1 symmetric matrix
 *                               - 2 non symmmetric matrix
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     i_face_cells  cell indexes of interior faces
 * \param[in]     b_face_cells  cell indexes of boundary faces
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (Impplicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (Impplicit part - 3x3 tensor array)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector(
                 int                       n_cells_ext,
                 int                       n_cells,
                 int                       n_i_faces,
                 int                       n_b_faces,
                 int                       iconvp,
                 int                       idiffp,
                 int                       ndircp,
                 int                       isym,
                 double                    thetap,
                 const cs_lnum_2_t         i_face_cells[],
                 const cs_int_t            b_face_cells[],
                 const cs_real_33_t        coefbu[],
                 const cs_real_33_t        cofbfu[],
                 const cs_real_33_t        fimp[],
                 const cs_real_t           i_massflux[],
                 const cs_real_t           b_massflux[],
                 const cs_real_t           i_visc[],
                 const cs_real_t           b_visc[],
                 cs_real_33_t    *restrict da,
                 cs_real_2_t     *restrict xa)
{

  /* Local variables */

  int face_id,ii,jj,cell_id, isou, jsou;
  double flui,fluj,epsi;

  /* 1. Initialization */

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  epsi = 1.e-7;

  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    for (isou = 0; isou < 3; isou++) {
      for (jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if(n_cells_ext > n_cells) {
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  if(isym == 2) {
    for (face_id = 0; face_id < n_i_faces; face_id++) {
      xa[face_id][0] = 0.;
      xa[face_id][1] = 0.;
    }
  } else {
    for (face_id = 0; face_id < n_i_faces; face_id++) {
      xa[face_id][0] = 0.;
    }
  }

  /* 2. Computation of extradiagonal terms */

  if(isym == 2) {

    for (face_id = 0; face_id <n_i_faces; face_id++) {
      flui = 0.5*( i_massflux[face_id] -fabs(i_massflux[face_id]) );
      fluj =-0.5*( i_massflux[face_id] +fabs(i_massflux[face_id]) );
      xa[face_id][0] = thetap*(iconvp*flui -idiffp*i_visc[face_id]);
      xa[face_id][1] = thetap*(iconvp*fluj -idiffp*i_visc[face_id]);
    }

  } else {

    for (face_id = 0; face_id <n_i_faces; face_id++) {
      flui = 0.5*( i_massflux[face_id] -fabs(i_massflux[face_id]) );
      xa[face_id][0] = thetap*(iconvp*flui -idiffp*i_visc[face_id]);
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  if(isym == 2) {

    for (face_id = 0; face_id <n_i_faces; face_id++) {
      ii = i_face_cells[face_id][0] - 1;
      jj = i_face_cells[face_id][1] - 1;
      for (isou = 0; isou < 3; isou++) {
        da[ii][isou][isou] -= xa[face_id][0];
        da[jj][isou][isou] -= xa[face_id][1];
      }
    }


  } else {

    for (face_id = 0; face_id <n_i_faces; face_id++) {
      ii = i_face_cells[face_id][0] - 1;
      jj = i_face_cells[face_id][1] - 1;
      for (isou = 0; isou < 3; isou++) {
        da[ii][isou][isou] -= xa[face_id][0];
        da[jj][isou][isou] -= xa[face_id][0];
      }
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (face_id = 0; face_id <n_b_faces; face_id++) {
    ii = b_face_cells[face_id] - 1;
    flui = 0.5*( b_massflux[face_id] -fabs(b_massflux[face_id]) );
    for (isou = 0; isou < 3; isou++) {
      for (jsou = 0; jsou < 3; jsou++) {
        if(isou == jsou) {
          da[ii][jsou][isou] += thetap*( iconvp*flui
                                         *(coefbu[face_id][jsou][isou]-1.)
                                        +idiffp*b_visc[face_id]
                                         *cofbfu[face_id][jsou][isou] );
        } else {
          da[ii][jsou][isou] += thetap*( iconvp*flui*coefbu[face_id][jsou][isou]
                                        +idiffp*b_visc[face_id]
                                         *cofbfu[face_id][jsou][isou] );
        }
      }
    }
  }


  /* 5. If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      for (isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

}

/*----------------------------------------------------------------------------*/

/*! \brief Construction of the diagonal of the advection/diffusion matrix
  for determining the variable time step, flow, Fourier.

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     isym          indicator
 *                               - 1 symmetric matrix
 *                               - 2 non symmmetric matrix
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_time_step(
                    int                       iconvp,
                    int                       idiffp,
                    int                       isym,
                    const cs_real_t           coefbp[],
                    const cs_real_t           cofbfp[],
                    const cs_real_t           i_massflux[],
                    const cs_real_t           b_massflux[],
                    const cs_real_t           i_visc[],
                    const cs_real_t           b_visc[],
                    cs_real_t       *restrict da)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Local variables */

  int face_id, ii, jj, cell_id, g_id, t_id;
  double flui, fluj, xaifa1, xaifa2;

  /* 1. Initialization */

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = 0.;
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > THR_MIN)
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }

  /* 2. Computation of extradiagonal terms unnecessary */

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  if (isym == 2) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, fluj, flui,             \
                                      xaifa2, xaifa1)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0] - 1;
          jj = i_face_cells[face_id][1] - 1;

          fluj =-0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
          flui = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

          xaifa2 = iconvp*fluj -idiffp*i_visc[face_id];
          xaifa1 = iconvp*flui -idiffp*i_visc[face_id];
          da[ii] -= xaifa2;
          da[jj] -= xaifa1;

        }
      }
    }

  } else {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, flui, xaifa1)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0] - 1;
          jj = i_face_cells[face_id][1] - 1;

          flui = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

          xaifa1 = iconvp*flui -idiffp*i_visc[face_id];
          da[ii] -= xaifa1;
          da[jj] -= xaifa1;

        }
      }
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for private(face_id, ii, flui, fluj) \
               if(m->n_b_faces > THR_MIN)
    for (t_id = 0; t_id < n_b_threads; t_id++) {
      for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        ii = b_face_cells[face_id] - 1;

        flui =  0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));
        fluj = -0.5*(b_massflux[face_id] + fabs(b_massflux[face_id]));

        da[ii] +=   iconvp*(-fluj + flui*coefbp[face_id])
                  + idiffp*b_visc[face_id]*cofbfp[face_id];
      }
    }
  }

}

/*----------------------------------------------------------------------------*/

/*! \brief This function builds the matrix of advection/diffusion for a vector
  field with a tensorial diffusivity.

  The advection is upwind, the diffusion is not reconstructed.
  The matrix is splitted into a diagonal block (3x3 times number of cells)
  and an extra diagonal part (of dimension 2 times 3x3 the number of internal
  faces).

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     n_cells_ext   number of extended (real + ghost) cells
 * \param[in]     n_cells       number of cells
 * \param[in]     n_i_faces     number of interior faces
 * \param[in]     n_b_faces     number of boundary faces
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     isym          indicator
 *                               - 1 symmetric matrix
 *                               - 2 non symmmetric matrix
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     i_face_cells  cell indexes of interior faces
 * \param[in]     b_face_cells  cell indexes of boundary faces
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (Impplicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (Impplicit part - 3x3 tensor array)
 * \param[in]     fimp
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_tensorial_diffusion(
                              int                       n_cells_ext,
                              int                       n_cells,
                              int                       n_i_faces,
                              int                       n_b_faces,
                              int                       iconvp,
                              int                       idiffp,
                              int                       ndircp,
                              int                       isym,
                              double                    thetap,
                              const cs_lnum_2_t         i_face_cells[],
                              const cs_int_t            b_face_cells[],
                              const cs_real_33_t        coefbu[],
                              const cs_real_33_t        cofbfu[],
                              const cs_real_33_t        fimp[],
                              const cs_real_t           i_massflux[],
                              const cs_real_t           b_massflux[],
                              const cs_real_33_t        i_visc[],
                              const cs_real_t           b_visc[],
                              cs_real_33_t    *restrict da,
                              cs_real_332_t   *restrict xa)
{

  /* Local variables */

  int ifac,ii,jj,cell_id, isou, jsou;
  double flui,fluj,epsi;

  /* 1. Initialization */

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  epsi = 1.e-7;

  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    for (isou = 0; isou < 3; isou++) {
      for (jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if(n_cells_ext > n_cells) {
    for (cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  if(isym == 2) {
    for (ifac = 0; ifac < n_i_faces; ifac++) {
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++) {
          xa[ifac][0][jsou][isou] = 0.;
          xa[ifac][1][jsou][isou] = 0.;
        }
      }
    }
  } else {
    for (ifac = 0; ifac < n_i_faces; ifac++) {
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++) {
          xa[ifac][0][jsou][isou] = 0.;
        }
      }
    }
  }

  /* 2. Computation of extradiagonal terms */

  if(isym == 2) {

    for (ifac = 0; ifac < n_i_faces; ifac++) {
      flui = 0.5*( i_massflux[ifac] -fabs(i_massflux[ifac]) );
      fluj =-0.5*( i_massflux[ifac] +fabs(i_massflux[ifac]) );
      for (isou = 0; isou < 3; isou++) {
        xa[ifac][0][isou][isou] = iconvp*flui;
        xa[ifac][1][isou][isou] = iconvp*fluj;
        for (jsou = 0; jsou < 3; jsou++) {
          xa[ifac][0][jsou][isou] = thetap*( xa[ifac][0][jsou][isou]
                                           - idiffp*i_visc[ifac][jsou][isou]);
          xa[ifac][1][jsou][isou] = thetap*( xa[ifac][1][jsou][isou]
                                           - idiffp*i_visc[ifac][jsou][isou]);
        }
      }
    }

  } else {

    for (ifac = 0; ifac < n_i_faces; ifac++) {
      flui = 0.5*(i_massflux[ifac] -fabs(i_massflux[ifac]));
      for (isou = 0; isou < 3; isou++) {
        xa[ifac][0][isou][isou] = iconvp*flui;
        for (jsou = 0; jsou < 3; jsou++) {
          xa[ifac][0][jsou][isou] = thetap*( xa[ifac][0][jsou][isou]
                                           - idiffp*i_visc[ifac][jsou][isou]);
        }
      }
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  if (isym == 2) {

    for (ifac = 0; ifac < n_i_faces; ifac++) {
      ii = i_face_cells[ifac][0] - 1;
      jj = i_face_cells[ifac][1] - 1;
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++) {
          da[ii][jsou][isou] -= xa[ifac][0][jsou][isou];
          da[jj][jsou][isou] -= xa[ifac][1][jsou][isou];
        }
      }
    }

  } else {

    for (ifac = 0; ifac < n_i_faces; ifac++) {
      ii = i_face_cells[ifac][0] - 1;
      jj = i_face_cells[ifac][1] - 1;
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++) {
          da[ii][jsou][isou] -= xa[ifac][0][jsou][isou];
          da[jj][jsou][isou] -= xa[ifac][0][jsou][isou];
        }
      }
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (ifac = 0; ifac < n_b_faces; ifac++) {
    ii = b_face_cells[ifac] - 1;
    flui = 0.5*( b_massflux[ifac] -fabs(b_massflux[ifac]) );
    for (isou = 0; isou < 3; isou++) {
      for (jsou = 0; jsou < 3; jsou++) {
        if(isou == jsou) {
          da[ii][jsou][isou] += thetap*( iconvp*flui
                                         *(coefbu[ifac][jsou][isou]-1.)
                                        +idiffp*b_visc[ifac]
                                         *cofbfu[ifac][jsou][isou] );
        } else {
          da[ii][jsou][isou] += thetap*( iconvp*flui*coefbu[ifac][jsou][isou]
                                        +idiffp*b_visc[ifac]
                                         *cofbfu[ifac][jsou][isou] );
        }
      }
    }
  }


  /* 5. If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum. */

  if ( ndircp <= 0 ) {
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      for (isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
