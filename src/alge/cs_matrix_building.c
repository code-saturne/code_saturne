/*============================================================================
 * Matrix building
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Wrapper to cs_matrix_scalar (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrix, MATRIX)
(
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   ndircp,
 const cs_int_t  *const   isym,
 const cs_real_t *const   thetap,
 const cs_int_t  *const   imucpp,
 const cs_real_t          coefbp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          rovsdt[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          xcpp[],
 cs_real_t                da[],
 cs_real_t                xa[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  if (*isym != 1 && *isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (*isym == 1) {
    cs_sym_matrix_scalar(m,
                         *idiffp,
                         *ndircp,
                         *thetap,
                         cofbfp,
                         rovsdt,
                         i_visc,
                         b_visc,
                         da,
                         xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_scalar(m,
                     *iconvp,
                     *idiffp,
                     *ndircp,
                     *thetap,
                     *imucpp,
                     coefbp,
                     cofbfp,
                     rovsdt,
                     i_massflux,
                     b_massflux,
                     i_visc,
                     b_visc,
                     xcpp,
                     da,
                     (cs_real_2_t*) xa);
  }

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_vector (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrxv, MATRXV)
(
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   ndircp,
 const cs_int_t  *const   isym,
 const cs_real_t *const   thetap,
 const cs_real_33_t       coefbu[],
 const cs_real_33_t       cofbfu[],
 const cs_real_33_t       fimp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_33_t             da[],
 cs_real_t                xa[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  if (*isym != 1 && *isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (*isym == 1) {
    cs_sym_matrix_vector(m,
                         *idiffp,
                         *ndircp,
                         *thetap,
                         cofbfu,
                         fimp,
                         i_visc,
                         b_visc,
                         da,
                         xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_vector(m,
                     *iconvp,
                     *idiffp,
                     *ndircp,
                     *thetap,
                     coefbu,
                     cofbfu,
                     fimp,
                     i_massflux,
                     b_massflux,
                     i_visc,
                     b_visc,
                     da,
                     (cs_real_2_t*) xa);
  }
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_time_step
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
 cs_real_t                da[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  cs_matrix_time_step(m,
                      *iconvp,
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
 * Wrapper to cs_matrix_anisotropic_diffusion (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrvv, MATRVV)
(
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   ndircp,
 const cs_int_t  *const   isym,
 const cs_real_t *const   thetap,
 const cs_real_33_t       coefbu[],
 const cs_real_33_t       cofbfu[],
 const cs_real_33_t       fimp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_33_t       i_visc[],
 const cs_real_t          b_visc[],
 cs_real_33_t             da[],
 cs_real_332_t            xa[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  if (*isym != 1 && *isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (*isym == 1) {
    cs_sym_matrix_anisotropic_diffusion(m,
                                        *idiffp,
                                        *ndircp,
                                        *thetap,
                                        cofbfu,
                                        fimp,
                                        i_visc,
                                        b_visc,
                                        da,
                                        (cs_real_33_t*) xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_anisotropic_diffusion(m,
                                    *iconvp,
                                    *idiffp,
                                    *ndircp,
                                    *thetap,
                                    coefbu,
                                    cofbfu,
                                    fimp,
                                    i_massflux,
                                    b_massflux,
                                    i_visc,
                                    b_visc,
                                    da,
                                    (cs_real_332_t*) xa);
  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a scalar field.
 * (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (Implicit part)
 * \param[in]     rovsdt        working array
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_sym_matrix_scalar(const cs_mesh_t          *m,
         int                       idiffp,
         int                       ndircp,
         double                    thetap,
         const cs_real_t           cofbfp[],
         const cs_real_t           rovsdt[],
         const cs_real_t           i_visc[],
         const cs_real_t           b_visc[],
         cs_real_t       *restrict da,
         cs_real_t       *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
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

  const double epsi = 1.e-7;

  /*===============================================================================*/

  /*===============================================================================
    1. Initialization
    ===============================================================================*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = rovsdt[cell_id];
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }

# pragma omp parallel for
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    xa[face_id] = 0.;
  }

  /* 2. Computation of extradiagonal terms */

# pragma omp parallel for firstprivate(thetap, iconvp, idiffp)
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    xa[face_id] = -thetap*idiffp*i_visc[face_id];

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
     face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
     face_id++) {

  cs_lnum_t ii = i_face_cells[face_id][0] - 1;
  cs_lnum_t jj = i_face_cells[face_id][1] - 1;

  da[ii] -= xa[face_id];
  da[jj] -= xa[face_id];

      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                        if(n_b_faces > THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
     face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
     face_id++) {

  cs_lnum_t ii = b_face_cells[face_id] - 1;

  da[ii] += thetap*idiffp*b_visc[face_id]*cofbfp[face_id];

      }
    }
  }

  /* 5. If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      da[cell_id] = (1.+epsi)*da[cell_id];
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a scalar field
 * (non-symmetric matrix).
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (number of cells)
 * and an extra diagonal part (of dimension 2 time the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convective term by Cp
 *                               - 1 do multiply the convective term by Cp
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (Implicit part)
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
/*----------------------------------------------------------------------------*/

void
cs_matrix_scalar(const cs_mesh_t          *m,
                 int                       iconvp,
                 int                       idiffp,
                 int                       ndircp,
                 double                    thetap,
                 int                       imucpp,
                 const cs_real_t           coefbp[],
                 const cs_real_t           cofbfp[],
                 const cs_real_t           rovsdt[],
                 const cs_real_t           i_massflux[],
                 const cs_real_t           b_massflux[],
                 const cs_real_t           i_visc[],
                 const cs_real_t           b_visc[],
                 const cs_real_t           xcpp[],
                 cs_real_t       *restrict da,
                 cs_real_2_t     *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
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

  const double epsi = 1.e-7;

  /*===============================================================================*/

  /*===============================================================================
    1. Initialization
    ===============================================================================*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = rovsdt[cell_id];
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }


# pragma omp parallel for
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    xa[face_id][0] = 0.;
    xa[face_id][1] = 0.;
  }

  /* When solving the temperature, the convective part is multiplied by Cp */
  if (imucpp == 0) {

    /* 2. Computation of extradiagonal terms */

#   pragma omp parallel for firstprivate(thetap, iconvp, idiffp)
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      double flui = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));
      double fluj =-0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));

      xa[face_id][0] = thetap*(iconvp*flui -idiffp*i_visc[face_id]);
      xa[face_id][1] = thetap*(iconvp*fluj -idiffp*i_visc[face_id]);

    }

    /* 3. Contribution of the extra-diagonal terms to the diagonal */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
  for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
       face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
       face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0] - 1;
    cs_lnum_t jj = i_face_cells[face_id][1] - 1;

    da[ii] -= xa[face_id][0];
    da[jj] -= xa[face_id][1];

  }
      }
    }

    /* 4. Contribution of border faces to the diagonal */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                          if(n_b_faces > THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id] - 1;

          double flui = 0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));

          da[ii] += thetap*(  iconvp*flui*(coefbp[face_id]-1.)
                            + idiffp*b_visc[face_id]*cofbfp[face_id]);

        }
      }
    }

    /* When solving the temperature, the convective part is multiplied by Cp */
  } else {

    /* 2. Computation of extradiagonal terms */

#   pragma omp parallel for firstprivate(thetap, iconvp, idiffp)
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      double flui = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));
      double fluj =-0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));

      xa[face_id][0] = thetap*( iconvp*xcpp[i_face_cells[face_id][0] - 1]*flui
                               -idiffp*i_visc[face_id]);
      xa[face_id][1] = thetap*( iconvp*xcpp[i_face_cells[face_id][1] - 1]*fluj
                               -idiffp*i_visc[face_id]);

    }

    /* 3. Contribution of the extra-diagonal terms to the diagonal */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
  for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
       face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
       face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0] - 1;
    cs_lnum_t jj = i_face_cells[face_id][1] - 1;

    da[ii] -= xa[face_id][0];
    da[jj] -= xa[face_id][1];

  }
      }
    }

    /* 4. Contribution of boundary faces to the diagonal */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                 if(n_b_faces > THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id] - 1;
          double flui = 0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));
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
#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      da[cell_id] = (1.+epsi)*da[cell_id];
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a vector field
 * (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (Implicit part - 3x3 tensor array)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_sym_matrix_vector(const cs_mesh_t          *m,
                     int                       idiffp,
                     int                       ndircp,
                     double                    thetap,
                     const cs_real_33_t        cofbfu[],
                     const cs_real_33_t        fimp[],
                     const cs_real_t           i_visc[],
                     const cs_real_t           b_visc[],
                     cs_real_33_t    *restrict da,
                     cs_real_t       *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const double epsi = 1.e-7;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }

  }

  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {

      for (int isou = 0; isou < 3; isou++) {
        for (int jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }

    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    xa[face_id] = 0.;
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t face_id = 0; face_id <n_i_faces; face_id++) {

    xa[face_id] = -thetap*idiffp*i_visc[face_id];

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_i_faces; face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0] - 1;
    cs_lnum_t jj = i_face_cells[face_id][1] - 1;

    for (int isou = 0; isou < 3; isou++) {
      da[ii][isou][isou] -= xa[face_id];
      da[jj][isou][isou] -= xa[face_id];
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id] - 1;

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        if(isou == jsou) {
          da[ii][jsou][isou] += thetap*idiffp*b_visc[face_id]
                                      *cofbfu[face_id][jsou][isou];
        } else {
          da[ii][jsou][isou] += thetap*idiffp*b_visc[face_id]
                                      *cofbfu[face_id][jsou][isou];
        }
      }
    }

  }


  /* 5. If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

      for (int isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a vector field
 * (non-symmetric matrix).
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 2 time the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (Implicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (Implicit part - 3x3 tensor array)
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
cs_matrix_vector(const cs_mesh_t          *m,
                 int                       iconvp,
                 int                       idiffp,
                 int                       ndircp,
                 double                    thetap,
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
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const double epsi = 1.e-7;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }

  }

  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {

      for (int isou = 0; isou < 3; isou++) {
        for (int jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }

    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    xa[face_id][0] = 0.;
    xa[face_id][1] = 0.;
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t face_id = 0; face_id <n_i_faces; face_id++) {

    double flui = 0.5*( i_massflux[face_id] -fabs(i_massflux[face_id]) );
    double fluj =-0.5*( i_massflux[face_id] +fabs(i_massflux[face_id]) );

    xa[face_id][0] = thetap*(iconvp*flui -idiffp*i_visc[face_id]);
    xa[face_id][1] = thetap*(iconvp*fluj -idiffp*i_visc[face_id]);

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_i_faces; face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0] - 1;
    cs_lnum_t jj = i_face_cells[face_id][1] - 1;

    for (int isou = 0; isou < 3; isou++) {
      da[ii][isou][isou] -= xa[face_id][0];
      da[jj][isou][isou] -= xa[face_id][1];
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id] - 1;
    double flui = 0.5*( b_massflux[face_id] -fabs(b_massflux[face_id]) );

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
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
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

      for (int isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \param[in]     m             pointer to mesh structure
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
cs_matrix_time_step(const cs_mesh_t          *m,
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

  /* 1. Initialization */

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = 0.;
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }

  /* 2. Computation of extradiagonal terms unnecessary */

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  if (isym == 2) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0] - 1;
          cs_lnum_t jj = i_face_cells[face_id][1] - 1;

          double fluj =-0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
          double flui = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

          double xaifa2 = iconvp*fluj -idiffp*i_visc[face_id];
          double xaifa1 = iconvp*flui -idiffp*i_visc[face_id];
          da[ii] -= xaifa2;
          da[jj] -= xaifa1;

        }
      }
    }

  } else {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0] - 1;
          cs_lnum_t jj = i_face_cells[face_id][1] - 1;

          double flui = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

          double xaifa1 = iconvp*flui -idiffp*i_visc[face_id];
          da[ii] -= xaifa1;
          da[jj] -= xaifa1;

        }
      }
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for if(n_b_faces > THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id] - 1;

        double flui =  0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));
        double fluj = -0.5*(b_massflux[face_id] + fabs(b_massflux[face_id]));

        da[ii] +=   iconvp*(-fluj + flui*coefbp[face_id])
                  + idiffp*b_visc[face_id]*cofbfp[face_id];

      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a vector field with a
 * tensorial diffusivity.
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 2 times 3x3 the number of internal
 * faces).
 *
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (Implicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (Implicit part - 3x3 tensor array)
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
cs_matrix_anisotropic_diffusion(const cs_mesh_t          *m,
                                int                       iconvp,
                                int                       idiffp,
                                int                       ndircp,
                                double                    thetap,
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
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const double epsi = 1.e-7;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        for (int jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {
    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        xa[ifac][0][jsou][isou] = 0.;
        xa[ifac][1][jsou][isou] = 0.;
      }
    }
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {

    double flui = 0.5*( i_massflux[ifac] -fabs(i_massflux[ifac]) );
    double fluj =-0.5*( i_massflux[ifac] +fabs(i_massflux[ifac]) );

    for (int isou = 0; isou < 3; isou++) {
      xa[ifac][0][isou][isou] = iconvp*flui;
      xa[ifac][1][isou][isou] = iconvp*fluj;
      for (int jsou = 0; jsou < 3; jsou++) {
        xa[ifac][0][jsou][isou] = thetap*( xa[ifac][0][jsou][isou]
                                         - idiffp*i_visc[ifac][jsou][isou]);
        xa[ifac][1][jsou][isou] = thetap*( xa[ifac][1][jsou][isou]
                                         - idiffp*i_visc[ifac][jsou][isou]);
      }
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {

    cs_lnum_t ii = i_face_cells[ifac][0] - 1;
    cs_lnum_t jj = i_face_cells[ifac][1] - 1;

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[ii][jsou][isou] -= xa[ifac][0][jsou][isou];
        da[jj][jsou][isou] -= xa[ifac][1][jsou][isou];
      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

    cs_lnum_t ii = b_face_cells[ifac] - 1;
    double flui = 0.5*( b_massflux[ifac] -fabs(b_massflux[ifac]) );

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
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
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a vector field with a
 * tensorial diffusivity (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 3x3 the number of internal
 * faces).
 *
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     ndircp        indicator
 *                               - 0 if the diagonal stepped aside
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (Implicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (Implicit part - 3x3 tensor array)
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
cs_sym_matrix_anisotropic_diffusion(const cs_mesh_t           *m,
                                    int                       idiffp,
                                    int                       ndircp,
                                    double                    thetap,
                                    const cs_real_33_t        cofbfu[],
                                    const cs_real_33_t        fimp[],
                                    const cs_real_33_t        i_visc[],
                                    const cs_real_t           b_visc[],
                                    cs_real_33_t    *restrict da,
                                    cs_real_33_t    *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const double epsi = 1.e-7;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        for (int jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {
    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        xa[ifac][jsou][isou] = 0.;
      }
    }
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        xa[ifac][jsou][isou] = -thetap*idiffp*i_visc[ifac][jsou][isou];
      }
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {

    cs_lnum_t ii = i_face_cells[ifac][0] - 1;
    cs_lnum_t jj = i_face_cells[ifac][1] - 1;

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[ii][jsou][isou] -= xa[ifac][jsou][isou];
        da[jj][jsou][isou] -= xa[ifac][jsou][isou];
      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

    cs_lnum_t ii = b_face_cells[ifac] - 1;

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        if(isou == jsou) {
          da[ii][jsou][isou] += thetap*idiffp*b_visc[ifac]
                                      *cofbfu[ifac][jsou][isou];
        } else {
          da[ii][jsou][isou] += thetap*idiffp*b_visc[ifac]
                                      *cofbfu[ifac][jsou][isou];
        }
      }
    }

  }


  /* 5. If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum. */

  if ( ndircp <= 0 ) {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
