/*============================================================================
 * Matrix building
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fortran wrapper to cs_matrix_scalar (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void CS_PROCF (matrix, MATRIX)
(
 const int        *iconvp,
 const int        *idiffp,
 const int        *ndircp,
 const int        *isym,
 const cs_real_t  *thetap,
 const int        *imucpp,
 const cs_real_t   coefbp[],
 const cs_real_t   cofbfp[],
 const cs_real_t   rovsdt[],
 const cs_real_t   i_massflux[],
 const cs_real_t   b_massflux[],
 const cs_real_t   i_visc[],
 const cs_real_t   b_visc[],
 const cs_real_t   xcpp[],
 cs_real_t         da[],
 cs_real_t         xa[]
)
{
  cs_matrix_wrapper_scalar(*iconvp,
                           *idiffp,
                           *ndircp,
                           *isym,
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
                           xa);
}

/*----------------------------------------------------------------------------
 * Fortran wrapper to cs_matrix_time_step
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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_scalar (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_scalar(int               iconvp,
                         int               idiffp,
                         int               ndircp,
                         int               isym,
                         double            thetap,
                         int               imucpp,
                         const cs_real_t   coefbp[],
                         const cs_real_t   cofbfp[],
                         const cs_real_t   rovsdt[],
                         const cs_real_t   i_massflux[],
                         const cs_real_t   b_massflux[],
                         const cs_real_t   i_visc[],
                         const cs_real_t   b_visc[],
                         const cs_real_t   xcpp[],
                         cs_real_t         da[],
                         cs_real_t         xa[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (isym == 1) {
    cs_sym_matrix_scalar(m,
                         idiffp,
                         thetap,
                         cofbfp,
                         rovsdt,
                         i_visc,
                         b_visc,
                         da,
                         xa);
  }

  /* Non-symmetric matrix */
  else {
    cs_matrix_scalar(m,
                     iconvp,
                     idiffp,
                     thetap,
                     imucpp,
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

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    const double epsi = 1.e-7;

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      da[cell_id] = (1.+epsi)*da[cell_id];
    }
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] += mq->c_solid_flag[CS_MIN(cs_glob_porous_model, 1)*cell_id];
  }

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_scalar for convection/diffusion multigrid
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_scalar_conv_diff(int               iconvp,
                                   int               idiffp,
                                   int               ndircp,
                                   double            thetap,
                                   int               imucpp,
                                   const cs_real_t   coefbp[],
                                   const cs_real_t   cofbfp[],
                                   const cs_real_t   rovsdt[],
                                   const cs_real_t   i_massflux[],
                                   const cs_real_t   b_massflux[],
                                   const cs_real_t   i_visc[],
                                   const cs_real_t   b_visc[],
                                   const cs_real_t   xcpp[],
                                   cs_real_t         da[],
                                   cs_real_t         xa[],
                                   cs_real_t         da_conv[],
                                   cs_real_t         xa_conv[],
                                   cs_real_t         da_diff[],
                                   cs_real_t         xa_diff[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  /* Diffusion matrix */
  cs_sym_matrix_scalar(m,
                       idiffp,
                       thetap,
                       cofbfp,
                       rovsdt,
                       i_visc,
                       b_visc,
                       da_diff,
                       xa_diff);

  /* Convection matrix */
  cs_matrix_scalar(m,
                   iconvp,
                   0.,
                   thetap,
                   imucpp,
                   coefbp,
                   cofbfp,
                   rovsdt,
                   i_massflux,
                   b_massflux,
                   i_visc,
                   b_visc,
                   xcpp,
                   da_conv,
                   (cs_real_2_t*)xa_conv);

  /* Global matrix */

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da_conv[cell_id] -= rovsdt[cell_id];
    da_diff[cell_id] -= rovsdt[cell_id];
    da[cell_id] = rovsdt[cell_id] + da_diff[cell_id] + da_conv[cell_id];
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }

# pragma omp parallel for
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    xa[2*face_id]    = xa_diff[face_id] + xa_conv[2*face_id];
    xa[2*face_id +1] = xa_diff[face_id] + xa_conv[2*face_id +1];
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    const double epsi = 1.e-7;

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      da[cell_id] = (1.+epsi)*da[cell_id];
    }
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] += mq->c_solid_flag[CS_MIN(cs_glob_porous_model, 1)*cell_id];
  }

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_vector (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_vector(int                  iconvp,
                         int                  idiffp,
                         int                  ndircp,
                         int                  isym,
                         double               thetap,
                         const cs_real_33_t   coefbu[],
                         const cs_real_33_t   cofbfu[],
                         const cs_real_33_t   fimp[],
                         const cs_real_t      i_massflux[],
                         const cs_real_t      b_massflux[],
                         const cs_real_t      i_visc[],
                         const cs_real_t      b_visc[],
                         cs_real_33_t         da[],
                         cs_real_t            xa[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (isym == 1) {
    cs_sym_matrix_vector(m,
                         idiffp,
                         thetap,
                         cofbfu,
                         fimp,
                         i_visc,
                         b_visc,
                         da,
                         xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_vector(m,
                     iconvp,
                     idiffp,
                     thetap,
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

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum. */

  if (ndircp <= 0) {
    const double epsi = 1.e-7;
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 3; isou++)
      da[cell_id][isou][isou]
        += mq->c_solid_flag[CS_MIN(cs_glob_porous_model, 1)*cell_id];
  }

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_tensor (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_tensor(int                  iconvp,
                         int                  idiffp,
                         int                  ndircp,
                         int                  isym,
                         double               thetap,
                         const cs_real_66_t   coefbts[],
                         const cs_real_66_t   cofbfts[],
                         const cs_real_66_t   fimp[],
                         const cs_real_t      i_massflux[],
                         const cs_real_t      b_massflux[],
                         const cs_real_t      i_visc[],
                         const cs_real_t      b_visc[],
                         cs_real_66_t         da[],
                         cs_real_t            xa[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (isym == 1) {
    cs_sym_matrix_tensor(m,
                         idiffp,
                         thetap,
                         cofbfts,
                         fimp,
                         i_visc,
                         b_visc,
                         da,
                         xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_tensor(m,
                     iconvp,
                     idiffp,
                     thetap,
                     coefbts,
                     cofbfts,
                     fimp,
                     i_massflux,
                     b_massflux,
                     i_visc,
                     b_visc,
                     da,
                     (cs_real_2_t*) xa);
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    const double epsi = 1.e-7;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 6; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_lnum_t corr_cell_id = CS_MIN(cs_glob_porous_model, 1)*cell_id;
    for (int isou = 0; isou < 6; isou++) {
      da[cell_id][isou][isou] += mq->c_solid_flag[corr_cell_id];
    }
  }

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_anisotropic_diffusion (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_anisotropic_diffusion_wrapper(int                  iconvp,
                                        int                  idiffp,
                                        int                  ndircp,
                                        int                  isym,
                                        double               thetap,
                                        const cs_real_33_t   coefbu[],
                                        const cs_real_33_t   cofbfu[],
                                        const cs_real_33_t   fimp[],
                                        const cs_real_t      i_massflux[],
                                        const cs_real_t      b_massflux[],
                                        const cs_real_33_t   i_visc[],
                                        const cs_real_t      b_visc[],
                                        cs_real_33_t         da[],
                                        cs_real_t            xa[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;


  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (isym == 1) {
    cs_sym_matrix_anisotropic_diffusion(m,
                                        idiffp,
                                        thetap,
                                        cofbfu,
                                        fimp,
                                        i_visc,
                                        b_visc,
                                        da,
                                        (cs_real_33_t *) xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_anisotropic_diffusion(m,
                                    iconvp,
                                    idiffp,
                                    thetap,
                                    coefbu,
                                    cofbfu,
                                    fimp,
                                    i_massflux,
                                    b_massflux,
                                    i_visc,
                                    b_visc,
                                    da,
                                    (cs_real_332_t *) xa);
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum. */

  if (ndircp <= 0) {
    const double epsi = 1.e-7;
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 3; isou++)
      da[cell_id][isou][isou] += mq->c_solid_flag[CS_MIN(cs_glob_porous_model, 1)*cell_id];
  }

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_anisotropic_diffusion_tensor (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_anisotropic_diffusion_wrapper_tensor(int                  iconvp,
                                               int                  idiffp,
                                               int                  ndircp,
                                               int                  isym,
                                               double               thetap,
                                               const cs_real_66_t   coefbts[],
                                               const cs_real_66_t   cofbfts[],
                                               const cs_real_66_t   fimp[],
                                               const cs_real_t      i_massflux[],
                                               const cs_real_t      b_massflux[],
                                               const cs_real_66_t   i_visc[],
                                               const cs_real_t      b_visc[],
                                               cs_real_66_t         da[],
                                               cs_real_662_t        xa[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;

  if (isym != 1 && isym != 2) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isym"));
  }

  /* Symmetric matrix */
  if (isym == 1) {
    cs_sym_matrix_anisotropic_diffusion_tensor(m,
                                               idiffp,
                                               thetap,
                                               cofbfts,
                                               fimp,
                                               i_visc,
                                               b_visc,
                                               da,
                                               (cs_real_66_t *)xa);

  /* Non-symmetric matrix */
  } else {
    cs_matrix_anisotropic_diffusion_tensor(m,
                                           iconvp,
                                           idiffp,
                                           thetap,
                                           coefbts,
                                           cofbfts,
                                           fimp,
                                           i_massflux,
                                           b_massflux,
                                           i_visc,
                                           b_visc,
                                           da,
                                           xa);
  }

  /* Penalization if non invertible matrix */

  /* If no Dirichlet condition, the diagonal is slightly increased in order
     to shift the eigenvalues spectrum (if IDIRCL=0, we force NDIRCP to be at
     least 1 in order not to shift the diagonal). */

  if (ndircp <= 0) {
    const double epsi = 1.e-7;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 6; isou++) {
        da[cell_id][isou][isou] = (1.+epsi)*da[cell_id][isou][isou];
      }
    }
  }

  /* If a whole line of the matrix is 0, the diagonal is set to 1 */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_lnum_t corr_cell_id = CS_MIN(cs_glob_porous_model, 1)*cell_id;
    for (int isou = 0; isou < 6; isou++) {
      da[cell_id][isou][isou] += mq->c_solid_flag[corr_cell_id];
    }
  }

}

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
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (implicit part)
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

  /*===============================================================================*/

  /*===============================================================================
    1. Initialization
    ===============================================================================*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = rovsdt[cell_id];
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      da[cell_id] = 0.;
    }
  }

# pragma omp parallel for
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    xa[face_id] = 0.;
  }

  /* 2. Computation of extradiagonal terms */

# pragma omp parallel for firstprivate(thetap, idiffp)
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

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        da[ii] -= xa[face_id];
        da[jj] -= xa[face_id];

      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (int g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for firstprivate(thetap, idiffp) \
                        if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
        face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
        face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        da[ii] += thetap*idiffp*b_visc[face_id]*cofbfp[face_id];

      }
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
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convective term by Cp
 *                               - 1 do multiply the convective term by Cp
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (implicit part)
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

  /*===============================================================================*/

  /*===============================================================================
    1. Initialization
    ===============================================================================*/

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    da[cell_id] = rovsdt[cell_id];
  }
  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if (n_cells_ext - n_cells > CS_THR_MIN)
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

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          /* D_ii =  theta (m_ij)^+ - m_ij
           *      = -X_ij - (1-theta)*m_ij
           * D_jj = -theta (m_ij)^- + m_ij
           *      = -X_ji + (1-theta)*m_ij
           */
          da[ii] -= xa[face_id][0] + iconvp*(1. - thetap)*i_massflux[face_id];
          da[jj] -= xa[face_id][1] - iconvp*(1. - thetap)*i_massflux[face_id];

        }
      }
    }

    /* 4. Contribution of border faces to the diagonal */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                          if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          double flui = 0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));

          /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
           *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
           *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
           */
          da[ii] += iconvp*(flui*thetap*(coefbp[face_id]-1.)
                           -(1.-thetap)*b_massflux[face_id])
                  + idiffp*thetap*b_visc[face_id]*cofbfp[face_id];
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

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      xa[face_id][0] = thetap*( iconvp*xcpp[ii]*flui
                               -idiffp*i_visc[face_id]);
      xa[face_id][1] = thetap*( iconvp*xcpp[jj]*fluj
                               -idiffp*i_visc[face_id]);

    }

    /* 3. Contribution of the extra-diagonal terms to the diagonal */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
            face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
            face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          /* D_ii =  theta (m_ij)^+ - m_ij
           *      = -X_ij - (1-theta)*m_ij
           * D_jj = -theta (m_ij)^- + m_ij
           *      = -X_ji + (1-theta)*m_ij
           */
          da[ii] -= xa[face_id][0] + iconvp*(1. - thetap)*xcpp[ii]*i_massflux[face_id];
          da[jj] -= xa[face_id][1] - iconvp*(1. - thetap)*xcpp[jj]*i_massflux[face_id];

        }
      }
    }

    /* 4. Contribution of boundary faces to the diagonal */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for firstprivate(thetap, iconvp, idiffp) \
                 if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          double flui = 0.5*(b_massflux[face_id] - fabs(b_massflux[face_id]));

          /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
           *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
           *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
           */
          da[ii] += iconvp*xcpp[ii]*(flui*thetap*(coefbp[face_id]-1.)
                                    -(1.-thetap)*b_massflux[face_id])
                  + idiffp*thetap*b_visc[face_id]*cofbfp[face_id];
        }
      }
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
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
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

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    for (int isou = 0; isou < 3; isou++) {
      da[ii][isou][isou] -= xa[face_id];
      da[jj][isou][isou] -= xa[face_id];
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[ii][jsou][isou] += thetap*idiffp*b_visc[face_id]
                                    *cofbfu[face_id][jsou][isou];
      }
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diffusion matrix for a tensor field
 * (symmetric matrix).
 *
 * The diffusion is not reconstructed.
 * The matrix is split into a diagonal block (6x6 times number of cells)
 * and an extra diagonal part (of dimension the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfts        boundary condition array for the variable flux
 *                               (Implicit part - 6x6 tensor array)
 * \param[in]     fimp          part of the diagonal
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_sym_matrix_tensor(const cs_mesh_t          *m,
                     int                       idiffp,
                     double                    thetap,
                     const cs_real_66_t        cofbfts[],
                     const cs_real_66_t        fimp[],
                     const cs_real_t           i_visc[],
                     const cs_real_t           b_visc[],
                     cs_real_66_t              *restrict da,
                     cs_real_t                 *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }

  }

  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {

      for (int isou = 0; isou < 6; isou++) {
        for (int jsou = 0; jsou < 6; jsou++) {
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

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    for (int isou = 0; isou < 6; isou++) {
      da[ii][isou][isou] -= xa[face_id];
      da[jj][isou][isou] -= xa[face_id];
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[ii][jsou][isou] += thetap*idiffp*b_visc[face_id]
                                    *cofbfts[face_id][jsou][isou];
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
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     fimp          part of the diagonal
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

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }

  }

  if (n_cells_ext > n_cells) {
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

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    /* D_ii =  theta (m_ij)^+ - m_ij
     *      = -X_ij - (1-theta)*m_ij
     * D_jj = -theta (m_ij)^- + m_ij
     *      = -X_ji + (1-theta)*m_ij
     */
    for (int isou = 0; isou < 3; isou++) {
      da[ii][isou][isou] -= xa[face_id][0]
                          + iconvp*(1. - thetap)*i_massflux[face_id];
      da[jj][isou][isou] -= xa[face_id][1]
                          - iconvp*(1. - thetap)*i_massflux[face_id];
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];
    double flui = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
         *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
         *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
         */
        if(isou == jsou) {
          da[ii][jsou][isou] += iconvp*( thetap*flui
                                        *(coefbu[face_id][jsou][isou]-1.)
                                       - (1. - thetap)*b_massflux[face_id])
                              + idiffp*thetap*b_visc[face_id]
                                      *cofbfu[face_id][jsou][isou];
        } else {
          da[ii][jsou][isou] += thetap*( iconvp*flui*coefbu[face_id][jsou][isou]
                                        +idiffp*b_visc[face_id]
                                         *cofbfu[face_id][jsou][isou] );
        }
      }
    }

  }

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a tensor field
 * (non-symmetric matrix).
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (6x6 times number of cells)
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
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbts        boundary condition array for the variable
 *                               (Implicit part - 6x6 tensor array)
 * \param[in]     cofbfts        boundary condition array for the variable flux
 *                               (Implicit part - 6x6 tensor array)
 * \param[in]     fimp          part of the diagonal
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
cs_matrix_tensor(const cs_mesh_t          *m,
                 int                       iconvp,
                 int                       idiffp,
                 double                    thetap,
                 const cs_real_66_t        coefbts[],
                 const cs_real_66_t        cofbfts[],
                 const cs_real_66_t        fimp[],
                 const cs_real_t           i_massflux[],
                 const cs_real_t           b_massflux[],
                 const cs_real_t           i_visc[],
                 const cs_real_t           b_visc[],
                 cs_real_66_t              *restrict da,
                 cs_real_2_t               *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }

  }

  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {

      for (int isou = 0; isou < 6; isou++) {
        for (int jsou = 0; jsou < 6; jsou++) {
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

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    /* D_ii =  theta (m_ij)^+ - m_ij
     *      = -X_ij - (1-theta)*m_ij
     * D_jj = -theta (m_ij)^- + m_ij
     *      = -X_ji + (1-theta)*m_ij
     */
    for (int isou = 0; isou < 6; isou++) {
      da[ii][isou][isou] -= xa[face_id][0]
                          + iconvp*(1. - thetap)*i_massflux[face_id];
      da[jj][isou][isou] -= xa[face_id][1]
                          - iconvp*(1. - thetap)*i_massflux[face_id];
    }

  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id <n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];
    double flui = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
         *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
         *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
         */
        if(isou == jsou) {
          da[ii][jsou][isou] += iconvp*( thetap*flui
                                        *(coefbts[face_id][jsou][isou]-1.)
                                       - (1. - thetap)*b_massflux[face_id])
                              + idiffp*thetap*b_visc[face_id]
                                      *cofbfts[face_id][jsou][isou];
        } else {
          da[ii][jsou][isou] += thetap*( iconvp*flui*coefbts[face_id][jsou][isou]
                                        +idiffp*b_visc[face_id]
                                         *cofbfts[face_id][jsou][isou] );
        }
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
 *                               (implicit part)
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (implicit part)
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
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
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

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

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

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

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
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

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
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     cofbfp        boundary condition array for the variable flux
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     fimp          part of the diagonal
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
                                double                    thetap,
                                const cs_real_33_t        coefbp[],
                                const cs_real_33_t        cofbfp[],
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

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if (n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        for (int jsou = 0; jsou < 3; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        xa[face_id][0][jsou][isou] = 0.;
        xa[face_id][1][jsou][isou] = 0.;
      }
    }
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    double flui = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));
    double fluj =-0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));

    for (int isou = 0; isou < 3; isou++) {
      xa[face_id][0][isou][isou] = iconvp*flui;
      xa[face_id][1][isou][isou] = iconvp*fluj;
      for (int jsou = 0; jsou < 3; jsou++) {
        xa[face_id][0][jsou][isou] = thetap*( xa[face_id][0][jsou][isou]
                                         - idiffp*i_visc[face_id][jsou][isou]);
        xa[face_id][1][jsou][isou] = thetap*( xa[face_id][1][jsou][isou]
                                         - idiffp*i_visc[face_id][jsou][isou]);
      }
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    /* D_ii =  theta (m_ij)^+ - m_ij
     *      = -X_ij - (1-theta)*m_ij
     * D_jj = -theta (m_ij)^- + m_ij
     *      = -X_ji + (1-theta)*m_ij
     */
    for (int isou = 0; isou < 3; isou++) {
      da[ii][isou][isou] -= iconvp*(1. - thetap)*i_massflux[face_id];
      da[jj][isou][isou] += iconvp*(1. - thetap)*i_massflux[face_id];
      for (int jsou = 0; jsou < 3; jsou++) {
        da[ii][jsou][isou] -= xa[face_id][0][jsou][isou];
        da[jj][jsou][isou] -= xa[face_id][1][jsou][isou];
      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];
    double flui = 0.5*( b_massflux[face_id] -fabs(b_massflux[face_id]) );

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
         *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
         *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
         */
        if(isou == jsou) {
          da[ii][jsou][isou] += iconvp*( thetap*flui
                                        *(coefbp[face_id][jsou][isou]-1.)
                                       - (1. - thetap)*b_massflux[face_id])
                              + idiffp*thetap*b_visc[face_id]
                                      *cofbfp[face_id][jsou][isou];
        } else {
          da[ii][jsou][isou] += thetap*( iconvp*flui*coefbp[face_id][jsou][isou]
                                        +idiffp*b_visc[face_id]
                                         *cofbfp[face_id][jsou][isou] );
        }
      }
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the advection/diffusion matrix for a tensor field with a
 * tensorial diffusivity.
 *
 * The advection is upwind, the diffusion is not reconstructed.
 * The matrix is split into a diagonal block (3x3 times number of cells)
 * and an extra diagonal part (of dimension 2 times 3x3 the number of internal
 * faces).
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     coefbu        boundary condition array for the variable
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     fimp          part of the diagonal
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
cs_matrix_anisotropic_diffusion_tensor(const cs_mesh_t          *m,
                                       int                       iconvp,
                                       int                       idiffp,
                                       double                    thetap,
                                       const cs_real_66_t        coefbu[],
                                       const cs_real_66_t        cofbfu[],
                                       const cs_real_66_t        fimp[],
                                       const cs_real_t           i_massflux[],
                                       const cs_real_t           b_massflux[],
                                       const cs_real_66_t        i_visc[],
                                       const cs_real_t           b_visc[],
                                       cs_real_66_t    *restrict da,
                                       cs_real_662_t   *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 6; isou++) {
        for (int jsou = 0; jsou < 6; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        xa[face_id][0][jsou][isou] = 0.;
        xa[face_id][1][jsou][isou] = 0.;
      }
    }
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    double flui = 0.5*( i_massflux[face_id] -fabs(i_massflux[face_id]) );
    double fluj =-0.5*( i_massflux[face_id] +fabs(i_massflux[face_id]) );

    for (int isou = 0; isou < 6; isou++) {
      xa[face_id][0][isou][isou] = iconvp*flui;
      xa[face_id][1][isou][isou] = iconvp*fluj;
      for (int jsou = 0; jsou < 6; jsou++) {
        xa[face_id][0][jsou][isou] = thetap*( xa[face_id][0][jsou][isou]
                                         - idiffp*i_visc[face_id][jsou][isou]);
        xa[face_id][1][jsou][isou] = thetap*( xa[face_id][1][jsou][isou]
                                         - idiffp*i_visc[face_id][jsou][isou]);
      }
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    /* D_ii =  theta (m_ij)^+ - m_ij
     *      = -X_ij - (1-theta)*m_ij
     * D_jj = -theta (m_ij)^- + m_ij
     *      = -X_ji + (1-theta)*m_ij
     */
    for (int isou = 0; isou < 6; isou++) {
      da[ii][isou][isou] -= iconvp*(1. - thetap)*i_massflux[face_id];
      da[jj][isou][isou] += iconvp*(1. - thetap)*i_massflux[face_id];

      for (int jsou = 0; jsou < 6; jsou++) {
        da[ii][jsou][isou] -= xa[face_id][0][jsou][isou];
        da[jj][jsou][isou] -= xa[face_id][1][jsou][isou];
      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];
    double flui = 0.5*( b_massflux[face_id] -fabs(b_massflux[face_id]) );

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        /* D_ii = theta (m_f)^+ + theta B (m_f)^- - m_f
         *      = (theta B -1)*(m_f)^- - (1-theta)*(m_f)^+
         *      = theta*(B -1)*(m_f)^- - (1-theta)*m_f
         */
        if(isou == jsou) {
          da[ii][jsou][isou] += iconvp*( thetap*flui
                                        *(coefbu[face_id][jsou][isou]-1.)
                                       - (1. - thetap)*b_massflux[face_id])
                              + idiffp*thetap*b_visc[face_id]
                                      *cofbfu[face_id][jsou][isou];
        } else {
          da[ii][jsou][isou] += thetap*( iconvp*flui*coefbu[face_id][jsou][isou]
                                        +idiffp*b_visc[face_id]
                                         *cofbfu[face_id][jsou][isou] );
        }
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
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
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

    cs_lnum_t ii = i_face_cells[ifac][0];
    cs_lnum_t jj = i_face_cells[ifac][1];

    for (int isou = 0; isou < 3; isou++) {
      for (int jsou = 0; jsou < 3; jsou++) {
        da[ii][jsou][isou] -= xa[ifac][jsou][isou];
        da[jj][jsou][isou] -= xa[ifac][jsou][isou];
      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

    cs_lnum_t ii = b_face_cells[ifac];

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
 * \param[in]     m             pointer to mesh structure
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     thetap        weighting coefficient for the theta-scheme,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centered
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     cofbfu        boundary condition array for the variable flux
 *                               (implicit part - 3x3 tensor array)
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 * \param[out]    xa            extra interleaved diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_sym_matrix_anisotropic_diffusion_tensor(const cs_mesh_t           *m,
                                           int                       idiffp,
                                           double                    thetap,
                                           const cs_real_66_t        cofbfu[],
                                           const cs_real_66_t        fimp[],
                                           const cs_real_66_t        i_visc[],
                                           const cs_real_t           b_visc[],
                                           cs_real_66_t    *restrict da,
                                           cs_real_66_t    *restrict xa)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* 1. Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[cell_id][jsou][isou] = fimp[cell_id][jsou][isou];
      }
    }
  }
  if(n_cells_ext > n_cells) {
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 6; isou++) {
        for (int jsou = 0; jsou < 6; jsou++) {
          da[cell_id][jsou][isou] = 0.;
        }
      }
    }
  }

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {
    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        xa[ifac][jsou][isou] = 0.;
      }
    }
  }

  /* 2. Computation of extradiagonal terms */

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        xa[ifac][jsou][isou] = -thetap*idiffp*i_visc[ifac][jsou][isou];
      }
    }

  }

  /* 3. Contribution of the extra-diagonal terms to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_i_faces; ifac++) {

    cs_lnum_t ii = i_face_cells[ifac][0];
    cs_lnum_t jj = i_face_cells[ifac][1];

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[ii][jsou][isou] -= xa[ifac][jsou][isou];
        da[jj][jsou][isou] -= xa[ifac][jsou][isou];
      }
    }
  }

  /* 4. Contribution of border faces to the diagonal */

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

    cs_lnum_t ii = b_face_cells[ifac];

    for (int isou = 0; isou < 6; isou++) {
      for (int jsou = 0; jsou < 6; jsou++) {
        da[ii][jsou][isou] += thetap*idiffp*b_visc[ifac]
                                    *cofbfu[ifac][jsou][isou];
      }
    }

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
