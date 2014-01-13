/*============================================================================
 * Finite volume operators
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

#include "cs_finite_volume_operators.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_finite_volume_operators.c

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
 * Wrapper to cs_finite_volume_divergence
 *----------------------------------------------------------------------------*/

void CS_PROCF (divmas, DIVMAS)
(
 const cs_int_t  *const   n_cells_ext,
 const cs_int_t  *const   n_cells,
 const cs_int_t  *const   n_i_faces,
 const cs_int_t  *const   n_b_faces,
 const cs_int_t  *const   init,
 const cs_int_t  *const   nfecra,
 const cs_lnum_2_t        i_face_cells[],
 const cs_int_t           b_face_cells[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 cs_real_t                diverg[])
{
  cs_finite_volume_divergence(*n_cells_ext,
                              *n_cells,
                              *n_i_faces,
                              *n_b_faces,
                              *init,
                              i_face_cells,
                              b_face_cells,
                              i_massflux,
                              b_massflux,
                              diverg);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_face_gradient_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmas, ITRMAS)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[])
{
  cs_finite_volume_face_gradient_scalar(*init,
                                        *inc,
                                        *imrgra,
                                        *iccocg,
                                        *nswrgp,
                                        *imligp,
                                        *iphydp,
                                        *iwarnp,
                                        *epsrgp,
                                        *climgp,
                                        *extrap,
                                        frcxt,
                                        pvar,
                                        coefap,
                                        coefbp,
                                        cofafp,
                                        cofbfp,
                                        i_visc,
                                        b_visc,
                                        viselx,
                                        visely,
                                        viselz,
                                        i_massflux,
                                        b_massflux);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_div_face_gradient_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrp, ITRGRP)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[],
 cs_real_t                diverg[])
{
  cs_finite_volume_div_face_gradient_scalar(*init,
                                                   *inc,
                                                   *imrgra,
                                                   *iccocg,
                                                   *nswrgp,
                                                   *imligp,
                                                   *iphydp,
                                                   *iwarnp,
                                                   *epsrgp,
                                                   *climgp,
                                                   *extrap,
                                                   frcxt,
                                                   pvar,
                                                   coefap,
                                                   coefbp,
                                                   cofafp,
                                                   cofbfp,
                                                   i_visc,
                                                   b_visc,
                                                   viselx,
                                                   visely,
                                                   viselz,
                                                   diverg);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_face_source_terms_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (projts, PROJTS)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_int_t  *const   nfecra,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[])
{
  cs_finite_volume_face_source_terms_scalar(*init,
                                            *nswrgu,
                                            frcxt,
                                            cofbfp,
                                            i_massflux,
                                            b_massflux,
                                            i_visc,
                                            b_visc,
                                            viselx,
                                            visely,
                                            viselz);

}

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_face_source_terms_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (projtv, PROJTV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   nfecra,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_6_t        viscel[],
 const cs_real_2_t        weighf[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[])
{
  cs_finite_volume_face_source_terms_vector(*init,
                                            *nswrgu,
                                            *ircflp,
                                            frcxt,
                                            cofbfp,
                                            i_visc,
                                            b_visc,
                                            viscel,
                                            weighf,
                                            i_massflux,
                                            b_massflux);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*! \brief This function adds the integrated mass flux on the cells.

  \f[
  \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
  \f]

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
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_face_cells  cell indexes of interior faces
 * \param[in]     b_face_cells  cell indexes of boundary faces
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_divergence(
                            int                       n_cells_ext,
                            int                       n_cells,
                            int                       n_i_faces,
                            int                       n_b_faces,
                            int                       init,
                            const cs_lnum_2_t         i_face_cells[],
                            const cs_int_t            b_face_cells[],
                            const cs_real_t           i_massflux[],
                            const cs_real_t           b_massflux[],
                            cs_real_t       *restrict diverg)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  /* Local variables */

  int cell_id, face_id, ii, jj, g_id, t_id;

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      diverg[cell_id] = 0.;
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > THR_MIN)
    for (cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      diverg[cell_id] = 0.;
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }


  /*==========================================================================
    2. Integration on internal faces
    ==========================================================================*/

  for (g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for private(face_id, ii, jj)
    for (t_id = 0; t_id < n_i_threads; t_id++) {
      for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        ii = i_face_cells[face_id][0] - 1;
        jj = i_face_cells[face_id][1] - 1;

        diverg[ii] += i_massflux[face_id];
        diverg[jj] -= i_massflux[face_id];

      }
    }
  }


  /*==========================================================================
    3. Integration on border faces
    ==========================================================================*/

  for (g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for private(face_id, ii) if(n_b_faces > THR_MIN)
    for (t_id = 0; t_id < n_b_threads; t_id++) {
      for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        ii = b_face_cells[face_id] - 1;
        diverg[ii] += b_massflux[face_id];

      }
    }
  }

}

/*----------------------------------------------------------------------------*/

/*! \brief This function updates the face mass flux with the face pressure
  (or pressure increment, or pressure double increment) gradient.

   \f[
   \dot{m}_\ij = \dot{m}_\ij
               - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
   \f]

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_face_gradient_scalar(
                                      int                       init,
                                      int                       inc,
                                      int                       imrgra,
                                      int                       iccocg,
                                      int                       nswrgp,
                                      int                       imligp,
                                      int                       iphydp,
                                      int                       iwarnp,
                                      double                    epsrgp,
                                      double                    climgp,
                                      double                    extrap,
                                      cs_real_3_t     *restrict frcxt,
                                      cs_real_t       *restrict pvar,
                                      const cs_real_t           coefap[],
                                      const cs_real_t           coefbp[],
                                      const cs_real_t           cofafp[],
                                      const cs_real_t           cofbfp[],
                                      const cs_real_t           i_visc[],
                                      const cs_real_t           b_visc[],
                                      const cs_real_t           viselx[],
                                      const cs_real_t           visely[],
                                      const cs_real_t           viselz[],
                                      cs_real_t       *restrict i_massflux,
                                      cs_real_t       *restrict b_massflux)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

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
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* Local variables */

  char var_name[32];
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  int face_id, ii, jj, g_id, t_id;
  double pfac,pip;
  double dpxf  , dpyf  , dpzf;
  double dijpfx, dijpfy, dijpfz;
  double diipbx, diipby, diipbz;
  double dijx  , dijy  , dijz;

  cs_real_3_t *grad;
  cs_real_3_t *visel;

  /*==========================================================================*/

  /* i_visc and visel carry similar information */

  /*============================================================================
    1. Initialization
    ==========================================================================*/

  BFT_MALLOC(visel, n_cells_ext, cs_real_3_t);
  for (ii = 0; ii < n_cells_ext; ii++) {
    visel[ii][0] = viselx[ii];
    visel[ii][1] = visely[ii];
    visel[ii][2] = viselz[ii];
  }

  if (init >= 1) {
#   pragma omp parallel for
    for (face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
#   pragma omp parallel for if(m->n_b_faces > THR_MIN)
    for (face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }
  } else if(init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  if (imrgra < 0)
    cs_gradient_type_by_imrgra(-imrgra,
                               &gradient_type,
                               &halo_type);

  snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  /* Handle parallelism and periodicity */

  if (halo != NULL)
    cs_halo_sync_var(halo, halo_type, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0] - 1;
          jj = i_face_cells[face_id][1] - 1;

          i_massflux[face_id] += i_visc[face_id]*(pvar[ii] -pvar[jj]);

        }
      }
    }

    /* Mass flow through boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, pfac)                       \
                 if(m->n_b_faces > THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id] - 1;
          pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

          b_massflux[face_id] += b_visc[face_id]*pfac;

        }
      }
    }

  }

  /*==========================================================================
    3. Update mass flux with reconstruction technics if the mesh is non
       orthogonal
    ==========================================================================*/

  if (nswrgp > 1) {

    /* Allocate a work array for the gradient calculation */
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    /* Compute gradient */

    cs_gradient_scalar(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       recompute_cocg,
                       nswrgp,
                       tr_dim,
                       iphydp,
                       iwarnp,
                       imligp,
                       epsrgp,
                       extrap,
                       climgp,
                       frcxt,
                       coefap,
                       coefbp,
                       pvar,
                       NULL, /* Weighted gradient */
                       grad);

    /* Handle parallelism and periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)visel, 3);
      if (cs_glob_mesh->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)visel, 3);
    }

    /* Mass flow through interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, dpxf, dpyf, dpzf, \
                                      dijpfx, dijpfy, dijpfz, dijx, dijy, dijz)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0] - 1;
          jj = i_face_cells[face_id][1] - 1;

          dpxf = 0.5*(visel[ii][0]*grad[ii][0] + visel[jj][0]*grad[jj][0]);
          dpyf = 0.5*(visel[ii][1]*grad[ii][1] + visel[jj][1]*grad[jj][1]);
          dpzf = 0.5*(visel[ii][2]*grad[ii][2] + visel[jj][2]*grad[jj][2]);

          dijpfx = dijpf[face_id][0];
          dijpfy = dijpf[face_id][1];
          dijpfz = dijpf[face_id][2];

          /*---> Dij = IJ - (IJ.N) N */
          dijx = (cell_cen[jj][0]-cell_cen[ii][0])-dijpfx;
          dijy = (cell_cen[jj][1]-cell_cen[ii][1])-dijpfy;
          dijz = (cell_cen[jj][2]-cell_cen[ii][2])-dijpfz;

          i_massflux[face_id] =  i_massflux[face_id]
                               + i_visc[face_id]*(pvar[ii] - pvar[jj])
                               +  (dpxf *dijx + dpyf*dijy + dpzf*dijz)
                                 * i_face_surf[face_id]/i_dist[face_id];

        }
      }
    }

    /* Mass flow through boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz,     \
                                      pip, pfac)                               \
                 if(m->n_b_faces > THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id] - 1;

          diipbx = diipb[face_id][0];
          diipby = diipb[face_id][1];
          diipbz = diipb[face_id][2];

          pip = pvar[ii] + grad[ii][0]*diipbx
                         + grad[ii][1]*diipby
                         + grad[ii][2]*diipbz;
          pfac = inc*cofafp[face_id] + cofbfp[face_id]*pip;

          b_massflux[face_id] += b_visc[face_id]*pfac;

        }
      }
    }

    /* Free memory */
    BFT_FREE(grad);
  }
  BFT_FREE(visel);
}

/*----------------------------------------------------------------------------*/

/*! \brief This function updates the cell mass flux divergence with the face
  pressure (or pressure increment, or pressure double increment) gradient.

  \f[
  \dot{m}_\ij = \dot{m}_\ij
              - \sum_j \Delta t \grad_\fij p \cdot \vect{S}_\ij
  \f]

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_div_face_gradient_scalar(
                                          int                       init,
                                          int                       inc,
                                          int                       imrgra,
                                          int                       iccocg,
                                          int                       nswrgp,
                                          int                       imligp,
                                          int                       iphydp,
                                          int                       iwarnp,
                                          double                    epsrgp,
                                          double                    climgp,
                                          double                    extrap,
                                          cs_real_3_t     *restrict frcxt,
                                          cs_real_t       *restrict pvar,
                                          const cs_real_t           coefap[],
                                          const cs_real_t           coefbp[],
                                          const cs_real_t           cofafp[],
                                          const cs_real_t           cofbfp[],
                                          const cs_real_t           i_visc[],
                                          const cs_real_t           b_visc[],
                                          const cs_real_t           viselx[],
                                          const cs_real_t           visely[],
                                          const cs_real_t           viselz[],
                                          cs_real_t       *restrict diverg)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

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
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* Local variables */

  char var_name[32];
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  int face_id, ii, jj, g_id, t_id;
  double pfac,pip;
  double dpxf  , dpyf  , dpzf  , i_massflux, b_massflux;
  double dijpfx, dijpfy, dijpfz;
  double diipbx, diipby, diipbz;
  double dijx  , dijy  , dijz;

  cs_real_3_t *grad;
  cs_real_3_t *visel;

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  BFT_MALLOC(visel, n_cells_ext, cs_real_3_t);
  for (ii = 0; ii < n_cells_ext; ii++) {
    visel[ii][0] = viselx[ii];
    visel[ii][1] = visely[ii];
    visel[ii][2] = viselz[ii];
  }

  if (init >= 1) {
#   pragma omp parallel for
    for (ii = 0; ii < n_cells_ext; ii++) {
      diverg[ii] = 0.;
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > THR_MIN)
    for (ii = n_cells; ii < n_cells_ext; ii++) {
      diverg[ii] = 0.;
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  if (imrgra < 0)
    cs_gradient_type_by_imrgra(imrgra,
                               &gradient_type,
                               &halo_type);

  snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  /* Handle parallelism and periodicity */

  if (halo != NULL)
    cs_halo_sync_var(halo, halo_type, pvar);


  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, i_massflux)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0] - 1;
          jj = i_face_cells[face_id][1] - 1;

          i_massflux = i_visc[face_id]*(pvar[ii] -pvar[jj]);
          diverg[ii] += i_massflux;
          diverg[jj] -= i_massflux;

        }
      }
    }

    /* Mass flow through boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, pfac, b_massflux)           \
                 if(m->n_b_faces > THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id] - 1;
          pfac = inc*cofafp[face_id] +cofbfp[face_id]*pvar[ii];

          b_massflux = b_visc[face_id]*pfac;
          diverg[ii] += b_massflux;

        }
      }
    }

  }


  /*==========================================================================
    3. Update mass flux with reconstruction technics if the mesh is non
       orthogonal
    ==========================================================================*/

  if (nswrgp > 1) {

    /* Allocate a work array for the gradient calculation */
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    /* Compute gradient */


    /* Compute gradient */

    cs_gradient_scalar(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       recompute_cocg,
                       nswrgp,
                       tr_dim,
                       iphydp,
                       iwarnp,
                       imligp,
                       epsrgp,
                       extrap,
                       climgp,
                       frcxt,
                       coefap,
                       coefbp,
                       pvar,
                       NULL, /* Weighted gradient */
                       grad);

    /* Handle parallelism and periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)visel, 3);
      if (cs_glob_mesh->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)visel, 3);
    }

    /* Mass flow through interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, dpxf, dpyf, dpzf,       \
                                      dijpfx, dijpfy, dijpfz,                  \
                                      dijx, dijy, dijz, i_massflux)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0] - 1;
          jj = i_face_cells[face_id][1] - 1;

          dijpfx = dijpf[face_id][0];
          dijpfy = dijpf[face_id][1];
          dijpfz = dijpf[face_id][2];

          /*---> Dij = IJ - (IJ.N) N */
          dijx = (cell_cen[jj][0]-cell_cen[ii][0])-dijpfx;
          dijy = (cell_cen[jj][1]-cell_cen[ii][1])-dijpfy;
          dijz = (cell_cen[jj][2]-cell_cen[ii][2])-dijpfz;

          dpxf = 0.5*(visel[ii][0]*grad[ii][0] + visel[jj][0]*grad[jj][0]);
          dpyf = 0.5*(visel[ii][1]*grad[ii][1] + visel[jj][1]*grad[jj][1]);
          dpzf = 0.5*(visel[ii][2]*grad[ii][2] + visel[jj][2]*grad[jj][2]);

          i_massflux =   i_visc[face_id]*(pvar[ii] -pvar[jj])
                       +  (dpxf*dijx + dpyf*dijy + dpzf*dijz)
                         *i_face_surf[face_id]/i_dist[face_id];

          diverg[ii] += i_massflux;
          diverg[jj] -= i_massflux;

        }
      }
    }

    /* Mass flow through boundary faces */

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz,     \
                                      pip, pfac, b_massflux)                   \
                 if(m->n_b_faces > THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id] - 1;

          diipbx = diipb[face_id][0];
          diipby = diipb[face_id][1];
          diipbz = diipb[face_id][2];

          pip = pvar[ii] + grad[ii][0]*diipbx
                         + grad[ii][1]*diipby
                         + grad[ii][2]*diipbz;
          pfac = inc*cofafp[face_id] +cofbfp[face_id]*pip;

          b_massflux = b_visc[face_id]*pfac;
          diverg[ii] += b_massflux;

        }
      }
    }

    /* Free memory */
    BFT_FREE(grad);
  }
  BFT_FREE(visel);
}

/*----------------------------------------------------------------------------*/

/*! \brief This function projects the external source terms to the faces
  in coherence with itrmas.f90 for the improved hydrostatic pressure
  algorithm (iphydr=1).

*/
/*----------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgu        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_face_source_terms_scalar(
                                          int                       init,
                                          int                       nswrgu,
                                          const cs_real_3_t         frcxt[],
                                          const cs_real_t           cofbfp[],
                                          cs_real_t       *restrict i_massflux,
                                          cs_real_t       *restrict b_massflux,
                                          const cs_real_t           i_visc[],
                                          const cs_real_t           b_visc[],
                                          const cs_real_t           viselx[],
                                          const cs_real_t           visely[],
                                          const cs_real_t           viselz[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;

  /* Local variables */

  int face_id, ii, jj;
  double dijpfx,dijpfy,dijpfz;
  double diipx,diipy,diipz;
  double djjpx,djjpy,djjpz;
  double distbf,surfn;

  double pnd;

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    for (face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
    for (face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Mass flow through interior faces */

    for (face_id = 0; face_id < m->n_i_faces; face_id++) {

      ii = i_face_cells[face_id][0] - 1;
      jj = i_face_cells[face_id][1] - 1;

      i_massflux[face_id] =  i_massflux[face_id]
                           + i_visc[face_id]*(
                                              ( i_face_cog[face_id][0]
                                               -cell_cen[ii][0] )*frcxt[ii][0]
                                             +( i_face_cog[face_id][1]
                                               -cell_cen[ii][1] )*frcxt[ii][1]
                                             +( i_face_cog[face_id][2]
                                               -cell_cen[ii][2] )*frcxt[ii][2]
                                             -( i_face_cog[face_id][0]
                                               -cell_cen[jj][0] )*frcxt[jj][0]
                                             -( i_face_cog[face_id][1]
                                               -cell_cen[jj][1] )*frcxt[jj][1]
                                             -( i_face_cog[face_id][2]
                                               -cell_cen[jj][2] )*frcxt[jj][2]
                                             );

    }

    /* Mass flux through boundary faces */

    for (face_id = 0; face_id < m->n_b_faces; face_id++) {

      ii = b_face_cells[face_id] - 1;
      surfn = b_face_surf[face_id];
      distbf = b_dist[face_id];

      b_massflux[face_id] =  b_massflux[face_id]
                           +  b_visc[face_id]*distbf/surfn
                             *cofbfp[face_id]
                             *( frcxt[ii][0]*b_face_normal[face_id][0]
                               +frcxt[ii][1]*b_face_normal[face_id][1]
                               +frcxt[ii][2]*b_face_normal[face_id][2] );

    }

  /*==========================================================================
    3. Update mass flux with reconstruction technics
    ==========================================================================*/

  } else {


    /* Mass flux through interior faces */

    for (face_id = 0; face_id < m->n_i_faces; face_id++) {

      ii = i_face_cells[face_id][0] - 1;
      jj = i_face_cells[face_id][1] - 1;

      pnd = weight[face_id];

      dijpfx = dijpf[face_id][0];
      dijpfy = dijpf[face_id][1];
      dijpfz = dijpf[face_id][2];

      surfn = i_face_surf[face_id];

      /* Recompute II' and JJ' at this level */

      diipx = i_face_cog[face_id][0]-cell_cen[ii][0]-(1.-pnd)*dijpfx;
      diipy = i_face_cog[face_id][1]-cell_cen[ii][1]-(1.-pnd)*dijpfy;
      diipz = i_face_cog[face_id][2]-cell_cen[ii][2]-(1.-pnd)*dijpfz;
      djjpx = i_face_cog[face_id][0]-cell_cen[jj][0]+pnd*dijpfx;
      djjpy = i_face_cog[face_id][1]-cell_cen[jj][1]+pnd*dijpfy;
      djjpz = i_face_cog[face_id][2]-cell_cen[jj][2]+pnd*dijpfz;

      i_massflux[face_id] =  i_massflux[face_id]
                           + i_visc[face_id]*(
                                               ( i_face_cog[face_id][0]
                                                -cell_cen[ii][0] )*frcxt[ii][0]
                                              +( i_face_cog[face_id][1]
                                                -cell_cen[ii][1] )*frcxt[ii][1]
                                              +( i_face_cog[face_id][2]
                                                -cell_cen[ii][2] )*frcxt[ii][2]
                                              -( i_face_cog[face_id][0]
                                                -cell_cen[jj][0] )*frcxt[jj][0]
                                              -( i_face_cog[face_id][1]
                                                -cell_cen[jj][1] )*frcxt[jj][1]
                                              -( i_face_cog[face_id][2]
                                                -cell_cen[jj][2] )*frcxt[jj][2]
                                              )
                            + surfn/i_dist[face_id]*0.5
                             *( (djjpx-diipx)*( viselx[ii]*frcxt[ii][0]
                                               +viselx[jj]*frcxt[jj][0] )
                               +(djjpy-diipy)*( visely[ii]*frcxt[ii][1]
                                               +visely[jj]*frcxt[jj][1] )
                               +(djjpz-diipz)*( viselz[ii]*frcxt[ii][2]
                                               +viselz[jj]*frcxt[jj][2] )
                              );

    }


    /* Mass flux through boundary faces */

    for (face_id = 0; face_id < m->n_b_faces; face_id++) {

      ii = b_face_cells[face_id] - 1;
      surfn = b_face_surf[face_id];
      distbf = b_dist[face_id];

      b_massflux[face_id] = b_massflux[face_id]
                           + b_visc[face_id]*distbf/surfn
                            *cofbfp[face_id]
                            *( frcxt[ii][0]*b_face_normal[face_id][0]
                              +frcxt[ii][1]*b_face_normal[face_id][1]
                              +frcxt[ii][2]*b_face_normal[face_id][2] );

    }
  }

}

/*----------------------------------------------------------------------------*/

/*! \brief This function projects the external source terms to the faces
  in coherence with itrmav.f90 for the improved hydrostatic pressure
  algorithm (iphydr=1).

*/
/*----------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_face_source_terms_vector(
                                          int                       init,
                                          int                       nswrgp,
                                          int                       ircflp,
                                          const cs_real_3_t         frcxt[],
                                          const cs_real_t           cofbfp[],
                                          const cs_real_t           i_visc[],
                                          const cs_real_t           b_visc[],
                                          const cs_real_6_t         viscel[],
                                          const cs_real_2_t         weighf[],
                                          cs_real_t       *restrict i_massflux,
                                          cs_real_t       *restrict b_massflux)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;

  /* Local variables */

  int face_id, ii, jj, i;
  double distbf,surfn;
  double diippf[3], djjppf[3];
  double visci[3][3], viscj[3][3];
  double fikdvi, fjkdvi;

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    for (face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
    for (face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* ---> Contribution from interior faces */

    for (face_id = 0; face_id < m->n_i_faces; face_id++) {

      ii = i_face_cells[face_id][0] - 1;
      jj = i_face_cells[face_id][1] - 1;

      i_massflux[face_id] =  i_massflux[face_id]+ i_visc[face_id]*(
                                               ( i_face_cog[face_id][0]
                                                -cell_cen[ii][0])*frcxt[ii][0]
                                              +( i_face_cog[face_id][1]
                                                -cell_cen[ii][1])*frcxt[ii][1]
                                              +( i_face_cog[face_id][2]
                                                -cell_cen[ii][2])*frcxt[ii][2]
                                              -( i_face_cog[face_id][0]
                                                -cell_cen[jj][0])*frcxt[jj][0]
                                              -( i_face_cog[face_id][1]
                                                -cell_cen[jj][1])*frcxt[jj][1]
                                              -( i_face_cog[face_id][2]
                                                -cell_cen[jj][2])*frcxt[jj][2]
                                              );

    }

    /* ---> Contribution from boundary faces */

    for (face_id = 0; face_id < m->n_b_faces; face_id++) {

      ii = b_face_cells[face_id] - 1;
      surfn = b_face_surf[face_id];
      distbf = b_dist[face_id];

      b_massflux[face_id] =  b_massflux[face_id]
                           + b_visc[face_id]*distbf/surfn
                            *cofbfp[face_id]
                            *( frcxt[ii][0]*b_face_normal[face_id][0]
                              +frcxt[ii][1]*b_face_normal[face_id][1]
                              +frcxt[ii][2]*b_face_normal[face_id][2] );

    }

    /*========================================================================
      3. Update mass flux with reconstruction technics
      ========================================================================*/

  } else {

    /* ---> Contribution from interior faces */

    for (face_id = 0; face_id < m->n_i_faces; face_id++) {

      ii = i_face_cells[face_id][0] - 1;
      jj = i_face_cells[face_id][1] - 1;

      /* Recompute II' and JJ' at this level */

      visci[0][0] = viscel[ii][0];
      visci[1][1] = viscel[ii][1];
      visci[2][2] = viscel[ii][2];
      visci[1][0] = viscel[ii][3];
      visci[0][1] = viscel[ii][3];
      visci[2][1] = viscel[ii][4];
      visci[1][2] = viscel[ii][4];
      visci[2][0] = viscel[ii][5];
      visci[0][2] = viscel[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      fikdvi = weighf[face_id][0];

      /* II" = IF + FI" */
      for (i = 0; i < 3; i++) {
        diippf[i] =  i_face_cog[face_id][i]-cell_cen[ii][i]
                   - fikdvi*(  visci[0][i]*i_face_normal[face_id][0]
                             + visci[1][i]*i_face_normal[face_id][1]
                             + visci[2][i]*i_face_normal[face_id][2] );
      }

      viscj[0][0] = viscel[jj][0];
      viscj[1][1] = viscel[jj][1];
      viscj[2][2] = viscel[jj][2];
      viscj[1][0] = viscel[jj][3];
      viscj[0][1] = viscel[jj][3];
      viscj[2][1] = viscel[jj][4];
      viscj[1][2] = viscel[jj][4];
      viscj[2][0] = viscel[jj][5];
      viscj[0][2] = viscel[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      fjkdvi = weighf[face_id][1];

      /* JJ" = JF + FJ" */
      for (i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi*(  viscj[0][i]*i_face_normal[face_id][0]
                              + viscj[1][i]*i_face_normal[face_id][1]
                              + viscj[2][i]*i_face_normal[face_id][2] );
      }

      i_massflux[face_id] =  i_massflux[face_id]
                            + i_visc[face_id]
                             *(  frcxt[ii][0]*( i_face_cog[face_id][0]
                                               -cell_cen[ii][0] )
                               + frcxt[ii][1]*( i_face_cog[face_id][1]
                                               -cell_cen[ii][1] )
                               + frcxt[ii][2]*( i_face_cog[face_id][2]
                                               -cell_cen[ii][2] )
                               - frcxt[jj][0]*( i_face_cog[face_id][0]
                                               -cell_cen[jj][0] )
                               - frcxt[jj][1]*( i_face_cog[face_id][1]
                                               -cell_cen[jj][1] )
                               - frcxt[jj][2]*( i_face_cog[face_id][2]
                                               -cell_cen[jj][2] )
                              )
                            + i_visc[face_id]*ircflp
                             *(- frcxt[ii][0]*diippf[0]
                               - frcxt[ii][1]*diippf[1]
                               - frcxt[ii][2]*diippf[2]
                               + frcxt[jj][0]*djjppf[0]
                               + frcxt[jj][1]*djjppf[1]
                               + frcxt[jj][2]*djjppf[2]
                              );

    }

    /* ---> Contribution from boundary faces */

    for (face_id = 0; face_id < m->n_b_faces; face_id++) {

      ii = b_face_cells[face_id] - 1;

      surfn = b_face_surf[face_id];
      distbf = b_dist[face_id];

      /* FIXME: wrong if dirichlet and viscel is really a tensor */
      b_massflux[face_id] =  b_massflux[face_id]
                            + b_visc[face_id]*distbf/surfn*cofbfp[face_id]
                             *(  frcxt[ii][0]*b_face_normal[face_id][0]
                               + frcxt[ii][1]*b_face_normal[face_id][1]
                               + frcxt[ii][2]*b_face_normal[face_id][2] );

    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
