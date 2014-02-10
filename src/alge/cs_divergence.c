/*============================================================================
 * Divergence operators.
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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_divergence.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_divergence.c

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
 * Wrapper to cs_divergence
 *----------------------------------------------------------------------------*/

void CS_PROCF (divmas, DIVMAS)
(
 const cs_int_t  *const   init,
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 cs_real_t                diverg[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  cs_divergence(m,
                *init,
                i_massflux,
                b_massflux,
                diverg);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_ext_force_flux
 *----------------------------------------------------------------------------*/

void CS_PROCF (projts, PROJTS)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_ext_force_flux(m,
                    fvq,
                    *init,
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
 * Wrapper to cs_ext_force_anisotropic_flux
 *----------------------------------------------------------------------------*/

void CS_PROCF (projtv, PROJTV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_int_t  *const   ircflp,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_6_t        viscel[],
 const cs_real_2_t        weighf[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_ext_force_anisotropic_flux(m,
                                fvq,
                                *init,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integrated mass flux on the cells.
 *
 * \f[
 * \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_divergence(const cs_mesh_t          *m,
              int                       init,
              const cs_real_t           i_massflux[],
              const cs_real_t           b_massflux[],
              cs_real_t       *restrict diverg)
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

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      diverg[cell_id] = 0.;
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > THR_MIN)
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      diverg[cell_id] = 0.;
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }


  /*==========================================================================
    2. Integration on internal faces
    ==========================================================================*/

  for (int g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0] - 1;
        cs_lnum_t jj = i_face_cells[face_id][1] - 1;

        diverg[ii] += i_massflux[face_id];
        diverg[jj] -= i_massflux[face_id];

      }
    }
  }


  /*==========================================================================
    3. Integration on border faces
    ==========================================================================*/

  for (int g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for if(m->n_b_faces > THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id] - 1;
        diverg[ii] += b_massflux[face_id];

      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Project the external source terms to the faces in coherence with
 * cs_face_diffusion_scalar for the improved hydrostatic pressure algorithm
 * (iphydr=1).
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
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
cs_ext_force_flux(const cs_mesh_t          *m,
                  cs_mesh_quantities_t     *fvq,
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

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
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

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0] - 1;
      cs_lnum_t jj = i_face_cells[face_id][1] - 1;

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

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id] - 1;
      double surfn = b_face_surf[face_id];
      double distbf = b_dist[face_id];

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

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0] - 1;
      cs_lnum_t jj = i_face_cells[face_id][1] - 1;

      double pnd = weight[face_id];

      double dijpfx = dijpf[face_id][0];
      double dijpfy = dijpf[face_id][1];
      double dijpfz = dijpf[face_id][2];

      double surfn = i_face_surf[face_id];

      /* Recompute II' and JJ' at this level */

      double diipx = i_face_cog[face_id][0]-cell_cen[ii][0]-(1.-pnd)*dijpfx;
      double diipy = i_face_cog[face_id][1]-cell_cen[ii][1]-(1.-pnd)*dijpfy;
      double diipz = i_face_cog[face_id][2]-cell_cen[ii][2]-(1.-pnd)*dijpfz;
      double djjpx = i_face_cog[face_id][0]-cell_cen[jj][0]+pnd*dijpfx;
      double djjpy = i_face_cog[face_id][1]-cell_cen[jj][1]+pnd*dijpfy;
      double djjpz = i_face_cog[face_id][2]-cell_cen[jj][2]+pnd*dijpfz;

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

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id] - 1;
      double surfn = b_face_surf[face_id];
      double distbf = b_dist[face_id];

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
/*!
 * \brief Project the external source terms to the faces in coherence with
 * cs_face_anisotropic_diffusion_scalar for the improved hydrostatic pressure
 * algorithm (iphydr=1).
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
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
cs_ext_force_anisotropic_flux(const cs_mesh_t          *m,
                              cs_mesh_quantities_t     *fvq,
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

  double diippf[3], djjppf[3];
  double visci[3][3], viscj[3][3];

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
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

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0] - 1;
      cs_lnum_t jj = i_face_cells[face_id][1] - 1;

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

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id] - 1;
      double surfn = b_face_surf[face_id];
      double distbf = b_dist[face_id];

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

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0] - 1;
      cs_lnum_t jj = i_face_cells[face_id][1] - 1;

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
      double fikdvi = weighf[face_id][0];

      /* II" = IF + FI" */
      for (int i = 0; i < 3; i++) {
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
      double fjkdvi = weighf[face_id][1];

      /* JJ" = JF + FJ" */
      for (int i = 0; i < 3; i++) {
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

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id] - 1;

      double surfn = b_face_surf[face_id];
      double distbf = b_dist[face_id];

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
