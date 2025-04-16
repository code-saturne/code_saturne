/*============================================================================
 * Divergence operators.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_dispatch.h"
#include "alge/cs_gradient_boundary.h"
#include "base/cs_porous_model.h"
#include "base/cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_divergence.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_divergence.cpp

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
 * Wrapper to cs_divergence
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add \f$ \rho \vect{u} \cdot \vect{s}_\ij\f$ to
 * the mass flux \f$ \dot{m}_\ij \f$.
 *
 * For the reconstruction, \f$ \gradt \left(\rho \vect{u} \right) \f$ is
 * computed with the following approximated boundary conditions:
 *  - \f$ \vect{a}_{\rho u} = \rho_\fib \vect{a}_u \f$
 *  - \f$ \tens{b}_{\rho u} = \tens{b}_u \f$
 *
 * For the mass flux at the boundary we have:
 * \f[
 * \dot{m}_\ib = \left[ \rho_\fib \vect{a}_u  + \rho_\fib \tens{b}_u \vect{u}
 * + \tens{b}_u \left(\gradt \vect{u} \cdot \vect{\centi \centip}\right)\right]
 * \cdot \vect{s}_\ij
 * \f]
 * The last equation uses some approximations detailed in the theory guide.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     f             pointer to field
 * \param[in]     itypfl        indicator (take rho into account or not)
 *                               - 1 compute \f$ \rho\vect{u}\cdot\vect{s} \f$
 *                               - 0 compute \f$ \vect{u}\cdot\vect{s} \f$
 * \param[in]     iflmb0        the mass flux is set to 0 on symmetries if = 1
 * \param[in]     init          the mass flux is initialized to 0 if > 0
 * \param[in]     inc           indicator
 *                               - 0 solve an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     nswrgu        number of sweeps for the reconstruction
 *                               of the gradients
 * \param[in]     imligu        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thanks to neighbooring gradients
 *                               - = 1 thanks to the mean gradient
 * \param[in]     iwarnu        verbosity
 * \param[in]     epsrgu        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgu        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     rom           cell density
 * \param[in]     romb          density at boundary faces
 * \param[in]     vel           vector variable
 * \param[in]     bc_coeffs_v   BC structure for the vector variable
 * \param[in,out] i_massflux    mass flux at interior faces \f$ \dot{m}_\fij \f$
 * \param[in,out] b_massflux    mass flux at boundary faces \f$ \dot{m}_\fib \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_mass_flux(const cs_mesh_t             *m,
             const cs_mesh_quantities_t  *fvq,
             int                          f_id,
             int                          itypfl,
             int                          iflmb0,
             int                          init,
             int                          inc,
             int                          imrgra,
             int                          nswrgu,
             cs_gradient_limit_t          imligu,
             int                          iwarnu,
             double                       epsrgu,
             double                       climgu,
             const cs_real_t              rom[],
             const cs_real_t              romb[],
             const cs_real_3_t            vel[],
             cs_field_bc_coeffs_t        *bc_coeffs_v,
             cs_real_t          *restrict i_massflux,
             cs_real_t          *restrict b_massflux)
{
  cs_real_3_t *coefav = (cs_real_3_t *)bc_coeffs_v->a;
  cs_real_33_t *coefbv = (cs_real_33_t *)bc_coeffs_v->b;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  cs_real_2_t *i_f_face_factor;
  cs_real_t *b_f_face_factor;

  /* Parallel or device dispatch */

  cs_dispatch_context ctx, ctx_c;
  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  ctx_c.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_c.use_gpu())
    ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  const bool on_device = ctx.use_gpu();

  /* Local variables */

  /* Discontinuous porous treatment */

  int is_p = 0; /* Is porous? */
  cs_real_2_t *_i_f_face_factor = nullptr;
  cs_real_t *_b_f_face_factor = nullptr;

  if (cs_glob_porous_model == 3) {
    i_f_face_factor = fvq->i_f_face_factor;
    b_f_face_factor = fvq->b_f_face_factor;
    is_p = 1;
  }
  else {
    CS_MALLOC_HD(_i_f_face_factor, 1, cs_real_2_t, amode);
    CS_MALLOC_HD(_b_f_face_factor, 1, cs_real_t, amode);
    i_f_face_factor = _i_f_face_factor;
    b_f_face_factor = _b_f_face_factor;
    i_f_face_factor[0][0] = i_f_face_factor[0][1] = 1.0;
    b_f_face_factor[0] = 1.0;
  }

  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_3_t *restrict dofij = fvq->dofij;

  char var_name[64];

  cs_field_t *f;

  cs_real_33_t *grdqdm = nullptr;
  cs_real_3_t *qdm, *f_momentum;

  CS_MALLOC_HD(qdm, n_cells_ext, cs_real_3_t, amode);
  CS_MALLOC_HD(f_momentum, n_b_faces, cs_real_3_t, amode);

  cs_field_bc_coeffs_t bc_coeffs_v_loc;
  cs_field_bc_coeffs_shallow_copy(bc_coeffs_v, &bc_coeffs_v_loc);

  CS_MALLOC_HD(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t, amode);
  cs_real_3_t *coefaq = (cs_real_3_t *)bc_coeffs_v_loc.a;

  /*==========================================================================
    1.  Initialization
    ==========================================================================*/

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 63, "%s", f->name);
  }
  else {
    strncpy(var_name, "[momentum]", 63);
  }
  var_name[63] = '\0';

  /* Momentum computation */

  if (init == 1) {
    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      i_massflux[face_id] = 0.;
    });
    ctx_c.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      b_massflux[face_id] = 0.;
    });

  }
  else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = rom[cell_id]*vel[cell_id][isou];
        }
      });
      /* With porosity */
    }
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = rom[cell_id]*vel[cell_id][isou]*porosi[cell_id];
        }
      });
      /* With anisotropic porosity */
    }
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        qdm[cell_id][0] = (  porosf[cell_id][0]*vel[cell_id][0]
                           + porosf[cell_id][3]*vel[cell_id][1]
                           + porosf[cell_id][5]*vel[cell_id][2])
                          * rom[cell_id];
        qdm[cell_id][1] = (  porosf[cell_id][3]*vel[cell_id][0]
                           + porosf[cell_id][1]*vel[cell_id][1]
                           + porosf[cell_id][4]*vel[cell_id][2])
                          * rom[cell_id];
        qdm[cell_id][2] = (  porosf[cell_id][5]*vel[cell_id][0]
                           + porosf[cell_id][4]*vel[cell_id][1]
                           + porosf[cell_id][2]*vel[cell_id][2])
                          * rom[cell_id];
      });
    }

    /* Velocity flux */
  }
  else {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = vel[cell_id][isou];
        }
      });
      /* With porosity */
    }
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = vel[cell_id][isou]*porosi[cell_id];
        }
      });
      /* With anisotropic porosity */
    }
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        qdm[cell_id][0] = porosf[cell_id][0]*vel[cell_id][0]
                        + porosf[cell_id][3]*vel[cell_id][1]
                        + porosf[cell_id][5]*vel[cell_id][2];
        qdm[cell_id][1] = porosf[cell_id][3]*vel[cell_id][0]
                        + porosf[cell_id][1]*vel[cell_id][1]
                        + porosf[cell_id][4]*vel[cell_id][2];
        qdm[cell_id][2] = porosf[cell_id][5]*vel[cell_id][0]
                        + porosf[cell_id][4]*vel[cell_id][1]
                        + porosf[cell_id][2]*vel[cell_id][2];
      });
    }
  }

  ctx_c.wait();

  /* Periodicity and parallelism treatment */

  if (halo != nullptr)
    cs_halo_sync_r(halo, halo_type, on_device, qdm);

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = romb[face_id]*coefav[face_id][isou];
          f_momentum[face_id][isou] = romb[face_id]*vel[cell_id][isou];
        }
      });
    } /* With porosity */
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = romb[face_id]
                                 *coefav[face_id][isou]*porosi[cell_id];
          f_momentum[face_id][isou] =  romb[face_id]*vel[cell_id][isou]
                                      *porosi[cell_id];
        }
      });

    } /* With anisotropic porosity */
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        coefaq[face_id][0] = (  porosf[cell_id][0]*coefav[face_id][0]
                              + porosf[cell_id][3]*coefav[face_id][1]
                              + porosf[cell_id][5]*coefav[face_id][2])
                             * romb[face_id];
        coefaq[face_id][1] = (  porosf[cell_id][3]*coefav[face_id][0]
                              + porosf[cell_id][1]*coefav[face_id][1]
                              + porosf[cell_id][4]*coefav[face_id][2])
                             * romb[face_id];
        coefaq[face_id][2] = (  porosf[cell_id][5]*coefav[face_id][0]
                              + porosf[cell_id][4]*coefav[face_id][1]
                              + porosf[cell_id][2]*coefav[face_id][2])
                             * romb[face_id];
        f_momentum[face_id][0] = (  porosf[cell_id][0]*vel[cell_id][0]
                                  + porosf[cell_id][3]*vel[cell_id][1]
                                  + porosf[cell_id][5]*vel[cell_id][2])
                                 * romb[face_id];
        f_momentum[face_id][1] = (  porosf[cell_id][3]*vel[cell_id][0]
                                  + porosf[cell_id][1]*vel[cell_id][1]
                                  + porosf[cell_id][4]*vel[cell_id][2])
                                 * romb[face_id];
        f_momentum[face_id][2] = (  porosf[cell_id][5]*vel[cell_id][0]
                                  + porosf[cell_id][4]*vel[cell_id][1]
                                  + porosf[cell_id][2]*vel[cell_id][2])
                                 * romb[face_id];
      });
    }

    /* Velocity flux */
  }
  else {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou];
          f_momentum[face_id][isou] = vel[cell_id][isou];
        }
      });
    } /* With porosity */
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou]*porosi[cell_id];
          f_momentum[face_id][isou] = vel[cell_id][isou]*porosi[cell_id];
        }
      });
    } /* With anisotropic porosity */
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        coefaq[face_id][0] =   porosf[cell_id][0]*coefav[face_id][0]
                             + porosf[cell_id][3]*coefav[face_id][1]
                             + porosf[cell_id][5]*coefav[face_id][2];
        coefaq[face_id][1] =   porosf[cell_id][3]*coefav[face_id][0]
                             + porosf[cell_id][1]*coefav[face_id][1]
                             + porosf[cell_id][4]*coefav[face_id][2];
        coefaq[face_id][2] =   porosf[cell_id][5]*coefav[face_id][0]
                             + porosf[cell_id][4]*coefav[face_id][1]
                             + porosf[cell_id][2]*coefav[face_id][2];
        f_momentum[face_id][0] = (  porosf[cell_id][0]*vel[cell_id][0]
                                  + porosf[cell_id][3]*vel[cell_id][1]
                                  + porosf[cell_id][5]*vel[cell_id][2]);
        f_momentum[face_id][1] = (  porosf[cell_id][3]*vel[cell_id][0]
                                  + porosf[cell_id][1]*vel[cell_id][1]
                                  + porosf[cell_id][4]*vel[cell_id][2]);
        f_momentum[face_id][2] = (  porosf[cell_id][5]*vel[cell_id][0]
                                  + porosf[cell_id][4]*vel[cell_id][1]
                                  + porosf[cell_id][2]*vel[cell_id][2]);
      });
    }

  }

  /*==========================================================================
    2. Compute mass flux without recontructions
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];
      cs_lnum_t _p = is_p*face_id;
      cs_real_t w_i = weight[face_id] * i_f_face_factor[_p][0];
      cs_real_t w_j = (1. - weight[face_id]) * i_f_face_factor[_p][1];
      /* u, v, w Components */
      cs_real_t q = 0;
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        q +=   (w_i * qdm[ii][isou] + w_j * qdm[jj][isou])
             * i_face_u_normal[face_id][isou];
      }
      i_massflux[face_id] += q * i_face_surf[face_id];
    });

    /* Boundary faces */

    ctx_c.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t _p = is_p*face_id;
      /* u, v, w Components */
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        double pfac = inc*coefaq[face_id][isou];

        /* coefbv is a matrix */
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
          pfac += coefbv[face_id][jsou][isou]*f_momentum[face_id][jsou];
        }
        pfac *= b_f_face_factor[_p];

        b_massflux[face_id] +=   pfac * b_face_surf[face_id]
                               * b_face_u_normal[face_id][isou];
      }
    });
  }

  ctx.wait();
  ctx_c.wait();

  /*==========================================================================
    4. Compute mass flux with reconstruction method if the mesh is
       non orthogonal
    ==========================================================================*/

  if (nswrgu > 1) {

    CS_MALLOC_HD(grdqdm, n_cells_ext, cs_real_33_t, amode);

    /* Computation of momentum gradient
       (vectorial gradient, the periodicity has already been treated) */

    cs_gradient_vector(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       nswrgu,
                       iwarnu,
                       imligu,
                       epsrgu,
                       climgu,
                       &bc_coeffs_v_loc,
                       qdm,
                       nullptr, /* weighted gradient */
                       nullptr, /* cpl */
                       grdqdm);

    /* Mass flow through interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];
      cs_lnum_t _p = is_p*face_id;

      const cs_real_t *f_dofij = dofij[face_id];

      cs_real_t w_i = weight[face_id] * i_f_face_factor[_p][0];
      cs_real_t w_j = (1. - weight[face_id]) * i_f_face_factor[_p][1];

      /* Terms along U, V, W */

      cs_real_t q = 0.;
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        q     /* Non-reconstructed term */
          += (w_i * qdm[ii][isou] + w_j * qdm[jj][isou]

                /*  --->     ->    -->      ->
                 *  (Grad(rho U ) . OFij ) . Sij
                 * FIXME for discontinuous porous modelling */
              + 0.5*(  (grdqdm[ii][isou][0] + grdqdm[jj][isou][0])*f_dofij[0]
                     + (grdqdm[ii][isou][1] + grdqdm[jj][isou][1])*f_dofij[1]
                     + (grdqdm[ii][isou][2] + grdqdm[jj][isou][2])*f_dofij[2]))
             * i_face_u_normal[face_id][isou];
      }
      i_massflux[face_id] += q * i_face_surf[face_id];

    });

     /* Mass flow through boundary faces */

    ctx_c.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = b_face_cells[face_id];
      cs_lnum_t _p = is_p*face_id;

      /* Terms along U, V, W */
      cs_real_t q = 0;
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        cs_real_t pfac = inc*coefaq[face_id][isou];

        /* coefu is a matrix */
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {

          double pip =   f_momentum[face_id][jsou]
                       + grdqdm[ii][jsou][0]*diipb[face_id][0]
                       + grdqdm[ii][jsou][1]*diipb[face_id][1]
                       + grdqdm[ii][jsou][2]*diipb[face_id][2];

          pfac += coefbv[face_id][jsou][isou]*pip;

        }

        q += pfac * b_face_u_normal[face_id][isou];
      }
      b_massflux[face_id] += q * b_f_face_factor[_p] * b_face_surf[face_id];
    });
  }

  /*==========================================================================
    6. Here, we make sure that the mass flux is null at the boundary faces of
       type symmetry and coupled walls.
    ==========================================================================*/

  if (iflmb0 == 1) {
    /* Force flumab to 0 for velocity */
    cs_lnum_t *b_sym_flag = fvq->b_sym_flag;
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      if (b_sym_flag[face_id] == 0) {
        b_massflux[face_id] = 0.;
      }
    });
  }

  ctx.wait();
  ctx_c.wait();

  CS_FREE_HD(grdqdm);
  CS_FREE_HD(qdm);
  CS_FREE_HD(f_momentum);
  CS_FREE_HD(_i_f_face_factor);
  CS_FREE_HD(_b_f_face_factor);

  coefaq = nullptr;
  cs_field_bc_coeffs_free_copy(bc_coeffs_v, &bc_coeffs_v_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add \f$ \rho \tens{r} \vect{s}_\ij\f$ to a flux.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     f_id          field id (or -1)
 * \param[in]     itypfl        indicator (take rho into account or not)
 *                               - 1 compute \f$ \rho\vect{u}\cdot\vect{s} \f$
 *                               - 0 compute \f$ \vect{u}\cdot\vect{s} \f$
 * \param[in]     iflmb0        the mass flux is set to 0 on symmetries if = 1
 * \param[in]     init          the mass flux is initialized to 0 if > 0
 * \param[in]     inc           indicator
 *                               - 0 solve an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     nswrgu        number of sweeps for the reconstruction
 *                               of the gradients
 * \param[in]     imligu        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thanks to neighbooring gradients
 *                               - = 1 thanks to the mean gradient
 * \param[in]     iwarnu        verbosity
 * \param[in]     epsrgu        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgu        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     c_rho         cell density
 * \param[in]     b_rho         density at boundary faces
 * \param[in]     c_var         variable
 * \param[in]     bc_coeffs_ts  boundary condition structure for the variable
 * \param[in,out] i_massflux    mass flux at interior faces \f$ \dot{m}_\fij \f$
 * \param[in,out] b_massflux    mass flux at boundary faces \f$ \dot{m}_\fib \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_tensor_face_flux(const cs_mesh_t             *m,
                    const cs_mesh_quantities_t  *fvq,
                    int                          f_id,
                    int                          itypfl,
                    int                          iflmb0,
                    int                          init,
                    int                          inc,
                    int                          imrgra,
                    int                          nswrgu,
                    cs_gradient_limit_t          imligu,
                    int                          iwarnu,
                    double                       epsrgu,
                    double                       climgu,
                    const cs_real_t              c_rho[],
                    const cs_real_t              b_rho[],
                    const cs_real_6_t            c_var[],
                    const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                    cs_real_3_t        *restrict i_massflux,
                    cs_real_3_t        *restrict b_massflux)
{
  cs_real_6_t  *coefav = (cs_real_6_t  *)bc_coeffs_ts->a;
  cs_real_66_t *coefbv = (cs_real_66_t *)bc_coeffs_ts->b;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;
  const cs_real_3_t *restrict dofij = fvq->dofij;

  /* Local variables */

  char var_name[64];

  cs_real_63_t *c_grad_mvar = nullptr;

  cs_field_t *f;

  /* Parallel or device dispatch */

  cs_dispatch_context ctx, ctx_c;
  cs_alloc_mode_t amode = ctx.alloc_mode(true);
  ctx_c.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_c.use_gpu())
    ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  const bool on_device = ctx.use_gpu();

  cs_real_6_t *c_mass_var, *b_mass_var;
  CS_MALLOC_HD(c_mass_var, n_cells_ext, cs_real_6_t, amode);
  CS_MALLOC_HD(b_mass_var, m->n_b_faces, cs_real_6_t, amode);

  cs_field_bc_coeffs_t bc_coeffs_ts_loc;
  cs_field_bc_coeffs_shallow_copy(bc_coeffs_ts, &bc_coeffs_ts_loc);

  CS_MALLOC_HD(bc_coeffs_ts_loc.a, 6*m->n_b_faces, cs_real_t, amode);
  cs_real_6_t *coefaq = (cs_real_6_t *)bc_coeffs_ts_loc.a;

  /*==========================================================================
    1.  Initialization
    ==========================================================================*/

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 63, "%s", f->name);
  }
  else {
    strncpy(var_name, "[tensor face flux]", 63);
  }
  var_name[63] = '\0';

  /* ---> Momentum computation */

  if (init == 1) {
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        i_massflux[face_id][i] = 0.;
    });
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        b_massflux[face_id][i] = 0.;
    });

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
         for (cs_lnum_t isou = 0; isou < 6; isou++) {
           c_mass_var[cell_id][isou] = c_rho[cell_id]*c_var[cell_id][isou];
         }
      });
    }
    /* With porosity */
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] =   c_rho[cell_id]*c_var[cell_id][isou]
                                      * porosi[cell_id];
        }
      });
    }
    /* With anisotropic porosity */
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               c_mass_var[cell_id]);

        for (cs_lnum_t isou = 0; isou < 6; isou++)
          c_mass_var[cell_id][isou] *= c_rho[cell_id];
      });
    }

  }

  /* Velocity flux */
  else {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] = c_var[cell_id][isou];
        }
      });
    }
    /* With porosity */
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] = c_var[cell_id][isou]*porosi[cell_id];
        }
      });
    }
    /* With anisotropic porosity */
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               c_mass_var[cell_id]);
      });
    }
  }

  ctx_c.wait();

  /* Periodicity and parallelism treatment */

  if (halo != nullptr)
    cs_halo_sync_r(halo, halo_type, on_device, c_mass_var);

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = b_rho[face_id]*coefav[face_id][isou];
          b_mass_var[face_id][isou] = b_rho[face_id]*c_var[cell_id][isou];
        }
      });
    }
    /* With porosity */
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = b_rho[face_id]
                                 *coefav[face_id][isou]*porosi[cell_id];
          b_mass_var[face_id][isou] =   b_rho[face_id]*c_var[cell_id][isou]
                                      * porosi[cell_id];
        }
      });
    }
    /* With anisotropic porosity */
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        cs_math_sym_33_product(porosf[cell_id],
                               coefav[face_id],
                               coefaq[face_id]);

        for (cs_lnum_t isou = 0; isou < 6; isou++)
          coefaq[face_id][isou] *= b_rho[face_id];

        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               b_mass_var[face_id]);

        for (cs_lnum_t isou = 0; isou < 6; isou++)
          b_mass_var[face_id][isou] *= b_rho[face_id];

      });
    }

  }

  /* Velocity flux */
  else {

    /* Without porosity */
    if (porosi == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou];
          b_mass_var[face_id][isou] = c_var[cell_id][isou];
        }
      });
    }
    /* With porosity */
    else if (porosi != nullptr && porosf == nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou]*porosi[cell_id];
          b_mass_var[face_id][isou] = c_var[cell_id][isou]*porosi[cell_id];
        }
      });
    }
    /* With anisotropic porosity */
    else if (porosi != nullptr && porosf != nullptr) {
      ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        cs_math_sym_33_product(porosf[cell_id],
                               coefav[face_id],
                               coefaq[face_id]);

        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               b_mass_var[face_id]);
      });
    }

  }

  /*==========================================================================
    2. Compute mass flux without recontructions
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t w_i = weight[face_id];
      cs_real_t w_j = (1. - weight[face_id]);

      cs_real_t f_mass_var[6];
      for (cs_lnum_t isou = 0; isou < 6; isou++)
        f_mass_var[isou] = w_i * c_mass_var[ii][isou] + w_j * c_mass_var[jj][isou];

      cs_real_t i_face_normal[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        i_face_normal[isou] =   i_face_u_normal[face_id][isou]
                              * i_face_surf[face_id];

      cs_math_sym_33_3_product_add(f_mass_var,
                                   i_face_normal,
                                   i_massflux[face_id]);

    });

    /* Boundary faces */

    ctx_c.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_real_t f_mass_var[6];

      /* var_f = a + b * var_i */
      for (cs_lnum_t isou = 0; isou < 6; isou++)
        f_mass_var[isou] = inc*coefaq[face_id][isou];

      cs_math_66_6_product_add(coefbv[face_id],
                               b_mass_var[face_id],
                               f_mass_var);

      cs_real_t b_face_normal[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        b_face_normal[isou] =   b_face_u_normal[face_id][isou]
                              * b_face_surf[face_id];

      cs_math_sym_33_3_product_add(f_mass_var,
                                   b_face_normal,
                                   b_massflux[face_id]);

    });
  }
  ctx.wait();
  ctx_c.wait();

  /*==========================================================================
    4. Compute mass flux with reconstruction technics if the mesh is
       non orthogonal
    ==========================================================================*/

  if (nswrgu > 1) {

    CS_MALLOC_HD(c_grad_mvar, n_cells_ext, cs_real_63_t, amode);

    /* As coefa has just been modified, face value for gradient
       computation need to be updated. coefb is the same */

    cs_real_6_t *val_f, *val_ip_g;
    CS_MALLOC_HD(val_f, n_b_faces, cs_real_6_t, amode);
    CS_MALLOC_HD(val_ip_g, n_b_faces, cs_real_6_t, amode);

    cs_gradient_boundary_iprime_lsq_strided<6>(ctx,
                                               m,
                                               fvq,
                                               n_b_faces,
                                               nullptr,
                                               halo_type,
                                               -1,
                                               nullptr,
                                               &bc_coeffs_ts_loc,
                                               nullptr, // gweight,
                                               (const cs_real_6_t *)c_mass_var,
                                               (cs_real_6_t *)val_ip_g,
                                               nullptr);

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      for (cs_lnum_t i = 0; i < 6; i++) {
        val_f[face_id][i] = coefaq[face_id][i];
        for (cs_lnum_t j = 0; j < 6; j++) {
          val_f[face_id][i] += coefbv[face_id][j][i]*val_ip_g[face_id][j];
        }
      }
    });

    ctx.wait();

    CS_FREE_HD(val_ip_g);

    /* Computation of c_mass_var gradient
       (tensor gradient, the periodicity has already been treated) */

    cs_gradient_tensor_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgu,
                                    iwarnu,
                                    imligu,
                                    epsrgu,
                                    climgu,
                                    &bc_coeffs_ts_loc,
                                    (const cs_real_6_t *)c_mass_var,
                                    val_f,
                                    c_grad_mvar);

    CS_FREE_HD(val_f);

    /* Mass flow through interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t w_i = weight[face_id];
      cs_real_t w_j = (1. - weight[face_id]);

      cs_real_t f_mass_var[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++)
        /* Non-reconstructed face value */
        f_mass_var[isou] = w_i * c_mass_var[ii][isou] + w_j * c_mass_var[jj][isou]
          /* Reconstruction: face gradient times OF */
          + 0.5*(c_grad_mvar[ii][isou][0] +c_grad_mvar[jj][isou][0])*dofij[face_id][0]
          + 0.5*(c_grad_mvar[ii][isou][1] +c_grad_mvar[jj][isou][1])*dofij[face_id][1]
          + 0.5*(c_grad_mvar[ii][isou][2] +c_grad_mvar[jj][isou][2])*dofij[face_id][2];

      cs_real_t i_face_normal[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        i_face_normal[isou] =   i_face_u_normal[face_id][isou]
                              * i_face_surf[face_id];

      cs_math_sym_33_3_product_add(f_mass_var,
                                   i_face_normal,
                                   i_massflux[face_id]);

    });

    /* Mass flow through boundary faces */
    ctx_c.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t f_mass_var[6];

      /* var_f = a + b * var_I' */
      for (cs_lnum_t isou = 0; isou < 6; isou++)
        f_mass_var[isou] = inc*coefaq[face_id][isou];

      /* Add the reconstruction to get value in I' */
      for (cs_lnum_t jsou = 0; jsou < 6; jsou++)
        b_mass_var[face_id][jsou] += c_grad_mvar[ii][jsou][0]*diipb[face_id][0]
                                   + c_grad_mvar[ii][jsou][1]*diipb[face_id][1]
                                   + c_grad_mvar[ii][jsou][2]*diipb[face_id][2];

      cs_math_66_6_product_add(coefbv[face_id],
                               b_mass_var[face_id],
                               f_mass_var);

      cs_real_t b_face_normal[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        b_face_normal[isou] =   b_face_u_normal[face_id][isou]
                              * b_face_surf[face_id];

      cs_math_sym_33_3_product_add(f_mass_var,
                                   b_face_normal,
                                   b_massflux[face_id]);

    });

  }

  /*==========================================================================
    6. Here, we make sure that the mass flux is null at the boundary faces of
       type symmetry and coupled walls.
    ==========================================================================*/

  if (iflmb0 == 1) {
    /* Force flumab to 0 for velocity */
    cs_lnum_t *b_sym_flag = fvq->b_sym_flag;
    ctx_c.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      if (b_sym_flag[face_id] == 0) {
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          b_massflux[face_id][isou] = 0.;
      }
    });
  }

  ctx.wait();
  ctx_c.wait();

  CS_FREE_HD(c_grad_mvar);
  CS_FREE_HD(c_mass_var);
  CS_FREE_HD(b_mass_var);

  coefaq = nullptr;
  cs_field_bc_coeffs_free_copy(bc_coeffs_ts, &bc_coeffs_ts_loc);
}

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
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  /* Parallel or device dispatch */

  cs_dispatch_context ctx;

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      diverg[cell_id] = 0.;
    });
  }
  else if (init < 0)
    bft_error(__FILE__, __LINE__, 0, _("invalid value of init"));


  /*==========================================================================
    2. Integration on internal faces
    ==========================================================================*/

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    if (ii < n_cells)
      cs_dispatch_sum(&diverg[ii], i_massflux[face_id], i_sum_type);
    if (jj < n_cells)
      cs_dispatch_sum(&diverg[jj],-i_massflux[face_id], i_sum_type);
  });

  /*==========================================================================
    3. Integration on border faces
    ==========================================================================*/

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    cs_lnum_t ii = b_face_cells[face_id];
    cs_dispatch_sum(&diverg[ii], b_massflux[face_id], b_sum_type);
  });
  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integrated mass flux on the cells for a tensor variable.
 *
 * \f[
 * \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_massflux    mass flux vector at interior faces
 * \param[in]     b_massflux    mass flux vector at boundary faces
 * \param[in,out] diverg        mass flux divergence vector
 */
/*----------------------------------------------------------------------------*/

void
cs_tensor_divergence(const cs_mesh_t            *m,
                     int                         init,
                     const cs_real_3_t           i_massflux[],
                     const cs_real_3_t           b_massflux[],
                     cs_real_3_t       *restrict diverg)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  cs_dispatch_context ctx;

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
     ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        diverg[cell_id][isou] = 0.;
      }
    });
  }
  else if (init < 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /*==========================================================================
    2. Integration on internal faces
    ==========================================================================*/

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];
    cs_real_3_t flux_p, flux_m;

    for (cs_lnum_t isou = 0; isou < 3; isou++) {
      flux_p[isou] =  i_massflux[face_id][isou];
      flux_m[isou] = -i_massflux[face_id][isou];
    }

    if (ii < n_cells)
      cs_dispatch_sum<3>(diverg[ii], flux_p, i_sum_type);
    if (jj < n_cells)
      cs_dispatch_sum<3>(diverg[jj], flux_m, i_sum_type);

  });

  /*==========================================================================
    3. Integration on border faces
    ==========================================================================*/

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    cs_lnum_t ii = b_face_cells[face_id];
    cs_dispatch_sum<3>(diverg[ii], b_massflux[face_id], b_sum_type);
  });

  ctx.wait();
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
 *                               of the variable (implicit part)
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
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_t *restrict i_f_face_surf = fvq->i_face_surf;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;

  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Parallel or device dispatch */

  cs_dispatch_context ctx, ctx_c;
  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  ctx_c.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_c.use_gpu())
    ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  /*Additional terms due to porosity */

  cs_field_t *f_i_poro_duq_0 = cs_field_by_name_try("i_poro_duq_0");

  cs_real_t *i_poro_duq_0;
  cs_real_t *i_poro_duq_1;
  cs_real_t *b_poro_duq;
  cs_real_t *_f_ext = nullptr;

  int is_p = 0; /* Is porous ? */
  if (f_i_poro_duq_0 != nullptr) {

    is_p = 1;
    i_poro_duq_0 = f_i_poro_duq_0->val;
    i_poro_duq_1 = cs_field_by_name("i_poro_duq_1")->val;
    b_poro_duq = cs_field_by_name("b_poro_duq")->val;

  }
  else {

    CS_MALLOC_HD(_f_ext, 1, cs_real_t, amode);
    _f_ext[0] = 0.;
    i_poro_duq_0 = _f_ext;
    i_poro_duq_1 = _f_ext;
    b_poro_duq = _f_ext;

  }

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      i_massflux[face_id] = 0.;
    });
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      b_massflux[face_id] = 0.;
    });
  }
  else if (init != 0)
    bft_error(__FILE__, __LINE__, 0, _("invalid value of init"));

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Mass flow through interior faces */
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];
      cs_lnum_t _p = is_p*face_id;

      cs_real_2_t poro = { i_poro_duq_0[_p],
                           i_poro_duq_1[_p] };

      i_massflux[face_id] +=
        i_visc[face_id]*(  (i_face_cog[face_id][0]-cell_cen[ii][0])*frcxt[ii][0]
                         + (i_face_cog[face_id][1]-cell_cen[ii][1])*frcxt[ii][1]
                         + (i_face_cog[face_id][2]-cell_cen[ii][2])*frcxt[ii][2]
                         + poro[0]
                         - (i_face_cog[face_id][0]-cell_cen[jj][0])*frcxt[jj][0]
                         - (i_face_cog[face_id][1]-cell_cen[jj][1])*frcxt[jj][1]
                         - (i_face_cog[face_id][2]-cell_cen[jj][2])*frcxt[jj][2]
                         - poro[1]
                        );

    });

    /* Mass flux through boundary faces */
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = b_face_cells[face_id];
      cs_lnum_t _p = is_p*face_id;

      const cs_nreal_t *normal = b_face_u_normal[face_id];
      const cs_real_t poro = b_poro_duq[_p];

      b_massflux[face_id] += b_visc[face_id] * cofbfp[face_id] *
        ( cs_math_3_dot_product(frcxt[ii], normal) * b_dist[face_id] + poro );

    });

  /*==========================================================================
    3. Update mass flux with reconstruction method
    ==========================================================================*/

  }
  else {

    /* Mass flux through interior faces */
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];
      cs_lnum_t _p = is_p*face_id;

      cs_real_2_t poro = { i_poro_duq_0[_p],
                           i_poro_duq_1[_p] };

      double surfn = i_f_face_surf[face_id];

      i_massflux[face_id] += i_visc[face_id]*
        ( ( i_face_cog[face_id][0] - cell_cen[ii][0] ) * frcxt[ii][0] +
          ( i_face_cog[face_id][1] - cell_cen[ii][1] ) * frcxt[ii][1] +
          ( i_face_cog[face_id][2] - cell_cen[ii][2] ) * frcxt[ii][2] +
          poro[0] -
          ( i_face_cog[face_id][0] - cell_cen[jj][0] ) * frcxt[jj][0] -
          ( i_face_cog[face_id][1] - cell_cen[jj][1] ) * frcxt[jj][1] -
          ( i_face_cog[face_id][2] - cell_cen[jj][2] ) * frcxt[jj][2] -
          poro[1]
        )
       + surfn/i_dist[face_id]*0.5
        *( (djjpf[face_id][0] - diipf[face_id][0])*
           (viselx[ii]*frcxt[ii][0] + viselx[jj]*frcxt[jj][0])
         + (djjpf[face_id][1] - diipf[face_id][1])*
           (visely[ii]*frcxt[ii][1] + visely[jj]*frcxt[jj][1])
         + (djjpf[face_id][2] - diipf[face_id][2])*
           (viselz[ii]*frcxt[ii][2] + viselz[jj]*frcxt[jj][2])
         );

    }); /* Loop on interior faces */

    /* Mass flux through boundary faces */
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      cs_lnum_t _p = is_p*face_id;

      const cs_nreal_t *normal = b_face_u_normal[face_id];
      const cs_real_t poro = b_poro_duq[_p];

      b_massflux[face_id] += b_visc[face_id] * cofbfp[face_id] *
        (cs_math_3_dot_product(frcxt[ii], normal) * b_dist[face_id] + poro);

    }); /* Loop on border faces */

  }

  ctx.wait();
  ctx_c.wait();

  CS_FREE_HD(_f_ext);
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
 *                               of the variable (implicit part)
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
                              cs_real_6_t               viscel[],
                              const cs_real_2_t         weighf[],
                              cs_real_t       *restrict i_massflux,
                              cs_real_t       *restrict b_massflux)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /*==========================================================================*/

  /* Parallel or device dispatch */

  cs_dispatch_context ctx, ctx_c;
  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  ctx_c.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_c.use_gpu())
    ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  const bool on_device = ctx.use_gpu();

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      i_massflux[face_id] = 0.;
    });
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      b_massflux[face_id] = 0.;
    });
  }
  else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  cs_real_6_t *w2 = nullptr;

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* ---> Contribution from interior faces */

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

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

    });

    /* ---> Contribution from boundary faces */

    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      const cs_nreal_t *normal = b_face_u_normal[face_id];

      b_massflux[face_id] += b_visc[face_id] * b_dist[face_id]
                            * cofbfp[face_id]
                            * cs_math_3_dot_product(frcxt[ii], normal);

    });

    /*========================================================================
      3. Update mass flux with reconstruction technique
      ========================================================================*/

  }
  else {

    cs_real_6_t *viscce = nullptr;

    /* Without porosity */
    if (porosi == nullptr) {
      viscce = viscel;

      /* With porosity */
    }
    else if (porosi != nullptr && porosf == nullptr) {
      CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, amode);
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
        }
      });
      viscce = w2;

      /* With tensorial porosity */
    } else if (porosi != nullptr && porosf != nullptr) {
      CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, amode);
      ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        cs_math_sym_33_product(porosf[cell_id],
                               viscel[cell_id],
                               w2[cell_id]);
      });
      viscce = w2;
    }

    ctx_c.wait();

    /* ---> Periodicity and parallelism treatment of symmetric tensors */

    if (halo != nullptr)
      cs_halo_sync_r(halo, on_device, viscce);

    /* ---> Contribution from interior faces */

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      /* Local variables */
      cs_real_t  diippf[3], djjppf[3];
      cs_real_t visci[3][3], viscj[3][3];
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      /* Recompute II' and JJ' at this level */

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      double fikdvi = weighf[face_id][0] * i_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =  i_face_cog[face_id][i]-cell_cen[ii][i]
                   - fikdvi*(  visci[0][i]*i_face_u_normal[face_id][0]
                             + visci[1][i]*i_face_u_normal[face_id][1]
                             + visci[2][i]*i_face_u_normal[face_id][2]);
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      double fjkdvi = weighf[face_id][1] * i_face_surf[face_id];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi*(  viscj[0][i]*i_face_u_normal[face_id][0]
                              + viscj[1][i]*i_face_u_normal[face_id][1]
                              + viscj[2][i]*i_face_u_normal[face_id][2]);
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

    });

    /* ---> Contribution from boundary faces */

    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      const cs_nreal_t *normal = b_face_u_normal[face_id];

      /* FIXME: wrong if dirichlet and viscce is really a tensor */
      b_massflux[face_id] += b_visc[face_id] * b_dist[face_id]
                            * cofbfp[face_id]
                            * cs_math_3_dot_product(frcxt[ii], normal);

    });

  }

  ctx.wait();
  ctx_c.wait();

  CS_FREE_HD(w2);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
