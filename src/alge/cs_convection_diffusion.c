/*============================================================================
 * Convection-diffusion operators.
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
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_convection_diffusion.h"

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
 * Return pointer to slope test indicator field values if active.
 *
 * parameters:
 *   f_id        <-- field id (or -1)
 *   var_cal_opt <-- variable calculation options
 *
 * return:
 *   pointer to local values array, or NULL;
 *----------------------------------------------------------------------------*/

static cs_real_t  *_get_v_slope_test(int                       f_id,
                                     const cs_var_cal_opt_t    var_cal_opt)
{
  const int iconvp = var_cal_opt.iconv;
  const int isstpp = var_cal_opt.isstpc;
  const double blencp = var_cal_opt.blencv;

  cs_real_t  *v_slope_test = NULL;

  if (f_id > -1 && iconvp > 0 && blencp > 0. && isstpp == 0) {

    static int _k_slope_test_f_id = -1;

    cs_field_t *f = cs_field_by_id(f_id);

    int f_track_slope_test_id = -1;

    if (_k_slope_test_f_id < 0)
      _k_slope_test_f_id = cs_field_key_id_try("slope_test_upwind_id");
    if (_k_slope_test_f_id > -1 && isstpp == 0)
      f_track_slope_test_id = cs_field_get_key_int(f, _k_slope_test_f_id);

    if (f_track_slope_test_id > -1)
    v_slope_test = (cs_field_by_id(f_track_slope_test_id))->val;

    if (v_slope_test != NULL) {
      const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
        v_slope_test[cell_id] = 0.;
    }

  }

  return v_slope_test;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc2, BILSC2)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   icvflb,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   iccocg,
 const cs_int_t          *const   ifaccp,
 cs_real_t                        pvar[],
 const cs_real_t                  pvara[],
 const cs_int_t                   bc_type[],
 const cs_int_t                   icvfli[],
 const cs_real_t                  coefap[],
 const cs_real_t                  coefbp[],
 const cs_real_t                  cofafp[],
 const cs_real_t                  cofbfp[],
 const cs_real_t                  i_massflux[],
 const cs_real_t                  b_massflux[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 cs_real_t                        rhs[]
)
{
  cs_convection_diffusion_scalar(*idtvar,
                                 *f_id,
                                 *var_cal_opt,
                                 *icvflb,
                                 *inc,
                                 *iccocg,
                                 *ifaccp,
                                 pvar,
                                 pvara,
                                 bc_type,
                                 icvfli,
                                 coefap,
                                 coefbp,
                                 cofafp,
                                 cofbfp,
                                 i_massflux,
                                 b_massflux,
                                 i_visc,
                                 b_visc,
                                 rhs);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc4, BILSC4)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   icvflb,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   ifaccp,
 const cs_int_t          *const   ivisep,
 cs_real_3_t                      pvar[],
 const cs_real_3_t                pvara[],
 const cs_int_t                   bc_type[],
 const cs_int_t                   icvfli[],
 const cs_real_3_t                coefav[],
 const cs_real_33_t               coefbv[],
 const cs_real_3_t                cofafv[],
 const cs_real_33_t               cofbfv[],
 const cs_real_t                  i_massflux[],
 const cs_real_t                  b_massflux[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 const cs_real_t                  secvif[],
 cs_real_3_t                      rhs[]
)
{
  cs_convection_diffusion_vector(*idtvar,
                                 *f_id,
                                 *var_cal_opt,
                                 *icvflb,
                                 *inc,
                                 *ifaccp,
                                 *ivisep,
                                 pvar,
                                 pvara,
                                 bc_type,
                                 icvfli,
                                 coefav,
                                 coefbv,
                                 cofafv,
                                 cofbfv,
                                 i_massflux,
                                 b_massflux,
                                 i_visc,
                                 b_visc,
                                 secvif,
                                 rhs);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_thermal
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsct, BILSCT)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   iccocg,
 const cs_int_t          *const   ifaccp,
 cs_real_t                        pvar[],
 const cs_real_t                  pvara[],
 const cs_int_t                   bc_type[],
 const cs_real_t                  coefap[],
 const cs_real_t                  coefbp[],
 const cs_real_t                  cofafp[],
 const cs_real_t                  cofbfp[],
 const cs_real_t                  i_massflux[],
 const cs_real_t                  b_massflux[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 const cs_real_t                  xcpp[],
 cs_real_t                        rhs[]
)
{
  cs_convection_diffusion_thermal(*idtvar,
                                  *f_id,
                                  *var_cal_opt,
                                  *inc,
                                  *iccocg,
                                  *ifaccp,
                                  pvar,
                                  pvara,
                                  bc_type,
                                  coefap,
                                  coefbp,
                                  cofafp,
                                  cofbfp,
                                  i_massflux,
                                  b_massflux,
                                  i_visc,
                                  b_visc,
                                  xcpp,
                                  rhs);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (diften, DIFTEN)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   iccocg,
 cs_real_t                        pvar[],
 const cs_real_t                  pvara[],
 const cs_real_t                  coefap[],
 const cs_real_t                  coefbp[],
 const cs_real_t                  cofafp[],
 const cs_real_t                  cofbfp[],
 const cs_real_t                  i_visc[],
 const cs_real_t                  b_visc[],
 cs_real_6_t                      viscel[],
 const cs_real_2_t                weighf[],
 const cs_real_t                  weighb[],
 cs_real_t                        rhs[]
)
{
  cs_anisotropic_diffusion_scalar(*idtvar,
                                  *f_id,
                                  *var_cal_opt,
                                  *inc,
                                  *iccocg,
                                  pvar,
                                  pvara,
                                  coefap,
                                  coefbp,
                                  cofafp,
                                  cofbfp,
                                  i_visc,
                                  b_visc,
                                  viscel,
                                  weighf,
                                  weighb,
                                  rhs);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (diftnv, DIFTNV)
(
 const cs_int_t          *const   idtvar,
 const cs_int_t          *const   f_id,
 const cs_var_cal_opt_t  *const   var_cal_opt,
 const cs_int_t          *const   inc,
 const cs_int_t          *const   ifaccp,
 const cs_int_t          *const   ivisep,
 cs_real_3_t                      pvar[],
 const cs_real_3_t                pvara[],
 const cs_int_t                   bc_type[],
 const cs_real_3_t                coefav[],
 const cs_real_33_t               coefbv[],
 const cs_real_3_t                cofafv[],
 const cs_real_33_t               cofbfv[],
 const cs_real_33_t               i_visc[],
 const cs_real_t                  b_visc[],
 const cs_real_t                  secvif[],
 cs_real_3_t                      rhs[]
)
{
  cs_anisotropic_diffusion_vector(*idtvar,
                                  *f_id,
                                  *var_cal_opt,
                                  *inc,
                                  *ifaccp,
                                  *ivisep,
                                  pvar,
                                  pvara,
                                  bc_type,
                                  coefav,
                                  coefbv,
                                  cofafv,
                                  cofbfv,
                                  i_visc,
                                  b_visc,
                                  secvif,
                                  rhs);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmas, ITRMAS)
(
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
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
 cs_real_t                b_massflux[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_face_diffusion_potential(*f_id,
                              m,
                              fvq,
                              *init,
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
 * Wrapper to cs_face_anisotropic_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmav, ITRMAV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
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
 cs_real_6_t              viscel[],
 const cs_real_2_t        weighf[],
 const cs_real_t          weighb[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_face_anisotropic_diffusion_potential(m,
                                          fvq,
                                          *init,
                                          *inc,
                                          *imrgra,
                                          *iccocg,
                                          *nswrgp,
                                          *imligp,
                                          *ircflp,
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
                                          viscel,
                                          weighf,
                                          weighb,
                                          i_massflux,
                                          b_massflux);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrp, ITRGRP)
(
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
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
 cs_real_t                diverg[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_diffusion_potential(*f_id,
                         m,
                         fvq,
                         *init,
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
 * Wrapper to cs_anisotropic_diffusion_potential
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrv, ITRGRV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
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
 cs_real_6_t              viscel[],
 const cs_real_2_t        weighf[],
 const cs_real_t          weighb[],
 cs_real_t                diverg[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_anisotropic_diffusion_potential(m,
                                     fvq,
                                     *init,
                                     *inc,
                                     *imrgra,
                                     *iccocg,
                                     *nswrgp,
                                     *imligp,
                                     *ircflp,
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
                                     viscel,
                                     weighf,
                                     weighb,
                                     diverg);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before calling bilsc2!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                   (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_scalar(int                       idtvar,
                               int                       f_id,
                               const cs_var_cal_opt_t    var_cal_opt,
                               int                       icvflb,
                               int                       inc,
                               int                       iccocg,
                               int                       ifaccp,
                               cs_real_t       *restrict pvar,
                               const cs_real_t *restrict pvara,
                               const cs_int_t            bc_type[],
                               const cs_int_t            icvfli[],
                               const cs_real_t           coefap[],
                               const cs_real_t           coefbp[],
                               const cs_real_t           cofafp[],
                               const cs_real_t           cofbfp[],
                               const cs_real_t           i_massflux[],
                               const cs_real_t           b_massflux[],
                               const cs_real_t           i_visc[],
                               const cs_real_t           b_visc[],
                               cs_real_t       *restrict rhs)
{
  const int iconvp = var_cal_opt.iconv;
  const int idiffp = var_cal_opt.idiff;
  const int nswrgp = var_cal_opt.nswrgr;
  const int imrgra = var_cal_opt.imrgra;
  const int imligp = var_cal_opt.imligr;
  const int ircflp = var_cal_opt.ircflu;
  const int ischcp = var_cal_opt.ischcv;
  const int isstpp = var_cal_opt.isstpc;
  const int iwarnp = var_cal_opt.iwarni;
  const double blencp = var_cal_opt.blencv;
  const double epsrgp = var_cal_opt.epsrgr;
  const double climgp = var_cal_opt.climgr;
  const double extrap = var_cal_opt.extrag;
  const double relaxp = var_cal_opt.relaxv;
  const double thetap = var_cal_opt.thetav;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* Local variables */

  char var_name[32];

  cs_gnum_t n_upwind;
  cs_lnum_t cell_id, face_id, ii, jj;
  int iupwin;
  int g_id, t_id;
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  double pfac,pfacd,flui,fluj,flux,fluxi,fluxj;
  double difx,dify,difz,djfx,djfy,djfz;
  double pi, pj, pia, pja;
  double pif,pjf,pip,pjp,pir,pjr,pipr,pjpr;
  double pifri,pifrj,pjfri,pjfrj;
  double testi,testj,testij;
  double dpxf,dpyf,dpzf;
  double dcc, ddi, ddj, tesqck;
  double dijpfx, dijpfy, dijpfz;
  double diipfx, diipfy, diipfz;
  double djjpfx, djjpfy, djjpfz;
  double diipbx, diipby, diipbz;
  double pnd, distf, srfan;
  double pfac1, pfac2, pfac3, unsvol;
  cs_real_t *coface, *cofbce;

  cs_real_3_t *grad;
  cs_real_3_t *grdpa;
  cs_field_t *f;

  cs_real_t *gweight = NULL;

  cs_real_t  *v_slope_test = _get_v_slope_test(f_id,  var_cal_opt);

  /* 1. Initialization */

  /* Allocate work arrays */

  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(grdpa, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    cs_gradient_perio_init_rij(f, &tr_dim, grad);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    } else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  iupwin = (blencp > 0.) ? 0 : 1;

  /* 2. Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test. */

  if (  (idiffp!=0 && ircflp==1)
     || (  iconvp!=0 && iupwin==0
        && (ischcp==0 || ircflp==1 || isstpp==0))) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && var_cal_opt.iwgrec == 1) {
        if (var_cal_opt.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
          }
        }
      }
    }

    cs_gradient_scalar(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       recompute_cocg,
                       nswrgp,
                       tr_dim,
                       0, /* hyd_p_flag */
                       iwarnp,
                       imligp,
                       epsrgp,
                       extrap,
                       climgp,
                       NULL, /* f_ext exterior force */
                       coefap,
                       coefbp,
                       pvar,
                       gweight, /* Weighted gradient */
                       grad);

  } else {
#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    }
  }

  /* 2.1 Compute the slope test gradient */

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    grdpa[cell_id][0] = 0.;
    grdpa[cell_id][1] = 0.;
    grdpa[cell_id][2] = 0.;
  }

  if (iconvp>0 && iupwin==0 && isstpp==0) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, difx, dify, difz, \
                                      djfx, djfy, djfz, pif, pjf, pfac,  \
                                      pfac1, pfac2, pfac3)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          difx = i_face_cog[face_id][0] - cell_cen[ii][0];
          dify = i_face_cog[face_id][1] - cell_cen[ii][1];
          difz = i_face_cog[face_id][2] - cell_cen[ii][2];
          djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
          djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
          djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

          pif = pvar[ii] + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
          pjf = pvar[jj] + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

          pfac = pjf;
          if (i_massflux[face_id] > 0.) pfac = pif;

          pfac1 = pfac*i_face_normal[face_id][0];
          pfac2 = pfac*i_face_normal[face_id][1];
          pfac3 = pfac*i_face_normal[face_id][2];

          grdpa[ii][0] = grdpa[ii][0] + pfac1;
          grdpa[ii][1] = grdpa[ii][1] + pfac2;
          grdpa[ii][2] = grdpa[ii][2] + pfac3;

          grdpa[jj][0] = grdpa[jj][0] - pfac1;
          grdpa[jj][1] = grdpa[jj][1] - pfac2;
          grdpa[jj][2] = grdpa[jj][2] - pfac3;

        }
      }
    }

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz, pfac) \
                              if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          diipbx = diipb[face_id][0];
          diipby = diipb[face_id][1];
          diipbz = diipb[face_id][2];
          pfac =   inc*coefap[face_id]
                 + coefbp[face_id] * (pvar[ii] + diipbx*grad[ii][0]
                                               + diipby*grad[ii][1]
                                               + diipbz*grad[ii][2]);
          grdpa[ii][0] = grdpa[ii][0] + pfac*b_face_normal[face_id][0];
          grdpa[ii][1] = grdpa[ii][1] + pfac*b_face_normal[face_id][1];
          grdpa[ii][2] = grdpa[ii][2] + pfac*b_face_normal[face_id][2];

        }
      }
    }

#   pragma omp parallel for private(unsvol)
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      unsvol = 1./cell_vol[cell_id];
      grdpa[cell_id][0] = grdpa[cell_id][0]*unsvol;
      grdpa[cell_id][1] = grdpa[cell_id][1]*unsvol;
      grdpa[cell_id][2] = grdpa[cell_id][2]*unsvol;
    }

    /* Synchronization for parallelism or periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)grdpa, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)grdpa, 3);

      // Gradient periodicity of rotation for Reynolds stress components
      if (m->have_rotation_perio > 0 && f_id != -1)
        cs_gradient_perio_process_rij(&f_id, grdpa);
    }

  }


  /* ======================================================================
    ---> Contribution from interior faces
    ======================================================================*/

  n_upwind = 0;

  if (n_cells_ext>n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      rhs[cell_id] = 0.;
    }
  }

  /* --> Pure upwind flux
    =====================*/

  if (iupwin==1) {

    /* Steady */
    if (idtvar<0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz,  \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,      \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf, pip,    \
                                        pjp, pipr, pjpr, flui, fluj, pif, pjf,    \
                                        fluxi, fluxj, pi, pj, pir, pjr, pia, pja) \
                            reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];
            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells) {
              n_upwind++;
            }

            pi = pvar[ii];
            pj = pvar[jj];
            pia = pvara[ii];
            pja = pvara[jj];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            /* reconstruction only if IRCFLP = 1 */
            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            pir = pi/relaxp - (1.-relaxp)/relaxp * pia;
            pjr = pj/relaxp - (1.-relaxp)/relaxp * pja;

            pipr = pir
                 + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjpr = pjr
                 + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            pif  = pi;
            pjf  = pj;

            fluxi = iconvp*(flui*pir + fluj*pjf - i_massflux[face_id]*pi)
                  + idiffp*i_visc[face_id]*(pipr - pjp);
            fluxj = iconvp*(flui*pif + fluj*pjr - i_massflux[face_id]*pj)
                  + idiffp*i_visc[face_id]*(pip - pjpr);

            rhs[ii] = rhs[ii] - fluxi;
            rhs[jj] = rhs[jj] + fluxj;

          }
        }
      }

    /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz,             \
                                        djjpfx, djjpfy, djjpfz, dpxf, dpyf, dpzf,\
                                        pip, pjp, flui, fluj,                    \
                                        pif, pjf, flux, pi, pj)                  \
                            reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];
            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells) {
              n_upwind++;
            }

            pi = pvar[ii];
            pj = pvar[jj];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */

            diipfx = i_face_cog[face_id][0] - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy = i_face_cog[face_id][1] - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz = i_face_cog[face_id][2] - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            pif = pi;
            pjf = pj;

            flux =   iconvp*(flui*pif +fluj*pjf)
                   + idiffp*i_visc[face_id]*(pip -pjp);

            rhs[ii] = rhs[ii] - thetap *(flux - iconvp*i_massflux[face_id]*pi);
            rhs[jj] = rhs[jj] + thetap *(flux - iconvp*i_massflux[face_id]*pj);

          }
        }
      }

    }

  /* --> Flux with no slope test
    ============================*/

  } else if (isstpp==1) {

    if (ischcp<0 || ischcp>1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    /* Steady */
    if (idtvar<0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,     \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf, pip,   \
                                        pjp, pipr, pjpr, flui, fluj, pir, pjr,   \
                                        pifri, pjfri, pifrj, pjfrj,              \
                                        difx, dify, difz, djfx, djfy, djfz,      \
                                        fluxi, fluxj, pi, pj, pia, pja)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            pi = pvar[ii];
            pj = pvar[jj];
            pia = pvara[ii];
            pja = pvara[jj];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            pir = pi/relaxp - (1. - relaxp)/relaxp*pia;
            pjr = pj/relaxp - (1. - relaxp)/relaxp*pja;

            pipr = pir
                 + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjpr = pjr
                 + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));


            /* Centered
              --------*/

            if (ischcp==1) {

              pifri = pnd*pipr +(1.-pnd)*pjp;
              pjfri = pifri;
              pifrj = pnd*pip  +(1.-pnd)*pjpr;
              pjfrj = pifrj;

            /* Second order
              ------------*/

            } else { /* if (ischcp==0)  */

              difx = i_face_cog[face_id][0] - cell_cen[ii][0];
              dify = i_face_cog[face_id][1] - cell_cen[ii][1];
              difz = i_face_cog[face_id][2] - cell_cen[ii][2];
              djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
              djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
              djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

              /* leave reconstruction of PIF and PJF even if IRCFLP=0
                otherwise, it is the same as using upwind */
              pifri = pir + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
              pifrj = pi  + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
              pjfrj = pjr + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];
              pjfri = pj  + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

            }

            /* Blending
              --------*/

            pifri = blencp*pifri+(1.-blencp)*pir;
            pifrj = blencp*pifrj+(1.-blencp)*pi;
            pjfri = blencp*pjfri+(1.-blencp)*pj;
            pjfrj = blencp*pjfrj+(1.-blencp)*pjr;

            /* Flux
              ----*/

            fluxi =   iconvp*(flui*pifri + fluj*pjfri - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pipr -pjp);
            fluxj =   iconvp*(flui*pifrj + fluj*pjfrj - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip -pjpr);

            /* Assembly
              --------*/

            rhs[ii] = rhs[ii] - fluxi;
            rhs[jj] = rhs[jj] + fluxj;

          }
        }
      }

    /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,     \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf,        \
                                        pip, pjp, flui, fluj, pif, pjf, difx,    \
                                        dify, difz, djfx, djfy, djfz, flux,      \
                                        pi, pj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            pi = pvar[ii];
            pj = pvar[jj];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            /* Centered
              --------*/

            if (ischcp==1) {

              pif = pnd*pip +(1.-pnd)*pjp;
              pjf = pif;

            /* Second order
              ------------*/

            } else { /* if (ischcp==0)  */

              difx = i_face_cog[face_id][0] - cell_cen[ii][0];
              dify = i_face_cog[face_id][1] - cell_cen[ii][1];
              difz = i_face_cog[face_id][2] - cell_cen[ii][2];
              djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
              djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
              djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

              /* leave reconstruction of PIF and PJF even if IRCFLP=0
                otherwise, it is the same as using upwind */
              pif = pi + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
              pjf = pj + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

            }

            /* Blending
              --------*/

            pif = blencp*pif+(1.-blencp)*pi;
            pjf = blencp*pjf+(1.-blencp)*pj;

            /* Flux
              ----*/

            flux =   iconvp*(flui*pif +fluj*pjf)
                   + idiffp*i_visc[face_id]*(pip -pjp);

            /* Assembly
              --------*/

            rhs[ii] = rhs[ii] - thetap *(flux - iconvp*i_massflux[face_id]*pi);
            rhs[jj] = rhs[jj] + thetap *(flux - iconvp*i_massflux[face_id]*pj);

          }
        }
      }

    }

  /* --> Flux with slope test
    =========================*/

  } else {

    if (ischcp<0 || ischcp>1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    /* Steady */
    if (idtvar<0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, distf, srfan, diipfx, diipfy,       \
                                        diipfz, djjpfx, djjpfy, djjpfz, dpxf,    \
                                        dpyf, dpzf, pip, pjp, pipr, pjpr, flui,  \
                                        fluj, pir, pjr, testi, testj, testij,    \
                                        dcc, ddi, ddj, tesqck, pifri, pjfri,     \
                                        pifrj, pjfrj, difx, dify, difz, djfx,    \
                                        djfy, djfz, fluxi, fluxj, pi, pj,        \
                                        pia, pja)                                \
                                reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd    = weight[face_id];
            distf  = i_dist[face_id];
            srfan  = i_face_surf[face_id];

            pi = pvar[ii];
            pj = pvar[jj];
            pia = pvara[ii];
            pja = pvara[jj];

            /* Recompute II' and JJ' at this level */
            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            pir = pi/relaxp - (1. - relaxp)/relaxp*pia;
            pjr = pj/relaxp - (1. - relaxp)/relaxp*pja;

            pipr = pir
                 + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjpr = pjr
                 + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));


            /* Slope test
              ----------*/

            testi =   grdpa[ii][0]*i_face_normal[face_id][0]
                    + grdpa[ii][1]*i_face_normal[face_id][1]
                    + grdpa[ii][2]*i_face_normal[face_id][2];
            testj =   grdpa[jj][0]*i_face_normal[face_id][0]
                    + grdpa[jj][1]*i_face_normal[face_id][1]
                    + grdpa[jj][2]*i_face_normal[face_id][2];
            testij=   grdpa[ii][0]*grdpa[jj][0]
                    + grdpa[ii][1]*grdpa[jj][1]
                    + grdpa[ii][2]*grdpa[jj][2];

            if (i_massflux[face_id]>0.) {
              dcc =   grad[ii][0]*i_face_normal[face_id][0]
                    + grad[ii][1]*i_face_normal[face_id][1]
                    + grad[ii][2]*i_face_normal[face_id][2];
              ddi = testi;
              ddj = (pj-pi)/distf *srfan;
            } else {
              dcc =   grad[jj][0]*i_face_normal[face_id][0]
                    + grad[jj][1]*i_face_normal[face_id][1]
                    + grad[jj][2]*i_face_normal[face_id][2];
              ddi = (pj-pi)/distf *srfan;
              ddj = testj;
            }
            tesqck = pow(dcc, 2.) - pow(ddi-ddj, 2.);

            /* Upwind
              ------*/

            if (tesqck<=0. || testij<=0.) {

              pifri = pir;
              pifrj = pi;
              pjfri = pj;
              pjfrj = pjr;
              /* in parallel, face will be counted by one and only one rank */
              if (ii < n_cells)
                n_upwind++;
              if (v_slope_test != NULL) {
                v_slope_test[ii] += fabs(i_massflux[face_id]) / cell_vol[ii];
                v_slope_test[jj] += fabs(i_massflux[face_id]) / cell_vol[jj];
              }

            } else {

              /* Centered
                --------*/

              if (ischcp==1) {

                pifri = pnd*pipr +(1.-pnd)*pjp;
                pjfri = pifri;
                pifrj = pnd*pip  +(1.-pnd)*pjpr;
                pjfrj = pifrj;

              /* Second order
                ------------*/

              } else { /* if (ischcp==0)  */

                difx = i_face_cog[face_id][0] - cell_cen[ii][0];
                dify = i_face_cog[face_id][1] - cell_cen[ii][1];
                difz = i_face_cog[face_id][2] - cell_cen[ii][2];
                djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
                djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
                djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

                /* leave reconstruction of PIF and PJF even if IRCFLP=0
                  otherwise, it is the same as using upwind */
                pifri = pir + difx*grad[ii][0]
                            + dify*grad[ii][1]
                            + difz*grad[ii][2];
                pifrj = pi  + difx*grad[ii][0]
                            + dify*grad[ii][1]
                            + difz*grad[ii][2];
                pjfrj = pjr + djfx*grad[jj][0]
                            + djfy*grad[jj][1]
                            + djfz*grad[jj][2];
                pjfri = pj  + djfx*grad[jj][0]
                            + djfy*grad[jj][1]
                            + djfz*grad[jj][2];

              }

            }

            /* Blending
              --------*/

            pifri = blencp*pifri+(1.-blencp)*pir;
            pifrj = blencp*pifrj+(1.-blencp)*pi;
            pjfri = blencp*pjfri+(1.-blencp)*pj;
            pjfrj = blencp*pjfrj+(1.-blencp)*pjr;

            /* Flux
              ----*/

            fluxi =   iconvp*(flui*pifri + fluj*pjfri - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pipr -pjp);
            fluxj =   iconvp*(flui*pifrj + fluj*pjfrj - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip -pjpr);

            /* Assembly
              --------*/

            rhs[ii] = rhs[ii] - fluxi;
            rhs[jj] = rhs[jj] + fluxj;

          }
        }
      }

    /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, pnd, \
                                        distf, srfan, diipfx, diipfy, diipfz, djjpfx, \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf, pip, pjp,   \
                                        flui, fluj, testi, testj, testij, dcc, ddi,   \
                                        ddj, tesqck, pif, pjf, difx, dify, difz,      \
                                        djfx, djfy, djfz, flux, pi, pj)               \
                            reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd    = weight[face_id];
            distf  = i_dist[face_id];
            srfan  = i_face_surf[face_id];

            pi = pvar[ii];
            pj = pvar[jj];

            /* Recompute II' and JJ' at this level */

            diipfx = i_face_cog[face_id][0] - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy = i_face_cog[face_id][1] - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz = i_face_cog[face_id][2] - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            /* Slope test
              ----------*/

            testi =   grdpa[ii][0]*i_face_normal[face_id][0]
                    + grdpa[ii][1]*i_face_normal[face_id][1]
                    + grdpa[ii][2]*i_face_normal[face_id][2];
            testj =   grdpa[jj][0]*i_face_normal[face_id][0]
                    + grdpa[jj][1]*i_face_normal[face_id][1]
                    + grdpa[jj][2]*i_face_normal[face_id][2];
            testij =   grdpa[ii][0]*grdpa[jj][0]
                     + grdpa[ii][1]*grdpa[jj][1]
                     + grdpa[ii][2]*grdpa[jj][2];

            if (i_massflux[face_id]>0.) {
              dcc =   grad[ii][0]*i_face_normal[face_id][0]
                    + grad[ii][1]*i_face_normal[face_id][1]
                    + grad[ii][2]*i_face_normal[face_id][2];
              ddi = testi;
              ddj = (pj-pi)/distf *srfan;
            } else {
              dcc =   grad[jj][0]*i_face_normal[face_id][0]
                    + grad[jj][1]*i_face_normal[face_id][1]
                    + grad[jj][2]*i_face_normal[face_id][2];
              ddi = (pj-pi)/distf *srfan;
              ddj = testj;
            }
            tesqck = pow(dcc, 2.) - pow(ddi-ddj, 2.);

            /* Upwind
              ------*/

            if (tesqck<=0. || testij<=0.) {

              pif = pi;
              pjf = pj;
              /* in parallel, face will be counted by one and only one rank */
              if (ii < n_cells)
                n_upwind++;
              if (v_slope_test != NULL) {
                v_slope_test[ii] += fabs(i_massflux[face_id]) / cell_vol[ii];
                v_slope_test[jj] += fabs(i_massflux[face_id]) / cell_vol[jj];
              }

            } else {

              /* Centered
                --------*/

              if (ischcp==1) {

                pif = pnd*pip +(1.-pnd)*pjp;
                pjf = pif;

              /* Second order
                ------------*/

              } else { /* if (ischcp==0)  */

                difx = i_face_cog[face_id][0] - cell_cen[ii][0];
                dify = i_face_cog[face_id][1] - cell_cen[ii][1];
                difz = i_face_cog[face_id][2] - cell_cen[ii][2];
                djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
                djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
                djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

                /* leave reconstruction of PIF and PJF even if IRCFLP=0
                  otherwise, it is the same as using upwind */
                pif = pi + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
                pjf = pj + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

              }

            }

            /* Blending
              --------*/

            pif = blencp*pif+(1.-blencp)*pi;
            pjf = blencp*pjf+(1.-blencp)*pj;

            /* Flux
              ----*/

            flux =   iconvp*(flui*pif +fluj*pjf)
                   + idiffp*i_visc[face_id]*(pip -pjp);

            /* Assembly
              --------*/

            rhs[ii] = rhs[ii] - thetap *(flux - iconvp*i_massflux[face_id]*pi);
            rhs[jj] = rhs[jj] + thetap *(flux - iconvp*i_massflux[face_id]*pj);

          }
        }
      }

    } /* idtvar */

  } /* iupwin */


  if (iwarnp >= 2) {

    /* Sum number of clippings */
    cs_parall_counter(&n_upwind, 1);

  bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces \n"),
             var_name, (unsigned long long)n_upwind,
             (unsigned long long)m->n_g_i_faces);
  }

  /* ======================================================================
    ---> Contribution from boundary faces
    ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

    /* Steady */
    if (idtvar<0) {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz, \
                                        flui, fluj, pir, pipr, pfac, pfacd,  \
                                        flux, pi, pia)                       \
                            if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            pi = pvar[ii];
            pia = pvara[ii];

            diipbx = diipb[face_id][0];
            diipby = diipb[face_id][1];
            diipbz = diipb[face_id][2];

            /* Remove decentering for coupled faces */
            if (ifaccp==1 && bc_type[face_id]==CS_COUPLED) {
              flui = 0.0;
              fluj = b_massflux[face_id];
            } else {
              flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
              fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
            }

            pir  = pi/relaxp - (1.-relaxp)/relaxp*pia;
            pipr =   pir
                   + ircflp*(  grad[ii][0]*diipbx
                             + grad[ii][1]*diipby
                             + grad[ii][2]*diipbz);

            pfac  = inc*coefap[face_id] +coefbp[face_id]*pipr;
            pfacd = inc*cofafp[face_id] +cofbfp[face_id]*pipr;

            flux =   iconvp*(flui*pir + fluj*pfac - b_massflux[face_id]*pi )
                   + idiffp*b_visc[face_id]*pfacd;
            rhs[ii] = rhs[ii] - flux;

          }
        }
      }

    /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz,    \
                                        flui, fluj, pip, pfac, pfacd, flux, pi) \
                            if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            pi = pvar[ii];

            diipbx = diipb[face_id][0];
            diipby = diipb[face_id][1];
            diipbz = diipb[face_id][2];

            /* Remove decentering for coupled faces */
            if (ifaccp==1 && bc_type[face_id]==CS_COUPLED) {
              flui = 0.0;
              fluj = b_massflux[face_id];
            } else {
              flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
              fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
            }

            pip = pi
                + ircflp*(  grad[ii][0]*diipbx
                          + grad[ii][1]*diipby
                          + grad[ii][2]*diipbz);

            pfac  = inc*coefap[face_id] + coefbp[face_id]*pip;
            pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pip;

            flux =   iconvp*((flui - b_massflux[face_id])*pi + fluj*pfac)
                   + idiffp*b_visc[face_id]*pfacd;
            rhs[ii] = rhs[ii] - thetap * flux;

          }
        }
      }

    }

  /* Boundary convective flux is imposed at some faces
     (tagged in icvfli array) */
  } else if (icvflb==1) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    } else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    /* Steady */
    if (idtvar<0) {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz, \
                                        flui, fluj, pir, pipr, pfac, pfacd,  \
                                        flux, pi, pia)                       \
                            if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            pi = pvar[ii];
            pia = pvara[ii];

            diipbx = diipb[face_id][0];
            diipby = diipb[face_id][1];
            diipbz = diipb[face_id][2];

            pir  = pi/relaxp - (1.-relaxp)/relaxp*pia;
            pipr =   pir
                   + ircflp*(  grad[ii][0]*diipbx
                             + grad[ii][1]*diipby
                             + grad[ii][2]*diipbz);

            /* Computed convective flux */
            if (icvfli[face_id] == 0) {
              /* Remove decentering for coupled faces */
              if (ifaccp==1 && bc_type[face_id]==CS_COUPLED) {
                flui = 0.0;
                fluj = b_massflux[face_id];
              } else {
                flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
                fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
              }
              pfac = inc*coefap[face_id] +coefbp[face_id]*pipr;
              flux = iconvp*(flui*pir + fluj*pfac - b_massflux[face_id]*pi);

            /* Imposed convective flux */
            } else {
              pfac = inc*coface[face_id] + cofbce[face_id]*pipr;
              flux = iconvp*(- b_massflux[face_id]*pi + pfac);
            }

            pfacd = inc*cofafp[face_id] +cofbfp[face_id]*pipr;

            flux += idiffp*b_visc[face_id]*pfacd;
            rhs[ii] = rhs[ii] - flux;

          }
        }
      }

    /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz,   \
                                        flui, fluj, pip, pfac, pfacd, flux, pi)\
                            if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            pi = pvar[ii];

            diipbx = diipb[face_id][0];
            diipby = diipb[face_id][1];
            diipbz = diipb[face_id][2];

            pip = pi
                + ircflp*(  grad[ii][0]*diipbx
                          + grad[ii][1]*diipby
                          + grad[ii][2]*diipbz);

            /* Computed convective flux */
            if (icvfli[face_id] == 0) {
              /* Remove decentering for coupled faces */
              if (ifaccp==1 && bc_type[face_id]==CS_COUPLED) {
                flui = 0.0;
                fluj = b_massflux[face_id];
              } else {
                flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
                fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
              }

              pfac  = inc*coefap[face_id] + coefbp[face_id]*pip;
              flux = iconvp*((flui- b_massflux[face_id])*pi + fluj*pfac);

            /* Imposed convective flux */
            } else {
              pfac = inc*coface[face_id] + cofbce[face_id]*pip;
              flux = iconvp*(- b_massflux[face_id]*pi + pfac);
            }

            pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pip;

            flux += idiffp*b_visc[face_id]*pfacd;
            rhs[ii] = rhs[ii] - thetap * flux;

          }
        }
      }

    }
  }

  /* Free memory */
  BFT_FREE(grad);
  BFT_FREE(grdpa);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     coefav        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_vector(int                         idtvar,
                               int                         f_id,
                               const cs_var_cal_opt_t      var_cal_opt,
                               int                         icvflb,
                               int                         inc,
                               int                         ifaccp,
                               int                         ivisep,
                               cs_real_3_t       *restrict pvar,
                               const cs_real_3_t *restrict pvara,
                               const cs_int_t              bc_type[],
                               const cs_int_t              icvfli[],
                               const cs_real_3_t           coefav[],
                               const cs_real_33_t          coefbv[],
                               const cs_real_3_t           cofafv[],
                               const cs_real_33_t          cofbfv[],
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               const cs_real_t             secvif[],
                               cs_real_3_t       *restrict rhs)
{
  const int iconvp = var_cal_opt.iconv;
  const int idiffp = var_cal_opt.idiff;
  const int nswrgp = var_cal_opt.nswrgr;
  const int imrgra = var_cal_opt.imrgra;
  const int imligp = var_cal_opt.imligr;
  const int ircflp = var_cal_opt.ircflu;
  const int ischcp = var_cal_opt.ischcv;
  const int isstpp = var_cal_opt.isstpc;
  const int iwarnp = var_cal_opt.iwarni;
  const double blencp = var_cal_opt.blencv;
  const double epsrgp = var_cal_opt.epsrgr;
  const double climgp = var_cal_opt.climgr;
  const double relaxp = var_cal_opt.relaxv;
  const double thetap = var_cal_opt.thetav;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* Local variables */

  char var_name[32];

  cs_gnum_t n_upwind;
  cs_lnum_t face_id, ii, jj, cell_id;
  int iupwin, g_id, t_id;
  int isou, jsou, ityp;
  double pfac,pfacd,flui,fluj,flux,fluxi,fluxj;
  double vfac[3];
  double difv[3], djfv[3];
  double pi , pj, pia, pja;
  double pif,pjf,pip,pjp,pir,pjr,pipr,pjpr;
  double pifr,pjfr,pifri,pifrj,pjfri,pjfrj;
  double testi,testj,testij;
  double dpvf[3];
  double dcc, ddi, ddj, tesqck;
  double dijpfv[3];
  double diipfv[3];
  double djjpfv[3];
  double diipbv[3];
  double pnd, distf, srfan;
  double unsvol, visco, grdtrv, tgrdfl, secvis;

  cs_real_33_t *gradv, *gradva;
  cs_real_t *bndcel;

  cs_real_3_t *cofacv;
  cs_real_33_t *cofbcv;

  cs_field_t *f;

  cs_real_t  *v_slope_test = _get_v_slope_test(f_id,  var_cal_opt);

  /*==========================================================================*/

  /* 1. Initialization */

  /* Allocate work arrays */

  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);
  BFT_MALLOC(gradva, n_cells_ext, cs_real_33_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf
        (
         _(" %s: Convection in centered blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    } else {
      bft_printf
        (
         _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
  }

  iupwin = (blencp > 0.) ? 0 : 1;

  /* 2. Compute the balance with reconstruction */

  /* Compute the gradient of the velocity

     The gradient (gradv) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test. */

  if (  (idiffp != 0 && ircflp == 1) || ivisep == 1
     || (   iconvp != 0 && iupwin == 0
         && (ischcp == 0 || ircflp == 1 || isstpp == 0))) {

    cs_gradient_vector(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       nswrgp,
                       iwarnp,
                       imligp,
                       epsrgp,
                       climgp,
                       coefav,
                       coefbv,
                       pvar,
                       gradv);

  } else {
#   pragma omp parallel for private(isou, jsou)
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++)
          gradv[cell_id][isou][jsou] = 0.;
      }
    }
  }

  /* ======================================================================
     ---> Compute uncentred gradient gradva for the slope test
     ======================================================================*/

# pragma omp parallel for private(isou, jsou)
  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    for (jsou = 0; jsou < 3; jsou++) {
      for (isou = 0; isou < 3; isou++)
        gradva[cell_id][isou][jsou] = 0.;
    }
  }

  if (iconvp > 0 && iupwin == 0 && isstpp == 0) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, isou, jsou,      \
                                      difv, djfv, pif, pjf, pfac, vfac)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          for (jsou = 0; jsou < 3; jsou++) {
            difv[jsou] = i_face_cog[face_id][jsou] - cell_cen[ii][jsou];
            djfv[jsou] = i_face_cog[face_id][jsou] - cell_cen[jj][jsou];
          }

          /*-----------------
            X-Y-Z component, p=u, v, w */

          for (isou = 0; isou < 3; isou++) {
            pif = pvar[ii][isou];
            pjf = pvar[jj][isou];
            for (jsou = 0; jsou < 3; jsou++) {
              pif = pif + gradv[ii][isou][jsou]*difv[jsou];
              pjf = pjf + gradv[jj][isou][jsou]*djfv[jsou];
            }

            pfac = pjf;
            if (i_massflux[face_id] > 0.) pfac = pif;

            /* U gradient */

            for (jsou = 0; jsou < 3; jsou++) {
              vfac[jsou] = pfac*i_face_normal[face_id][jsou];

              gradva[ii][isou][jsou] = gradva[ii][isou][jsou] + vfac[jsou];
              gradva[jj][isou][jsou] = gradva[jj][isou][jsou] - vfac[jsou];
            }
          }

        }
      }
    }

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, isou, jsou, diipbv, pfac) \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          for (jsou = 0; jsou < 3; jsou++)
            diipbv[jsou] = diipb[face_id][jsou];

          /*-----------------
            X-Y-Z components, p=u, v, w */
          for (isou = 0; isou <3; isou++) {
            pfac = inc*coefav[face_id][isou];
            /*coefu is a matrix */
            for (jsou =  0; jsou < 3; jsou++)
              pfac += coefbv[face_id][jsou][isou]*(  pvar[ii][jsou]
                                                   + gradv[ii][jsou][0]*diipbv[0]
                                                   + gradv[ii][jsou][1]*diipbv[1]
                                                   + gradv[ii][jsou][2]*diipbv[2]);

            for (jsou = 0; jsou < 3; jsou++)
              gradva[ii][isou][jsou] += pfac*b_face_normal[face_id][jsou];
          }

        }
      }
    }

#   pragma omp parallel for private(isou, jsou, unsvol)
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      unsvol = 1./cell_vol[cell_id];
      for (isou = 0; isou < 3; isou++) {
        for (jsou = 0; jsou < 3; jsou++)
          gradva[cell_id][isou][jsou] = gradva[cell_id][isou][jsou]*unsvol;
      }
    }

    /* Handle parallelism and periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)gradva, 9);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_tens(halo, halo_type, (cs_real_t *)gradva);
    }

  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  n_upwind = 0;

  if (n_cells_ext > n_cells) {
#   pragma omp parallel for private(isou) if(n_cells_ext -n_cells > CS_THR_MIN)
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      for (isou = 0; isou < 3; isou++)
        rhs[cell_id][isou] = 0.;
    }
  }

  /* --> Pure upwind flux
     =====================*/

  if (iupwin == 1) {

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,   \
                                        dijpfv, diipfv, djjpfv,             \
                                        flui, fluj, dpvf, pi, pj, pia, pja, \
                                        pip, pjp, pipr, pjpr, pifr, pjfr,   \
                                        fluxi, fluxj)                       \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells) {
              n_upwind++;
            }

            for (jsou = 0; jsou < 3; jsou++) {
              dijpfv[jsou] = dijpf[face_id][jsou];
            }

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */
            for (jsou = 0; jsou < 3; jsou++) {
              diipfv[jsou] =   i_face_cog[face_id][jsou]
                             - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
              djjpfv[jsou] =   i_face_cog[face_id][jsou]
                             - cell_cen[jj][jsou] + pnd  * dijpfv[jsou];
            }

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              for (jsou = 0; jsou < 3; jsou++)
                dpvf[jsou] = 0.5*(  gradv[ii][isou][jsou]
                                  + gradv[jj][isou][jsou]);

              /* reconstruction only if IRCFLP = 1 */
              pi  = pvar[ii][isou];
              pj  = pvar[jj][isou];

              pia = pvara[ii][isou];
              pja = pvara[jj][isou];

              pip = pi + ircflp*(  dpvf[0]*diipfv[0]
                                 + dpvf[1]*diipfv[1]
                                 + dpvf[2]*diipfv[2]);
              pjp = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                 + dpvf[1]*djjpfv[1]
                                 + dpvf[2]*djjpfv[2]);

              pipr = pi /relaxp - (1.-relaxp)/relaxp * pia
                     + ircflp*(  dpvf[0]*diipfv[0]
                               + dpvf[1]*diipfv[1]
                               + dpvf[2]*diipfv[2]);
              pjpr = pj /relaxp - (1.-relaxp)/relaxp * pja
                     + ircflp*(  dpvf[0]*djjpfv[0]
                               + dpvf[1]*djjpfv[1]
                               + dpvf[2]*djjpfv[2]);

              pifr = pi /relaxp - (1.-relaxp)/relaxp * pia;
              pjfr = pj /relaxp - (1.-relaxp)/relaxp * pja;

              fluxi =   iconvp*(flui*pifr + fluj*pj - i_massflux[face_id]*pi)
                      + idiffp*i_visc[face_id]*(pipr -pjp);
              fluxj =   iconvp*(flui*pi + fluj*pjfr - i_massflux[face_id]*pj)
                      + idiffp*i_visc[face_id]*(pip -pjpr);

              rhs[ii][isou] = rhs[ii][isou] - fluxi;
              rhs[jj][isou] = rhs[jj][isou] + fluxj;

            }

          }
        }
      }

    /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,   \
                                        dijpfv, diipfv, djjpfv, flui, fluj, \
                                        dpvf, pi, pj, pip, pjp, flux)       \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells) {
              n_upwind++;
            }

            for (jsou = 0; jsou < 3; jsou++)
              dijpfv[jsou] = dijpf[face_id][jsou];

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */
            for (jsou = 0; jsou < 3; jsou++) {
              diipfv[jsou] =   i_face_cog[face_id][jsou]
                             - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
              djjpfv[jsou] =   i_face_cog[face_id][jsou]
                             - cell_cen[jj][jsou] +  pnd * dijpfv[jsou];
            }

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              for (jsou = 0; jsou < 3; jsou++)
                dpvf[jsou] = 0.5*(  gradv[ii][isou][jsou]
                                  + gradv[jj][isou][jsou]);

              pi = pvar[ii][isou];
              pj = pvar[jj][isou];

              pip = pi + ircflp*(  dpvf[0]*diipfv[0]
                                 + dpvf[1]*diipfv[1]
                                 + dpvf[2]*diipfv[2]);
              pjp = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                 + dpvf[1]*djjpfv[1]
                                 + dpvf[2]*djjpfv[2]);

              flux =   iconvp*(flui*pi +fluj*pj)
                + idiffp*i_visc[face_id]*(pip -pjp);

              rhs[ii][isou] -= thetap*(flux - iconvp*i_massflux[face_id]*pi);
              rhs[jj][isou] += thetap*(flux - iconvp*i_massflux[face_id]*pj);
            }

          }
        }
      }

    }

    /* --> Flux with no slope test
       ============================*/

  } else if (isstpp == 1) {

    if (ischcp<0 || ischcp>1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,    \
                                        dijpfv, diipfv, djjpfv, flui, fluj,  \
                                        difv, djfv, dpvf, pi, pj, pia, pja,  \
                                        pip, pjp, pipr, pjpr, pir, pjr,      \
                                        pifri, pjfri, pifrj, pjfrj, fluxi, fluxj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            for (jsou = 0; jsou < 3; jsou++)
              dijpfv[jsou] = dijpf[face_id][jsou];

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */
            for (jsou = 0; jsou < 3; jsou++) {
              diipfv[jsou] =   i_face_cog[face_id][jsou]
                             - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
              djjpfv[jsou] =   i_face_cog[face_id][jsou]
                             - cell_cen[jj][jsou] + pnd  * dijpfv[jsou];
            }

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            /* For second order, define IF */
            if (ischcp == 0) {
              for (jsou = 0; jsou < 3; jsou++) {
                difv[jsou] = i_face_cog[face_id][jsou] - cell_cen[ii][jsou];
                djfv[jsou] = i_face_cog[face_id][jsou] - cell_cen[jj][jsou];
              }
            }

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              for (jsou = 0; jsou < 3; jsou++)
                dpvf[jsou] = 0.5*(  gradv[ii][isou][jsou]
                                  + gradv[jj][isou][jsou]);

              pi = pvar[ii][isou];
              pj = pvar[jj][isou];

              pia = pvara[ii][isou];
              pja = pvara[jj][isou];

              pip = pi + ircflp*(  dpvf[0]*diipfv[0]
                                 + dpvf[1]*diipfv[1]
                                 + dpvf[2]*diipfv[2]);
              pjp = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                 + dpvf[1]*djjpfv[1]
                                 + dpvf[2]*djjpfv[2]);

              pipr = pi /relaxp - (1.-relaxp)/relaxp * pia
                     + ircflp*(  dpvf[0]*diipfv[0]
                               + dpvf[1]*diipfv[1]
                               + dpvf[2]*diipfv[2]);
              pjpr = pj /relaxp - (1.-relaxp)/relaxp * pja
                     + ircflp*(  dpvf[0]*djjpfv[0]
                               + dpvf[1]*djjpfv[1]
                               + dpvf[2]*djjpfv[2]);

              pir = pi /relaxp - (1. - relaxp)/relaxp* pia;
              pjr = pj /relaxp - (1. - relaxp)/relaxp* pja;

              /* Centered
                 --------*/

              if (ischcp == 1) {

                pifri = pnd*pipr + (1.-pnd)*pjp;
                pjfri = pifri;
                pifrj = pnd*pip  + (1.-pnd)*pjpr;
                pjfrj = pifrj;

              /* Second order
                 ------------*/

              } else {/* if (ischcp.eq.0) then
                       dif* is already defined */

                /* leave reconstruction of PIF and PJF even if IRCFLP=0
                   otherwise, it is the same as using upwind */
                pifri = pir + difv[0]*gradv[ii][isou][0]
                            + difv[1]*gradv[ii][isou][1]
                            + difv[2]*gradv[ii][isou][2];
                pifrj = pi  + difv[0]*gradv[ii][isou][0]
                            + difv[1]*gradv[ii][isou][1]
                            + difv[2]*gradv[ii][isou][2];

                pjfrj = pjr + djfv[0]*gradv[jj][isou][0]
                            + djfv[1]*gradv[jj][isou][1]
                            + djfv[2]*gradv[jj][isou][2];
                pjfri = pj  + djfv[0]*gradv[jj][isou][0]
                            + djfv[1]*gradv[jj][isou][1]
                            + djfv[2]*gradv[jj][isou][2];

              }

              /* Blending
                 --------*/

              pifri = blencp*pifri+(1.-blencp)*pir;
              pifrj = blencp*pifrj+(1.-blencp)*pi;
              pjfri = blencp*pjfri+(1.-blencp)*pj;
              pjfrj = blencp*pjfrj+(1.-blencp)*pjr;

              /* Flux
                 ----*/

              fluxi =   iconvp*(  flui*pifri + fluj*pjfri
                                - i_massflux[face_id]*pi)
                      + idiffp*i_visc[face_id]*(pipr -pjp);
              fluxj =   iconvp*(  flui*pifrj + fluj*pjfrj
                                - i_massflux[face_id]*pj)
                      + idiffp*i_visc[face_id]*(pip -pjpr);

              /* Assembly
                 --------*/

              rhs[ii][isou] = rhs[ii][isou] - fluxi;
              rhs[jj][isou] = rhs[jj][isou] + fluxj;

            } /* isou */

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,     \
                                        dijpfv, diipfv, djjpfv, flui, fluj,   \
                                        difv, djfv, dpvf,                     \
                                        pi, pj, pip, pjp, pif, pjf, flux)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            for (jsou = 0; jsou < 3; jsou++) {
              dijpfv[jsou] = dijpf[face_id][jsou];
            }

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */
            for (jsou = 0; jsou < 3; jsou++) {
              diipfv[jsou] =   i_face_cog[face_id][jsou]
                             - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
              djjpfv[jsou] =   i_face_cog[face_id][jsou]
                             - cell_cen[jj][jsou] + pnd  * dijpfv[jsou];
            }

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            /* For second order, define IF */
            if (ischcp == 0) {
              for (jsou = 0; jsou < 3; jsou++) {
                difv[jsou] = i_face_cog[face_id][jsou] - cell_cen[ii][jsou];
                djfv[jsou] = i_face_cog[face_id][jsou] - cell_cen[jj][jsou];
              }
            }

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              for (jsou = 0; jsou < 3; jsou++) {
                dpvf[jsou] = 0.5*(  gradv[ii][isou][jsou]
                                  + gradv[jj][isou][jsou]);
              }

              pi = pvar[ii][isou];
              pj = pvar[jj][isou];

              pip = pi + ircflp*(  dpvf[0]*diipfv[0]
                                 + dpvf[1]*diipfv[1]
                                 + dpvf[2]*diipfv[2]);
              pjp = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                 + dpvf[1]*djjpfv[1]
                                 + dpvf[2]*djjpfv[2]);

              /* Centered
                 --------*/

              if (ischcp == 1) {

                pif = pnd*pip +(1.-pnd)*pjp;
                pjf = pif;

              /* Second order
                 ------------*/

              } else {/* if (ischcp.eq.0) then
                       dif* is already defined */

                /* leave reconstruction of PIF and PJF even if IRCFLP=0
                   otherwise, it is the same as using upwind */
                pif = pi;
                pjf = pj;
                for (jsou = 0; jsou < 3; jsou++) {
                  pif = pif + gradv[ii][isou][jsou]*difv[jsou];
                  pjf = pjf + gradv[jj][isou][jsou]*djfv[jsou];
                }

              }

              /* Blending
                 --------*/

              pif = blencp*pif+(1.-blencp)*pi;
              pjf = blencp*pjf+(1.-blencp)*pj;

              /* Flux
                 ----*/

              flux =   iconvp*(flui*pif +fluj*pjf)
                     + idiffp*i_visc[face_id]*(pip -pjp);

              /* Assembly
                 --------*/

              rhs[ii][isou] -= thetap*(flux - iconvp*i_massflux[face_id]*pi);
              rhs[jj][isou] += thetap*(flux - iconvp*i_massflux[face_id]*pj);

            } /* isou */

          }
        }
      }

    }

    /* --> Flux with slope test
       =========================*/

  } else {

    if (ischcp<0 || ischcp>1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,      \
                                        distf, srfan, dijpfv, diipfv, djjpfv,  \
                                        flui, fluj, difv, djfv, dpvf, pi, pj,  \
                                        pia, pja, pip, pjp, pipr, pjpr,        \
                                        pir, pjr, testij, testi, testj, dcc,   \
                                        ddi, ddj, tesqck, pifri, pifrj, pjfri, \
                                        pjfrj, fluxi, fluxj)                   \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            for (jsou = 0; jsou < 3; jsou++) {
              dijpfv[jsou] = dijpf[face_id][jsou];
            }

            pnd   = weight[face_id];
            distf = i_dist[face_id];
            srfan = i_face_surf[face_id];

            /* Recompute II' and JJ' at this level */
            for (jsou = 0; jsou < 3; jsou++) {
              diipfv[jsou] =   i_face_cog[face_id][jsou]
                             - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
              djjpfv[jsou] =    i_face_cog[face_id][jsou]
                             - cell_cen[jj][jsou] + pnd * dijpfv[jsou];
            }

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            /* For second order, define IF */
            if (ischcp == 0) {
              for (jsou = 0; jsou < 3; jsou++) {
                difv[jsou] = i_face_cog[face_id][jsou] - cell_cen[ii][jsou];
                djfv[jsou] = i_face_cog[face_id][jsou] - cell_cen[jj][jsou];
              }
            }

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              for (jsou = 0; jsou < 3; jsou++)
                dpvf[jsou] = 0.5*(  gradv[ii][isou][jsou]
                                  + gradv[jj][isou][jsou]);

              pi  = pvar[ii][isou];
              pj  = pvar[jj][isou];

              pia = pvara[ii][isou];
              pja = pvara[jj][isou];

              pip = pi + ircflp*(  dpvf[0]*diipfv[0]
                                 + dpvf[1]*diipfv[1]
                                 + dpvf[2]*diipfv[2]);
              pjp = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                 + dpvf[1]*djjpfv[1]
                                 + dpvf[2]*djjpfv[2]);

              pipr = pi /relaxp - (1.-relaxp)/relaxp * pia
                     + ircflp*(  dpvf[0]*diipfv[0]
                               + dpvf[1]*diipfv[1]
                               + dpvf[2]*diipfv[2]);
              pjpr = pj /relaxp - (1.-relaxp)/relaxp * pja
                     + ircflp*(  dpvf[0]*djjpfv[0]
                               + dpvf[1]*djjpfv[1]
                               + dpvf[2]*djjpfv[2]);

              pir = pi /relaxp - (1. - relaxp)/relaxp*pia;
              pjr = pj /relaxp - (1. - relaxp)/relaxp*pja;

              /* Slope test
                 ----------*/

              testi =   gradva[ii][isou][0]*i_face_normal[face_id][0]
                      + gradva[ii][isou][1]*i_face_normal[face_id][1]
                      + gradva[ii][isou][2]*i_face_normal[face_id][2];
              testj =   gradva[jj][isou][0]*i_face_normal[face_id][0]
                      + gradva[jj][isou][1]*i_face_normal[face_id][1]
                      + gradva[jj][isou][2]*i_face_normal[face_id][2];
              testij =   gradva[ii][isou][0]*gradva[jj][isou][0]
                       + gradva[ii][isou][1]*gradva[jj][isou][1]
                       + gradva[ii][isou][2]*gradva[jj][isou][2];

              if (i_massflux[face_id] > 0.) {
                dcc =   gradv[ii][isou][0]*i_face_normal[face_id][0]
                      + gradv[ii][isou][1]*i_face_normal[face_id][1]
                      + gradv[ii][isou][2]*i_face_normal[face_id][2];
                ddi = testi;
                ddj = (pj - pi)/distf *srfan;
              }
              else {
                dcc =   gradv[jj][isou][0]*i_face_normal[face_id][0]
                      + gradv[jj][isou][1]*i_face_normal[face_id][1]
                      + gradv[jj][isou][2]*i_face_normal[face_id][2];
                ddi = (pj - pi)/distf *srfan;
                ddj = testj;
              }
              tesqck = pow(dcc,2) -pow((ddi-ddj),2);

              /* Upwind
                 ------*/

              if (tesqck <= 0. || testij <= 0.) {

                pifri = pir;
                pifrj = pi;
                pjfri = pj;
                pjfrj = pjr;
                /* in parallel, face will be counted by one and only one rank */
                if (ii < n_cells)
                  n_upwind++;
                if (v_slope_test != NULL) {
                  v_slope_test[ii] += fabs(i_massflux[face_id]) / cell_vol[ii];
                  v_slope_test[jj] += fabs(i_massflux[face_id]) / cell_vol[jj];
                }
              }
              else {

                /* Centered
                   --------*/

                if (ischcp == 1) {

                  pifri = pnd*pipr +(1.-pnd)*pjp;
                  pjfri = pifri;
                  pifrj = pnd*pip  +(1.-pnd)*pjpr;
                  pjfrj = pifrj;

                /* Second order
                   ------------*/

                }
                else { /* if (ischcp.eq.0) then difv already defined */

                  /* leave reconstruction of PIF and PJF even if IRCFLP=0
                     otherwise, it is the same as using upwind */
                  pifri = pir + difv[0]*gradv[ii][isou][0]
                              + difv[1]*gradv[ii][isou][1]
                              + difv[2]*gradv[ii][isou][2];
                  pifrj = pi  + difv[0]*gradv[ii][isou][0]
                              + difv[1]*gradv[ii][isou][1]
                              + difv[2]*gradv[ii][isou][2];

                  pjfrj = pjr + djfv[0]*gradv[jj][isou][0]
                              + djfv[1]*gradv[jj][isou][1]
                              + djfv[2]*gradv[jj][isou][2];
                  pjfri = pj  + djfv[0]*gradv[jj][isou][0]
                              + djfv[1]*gradv[jj][isou][1]
                              + djfv[2]*gradv[jj][isou][2];

                }

              }

              /* Blending
                 --------*/

              pifri = blencp*pifri+(1.-blencp)*pir;
              pifrj = blencp*pifrj+(1.-blencp)*pi;
              pjfri = blencp*pjfri+(1.-blencp)*pj;
              pjfrj = blencp*pjfrj+(1.-blencp)*pjr;

              /* Flux
                 ----*/

              fluxi =    iconvp*(  flui*pifri + fluj*pjfri
                                 - i_massflux[face_id]*pi)
                + idiffp*i_visc[face_id]*(pipr -pjp);
              fluxj =    iconvp*(  flui*pifrj + fluj*pjfrj
                                 - i_massflux[face_id]*pj)
                + idiffp*i_visc[face_id]*(pip -pjpr);

              /* Assembly
                 --------*/

              rhs[ii][isou] = rhs[ii][isou] - fluxi;
              rhs[jj][isou] = rhs[jj][isou] + fluxj;

            } /* isou */

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,     \
                                        distf, srfan, dijpfv, diipfv, djjpfv, \
                                        flui, fluj, difv, djfv, dpvf, pi, pj, \
                                        pip, pjp, testi, testj, testij, dcc,  \
                                        ddi, ddj, tesqck, pif, pjf, flux)     \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            for (jsou = 0; jsou < 3; jsou++) {
              dijpfv[jsou] = dijpf[face_id][jsou];
            }

            pnd   = weight[face_id];
            distf = i_dist[face_id];
            srfan = i_face_surf[face_id];

            /* Recompute II' and JJ' at this level */
            for (jsou = 0; jsou < 3; jsou++) {
              diipfv[jsou] =   i_face_cog[face_id][jsou]
                             - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
              djjpfv[jsou] =   i_face_cog[face_id][jsou]
                             - cell_cen[jj][jsou] + pnd * dijpfv[jsou];
            }

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            /* For second order, define IF */
            if (ischcp == 0) {
              for (jsou = 0; jsou < 3; jsou++) {
                difv[jsou] = i_face_cog[face_id][jsou] - cell_cen[ii][jsou];
                djfv[jsou] = i_face_cog[face_id][jsou] - cell_cen[jj][jsou];
              }
            }

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              for (jsou = 0; jsou < 3; jsou++) {
                dpvf[jsou] = 0.5*(  gradv[ii][isou][jsou]
                                  + gradv[jj][isou][jsou]);
              }

              pi = pvar[ii][isou];
              pj = pvar[jj][isou];

              pip = pi + ircflp*(dpvf[0]*diipfv[0]
                                 +dpvf[1]*diipfv[1]
                                 +dpvf[2]*diipfv[2]);
              pjp = pj + ircflp*(dpvf[0]*djjpfv[0]
                                 +dpvf[1]*djjpfv[1]
                                 +dpvf[2]*djjpfv[2]);

              /* Slope test
                 ----------*/

              testi = gradva[ii][isou][0]*i_face_normal[face_id][0]
                    + gradva[ii][isou][1]*i_face_normal[face_id][1]
                    + gradva[ii][isou][2]*i_face_normal[face_id][2];
              testj = gradva[jj][isou][0]*i_face_normal[face_id][0]
                    + gradva[jj][isou][1]*i_face_normal[face_id][1]
                    + gradva[jj][isou][2]*i_face_normal[face_id][2];
              testij=   gradva[ii][isou][0]*gradva[jj][isou][0]
                      + gradva[ii][isou][1]*gradva[jj][isou][1]
                      + gradva[ii][isou][2]*gradva[jj][isou][2];

              if (i_massflux[face_id] > 0.) {
                dcc =   gradv[ii][isou][0]*i_face_normal[face_id][0]
                      + gradv[ii][isou][1]*i_face_normal[face_id][1]
                      + gradv[ii][isou][2]*i_face_normal[face_id][2];
                ddi = testi;
                ddj = (pj - pi)/distf *srfan;
              } else {
                dcc =   gradv[jj][isou][0]*i_face_normal[face_id][0]
                      + gradv[jj][isou][1]*i_face_normal[face_id][1]
                      + gradv[jj][isou][2]*i_face_normal[face_id][2];
                ddi = (pj - pi)/distf *srfan;
                ddj = testj;
              }
              tesqck = pow(dcc,2) -pow((ddi-ddj),2);

              /* Upwind
                 ------*/

              if (tesqck <= 0. || testij <= 0.) {

                pif = pi;
                pjf = pj;
                /* in parallel, face will be counted by one and only one rank */
                if (ii < n_cells)
                  n_upwind++;
                if (v_slope_test != NULL) {
                  v_slope_test[ii] += fabs(i_massflux[face_id]) / cell_vol[ii];
                  v_slope_test[jj] += fabs(i_massflux[face_id]) / cell_vol[jj];
                }


              } else {

                /* Centered
                   --------*/

                if (ischcp == 1) {

                  pif = pnd*pip +(1.-pnd)*pjp;
                  pjf = pif;

                  /* Second order
                     ------------*/

                } else { /* if (ischcp.eq.0) then
                            dif* is already defined */

                  pif = pi;
                  pjf = pj;
                  for (jsou = 0; jsou < 3; jsou++) {
                    pif = pif + gradv[ii][isou][jsou]*difv[jsou];
                    pjf = pjf + gradv[jj][isou][jsou]*djfv[jsou];
                  }

                  /* leave reconstruction of PIF and PJF even if IRCFLP=0
                     otherwise, it is the same as using upwind */

                }

              }

              /* Blending
                 --------*/

              pif = blencp*pif+(1.-blencp)*pi;
              pjf = blencp*pjf+(1.-blencp)*pj;

              /* Flux
                 ----*/

              flux =   iconvp*(flui*pif +fluj*pjf)
                     + idiffp*i_visc[face_id]*(pip -pjp);

              /* Assembly
                 --------*/

              rhs[ii][isou] -= thetap*(flux - iconvp*i_massflux[face_id]*pi);
              rhs[jj][isou] += thetap*(flux - iconvp*i_massflux[face_id]*pj);

            } /* isou */

          }
        }
      }

    } /* idtvar */

  } /* iupwin */

  if (iwarnp >= 2) {

    /* Sum number of clippings */
    cs_parall_counter(&n_upwind, 1);

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces \n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_faces);
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, isou, jsou, diipbv, \
                                        flui, fluj, pfac, pfacd, pir,    \
                                        pipr, flux, pi, pia)            \
                   if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            for (jsou = 0; jsou < 3; jsou++)
              diipbv[jsou] = diipb[face_id][jsou];

            /* Remove decentering for coupled faces */
            if (ifaccp == 1 && bc_type[face_id] == CS_COUPLED) {
              flui = 0.0;
              fluj = b_massflux[face_id];
            } else {
              flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
              fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
            }

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              pfac  = inc*coefav[face_id][isou];
              pfacd = inc*cofafv[face_id][isou];

              /*coefu and cofuf are matrices */
              for (jsou = 0; jsou < 3; jsou++) {
                pir  =   pvar[ii][jsou]/relaxp
                       - (1.-relaxp)/relaxp*pvara[ii][jsou];

                pipr = pir +ircflp*(   gradv[ii][jsou][0]*diipbv[0]
                                     + gradv[ii][jsou][1]*diipbv[1]
                                     + gradv[ii][jsou][2]*diipbv[2]);
                pfac  = pfac  + coefbv[face_id][jsou][isou]*pipr;
                pfacd = pfacd + cofbfv[face_id][jsou][isou]*pipr;
              }

              pi  = pvar[ii][isou];
              pia = pvara[ii][isou];

              pir  = pi/relaxp - (1.-relaxp)/relaxp*pia;

              flux = iconvp*(flui*pir + fluj*pfac - b_massflux[face_id]*pi)
                + idiffp*b_visc[face_id]*pfacd;
              rhs[ii][isou] = rhs[ii][isou] - flux;

            } /* isou */

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, isou, jsou, diipbv,         \
                                        flui, fluj, pfac, pfacd, pip, flux, pi)  \
                   if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            for (jsou = 0; jsou < 3; jsou++) {
              diipbv[jsou] = diipb[face_id][jsou];
            }

            /* Remove decentering for coupled faces */
            if (ifaccp == 1 && bc_type[face_id] == CS_COUPLED) {
              flui = 0.0;
              fluj = b_massflux[face_id];
            } else {
              flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
              fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
            }

            /*-----------------
              X-Y-Z components, p=u, v, w */
            for (isou = 0; isou < 3; isou++) {

              pfac  = inc*coefav[face_id][isou];
              pfacd = inc*cofafv[face_id][isou];

              /*coefu and cofuf are matrices */
              for (jsou = 0; jsou < 3; jsou++) {
                pip = pvar[ii][jsou] + ircflp*(   gradv[ii][jsou][0]*diipbv[0]
                                                + gradv[ii][jsou][1]*diipbv[1]
                                                + gradv[ii][jsou][2]*diipbv[2]);
                pfac  = pfac  + coefbv[face_id][jsou][isou]*pip;
                pfacd = pfacd + cofbfv[face_id][jsou][isou]*pip;
              }

              pi = pvar[ii][isou];

              flux = iconvp*((flui-b_massflux[face_id])*pi + fluj*pfac)
                + idiffp*b_visc[face_id]*pfacd;
              rhs[ii][isou] = rhs[ii][isou] - thetap * flux;

            } /* isou */

          }
        }
      }

    } /* idtvar */

    /* Boundary convective flux imposed at some faces (tags in icvfli array) */
  } else if (icvflb==1) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      cofacv = (cs_real_3_t *)(f->bc_coeffs->ac);
      cofbcv = (cs_real_33_t *)(f->bc_coeffs->bc);
    } else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, isou, jsou, diipbv,    \
                                        flui, fluj, pfac, pfacd, pir, pipr, \
                                        flux, pi, pia)                      \
                   if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            for (jsou = 0; jsou < 3; jsou++) {
              diipbv[jsou] = diipb[face_id][jsou];
            }

            /* Computed convective flux */
            if (icvfli[face_id] == 0) {
              /* Remove decentering for coupled faces */
              if (ifaccp == 1 && bc_type[face_id] == CS_COUPLED) {
                flui = 0.0;
                fluj = b_massflux[face_id];
              } else {
                flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
                fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
              }

              /*-----------------
                X-Y-Z components, p=u, v, w */
              for (isou = 0; isou < 3; isou++) {

                pfac  = inc*coefav[face_id][isou];
                pfacd = inc*cofafv[face_id][isou];

                /*coefu and cofuf are matrices */
                for (jsou = 0; jsou < 3; jsou++) {
                  pir  =   pvar[ii][jsou]/relaxp
                         - (1.-relaxp)/relaxp*pvara[ii][jsou];

                  pipr = pir +ircflp*( gradv[ii][jsou][0]*diipbv[0]
                                       + gradv[ii][jsou][1]*diipbv[1]
                                       + gradv[ii][jsou][2]*diipbv[2]);
                  pfac  = pfac  + coefbv[face_id][jsou][isou]*pipr;
                  pfacd = pfacd + cofbfv[face_id][jsou][isou]*pipr;
                }

                pi  = pvar[ii][isou];
                pia = pvara[ii][isou];

                pir  = pi/relaxp - (1.-relaxp)/relaxp*pia;

                flux = iconvp*(flui*pir + fluj*pfac - b_massflux[face_id]*pi)
                  + idiffp*b_visc[face_id]*pfacd;
                rhs[ii][isou] = rhs[ii][isou] - flux;

              } /* isou */

              /* Imposed convective flux */
            } else {
              /*-----------------
                X-Y-Z components, p=u, v, w */
              for (isou = 0; isou < 3; isou++) {

                pfac  = inc*cofacv[face_id][isou];
                pfacd = inc*cofafv[face_id][isou];

                /*coefu and cofuf are matrices */
                for (jsou = 0; jsou < 3; jsou++) {
                  pir  =   pvar[ii][jsou]/relaxp
                         - (1.-relaxp)/relaxp*pvara[ii][jsou];

                  pipr = pir +ircflp*(   gradv[ii][jsou][0]*diipbv[0]
                                       + gradv[ii][jsou][1]*diipbv[1]
                                       + gradv[ii][jsou][2]*diipbv[2]);
                  pfac  = pfac + cofbcv[face_id][jsou][isou]*pipr;
                  pfacd = pfacd + cofbfv[face_id][jsou][isou]*pipr;
                }

                pi  = pvar[ii][isou];
                pia = pvara[ii][isou];

                flux =   iconvp*( - b_massflux[face_id]*pi
                                  + pfac*icvfli[face_id])
                       + idiffp*b_visc[face_id]*pfacd;
                rhs[ii][isou] = rhs[ii][isou] - flux;

              } /* isou */

            }

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_b_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, isou, jsou, diipbv,        \
                                        flui, fluj, pfac, pfacd, pip, flux, pi) \
                   if(m->n_b_faces > CS_THR_MIN)
        for (t_id = 0; t_id < n_b_threads; t_id++) {
          for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
               face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
               face_id++) {

            ii = b_face_cells[face_id];

            for (jsou = 0; jsou < 3; jsou++) {
              diipbv[jsou] = diipb[face_id][jsou];
            }

            /* Computed convective flux */
            if (icvfli[face_id] == 0) {
              /* Remove decentering for coupled faces */
              if (ifaccp == 1 && bc_type[face_id] == CS_COUPLED) {
                flui = 0.0;
                fluj = b_massflux[face_id];
              } else {
                flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
                fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
              }

              /*-----------------
                X-Y-Z components, p=u, v, w */
              for (isou = 0; isou < 3; isou++) {

                pfac  = inc*coefav[face_id][isou];
                pfacd = inc*cofafv[face_id][isou];

                /*coefu and cofuf are matrices */
                for (jsou = 0; jsou < 3; jsou++) {
                  pip =   pvar[ii][jsou]
                        + ircflp*(  gradv[ii][jsou][0]*diipbv[0]
                                  + gradv[ii][jsou][1]*diipbv[1]
                                  + gradv[ii][jsou][2]*diipbv[2]);
                  pfac  = pfac  + coefbv[face_id][jsou][isou]*pip;
                  pfacd = pfacd + cofbfv[face_id][jsou][isou]*pip;
                }

                pi = pvar[ii][isou];

                flux = iconvp*((flui-b_massflux[face_id])*pi + fluj*pfac)
                  + idiffp*b_visc[face_id]*pfacd;
                rhs[ii][isou] = rhs[ii][isou] - thetap * flux;

              } /* isou */

              /* Imposed convective flux */
            } else {

              /*-----------------
                X-Y-Z components, p=u, v, w */
              for (isou = 0; isou < 3; isou++) {

                pfac = inc*cofacv[face_id][isou];
                pfacd = inc*cofafv[face_id][isou];

                /*coefu and cofuf are matrices */
                for (jsou = 0; jsou < 3; jsou++) {
                  pip = pvar[ii][jsou]
                        + ircflp*(  gradv[ii][jsou][0]*diipbv[0]
                                  + gradv[ii][jsou][1]*diipbv[1]
                                  + gradv[ii][jsou][2]*diipbv[2]);
                  pfac = pfac + cofbcv[face_id][jsou][isou]*pip;
                  pfacd = pfacd + cofbfv[face_id][jsou][isou]*pip;
                }

                pi = pvar[ii][isou];

                flux = iconvp*( -b_massflux[face_id]*pi + pfac*icvfli[face_id])
                  + idiffp*b_visc[face_id]*pfacd;
                rhs[ii][isou] = rhs[ii][isou] - thetap * flux;

              } /* isou */

            }


          }
        }
      }

    } /* idtvar */

  }

  /* 3. Computation of the transpose grad(vel) term and grad(-2/3 div(vel)) */

  if (ivisep == 1) {

    /* We do not know what condition to put in the inlets and the outlets, so we
       assume that there is an equilibrium. Moreover, cells containing a coupled
       are removed. */

    /* Allocate a temporary array */
    BFT_MALLOC(bndcel, n_cells_ext, cs_real_t);

#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++)
      bndcel[cell_id] = 1.;

#   pragma omp parallel for private(ityp) if(m->n_b_faces > CS_THR_MIN)
    for (face_id = 0; face_id < m->n_b_faces; face_id++) {
      ityp = bc_type[face_id];
      if (   ityp == CS_OUTLET
          || ityp == CS_INLET
          || (ifaccp == 1 && ityp == CS_COUPLED))
        bndcel[b_face_cells[face_id]] = 0.;
    }

    if (halo != NULL)
      cs_halo_sync_var(halo, halo_type, bndcel);

    /* ---> Interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd, secvis, \
                                      visco, grdtrv, tgrdfl, flux)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          pnd = weight[face_id];
          secvis = secvif[face_id];
          visco = i_visc[face_id];

          grdtrv =      pnd*(gradv[ii][0][0]+gradv[ii][1][1]+gradv[ii][2][2])
                 + (1.-pnd)*(gradv[jj][0][0]+gradv[jj][1][1]+gradv[jj][2][2]);

          /* We need to compute trans_grad(u).IJ which is equal to IJ.grad(u) */

          for (isou = 0; isou < 3; isou++) {

            tgrdfl = dijpf[face_id][0] * (      pnd*gradv[ii][0][isou]
                                         + (1.-pnd)*gradv[jj][0][isou])
                   + dijpf[face_id][1] * (      pnd*gradv[ii][1][isou]
                                         + (1.-pnd)*gradv[jj][1][isou])
                   + dijpf[face_id][2] * (      pnd*gradv[ii][2][isou]
                                         + (1.-pnd)*gradv[jj][2][isou]);

            flux = visco*tgrdfl + secvis*grdtrv*i_face_normal[face_id][isou];

            rhs[ii][isou] = rhs[ii][isou] + flux*bndcel[ii];
            rhs[jj][isou] = rhs[jj][isou] - flux*bndcel[jj];

          }

        }
      }
    }

    /* ---> Boundary FACES
       the whole flux term of the stress tensor is already taken into account
       (so, no corresponding term in forbr)
       TODO in theory we should take the normal component into account (the
       tangential one is modeled by the wall law) */

    /*Free memory */
    BFT_FREE(bndcel);

  }

  /* Free memory */
  BFT_FREE(gradva);
  BFT_FREE(gradv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 * equation of a scalar field \f$ \varia \f$ such as the temperature.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs + \sum_{\fij \in \Facei{\celli}}      \left(
 *        C_p\dot{m}_\ij \varia_\fij
 *      - \lambda_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * \f$ Rhs \f$ has already been initialized before calling bilsct!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$)
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_thermal(int                       idtvar,
                                int                       f_id,
                                const cs_var_cal_opt_t    var_cal_opt,
                                int                       inc,
                                int                       iccocg,
                                int                       ifaccp,
                                cs_real_t       *restrict pvar,
                                const cs_real_t *restrict pvara,
                                const cs_int_t            bc_type[],
                                const cs_real_t           coefap[],
                                const cs_real_t           coefbp[],
                                const cs_real_t           cofafp[],
                                const cs_real_t           cofbfp[],
                                const cs_real_t           i_massflux[],
                                const cs_real_t           b_massflux[],
                                const cs_real_t           i_visc[],
                                const cs_real_t           b_visc[],
                                const cs_real_t           xcpp[],
                                cs_real_t       *restrict rhs)
{
  const int iconvp = var_cal_opt.iconv ;
  const int idiffp = var_cal_opt.idiff ;
  const int nswrgp = var_cal_opt.nswrgr;
  const int imrgra = var_cal_opt.imrgra;
  const int imligp = var_cal_opt.imligr;
  const int ircflp = var_cal_opt.ircflu;
  const int ischcp = var_cal_opt.ischcv;
  const int isstpp = var_cal_opt.isstpc;
  const int iwarnp = var_cal_opt.iwarni;
  const double blencp = var_cal_opt.blencv;
  const double epsrgp = var_cal_opt.epsrgr;
  const double climgp = var_cal_opt.climgr;
  const double extrap = var_cal_opt.extrag;
  const double relaxp = var_cal_opt.relaxv;
  const double thetap = var_cal_opt.thetav;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* Local variables */

  char var_name[32];

  cs_gnum_t n_upwind;
  cs_lnum_t face_id, ii, jj, cell_id;
  int iupwin, g_id, t_id;
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  double pfac,pfacd,flui,fluj,flux,fluxi,fluxj;
  double difx,dify,difz,djfx,djfy,djfz;
  double pi, pj, pia, pja;
  double pif,pjf,pip,pjp,pir,pjr,pipr,pjpr;
  double pifri,pifrj,pjfri,pjfrj;
  double testi,testj,testij;
  double dpxf,dpyf,dpzf;
  double dcc, ddi, ddj, tesqck;
  double dijpfx, dijpfy, dijpfz;
  double diipfx, diipfy, diipfz;
  double djjpfx, djjpfy, djjpfz;
  double diipbx, diipby, diipbz;
  double pnd, distf, srfan;
  double pfac1, pfac2, pfac3, unsvol;

  cs_real_3_t *grad;
  cs_real_3_t *grdpa;
  cs_field_t *f;

  int f_track_slope_test_id = -1;
  cs_real_t *gweight = NULL;

  cs_real_t  *v_slope_test = _get_v_slope_test(f_id,  var_cal_opt);

  /*==========================================================================*/

  /* 1. Initialization */

  /* Allocate work arrays */
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(grdpa, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf
        (
         _(" %s: Convection in centered blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    } else {
      bft_printf
        (
         _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
  }

  iupwin = (blencp > 0.) ? 0 : 1;

  /* 2. Compute the balance with reconstruction technics */

  /* Compute the gradient of the variable

     The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test. */

  if (  (idiffp!=0 && ircflp==1)
        || (  iconvp!=0 && iupwin==0
              && (ischcp==0 || ircflp==1 || isstpp==0))) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && var_cal_opt.iwgrec == 1) {
        if (var_cal_opt.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
          }
        }
      }
    }
    cs_gradient_scalar(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       recompute_cocg,
                       nswrgp,
                       tr_dim,
                       0, /* hyd_p_flag */
                       iwarnp,
                       imligp,
                       epsrgp,
                       extrap,
                       climgp,
                       NULL, /* f_ext exterior force */
                       coefap,
                       coefbp,
                       pvar,
                       gweight, /* Weighted gradient */
                       grad);

  } else {
#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    }
  }

  /* ======================================================================
     ---> Compute uncentred gradient gradpa for the slope test
     ======================================================================*/

# pragma omp parallel for
  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    grdpa[cell_id][0] = 0.;
    grdpa[cell_id][1] = 0.;
    grdpa[cell_id][2] = 0.;
  }

  if (iconvp > 0 && iupwin == 0 && isstpp == 0) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, difx, dify, difz,  \
                                      djfx, djfy, djfz, pif, pjf, pfac,   \
                                      pfac1, pfac2, pfac3)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          difx = i_face_cog[face_id][0] - cell_cen[ii][0];
          dify = i_face_cog[face_id][1] - cell_cen[ii][1];
          difz = i_face_cog[face_id][2] - cell_cen[ii][2];
          djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
          djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
          djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

          pif = pvar[ii] + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
          pjf = pvar[jj] + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

          pfac = pjf;
          if (i_massflux[face_id] > 0.) pfac = pif;

          pfac1 = pfac*i_face_normal[face_id][0];
          pfac2 = pfac*i_face_normal[face_id][1];
          pfac3 = pfac*i_face_normal[face_id][2];

          grdpa[ii][0] = grdpa[ii][0] + pfac1;
          grdpa[ii][1] = grdpa[ii][1] + pfac2;
          grdpa[ii][2] = grdpa[ii][2] + pfac3;

          grdpa[jj][0] = grdpa[jj][0] - pfac1;
          grdpa[jj][1] = grdpa[jj][1] - pfac2;
          grdpa[jj][2] = grdpa[jj][2] - pfac3;

        }
      }
    }

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz, pfac) \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];
          diipbx = diipb[face_id][0];
          diipby = diipb[face_id][1];
          diipbz = diipb[face_id][2];
          pfac =   inc*coefap[face_id]
            + coefbp[face_id] * (  pvar[ii]          + diipbx*grad[ii][0]
                                   + diipby*grad[ii][1] + diipbz*grad[ii][2]);
          grdpa[ii][0] = grdpa[ii][0] + pfac*b_face_normal[face_id][0];
          grdpa[ii][1] = grdpa[ii][1] + pfac*b_face_normal[face_id][1];
          grdpa[ii][2] = grdpa[ii][2] + pfac*b_face_normal[face_id][2];

        }
      }
    }

#   pragma omp parallel for private(unsvol)
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      unsvol = 1./cell_vol[cell_id];
      grdpa[cell_id][0] = grdpa[cell_id][0]*unsvol;
      grdpa[cell_id][1] = grdpa[cell_id][1]*unsvol;
      grdpa[cell_id][2] = grdpa[cell_id][2]*unsvol;
    }

    /* Synchronization for parallelism or periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)grdpa, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)grdpa, 3);
    }

  }


  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  n_upwind = 0;

  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++)
      rhs[cell_id] = 0.;
  }

  /* --> Pure upwind flux
     =====================*/

  if (iupwin == 1) {

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,     \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf,        \
                                        pip, pjp, pipr, pjpr, flui, fluj,        \
                                        pif, pjf, fluxi, fluxj,                  \
                                        pi, pj, pir, pjr, pia, pja)              \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];
            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells) {
              n_upwind = n_upwind+1;
            }

            pi = pvar[ii];
            pj = pvar[jj];
            pia = pvara[ii];
            pja = pvara[jj];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            /* reconstruction only if IRCFLP = 1 */
            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            pir = pi/relaxp - (1.-relaxp)/relaxp * pia;
            pjr = pj/relaxp - (1.-relaxp)/relaxp * pja;

            pipr = pir
              + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjpr = pjr
              + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            pif  = pi;
            pjf  = pj;

            fluxi =   iconvp*xcpp[ii]*(  flui*pir + fluj*pjf
                                       - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pipr - pjp);
            fluxj =   iconvp*xcpp[jj]*(flui*pif + fluj*pjr
                                       - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip - pjpr);

            rhs[ii] = rhs[ii] - fluxi;
            rhs[jj] = rhs[jj] + fluxj;

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,     \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf,        \
                                        pip, pjp, flui, fluj, pif, pjf,          \
                                        fluxi, fluxj, pi, pj)                    \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];
            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells) {
              n_upwind = n_upwind+1;
            }

            pi = pvar[ii];
            pj = pvar[jj];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            pif = pi;
            pjf = pj;

            fluxi =   iconvp*xcpp[ii]*(  flui*pif +fluj*pjf
                                       - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pip -pjp);
            fluxj =   iconvp*xcpp[jj]*(  flui*pif +fluj*pjf
                                       - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip -pjp);

            rhs[ii] = rhs[ii] - thetap * fluxi;
            rhs[jj] = rhs[jj] + thetap * fluxj;

          }
        }
      }

    }


    /* --> Flux with no slope test
       ============================*/

  } else if (isstpp==1) {

    if (ischcp<0 || ischcp>1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,     \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf,        \
                                        pip, pjp, pipr, pjpr, flui,              \
                                        fluj, pir, pjr, pifri, pjfri,            \
                                        pifrj, pjfrj,                            \
                                        difx, dify, difz, djfx, djfy, djfz,      \
                                        fluxi, fluxj, pi, pj, pia, pja)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            pi = pvar[ii];
            pj = pvar[jj];
            pia = pvara[ii];
            pja = pvara[jj];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =  i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            pir = pi/relaxp - (1. - relaxp)/relaxp*pia;
            pjr = pj/relaxp - (1. - relaxp)/relaxp*pja;

            pipr = pir
              + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjpr = pjr
              + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));


            /* Centered
               --------*/

            if (ischcp == 1) {

              pifri = pnd*pipr +(1.-pnd)*pjp;
              pjfri = pifri;
              pifrj = pnd*pip  +(1.-pnd)*pjpr;
              pjfrj = pifrj;

              /* Second order
                 ------------*/

            } else { /* if (ischcp.eq.0) then */

              difx = i_face_cog[face_id][0] - cell_cen[ii][0];
              dify = i_face_cog[face_id][1] - cell_cen[ii][1];
              difz = i_face_cog[face_id][2] - cell_cen[ii][2];
              djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
              djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
              djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

              /* leave reconstruction of PIF and PJF even if IRCFLP=0
                 otherwise, it is the same as using upwind */
              pifri = pir + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
              pifrj = pi + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
              pjfrj = pjr + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];
              pjfri = pj + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

            }

            /* Blending
               --------*/

            pifri = blencp*pifri+(1.-blencp)*pir;
            pifrj = blencp*pifrj+(1.-blencp)*pi;
            pjfri = blencp*pjfri+(1.-blencp)*pj;
            pjfrj = blencp*pjfrj+(1.-blencp)*pjr;

            /* Flux
               ----*/

            fluxi =   iconvp*xcpp[ii]*(  flui*pifri + fluj*pjfri
                                       - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pipr -pjp);
            fluxj =   iconvp*xcpp[jj]*(  flui*pifrj + fluj*pjfrj
                                       - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip -pjpr);

            /* Assembly
               --------*/

            rhs[ii] = rhs[ii] - fluxi;
            rhs[jj] = rhs[jj] + fluxj;

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, diipfx, diipfy, diipfz, djjpfx,     \
                                        djjpfy, djjpfz, dpxf, dpyf, dpzf,        \
                                        pip, pjp, flui, fluj, pif, pjf,          \
                                        difx, dify, difz, djfx, djfy, djfz,      \
                                        fluxi, fluxj, pi, pj)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd = weight[face_id];

            pi = pvar[ii];
            pj = pvar[jj];

            /* Recompute II' and JJ' at this level */

            diipfx = i_face_cog[face_id][0] - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy = i_face_cog[face_id][1] - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz = i_face_cog[face_id][2] - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] + fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] - fabs(i_massflux[face_id]));

            /* Centered
               --------*/

            if (ischcp == 1) {

              pif = pnd*pip +(1.-pnd)*pjp;
              pjf = pif;

              /* Second order
                 ------------*/

            } else { /* if (ischcp.eq.0) then */

              difx = i_face_cog[face_id][0] - cell_cen[ii][0];
              dify = i_face_cog[face_id][1] - cell_cen[ii][1];
              difz = i_face_cog[face_id][2] - cell_cen[ii][2];
              djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
              djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
              djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

              /* leave reconstruction of PIF and PJF even if IRCFLP=0
                 otherwise, it is the same as using upwind */
              pif = pi + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
              pjf = pj + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

            }

            /* Blending
               --------*/

            pif = blencp*pif+(1.-blencp)*pi;
            pjf = blencp*pjf+(1.-blencp)*pj;

            /* Flux
               ----*/

            fluxi =   iconvp*xcpp[ii]*(  flui*pif +fluj*pjf
                                      - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pip -pjp);
            fluxj =   iconvp*xcpp[jj]*(  flui*pif +fluj*pjf
                                       - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip -pjp);

            /* Assembly
               --------*/

            rhs[ii] = rhs[ii] - thetap * fluxi;
            rhs[jj] = rhs[jj] + thetap * fluxj;

          }
        }
      }

    }

    /* --> Flux with slope test
       =========================*/

  } else {

    if (ischcp<0 || ischcp>1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    /* Steady */
    if (idtvar < 0) {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy,        \
                                        dijpfz, pnd, distf, srfan, diipfx,      \
                                        diipfy, diipfz, djjpfx, djjpfy, djjpfz, \
                                        dpxf, dpyf, dpzf, pip, pjp, pipr, pjpr, \
                                        flui, fluj, pir, pjr, testi, testj,     \
                                        testij, dcc, ddi, ddj, tesqck, pifri,   \
                                        pjfri, pifrj, pjfrj, difx, dify, difz,  \
                                        djfx, djfy, djfz, fluxi, fluxj, pi, pj, \
                                        pia, pja)                               \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd    = weight[face_id];
            distf  = i_dist[face_id];
            srfan  = i_face_surf[face_id];

            pi = pvar[ii];
            pj = pvar[jj];
            pia = pvara[ii];
            pja = pvara[jj];

            /* Recompute II' and JJ' at this level */
            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            pir = pi/relaxp - (1. - relaxp)/relaxp*pia;
            pjr = pj/relaxp - (1. - relaxp)/relaxp*pja;

            pipr = pir
              + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjpr = pjr
              + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));


            /* Slope test
               ----------*/

            testi =   grdpa[ii][0]*i_face_normal[face_id][0]
                    + grdpa[ii][1]*i_face_normal[face_id][1]
                    + grdpa[ii][2]*i_face_normal[face_id][2];
            testj =   grdpa[jj][0]*i_face_normal[face_id][0]
                    + grdpa[jj][1]*i_face_normal[face_id][1]
                    + grdpa[jj][2]*i_face_normal[face_id][2];
            testij=   grdpa[ii][0]*grdpa[jj][0]
                    + grdpa[ii][1]*grdpa[jj][1]
                    + grdpa[ii][2]*grdpa[jj][2];

            if (i_massflux[face_id] > 0.) {
              dcc =   grad[ii][0]*i_face_normal[face_id][0]
                    + grad[ii][1]*i_face_normal[face_id][1]
                    + grad[ii][2]*i_face_normal[face_id][2];
              ddi = testi;
              ddj = (pj-pi)/distf *srfan;
            } else {
              dcc =   grad[jj][0]*i_face_normal[face_id][0]
                    + grad[jj][1]*i_face_normal[face_id][1]
                    + grad[jj][2]*i_face_normal[face_id][2];
              ddi = (pj-pi)/distf *srfan;
              ddj = testj;
            }
            tesqck = pow(dcc,2) -pow((ddi-ddj),2);

            /* Upwind
               ------*/

            if (tesqck <= 0. || testij <= 0.) {

              pifri = pir;
              pifrj = pi;
              pjfri = pj;
              pjfrj = pjr;
              /* in parallel, face will be counted by one and only one rank */
              if (ii < n_cells)
                n_upwind = n_upwind+1;
              if (v_slope_test != NULL) {
                v_slope_test[ii] += fabs(i_massflux[face_id]) / cell_vol[ii];
                v_slope_test[jj] += fabs(i_massflux[face_id]) / cell_vol[jj];
              }

            } else {

              /* Centered
                 --------*/

              if (ischcp == 1) {

                pifri = pnd*pipr +(1.-pnd)*pjp;
                pjfri = pifri;
                pifrj = pnd*pip  +(1.-pnd)*pjpr;
                pjfrj = pifrj;

                /* Second order
                   ------------*/

              } else { /* if (ischcp.eq.0) then */

                difx = i_face_cog[face_id][0] - cell_cen[ii][0];
                dify = i_face_cog[face_id][1] - cell_cen[ii][1];
                difz = i_face_cog[face_id][2] - cell_cen[ii][2];
                djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
                djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
                djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

                /* leave reconstruction of PIF and PJF even if IRCFLP=0
                   otherwise, it is the same as using upwind */
                pifri = pir + difx*grad[ii][0]
                            + dify*grad[ii][1]
                            + difz*grad[ii][2];
                pifrj = pi  + difx*grad[ii][0]
                            + dify*grad[ii][1]
                            + difz*grad[ii][2];
                pjfrj = pjr + djfx*grad[jj][0]
                            + djfy*grad[jj][1]
                            + djfz*grad[jj][2];
                pjfri = pj  + djfx*grad[jj][0]
                            + djfy*grad[jj][1]
                            + djfz*grad[jj][2];

              }

            }

            /* Blending
               --------*/

            pifri = blencp*pifri+(1.-blencp)*pir;
            pifrj = blencp*pifrj+(1.-blencp)*pi;
            pjfri = blencp*pjfri+(1.-blencp)*pj;
            pjfrj = blencp*pjfrj+(1.-blencp)*pjr;

            /* Flux
               ----*/

            fluxi =   iconvp*xcpp[ii]*(  flui*pifri + fluj*pjfri
                                       - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pipr -pjp);
            fluxj =   iconvp*xcpp[jj]*(  flui*pifrj + fluj*pjfrj
                                       - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip -pjpr);

            /* Assembly
               --------*/

            rhs[ii] = rhs[ii] - fluxi;
            rhs[jj] = rhs[jj] + fluxj;

          }
        }
      }

      /* Unsteady */
    } else {

      for (g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for private(face_id, ii, jj, dijpfx, dijpfy, dijpfz, \
                                        pnd, distf, srfan, diipfx, diipfy,       \
                                        diipfz, djjpfx, djjpfy, djjpfz,          \
                                        dpxf, dpyf, dpzf, pip, pjp, flui, fluj,  \
                                        testi, testj, testij, dcc, ddi, ddj,     \
                                        tesqck, pif, pjf, difx, dify, difz,      \
                                        djfx, djfy, djfz, fluxi, fluxj, pi, pj)  \
                   reduction(+:n_upwind)
        for (t_id = 0; t_id < n_i_threads; t_id++) {
          for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
               face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
               face_id++) {

            ii = i_face_cells[face_id][0];
            jj = i_face_cells[face_id][1];

            dijpfx = dijpf[face_id][0];
            dijpfy = dijpf[face_id][1];
            dijpfz = dijpf[face_id][2];

            pnd    = weight[face_id];
            distf  = i_dist[face_id];
            srfan  = i_face_surf[face_id];

            pi = pvar[ii];
            pj = pvar[jj];

            /* Recompute II' and JJ' at this level */

            diipfx =   i_face_cog[face_id][0]
                     - (cell_cen[ii][0] + (1.-pnd) * dijpfx);
            diipfy =   i_face_cog[face_id][1]
                     - (cell_cen[ii][1] + (1.-pnd) * dijpfy);
            diipfz =   i_face_cog[face_id][2]
                     - (cell_cen[ii][2] + (1.-pnd) * dijpfz);
            djjpfx = i_face_cog[face_id][0] -  cell_cen[jj][0] + pnd * dijpfx;
            djjpfy = i_face_cog[face_id][1] -  cell_cen[jj][1] + pnd * dijpfy;
            djjpfz = i_face_cog[face_id][2] -  cell_cen[jj][2] + pnd * dijpfz;

            dpxf = 0.5*(grad[ii][0] + grad[jj][0]);
            dpyf = 0.5*(grad[ii][1] + grad[jj][1]);
            dpzf = 0.5*(grad[ii][2] + grad[jj][2]);

            pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz);
            pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz);

            flui = 0.5*(i_massflux[face_id] +fabs(i_massflux[face_id]));
            fluj = 0.5*(i_massflux[face_id] -fabs(i_massflux[face_id]));

            /* Slope test
               ----------*/

            testi =   grdpa[ii][0]*i_face_normal[face_id][0]
                    + grdpa[ii][1]*i_face_normal[face_id][1]
                    + grdpa[ii][2]*i_face_normal[face_id][2];
            testj =   grdpa[jj][0]*i_face_normal[face_id][0]
                    + grdpa[jj][1]*i_face_normal[face_id][1]
                    + grdpa[jj][2]*i_face_normal[face_id][2];
            testij =   grdpa[ii][0]*grdpa[jj][0]
                     + grdpa[ii][1]*grdpa[jj][1]
                     + grdpa[ii][2]*grdpa[jj][2];

            if (i_massflux[face_id] > 0.) {
              dcc =   grad[ii][0]*i_face_normal[face_id][0]
                    + grad[ii][1]*i_face_normal[face_id][1]
                    + grad[ii][2]*i_face_normal[face_id][2];
              ddi = testi;
              ddj = (pj-pi)/distf *srfan;
            } else {
              dcc =   grad[jj][0]*i_face_normal[face_id][0]
                    + grad[jj][1]*i_face_normal[face_id][1]
                    + grad[jj][2]*i_face_normal[face_id][2];
              ddi = (pj-pi)/distf *srfan;
              ddj = testj;
            }
            tesqck = pow(dcc,2) -pow((ddi-ddj),2);

            /* Upwind
               ------*/

            if (tesqck <= 0. || testij <= 0.) {

              pif = pi;
              pjf = pj;
              /* in parallel, face will be counted by one and only one rank */
              if (ii < n_cells)
                n_upwind = n_upwind+1;
              if (v_slope_test != NULL) {
                v_slope_test[ii] += fabs(i_massflux[face_id]) / cell_vol[ii];
                v_slope_test[jj] += fabs(i_massflux[face_id]) / cell_vol[jj];
              }

            } else {

              /* Centered
                 --------*/

              if (ischcp == 1) {

                pif = pnd*pip +(1.-pnd)*pjp;
                pjf = pif;

                /* Second order
                   ------------*/

              } else { /* if (ischcp.eq.0) then */

                difx = i_face_cog[face_id][0] - cell_cen[ii][0];
                dify = i_face_cog[face_id][1] - cell_cen[ii][1];
                difz = i_face_cog[face_id][2] - cell_cen[ii][2];
                djfx = i_face_cog[face_id][0] - cell_cen[jj][0];
                djfy = i_face_cog[face_id][1] - cell_cen[jj][1];
                djfz = i_face_cog[face_id][2] - cell_cen[jj][2];

                /* leave reconstruction of PIF and PJF even if IRCFLP=0
                   otherwise, it is the same as using upwind */
                pif = pi + difx*grad[ii][0]+dify*grad[ii][1]+difz*grad[ii][2];
                pjf = pj + djfx*grad[jj][0]+djfy*grad[jj][1]+djfz*grad[jj][2];

              }

            }

            /* Blending
               --------*/

            pif = blencp*pif+(1.-blencp)*pi;
            pjf = blencp*pjf+(1.-blencp)*pj;

            /* Flux
               ----*/

            fluxi =   iconvp*xcpp[ii]*(  flui*pif + fluj*pjf
                                       - i_massflux[face_id]*pi)
                    + idiffp*i_visc[face_id]*(pip-pjp);
            fluxj =   iconvp*xcpp[jj]*(  flui*pif + fluj*pjf
                                       - i_massflux[face_id]*pj)
                    + idiffp*i_visc[face_id]*(pip-pjp);

            /* Assembly
               --------*/

            rhs[ii] = rhs[ii] - thetap *fluxi;
            rhs[jj] = rhs[jj] + thetap *fluxj;

          }
        }
      }

    } /* idtvar */

  } /* iupwin */


  if (iwarnp >= 2) {

    /* Sum number of clippings */
    cs_parall_counter(&n_upwind, 1);

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces \n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_faces);
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz,     \
                                      flui, fluj, pir, pipr, pfac, pfacd,      \
                                      flux, pi, pia)                           \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          pi = pvar[ii];
          pia = pvara[ii];

          diipbx = diipb[face_id][0];
          diipby = diipb[face_id][1];
          diipbz = diipb[face_id][2];

          /* Remove decentering for coupled faces */
          if (ifaccp == 1 && bc_type[face_id] == CS_COUPLED) {
            flui = 0.0;
            fluj = b_massflux[face_id];
          } else {
            flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
            fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
          }

          pir  = pi/relaxp - (1.-relaxp)/relaxp*pia;
          pipr = pir
            + ircflp*(grad[ii][0]*diipbx+grad[ii][1]*diipby+grad[ii][2]*diipbz);

          pfac  = inc*coefap[face_id] +coefbp[face_id]*pipr;
          pfacd = inc*cofafp[face_id] +cofbfp[face_id]*pipr;

          flux = iconvp*xcpp[ii]*(flui*pir + fluj*pfac - b_massflux[face_id]*pi )
            + idiffp*b_visc[face_id]*pfacd;
          rhs[ii] = rhs[ii] - flux;

        }
      }
    }

    /* Unsteady */
  } else {

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, diipbx, diipby, diipbz,     \
                                      flui, fluj, pip, pfac, pfacd, flux, pi)  \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          pi = pvar[ii];

          diipbx = diipb[face_id][0];
          diipby = diipb[face_id][1];
          diipbz = diipb[face_id][2];

          /* Remove decentering for coupled faces */
          if (ifaccp == 1 && bc_type[face_id] == CS_COUPLED) {
            flui = 0.0;
            fluj = b_massflux[face_id];
          } else {
            flui = 0.5*(b_massflux[face_id] +fabs(b_massflux[face_id]));
            fluj = 0.5*(b_massflux[face_id] -fabs(b_massflux[face_id]));
          }

          pip = pi
            + ircflp*(grad[ii][0]*diipbx+grad[ii][1]*diipby+grad[ii][2]*diipbz);

          pfac  = inc*coefap[face_id] + coefbp[face_id]*pip;
          pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pip;

          flux = iconvp*xcpp[ii]*((flui - b_massflux[face_id])*pi + fluj*pfac)
            + idiffp*b_visc[face_id]*pfacd;
          rhs[ii] = rhs[ii] - thetap * flux;

        }
      }
    }

  }

  /* Free memory */
  BFT_FREE(grad);
  BFT_FREE(grdpa);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before
 *   calling cs_anisotropic_diffusion_scalar!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
                                (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_scalar(int                       idtvar,
                                int                       f_id,
                                const cs_var_cal_opt_t    var_cal_opt,
                                int                       inc,
                                int                       iccocg,
                                cs_real_t       *restrict pvar,
                                const cs_real_t *restrict pvara,
                                const cs_real_t           coefap[],
                                const cs_real_t           coefbp[],
                                const cs_real_t           cofafp[],
                                const cs_real_t           cofbfp[],
                                const cs_real_t           i_visc[],
                                const cs_real_t           b_visc[],
                                cs_real_6_t     *restrict viscel,
                                const cs_real_2_t         weighf[],
                                const cs_real_t           weighb[],
                                cs_real_t       *restrict rhs)
{
  const int nswrgp = var_cal_opt.nswrgr;
  const int imrgra = var_cal_opt.imrgra;
  const int imligp = var_cal_opt.imligr;
  const int ircflp = var_cal_opt.ircflu;
  const int iwarnp = var_cal_opt.iwarni;
  const double epsrgp = var_cal_opt.epsrgr;
  const double climgp = var_cal_opt.climgr;
  const double extrap = var_cal_opt.extrag;
  const double relaxp = var_cal_opt.relaxv;
  const double thetap = var_cal_opt.thetav;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  /* Local variables */

  char var_name[32];

  cs_lnum_t face_id, ii, jj, cell_id;
  int isou, n_upwind, g_id, t_id;
  double pfacd,flux,fluxi,fluxj;
  double pi, pj, pia, pja;
  double pir,pjr,pippr,pjppr;
  double diippf[3], djjppf[3], pipp, pjpp;
  double visci[3][3], viscj[3][3];
  double fikdvi, fjkdvi;
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  cs_real_6_t *viscce;
  cs_real_6_t *w2;
  cs_real_3_t *grad;

  cs_field_t *f;

  cs_real_t *gweight = NULL;

  /* 1. Initialization */

  viscce = NULL;
  w2 = NULL;

  /* Allocate work arrays */
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  /* Porosity fields */
  cs_field_t *fporo = CS_F_(poro);
  cs_field_t *ftporo = CS_F_(t_poro);

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (fporo != NULL) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Without porosity */
  if (porosi == NULL) {
    viscce = viscel;

    /* With porosity */
  } else if (porosi != NULL && porosf == NULL) {
    BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      for (isou = 0; isou < 6; isou++) {
        w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
      }
    }
    viscce = w2;

    /* With tensorial porosity */
  } else if (porosi != NULL && porosf != NULL) {
    BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      cs_math_sym_33_product(porosf[cell_id],
                             viscel[cell_id],
                             w2[cell_id]);
    }
    viscce = w2;
  }

  /* ---> Periodicity and parallelism treatment of symmetric tensors */
  if (halo != NULL) {
    cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)viscce, 6);
    if (m->n_init_perio > 0)
      cs_halo_perio_sync_var_sym_tens(halo, halo_type, (cs_real_t *)viscce);
  }

  /* 2. Compute the diffusive part with reconstruction technics */

  /* ======================================================================
     ---> Compute the gradient of the current variable if needed
     ======================================================================*/

  if (ircflp == 1) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && var_cal_opt.iwgrec == 1) {
        if (var_cal_opt.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
          }
        }
      }
    }
    cs_gradient_scalar(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       recompute_cocg,
                       nswrgp,
                       tr_dim,
                       0, /* hyd_p_flag */
                       iwarnp,
                       imligp,
                       epsrgp,
                       extrap,
                       climgp,
                       NULL, /* f_ext exterior force */
                       coefap,
                       coefbp,
                       pvar,
                       gweight, /* Weighted gradient */
                       grad);

  } else {
#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    }
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  n_upwind = 0;

  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      rhs[cell_id] = 0.;
    }
  }

  /* Steady */
  if (idtvar < 0) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, visci, viscj, \
                                      pipp, pjpp, pippr, pjppr,      \
                                      fluxi, fluxj, fikdvi,          \
                                      pi, pj, pir, pjr, pia, pja)    \
                 reduction(+:n_upwind)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];
          /* in parallel, face will be counted by one and only one rank */
          if (ii < n_cells) {
            n_upwind++;
          }

          pi = pvar[ii];
          pj = pvar[jj];
          pia = pvara[ii];
          pja = pvara[jj];

          /* Recompute II" and JJ"
             ----------------------*/

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
          fikdvi = weighf[face_id][0];

          /* II" = IF + FI" */
          for (int i = 0; i < 3; i++) {
            diippf[i] = i_face_cog[face_id][i]-cell_cen[ii][i]
                      - fikdvi*( visci[0][i]*i_face_normal[face_id][0]
                               + visci[1][i]*i_face_normal[face_id][1]
                               + visci[2][i]*i_face_normal[face_id][2] );
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
          fjkdvi = weighf[face_id][1];

          /* JJ" = JF + FJ" */
          for (int i = 0; i < 3; i++) {
            djjppf[i] = i_face_cog[face_id][i]-cell_cen[jj][i]
                      + fjkdvi*( viscj[0][i]*i_face_normal[face_id][0]
                               + viscj[1][i]*i_face_normal[face_id][1]
                               + viscj[2][i]*i_face_normal[face_id][2] );
          }

          /* p in I" and J" */
          pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                             + grad[ii][1]*diippf[1]
                             + grad[ii][2]*diippf[2]);
          pjpp = pj + ircflp*( grad[jj][0]*djjppf[0]
                             + grad[jj][1]*djjppf[1]
                             + grad[jj][2]*djjppf[2]);

          pir = pi/relaxp - (1.-relaxp)/relaxp * pia;
          pjr = pj/relaxp - (1.-relaxp)/relaxp * pja;

          /* pr in I" and J" */
          pippr = pir + ircflp*( grad[ii][0]*diippf[0]
                               + grad[ii][1]*diippf[1]
                               + grad[ii][2]*diippf[2]);
          pjppr = pjr + ircflp*( grad[jj][0]*djjppf[0]
                               + grad[jj][1]*djjppf[1]
                               + grad[jj][2]*djjppf[2]);


          fluxi = i_visc[face_id]*(pippr - pjpp);
          fluxj = i_visc[face_id]*(pipp - pjppr);

          rhs[ii] -= fluxi;
          rhs[jj] += fluxj;

        }
      }
    }

    /* Unsteady */
  } else {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj,              \
                                      visci, viscj, fikdvi, fjkdvi, \
                                      pipp, pjpp, diippf, djjppf,   \
                                      flux, pi, pj)                 \
                 reduction(+:n_upwind)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];
          /* in parallel, face will be counted by one and only one rank */
          if (ii < n_cells) {
            n_upwind++;
          }

          pi = pvar[ii];
          pj = pvar[jj];

          /* Recompute II" and JJ"
             ----------------------*/

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
          fikdvi = weighf[face_id][0];

          /* II" = IF + FI" */
          for (int i = 0; i < 3; i++) {
            diippf[i] = i_face_cog[face_id][i]-cell_cen[ii][i]
                      - fikdvi*( visci[0][i]*i_face_normal[face_id][0]
                               + visci[1][i]*i_face_normal[face_id][1]
                               + visci[2][i]*i_face_normal[face_id][2] );
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
          fjkdvi = weighf[face_id][1];

          /* JJ" = JF + FJ" */
          for (int i = 0; i < 3; i++) {
            djjppf[i] = i_face_cog[face_id][i]-cell_cen[jj][i]
                      + fjkdvi*( viscj[0][i]*i_face_normal[face_id][0]
                               + viscj[1][i]*i_face_normal[face_id][1]
                               + viscj[2][i]*i_face_normal[face_id][2] );
          }

          /* p in I" and J" */
          pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                             + grad[ii][1]*diippf[1]
                             + grad[ii][2]*diippf[2]);
          pjpp = pj + ircflp*( grad[jj][0]*djjppf[0]
                             + grad[jj][1]*djjppf[1]
                             + grad[jj][2]*djjppf[2]);

          flux = i_visc[face_id]*(pipp -pjpp);

          rhs[ii] -= thetap*flux;
          rhs[jj] += thetap*flux;

        }
      }
    }
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, visci, fikdvi,       \
                                      pir, pippr, pfacd, flux, pi, pia) \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          pi = pvar[ii];
          pia = pvara[ii];

          pir = pi/relaxp - (1.-relaxp)/relaxp*pia;

          /* Recompute II"
             --------------*/

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
          fikdvi = weighb[face_id];

          /* II" = IF + FI" */
          for (int i = 0; i < 3; i++) {
            diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                      - fikdvi*( visci[0][i]*b_face_normal[face_id][0]
                               + visci[1][i]*b_face_normal[face_id][1]
                               + visci[2][i]*b_face_normal[face_id][2] );
          }

          pippr = pir + ircflp*( grad[ii][0]*diippf[0]
                               + grad[ii][1]*diippf[1]
                               + grad[ii][2]*diippf[2]);

          pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pippr;

          flux = b_visc[face_id]*pfacd;
          rhs[ii] -= flux;

        }
      }
    }

    /* Unsteady */
  } else {

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, visci, fikdvi, \
                                      diippf, pipp, pfacd, flux, pi)      \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          pi = pvar[ii];

          /* Recompute II"
             --------------*/

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
          fikdvi = weighb[face_id];

          /* II" = IF + FI" */
          for (int i = 0; i < 3; i++) {
            diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                      - fikdvi*( visci[0][i]*b_face_normal[face_id][0]
                               + visci[1][i]*b_face_normal[face_id][1]
                               + visci[2][i]*b_face_normal[face_id][2]);
          }

          pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                             + grad[ii][1]*diippf[1]
                             + grad[ii][2]*diippf[2]);

          pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pipp;

          flux = b_visc[face_id]*pfacd;
          rhs[ii] -= thetap*flux;

        }
      }
    }

  }

  /* Free memory */
  BFT_FREE(grad);
  BFT_FREE(w2);
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensorial
 * diffusivity for a transport equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling diftnv!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     var_cal_opt   variable calculation options
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     coefav        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf\jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_vector(int                         idtvar,
                                int                         f_id,
                                const cs_var_cal_opt_t      var_cal_opt,
                                int                         inc,
                                int                         ifaccp,
                                int                         ivisep,
                                cs_real_3_t       *restrict pvar,
                                const cs_real_3_t *restrict pvara,
                                const cs_int_t              bc_type[],
                                const cs_real_3_t           coefav[],
                                const cs_real_33_t          coefbv[],
                                const cs_real_3_t           cofafv[],
                                const cs_real_33_t          cofbfv[],
                                const cs_real_33_t          i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t             secvif[],
                                cs_real_3_t       *restrict rhs)
{
  const int nswrgp = var_cal_opt.nswrgr;
  const int imrgra = var_cal_opt.imrgra;
  const int imligp = var_cal_opt.imligr;
  const int ircflp = var_cal_opt.ircflu;
  const int iwarnp = var_cal_opt.iwarni;
  const double epsrgp = var_cal_opt.epsrgr;
  const double climgp = var_cal_opt.climgr;
  const double relaxp = var_cal_opt.relaxv;
  const double thetap = var_cal_opt.thetav;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* Local variables */

  char var_name[32];

  cs_gnum_t n_upwind;
  cs_lnum_t face_id, ii, jj, cell_id;
  int g_id, t_id, isou, jsou;
  int ityp;
  int i, j, k;
  double pfacd, flux, fluxi, fluxj, pnd;
  double pi, pj, pia, pja;
  double pir, pipr[3], pjpr[3];
  double dpvf[3];
  double dijpfv[3];
  double diipfv[3], djjpfv[3], pip[3], pjp[3];
  double diipbv[3];
  double grdtrv, secvis;

  cs_real_33_t *gradv;
  cs_real_t *bndcel;

  cs_field_t *f;

  /* 1. Initialization */

  /* Allocate work arrays */
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  /* 2. Compute the diffusive part with reconstruction technics */

  /* Compute the gradient of the current variable if needed */

  if (ircflp == 1 || ivisep == 1) {

    cs_gradient_vector(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       nswrgp,
                       iwarnp,
                       imligp,
                       epsrgp,
                       climgp,
                       coefav,
                       coefbv,
                       pvar,
                       gradv);

  } else {
#   pragma omp parallel for private(isou, jsou)
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      for (isou = 0; isou < 3 ; isou++) {
        for (jsou = 0; jsou < 3; jsou++)
          gradv[cell_id][isou][jsou] = 0.;
      }
    }
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  n_upwind = 0;

  if (n_cells_ext > n_cells) {
#   pragma omp parallel for private(isou) if(n_cells_ext -n_cells > CS_THR_MIN)
    for (cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      for (isou = 0; isou < 3; isou++) {
        rhs[cell_id][isou] = 0.;
      }
    }
  }

  /* Steady */
  if (idtvar < 0) {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd,     \
                                      dijpfv, diipfv, djjpfv, dpvf, pi, pj, \
                                      pia, pja, pip, pjp, pipr, pjpr,       \
                                      fluxi, fluxj)                         \
                 reduction(+:n_upwind)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          /* in parallel, face will be counted by one and only one rank */
          if (ii < n_cells) {
            n_upwind++;
          }

          for (jsou = 0; jsou < 3; jsou++) {
            dijpfv[jsou] = dijpf[face_id][jsou];
          }

          pnd = weight[face_id];

          /* Recompute II' and JJ' at this level */
          for (jsou = 0; jsou < 3; jsou++) {
            diipfv[jsou] =   i_face_cog[face_id][jsou]
                           - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
            djjpfv[jsou] =   i_face_cog[face_id][jsou]
                           - cell_cen[jj][jsou]  + pnd  * dijpfv[jsou];
          }

          /*-----------------
            X-Y-Z components, p=u, v, w */
          for (isou = 0; isou < 3; isou++) {

            for (jsou = 0; jsou < 3; jsou++) {
              dpvf[jsou] = 0.5*(gradv[ii][isou][jsou] + gradv[jj][isou][jsou]);
            }

            pi  = pvar [ii][isou];
            pj  = pvar [jj][isou];

            pia = pvara[ii][isou];
            pja = pvara[jj][isou];

            /* reconstruction only if IRCFLP = 1 */
            pip[isou] = pi + ircflp*(  dpvf[0]*diipfv[0]
                                     + dpvf[1]*diipfv[1]
                                     + dpvf[2]*diipfv[2]);
            pjp[isou] = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                     + dpvf[1]*djjpfv[1]
                                     + dpvf[2]*djjpfv[2]);

            pipr[isou] = pi /relaxp - (1.-relaxp)/relaxp * pia
                         + ircflp*(  dpvf[0]*diipfv[0]
                                   + dpvf[1]*diipfv[1]
                                   + dpvf[2]*diipfv[2]);
            pjpr[isou] = pj /relaxp - (1.-relaxp)/relaxp * pja
                         + ircflp*(  dpvf[0]*djjpfv[0]
                                   + dpvf[1]*djjpfv[1]
                                   + dpvf[2]*djjpfv[2]);

          }

          for (isou = 0; isou < 3; isou++) {

            fluxi = i_visc[face_id][0][isou]*(pipr[0] - pjp[0])
                  + i_visc[face_id][1][isou]*(pipr[1] - pjp[1])
                  + i_visc[face_id][2][isou]*(pipr[2] - pjp[2]);
            fluxj = i_visc[face_id][0][isou]*(pip[0] - pjpr[0])
                  + i_visc[face_id][1][isou]*(pip[1] - pjpr[1])
                  + i_visc[face_id][2][isou]*(pip[2] - pjpr[2]);

            rhs[ii][isou] = rhs[ii][isou] - fluxi;
            rhs[jj][isou] = rhs[jj][isou] + fluxj;

          }

        }
      }
    }

    /* Unsteady */
  } else {

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd, dijpfv,\
                                      diipfv, djjpfv, dpvf, pi, pj,            \
                                      pip, pjp, flux)                          \
                 reduction(+:n_upwind)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          /* in parallel, face will be counted by one and only one rank */
          if (ii < n_cells) {
            n_upwind++;
          }

          for (jsou = 0; jsou < 3; jsou++) {
            dijpfv[jsou] = dijpf[face_id][jsou];
          }

          pnd = weight[face_id];

          /* Recompute II' and JJ' at this level */
          for (jsou = 0; jsou < 3; jsou++) {
            diipfv[jsou] =   i_face_cog[face_id][jsou]
                           - (cell_cen[ii][jsou] + (1.-pnd) * dijpfv[jsou]);
            djjpfv[jsou] =   i_face_cog[face_id][jsou]
                           - cell_cen[jj][jsou]  +  pnd * dijpfv[jsou];
          }

          /*-----------------
            X-Y-Z components, p=u, v, w */
          for (isou = 0; isou < 3; isou++) {

            for (jsou = 0; jsou < 3; jsou++) {
              dpvf[jsou] = 0.5*(gradv[ii][isou][jsou] + gradv[jj][isou][jsou]);
            }

            pi = pvar[ii][isou];
            pj = pvar[jj][isou];

            pip[isou] = pi + ircflp*(  dpvf[0]*diipfv[0]
                                     + dpvf[1]*diipfv[1]
                                     + dpvf[2]*diipfv[2]);
            pjp[isou] = pj + ircflp*(  dpvf[0]*djjpfv[0]
                                     + dpvf[1]*djjpfv[1]
                                     + dpvf[2]*djjpfv[2]);

          }

          for (isou = 0; isou < 3; isou++) {

            flux = i_visc[face_id][0][isou]*(pip[0] - pjp[0])
                 + i_visc[face_id][1][isou]*(pip[1] - pjp[1])
                 + i_visc[face_id][2][isou]*(pip[2] - pjp[2]);

            rhs[ii][isou] = rhs[ii][isou] - thetap*flux;
            rhs[jj][isou] = rhs[jj][isou] + thetap*flux;

          }

        }
      }
    }

  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, isou, jsou, diipbv, \
                                      pfacd, pir, pipr, flux)          \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          for (jsou = 0; jsou < 3; jsou++)
            diipbv[jsou] = diipb[face_id][jsou];

          /*-----------------
            X-Y-Z components, p=u, v, w */
          for (isou = 0; isou < 3; isou++) {

            pfacd = inc*cofafv[face_id][isou];

            /*coefu and cofuf are matrices */
            for (jsou = 0; jsou < 3; jsou++) {
              pir  = pvar[ii][jsou]/relaxp - (1.-relaxp)/relaxp*pvara[ii][jsou];

              pipr[jsou] = pir +ircflp*(  gradv[ii][jsou][0]*diipbv[0]
                                        + gradv[ii][jsou][1]*diipbv[1]
                                        + gradv[ii][jsou][2]*diipbv[2]);
              pfacd = pfacd + cofbfv[face_id][jsou][isou]*pipr[jsou];
            }

            flux = b_visc[face_id]*pfacd;
            rhs[ii][isou] = rhs[ii][isou] - flux;

          } /* isou */

        }
      }
    }

    /* Unsteady */
  } else {

    for (g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, isou, jsou, diipbv, \
                                      pfacd, pir, flux, pi)            \
                 if(m->n_b_faces > CS_THR_MIN)
      for (t_id = 0; t_id < n_b_threads; t_id++) {
        for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          ii = b_face_cells[face_id];

          for (jsou = 0; jsou < 3; jsou++)
            diipbv[jsou] = diipb[face_id][jsou];

          /*-----------------
            X-Y-Z components, p=u, v, w */
          for (isou = 0; isou < 3; isou++) {

            pfacd = inc*cofafv[face_id][isou];

            /*coefu and cofuf are matrices */
            for (jsou = 0; jsou < 3; jsou++) {
              pir = pvar[ii][jsou] + ircflp*(  gradv[ii][jsou][0]*diipbv[0]
                                             + gradv[ii][jsou][1]*diipbv[1]
                                             + gradv[ii][jsou][2]*diipbv[2]);
              pfacd = pfacd + cofbfv[face_id][jsou][isou]*pir;
            }

            flux = b_visc[face_id]*pfacd;
            rhs[ii][isou] = rhs[ii][isou] - thetap * flux;

          } /* isou */

        }
      }
    }

  } /* idtvar */

  /* 3. Computation of the transpose grad(vel) term and grad(-2/3 div(vel)) */

  if (ivisep == 1) {

    /* We do not know what condition to put in the inlets and the outlets, so we
       assume that there is an equilibrium. Moreover, cells containing a coupled
       are removed. */

    /* Allocate a temporary array */
    BFT_MALLOC(bndcel, n_cells_ext, cs_real_t);

#   pragma omp parallel for
    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      bndcel[cell_id] = 1.;
    }

#   pragma omp parallel for private(ityp) if(m->n_b_faces > CS_THR_MIN)
    for (face_id = 0; face_id < m->n_b_faces; face_id++) {
      ityp = bc_type[face_id];
      if (   ityp == CS_OUTLET
          || ityp == CS_INLET
          || (ifaccp == 1 && ityp == CS_COUPLED)) {
        bndcel[b_face_cells[face_id]] = 0.;
      }
    }

    if (halo != NULL)
      cs_halo_sync_var(halo, halo_type, bndcel);

    /* ---> Interior faces */

    for (g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for private(face_id, ii, jj, isou, jsou, pnd, secvis,\
                                      grdtrv, flux, i, j, k)
      for (t_id = 0; t_id < n_i_threads; t_id++) {
        for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          ii = i_face_cells[face_id][0];
          jj = i_face_cells[face_id][1];

          pnd = weight[face_id];
          secvis = secvif[face_id];

          grdtrv =      pnd*(gradv[ii][0][0]+gradv[ii][1][1]+gradv[ii][2][2])
                 + (1.-pnd)*(gradv[jj][0][0]+gradv[jj][1][1]+gradv[jj][2][2]);

          for (i = 0; i < 3; i++) {

            flux = secvis*grdtrv*i_face_normal[ face_id][i];

            /* We need to compute (K grad(u)^T) .IJ
               which is equal to IJ . (grad(u) . K^T)
               But: (IJ . (grad(u) . K^T))_i = IJ_k grad(u)_kj K_ij
               Warning, in FORTRAN K_ij = K(j, i) */

            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                flux = flux + dijpf[face_id][k]
                            *(pnd*gradv[ii][k][j]+(1-pnd)*gradv[jj][k][j])
                            *i_visc[face_id][i][j];
              }
            }

            rhs[ii][i] = rhs[ii][i] + flux*bndcel[ii];
            rhs[jj][i] = rhs[jj][i] - flux*bndcel[jj];

          }

        }
      }
    }

    /* ---> Boundary FACES
       the whole flux term of the stress tensor is already taken into account
       (so, no corresponding term in forbr)
       TODO in theory we should take the normal component into account (the
       tangential one is modeled by the wall law) */

    /*Free memory */
    BFT_FREE(bndcel);

  }

  /* Free memory */
  BFT_FREE(gradv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the face mass flux with the face pressure (or pressure
 * increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
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
 *                               (Implicit part)
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
cs_face_diffusion_potential(const int                 f_id,
                            const cs_mesh_t          *m,
                            cs_mesh_quantities_t     *fvq,
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
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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

  cs_real_3_t *grad;
  cs_real_3_t *visel;
  cs_field_t *f;

  cs_real_t *gweight = NULL;

  /*==========================================================================*/

  /* i_visc and visel carry similar information */

  /*============================================================================
    1. Initialization
    ==========================================================================*/

  BFT_MALLOC(visel, n_cells_ext, cs_real_3_t);
  for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++) {
    visel[ii][0] = viselx[ii];
    visel[ii][1] = visely[ii];
    visel[ii][2] = viselz[ii];
  }

  if (init >= 1) {
#   pragma omp parallel for
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }
  } else if(init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  if (imrgra >= 10)
    imrgra = 10;
  else if (imrgra > 0)
    imrgra = 0;
  if (imrgra < 0)
    imrgra = -imrgra;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else {
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';
  }

  /* Handle parallelism and periodicity */

  if (halo != NULL)
    cs_halo_sync_var(halo, halo_type, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          i_massflux[face_id] += i_visc[face_id]*(pvar[ii] -pvar[jj]);

        }
      }
    }

    /* Mass flow through boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];
          double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

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
    if (f_id > -1) {
      /* Get the calculation option from the field */
      int key_cal_opt_id = cs_field_key_id("var_cal_opt");
      cs_var_cal_opt_t var_cal_opt;
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      if (f->type & CS_FIELD_VARIABLE && var_cal_opt.iwgrec == 1) {
        if (var_cal_opt.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
          }
        }
      }
    }
    else if (f_id == -2) {
      gweight = viselx;
    }
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
                       gweight, /* Weighted gradient */
                       grad);

    /* Handle parallelism and periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)visel, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)visel, 3);
    }

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double dpxf = 0.5*(  visel[ii][0]*grad[ii][0]
                             + visel[jj][0]*grad[jj][0]);
          double dpyf = 0.5*(  visel[ii][1]*grad[ii][1]
                             + visel[jj][1]*grad[jj][1]);
          double dpzf = 0.5*(  visel[ii][2]*grad[ii][2]
                             + visel[jj][2]*grad[jj][2]);

          double dijpfx = dijpf[face_id][0];
          double dijpfy = dijpf[face_id][1];
          double dijpfz = dijpf[face_id][2];

          /*---> Dij = IJ - (IJ.N) N */
          double dijx = (cell_cen[jj][0]-cell_cen[ii][0])-dijpfx;
          double dijy = (cell_cen[jj][1]-cell_cen[ii][1])-dijpfy;
          double dijz = (cell_cen[jj][2]-cell_cen[ii][2])-dijpfz;

          i_massflux[face_id] =  i_massflux[face_id]
                               + i_visc[face_id]*(pvar[ii] - pvar[jj])
                               +  (dpxf *dijx + dpyf*dijy + dpzf*dijz)
                                 * i_face_surf[face_id]/i_dist[face_id];

        }
      }
    }

    /* Mass flow through boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          double diipbx = diipb[face_id][0];
          double diipby = diipb[face_id][1];
          double diipbz = diipb[face_id][2];

          double pip = pvar[ii] + grad[ii][0]*diipbx
                                + grad[ii][1]*diipby
                                + grad[ii][2]*diipbz;
          double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pip;

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
/*!
 * \brief Add the explicit part of the pressure gradient term to the mass flux
 * in case of anisotropic diffusion of the pressure field \f$ P \f$.
 *
 * More precisely, the mass flux side \f$ \dot{m}_\fij \f$ is updated as
 * follows:
 * \f[
 * \dot{m}_\fij = \dot{m}_\fij -
 *              \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
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
                                    (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_diffusion_potential(const cs_mesh_t          *m,
                                        cs_mesh_quantities_t     *fvq,
                                        int                       init,
                                        int                       inc,
                                        int                       imrgra,
                                        int                       iccocg,
                                        int                       nswrgp,
                                        int                       imligp,
                                        int                       ircflp,
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
                                        cs_real_6_t     *restrict viscel,
                                        const cs_real_2_t         weighf[],
                                        const cs_real_t           weighb[],
                                        cs_real_t       *restrict i_massflux,
                                        cs_real_t       *restrict b_massflux)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  /* Local variables */

  char var_name[32];
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  double diippf[3], djjppf[3];
  double visci[3][3], viscj[3][3];

  cs_real_6_t *viscce;
  cs_real_6_t *w2;
  cs_real_3_t *grad;

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
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

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  if (imrgra >= 10)
    imrgra = 10;
  else if (imrgra > 0)
    imrgra = 0;
  if (imrgra < 0)
    imrgra = -imrgra;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  /* Porosity fields */

  cs_field_t *fporo = CS_F_(poro);
  cs_field_t *ftporo = CS_F_(t_poro);

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (fporo != NULL) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Handle parallelism and periodicity */

  if (halo != NULL)
    cs_halo_sync_var(halo, halo_type, pvar);


  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* ---> Contribution from interior faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      i_massflux[face_id] += i_visc[face_id]*(pvar[ii] - pvar[jj]);

    }

    /* ---> Contribution from boundary faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];
      double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

      b_massflux[face_id] += b_visc[face_id]*pfac;

    }

  }

  /*==========================================================================
    3. Update mass flux with reconstruction technics
    ==========================================================================*/

  if (nswrgp > 1) {

    viscce = NULL;
    w2 = NULL;

    /* Without porosity */
    if (porosi == NULL) {
      viscce = viscel;

      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
      BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 6; isou++) {
          w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
        }
      }
      viscce = w2;

      /* With tensorial porosity */
    } else if (porosi != NULL && porosf != NULL) {
      BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cs_math_sym_33_product(porosf[cell_id],
                               viscel[cell_id],
                               w2[cell_id]);
      }
      viscce = w2;
    }

    /* ---> Periodicity and parallelism treatment of symmetric tensors */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, (cs_real_t *)viscce, 6);

      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_sym_tens(halo,
                                        CS_HALO_STANDARD,
                                        (cs_real_t *)viscce);
    }

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

    /* Mass flow through interior faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      double pi = pvar[ii];
      double pj = pvar[jj];

      /* Recompute II" and JJ"
         ----------------------*/

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
      double fikdvi = weighf[face_id][0];

      /* II" = IF + FI" */
      for (int i = 0; i < 3; i++) {
        diippf[i] = i_face_cog[face_id][i]-cell_cen[ii][i]
                  - fikdvi*( visci[0][i]*i_face_normal[face_id][0]
                           + visci[1][i]*i_face_normal[face_id][1]
                           + visci[2][i]*i_face_normal[face_id][2] );
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
      double fjkdvi = weighf[face_id][1];

      /* JJ" = JF + FJ" */
      for (int i = 0; i < 3; i++) {
        djjppf[i] = i_face_cog[face_id][i]-cell_cen[jj][i]
                  + fjkdvi*( viscj[0][i]*i_face_normal[face_id][0]
                           + viscj[1][i]*i_face_normal[face_id][1]
                           + viscj[2][i]*i_face_normal[face_id][2] );
      }

      /* p in I" and J" */
      double pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                         + grad[ii][1]*diippf[1]
                         + grad[ii][2]*diippf[2]);
      double pjpp = pj + ircflp*( grad[jj][0]*djjppf[0]
                         + grad[jj][1]*djjppf[1]
                         + grad[jj][2]*djjppf[2]);

      i_massflux[face_id] += i_visc[face_id]*(pipp - pjpp);

    }

    /* ---> Contribution from boundary faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      double pi = pvar[ii];

      /* Recompute II"
         --------------*/

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
      double fikdvi = weighb[face_id];

      /* II" = IF + FI" */
      for (int i = 0; i < 3; i++) {
        diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                  - fikdvi*( visci[0][i]*b_face_normal[face_id][0]
                           + visci[1][i]*b_face_normal[face_id][1]
                           + visci[2][i]*b_face_normal[face_id][2] );
      }

      double pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                         + grad[ii][1]*diippf[1]
                         + grad[ii][2]*diippf[2]);


      double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pipp;

      b_massflux[face_id] += b_visc[face_id]*pfac;

    }

    /* Free memory */
    BFT_FREE(grad);
    BFT_FREE(w2);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the cell mass flux divergence with the face pressure (or
 * pressure increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \sum_j \Delta t \grad_\fij p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
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
 *                               (Implicit part)
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
cs_diffusion_potential(const int                 f_id,
                       const cs_mesh_t          *m,
                       cs_mesh_quantities_t     *fvq,
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
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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

  cs_real_3_t *grad;
  cs_real_3_t *visel;
  cs_field_t *f;

  const cs_real_t *gweight = NULL;

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  BFT_MALLOC(visel, n_cells_ext, cs_real_3_t);
  for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++) {
    visel[ii][0] = viselx[ii];
    visel[ii][1] = visely[ii];
    visel[ii][2] = viselz[ii];
  }

  if (init >= 1) {
#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++) {
      diverg[ii] = 0.;
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t ii = n_cells; ii < n_cells_ext; ii++) {
      diverg[ii] = 0.;
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  if (imrgra >= 10)
    imrgra = 10;
  else if (imrgra > 0)
    imrgra = 0;
  if (imrgra < 0)
    imrgra = -imrgra;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);
  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else {
    snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';
  }

  /* Handle parallelism and periodicity */

  if (halo != NULL)
    cs_halo_sync_var(halo, halo_type, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double i_massflux = i_visc[face_id]*(pvar[ii] -pvar[jj]);
          diverg[ii] += i_massflux;
          diverg[jj] -= i_massflux;

        }
      }
    }

    /* Mass flow through boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];
          double pfac = inc*cofafp[face_id] +cofbfp[face_id]*pvar[ii];

          double b_massflux = b_visc[face_id]*pfac;
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
    if (f_id != -1) {
      /* Get the calculation option from the field */
      int key_cal_opt_id = cs_field_key_id("var_cal_opt");
      cs_var_cal_opt_t var_cal_opt;
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      if (f->type & CS_FIELD_VARIABLE && var_cal_opt.iwgrec == 1) {
        if (var_cal_opt.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
          }
        }
      }
    }

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
                       gweight, /* Weighted gradient */
                       grad);

    /* Handle parallelism and periodicity */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)visel, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)visel, 3);
    }

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double dijpfx = dijpf[face_id][0];
          double dijpfy = dijpf[face_id][1];
          double dijpfz = dijpf[face_id][2];

          /*---> Dij = IJ - (IJ.N) N */
          double dijx = (cell_cen[jj][0]-cell_cen[ii][0])-dijpfx;
          double dijy = (cell_cen[jj][1]-cell_cen[ii][1])-dijpfy;
          double dijz = (cell_cen[jj][2]-cell_cen[ii][2])-dijpfz;

          double dpxf = 0.5*(  visel[ii][0]*grad[ii][0]
                             + visel[jj][0]*grad[jj][0]);
          double dpyf = 0.5*(  visel[ii][1]*grad[ii][1]
                             + visel[jj][1]*grad[jj][1]);
          double dpzf = 0.5*(  visel[ii][2]*grad[ii][2]
                             + visel[jj][2]*grad[jj][2]);

          double i_massflux =   i_visc[face_id]*(pvar[ii] -pvar[jj])
                              +  (dpxf*dijx + dpyf*dijy + dpzf*dijz)
                                *i_face_surf[face_id]/i_dist[face_id];

          diverg[ii] += i_massflux;
          diverg[jj] -= i_massflux;

        }
      }
    }

    /* Mass flow through boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          double diipbx = diipb[face_id][0];
          double diipby = diipb[face_id][1];
          double diipbz = diipb[face_id][2];

          double pip = pvar[ii] + grad[ii][0]*diipbx
                                + grad[ii][1]*diipby
                                + grad[ii][2]*diipbz;
          double pfac = inc*cofafp[face_id] +cofbfp[face_id]*pip;

          double b_massflux = b_visc[face_id]*pfac;
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
/*!
 * \brief Add the explicit part of the divergence of the mass flux due to the
 * pressure gradient (routine analog to diften).
 *
 * More precisely, the divergence of the mass flux side
 * \f$ \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij \f$ is updated as follows:
 * \f[
 * \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  = \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  - \sum_{\fij \in \Facei{\celli}}
 *    \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
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
                                     (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] diverg        divergence of the mass flux
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_potential(const cs_mesh_t          *m,
                                   cs_mesh_quantities_t     *fvq,
                                   int                       init,
                                   int                       inc,
                                   int                       imrgra,
                                   int                       iccocg,
                                   int                       nswrgp,
                                   int                       imligp,
                                   int                       ircflp,
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
                                   cs_real_6_t     *restrict viscel,
                                   const cs_real_2_t         weighf[],
                                   const cs_real_t           weighb[],
                                   cs_real_t       *restrict diverg)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
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
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  /* Local variables */

  char var_name[32];
  int tr_dim = 0;

  bool recompute_cocg = (iccocg) ? true : false;

  double diippf[3], djjppf[3];
  double visci[3][3], viscj[3][3];

  cs_real_6_t *viscce;
  cs_real_6_t *w2;
  cs_real_3_t *grad;

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++) {
      diverg[ii] = 0.;
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t ii = n_cells+0; ii < n_cells_ext; ii++) {
      diverg[ii] = 0.;
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  if (imrgra >= 10)
    imrgra = 10;
  else if (imrgra > 0)
    imrgra = 0;
  if (imrgra < 0)
    imrgra = -imrgra;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  snprintf(var_name, 31, "Var. 0"); var_name[31] = '\0';

  /* Porosity fields */

  cs_field_t *fporo = CS_F_(poro);
  cs_field_t *ftporo = CS_F_(t_poro);

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (fporo != NULL) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Handle parallelism and periodicity */

  if (halo != NULL)
    cs_halo_sync_var(halo, halo_type, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double flux = i_visc[face_id]*(pvar[ii] - pvar[jj]);
          diverg[ii] += flux;
          diverg[jj] -= flux;

        }
      }
    }

    /* Mass flow though boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];
          double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

          double flux = b_visc[face_id]*pfac;
          diverg[ii] += flux;

        }
      }
    }

  }

  /*==========================================================================
    3. Update mass flux with reconstruction technics
    ==========================================================================*/

  if (nswrgp > 1) {

    viscce = NULL;
    w2 = NULL;

    /* Without porosity */
    if (porosi == NULL) {
      viscce = viscel;

      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
      BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 6; isou++) {
          w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
        }
      }
      viscce = w2;

      /* With tensorial porosity */
    } else if (porosi != NULL && porosf != NULL) {
      BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cs_math_sym_33_product(porosf[cell_id],
                               viscel[cell_id],
                               w2[cell_id]);
      }
      viscce = w2;
    }

    /* ---> Periodicity and parallelism treatment of symmetric tensors */

    if (halo != NULL) {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, (cs_real_t *)viscce, 6);

      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_sym_tens(halo,
                                        CS_HALO_STANDARD,
                                        (cs_real_t *)viscce);
    }

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

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double pi = pvar[ii];
          double pj = pvar[jj];

          /* Recompute II" and JJ"
             ----------------------*/

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
          double fikdvi = weighf[face_id][0];

          /* II" = IF + FI" */
          for (int i = 0; i < 3; i++) {
            diippf[i] = i_face_cog[face_id][i]-cell_cen[ii][i]
                      - fikdvi*( visci[0][i]*i_face_normal[face_id][0]
                               + visci[1][i]*i_face_normal[face_id][1]
                               + visci[2][i]*i_face_normal[face_id][2] );
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
          double fjkdvi = weighf[face_id][1];

          /* JJ" = JF + FJ" */
          for (int i = 0; i < 3; i++) {
            djjppf[i] = i_face_cog[face_id][i]-cell_cen[jj][i]
                      + fjkdvi*( viscj[0][i]*i_face_normal[face_id][0]
                               + viscj[1][i]*i_face_normal[face_id][1]
                               + viscj[2][i]*i_face_normal[face_id][2] );
          }

          /* p in I" and J" */
          double pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                             + grad[ii][1]*diippf[1]
                             + grad[ii][2]*diippf[2]);
          double pjpp = pj + ircflp*( grad[jj][0]*djjppf[0]
                             + grad[jj][1]*djjppf[1]
                             + grad[jj][2]*djjppf[2]);

          double flux = i_visc[face_id]*(pipp - pjpp);

          diverg[ii] += flux;
          diverg[jj] -= flux;

        }
      }
    }

    /* Mass flow though boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
             face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          double pi = pvar[ii];

          /* Recompute II"
             --------------*/

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
          double fikdvi = weighb[face_id];

          /* II" = IF + FI" */
          for (int i = 0; i < 3; i++) {
            diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                      - fikdvi*( visci[0][i]*b_face_normal[face_id][0]
                               + visci[1][i]*b_face_normal[face_id][1]
                               + visci[2][i]*b_face_normal[face_id][2] );
          }

          double pipp = pi + ircflp*( grad[ii][0]*diippf[0]
                             + grad[ii][1]*diippf[1]
                             + grad[ii][2]*diippf[2]);


          double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pipp;

          double flux = b_visc[face_id]*pfac;
          diverg[ii] += flux;

        }
      }
    }

    /* Free memory */
    BFT_FREE(grad);
    BFT_FREE(w2);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
