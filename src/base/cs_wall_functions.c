/*============================================================================
 * Wall functions
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
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
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
#include "cs_turbulence_model.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_field_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_wall_functions.c
        Wall functions descriptor and computation.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_wall_functions_t

  \brief wall functions descriptor.

  Members of this wall functions descriptor are publicly accessible, to allow
  for concise syntax, as it is expected to be used in many places.

  \var  cs_wall_functions_t::iwallf
        Indicates the type of wall function used for the velocity
        boundary conditions on a frictional wall.\n
        - CS_WALL_F_DISABLED: no wall functions
        - CS_WALL_F_1SCALE_POWER: one scale of friction velocities (power law)
        - CS_WALL_F_1SCALE_LOG: one scale of friction velocities (log law)
        - CS_WALL_F_2SCALES_LOG: two scales of friction velocities (log law)
        - CS_WALL_F_SCALABLE_2SCALES_LOG: two scales of friction velocities
          (log law) (scalable wall functions)
        - CS_WALL_F_2SCALES_VDRIEST: two scales of friction velocities
          (mixing length based on V. Driest analysis)
        - CS_WALL_F_2SCALES_SMOOTH_ROUGH: wall function unifying rough and smooth
          friction regimes
        - CS_WALL_F_2SCALES_CONTINUOUS: All \f$ y^+ \f$  for low Reynolds models\n
        \ref iwallf is initialised to CS_WALL_F_1SCALE_LOG for \ref iturb = 10,
          40, 41 or 70
        (mixing length, LES and Spalart Allmaras).\n
        \ref iwallf is initialised to CS_WALL_F_DISABLED for \ref iturb = 0, 32,
          50 or 51\n
        \ref iwallf is initialised to CS_WALL_F_2SCALES_LOG for \ref iturb = 20,
          21, 30, 31 or 60
        (\f$k-\epsilon\f$, \f$R_{ij}-\epsilon\f$ LRR, \f$R_{ij}-\epsilon\f$ SSG
          and \f$k-\omega\f$ SST models).\n
        The v2f model (\ref iturb=50) is not designed to use wall functions
        (the mesh must be low Reynolds).\n
        The value \ref iwallf = CS_WALL_F_2SCALES_LOG is not compatible with
          \ref iturb=0, 10, 40
        or 41 (laminar, mixing length and LES).\n
        Concerning the \f$k-\epsilon\f$ and \f$R_{ij}-\epsilon\f$ models, the
        two-scales model is usually at least as satisfactory as the one-scale
        model.\n
        The scalable wall function allows to virtually shift the wall when
        necessary in order to be always in a logarithmic layer. It is used to make up for
        the problems related to the use of High-Reynolds models on very refined
        meshes.\n
        Useful if \ref iturb is different from 50.
  \var  cs_wall_functions_t::iwalfs
        wall functions for scalar
        - CS_WALL_F_S_ARPACI_LARSEN: three layers (Arpaci and Larsen) or two layers (Prandtl-Taylor) for
             Prandtl number smaller than 0.1
        - CS_WALL_F_S_VDRIEST: consistant with the 2 scales wall function for velocity using Van
             Driest mixing length
        - CS_WALL_F_S_LOUIS: default wall function for atmospheric flows for
             potential temperature. This has influence on the dynamic.
        - CS_WALL_F_S_MONIN_OBUKHOV: Monin Obukhov wall function for atmospheric flows for
             potential temperature. This has influence on the dynamic.

  \var  cs_wall_functions_t::ypluli
        limit value of \f$y^+\f$ for the viscous sublayer

        \ref ypluli depends on the chosen wall function: it is initialized to
        10.88 for the scalable wall function (\ref iwallf=CS_WALL_F_SCALABLE_2SCALES_LOG), otherwise it is
        initialized to \f$1/\kappa\approx 2,38\f$. In LES, \ref ypluli is taken
        by default to be 10.88. Always useful.

*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Wall functions structure and associated pointers

   ypluli is set to -grand*10. If the user has not changed this value,
   its value is modified in  modini (10.88 with invariant wall laws,
   1/kappa otherwise). */

static cs_wall_functions_t  _wall_functions =
{
  .iwallf = -999,
  .iwalfs = -999,
  .ypluli = -1e13
};

const cs_wall_functions_t  * cs_glob_wall_functions = &_wall_functions;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_wall_functions_get_pointers(int     **iwallf,
                                 int     **iwalfs,
                                 double  **ypluli);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the wall functions structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iwallf  --> pointer to cs_glob_wall_functions->iwallf
 *   iwalfs  --> pointer to cs_glob_wall_functions->iwalfs
 *   ypluli  --> pointer to cs_glob_wall_functions->ypluli
 *----------------------------------------------------------------------------*/

void
cs_f_wall_functions_get_pointers(int     **iwallf,
                                 int     **iwalfs,
                                 double  **ypluli)
{
  *iwallf  = (int *)&(_wall_functions.iwallf);
  *iwalfs = (int *)&(_wall_functions.iwalfs);
  *ypluli  = &(_wall_functions.ypluli);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_wall_functions_velocity
 *----------------------------------------------------------------------------*/

void CS_PROCF (wallfunctions, WALLFUNCTIONS)
(
 const int        *const iwallf,
 const cs_lnum_t  *const ifac,
 const cs_real_t  *const l_visc,
 const cs_real_t  *const t_visc,
 const cs_real_t  *const vel,
 const cs_real_t  *const y,
 const cs_real_t  *const rough_d,
 const cs_real_t  *const rnnb,
 const cs_real_t  *const kinetic_en,
       int              *iuntur,
       cs_lnum_t        *nsubla,
       cs_lnum_t        *nlogla,
       cs_real_t        *ustar,
       cs_real_t        *uk,
       cs_real_t        *yplus,
       cs_real_t        *ypup,
       cs_real_t        *cofimp,
       cs_real_t        *dplus
)
{
  assert(*iwallf >= 0 && *iwallf <= 7);

  cs_gnum_t nsubla_loc = *nsubla;
  cs_gnum_t nlogla_loc = *nlogla;

  /* Get the adjacent border cell and its fluid/solid tag */
  cs_lnum_t cell_id = cs_glob_mesh->b_face_cells[*ifac-1] ;

  /* If the cell is a solid cell, disable wall functions */

  cs_wall_f_type_t iwallf_loc = (cs_wall_f_type_t)*iwallf;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  if (mq->has_disable_flag) {
    if (mq->c_disable_flag[cell_id])
      iwallf_loc = CS_WALL_F_DISABLED;
  }

  cs_wall_functions_velocity(iwallf_loc,
                             *l_visc,
                             *t_visc,
                             *vel,
                             *y,
                             *rough_d,
                             *rnnb,
                             *kinetic_en,
                             iuntur,
                             &nsubla_loc,
                             &nlogla_loc,
                             ustar,
                             uk,
                             yplus,
                             ypup,
                             cofimp,
                             dplus);

  *nsubla = (cs_lnum_t)nsubla_loc;
  *nlogla = (cs_lnum_t)nlogla_loc;
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_wall_functions_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (hturbp, HTURBP)
(
 const int        *const iwalfs,
 const cs_real_t  *const l_visc,
 const cs_real_t  *const prl,
 const cs_real_t  *const prt,
 const cs_real_t  *const rough_t,
 const cs_real_t  *const uk,
 const cs_real_t  *const yplus,
 const cs_real_t  *const dplus,
       cs_real_t        *htur,
       cs_real_t        *yplim
)
{
  cs_wall_functions_scalar((cs_wall_f_s_type_t)*iwalfs,
                           *l_visc,
                           *prl,
                           *prt,
                           *rough_t,
                           *uk,
                           *yplus,
                           *dplus,
                           htur,
                           yplim);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *! \brief Provide access to cs_glob_wall_functions
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_wall_functions_t *
cs_get_glob_wall_functions(void)
{
  return &_wall_functions;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the friction velocity and \f$y^+\f$ / \f$u^+\f$.
 *
 * \param[in]     iwallf        wall function type
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     t_visc        turbulent kinematic viscosity
 * \param[in]     vel           wall projected cell center velocity
 * \param[in]     y             wall distance
 * \param[in]     rough_d       roughness length scale
 *                              (not sand grain roughness)
 * \param[in]     rnnb          \f$\vec{n}.(\tens{R}\vec{n})\f$
 * \param[in]     kinetic_en    turbulent kinetic energy (cell center)
 * \param[in,out] iuntur        indicator: 0 in the viscous sublayer
 * \param[in,out] nsubla        counter of cell in the viscous sublayer
 * \param[in,out] nlogla        counter of cell in the log-layer
 * \param[out]    ustar         friction velocity
 * \param[out]    uk            friction velocity
 * \param[out]    yplus         dimensionless distance to the wall
 * \param[out]    ypup          yplus projected vel ratio
 * \param[out]    cofimp        \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                              turbulence production
 * \param[out]    dplus         dimensionless shift to the wall for scalable
 *                              wall functions
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_velocity(cs_wall_f_type_t  iwallf,
                           cs_real_t         l_visc,
                           cs_real_t         t_visc,
                           cs_real_t         vel,
                           cs_real_t         y,
                           cs_real_t         rough_d,
                           cs_real_t         rnnb,
                           cs_real_t         kinetic_en,
                           int              *iuntur,
                           cs_gnum_t        *nsubla,
                           cs_gnum_t        *nlogla,
                           cs_real_t        *ustar,
                           cs_real_t        *uk,
                           cs_real_t        *yplus,
                           cs_real_t        *ypup,
                           cs_real_t        *cofimp,
                           cs_real_t        *dplus)
{
  cs_real_t lmk;
  bool wf = true;

  /* Pseudo shift of the wall, 0 by default */
  *dplus = 0.;

  /* Activation of wall function by default */
  *iuntur = 1;

  /* Sand Grain roughness */
  cs_real_t sg_rough = rough_d * exp(cs_turb_xkappa*cs_turb_cstlog_rough);

  switch (iwallf) {
  case CS_WALL_F_DISABLED:
    cs_wall_functions_disabled(l_visc,
                               t_visc,
                               vel,
                               y,
                               iuntur,
                               nsubla,
                               nlogla,
                               ustar,
                               uk,
                               yplus,
                               dplus,
                               ypup,
                               cofimp);
    break;
  case CS_WALL_F_1SCALE_POWER:
    cs_wall_functions_1scale_power(l_visc,
                                   vel,
                                   y,
                                   iuntur,
                                   nsubla,
                                   nlogla,
                                   ustar,
                                   uk,
                                   yplus,
                                   ypup,
                                   cofimp);
    break;
  case CS_WALL_F_1SCALE_LOG:
    cs_wall_functions_1scale_log(l_visc,
                                 vel,
                                 y,
                                 iuntur,
                                 nsubla,
                                 nlogla,
                                 ustar,
                                 uk,
                                 yplus,
                                 ypup,
                                 cofimp);
    break;
  case CS_WALL_F_2SCALES_LOG:
    cs_wall_functions_2scales_log(l_visc,
                                  t_visc,
                                  vel,
                                  y,
                                  kinetic_en,
                                  iuntur,
                                  nsubla,
                                  nlogla,
                                  ustar,
                                  uk,
                                  yplus,
                                  ypup,
                                  cofimp);
    break;
  case CS_WALL_F_SCALABLE_2SCALES_LOG:
    cs_wall_functions_2scales_scalable(l_visc,
                                       t_visc,
                                       vel,
                                       y,
                                       kinetic_en,
                                       iuntur,
                                       nsubla,
                                       nlogla,
                                       ustar,
                                       uk,
                                       yplus,
                                       dplus,
                                       ypup,
                                       cofimp);
    break;
  case CS_WALL_F_2SCALES_VDRIEST:
   cs_wall_functions_2scales_vdriest(rnnb,
                                      l_visc,
                                      vel,
                                      y,
                                      kinetic_en,
                                      iuntur,
                                      nsubla,
                                      nlogla,
                                      ustar,
                                      uk,
                                      yplus,
                                      ypup,
                                      cofimp,
                                      &lmk,
                                      sg_rough,
                                      wf);
    break;
  case CS_WALL_F_2SCALES_SMOOTH_ROUGH:
    cs_wall_functions_2scales_smooth_rough(l_visc,
                                           t_visc,
                                           vel,
                                           y,
                                           rough_d,
                                           kinetic_en,
                                           iuntur,
                                           nsubla,
                                           nlogla,
                                           ustar,
                                           uk,
                                           yplus,
                                           dplus,
                                           ypup,
                                           cofimp);
    break;
  case CS_WALL_F_2SCALES_CONTINUOUS:
    cs_wall_functions_2scales_continuous(rnnb,
                                         l_visc,
                                         t_visc,
                                         vel,
                                         y,
                                         kinetic_en,
                                         iuntur,
                                         nsubla,
                                         nlogla,
                                         ustar,
                                         uk,
                                         yplus,
                                         ypup,
                                         cofimp);
    break;
  default:
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Compute the correction of the exchange coefficient between the
 *         fluid and the wall for a turbulent flow.
 *
 *  This is function of the dimensionless
 *  distance to the wall \f$ y^+ = \dfrac{\centip \centf u_\star}{\nu}\f$.
 *
 *  Then the return coefficient reads:
 *  \f[
 *  h_{tur} = Pr \dfrac{y^+}{T^+}
 *  \f]
 *
 * \param[in]     iwalfs        type of wall functions for scalar
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     prl           laminar Prandtl number
 * \param[in]     prt           turbulent Prandtl number
 * \param[in]     rough_t       scalar roughness lenghth scale
 * \param[in]     uk            velocity scale based on TKE
 * \param[in]     yplus         dimensionless distance to the wall
 * \param[in]     dplus         dimensionless distance for scalable
 *                              wall functions
 * \param[out]    htur          corrected exchange coefficient
 * \param[out]    yplim         value of the limit for \f$ y^+ \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_scalar(cs_wall_f_s_type_t  iwalfs,
                         cs_real_t           l_visc,
                         cs_real_t           prl,
                         cs_real_t           prt,
                         cs_real_t           rough_t,
                         cs_real_t           uk,
                         cs_real_t           yplus,
                         cs_real_t           dplus,
                         cs_real_t          *htur,
                         cs_real_t          *yplim)
{
  switch (iwalfs) {
  case CS_WALL_F_S_ARPACI_LARSEN:
    cs_wall_functions_s_arpaci_larsen(l_visc,
                                      prl,
                                      prt,
                                      rough_t,
                                      uk,
                                      yplus,
                                      dplus,
                                      htur,
                                      yplim);
    break;
  case CS_WALL_F_S_VDRIEST:
    cs_wall_functions_s_vdriest(prl,
                                prt,
                                yplus,
                                htur);
    break;
  case CS_WALL_F_S_SMOOTH_ROUGH:
    cs_wall_functions_s_smooth_rough(l_visc,
                                     prl,
                                     prt,
                                     rough_t,
                                     uk,
                                     yplus,
                                     dplus,
                                     htur);
    break;
  default:
    /* TODO Monin Obukhov or Louis atmospheric wall function
     * must be adapted to smooth wall functions.
     * Arpaci and Larsen wall functions are put as in previous versions of
     * code_Saturne.
     * */
    cs_wall_functions_s_arpaci_larsen(l_visc,
                                      prl,
                                      prt,
                                      rough_t,
                                      uk,
                                      yplus,
                                      dplus,
                                      htur,
                                      yplim);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute boundary contributions for all immersed boundaries.
 *
 * This is similar to clptur.
 *
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_immersed_boundary_wall_functions(int         f_id,
                                    cs_real_t  *st_exp,
                                    cs_real_t  *st_imp)
{
  const cs_domain_t *domain = cs_glob_domain;

  const cs_field_t  *f = cs_field_by_id(f_id);

  /* mesh quantities */
  cs_mesh_t *m = domain->mesh;
  cs_mesh_quantities_t *mq = domain->mesh_quantities;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_real_t  *cell_f_vol = mq->cell_f_vol;

  /*  Wall */
  const cs_real_t *c_w_face_surf
    = (const cs_real_t *restrict)mq->c_w_face_surf;
  const cs_real_3_t *c_w_face_normal
    = (const cs_real_3_t *restrict)mq->c_w_face_normal;
  const cs_real_t *c_w_dist_inv
    = (const cs_real_t *restrict)mq->c_w_dist_inv;

  cs_real_t *cpro_roughness = NULL;
  if (cs_field_by_name_try("immersed_boundary_roughness") != NULL)
    cpro_roughness = cs_field_by_name_try("immersed_boundary_roughness")->val;

  /* Dynamic viscosity */
  const cs_real_t  *cpro_mu = CS_F_(mu)->val;
  const cs_real_t  *cpro_mut = CS_F_(mu_t)->val;

  cs_field_t *ib_uk = cs_field_by_name_try("immersed_boundary_uk");
  cs_field_t *ib_yplus = cs_field_by_name_try("immersed_boundary_yplus");
  cs_field_t *ib_dplus = cs_field_by_name_try("immersed_boundary_dplus");

  /* Density */
  const cs_real_t *rho = (cs_real_t *)CS_F_(rho)->val;

  /* Velocity */
  const cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

  cs_wall_functions_t *wall_functions = cs_get_glob_wall_functions();
  const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();

  cs_real_6_t *cvar_rij = (CS_F_(rij) != NULL) ?
                          (cs_real_6_t *)CS_F_(rij)->val :
                          NULL;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  /* velocity */
  if (f == CS_F_(vel) && eqp->idiff > 0) {
    /* cast to 3D vectors for readability */
    cs_real_3_t   *_st_exp = (cs_real_3_t *)st_exp;
    cs_real_33_t  *_st_imp = (cs_real_33_t *)st_imp;

    cs_gnum_t nsubla = 0;
    cs_gnum_t nlogla = 0;
    cs_real_t rnnb = 0.;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      //Geometric quantities
      cs_real_t solid_surf = c_w_face_surf[c_id];

      if (solid_surf > cs_math_epzero*pow(cell_f_vol[c_id], 2./3.)) {

        cs_real_t wall_dist  = (c_w_dist_inv[c_id] < DBL_MIN) ?
                                0.:
                                1. / c_w_dist_inv[c_id];

        /* 1. Velocity components on the solid surface */

        // Unit normal
        cs_real_3_t nw;
        cs_math_3_normalize(c_w_face_normal[c_id], nw);

        // Tangential velocity
        cs_real_3_t nt;
        cs_real_t utau;
        cs_math_3_orthogonal_projection(nw, vel[c_id], nt);
        cs_real_t ek = 0.;
        utau = cs_math_3_norm(nt);

        // Unit tangent
        cs_math_3_normalize(nt, nt);

        /* 2. Friction velocities */

        if (utau <= cs_math_epzero) utau = cs_math_epzero;

        cs_real_t l_visc = cpro_mu[c_id]/rho[c_id];
        cs_real_t t_visc = cpro_mut[c_id]/rho[c_id];

        if (   turb_model->order == CS_TURB_FIRST_ORDER
            && turb_model->type  == CS_TURB_RANS) {
          ek = CS_F_(k)->val[c_id];
        }
        else if (   turb_model->order == CS_TURB_SECOND_ORDER
                 && turb_model->type  == CS_TURB_RANS) {

          ek = 0.5*(cvar_rij[c_id][0]+cvar_rij[c_id][1]+cvar_rij[c_id][2]);
          rnnb = cs_math_3_sym_33_3_dot_product(nw, cvar_rij[c_id], nw);
        }
        /* Cell wall roughness for immersed boundary */
        cs_real_t w_roughness = 0;
        if (cpro_roughness != NULL)
          w_roughness = cpro_roughness[c_id];

        int iuntur;
        cs_real_t ustar;
        cs_real_t uk;
        cs_real_t yplus;
        cs_real_t ypup;
        cs_real_t cofimp;
        cs_real_t dplus;

        /* If the cell is a solid cell, disable wall functions */

        cs_wall_f_type_t iwallf_loc = wall_functions->iwallf;
        if (mq->has_disable_flag) {
          if (mq->c_disable_flag[c_id])
            iwallf_loc = CS_WALL_F_DISABLED;
        }

        cs_wall_functions_velocity(iwallf_loc,
                                   l_visc,
                                   t_visc,
                                   utau,
                                   wall_dist,
                                   w_roughness,
                                   rnnb,
                                   ek,
                                   &iuntur,
                                   &nsubla,
                                   &nlogla,
                                   &ustar,
                                   &uk,
                                   &yplus,
                                   &ypup,
                                   &cofimp,
                                   &dplus);

        if (ib_uk != NULL)
          ib_uk->val[c_id] = uk;
        if (ib_yplus != NULL)
          ib_yplus->val[c_id] = yplus;
        if (ib_dplus != NULL)
          ib_dplus->val[c_id] = dplus;

        cs_real_t hint = (   turb_model->order == CS_TURB_SECOND_ORDER
                          && turb_model->type  == CS_TURB_RANS) ?
                         cpro_mu[c_id]*c_w_dist_inv[c_id] :
                         (cpro_mu[c_id]+cpro_mut[c_id])*c_w_dist_inv[c_id];

        cs_real_t hflui = cpro_mu[c_id]*c_w_dist_inv[c_id]*ypup;

        for (cs_lnum_t i = 0; i < 3; i++) {
          //TODO moving walls
          _st_exp[c_id][i] = 0.;
          for (cs_lnum_t j = 0; j < 3; j++) {
            _st_imp[c_id][i][j] = - ( hflui * cs_math_33_identity[i][j]
                                    + (hint - hflui) * nw[i] * nw[j] )
                                    * solid_surf;
          }
        }
      }
    }
    cs_parall_counter(&nsubla, 1);
    cs_parall_counter(&nlogla, 1);
  }
  /* Turbulent kinetic energy dissipation */
  if (f == CS_F_(eps)) {

    cs_real_t sigmae =
      cs_field_get_key_double(f, cs_field_key_id("turbulent_schmidt"));

    /* If turbulent Schmidt is variable, id of the corresponding field */
    int sigmae_id =
      cs_field_get_key_int(f, cs_field_key_id("turbulent_schmidt_id"));

    cs_real_t *cpro_sigmae = NULL;
    int sigmae_step = 0;
    if (sigmae_id >= 0) {
      cpro_sigmae = cs_field_by_id(sigmae_id)->val;
      sigmae_step = 1;
    }
    else {
      cpro_sigmae = &sigmae;
    }

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      //Geometric quantities
      cs_real_t solid_surf = c_w_face_surf[c_id];

      /* TODO loop over selection only */
      if (solid_surf > cs_math_epzero*pow(cell_f_vol[c_id], 2./3.) || true) {


        cs_real_t hint = ( cpro_mu[c_id]
                         + cpro_mut[c_id]/cpro_sigmae[c_id*sigmae_step])
                         * c_w_dist_inv[c_id];

        cs_real_t wall_dist  = (c_w_dist_inv[c_id] < DBL_MIN) ?
                                0.:
                                1. / c_w_dist_inv[c_id];

        /* Target imposed value at the wall */
        cs_real_t pimp = 0.;

        cs_real_t uk = ib_uk->val[c_id];
        cs_real_t yplus = ib_yplus->val[c_id];
        cs_real_t dplus = ib_dplus->val[c_id];

        /* If not in the viscous sublayer */
        const double ypluli = cs_glob_wall_functions->ypluli;
        if (yplus > ypluli) {
          cs_real_t l_visc = cpro_mu[c_id]/rho[c_id];
          pimp = 4. * wall_dist * cs_math_pow5(uk) /
              (cs_turb_xkappa * cs_math_pow2(l_visc * (yplus + 2.*dplus)));
        }

        st_exp[c_id] = hint * solid_surf * pimp;
        st_imp[c_id] = - hint * solid_surf;

      }
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
