/*============================================================================
 * Wall boundary condition management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_ale.h"
#include "cs_atmo.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_internal_coupling.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_rad_transfer.h"
#include "cs_thermal_model.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_set_coeffs_turb.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_lnum_t _ntlast = -1;

static cs_lnum_t _iaff = 0;

static int _kbfid = -1;

static const cs_lnum_t _iv2t[6] = {0, 1, 2, 0, 1, 0};
static const cs_lnum_t _jv2t[6] = {0, 1, 2, 1, 2, 2};

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_mo_compute_from_thermal_diff(const cs_real_t  *z,
                                  const cs_real_t  *z0,
                                  const cs_real_t  *du,
                                  const cs_real_t  *dt,
                                  const cs_real_t  *tm,
                                  const cs_real_t  *gredu,
                                  cs_real_t        *dlmo,
                                  cs_real_t        *ustar);

void
cs_f_mo_compute_from_thermal_flux(const cs_real_t  *z,
                                  const cs_real_t  *z0,
                                  const cs_real_t  *du,
                                  const cs_real_t  *flux,
                                  const cs_real_t  *tm,
                                  const cs_real_t  *gredu,
                                  cs_real_t        *dlmo,
                                  cs_real_t        *ustar);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute boundary coefficients for smooth/rough walls for scalar.
 *
 * \param[in]     f_sc          scalar field
 * \param[in]     isvhb         indicator to save exchange coeffient
 * \param[in]     byplus        dimensionless distance to the wall
 * \param[in]     bdplus        dimensionless shift to the wall
 *                              for scalable wall functions
 * \param[in]     buk           dimensionless velocity
 * \param[in]     buet          boundary ustar value
 * \param[in]     bcfnns        boundary correction factor
 * \param[in]     bdlmo         boundary Monin Obukhov length inverse,
 *                              or NULL
 * \param[in,out] hbord         exchange coefficient at boundary
 * \param[in]     theipb        value of thermal scalar at \f$ \centip \f$
 *                              of boundary cells
 * \param[out]    tetmax        maximum local ustar value
 * \param[out]    tetmin        minimum local ustar value
 * \param[out]    tplumx        maximum local tplus value
 * \param[out]    tplumn        minimum local tplus value
 */
/*----------------------------------------------------------------------------*/

static void
_cs_boundary_conditions_set_coeffs_turb_scalar(cs_field_t  *f_sc,
                                               int          isvhb,
                                               cs_real_t    byplus[],
                                               cs_real_t    bdplus[],
                                               cs_real_t    buk[],
                                               cs_real_t    buet[],
                                               cs_real_t    bcfnns[],
                                               cs_real_t    bdlmo[],
                                               cs_real_t    hbord[],
                                               cs_real_t    theipb[],
                                               cs_real_t   *tetmax,
                                               cs_real_t   *tetmin,
                                               cs_real_t   *tplumx,
                                               cs_real_t   *tplumn)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const cs_turb_model_type_t iturb = cs_glob_turb_model->iturb;
  const cs_real_t xkappa = cs_turb_xkappa;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)fvq->b_face_u_normal;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)fvq->b_face_cog;

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  const int kivisl  = cs_field_key_id("diffusivity_id");
  const int kturt   = cs_field_key_id("turbulent_flux_model");
  const int kscacp  = cs_field_key_id("is_temperature");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");

  const int ifcvsl = cs_field_get_key_int(f_sc, kivisl);
  const int thermal_variable = cs_glob_thermal_model->thermal_variable;
  const cs_field_t *f_th = cs_thermal_model_field();

  const cs_real_t cp0 = fluid_props->cp0;
  const cs_real_t cv0 = fluid_props->cv0;
  const int icp = fluid_props->icp;
  const int icv = fluid_props->icv;
  const cs_real_t rair = fluid_props->r_pg_cnst;

  const cs_wall_f_s_type_t iwalfs = cs_glob_wall_functions->iwalfs;

  const cs_real_t *viscls = NULL;
  if (ifcvsl >= 0)
    viscls = cs_field_by_id(ifcvsl)->val;

  cs_real_t *val_s = f_sc->val;
  cs_equation_param_t *eqp_sc = cs_field_get_equation_param(f_sc);

  /* If we have no diffusion, no boundary face should have a wall BC type
     (this is ensured in cs_boundary_conditions_type) */

  if (eqp_sc->idiff == 0) {
    *tetmax = 0.;
    *tetmin = 0.;
    *tplumx = 0.;
    *tplumn = 0.;
    return;
  }

  /* Get the turbulent flux model for the scalar */

  const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
  const cs_real_t ctheta = cs_field_get_key_double(f_sc, kctheta);

  const int turb_flux_model = cs_field_get_key_int(f_sc, kturt);
  const int turb_flux_model_type = turb_flux_model / 10;

  cs_real_6_t *visten = NULL;

  if (   eqp_sc->idften & CS_ANISOTROPIC_DIFFUSION
      || turb_flux_model_type == CS_TURB_HYBRID) {

    if (   iturb != CS_TURB_RIJ_EPSILON_EBRSM
        || turb_flux_model_type == CS_TURB_HYBRID) {
      cs_field_t *f_a_t_visc
        = cs_field_by_name("anisotropic_turbulent_viscosity");
      visten = (cs_real_6_t *)f_a_t_visc->val;
    }
    else { /* EBRSM and (GGDH or AFM) */
      cs_field_t *f_vis
        = cs_field_by_name("anisotropic_turbulent_viscosity_scalar");
      visten = (cs_real_6_t *)f_vis->val;
    }

  }

  cs_real_t *coefa_sc = f_sc->bc_coeffs->a;
  cs_real_t *coefb_sc = f_sc->bc_coeffs->b;
  cs_real_t *cofaf_sc = f_sc->bc_coeffs->af;
  cs_real_t *cofbf_sc = f_sc->bc_coeffs->bf;

  cs_real_t *crom = CS_F_(rho)->val;

  const cs_real_t *cpro_cp = NULL, *cpro_cv = NULL;
  if (icp >= 0)
    cpro_cp = CS_F_(cp)->val;

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0 && icv >= 0)
    cpro_cv = cs_field_by_id(icv)->val;

  const int keysca = cs_field_key_id("scalar_id");
  int scal_id = cs_field_get_key_int(f_sc, keysca);

  int isvhbl = 0;
  if (scal_id == isvhb)
    isvhbl = isvhb;

  if (f_sc == f_th) {
    /* min. and max. of wall friction of the thermal scalar */
    *tetmax = -cs_math_big_r;
    *tetmin =  cs_math_big_r;
    /* min. and max. of T+ */
    *tplumx = -cs_math_big_r;
    *tplumn =  cs_math_big_r;
  }

  cs_real_t rinfiv[3] = {cs_math_infinite_r,
                         cs_math_infinite_r,
                         cs_math_infinite_r};

  /* pointers to T+ and T* if saved */

  cs_real_t *tplusp = NULL, *tstarp = NULL;
  cs_real_t *dist_theipb = NULL;

  if (f_sc == f_th) {
    cs_field_t *itplus = cs_field_by_name_try("tplus");
    cs_field_t *itstar = cs_field_by_name_try("tstar");

    if (itplus != NULL)
      tplusp = itplus->val;

    if (itstar != NULL)
      tstarp = itstar->val;

    if (eqp_sc->icoupl > 0) {
      BFT_MALLOC(dist_theipb, n_b_faces, cs_real_t);
      cs_ic_field_dist_data_by_face_id(f_sc->id, 1, theipb, dist_theipb);
    }
  }

  cs_real_t *bpro_rough_t = NULL;
  cs_field_t *f_rough = cs_field_by_name_try("boundary_roughness");
  cs_field_t *f_rough_t = cs_field_by_name_try("boundary_thermal_roughness");

  if (f_rough != NULL)
    /* same thermal roughness if not specified */
    bpro_rough_t = f_rough->val;

  if (f_rough_t != NULL)
    bpro_rough_t = f_rough_t->val;

  bool *cpl_faces = NULL;
  if (eqp_sc->icoupl > 0) {
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    const int coupling_id = cs_field_get_key_int(f_sc, coupling_key_id);
    const cs_internal_coupling_t *cpl = cs_internal_coupling_by_id(coupling_id);

    cpl_faces = cpl->coupled_faces;
  }

  /* Pointers to specific fields */

  cs_real_t *bfconv  = NULL, *bhconv = NULL;

  if (cs_glob_rad_transfer_params->type >= 1 && f_sc == f_th) {
    bfconv = cs_field_by_name("rad_convective_flux")->val;
    bhconv = cs_field_by_name("rad_exchange_coefficient")->val;
  }

  /* FIXME not really the BC value */
  if (_kbfid < 0)
    _kbfid = cs_field_key_id("boundary_value_id");

  cs_field_t *f_scal_b = NULL;
  cs_real_t *bvar_s = NULL;

  int b_f_id = cs_field_get_key_int(f_sc, _kbfid);

  if (b_f_id > -1)
    f_scal_b = cs_field_by_id(b_f_id);
  else {
    /* if thermal variable has no boundary but temperature does, use it */
    if (f_sc == f_th && f_sc == CS_F_(h))
      f_scal_b = cs_field_by_name_try("boundary_temperature");
  }

  if (f_scal_b != NULL)
    bvar_s = f_scal_b->val;

  /* Does the scalar behave as a temperature ? */
  const int iscacp = cs_field_get_key_int(f_sc, kscacp);

  /* Retrieve turbulent Schmidt value for current scalar */
  const cs_real_t turb_schmidt = cs_field_get_key_double(f_sc, ksigmas);

  /* Reference diffusivity */
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  cs_real_t visls_0 = cs_field_get_key_double(f_sc, kvisl0);

  cs_field_t *f_id_cv = cs_field_by_name_try("isobaric_heat_capacity");
  if (f_id_cv != NULL)
    cpro_cv = f_id_cv->val;

  int *icodcl_vel = CS_F_(vel)->bc_coeffs->icodcl;
  int *icodcl_sc = f_sc->bc_coeffs->icodcl;
  cs_real_t *rcodcl1_sc = f_sc->bc_coeffs->rcodcl1;
  cs_real_t *rcodcl2_sc = f_sc->bc_coeffs->rcodcl2;
  cs_real_t *rcodcl3_sc = f_sc->bc_coeffs->rcodcl3;

  cs_real_t ypth = 0.0;

  cs_real_t *hbnd, *hint, *yptp;
  BFT_MALLOC(hbnd, n_b_faces, cs_real_t);
  BFT_MALLOC(hint, n_b_faces, cs_real_t);
  BFT_MALLOC(yptp, n_b_faces, cs_real_t);

  /* Loop on boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test on the presence of a smooth/rough wall condition (start) */
    if (icodcl_vel[f_id] != 5 && icodcl_vel[f_id] != 6)
      continue;

    const cs_lnum_t c_id = b_face_cells[f_id];

    /* Physical quantities */
    const cs_real_t yplus = byplus[f_id];
    const cs_real_t dplus = bdplus[f_id];
    const cs_real_t uk = buk[f_id];

    const cs_real_t visclc = viscl[c_id];
    const cs_real_t visctc = visct[c_id];
    const cs_real_t romc = crom[c_id];
    const cs_real_t xnuii = visclc / romc;

    /* Geometric quantities */
    const cs_real_t *rnxyz = b_face_u_normal[f_id];
    const cs_real_t distbf = b_dist[f_id];

    cs_real_t cpp = 1.;
    if (iscacp == 1)
      cpp = (icp >= 0) ? cpro_cp[c_id] : cp0;
    else if (iscacp == 2)
      cpp = (icp >= 0) ? cpro_cv[c_id] : cp0;

    const cs_real_t rkl = (ifcvsl < 0) ? visls_0 : viscls[c_id];
    cs_real_t prdtl = cpp * visclc / rkl;

    /* Compressible module:
       We assume that the Prandlt number should be defined in the same manner
       whether we solve for enthalpy or energy, that is Mu*Cp/Lambda.
       If we solve in energy we have computed Mu*Cv/Lambda above. */

    if (f_sc == f_th && thermal_variable == CS_THERMAL_MODEL_TOTAL_ENERGY) {
      prdtl = (icp >= 0) ? prdtl * cpro_cp[c_id] : prdtl * cp0;
      prdtl = (icv >= 0) ? prdtl / cpro_cv[c_id] : prdtl / cv0;
    }

    /* Scalar diffusivity */
    if (eqp_sc->idften & CS_ISOTROPIC_DIFFUSION) {

      /* In compressible, for energy Lambda/Cv+Cp/Cv*(mu_t/turb_schmidt) */
      if (f_sc == f_th && thermal_variable == CS_THERMAL_MODEL_TOTAL_ENERGY) {
        cs_real_t cpscv = (icp >= 0) ? cpro_cp[c_id] : cp0;
        cpscv = (icv >= 0) ? cpscv / cpro_cv[c_id] : cpscv / cv0;
        hint[f_id]
          = (rkl + eqp_sc->idifft * cpscv * visctc / turb_schmidt) / distbf;
      }
      else
        hint[f_id]
          = (rkl + eqp_sc->idifft * cpp * visctc / turb_schmidt) / distbf;
    }

    /* Symmetric tensor diffusivity (GGDH or AFM) */
    else if (eqp_sc->idften & CS_ANISOTROPIC_DIFFUSION) {

      cs_real_t temp = 0.0;
      cs_real_t visci[3][3], dist[3];

      dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
      dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
      dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

      /* In compressible, for energy Lambda/Cv+Cp/Cv*(mu_t/sigmas) */
      if (f_sc == f_th && thermal_variable == CS_THERMAL_MODEL_TOTAL_ENERGY) {
        cs_real_t cpscv = (icp >= 0) ? cpro_cp[c_id] : cp0;
        cpscv = (icv >= 0) ? cpscv / cpro_cv[c_id] : cpscv / cv0;
        temp = eqp_sc->idifft * cpscv * ctheta / cs_turb_csrij;
      }
      else
        temp = eqp_sc->idifft * cpp * ctheta / cs_turb_csrij;

      visci[0][0] = temp * visten[c_id][0] + rkl;
      visci[1][1] = temp * visten[c_id][1] + rkl;
      visci[2][2] = temp * visten[c_id][2] + rkl;
      visci[0][1] = temp * visten[c_id][3];
      visci[1][0] = temp * visten[c_id][3];
      visci[1][2] = temp * visten[c_id][4];
      visci[2][1] = temp * visten[c_id][4];
      visci[0][2] = temp * visten[c_id][5];
      visci[2][0] = temp * visten[c_id][5];

      /* ||Ki.n||^2 */
      const cs_real_t viscis =   cs_math_pow2(  visci[0][0]*rnxyz[0]
                                              + visci[1][0]*rnxyz[1]
                                              + visci[2][0]*rnxyz[2])
                               + cs_math_pow2(  visci[0][1]*rnxyz[0]
                                              + visci[1][1]*rnxyz[1]
                                              + visci[2][1]*rnxyz[2])
                               + cs_math_pow2(  visci[0][2]*rnxyz[0]
                                              + visci[1][2]*rnxyz[1]
                                              + visci[2][2]*rnxyz[2]);

      /* IF.Ki.n */
      cs_real_t fikis
        = (  cs_math_3_dot_product(dist, visci[0]) * rnxyz[0]
           + cs_math_3_dot_product(dist, visci[1]) * rnxyz[1]
           + cs_math_3_dot_product(dist, visci[2]) * rnxyz[2]);

      /* Take I so that I"F= eps*||FI||*Ki.n when I" is not in cell i
         NB: eps =1.d-1 must be consistent with vitens.f90 */
      fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distbf);

      hint[f_id] = viscis / fikis;
    }

    if (icodcl_vel[f_id] == 6)
      continue;

    cs_real_t hflui = 0.0;

    /* Wall function and Dirichlet or Neumann on the scalar */
    if (   iturb != CS_TURB_NONE
        && (   icodcl_sc[f_id] == 5
            || icodcl_sc[f_id] == 6
            || icodcl_sc[f_id] == 15
            || icodcl_sc[f_id] == 3)) {
      /* Note: to make things clearer yplus is always
         "y uk /nu" even for rough modelling. And the roughness correction is
         multiplied afterwards where needed. */
      const cs_real_t rough_t = (f_rough != NULL) ? bpro_rough_t[f_id] : 0;

      cs_wall_functions_scalar(cs_glob_wall_functions->iwalfs,
                               xnuii,
                               prdtl,
                               turb_schmidt,
                               rough_t,
                               uk,
                               yplus,
                               dplus,
                               &hflui,
                               &ypth);

      /* Correction for non-neutral condition in atmospheric flows */
      hflui *= bcfnns[f_id];

      /* Compute yk/T+ *PrT, take stability into account */
      yptp[f_id] = hflui / prdtl;

      /* Compute
         lambda/y * Pr_l *yk/T+ = lambda / nu * Pr_l * uk / T+ = rho cp uk / T+
         so "Pr_l * yk/T+" is the correction factor compared to a
         laminar profile */
      hflui *= rkl / distbf;

      /* User exchange coefficient */
      if (icodcl_sc[f_id] == 15) {
        hflui = rcodcl2_sc[f_id];
        yptp[f_id] = hflui / prdtl * distbf / rkl;
      }

    }
    else {
      /* y+/T+ *PrT */
      yptp[f_id] = 1.0 / prdtl;
      hflui = hint[f_id];
    }

    hbnd[f_id] = hflui; // = exchange_coeff, to save un new bc_coeffs structure.

  } /* End loop en boundary faces */

  /* internal coupling */
  if (eqp_sc->icoupl > 0) {
    /* Update exchange coef. in coupling entity of current scalar */
    cs_ic_field_set_exchcoeff(f_sc, hbnd);
  }

  /* Model-dependent fields */
  cs_field_t *f_tf = cs_field_by_composite_name_try(f_sc->name,
                                                    "turbulent_flux");
  cs_field_t *f_al = cs_field_by_composite_name_try(f_sc->name, "alpha");

  /* Loop on boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    if (icodcl_vel[f_id] != 5 && icodcl_vel[f_id] != 6)
      continue;

    const cs_real_t yplus = byplus[f_id];
    const cs_real_t dplus = bdplus[f_id];
    const cs_real_t uk = buk[f_id];

    /* Geometric quantities */
    const cs_lnum_t c_id = b_face_cells[f_id];
    const cs_real_t distbf = b_dist[f_id];
    const cs_real_t *rnxyz = b_face_u_normal[f_id];

    /* Physical quantities */
    const cs_real_t visclc = viscl[c_id];
    const cs_real_t visctc = visct[c_id];
    const cs_real_t romc = crom[c_id];

    cs_real_t cpp = 1.0;
    if (iscacp == 1)
      cpp = (icp >= 0) ? cpro_cp[c_id] : cp0;
    else if (iscacp == 2) {
      cpp = (icp >= 0) ? cpro_cp[c_id] - rair : cp0 - rair;
      /* FIXME: this formula does not seem consistent with that in diffst,
         but was present in clptrg.f90 */
      if (icodcl_vel[f_id] == 6) {
        cpp = (icp >= 0) ? cpro_cv[c_id] : cp0;
      }
    }

    const cs_real_t rkl = (ifcvsl < 0) ? visls_0 : viscls[c_id];

    const cs_real_t pimp = rcodcl1_sc[f_id];
    const cs_real_t hext = rcodcl2_sc[f_id];
    cs_real_t heq = 0.0, cofimp = 0.0, hflui = 0.0, tplus = 0.0;

    if (icodcl_vel[f_id] == 5) {
      hflui = hbnd[f_id];

      /* T+ = (T_I - T_w) / Tet */
      tplus = cs_math_fmax(yplus, cs_math_epzero) / yptp[f_id];
    }

    else if (icodcl_vel[f_id] == 6) {  /* rough wall (legacy) */

      /* Note: for Neumann, Tplus is chosen for post-processing */
      const cs_real_t rough_t = bpro_rough_t[f_id];

      /* Modified wall function from Louis */
      if (iwalfs != CS_WALL_F_S_MONIN_OBUKHOV) {

        /* T+ = (T_I - T_w) / Tet */
        tplus = log((distbf + rough_t) / rough_t) / (xkappa * bcfnns[f_id]);
      }
      else {
        /* Dry atmosphere, Monin Obukhov */
        const cs_real_t coef_moh
          = cs_mo_psih(distbf + rough_t, rough_t, bdlmo[f_id]);
        /* T+ */
        tplus = coef_moh / xkappa;
      }

      /* Dirichlet on the scalar, with wall function */
      if (iturb != CS_TURB_NONE && icodcl_sc[f_id] == 6) {
        /* 1/T+ */
        const cs_real_t dtplus = 1.0 / tplus;
        /* FIXME apparently buet should be buk */
        hflui = romc * cpp * buet[f_id] * dtplus;

        /* Neumann on the scalar, with wall function (for post-processing) */
      }
      else
        hflui = hint[f_id];

    } /* End hflui computation */

    /* Compute heq for smooth and rough wall */
    if (   cs_math_fabs(hext) > 0.5*cs_math_infinite_r
        || (icodcl_sc[f_id] == 15 && icodcl_vel[f_id] == 5)) {
      heq = hflui;
      if (eqp_sc->icoupl > 0 && icodcl_vel[f_id] == 5) {
        /* ensure correct saving of flux in case of rad coupling */
        if (cpl_faces[f_id])
          heq = hflui * hext / (hflui + hext);
      }
    }
    else
      heq = hflui * hext / (hflui + hext);

    /* Dirichlet Boundary condition with a wall function correction
       with or without an additional exchange coefficient hext */

    bool is_wall_scalar_std = (  icodcl_vel[f_id] == 5
                               && (   icodcl_sc[f_id] == 5
                                   || icodcl_sc[f_id] == 6
                                   || icodcl_sc[f_id] == 15));

    bool is_wall_scalar_rough_legacy = (    (icodcl_vel[f_id] == 6
                                        &&  icodcl_sc[f_id] == 6));

    if (is_wall_scalar_std || is_wall_scalar_rough_legacy) {

      if (is_wall_scalar_std) {

        /* DFM: the gradient BCs are so that the production term
           of u'T' is correcty computed */

        if (turb_flux_model_type >= 1) {
          /* In the log layer */
          if (yplus >= ypth && iturb != CS_TURB_NONE) {
            const cs_real_t xmutlm = xkappa * visclc * yplus;

            const cs_real_t mut_lm_dmut
              = (cs_mesh_quantities_cell_is_active(fvq, c_id) == 1) ?
                 xmutlm/cs_math_fmax(visctc,1.e-12*visclc) : 0.0;

            const cs_real_t rcprod
              = cs_math_fmin(xkappa,
                             cs_math_fmax(1.0, sqrt(mut_lm_dmut))
                             / (yplus + dplus));

            cofimp = 1.0 -   yptp[f_id] * turb_schmidt
              / xkappa * (2.0 * rcprod - 1.0 / (2.0 * yplus + dplus));
          }
          /* In the viscous sub-layer */
          else
            cofimp = 0.0;
        }
        else
          cofimp = 1.0 - heq / hint[f_id];

      }

      /* Rough wall (legacy) */
      else if (is_wall_scalar_rough_legacy) {

        /* FIXME this should also be done for Neumann, but overwritten in
           cs_boundary_condition_set_coeffs for now
           Same remark for smooth wall...
           if ((icodcl(ifac,ivar).eq.6).or.(icodcl(ifac,ivar).eq.3)) then */

        /* Modified wall function from Louis */
        if (iwalfs != CS_WALL_F_S_MONIN_OBUKHOV) {
          cofimp = 1.0 - heq / hint[f_id];
        }
        /* Monin obukhov */
        else {
          const cs_real_t rough_t = bpro_rough_t[f_id];

          /* To approximately respect thermal turbulent
             production with 2 hypothesis */
          //FIXME should be dynamic rhoughness
          const cs_real_t coef_mom = cs_mo_phim(distbf + rough_t, bdlmo[f_id]);
          const cs_real_t coef_mohh
            = cs_mo_phih (2.0 * distbf + rough_t, bdlmo[f_id]);

          const cs_real_t rcprod
            =   2.0 * romc / visctc * distbf * uk * tplus / coef_mom
              - coef_mohh / (2.0 + rough_t / distbf);

          cofimp = 1.0 - rcprod / (xkappa * tplus);
        }

      }

      /* To be coherent with a wall function, clip it to 0 */
      cofimp = cs_math_fmax(cofimp, 0.0);

      /* Gradient BCs */
      coefa_sc[f_id] = (1.0 - cofimp) * pimp;
      coefb_sc[f_id] = cofimp;

      /* Flux BCs */
      cofaf_sc[f_id] = - heq * pimp;
      cofbf_sc[f_id] =   heq;

      /* Set coef for coupled face just to ensure relevant saving
         of bfconv if rad transfer activated */
      if (eqp_sc->icoupl > 0 && icodcl_vel[f_id] == 5) {
        if (cpl_faces[f_id]) {
          /* Flux BCs */
          cofaf_sc[f_id] = - heq * dist_theipb[f_id];
          cofbf_sc[f_id] =   heq;
        }
      }

      /* Storage of the thermal exchange coefficient
         (conversion in case of energy or enthalpy)
         the exchange coefficient is in W/(m2 K)
         Useful for thermal coupling or radiative transfer */

      cs_real_t exchange_coef = 0.0;
      if (   (cs_glob_rad_transfer_params->type >= 1 && f_sc == f_th)
          || isvhbl > 0) {

        /* Enthalpy */
        if (thermal_variable == CS_THERMAL_MODEL_ENTHALPY) {
          /* If Cp is variable */
          exchange_coef = (icp >= 0) ? hflui*cpro_cp[c_id] : hflui * cp0;
        }
        /* Total energy (compressible module) */
        else if (thermal_variable == CS_THERMAL_MODEL_TOTAL_ENERGY) {
          /* If Cv is variable */
          exchange_coef = (icv >= 0) ? hflui*cpro_cv[c_id] : hflui * cv0;
        }
        /* Temperature */
        else if (iscacp > 0) {
          exchange_coef = hflui;
        }

      }

      /* Thermal coupling, store h = lambda/d */
      if (isvhbl > 0)
        hbord[f_id] = exchange_coef;

      /* Radiative transfer */
      if (cs_glob_rad_transfer_params->type >= 1 && f_sc == f_th) {
        bhconv[f_id] = exchange_coef;

        /* The outgoing flux is stored (Q = h(Ti'-Tp): negative if
           gain for the fluid) in W/m2 */
        bfconv[f_id] = cofaf_sc[f_id] + cofbf_sc[f_id] * theipb[f_id];
      }

      /* For the coupled faces with h_user (ie icodcl_sc[f_id]=15)
         reset to zero af/bf coeff.
         By default icodcl(f_id,ivar)=3) for coupled faces */
      if (eqp_sc->icoupl > 0 && icodcl_vel[f_id] == 5) {
        if (cpl_faces[f_id]) {
          /* Flux BCs */
          cofaf_sc[f_id] = 0.0;
          cofbf_sc[f_id] = 0.0;
        }
      }

    } /* End if icodcl == 5 or 6 or 15 */

    /* Turbulent heat flux */

    if (turb_flux_model_type == 3) {

      cs_real_3_t  *coefa_tf = (cs_real_3_t  *)f_tf->bc_coeffs->a;
      cs_real_33_t *coefb_tf = (cs_real_33_t *)f_tf->bc_coeffs->b;
      cs_real_3_t  *cofaf_tf = (cs_real_3_t  *)f_tf->bc_coeffs->af;
      cs_real_33_t *cofbf_tf = (cs_real_33_t *)f_tf->bc_coeffs->bf;
      cs_real_3_t  *cofar_tf = (cs_real_3_t  *)f_tf->bc_coeffs->ad;
      cs_real_33_t *cofbr_tf = (cs_real_33_t *)f_tf->bc_coeffs->bd;

      /* Turbulent diffusive flux of the scalar T
         (blending factor so that the component v'T' have only
         mu_T/(mu+mu_T)* Phi_T) */

      cs_real_t phit = 0.0, hintt[6];

      if (icodcl_vel[f_id] == 5) {
        if (   icodcl_sc[f_id] == 5
            || icodcl_sc[f_id] == 6
            || icodcl_sc[f_id] == 15)
          phit = cofaf_sc[f_id] + cofbf_sc[f_id] * val_s[c_id];

        else if (icodcl_sc[f_id] == 3)
          phit = rcodcl3_sc[f_id];

        else if (icodcl_sc[f_id] == 1)
          phit = heq * (val_s[c_id] - pimp);

        else
          phit = 0.0;
      }
      else if (icodcl_vel[f_id] == 6)
        phit = cofaf_sc[f_id] + cofbf_sc[f_id] * val_s[c_id];

      hintt[0] =   0.5*(visclc+rkl)/distbf
                 + visten[c_id][0]*ctheta/distbf/cs_turb_csrij;

      hintt[1] =   0.5*(visclc+rkl)/distbf
                 + visten[c_id][1]*ctheta/distbf/cs_turb_csrij;

      hintt[2] =   0.5*(visclc+rkl)/distbf
                 + visten[c_id][2]*ctheta/distbf/cs_turb_csrij;

      hintt[3] = visten[c_id][3]*ctheta/distbf/cs_turb_csrij;
      hintt[4] = visten[c_id][4]*ctheta/distbf/cs_turb_csrij;
      hintt[5] = visten[c_id][5]*ctheta/distbf/cs_turb_csrij;

      /* Dirichlet Boundary Condition
         ---------------------------- */

      /* Add rho*uk*Tet to T'v' in High Reynolds */

      cs_real_t pimpv[3];
      if (yplus >= ypth || icodcl_vel[f_id] == 6) {
        for (int i = 0; i < 3; i++)
          pimpv[i] = rnxyz[i] * phit / (cpp * romc);
      }
      else {
        for (int isou = 0; isou < 3; isou++)
          pimpv[isou] = 0.0;
      }

      /* Turbulent flux */
      cs_boundary_conditions_set_dirichlet_vector_aniso
        (coefa_tf[f_id], cofaf_tf[f_id],
         coefb_tf[f_id], cofbf_tf[f_id],
         pimpv, hintt, rinfiv);

      /* Boundary conditions used in the temperature equation */
      for (int isou = 0; isou < 3; isou++) {
        cofar_tf[f_id][isou] = 0.;
        for (int jsou = 0; jsou < 3; jsou++)
          cofbr_tf[f_id][isou][jsou] = 0.;
      }

    }

    /* EB-GGDH/AFM/DFM alpha boundary conditions */

    if (f_al != NULL && icodcl_vel[f_id] == 5) {

      cs_real_t *coefa_al = f_al->bc_coeffs->a;
      cs_real_t *coefb_al = f_al->bc_coeffs->b;
      cs_real_t *cofaf_al = f_al->bc_coeffs->af;
      cs_real_t *cofbf_al = f_al->bc_coeffs->bf;

      /* Dirichlet Boundary Condition
         ---------------------------- */

      const cs_real_t pimp_al = 0.;
      const cs_real_t hint_al = 1.0 / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_al[f_id],
                                                  &cofaf_al[f_id],
                                                  &coefb_al[f_id],
                                                  &cofbf_al[f_id],
                                                  pimp_al,
                                                  hint_al,
                                                  cs_math_infinite_r);
    }

    /* Save the values of T^star and T^+ for post-processing */

    if (b_f_id >= 0 || f_sc == f_th) {
      cs_real_t phit = 0.0;

      /* Wall function */
      if (   icodcl_vel[f_id] == 5
          && (   icodcl_sc[f_id] == 5
              || icodcl_sc[f_id] == 6
              || icodcl_sc[f_id] == 15)) {

        phit = (f_sc == f_th) ?
          cofaf_sc[f_id] + cofbf_sc[f_id] * theipb[f_id] :
          cofaf_sc[f_id] + cofbf_sc[f_id] * bvar_s[f_id];
      }
      else if (icodcl_sc[f_id] == 1 && icodcl_vel[f_id] == 5) {
        phit = (f_sc == f_th) ?
                heq * (theipb[f_id] - pimp):
                heq * (bvar_s[f_id] - pimp);
      }
      else if (icodcl_vel[f_id] == 6 && icodcl_sc[f_id] == 6) {
        phit = cofaf_sc[f_id] + cofbf_sc[f_id] * theipb[f_id];
      }
      /* Imposed flux with wall function for post-processing */
      else if (icodcl_sc[f_id] == 3)
        phit = rcodcl3_sc[f_id]; /* = 0 if current face is coupled */
      else
        phit = 0.0;

      /* if face is coupled */
      if (eqp_sc->icoupl > 0 && icodcl_vel[f_id] == 5) {
        if (cpl_faces[f_id])
          phit = heq * (theipb[f_id] - dist_theipb[f_id]);
      }

      const cs_real_t var = (icodcl_vel[f_id] == 5) ? uk : buet[f_id];
      /* FIXME Should be uk rather than ustar? */

      const cs_real_t tet
        = phit / (romc * cpp *cs_math_fmax(var, cs_math_epzero));

      if (b_f_id >= 0)
        bvar_s[f_id] -= tplus * tet;

      if (tplusp != NULL)
        tplusp[f_id] = tplus;

      if (tstarp != NULL)
        tstarp[f_id] = tet;

      if (f_sc == f_th) {
        *tetmax = cs_math_fmax(tet, *tetmax);
        *tetmin = cs_math_fmin(tet, *tetmin);
        *tplumx = cs_math_fmax(tplus, *tplumx);
        *tplumn = cs_math_fmin(tplus, *tplumn);
      }

    }

  } /* End loop on faces */

  BFT_FREE(hbnd);
  BFT_FREE(hint);
  BFT_FREE(yptp);
  BFT_FREE(dist_theipb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute boundary coefficients for smooth walls for vector.
 *
 * \param[in]     f_v           vector field
 * \param[in]     byplus        dimensionless distance to the wall
 * \param[in]     bdplus        dimensionless shift to the wall
 *                              for scalable wall functions
 * \param[in]     buk           dimensionless velocity
 */
/*----------------------------------------------------------------------------*/

static void
_cs_boundary_conditions_set_coeffs_turb_vector(cs_field_t  *f_v,
                                               cs_real_t    byplus[],
                                               cs_real_t    bdplus[],
                                               cs_real_t    buk[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const cs_turb_model_type_t iturb = cs_glob_turb_model->iturb;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)fvq->b_face_u_normal;

  const int kscacp  = cs_field_key_id("is_temperature");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const int kturt   = cs_field_key_id("turbulent_flux_model");
  const int kivisl  = cs_field_key_id("diffusivity_id");

  const cs_real_t cp0 = fluid_props->cp0;
  const int icp = fluid_props->icp;
  const cs_real_t rair = fluid_props->r_pg_cnst;

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  const int ifcvsl = cs_field_get_key_int(f_v, kivisl);

  const cs_real_t *viscls = NULL;
  if (ifcvsl >= 0)
    viscls = cs_field_by_id(ifcvsl)->val;

  cs_equation_param_t *eqp_v = cs_field_get_equation_param(f_v);
  cs_real_3_t  *coefa_v = (cs_real_3_t  *)f_v->bc_coeffs->a;
  cs_real_33_t *coefb_v = (cs_real_33_t *)f_v->bc_coeffs->b;
  cs_real_3_t  *cofaf_v = (cs_real_3_t  *)f_v->bc_coeffs->af;
  cs_real_33_t *cofbf_v = (cs_real_33_t *)f_v->bc_coeffs->bf;

  const cs_real_t *crom = crom = CS_F_(rho)->val;

  const cs_real_t *cpro_cp = NULL;
  if (icp >= 0)
    cpro_cp = CS_F_(cp)->val;

  /* Does the vector behave as a temperature ? */
  const int iscacp = cs_field_get_key_int(f_v, kscacp);

  /* Retrieve turbulent Schmidt value for current vector */
  const cs_real_t turb_schmidt = cs_field_get_key_double(f_v, ksigmas);

  /* Reference diffusivity */
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  cs_real_t visls_0 = cs_field_get_key_double(f_v, kvisl0);

  /* Get the turbulent flux model for the vector */
  const int turb_flux_model = cs_field_get_key_int(f_v, kturt);
  const int turb_flux_model_type = turb_flux_model / 10;

  cs_real_t *bpro_rough_t = NULL;

  cs_field_t *f_rough = cs_field_by_name_try("boundary_roughness");
  cs_field_t *f_rough_t = cs_field_by_name_try("boundary_thermal_roughness");

  if (f_rough != NULL)
    /* same thermal roughness if not specified */
    bpro_rough_t = f_rough->val;

  if (f_rough_t != NULL)
    bpro_rough_t = f_rough_t->val;

  cs_real_t *hbnd;
  BFT_MALLOC(hbnd, n_b_faces, cs_real_t);

  cs_real_t *hint;
  BFT_MALLOC(hint, n_b_faces, cs_real_t);

  const int *icodcl_vel = CS_F_(vel)->bc_coeffs->icodcl;
  const int *icodcl_v = f_v->bc_coeffs->icodcl;
  const cs_real_t *rcodcl1_v = f_v->bc_coeffs->rcodcl1;
  const cs_real_t *rcodcl2_v = f_v->bc_coeffs->rcodcl2;

  cs_real_t yptp = 0.0, ypth = 0.0;;

  /* Loop on boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test on the presence of a smooth wall condition (start) */
    if (icodcl_vel[f_id] != 5)
      continue;

    /* Geometric quantities */
    const cs_lnum_t c_id = b_face_cells[f_id];
    const cs_real_t distbf = b_dist[f_id];

    const cs_real_t yplus = byplus[f_id];
    const cs_real_t dplus = bdplus[f_id];
    const cs_real_t uk = buk[f_id];

    /* Physical quantities */
    const cs_real_t visclc = viscl[c_id];
    const cs_real_t visctc = visct[c_id];
    const cs_real_t romc = crom[c_id];
    const cs_real_t xnuii = visclc / romc;

    cs_real_t cpp = 1.;
    if (iscacp == 1)
      cpp = (icp >= 0) ? cpro_cp[c_id] : cp0;
    else if (iscacp == 2)
      cpp = (icp >= 0) ? cpro_cp[c_id] - rair : cp0 - rair;

    const cs_real_t rkl = (ifcvsl < 0) ? visls_0 : viscls[c_id];
    cs_real_t prdtl = cpp * visclc / rkl;

    /* Scalar diffusivity */
    if (eqp_v->idften & CS_ISOTROPIC_DIFFUSION)
      hint[f_id] = (rkl + eqp_v->idifft * cpp * visctc / turb_schmidt) / distbf;
    else {
      /* TODO if (vcopt%idften == 6) */
      bft_error(__FILE__, __LINE__, 0,
                "%s: case with anisotropic diffusion not handled.",
                __func__);
    }

    cs_real_t hflui = 0.0;

    /* Wall function and Dirichlet or Neumann on the scalar */
    if (iturb != CS_TURB_NONE && (icodcl_v[f_id] == 5 || icodcl_v[f_id] == 3)) {

      const cs_real_t rough_t = (f_rough != NULL) ? bpro_rough_t[f_id] : 0;

      /* FIXME use Re* = rough_t * uk / nu ? * PrT ? */
      cs_wall_functions_scalar(cs_glob_wall_functions->iwalfs,
                               xnuii,
                               prdtl,
                               turb_schmidt,
                               rough_t,
                               uk,
                               yplus,
                               dplus,
                               &hflui,
                               &ypth);

      /* Compute (y+-d+)/T+ *PrT */
      yptp = hflui / prdtl;

      /* Compute lambda/y * (y+-d+)/T+ */
      hflui = rkl/distbf *hflui;

    }
    /* User exchange coefficient */
    else if (icodcl_v[f_id] == 15)
      hflui = rcodcl2_v[f_id];

    else {
      /* y+/T+ *PrT */
      yptp = 1.0 / prdtl;
      hflui = hint[f_id];
    }

    hbnd[f_id] = hflui;

  } /* End loop on boundary faces */

  /* internal coupling */
  if (eqp_v->icoupl > 0)
    /* Update exchange coef. in coupling entity of current scalar */
    cs_ic_field_set_exchcoeff(f_v, hbnd);

  /* Loop on boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test on the presence of a smooth wall condition (start) */
    if (icodcl_vel[f_id] != 5)
      continue;

    const cs_real_t yplus = byplus[f_id];
    const cs_real_t dplus = bdplus[f_id];

    /* Geometric quantities */
    const cs_lnum_t c_id = b_face_cells[f_id];
    const cs_real_t *rnxyz = b_face_u_normal[f_id];

    /* Physical quantities */
    const cs_real_t visclc = viscl[c_id];
    const cs_real_t visctc = visct[c_id];
    const cs_real_t hext = rcodcl2_v[f_id];
    const cs_real_t hflui = hbnd[f_id];

    /* Local framework
       --------------- */

    /* Handle Dirichlet vector values */

    cs_real_t rcodcxyz[3] = {rcodcl1_v[n_b_faces*0 + f_id],
                             rcodcl1_v[n_b_faces*1 + f_id],
                             rcodcl1_v[n_b_faces*2 + f_id]};

    /* Keep tangential part */

    cs_real_t rcodcn = cs_math_3_dot_product(rcodcxyz, rnxyz);
    rcodcxyz[0] = rcodcxyz[0] - rcodcn * rnxyz[0];
    rcodcxyz[1] = rcodcxyz[1] - rcodcn * rnxyz[1];
    rcodcxyz[2] = rcodcxyz[2] - rcodcn * rnxyz[2];

    rcodcn = cs_math_3_dot_product(rcodcxyz, rnxyz);

    cs_real_t heq = 0.0;
    if (cs_math_fabs(hext) > 0.5*cs_math_infinite_r || icodcl_v[f_id] == 15)
      heq = hflui;
    else
      heq = hflui * hext / (hflui + hext);

    /* Dirichlet Boundary condition with a wall function correction
       with or without an additional exchange coefficient hext */

    cs_real_t cofimp = 0.0;
    if (icodcl_v[f_id] == 5 || icodcl_v[f_id] == 15) {
      /* DFM: the gradient BCs are so that the production term
         of u'T' is correcty computed */

      if (turb_flux_model_type >= 1) {
        /* In the log layer */
        if (yplus >= ypth && iturb != CS_TURB_NONE) {
          const cs_real_t xmutlm = cs_turb_xkappa * visclc * (yplus + dplus);
          const cs_real_t rcprod
            = cs_math_fmin(cs_turb_xkappa,
                             cs_math_fmax(1.0, sqrt(xmutlm/visctc))
                           / (yplus+dplus));

          cofimp =   1.0 - yptp * turb_schmidt / cs_turb_xkappa
                   * (2.0 * rcprod - 1.0 / (2.0 * yplus + dplus));
        }
        /* In the viscous sub-layer */
        else
          cofimp = 0.0;
      }
      else
        cofimp = 1.0 - heq / hint[f_id];

      /* To be coherent with a wall function, clip it to 0 */
      cofimp = cs_math_fmax(cofimp, 0.0);

      /* Coupled solving of the velocity components */

      /* Gradient boundary conditions
         ---------------------------- */

      coefa_v[f_id][0] =   (1.0 - cofimp) * (rcodcxyz[0] - rcodcn*rnxyz[0])
                         + rcodcn*rnxyz[0];

      coefa_v[f_id][1] =   (1.0 - cofimp) * (rcodcxyz[1] - rcodcn*rnxyz[1])
                         + rcodcn*rnxyz[1];

      coefa_v[f_id][2] =   (1.0 - cofimp) * (rcodcxyz[2] - rcodcn*rnxyz[2])
                         + rcodcn*rnxyz[2];

      /* Projection in order to have the vector parallel to the wall
         B = cofimp * ( IDENTITY - n x n ) */

      coefb_v[f_id][0][0] =   cofimp * (1.0 - rnxyz[0] * rnxyz[0]);
      coefb_v[f_id][1][1] =   cofimp * (1.0 - rnxyz[1] * rnxyz[1]);
      coefb_v[f_id][2][2] =   cofimp * (1.0 - rnxyz[2] * rnxyz[2]);
      coefb_v[f_id][0][1] = - cofimp * rnxyz[0] * rnxyz[1];
      coefb_v[f_id][0][2] = - cofimp * rnxyz[0] * rnxyz[2];
      coefb_v[f_id][1][2] = - cofimp * rnxyz[1] * rnxyz[2];
      coefb_v[f_id][1][0] =   coefb_v[f_id][0][1];
      coefb_v[f_id][2][1] =   coefb_v[f_id][1][2];
      coefb_v[f_id][2][0] =   coefb_v[f_id][0][2];

      /* Flux boundary conditions
         ------------------------ */

      cofaf_v[f_id][0] = - heq * (rcodcxyz[0] - rcodcn * rnxyz[0])
                         - hint[f_id] * rcodcn * rnxyz[0];

      cofaf_v[f_id][1] = - heq * (rcodcxyz[1] - rcodcn * rnxyz[1])
                         - hint[f_id] * rcodcn * rnxyz[1];

      cofaf_v[f_id][2] = - heq * (rcodcxyz[2] - rcodcn * rnxyz[2])
                         - hint[f_id] * rcodcn * rnxyz[2];

      /* Projection
         B = heq*( IDENTITY - n x n ) */

      cofbf_v[f_id][0][0] =   heq*(1.-rnxyz[0]*rnxyz[0])
                            + hint[f_id]*rnxyz[0]*rnxyz[0];

      cofbf_v[f_id][1][1] =   heq*(1.-rnxyz[1]*rnxyz[1])
                            + hint[f_id]*rnxyz[1]*rnxyz[1];

      cofbf_v[f_id][2][2] =   heq*(1.-rnxyz[2]*rnxyz[2])
                            + hint[f_id]*rnxyz[2]*rnxyz[2];

      cofbf_v[f_id][0][1] = (hint[f_id] - heq) * rnxyz[0] * rnxyz[1];
      cofbf_v[f_id][0][2] = (hint[f_id] - heq) * rnxyz[0] * rnxyz[2];
      cofbf_v[f_id][1][0] = (hint[f_id] - heq) * rnxyz[1] * rnxyz[0];
      cofbf_v[f_id][1][2] = (hint[f_id] - heq) * rnxyz[1] * rnxyz[2];
      cofbf_v[f_id][2][0] = (hint[f_id] - heq) * rnxyz[2] * rnxyz[0];
      cofbf_v[f_id][2][1] = (hint[f_id] - heq) * rnxyz[2] * rnxyz[1];

      /* TODO: postprocessing at the boundary */

    } /* End if icodcl 5 or 15 */

  } /* End loop on boudary faces */

  BFT_FREE(hbnd);
  BFT_FREE(hint);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute hflui and uiptn for smooth wall
 *
 * \param[in]        c_id         cell id
 * \param[in]        visclc       kinematic viscosity
 * \param[in]        visctc       turbulent kinematic viscosity
 * \param[in]        romc         density evaluated at cell_id
 * \param[in]        distbf       distance between the cell center and
                                  the center of gravity of the border face
 * \param[in]        utau         tangential mean velocity
 * \param[in]        uet          boundary ustar value
 * \param[in]        uk           dimensionless velocity
 * \param[in]        yplus        dimensionless distance to the wall
 * \param[in]        ypup         yplus projected vel ratio
 * \param[in]        dplus        dimensionless shift to the wall for scalable
 *                                wall functions
 * \param[inout]     cofimp       \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                                turbulence production
 * \param[inout]     rcprod       ?
 * \param[inout]     hflui        internal exchange coefficient
 * \param[inout]     uiptn        counter of reversal layer
 *
 */
/*----------------------------------------------------------------------------*/

static void
_update_physical_quantities_smooth_wall(const cs_lnum_t  c_id,
                                        const cs_real_t  visclc,
                                        const cs_real_t  visctc,
                                        const cs_real_t  romc,
                                        const cs_real_t  distbf,
                                        const cs_real_t  utau,
                                        const cs_real_t  uet,
                                        const cs_real_t  uk,
                                        const cs_real_t  yplus,
                                        const cs_real_t  ypup,
                                        const cs_real_t  dplus,
                                        cs_real_t       *cofimp,
                                        cs_real_t       *rcprod,
                                        cs_real_t       *hflui,
                                        cs_real_t       *uiptn)

{
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t xkappa = cs_turb_xkappa;
  const cs_turb_model_type_t iturb  = cs_glob_turb_model->iturb;
  const int itytur = cs_glob_turb_model->itytur;

  /* Deprecated power law (Werner & Wengle) */
  if (cs_glob_wall_functions->iwallf == 1) {
    *uiptn  =    utau + uet * cs_turb_apow * cs_turb_bpow
              * pow(yplus, cs_turb_bpow)*(pow(2, cs_turb_bpow - 1) - 2);
  }

  /* Dependent on the turbulence Model*/
  else {

    /* uiptn respects the production of k
       in a conditional manner --> rcprod coefficient */

    /* k-epsilon and k-omega
       --------------------- */

    if (itytur == 2 || iturb == CS_TURB_K_OMEGA) {

      const cs_real_t xmutlm = xkappa * visclc * (yplus + dplus);
      /* FIXME should be efvisc... */

      const cs_real_t mut_lm_dmut
        = (cs_mesh_quantities_cell_is_active(fvq, c_id) == 1) ?
        (xmutlm / visctc) : 0;

      /* If yplus=0, uiptn is set to 0 to avoid division by 0.
         By the way, in this case: iuntur=0 */

      if (yplus > cs_math_epzero) { /* TODO use iuntur == 1 */
        /*FIXME not valid for rough */
        *rcprod = cs_math_fmin(xkappa,
                                 cs_math_fmax(1.0,sqrt(mut_lm_dmut))
                               / (yplus+dplus));

        *uiptn =   utau - distbf * uet * uk * romc / xkappa / visclc
                 * (2.0*(*rcprod) - 1.0 / (2.0 * yplus + dplus));
      }
      else {
        *uiptn = 0.;
      }

    }

    /* No turbulence, mixing length or Rij-espilon
       -------------------------------------------*/

    else if (   iturb == CS_TURB_NONE || iturb == CS_TURB_MIXING_LENGTH
             || itytur == 3) {

      /* In the case of elliptic weighting, we should ignore the wall laws.
         So we use a test on the turbulence model:
         With LRR or SSG use wall laws, with EBRSM, use no-slip condition. */

      if (iturb == CS_TURB_RIJ_EPSILON_EBRSM || iturb ==  CS_TURB_NONE) {
        *uiptn = 0.;

      }
      else {

        /* If yplus=0, uiptn is set to 0 to avoid division by 0.
           By the way, in this case: iuntur = 0 */
        if (yplus > cs_math_epzero) /* FIXME use iuntur */
          *uiptn =   utau - distbf * romc * uet * uk / xkappa / visclc
                   * (2.0 / (yplus + dplus) - 1.0 / (2.0 * yplus + dplus));
        else
          *uiptn = 0.;

      }
    }

    /* LES and Spalart Allmaras
       ------------------------ */

    else if (itytur == 4 || iturb == CS_TURB_SPALART_ALLMARAS) {

      *uiptn  = utau - 1.5 * uet / xkappa;

      /* If (mu+mut) becomes zero (dynamic models),
         an arbitrary value is set
         (zero flux) but without any problems because the flux
         is really zero at this face. */
      if (visctc+visclc <= 0)
        *hflui = 0.; /* TODO check if this is useful as it is implicited
                       and overwritten below */

    }

    /* v2f
       --- */

    else if (itytur == 5) {

      /* With these conditions, no need to compute uiptmx, uiptmn
         and iuiptn which are 0 (initialization value) */
      *uiptn  = 0.;

    }
  }

  /* Implicit the term (rho*uet*uk) */
  *hflui = visclc / distbf * ypup;

  /* To be coherent with a wall function, clip it to 0 */
  *cofimp = cs_math_fmax(*cofimp, 0.);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update physical quantities for rough wall
 *
 * \param[in]         visclc       kinematic viscosity
 * \param[in]         visctc       turbulent kinematic viscosity
 * \param[in]         romc         density
 * \param[in]         distbf       distance between the cell center and
                                   the center of gravity of the border face
 * \param[in]         utau         tangential mean velocity
 * \param[in]         uet          friction velocity
 * \param[in]         uk           friction velocity
 * \param[in]         uplus        dimensionless velocity
 * \param[in]         dplus        dimensionless shift to the wall for scalable
 *                                 wall functions
 * \param[in]         rough_d      roughness length scale (not sand grain)
 * \param[in]         dlmo         inverse Monin Obukhov length (for log only)
 * \param[in,out]     iuntur       indicator: 0 in the viscous sublayer
 * \param[in,out]     nlogla       counter of cell in the log-layer
 * \param[in,out]     nsubla       counter of cell in the viscous sublayer
 * \param[in,out]     cofimp       \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                                 turbulence production
 * \param[in,out]     rcprod       ?
 * \param[in,out]     hflui        internal exchange coefficient
 * \param[in,out]     uiptn        counter of reversal layer
 *
 */
/*----------------------------------------------------------------------------*/

static void
_update_physical_quantities_rough_wall(const cs_real_t  visclc,
                                       const cs_real_t  visctc,
                                       const cs_real_t  romc,
                                       const cs_real_t  distbf,
                                       const cs_real_t  utau,
                                       const cs_real_t  uet,
                                       const cs_real_t  uk,
                                       const cs_real_t  uplus,
                                       const cs_real_t  rough_d,
                                       const cs_real_t  dlmo,
                                       int             *iuntur,
                                       cs_gnum_t       *nlogla,
                                       cs_gnum_t       *nsubla,
                                       cs_real_t       *cofimp,
                                       cs_real_t       *rcprod,
                                       cs_real_t       *hflui,
                                       cs_real_t       *uiptn)
{
  const cs_real_t xkappa = cs_turb_xkappa;
  const cs_turb_model_type_t iturb  = cs_glob_turb_model->iturb;
  const int itytur = cs_glob_turb_model->itytur;
  const cs_wall_f_s_type_t iwalfs = cs_glob_wall_functions->iwalfs;

  /* uiptn respecte la production de k
     de facon conditionnelle --> Coef RCPROD

     All turbulence models (except v2f and EBRSM)
    -------------------------------------------- */

  if (   iturb == CS_TURB_NONE || itytur == 2 || itytur == 4
      || iturb == CS_TURB_K_OMEGA
      || iturb == CS_TURB_MIXING_LENGTH
      || iturb == CS_TURB_RIJ_EPSILON_LRR
      || iturb == CS_TURB_RIJ_EPSILON_SSG
      || iturb == CS_TURB_SPALART_ALLMARAS) {

    if (visctc > cs_math_epzero) {

      /* Pseudo shift of wall by rough_d ((distbf+rough_d)/rough_d) */
      const cs_real_t distb0 = distbf + rough_d;

      /* FIXME uk not modified for Louis yet.... */
      const cs_real_t xmutlm = xkappa * uk * distb0 * romc;

      if (iwalfs != CS_WALL_F_S_MONIN_OBUKHOV) {

        const cs_real_t var
          =  2.0 * sqrt(xmutlm/visctc) - distb0/distbf/(2.0 + rough_d / distb0);

        *rcprod = distbf / distb0 * cs_math_fmax(1.0, var);

        /* Ground apparent velocity (for log only) */
        *uiptn  = cs_math_fmax(utau - uet/xkappa* (*rcprod), 0.0);
        *iuntur = 1;
        *nlogla = *nlogla + 1;

        /* Coupled solving of the velocity components
           The boundary term for velocity gradient is implicit
           modified for non neutral boundary layer (in uplus) */

        *cofimp  = cs_math_fmax(1.0 - 1.0/(xkappa*uplus)*(*rcprod), 0.0);

        /*The term (rho*uet*uk) is implicit */

        /* TODO merge with MO without this max */
        const cs_real_t rcflux = cs_math_fmax(xmutlm, visctc) / distb0;

        *hflui = rcflux / (xkappa * uplus);
      }
      /* Monin Obukhov */
      else {

        /* Boundary condition on the velocity to have approximately the good
           turbulence production */

        const cs_real_t coef_mom = cs_mo_phim(distbf+rough_d, dlmo);
        const cs_real_t coef_momm = cs_mo_phim(2.0*distbf+rough_d, dlmo);

        *rcprod =   2*distbf*sqrt( xkappa*uk*romc*coef_mom/visctc/distb0 )
                 - coef_momm / (2.0 + rough_d / distbf);

        /* Ground apparent velocity (for log only) */
        *uiptn  = cs_math_fmax(utau - uet/xkappa*(*rcprod), 0.0);
        *iuntur = 1;
        *nlogla = *nlogla + 1;

        /* Coupled solving of the velocity components
          The boundary term for velocity gradient is implicit
          modified for non neutral boundary layer (in uplus) */

        *cofimp  = cs_math_fmin(cs_math_fmax(1 - 1/(xkappa*uplus)*(*rcprod), 0),
                                1);

        /* The term (rho*uet*uk) is implicit */
        *hflui = romc * uk / uplus;

      }

    }
    /* In the viscous sub-layer */
    else {
      *uiptn  = 0.0;
      *iuntur = 0;
      *nsubla = *nsubla + 1;

      /* Coupled solving of the velocity components */
      *cofimp  = 0.0;
      *hflui = visclc / distbf;
    }

  }

  /* Clipping :
     We bound U_f, grad between 0 and Utau (we could probably do better...)
     - 0    : forbid inversion at boundary, which is in contradiction
              with the log law hypothesis.
     - Utau : the turbulent production cannot be zero.
              We prevent U_f, flux from being negative */

  /* v2f and EBRSM !FIXME EBRSM
     --------------------------*/

  else if (itytur == 5) {

    /* With these conditions, no need to compute uiptmx, uiptmn
       and iuiptn which are zero (initialization value) */
    *iuntur = 0;
    *uiptn = 0.0;

    /* Coupled solving of the velocity components */
    *hflui = (visclc + visctc) / distbf;
    *cofimp = 0.0;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute friction velocity u* and surface sensible heat flux q0
 * for a non neutral atmospheric surface layer using the explicit formula
 * developed for the ECMWF by Louis (1982)

 * \param[in]     ifac        treated boundary face
 * \param[in]     utau        tangential mean velocity
 * \param[in]     rough_d     roughness z0
 * \param[in]     duplus      1 over dimensionless velocity in neutral
 *                            conditions
 * \param[in]     dtplus      1 over dimensionless temperature in neutral
 *                            conditions
 * \param[in]     yplus_t     thermal dimensionless wall distance
 * \param[out]    uet         friction velocity
 * \param[out]    gredu       reduced gravity for non horizontal wall
 * \param[out]    cfnns       non neutral correction coefficients for profiles
                              of scalar
 * \param[out]    cfnnk       non neutral correction coefficients
                              for profiles of k
 * \param[out]    cfnne       non neutral correction coefficients
                              for profiles of eps
 * \param[out]    dlmo        inverse Monin Obukhov length (for log only)
 * \param[in]     temp        potential temperature in boundary cell
 * \param[in]     totwt       total water content in boundary cell
 * \param[in]     liqwt       liquid water content in boundary cell
 */
/*----------------------------------------------------------------------------*/

static void
_atmo_cls(const cs_lnum_t  f_id,
          const cs_real_t  utau,
          const cs_real_t  rough_d,
          const cs_real_t  duplus,
          const cs_real_t  dtplus,
          const cs_real_t  yplus_t,
          cs_real_t       *uet,
          const cs_real_t  gredu,
          cs_real_t       *cfnns,
          cs_real_t       *cfnnk,
          cs_real_t       *cfnne,
          cs_real_t       *dlmo,
          const cs_real_t  temp,
          const cs_real_t  totwt,
          const cs_real_t  liqwt)
{
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t nt_cur = cs_glob_time_step->nt_cur;

  cs_field_t *f_th = cs_thermal_model_field();
  const cs_real_t *rcodcl1_th = f_th->bc_coeffs->rcodcl1;
  const int *icodcl_th = f_th->bc_coeffs->icodcl;

  cs_field_t *ym_water = cs_field_by_name_try("ym_water");
  cs_real_t *rcodcl1_ymw = (ym_water != NULL) ?
    ym_water->bc_coeffs->rcodcl1 : NULL;

  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_t distbf = b_dist[f_id];

  const cs_real_t rvsra = cs_glob_fluid_properties->rvsra;

  /* Initializations
     --------------- */

  const cs_real_t b = 5.0;
  const cs_real_t c = b;
  const cs_real_t d = b;

  /* TODO: Take into account humidity in ratio r/cp */
  /*
  const cs_real_t cpvcpa = cs_glob_air_props->cp_v / cs_glob_air_props->cp_a;

  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const cs_real_t cp0 = fluid_props->cp0;
  const cs_real_t rair = fluid_props->r_pg_cnst;

  cs_real_t rscp1, rscp2;
  if (ym_water != NULL) {
    rscp1 = (rair / cp0) * (1.0 + (rvsra-cpvcpa) * rcodcl1_ymw[f_id]);
    // Bouzerau PhD
    rscp2 = (rair / cp0) * (1.0 + (rvsra-cpvcpa) * (totwt-liqwt));
  }
  else {
    rscp1 = rair / cp0;
    rscp2 = rair / cp0;
  }
  */

  const cs_real_t tpot1 = rcodcl1_th[f_id];
  const cs_real_t tpot2 = temp;

  /* Compute virtual potential temperature at two levels */
  cs_real_t tpotv1, tpotv2;
  if (ym_water != NULL) {
    tpotv1 = tpot1 * (1.0 + (rvsra - 1.0) * rcodcl1_ymw[f_id]);
      /* Bouzerau PhD */
    tpotv2 = tpot2 * (1.0 + (rvsra - 1.0) * (totwt-liqwt));
  }
  else {
    tpotv1 = tpot1;
    tpotv2 = tpot2;
  }

  /* Patch for the initial time step when thermal field is not initalized */
  if (nt_cur == 1)
    tpotv2 = tpotv1;

  /* Compute layer average Richardson number */

  /* NB: rib = 0 if thermal flux conditions are imposed and tpot1 not defined */
  cs_real_t rib;
  if (cs_math_fabs(utau) < cs_math_epzero || icodcl_th[f_id] == 3)
    rib = 0.0;
  else
    rib = 2*gredu*distbf*(tpotv2 - tpotv1)/(tpotv1 + tpotv2)/utau/utau;

  /* Compute correction factors based on ECMWF parametrisation
     Louis (1982) */

  cs_real_t fm, fh, fmden1, fmden2, fhden;
  if (rib >= cs_math_epzero) {
    /* Stable case */
    fm = 1.0 / (1.0 + 2.0*b*rib/sqrt(1.0 + d*rib));
    fh = 1.0 / (1.0 + 3.0*b*rib*sqrt(1.0 + d*rib));
  }
  else {
    /* Unstable case */
    fmden1 = (yplus_t + 1.0) * cs_math_fabs(rib);
    fmden2 = 1.0 + 3.0 * b * c * duplus *dtplus * sqrt(fmden1);
    fm = 1.0 - 2.0 * b * rib / fmden2;
    fhden = 3.0 * b * c * duplus * dtplus * sqrt(yplus_t + 1.0);
    fh = 1.0 - (3.0*b*rib)/(1.0 + fhden * sqrt(cs_math_fabs(rib)));
  }

  if (fm <= cs_math_epzero)
    fm = cs_math_epzero;

  if (cs_math_fabs(fh) <= cs_math_epzero)
    fh = cs_math_epzero;

  if ((1.0-rib) > cs_math_epzero) {
    *cfnnk = sqrt(1.0 - rib); /* +correction with turbulent Prandtl */
    *cfnne = (1.0 - rib) / sqrt(fm);
  }
  else {
    *cfnnk = 1.0;
    *cfnne = 1.0;
  }

  /* Note: non neutral correction coefficients for profiles of wind
     (Re) compute friction velocity uet (for non neutral)
     uet = U/U^+ = U / U^{+,n} * sqrt(fm) */
  *uet = duplus * utau * sqrt(fm);

  /* Compute surface sensible heat flux q0 (can be useful for post-processing)
     Note: non-consistent with two velocity scales */
  *cfnns = fh / sqrt(fm);
  //const cs_real_t q0 = (tpot1 - tpot2) * (*uet) * dtplus * (*cfnns);
  /* FIXME tet should be output as uet is... */

  /* Compute local Monin-Obukhov inverse length for log
     1/L =  Ri / (z Phim) */
  *dlmo = rib * sqrt(fm) / (distbf + rough_d);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Boundary conditions for smooth walls (icodcl = 5).
 *
 * The wall functions may change the value of the diffusive flux.

 * The values at a boundary face \f$ \fib \f$ stored in the face center
 * \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
 * are written as:
 * \f[
 * P_{\face} = A_P^g + B_P^g P_{\centi}
 * \f]
 * and
 * \f[
 * Q_{\face} = A_P^f + B_P^f P_{\centi}
 * \f]
 * where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
 * neighboring cell.

 * Warning:

 * - For a vector field such as the velocity \f$ \vect{u} \f$ the boundary
 *   conditions may read:
 *   \f[
 *   \vect{u}_{\face} = \vect{A}_u^g + \tens{B}_u^g \vect{u}_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \vect{Q}_{\face} = \vect{A}_u^f + \tens{B}_u^f \vect{u}_{\centi}
 *   \f]
 *   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
 *   which coupled velocity components next to a boundary.

 * Please refer to the
 * <a href="../../theory.pdf#wallboundary"><b>wall boundary conditions</b></a>
 * section of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#clptur"><b>clptur</b></a> section.

 * \param[in]     isvhb         indicator to save exchange coeffient
 * \param[in]     velipb        value of the velocity at \f$ \centip \f$
 *                              of boundary cells
 * \param[in]     rijipb        value of \f$ R_{ij} \f$ at \f$ \centip \f$
 *                              of boundary cells
 * \param[out]    visvdr        dynamic viscosity after V. Driest damping in
 *                              boundary cells
 * \param[out]    hbord         exchange coefficient at boundary
 * \param[in]     theipb        value of thermal scalar at \f$ \centip \f$
 *                              of boundary cells
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs_turb(int        isvhb,
                                       cs_real_t  velipb[][3],
                                       cs_real_t  rijipb[][6],
                                       cs_real_t  visvdr[],
                                       cs_real_t  hbord[],
                                       cs_real_t  theipb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)fvq->b_face_u_normal;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)fvq->b_face_cog;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;

  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
  cs_field_t *f_th = cs_thermal_model_field();

  const cs_real_t pref = cs_glob_atmo_constants->ps;
  const cs_real_t p0 = fluid_props->p0;
  const cs_real_t r_pg_cnst = fluid_props->r_pg_cnst;
  const cs_real_t cp0 = fluid_props->cp0;
  const cs_real_t tref = fluid_props->t0;
  const int icp = fluid_props->icp;

  const int keysca  = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const int kdflim = cs_field_key_id("diffusion_limiter_id");

  /* Type of wall functions for scalar */
  const cs_wall_f_s_type_t iwalfs = cs_glob_wall_functions->iwalfs;

  const cs_real_t xkappa = cs_turb_xkappa;

  const cs_turb_model_t *turb_model = cs_get_glob_turb_model();
  const cs_turb_model_type_t iturb  = cs_glob_turb_model->iturb;
  const int n_fields = cs_field_n_fields();

  const cs_lnum_t nt_cur = cs_glob_time_step->nt_cur;
  const cs_lnum_t nt_max = cs_glob_time_step->nt_max;

  /* Initializations
     =============== */

  /* Initialize variables to avoid compiler warnings */

  cs_real_t cofimp = 0;
  cs_real_t ek = 0.;
  cs_real_t rcprod = 0;
  cs_real_t uiptn = 0;
  cs_real_t rnnb = 0;

  cs_real_t uet = 1;
  cs_real_t utau = 1;

  /* Constants */
  const cs_real_t sqrcmu = sqrt(cs_turb_cmu);

  /* Correction factors for stratification (used in atmospheric models) */
  cs_real_t cfnns = 1;
  cs_real_t cfnnk = 1;
  cs_real_t cfnne = 1;

  cs_field_t *rough = cs_field_by_name_try("boundary_roughness");
  cs_field_t *rough_t
    = cs_field_by_name_try("boundary_thermal_roughness");

  cs_real_t *bpro_rough = NULL;
  cs_real_t *bpro_rough_t = NULL;

  if (rough != NULL) {
    bpro_rough = rough->val;
    bpro_rough_t = rough->val;
  }

  if (rough_t != NULL)
    bpro_rough_t = rough_t->val;

  cs_field_t *boundary_ustar = cs_field_by_name_try("boundary_ustar");
  cs_field_t *boundary_uk = cs_field_by_name_try("boundary_uk");

  cs_real_t *bpro_ustar = NULL, *buet = NULL;

  /* Save wall friction velocity */

  if (boundary_ustar != NULL) {
    bpro_ustar = boundary_ustar->val;
  }
  else {
    BFT_MALLOC(buet, n_b_faces, cs_real_t);
    bpro_ustar = buet;
  }

  cs_real_t *bpro_uk = NULL, *_buk = NULL;

  if (boundary_uk != NULL) {
    bpro_uk = boundary_uk->val;
  }
  else {
    BFT_MALLOC(_buk, n_b_faces, cs_real_t);
    bpro_uk = _buk;
  }

  /* Pointers to y+ if saved */
  cs_field_t *f_yplus = cs_field_by_name_try("yplus");
  cs_real_t *yplbr = NULL;
  if (f_yplus != NULL)
    yplbr = f_yplus->val;

  const int itytur = cs_glob_turb_model->itytur;
  const int idirsm = cs_glob_turb_rans_model->idirsm;

  cs_field_t *f_a_t_visc = NULL;
  cs_real_6_t *visten = NULL;

  if (itytur == 3 && idirsm == 1) {
    f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
    visten = (cs_real_6_t *)f_a_t_visc->val;
  }

  /* Diffusion limiter for rough wall */

  cs_real_t *df_limiter_eps = NULL;
  cs_real_t *df_limiter_k = NULL;
  cs_real_t *df_limiter_rij = NULL;

  /* Gradient and flux boundary conditions */
  cs_field_t *vel = CS_F_(vel);
  cs_real_3_t  *coefa_vel = (cs_real_3_t  *)vel->bc_coeffs->a;
  cs_real_33_t *coefb_vel = (cs_real_33_t *)vel->bc_coeffs->b;
  cs_real_3_t  *cofaf_vel = (cs_real_3_t  *)vel->bc_coeffs->af;
  cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

  /* Lagrangian time scale */
  cs_field_t *f_tlag = cs_field_by_name_try("lagr_time");

  /* Physical quantities */

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;
  cs_real_t *visct = CS_F_(mu_t)->val;

  const cs_real_t *cpro_cp = NULL;
  if (icp >= 0)
    cpro_cp = (const cs_real_t *)CS_F_(cp)->val;

  cs_field_t *f_k = NULL, *f_eps = NULL, *f_rij = NULL, *f_alpha = NULL;
  cs_field_t *f_phi = NULL, *f_f_bar = NULL, *f_omg = NULL, *f_nusa = NULL;
  cs_equation_param_t *eqp_rij = NULL, *eqp_eps = NULL, *eqp_nusa = NULL;

  /* Turbulence variables */

  if (itytur == 2 || itytur == 5) {
    f_eps = CS_F_(eps);
    f_k = CS_F_(k);
    if (iturb == CS_TURB_V2F_PHI) {
      f_phi = CS_F_(phi);
      f_f_bar = CS_F_(f_bar);
    }
    else if (iturb == CS_TURB_V2F_BL_V2K) {
      f_phi = CS_F_(phi);
      f_alpha = CS_F_(alp_bl);
    }
  }
  else if (itytur == 3) {
    f_eps = CS_F_(eps);
    f_rij = CS_F_(rij);
    if (iturb == CS_TURB_RIJ_EPSILON_EBRSM)
      f_alpha = CS_F_(alp_bl);
    eqp_eps = cs_field_get_equation_param(f_eps);
    eqp_rij = cs_field_get_equation_param(f_rij);
  }
  else if (iturb == CS_TURB_K_OMEGA) {
    f_k = CS_F_(k);
    f_omg = CS_F_(omg);
  }
  else if (iturb == CS_TURB_SPALART_ALLMARAS) {
    f_nusa = CS_F_(nusa);
    eqp_nusa = cs_field_get_equation_param(f_nusa);
  }

  const cs_real_t sigmak = (f_k != NULL) ?
    cs_field_get_key_double(f_k, ksigmas) : 0;
  const cs_real_t sigmae = (f_eps != NULL) ?
    cs_field_get_key_double(f_eps, ksigmas) : 0;

  if (f_eps != NULL) {
    if (    (f_eps->type & CS_FIELD_VARIABLE)
        && !(f_eps->type & CS_FIELD_CDO)) {
      int df_limiter_id = cs_field_get_key_int(f_eps, kdflim);
      if (df_limiter_id > -1)
        df_limiter_k = cs_field_by_id(df_limiter_id)->val;
    }
  }

  cs_real_t *cvar_k = (f_k != NULL) ? f_k->val : NULL;
  if (f_k != NULL) {
    if (    (f_k->type & CS_FIELD_VARIABLE)
        && !(f_k->type & CS_FIELD_CDO)) {
      int df_limiter_id = cs_field_get_key_int(f_k, kdflim);
      if (df_limiter_id > -1)
        df_limiter_k = cs_field_by_id(df_limiter_id)->val;
    }
  }

  cs_real_6_t *cvar_rij = (f_rij != NULL) ? (cs_real_6_t *)f_rij->val : NULL;

  if (f_rij != NULL) {
    if (    (f_rij->type & CS_FIELD_VARIABLE)
        && !(f_rij->type & CS_FIELD_CDO)) {
      int df_limiter_id = cs_field_get_key_int(f_rij, kdflim);
      if (df_limiter_id > -1)
        df_limiter_rij = cs_field_by_id(df_limiter_id)->val;
    }
  }

  /* min. and max. of wall tangential velocity */
  cs_real_t uiptmx = -cs_math_big_r;
  cs_real_t uiptmn = cs_math_big_r;

  /* min. and max. of wall friction velocity */
  cs_real_t uetmax = -cs_math_big_r;
  cs_real_t uetmin =  cs_math_big_r;
  cs_real_t ukmax  = -cs_math_big_r;
  cs_real_t ukmin  =  cs_math_big_r;

  /* min. and max. of y+ */
  cs_real_t yplumx = -cs_math_big_r;
  cs_real_t yplumn =  cs_math_big_r;

  /* min. and max. of wall friction of the thermal scalar */
  cs_real_t tetmax = -cs_math_big_r;
  cs_real_t tetmin =  cs_math_big_r;

  /* min. and max. of inverse of MO length */
  cs_real_t dlmomax = -cs_math_big_r;
  cs_real_t dlmomin =  cs_math_big_r;

  /* min. and max. of T+ */
  cs_real_t tplumx = -cs_math_big_r;
  cs_real_t tplumn =  cs_math_big_r;

 /* Counters (turbulent, laminar, reversal, scale correction) */
  cs_gnum_t nlogla = 0;
  cs_gnum_t nsubla = 0;
  cs_lnum_t iuiptn = 0;

  cs_real_t alpha_rnn;
  if (   iturb == CS_TURB_RIJ_EPSILON_LRR
      && cs_math_fabs(cs_turb_crij2) <= cs_math_epzero
      && cs_turb_crij1 > 1.0) {
    /* Alpha constant for a realisable BC for R12 with the Rotta model */
    alpha_rnn = 1.0 / sqrt(cs_turb_crij_c0 + 2.0);
  }
  else {
    /* FIXME should be dereve from the algebraic model */
    /* Alpha constant for a realisable BC for R12 with the SSG model */
    alpha_rnn = 0.47;
  }

  /* See the different model */
  const cs_real_t cl = 1.0 / (0.5 + 0.75 * cs_turb_crij_c0);

  /* With v2f type model, (phi-fbar et BL-v2/k) u=0 is set directly, so
     uiptmx and uiptmn are necessarily 0 */
  if (itytur == 5) {
    uiptmx = 0;
    uiptmn = 0;
  }

  /* Pointers to specific fields */
  cs_real_t *byplus = NULL, *bdplus = NULL, *bdlmo = NULL; //*buk = NULL
  BFT_MALLOC(byplus, n_b_faces, cs_real_t);
  BFT_MALLOC(bdplus, n_b_faces, cs_real_t);
  BFT_MALLOC(bdlmo, n_b_faces, cs_real_t);

  // BFT_MALLOC(buk, n_b_faces, cs_real_t);

  /* Correction for atmospheric wall functions */

  cs_field_t *non_neutral_scalar_correction
    = cs_field_by_name_try("non_neutral_scalar_correction");
  cs_real_t *bcfnns = NULL, *bcfnns_loc = NULL;

  if (non_neutral_scalar_correction != NULL) {
    bcfnns = non_neutral_scalar_correction->val;
  }
  else {
    BFT_MALLOC(bcfnns_loc, n_b_faces, cs_real_t);
    bcfnns = bcfnns_loc;
  }

  cs_real_t *cvar_t = NULL;
  cs_real_t *cvar_totwt = NULL;
  cs_real_t *cpro_liqwt = NULL;

  const cs_real_t theta0 = (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) ?
                            tref * pow(pref/p0, r_pg_cnst/cp0) : 0.0;

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {
    cvar_t = f_th->val;

    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2) {
      cvar_totwt = CS_F_(ym_w)->val;
      cpro_liqwt = cs_field_by_name("liquid_water")->val;
    }

  }

  const int *icodcl_vel = vel->bc_coeffs->icodcl;
  cs_real_t *rcodcl1_vel = vel->bc_coeffs->rcodcl1;

  cs_real_t *coftur = NULL, *hfltur = NULL;
  if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
    cs_turbomachinery_get_wall_bc_coeffs(&coftur, &hfltur);
  }

  /* Loop on boundary faces
     ----------------------*/

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test on the presence of a smooth/rough wall condition (start) */
    if (icodcl_vel[f_id] != 5 && icodcl_vel[f_id] != 6)
      continue;

    const cs_lnum_t c_id = b_face_cells[f_id];

    /* Physical properties */
    const cs_real_t visclc = viscl[c_id];
    cs_real_t visctc = visct[c_id];
    const cs_real_t romc = crom[c_id];

    /* Geometric quantities */
    const cs_real_t distbf = b_dist[f_id];
    const cs_real_t *rnxyz = b_face_u_normal[f_id];
    const cs_real_t distfi = b_dist[f_id];

    /* Local reference frame
       --------------------- */

    /* Handle displacement velocity */

    cs_real_t rcodcxyz[3] = {rcodcl1_vel[n_b_faces*0 + f_id],
                             rcodcl1_vel[n_b_faces*1 + f_id],
                             rcodcl1_vel[n_b_faces*2 + f_id]};

    /* If we are not using ALE, force the displacement velocity for the face
       to be tangential (and update rcodcl for possible use)
       In frozen rotor (iturbo = 1), the velocity is neither tangential to the
       wall (absolute velocity solved in a relative frame of reference) */

    if (   cs_glob_ale == CS_ALE_NONE
        && cs_turbomachinery_get_model() == CS_TURBOMACHINERY_NONE) {

      const cs_real_t rcodcn = cs_math_3_dot_product(rcodcxyz, rnxyz);
      rcodcxyz[0] = rcodcxyz[0] - rcodcn * rnxyz[0];
      rcodcxyz[1] = rcodcxyz[1] - rcodcn * rnxyz[1];
      rcodcxyz[2] = rcodcxyz[2] - rcodcn * rnxyz[2];

      rcodcl1_vel[n_b_faces*0 + f_id] = rcodcxyz[0];
      rcodcl1_vel[n_b_faces*1 + f_id] = rcodcxyz[1];
      rcodcl1_vel[n_b_faces*2 + f_id] = rcodcxyz[2];

    }

    /* Relative tangential velocity */
    const cs_real_t upxyz[3] = {velipb[f_id][0] - rcodcxyz[0],
                                velipb[f_id][1] - rcodcxyz[1],
                                velipb[f_id][2] - rcodcxyz[2]};

    const cs_real_t usn = cs_math_3_dot_product(upxyz, rnxyz);

    cs_real_t txyz[3] = {upxyz[0] - usn*rnxyz[0],
                         upxyz[1] - usn*rnxyz[1],
                         upxyz[2] - usn*rnxyz[2]};

    /* Unit tangent (if the velocity is zero, Tx, Ty, Tz is not
       used (we cancel the velocity),
       so we assign any value (zero for example) */
    utau = cs_math_3_norm(txyz);
    cs_math_3_normalize(txyz, txyz);

    /* Complete if necessary for Rij-Epsilon */

    cs_real_t eloglo[3][3], alpha[6][6];

    if (itytur == 3) {

      /* --> T2 = RN X T (where X is the cross product) */

      const cs_real_t t2xyz[3] = {rnxyz[1]*txyz[2] - rnxyz[2]*txyz[1],
                                  rnxyz[2]*txyz[0] - rnxyz[0]*txyz[2],
                                  rnxyz[0]*txyz[1] - rnxyz[1]*txyz[0]};

      /* Orthogonal matrix for change of reference frame ELOGLOij
         (from local to global reference frame)

                  | TX    TY    TZ |
         ELOGLO = |-RNX  -RNY  -RNZ|
                  | T2X   T2Y   T2Z|

         Its transpose ELOGLOt is its inverse */

      eloglo[0][0] =  txyz[0];
      eloglo[1][0] = -rnxyz[0];
      eloglo[2][0] =  t2xyz[0];
      eloglo[0][1] =  txyz[1];
      eloglo[1][1] = -rnxyz[1];
      eloglo[2][1] =  t2xyz[1];
      eloglo[0][2] =  txyz[2];
      eloglo[1][2] = -rnxyz[2];
      eloglo[2][2] =  t2xyz[2];

      /* Compute Reynolds stress transformation matrix */

      int clsyme = 0;
      cs_turbulence_bc_rij_transform(clsyme, eloglo, alpha);

    }

    /* Friction velocities
       =================== */

    /* Compute Uet depending if we are in the log zone or not
       in 1 or 2 velocity scales
       and uk based on ek */

    if (cs_math_fabs(utau) < cs_math_epzero)
      utau = cs_math_epzero;

    const cs_real_t xnuii  = visclc / romc;
    const cs_real_t xnuit  = visctc / romc;

    cs_real_t rttb;
    if (cvar_k != NULL) {
      ek = cvar_k[c_id];
      /* TODO: we could add 2*nu_T dv/dy to rnnb */
      if (icodcl_vel[f_id] == 5)
        rnnb = (2./3.) * ek;
    }
    else if (   turb_model->order == CS_TURB_SECOND_ORDER
             && turb_model->type  == CS_TURB_RANS) {
      ek = 0.5 * (cvar_rij[c_id][0] + cvar_rij[c_id][1] + cvar_rij[c_id][2]);

      rnnb = cs_math_3_sym_33_3_dot_product(rnxyz, cvar_rij[c_id], rnxyz);
      rttb = cs_math_3_sym_33_3_dot_product(txyz, cvar_rij[c_id], txyz);
    }

    const cs_real_t rough_d = (rough != NULL) ? bpro_rough[f_id] : 0;

    int iuntur;
    cs_real_t uk, ypup, dplus, yplus;

    if (icodcl_vel[f_id] == 5) {

      cs_wall_f_type_t iwallf_loc = cs_glob_wall_functions->iwallf;
      if (fvq->has_disable_flag) {
        cs_lnum_t cell_id = cs_glob_mesh->b_face_cells[f_id];
        if (fvq->c_disable_flag[cell_id])
          iwallf_loc = CS_WALL_F_DISABLED;
      }

      cs_wall_functions_velocity(iwallf_loc,
                                 xnuii,
                                 xnuit,
                                 utau,
                                 distbf,
                                 rough_d,
                                 rnnb,
                                 ek,
                                 &iuntur,
                                 &nsubla,
                                 &nlogla,
                                 &uet,
                                 &uk,
                                 &yplus,
                                 &ypup,
                                 &cofimp,
                                 &dplus);

    }
    else if (icodcl_vel[f_id] == 6) {

      /* Neutral value, might be overwritten after */
      uk = cs_turb_cmu025*sqrt(ek);

      /* NB: for rough walls, yplus is computed from the roughness and not uk */
      assert(rough != NULL);
      yplus = distbf / rough_d;

    }

    /* Louis or Monin Obukhov wall function for atmospheric flows */

    cs_real_t dlmo = 0, yk = 0;

    if (iwalfs != CS_WALL_F_S_MONIN_OBUKHOV) {

      /* rough wall */
      if (icodcl_vel[f_id] == 6) {

        /* ustar for neutral, may be modified after */
        uet = utau / log(yplus+1.0) * xkappa;

        /* Dimensionless velocity, neutral wall function,
           may be modified after */
        const cs_real_t _uplus = log(yplus+1.0) / xkappa;

        /* Atmospheric Louis wall functions for rough wall */
        if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

          const cs_real_t gredu = cs_math_3_dot_product(gxyz, rnxyz);
          const cs_real_t temp = cvar_t[c_id];
          cs_real_t totwt = 0.;
          cs_real_t liqwt = 0.;

          if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2) {
            totwt = cvar_totwt[c_id];
            liqwt = cpro_liqwt[c_id];
          }

          /* 1/U+ for neutral */
          const cs_real_t duplus = 1.0 / _uplus;

          const cs_real_t brough_t = bpro_rough_t[f_id];
          const cs_real_t yplus_t = distbf / brough_t;

          /* 1/T+ for neutral */
          const cs_real_t dtplus
            = xkappa / log((distbf + brough_t) / brough_t);

          _atmo_cls(f_id,
                    utau,
                    rough_d,
                    duplus,
                    dtplus,
                    yplus_t,
                    &uet,
                    gredu,
                    &cfnns,
                    &cfnnk,
                    &cfnne,
                    &dlmo,
                    temp,
                    totwt,
                    liqwt);

        }

      }
      /* Louis for the smooth wall case */
      else if (   iwalfs == CS_WALL_F_S_LOUIS && icodcl_vel[f_id] == 5
               && cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

        /* Compute reduced gravity for non horizontal walls */
        const cs_real_t gredu = cs_math_3_dot_product(gxyz, rnxyz);
        const cs_real_t temp = cvar_t[c_id];
        cs_real_t totwt = 0.;
        cs_real_t liqwt = 0.;

        if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2) {
          totwt = cvar_totwt[c_id];
          liqwt = cpro_liqwt[c_id];
        }

        yk = distbf * uk / xnuii;
        /* 1/U+ for neutral */
        const cs_real_t duplus = ypup / yk;
        const cs_real_t brough_t = bpro_rough_t[f_id];

        /* 1/T+
         * "y+_t" tends to "y/rough_t" for rough regime and to "y+k"
         * times a shift for smooth regime
         *
         * Rough regime reads:
         *   T+ = 1/kappa ln(y/rough_t) = 1/kappa ln(y/zeta) + 8.5
         *                              = 1/kappa ln[y/zeta * exp(8.5 kappa)]
         *
         * Note zeta_t = rough_t * exp(8.5 kappa)
         *
         * Smooth regime reads:
         * T+ = 1/kappa ln(y uk/nu) + 5.2 = 1/kappa ln[y uk*exp(5.2 kappa) / nu]
         *
         * Mixed regime reads:
         *   T+ = 1/kappa ln[y uk*exp(5.2 kappa)/(nu + alpha uk zeta)]
         *      = 1/kappa ln[  y uk*exp(5.2 kappa)
         *                   / (nu + alpha uk rough_t * exp(8.5 kappa))]
         *      = 1/kappa ln[  y uk*exp(5.2 kappa)
         *                   / (nu + alpha uk rough_t * exp(8.5 kappa))]
         * with
         *   alpha * exp(8.5 kappa) / exp(5.2 kappa) = 1
         * ie
         *   alpha = exp(-(8.5-5.2) kappa) = 0.25
         * so
         *   T+ = 1/kappa ln[  y uk*exp(5.2 kappa)
         *                   / (nu + uk rough_t * exp(5.2 kappa))]
         *      = 1/kappa ln[y+k / (exp(-5.2 kappa) + uk rough_t/nu)]
         */

        /* shifted y+ */
        /* FIXME use log constant */
        const cs_real_t yplus_t
          = yk / (exp(-xkappa * 5.2) + uk *brough_t / xnuii);
        /* 1/T+ for neutral */
        const cs_real_t dtplus = xkappa / log(yplus_t);

        _atmo_cls(f_id,
                  utau,
                  rough_d,
                  duplus,
                  dtplus,
                  yplus_t,
                  &uet,
                  gredu,
                  &cfnns,
                  &cfnnk,
                  &cfnne,
                  &dlmo,
                  temp,
                  totwt,
                  liqwt);

      }

    }
    /* Monin Obukhov wall function for smooth and rough wall */
    else if (iwalfs == CS_WALL_F_S_MONIN_OBUKHOV) {

      /* Compute local LMO */
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

        const cs_real_t gredu = cs_math_3_dot_product(gxyz, rnxyz);

        const int *icodcl_th = f_th->bc_coeffs->icodcl;

        if (   icodcl_th[f_id] == 6
            || (icodcl_th[f_id] == 5 && icodcl_vel[f_id] == 5)) {

          const cs_real_t *rcodcl1_th = f_th->bc_coeffs->rcodcl1;
          const cs_real_t dt = theipb[f_id]-rcodcl1_th[f_id];

          cs_f_mo_compute_from_thermal_diff(&distbf,
                                            &rough_d,
                                            &utau,
                                            &dt,
                                            &theta0,
                                            &gredu,
                                            &dlmo,
                                            &uet);

        }
        else if (icodcl_th[f_id] == 3) {

          const cs_real_t *rcodcl3_th = f_th->bc_coeffs->rcodcl3;
          const cs_real_t cpp = (icp >= 0) ? cpro_cp[c_id] : cp0;
          const cs_real_t flux = rcodcl3_th[f_id] / romc / cpp;

          cs_f_mo_compute_from_thermal_flux(&distbf,
                                            &rough_d,
                                            &utau,
                                            &flux,
                                            &theta0,
                                            &gredu,
                                            &dlmo,
                                            &uet);
        }

      }
      else {

        /* No temperature delta: neutral */
        const cs_real_t dt = 0., _theta0 = 0., gredu = 0.;

        cs_f_mo_compute_from_thermal_diff(&distbf,
                                          &rough_d,
                                          &utau,
                                          &dt,
                                          &_theta0,
                                          &gredu,
                                          &dlmo,
                                          &uet);

      }

      /* Take stability into account for the turbulent velocity scale */
      cs_real_t coef_mom = cs_mo_phim(distbf + rough_d, dlmo);
      const cs_real_t one_minus_ri
        = 1 - (distbf + rough_d) * dlmo / coef_mom;

      if (one_minus_ri > 0) {
        /* Warning: overwritting uk, yplus should be recomputed */
        uk = uk / pow(one_minus_ri, 0.25);
        yplus = distbf * uk / xnuii;

        /* Epsilon should be modified as well to get
           P+G = P(1-Ri) = epsilon
           P = -R_tn dU/dn = uk^2 uet Phi_m / (kappa z) */
        cfnne = one_minus_ri * coef_mom;
        /* Nothing done for the moment for really high stability */
      }
      else {
        cfnne = 1.;
      }

      if (icodcl_vel[f_id] == 5) {
        /* Boundary condition on the velocity to have approximately
           the correct turbulence production */
        coef_mom = cs_mo_phim(distbf+rough_d, dlmo);
        const cs_real_t coef_momm = cs_mo_phim(2 * distbf + rough_d, dlmo);
        rcprod =   2*distbf*sqrt(  xkappa*uk*romc*coef_mom/visctc
                                   / (distbf+rough_d))
                 - coef_momm / (2.0 + rough_d / distbf);

        iuntur = 1;

        const cs_real_t _uplus = utau / uet;
        /* Coupled solving of the velocity components
           The boundary term for velocity gradient is implicit
           modified for non neutral boundary layer (in uplus) */
        cofimp  = cs_math_fmin(cs_math_fmax(1-1/(xkappa*_uplus)*rcprod, 0),
                               1);
        yk = distbf * uk / xnuii;

      }

    } /* End Monin Obukhov */

    /* Dimensionless velocity, recomputed and therefore may
       take stability into account */

    cs_real_t uplus = 0.0;
    if (   cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1
        && (iwalfs == 2 || iwalfs == 3) && icodcl_vel[f_id] == 5) {

      uplus = utau / uet;

      /* y+/U+ for non neutral is recomputed */
      ypup = yk / cs_math_fmax(uplus, cs_math_epzero);
    }
    else if (icodcl_vel[f_id] == 6)
      uplus = utau / uet;

    /* Store u_star and uk if needed */
    if (boundary_ustar != NULL && icodcl_vel[f_id] == 5)
      bpro_ustar[f_id] = uet;

    /* Rough wall: one velocity scale: set uk to uet */
    if (cs_glob_wall_functions->iwallf <= 2 && icodcl_vel[f_id] == 6)
      uk = uet;

    if(boundary_uk != NULL && icodcl_vel[f_id] == 5)
      bpro_uk[f_id] = uk;

    /* Save yplus if post-processed or condensation modelling */
    if (f_yplus != NULL)
      yplbr[f_id] = yplus;

    uetmax  = cs_math_fmax(uet, uetmax);
    uetmin  = cs_math_fmin(uet, uetmin);
    ukmax   = cs_math_fmax(uk, ukmax);
    ukmin   = cs_math_fmin(uk, ukmin);
    yplumx  = cs_math_fmax(yplus, yplumx);
    yplumn  = cs_math_fmin(yplus, yplumn);
    dlmomin = cs_math_fmin(dlmo, dlmomin);
    dlmomax = cs_math_fmax(dlmo, dlmomax);

    /* Save turbulent subgrid viscosity after van Driest damping in LES
       care is taken to not dampen it twice at boundary cells having more
       than one boundary face */
    if (itytur == 4 && cs_glob_turb_les_model->idries == 1) {
      if (visvdr[c_id] < -900.) {
        if (icodcl_vel[f_id] == 5)
          visct[c_id] =   visct[c_id]
                        * cs_math_pow2(1.0 - exp(-yplus/cs_turb_cdries));

        visvdr[c_id] = visct[c_id];
        visctc = visct[c_id];
      }
    }

    /* Velocity boundary conditions
       ============================ */

    cs_real_t hflui = 0;
    if (icodcl_vel[f_id] == 5)
      _update_physical_quantities_smooth_wall(c_id,
                                              visclc,
                                              visctc,
                                              romc,
                                              distbf,
                                              utau,
                                              uet,
                                              uk,
                                              yplus,
                                              ypup,
                                              dplus,
                                              &cofimp,
                                              &rcprod,
                                              &hflui,
                                              &uiptn);
    else if (icodcl_vel[f_id] == 6)
      _update_physical_quantities_rough_wall(visclc,
                                             visctc,
                                             romc,
                                             distbf,
                                             utau,
                                             uet,
                                             uk,
                                             uplus,
                                             rough_d,
                                             dlmo,
                                             &iuntur,
                                             &nlogla,
                                             &nsubla,
                                             &cofimp,
                                             &rcprod,
                                             &hflui,
                                             &uiptn);

    /* Min and Max and counter of reversal layer */
    uiptmn = cs_math_fmin(uiptn * iuntur, uiptmn);
    uiptmx = cs_math_fmax(uiptn * iuntur, uiptmx);

    if (uiptn * iuntur < - cs_math_epzero)
      iuiptn = iuiptn + 1;

    const cs_real_t hintv = (itytur == 3) ?
                            visclc / distbf:
                            (visclc + visctc) / distbf;

    /* Gradient boundary conditions
       ---------------------------- */

    const cs_real_t rcodcn = cs_math_3_dot_product(rcodcxyz, rnxyz);

    coefa_vel[f_id][0] =   (1.0 - cofimp) * (rcodcxyz[0] - rcodcn*rnxyz[0])
                         + rcodcn*rnxyz[0];

    coefa_vel[f_id][1] =   (1.0 - cofimp) * (rcodcxyz[1] - rcodcn*rnxyz[1])
                         + rcodcn*rnxyz[1];

    coefa_vel[f_id][2] =   (1.0 - cofimp) * (rcodcxyz[2] - rcodcn*rnxyz[2])
                         + rcodcn*rnxyz[2];

    /* Projection in order to have the velocity parallel to the wall
       B = cofimp * ( IDENTITY - n x n ) */

    coefb_vel[f_id][0][0] =   cofimp * (1.0 - rnxyz[0] * rnxyz[0]);
    coefb_vel[f_id][1][1] =   cofimp * (1.0 - rnxyz[1] * rnxyz[1]);
    coefb_vel[f_id][2][2] =   cofimp * (1.0 - rnxyz[2] * rnxyz[2]);
    coefb_vel[f_id][0][1] = - cofimp * rnxyz[0] * rnxyz[1];
    coefb_vel[f_id][0][2] = - cofimp * rnxyz[0] * rnxyz[2];
    coefb_vel[f_id][1][2] = - cofimp * rnxyz[1] * rnxyz[2];
    coefb_vel[f_id][1][0] =   coefb_vel[f_id][0][1];
    coefb_vel[f_id][2][1] =   coefb_vel[f_id][1][2];
    coefb_vel[f_id][2][0] =   coefb_vel[f_id][0][2];

    /* Flux boundary conditions
       ------------------------ */

    cofaf_vel[f_id][0] = - hflui*(rcodcxyz[0] - rcodcn*rnxyz[0])
                         - hintv*rcodcn*rnxyz[0];

    cofaf_vel[f_id][1] = - hflui*(rcodcxyz[1] - rcodcn*rnxyz[1])
                         - hintv*rcodcn*rnxyz[1];

    cofaf_vel[f_id][2] = - hflui*(rcodcxyz[2] - rcodcn*rnxyz[2])
                         - hintv*rcodcn*rnxyz[2];

    /* Projection in order to have the shear stress parallel to the wall
       B = hflui*( IDENTITY - n x n ) */

    cofbf_vel[f_id][0][0] =   hflui*(1.-rnxyz[0]*rnxyz[0])
                            + hintv*rnxyz[0]*rnxyz[0];
    cofbf_vel[f_id][1][1] =   hflui*(1.-rnxyz[1]*rnxyz[1])
                            + hintv*rnxyz[1]*rnxyz[1];
    cofbf_vel[f_id][2][2] =   hflui*(1.-rnxyz[2]*rnxyz[2])
                            + hintv*rnxyz[2]*rnxyz[2];

    cofbf_vel[f_id][0][1] = (hintv - hflui)*rnxyz[0]*rnxyz[1];
    cofbf_vel[f_id][0][2] = (hintv - hflui)*rnxyz[0]*rnxyz[2];
    cofbf_vel[f_id][1][2] = (hintv - hflui)*rnxyz[1]*rnxyz[2];

    cofbf_vel[f_id][1][0] = cofbf_vel[f_id][0][1];
    cofbf_vel[f_id][2][0] = cofbf_vel[f_id][0][2];
    cofbf_vel[f_id][2][1] = cofbf_vel[f_id][1][2];

    /* In case of transient turbomachinery computations, save the coefficents
       associated to turbulent wall velocity BC, in order to update the wall
       velocity after the geometry update
       (between prediction and correction step) */

    if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
      const int *irotce = cs_turbomachinery_get_cell_rotor_num();
      if (irotce[c_id] != 0) {
        coftur[f_id] = cofimp;
        hfltur[f_id] = hflui;
      }
    }

    /* Boundary conditions on k and epsilon
       ==================================== */

    const cs_real_t ydep = 0.5 * distbf + rough_d;

    if (itytur == 2) {

      /* Launder Sharma boundary conditions
         ================================== */

      if (iturb == CS_TURB_K_EPSILON_LS && icodcl_vel[f_id] == 5) {

        cs_real_t *coefa_ep = f_eps->bc_coeffs->a;
        cs_real_t *coefb_ep = f_eps->bc_coeffs->b;
        cs_real_t *cofaf_ep = f_eps->bc_coeffs->af;
        cs_real_t *cofbf_ep = f_eps->bc_coeffs->bf;

        cs_real_t *coefa_k = f_k->bc_coeffs->a;
        cs_real_t *coefb_k = f_k->bc_coeffs->b;
        cs_real_t *cofaf_k = f_k->bc_coeffs->af;
        cs_real_t *cofbf_k = f_k->bc_coeffs->bf;

        /* Dirichlet Boundary Condition on k
           --------------------------------- */

        cs_real_t pimp = 0.0, hint = 0.0;
        if (cs_glob_wall_functions->iwallf == 0)
          /* No wall functions forced by user */
          pimp = 0.;
        else
          /* Use of wall functions */
          pimp = (iuntur == 1) ? uk*uk/sqrcmu : 0.0;

        pimp = pimp * cfnnk;
        hint = (visclc + visctc / sigmak) / distbf;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_k[f_id],
                                                    &cofaf_k[f_id],
                                                    &coefb_k[f_id],
                                                    &cofbf_k[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        /* Dirichlet Boundary Condition on epsilon tilda
           --------------------------------------------- */
        const cs_real_t pimp_lam = 0.;

        if (cs_glob_wall_functions->iwallf == 0) {
          /* No wall functions forced by user */
          pimp = pimp_lam;
        }
        else {
          /* Use of wall functions */
          if (yplus > cs_math_epzero) {
            const cs_real_t pimp_turb
              = 5 * pow(uk, 4) * romc / (xkappa * visclc * yplus);

            /* Blending function, from JF Wald PhD (2016) */
            const cs_real_t fct_bl = exp( - 0.674e-3 * pow(yplus, 3) );
            const cs_real_t fep    = exp( - pow(0.25 * (yplus + dplus), 1.5) );
            const cs_real_t dep    = 1.0 - exp(-pow((yplus + dplus)/9.0, 2.1));

            /* Je comprend pas: pimp est calculé à partir de fct_bl
               puis re-calculé differemment avec fep et dep */
            pimp = pimp_lam * fct_bl + pimp_turb * (1 - fct_bl);
            pimp = fep * pimp_lam + (1.0 - fep) * dep * pimp_turb;
          }
          else
            pimp = pimp_lam;

        }

        hint = (visclc+visctc/sigmae)/distbf;
        pimp = pimp * cfnne;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_ep[f_id],
                                                    &cofaf_ep[f_id],
                                                    &coefb_ep[f_id],
                                                    &cofbf_ep[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        /* if defined set Dirichlet condition for the Lagrangian time scale */

        if (f_tlag != NULL) {

          cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
          cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
          cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
          cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

          if (cs_glob_wall_functions->iwallf == 0) {
            /* No wall functions forced by user */
            pimp = 0.;
          }
          else {
            /* Use of wall functions */
            if (iuntur == 1)
              pimp =   cfnnk / (cfnne * uk) * cl / sqrcmu * xkappa
                     * (dplus * visclc / (romc * uk) + rough_d);
            else
              pimp = 0.;
          }

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                      &cofaf_tlag[f_id],
                                                      &coefb_tlag[f_id],
                                                      &cofbf_tlag[f_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);
        }

      }

      /* Quadratic Baglietto k-epsilon model
         =================================== */

      else if (iturb == CS_TURB_K_EPSILON_QUAD && icodcl_vel[f_id] == 5) {

        cs_real_t *coefa_ep = f_eps->bc_coeffs->a;
        cs_real_t *coefb_ep = f_eps->bc_coeffs->b;
        cs_real_t *cofaf_ep = f_eps->bc_coeffs->af;
        cs_real_t *cofbf_ep = f_eps->bc_coeffs->bf;

        cs_real_t *coefa_k = f_k->bc_coeffs->a;
        cs_real_t *coefb_k = f_k->bc_coeffs->b;
        cs_real_t *cofaf_k = f_k->bc_coeffs->af;
        cs_real_t *cofbf_k = f_k->bc_coeffs->bf;

        /* Dirichlet Boundary Condition on k
           ---------------------------------- */

        cs_real_t pimp = 0.0, hint = 0.0;
        if (cs_glob_wall_functions->iwallf == 0)
          /* No wall functions forces by user */
          pimp = 0.;
        else
          /* Use of wall functions */
          pimp = (iuntur == 1) ? uk*uk/sqrcmu : 0.0;

        hint = (visclc + visctc / sigmak) / distbf;
        pimp = pimp * cfnnk;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_k[f_id],
                                                    &cofaf_k[f_id],
                                                    &coefb_k[f_id],
                                                    &cofbf_k[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        /* Dirichlet Boundary Condition on epsilon
           --------------------------------------- */

        if (cs_glob_wall_functions->iwallf != 0) {

          const cs_real_t pimp_lam
            = 2 * visclc / romc * cvar_k[c_id] / (distbf * distbf);

          if (yplus > cs_math_epzero) {
            const cs_real_t pimp_turb = 5 * pow(uk, 4) * romc
                                        / (xkappa * visclc * yplus);

            /* Blending between wall and homogeneous layer */
            const cs_real_t fep  = exp(- pow(0.25 * (yplus + dplus), 1.5));
            const cs_real_t dep  = 1.0 - exp(- pow((yplus + dplus) / 9.0, 2.1));
            pimp = fep * pimp_lam + (1.0 - fep) * dep * pimp_turb;
          }
          else
            pimp = pimp_lam;
          }
          else
            pimp = 2 * visclc / romc * cvar_k[c_id] / (distbf * distbf);

        pimp = pimp * cfnne;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_ep[f_id],
                                                    &cofaf_ep[f_id],
                                                    &coefb_ep[f_id],
                                                    &cofbf_ep[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        /* if defined set Dirichlet condition for the Lagrangian time scale */

        if (f_tlag != NULL) {

          cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
          cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
          cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
          cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

          if (cs_glob_wall_functions->iwallf == 0)
            /* No wall functions forced by user*/
            pimp = 0.;
          else {
            /* Use of wall functions */
            if (iuntur == 1)
              pimp = cfnnk / (cfnne * uk) * cl / sqrcmu * xkappa
                * (dplus * visclc / (romc * uk) + rough_d);
            else
              pimp = 0.;
          }

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                      &cofaf_tlag[f_id],
                                                      &coefb_tlag[f_id],
                                                      &cofbf_tlag[f_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);
        }

      }

      /* k-epsilon and k-epsilon LP boundary conditions
         ============================================== */

      else {

        /* Dirichlet Boundary Condition on k
           --------------------------------- */

        cs_real_t *coefa_ep = f_eps->bc_coeffs->a;
        cs_real_t *coefb_ep = f_eps->bc_coeffs->b;
        cs_real_t *cofaf_ep = f_eps->bc_coeffs->af;
        cs_real_t *cofbf_ep = f_eps->bc_coeffs->bf;

        cs_real_t *coefa_k = f_k->bc_coeffs->a;
        cs_real_t *coefb_k = f_k->bc_coeffs->b;
        cs_real_t *cofaf_k = f_k->bc_coeffs->af;
        cs_real_t *cofbf_k = f_k->bc_coeffs->bf;

        cs_real_t qimp = 0.0;
        cs_real_t pimp = (iuntur == 1 || icodcl_vel[f_id] == 6) ?
                         uk*uk*cfnnk/sqrcmu : 0.0;

        cs_real_t hint = (visclc + visctc / sigmak) / distbf;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_k[f_id],
                                                    &cofaf_k[f_id],
                                                    &coefb_k[f_id],
                                                    &cofbf_k[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        if (icodcl_vel[f_id] == 6 && df_limiter_k != NULL)
          df_limiter_k[c_id] = 0.0;

        /* Neumann Boundary Condition on epsilon
           ------------------------------------- */

        hint = (visclc + visctc / sigmae) / distbf;

        /* If yplus=0, uiptn is set to 0 to avoid division by 0.
           By the way, in this case: iuntur=0 */
        if (yplus > cs_math_epzero && iuntur == 1) { /* FIXME use only iuntur */
          pimp =   distbf * 4 * pow(uk, 5)
                 / (xkappa * xnuii *xnuii * cs_math_pow2(yplus+2*dplus));

          qimp = - pimp * hint; /* TODO transform it,
                                   it is only to be fully equivalent */
        }
        else {
          qimp = 0.;
        }

        pimp = pimp * cfnne;

        if (icodcl_vel[f_id] == 6) {

          pimp =   cs_math_pow3(uk) / (xkappa * ydep * ydep) * distbf * cfnne;
          qimp = - pimp * hint;
          /*TODO transform it to use d eps / d y directly */

        }

        cs_boundary_conditions_set_neumann_scalar(&coefa_ep[f_id],
                                                  &cofaf_ep[f_id],
                                                  &coefb_ep[f_id],
                                                  &cofbf_ep[f_id],
                                                  qimp,
                                                  hint);

        /* If defined set Dirichlet condition for the Lagrangian time scale */

        if (f_tlag != NULL) {

          cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
          cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
          cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
          cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

          if (cs_glob_wall_functions->iwallf == 0)
            /* No wall functions forced by user */
            pimp = 0.;
          else {
            /* Use of wall functions */
            if (iuntur == 1) {

              if (icodcl_vel[f_id] == 5) {

                pimp = cfnnk / (cfnne * uk) * cl / sqrcmu * xkappa
                  * (dplus * visclc / (romc * uk) + rough_d);
              }
              else if (icodcl_vel[f_id] == 6)
                pimp = cfnnk / (cfnne * uk) * cl / sqrcmu * xkappa * rough_d;
            }
            else {
              pimp = 0.;
            }
          }

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                      &cofaf_tlag[f_id],
                                                      &coefb_tlag[f_id],
                                                      &cofbf_tlag[f_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);
        }

        if (icodcl_vel[f_id] == 6 && df_limiter_eps != NULL)
          df_limiter_eps[c_id] = 0.0;
      }
    }

    /* Boundary conditions on Rij-epsilon
       ================================== */

    else if (itytur == 3) {

      cs_real_t visci[3][3], dist[3], hint = 0.0;

      cs_real_t *coefa_ep = f_eps->bc_coeffs->a;
      cs_real_t *coefb_ep = f_eps->bc_coeffs->b;
      cs_real_t *cofaf_ep = f_eps->bc_coeffs->af;
      cs_real_t *cofbf_ep = f_eps->bc_coeffs->bf;

      cs_real_6_t  *coefa_rij = (cs_real_6_t  *)f_rij->bc_coeffs->a;
      cs_real_66_t *coefb_rij = (cs_real_66_t *)f_rij->bc_coeffs->b;
      cs_real_6_t  *cofaf_rij = (cs_real_6_t  *)f_rij->bc_coeffs->af;
      cs_real_66_t *cofbf_rij = (cs_real_66_t *)f_rij->bc_coeffs->bf;
      cs_real_6_t  *cofad_rij = (cs_real_6_t  *)f_rij->bc_coeffs->ad;
      cs_real_66_t *cofbd_rij = (cs_real_66_t *)f_rij->bc_coeffs->bd;

      /* Exchange coefficient */

      /* Symmetric tensor diffusivity (Daly Harlow -- GGDH) */
      if (eqp_rij->idften & CS_ANISOTROPIC_RIGHT_DIFFUSION) {

        visci[0][0] = visclc + visten[c_id][0];
        visci[1][1] = visclc + visten[c_id][1];
        visci[2][2] = visclc + visten[c_id][2];
        visci[0][1] =          visten[c_id][3];
        visci[1][0] =          visten[c_id][3];
        visci[1][2] =          visten[c_id][4];
        visci[2][1] =          visten[c_id][4];
        visci[0][2] =          visten[c_id][5];
        visci[2][0] =          visten[c_id][5];

        dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
        dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
        dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

        /* ||Ki.n||^2 */
        const cs_real_t viscis = cs_math_pow2(  visci[0][0]*rnxyz[0]
                                              + visci[1][0]*rnxyz[1]
                                              + visci[2][0]*rnxyz[2])
                               + cs_math_pow2(  visci[0][1]*rnxyz[0]
                                              + visci[1][1]*rnxyz[1]
                                              + visci[2][1]*rnxyz[2])
                               + cs_math_pow2(  visci[0][2]*rnxyz[0]
                                              + visci[1][2]*rnxyz[1]
                                              + visci[2][2]*rnxyz[2]);

        /* IF.Ki.n */
        cs_real_t fikis
          = (  cs_math_3_dot_product(dist, visci[0]) * rnxyz[0]
             + cs_math_3_dot_product(dist, visci[1]) * rnxyz[1]
             + cs_math_3_dot_product(dist, visci[2]) * rnxyz[2]);

        /* Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
           NB: eps =1.d-1 must be consistent with vitens.f90 */
        fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distfi);

        hint = viscis / fikis;
      }
      /* Scalar diffusivity */
      else
        hint = (visclc + visctc * cs_turb_csrij / cs_turb_cmu) / distbf;

      /* ---> Tensor Rij (Partially or totally implicited) */

      cs_real_t fcoefa[6], fcoefb[6], fcofad[6], fcofbd[6];
      cs_real_t fcofaf[6], fcofbf[6];
      for (int isou = 0; isou < 6; isou++) {
        fcoefa[isou] = 0.;
        fcoefb[isou] = 0.;
        fcofad[isou] = 0.;
        fcofbd[isou] = 0.;
      }

      /* blending factor so that the component R(n,tau) have only
         -mu_T/(mu+mu_T)*uet*uk */
      const cs_real_t bldr12 = (icodcl_vel[f_id] == 5) ?
                                visctc / (visclc + visctc) : 1.0;

      for (int isou = 0; isou < 6; isou++) {

        cs_lnum_t jj = _iv2t[isou];
        cs_lnum_t kk = _jv2t[isou];

        /* LRR and the Standard SGG or EB-RSM + wall functions */
        if (      ((iuntur == 1)
               && (   iturb == CS_TURB_RIJ_EPSILON_LRR
                   || iturb == CS_TURB_RIJ_EPSILON_SSG))
            || (   iturb == CS_TURB_RIJ_EPSILON_EBRSM
                && cs_glob_wall_functions->iwallf != 0
                && yplus > cs_math_epzero)
            || icodcl_vel[f_id] == 6) {

          if (cs_glob_turb_rans_model->irijco == 1) {

            coefa_rij[f_id][isou] = - (  eloglo[0][jj]*eloglo[1][kk]
                                       + eloglo[1][jj]*eloglo[0][kk])
                                  * alpha_rnn * sqrt(rnnb * rttb) * cfnnk;

            cofaf_rij[f_id][isou] = - hint * coefa_rij[f_id][isou];
            cofad_rij[f_id][isou] = 0.0;

            for (int ii = 0; ii < 6; ii++) {
              coefb_rij[f_id][isou][ii] = alpha[ii][isou];

              if (ii == isou)
                cofbf_rij[f_id][isou][ii]
                  = hint * (1.0 - coefb_rij[f_id][isou][ii]);
              else
                cofbf_rij[f_id][isou][ii]
                  = - hint * coefb_rij[f_id][isou][ii];

              cofbd_rij[f_id][isou][ii] = coefb_rij[f_id][isou][ii];
            }
          }
          else if ((cs_glob_turb_rans_model->iclptr == 1)) {

            for (int ii = 0; ii < 6; ii++) {
              if (ii != isou)
                fcoefa[isou] += alpha[ii][isou] * rijipb[f_id][ii];
            }
            fcoefb[isou] = alpha[isou][isou];

          }
          else {
            for (int ii = 0; ii < 6; ii++)
              fcoefa[isou] = fcoefa[isou] + alpha[ii][isou] * rijipb[f_id][ii];

            fcoefb[isou] = 0.0;
          }

          /* Boundary conditions for the momentum equation */
          fcofad[isou] = fcoefa[isou];
          fcofbd[isou] = fcoefb[isou];

          /*fcoefa[isou] =   fcoefa[isou] - (eloglo[jj][0]*eloglo[kk][1]
            + eloglo[jj][1]*eloglo[kk][0]) * bldr12*uet*uk*cfnnk;*/

          fcoefa[isou] =   fcoefa[isou] - (eloglo[0][jj]*eloglo[1][kk]
            + eloglo[1][jj]*eloglo[0][kk]) * bldr12*uet*uk*cfnnk;

          /* Translate into Diffusive flux BCs for rough wall */
          if (icodcl_vel[f_id] == 6) {
            fcofaf[isou] = - hint * fcoefa[isou];
            fcofbf[isou] =   hint * (1.0 - fcoefb[isou]);
          }

        }

        /* In the viscous sublayer or for EBRSM: zero Reynolds' stresses
           (only for smooth wall) */
        else {
          if (cs_glob_turb_rans_model->irijco == 1) {
            coefa_rij[f_id][isou] = 0.;
            cofaf_rij[f_id][isou] = 0.;
            cofad_rij[f_id][isou] = 0.;
            for (int ii = 0; ii < 6; ii++) {
              coefb_rij[f_id][isou][ii]  = 0.;

              if (ii == isou)
                cofbf_rij[f_id][isou][ii] = hint;
              else
                cofbf_rij[f_id][isou][ii] = 0.;

              cofbd_rij[f_id][isou][ii] = 0.;
            }
          }
          else {
            fcoefa[isou] = 0.;
            fcofad[isou] = 0.;
            fcoefb[isou] = 0.;
            fcofbd[isou] = 0.;
          }
        }

        /* Translate into Diffusive flux BCs */
        fcofaf[isou] = - hint * fcoefa[isou];
        fcofbf[isou] =   hint * (1.0 - fcoefb[isou]);

      } /* End loop on isou */

      if (cs_glob_turb_rans_model->irijco != 1) {

        for (int isou = 0; isou < 6; isou++) {
          coefa_rij[f_id][isou] = fcoefa[isou];
          cofaf_rij[f_id][isou] = fcofaf[isou];
          cofad_rij[f_id][isou] = fcofad[isou];

          for (int ii = 0; ii < 6; ii++) {
            coefb_rij[f_id][isou][ii] = 0;
            cofbd_rij[f_id][isou][ii] = 0;
          }

          coefb_rij[f_id][isou][isou] = fcoefb[isou];
          cofbf_rij[f_id][isou][isou] = fcofbf[isou];
          cofbd_rij[f_id][isou][isou] = fcofbd[isou];
        }

      }

      if (icodcl_vel[f_id] == 6 && df_limiter_rij != NULL)
        df_limiter_rij[c_id] = 0.0;

      /* Epsilon
         NB: no reconstruction, possibility of partial implicitation */

      /* Symmetric tensor diffusivity (Daly Harlow -- GGDH) */
      if (eqp_eps->idften & CS_ANISOTROPIC_DIFFUSION) {

        visci[0][0] = visclc + visten[c_id][0] / sigmae;
        visci[1][1] = visclc + visten[c_id][1] / sigmae;
        visci[2][2] = visclc + visten[c_id][2] / sigmae;
        visci[0][1] =          visten[c_id][3] / sigmae;
        visci[1][0] =          visten[c_id][3] / sigmae;
        visci[1][2] =          visten[c_id][4] / sigmae;
        visci[2][1] =          visten[c_id][4] / sigmae;
        visci[0][2] =          visten[c_id][5] / sigmae;
        visci[2][0] =          visten[c_id][5] / sigmae;

        /* ||Ki.S||^2 */
        const cs_real_t viscis = cs_math_pow2(  visci[0][0]*rnxyz[0]
                                              + visci[1][0]*rnxyz[1]
                                              + visci[2][0]*rnxyz[2])
                               + cs_math_pow2(  visci[0][1]*rnxyz[0]
                                              + visci[1][1]*rnxyz[1]
                                              + visci[2][1]*rnxyz[2])
                               + cs_math_pow2(  visci[0][2]*rnxyz[0]
                                              + visci[1][2]*rnxyz[1]
                                              + visci[2][2]*rnxyz[2]);

        /* if.ki.s */
        cs_real_t fikis
          = (  cs_math_3_dot_product(dist, visci[0]) * rnxyz[0]
             + cs_math_3_dot_product(dist, visci[1]) * rnxyz[1]
             + cs_math_3_dot_product(dist, visci[2]) * rnxyz[2]);

        /* take i" so that i"f= eps*||fi||*ki.n when j" is in cell rji
             nb: eps =1.d-1 must be consistent with vitens.f90 */
          fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distbf);

          hint = viscis / fikis;
      }

      /* Scalar diffusivity */
      else {
        hint = (visclc + visctc / sigmae) / distbf;
      }

      if (iturb == CS_TURB_RIJ_EPSILON_LRR || iturb == CS_TURB_RIJ_EPSILON_SSG
          || (itytur == 3 && icodcl_vel[f_id] == 6)) {

        /* Si yplus=0, on met coefa a 0 directement pour eviter une division
           par 0 */
        /* Compute pimp for smooth wall */
        cs_real_t pimp = 0.;
        if (yplus > cs_math_epzero && iuntur == 1)
          pimp =   distbf * 4 * pow(uk, 5)
                 / (xkappa * xnuii * xnuii * cs_math_pow2(yplus + 2*dplus));
        else
          pimp = 0.;

        /* Neumann Boundary Condition
           -------------------------- */

        if (cs_glob_turb_rans_model->iclptr == 1 || icodcl_vel[f_id] == 6) {
          /* TODO not available for k-eps */

          /* TODO transform it, it is only to be fully equivalent */
          cs_real_t qimp = - pimp * hint;
          pimp = pimp * cfnne;

          if (icodcl_vel[f_id] == 6) {
            pimp = cs_math_pow3(uk)/(xkappa*ydep*ydep)*distbf*cfnne;

            /* TODO transform it to use d eps / d y directly */
            qimp = - pimp * hint;
          }

          cs_boundary_conditions_set_neumann_scalar(&coefa_ep[f_id],
                                                    &cofaf_ep[f_id],
                                                    &coefb_ep[f_id],
                                                    &cofbf_ep[f_id],
                                                    qimp,
                                                    hint);

        }

        /* Dirichlet Boundary Condition
           ---------------------------- */

        else { /* Only for smooth wall */

          const cs_real_t *cvar_ep = (const cs_real_t *)f_eps->val;

          pimp = pimp + cvar_ep[c_id];
          pimp = pimp * cfnne;

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_ep[f_id],
                                                      &cofaf_ep[f_id],
                                                      &coefb_ep[f_id],
                                                      &cofbf_ep[f_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);

        }

        /* If defined set Dirichlet condition for the Lagrangian time scale */

        if (f_tlag != NULL) {

          cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
          cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
          cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
          cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

          if (cs_glob_wall_functions->iwallf == 0) {
            /* No wall functions forced by user */
            pimp = 0.;
          }
          else {
            /* Use of wall functions */
            if (iuntur == 1) {

              if (icodcl_vel[f_id] == 5) {

                pimp = 0.5 * cfnnk / (cfnne * cs_math_pow3(uk)) * cl * xkappa
                * (  coefa_rij[f_id][0] + coefb_rij[f_id][0][0] *rijipb[f_id][0]
                   + coefa_rij[f_id][1] + coefb_rij[f_id][1][1] *rijipb[f_id][1]
                   + coefa_rij[f_id][2] + coefb_rij[f_id][2][2] *rijipb[f_id][2])
                * (dplus * visclc / (romc * uk) + rough_d);
              }
              else if (icodcl_vel[f_id] == 6) {
                pimp = 0.5 * cfnnk / (cfnne * cs_math_pow3(uk)) * cl * xkappa
                * (  coefa_rij[f_id][0] + coefb_rij[f_id][0][0] *rijipb[f_id][0]
                   + coefa_rij[f_id][1] + coefb_rij[f_id][1][1] *rijipb[f_id][1]
                   + coefa_rij[f_id][2] + coefb_rij[f_id][2][2] *rijipb[f_id][2])
                * rough_d;
              }

            }
            else
              pimp = 0.;
          }

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                      &cofaf_tlag[f_id],
                                                      &coefb_tlag[f_id],
                                                      &cofbf_tlag[f_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);
        }
      }

      /* process only for smooth wall here after */
      else if (iturb == CS_TURB_RIJ_EPSILON_EBRSM && icodcl_vel[f_id] == 5) {

        cs_real_t pimp = 0.0;

        if (cs_glob_wall_functions->iwallf != 0) {
          /* Use k at I' */
          const cs_real_t xkip
            = 0.5 * (rijipb[f_id][0] + rijipb[f_id][1] + rijipb[f_id][2]);

          const cs_real_t pimp_lam =   2 * visclc * xkip
                                     / (distbf * distbf * romc);

          if (yplus > cs_math_epzero) {
            const cs_real_t pimp_turb
              = 5 * pow(uk, 4) * romc / (xkappa * visclc * (yplus + 2 * dplus));

            /* Blending between wall and homogeneous layer
               from JF Wald PhD (2016) */
            const cs_real_t fep  = exp(- pow(0.25 * (yplus + dplus), 1.5));
            const cs_real_t dep  = 1.0 - exp(- pow((yplus + dplus) / 9.0, 2.1));
            pimp = fep * pimp_lam + (1.0 - fep) * dep * pimp_turb;
          }
          else
            pimp = pimp_lam;

        }
        else {
          /* Use k at I' */
          const cs_real_t xkip
            = 0.5 * (rijipb[f_id][0] + rijipb[f_id][1] + rijipb[f_id][2]);

          pimp = 2 * visclc * xkip / (distbf * distbf * romc);
        }

        pimp = pimp * cfnne;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_ep[f_id],
                                                    &cofaf_ep[f_id],
                                                    &coefb_ep[f_id],
                                                    &cofbf_ep[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        /* If defined set Dirichlet condition for the Lagrangian time scale */

        if (f_tlag != NULL) {

          cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
          cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
          cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
          cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

          if (cs_glob_wall_functions->iwallf == 0) {
            /* No wall functions forced by user */
            pimp = 0;
          }
          else {
            /* Use of wall functions */
            if (iuntur == 1)
              pimp =   0.5 * cfnnk / (cfnne * pow(uk, 3)) * cl * xkappa
               * (  coefa_rij[f_id][0] + coefb_rij[f_id][0][0] * rijipb[f_id][0]
                  + coefa_rij[f_id][1] + coefb_rij[f_id][1][1] * rijipb[f_id][1]
                  + coefa_rij[f_id][2] + coefb_rij[f_id][2][2] * rijipb[f_id][2])
               * (dplus * visclc / (romc * uk) + rough_d);
            else
              pimp = 0.;
          }

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                      &cofaf_tlag[f_id],
                                                      &coefb_tlag[f_id],
                                                      &cofbf_tlag[f_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);

        }

        /* Alpha */

        cs_real_t *coefa_al = f_alpha->bc_coeffs->a;
        cs_real_t *coefb_al = f_alpha->bc_coeffs->b;
        cs_real_t *cofaf_al = f_alpha->bc_coeffs->af;
        cs_real_t *cofbf_al = f_alpha->bc_coeffs->bf;

        /* Dirichlet Boundary Condition
           ---------------------------- */

        if (cs_glob_wall_functions->iwallf != 0) {

          if (yplus > cs_math_epzero) {
            const cs_real_t ypsd  = 0.5 * (yplus + dplus);

            const cs_real_t falpg
              = 16. /cs_math_pow2(16 + 4.e-2 * ypsd)
              * exp(- ypsd / (16. + 4.e-2*ypsd) );

            const cs_real_t falpv
              = 1.0 - exp(-(yplus + dplus) / (16. + 4.e-2*(yplus+dplus)));

            pimp  = falpv - (yplus + dplus) * falpg;
          }
          else {
            pimp = 0.;
          }
        }
        else {
          pimp = 0.;
        }

        hint = 1.0 / distbf;
        pimp = pimp * cfnne;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_al[f_id],
                                                    &cofaf_al[f_id],
                                                    &coefb_al[f_id],
                                                    &cofbf_al[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

      }

      if (icodcl_vel[f_id] == 6 && df_limiter_eps != NULL)
        df_limiter_eps[c_id] = 0.0;

    } /* End if itytur == 3 */

    /* Boundary conditions on k, epsilon, f_bar and phi in the phi_Fbar model
       ====================================================================== */

    else if (iturb == CS_TURB_V2F_PHI) {

      cs_real_t *coefa_ep = f_eps->bc_coeffs->a;
      cs_real_t *coefb_ep = f_eps->bc_coeffs->b;
      cs_real_t *cofaf_ep = f_eps->bc_coeffs->af;
      cs_real_t *cofbf_ep = f_eps->bc_coeffs->bf;

      cs_real_t *coefa_phi = f_phi->bc_coeffs->a;
      cs_real_t *coefb_phi = f_phi->bc_coeffs->b;
      cs_real_t *cofaf_phi = f_phi->bc_coeffs->af;
      cs_real_t *cofbf_phi = f_phi->bc_coeffs->bf;

      cs_real_t *coefa_fb = f_f_bar->bc_coeffs->a;
      cs_real_t *coefb_fb = f_f_bar->bc_coeffs->b;
      cs_real_t *cofaf_fb = f_f_bar->bc_coeffs->af;
      cs_real_t *cofbf_fb = f_f_bar->bc_coeffs->bf;

      cs_real_t *coefa_k = f_k->bc_coeffs->a;
      cs_real_t *coefb_k = f_k->bc_coeffs->b;
      cs_real_t *cofaf_k = f_k->bc_coeffs->af;
      cs_real_t *cofbf_k = f_k->bc_coeffs->bf;

      /* Dirichlet Boundary Condition on k
         --------------------------------- */

      cs_real_t pimp = 0.;
      cs_real_t hint = (visclc + visctc / sigmak) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_k[f_id],
                                                  &cofaf_k[f_id],
                                                  &coefb_k[f_id],
                                                  &cofbf_k[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on epsilon
         --------------------------------------- */

      pimp = 2 * visclc / romc * cvar_k[c_id] / (distbf * distbf);
      hint = (visclc + visctc / sigmae) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_ep[f_id],
                                                  &cofaf_ep[f_id],
                                                  &coefb_ep[f_id],
                                                  &cofbf_ep[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on Lagrangian time scale
         ----------------------------------------------------- */

      if (f_tlag != NULL) {

        cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
        cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
        cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
        cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

        pimp = 0.;
        hint = (visclc + visctc / sigmak) / distbf;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                    &cofaf_tlag[f_id],
                                                    &coefb_tlag[f_id],
                                                    &cofbf_tlag[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

      }

      /* Dirichlet Boundary Condition on Phi
         ------------------------------------ */

      pimp = 0.;
      hint = (visclc + visctc / sigmak) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_phi[f_id],
                                                  &cofaf_phi[f_id],
                                                  &coefb_phi[f_id],
                                                  &cofbf_phi[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on Fb
         ---------------------------------- */

      pimp = 0.;
      hint = 1./distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_fb[f_id],
                                                  &cofaf_fb[f_id],
                                                  &coefb_fb[f_id],
                                                  &cofbf_fb[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);
    }

    /* Boundary conditions on k, epsilon, phi and alpha in the Bl-v2/k model
       ===================================================================== */

    else if (iturb == CS_TURB_V2F_BL_V2K) {

      cs_real_t *coefa_ep = f_eps->bc_coeffs->a;
      cs_real_t *coefb_ep = f_eps->bc_coeffs->b;
      cs_real_t *cofaf_ep = f_eps->bc_coeffs->af;
      cs_real_t *cofbf_ep = f_eps->bc_coeffs->bf;

      cs_real_t *coefa_phi = f_phi->bc_coeffs->a;
      cs_real_t *coefb_phi = f_phi->bc_coeffs->b;
      cs_real_t *cofaf_phi = f_phi->bc_coeffs->af;
      cs_real_t *cofbf_phi = f_phi->bc_coeffs->bf;

      cs_real_t *coefa_al = f_alpha->bc_coeffs->a;
      cs_real_t *coefb_al = f_alpha->bc_coeffs->b;
      cs_real_t *cofaf_al = f_alpha->bc_coeffs->af;
      cs_real_t *cofbf_al = f_alpha->bc_coeffs->bf;

      cs_real_t *coefa_k = f_k->bc_coeffs->a;
      cs_real_t *coefb_k = f_k->bc_coeffs->b;
      cs_real_t *cofaf_k = f_k->bc_coeffs->af;
      cs_real_t *cofbf_k = f_k->bc_coeffs->bf;

      /* Dirichlet Boundary Condition on k
         --------------------------------- */

      cs_real_t pimp = 0.;
      cs_real_t hint = (visclc + visctc / sigmak) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_k[f_id],
                                                  &cofaf_k[f_id],
                                                  &coefb_k[f_id],
                                                  &cofbf_k[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on epsilon
         --------------------------------------- */

      pimp = visclc / romc * cvar_k[c_id] / (distbf * distbf);
      hint = (visclc + visctc / sigmae) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_ep[f_id],
                                                  &cofaf_ep[f_id],
                                                  &coefb_ep[f_id],
                                                  &cofbf_ep[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on Lagrangian time scale
         ----------------------------------------------------- */

      if (f_tlag != NULL) {

        cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
        cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
        cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
        cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

        pimp = 0.;
        hint = (visclc + visctc / sigmak) / distbf;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                    &cofaf_tlag[f_id],
                                                    &coefb_tlag[f_id],
                                                    &cofbf_tlag[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

      }

      /* Dirichlet Boundary Condition on Phi
         ----------------------------------- */

      pimp = 0.;
      hint = (visclc + visctc / sigmak) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_phi[f_id],
                                                  &cofaf_phi[f_id],
                                                  &coefb_phi[f_id],
                                                  &cofbf_phi[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on alpha
         ------------------------------------- */

      pimp = 0.;
      hint = 1.0 / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_al[f_id],
                                                  &cofaf_al[f_id],
                                                  &coefb_al[f_id],
                                                  &cofbf_al[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

    }

    /* Boundary conditions on k and omega
       ================================== */

    else if (iturb == CS_TURB_K_OMEGA) {

      cs_real_t *coefa_omg = f_omg->bc_coeffs->a;
      cs_real_t *coefb_omg = f_omg->bc_coeffs->b;
      cs_real_t *cofaf_omg = f_omg->bc_coeffs->af;
      cs_real_t *cofbf_omg = f_omg->bc_coeffs->bf;

      cs_real_t *coefa_k = f_k->bc_coeffs->a;
      cs_real_t *coefb_k = f_k->bc_coeffs->b;
      cs_real_t *cofaf_k = f_k->bc_coeffs->af;
      cs_real_t *cofbf_k = f_k->bc_coeffs->bf;

      /* Dirichlet Boundary Condition on k
         --------------------------------- */

      /* pimp > 0 if we are outside the visous sub-layer (really or through
         the scalable wall functions).
         pimp = 0 if we are in the viscous sub-layer */
      cs_real_t pimp = (iuntur == 1 || icodcl_vel[f_id] == 6) ?
                        uk*uk/sqrcmu : 0.0;

      /* FIXME it is wrong because sigma is computed within the model
         see cs_turbulence_kw.c */
      cs_real_t hint = (visclc + visctc / cs_turb_ckwsk2) / distbf;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_k[f_id],
                                                  &cofaf_k[f_id],
                                                  &coefb_k[f_id],
                                                  &cofbf_k[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition on omega
         ------------------------------------- */

      /* FIXME it is wrong because sigma is computed within the model
         see cs_turbulence_kw.c (so the flux is not the one we impose!) */
      hint = (visclc + visctc / cs_turb_ckwsw2) / distbf;

      if (cs_glob_turb_rans_model->ikwcln == 1 && icodcl_vel[f_id] == 5) {
        /* In viscous sub layer */
        const cs_real_t pimp_lam
          = 60 * visclc / (romc * cs_turb_ckwbt1 * distbf * distbf);

        /* If we are outside the viscous sub-layer (either naturally, or
           artificialy using scalable wall functions) */

        if (yplus > cs_math_epzero) {
          const cs_real_t pimp_turb
            = 5*uk*uk*romc / (sqrcmu*xkappa*visclc*(yplus+2*dplus));

          /* Use gamma function of Kader to weight
             between high and low Reynolds meshes */

          const cs_real_t gammap
            = - 0.01 * pow(yplus+2*dplus, 4) / (1.0 + 5.0 * (yplus+2*dplus));

          pimp = pimp_lam * exp(gammap) + exp(1.0/gammap) * pimp_turb;
        }
        else
          pimp = pimp_lam;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_omg[f_id],
                                                    &cofaf_omg[f_id],
                                                    &coefb_omg[f_id],
                                                    &cofbf_omg[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

      }
      /* If ikwcln is equal to 0, switch to deprecated Neumann
         condition on omega. */
      else {

        /* In viscous sub layer */
        const cs_real_t pimp_lam
          = 120 * 8 * visclc / (romc * cs_turb_ckwbt1 *distbf * distbf);

        /* If we are outside the viscous sub-layer (either naturally, or
           artificialy using scalable wall functions) */

        if (yplus > cs_math_epzero) {

          const cs_real_t pimp_turb = distbf*4*cs_math_pow3(uk)*romc*romc /
            (sqrcmu*xkappa*visclc*visclc*pow(yplus+2*dplus, 2));

          /* Use gamma function of Kader to weight
             between high and low Reynolds meshes */

          const cs_real_t gammap
            = -0.01 * pow(yplus+2*dplus, 4) / (1 + 5 *(yplus+2*dplus));

          pimp = pimp_lam * exp(gammap) + exp(1.0/gammap) * pimp_turb;
        }
        else
          pimp = pimp_lam;

        /* Compute pimp for rough wall */
        if (icodcl_vel[f_id] == 6)
          pimp = distbf*uk/(sqrcmu*xkappa*ydep*ydep) * cfnne / cfnnk;

        /* TODO transform it, it is only to be fully equivalent */
        const cs_real_t qimp = - pimp * hint;

        cs_boundary_conditions_set_neumann_scalar(&coefa_omg[f_id],
                                                  &cofaf_omg[f_id],
                                                  &coefb_omg[f_id],
                                                  &cofbf_omg[f_id],
                                                  qimp,
                                                  hint);

      }

      /* If defined set Dirichlet condition for the Lagrangian time scale */

      if (f_tlag != NULL) {

        cs_real_t *coefa_tlag = f_tlag->bc_coeffs->a;
        cs_real_t *coefb_tlag = f_tlag->bc_coeffs->b;
        cs_real_t *cofaf_tlag = f_tlag->bc_coeffs->af;
        cs_real_t *cofbf_tlag = f_tlag->bc_coeffs->af;

        if (cs_glob_wall_functions->iwallf == 0)
          /* No wall functions forced by user */
          pimp = 0.;
        else {
          /* Use of wall functions */
          if (iuntur == 1) {

            if (icodcl_vel[f_id] == 5)
              pimp = cfnnk / (cfnne * uk) * cl / sqrcmu * xkappa
                * (dplus * visclc / (romc * uk) + rough_d);
            else if (icodcl_vel[f_id] == 6)
              pimp = cfnnk / (cfnne * uk) * cl / sqrcmu * xkappa * rough_d;

          }
          else
            pimp = 0.;
        }

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_tlag[f_id],
                                                    &cofaf_tlag[f_id],
                                                    &coefb_tlag[f_id],
                                                    &cofbf_tlag[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);
      }
    }

    /* Boundary conditions on the Spalart Allmaras turbulence model
       ============================================================ */

    else if (iturb == CS_TURB_SPALART_ALLMARAS) {

      /* Dirichlet Boundary Condition on nusa
         ------------------------------------ */

      cs_real_t *coefa_nusa = f_nusa->bc_coeffs->a;
      cs_real_t *coefb_nusa = f_nusa->bc_coeffs->b;
      cs_real_t *cofaf_nusa = f_nusa->bc_coeffs->af;
      cs_real_t *cofbf_nusa = f_nusa->bc_coeffs->bf;

      if (icodcl_vel[f_id] == 5) {
        cs_real_t pimp = 0.;
        /* Note: nusa is zero at the wall */
        cs_real_t hint = visclc / distbf / cs_turb_csasig;

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_nusa[f_id],
                                                    &cofaf_nusa[f_id],
                                                    &coefb_nusa[f_id],
                                                    &cofbf_nusa[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);
      }
      else if (icodcl_vel[f_id] == 6) {

        const cs_real_t *cvara_nusa = f_nusa->val_pre;

        /* FIXME is it the sand grain roughness or the length scale as here? */
        const cs_real_t dsa0 = rough_d;
        const cs_real_t hint
          = (visclc + eqp_nusa->idifft * cvara_nusa[c_id] * romc * dsa0
             / (distbf + dsa0) ) / distbf / cs_turb_csasig;

        /* If we have a rough wall then:
           nusa_wall*(1- I'F/d0)=nusa_I'
           which is a Robin type BC */

        coefa_nusa[f_id] = 0.0;
        coefb_nusa[f_id] = dsa0 / (dsa0 + distbf);

        cofaf_nusa[f_id] = 0.0;
        cofbf_nusa[f_id] = hint * distbf / (dsa0 + distbf);

      }

    }

    byplus[f_id] = yplus;
    bdplus[f_id] = dplus;
    bpro_ustar[f_id] = uet;

    /* FIXME note taken into account yet in hturbp cfnns */
    bcfnns[f_id] = (icodcl_vel[f_id] == 5) ? 1.0 : cfnns;
    bdlmo[f_id] = dlmo;
    bpro_uk[f_id] = uk;

  } /* End of loop over faces */

  /* Boundary conditions on the other scalars
     (Specific treatment for the variances of the scalars next to walls:
     see cs_boundary_condition_set_coeffs)
     =================================================================== */

  for (int fld_id = 0; fld_id < n_fields; fld_id++) {

    cs_field_t *f_scal = cs_field_by_id(fld_id);

    if (!(f_scal->type & CS_FIELD_VARIABLE))
      continue;
    if (cs_field_get_key_int(f_scal, keysca) <= 0)
      continue;

    const int iscavr = cs_field_get_key_int(f_scal, kscavr);

    if (iscavr <= 0) {

      if (f_scal->dim == 1) {
        _cs_boundary_conditions_set_coeffs_turb_scalar
          (f_scal, isvhb, byplus, bdplus,
           bpro_uk, bpro_ustar, bcfnns, bdlmo, hbord,
          theipb, &tetmax, &tetmin, &tplumx, &tplumn);

      }

      /* Vector field */
      else {
        _cs_boundary_conditions_set_coeffs_turb_vector(f_scal,
                                                       byplus,
                                                       bdplus,
                                                       bpro_uk);
      }
    }

  }

  cs_gnum_t n_per_layer[3] = {nlogla, nsubla, iuiptn};
  cs_parall_counter(n_per_layer, 3);

  if (cs_glob_rank_id > -1) {

    int n_minmax = (f_th != NULL) ? 7 : 4;

    cs_real_t umin[7]
      = {uiptmn, uetmin, ukmin, yplumn, tetmin, tplumn, dlmomin};
    cs_parall_min(n_minmax, CS_DOUBLE, umin);

    uiptmn = umin[0];
    uetmin = umin[1];
    ukmin  = umin[2];
    yplumn = umin[3];

    cs_real_t umax[7]
      = {uiptmx, uetmax, ukmax, yplumx, tetmax, tplumx, dlmomax};
    cs_parall_max(n_minmax, CS_DOUBLE, umax);

    uiptmx = umax[0];
    uetmax = umax[1];
    ukmax  = umax[2];
    yplumx = umax[3];

    if (f_th != NULL) {
      tetmin  = umin[4];
      tplumn  = umin[5];
      dlmomin = umin[6];

      tetmax  = umax[4];
      tplumx  = umax[5];
      dlmomax = umax[6];
    }
  }

  BFT_FREE(byplus);
  BFT_FREE(_buk);
  BFT_FREE(buet);
  BFT_FREE(bcfnns_loc);
  BFT_FREE(bdplus);
  BFT_FREE(bdlmo);

  /* Logging
     ======= */

  /* Remark: so as not to encumber logs when only a few y+ values are not
     correct, the message is produced only at the 2 first and last time steps,
     or if the verbosity is  >= 2.
     We also indicate the number of the last time step at which y+ outside
     admissible bounds was encountered. */

  cs_real_t ypluli = cs_glob_wall_functions->ypluli;

  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(vel);

  if (eqp_vel->verbosity >= 0) {

    bool log_is_active = cs_log_default_is_active();
    if (eqp_vel->verbosity >= 2)
      log_is_active = true;

    bool warn_refine = false;

    if (   (iturb == CS_TURB_NONE && n_per_layer[0] != 0)
        ||  (itytur == 5 && n_per_layer[0] != 0)
        || ((itytur == 2 || itytur == 3) && n_per_layer[1] > 0))
      _ntlast = nt_cur;

    if (  (_ntlast == nt_cur && _iaff < 2)
        || (_ntlast >= 0 && nt_cur >= nt_max-1)
        || (_ntlast == nt_cur && eqp_vel->verbosity >= 2)) {
      _iaff = _iaff + 1;
      warn_refine = true;
    }

    if (log_is_active || warn_refine) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "   ** Boundary conditions for walls\n"
           "      -----------------------------\n\n"));
      cs_log_separator(CS_LOG_DEFAULT);
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("                                         Minimum     Maximum\n"));
      cs_log_separator(CS_LOG_DEFAULT);
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("   Rel velocity at the wall uiptn : %12.5e %12.5e\n"
           "   Friction velocity        uet   : %12.5e %12.5e\n"
           "   Friction velocity        uk    : %12.5e %12.5e\n"
           "   Dimensionless distance   yplus : %12.5e %12.5e\n"),
         uiptmn, uiptmx, uetmin, uetmax, ukmin, ukmax, yplumn, yplumx);

      if (f_th != NULL) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("   Friction thermal sca.    tstar : %12.5e %12.5e\n"
             "   Dim-less thermal sca.    tplus : %12.5e %12.5e\n"),
           tetmin, tetmax, tplumn, tplumx);
        if (iwalfs == CS_WALL_F_S_MONIN_OBUKHOV)
          cs_log_printf
            (CS_LOG_DEFAULT,
             _("   Inverse Monin-Ob. length dlmo  : %12.5e %12.5e\n"),
             dlmomin, dlmomax);
      }

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("   ------------------------------------------------------\n"
           "   Nb of reversals of the velocity at the wall: %llu\n"
           "   Nb of faces within the viscous sub-layer   : %llu\n"
           "   Total number of wall faces                 : %llu\n"
           "------------------------------------------------------------\n"),
         (unsigned long long)n_per_layer[2],
         (unsigned long long)n_per_layer[1],
         (unsigned long long)(n_per_layer[1] + n_per_layer[0]));

    }

    if (warn_refine) {

      bool need_close = false;

      if (iturb == CS_TURB_NONE) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ Warning: mesh too coarse at the wall\n"
             "@    ========\n"
             "@    The mesh does not seem to be refined enough at the wall\n"
             "@      to be able to run a laminar simulation.\n"
             "@\n"
             "@    The last time step at which too large values for the\n"
             "@      dimensionless distance to the wall (yplus) have been\n"
             "@      observed is the time step %d\n"
             "@\n"
             "@    The minimum value for yplus must be lower than the\n"
             "@      limit value 'ypluli' = %14.5e\n"
             "@\n"
             "@    Visualize the distribution of yplus at the wall\n"
             "@      (with ParaView for example) to conclude on\n"
             "@      the way the results quality might be affected.\n"),
           _ntlast, ypluli);
        need_close = true;
      }

      if (itytur == 5) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ Warning: mesh too coarse at the wall\n"
             "@    ========\n"
             "@    The mesh does not seem to be refined enough at the wall\n"
             "@      to be able to run a v2f simulation\n"
             "@      (phi-fbar or BL-v2/k)\n"
             "@\n"
             "@    The last time step at which too large values for the\n"
             "@      dimensionless distance to the wall (yplus) have been\n"
             "@      observed is the time step %10d\n"
             "@\n"
             "@    The minimum value for yplus must be lower than the\n"
             "@      limit value 'ypluli' = %14.5e\n"
             "@\n"
             "@    Visualize the distribution of yplus at the wall\n"
             "@      (with ParaView for example) to conclude on\n"
             "@      the way the results quality might be affected.\n"),
           _ntlast, ypluli);
        need_close = true;
      }

      /* No warnings in EBRSM */
      if (   (itytur ==  2 && iturb != CS_TURB_K_EPSILON_LS)
          || iturb == CS_TURB_RIJ_EPSILON_LRR
          || iturb == CS_TURB_RIJ_EPSILON_SSG) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ Warning: mesh too fine at the wall\n"
             "@    ========\n"
             "@    The mesh seems to be too fine at the wall to use\n"
             "@      a high-Reynolds turbulence model.\n"
             "@\n"
             "@    The last time step at which too small values for the\n"
             "@      dimensionless distance to the wall (yplus) have been\n"
             "@      observed is the time step %10d\n"
             "@\n"
             "@    The minimum value for yplus must be greater than the\n"
             "@      limit value 'ypluli' = %14.5e\n"
             "@\n"
             "@    Visualize the distribution of yplus at the wall\n"
             "@      (with ParaView for example) to conclude on\n"
             "@      the way the results quality might be affected.\n"),
           _ntlast, ypluli);
        need_close = true;
      }

      if (   eqp_vel->verbosity < 2
          && (   iturb !=  CS_TURB_RIJ_EPSILON_EBRSM
              && iturb !=  CS_TURB_K_EPSILON_LS)) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@\n"
             "@    This warning is only printed at the first two\n"
             "@    occurences of the problem and at the last time step\n"
             "@    of the calculation. The vanishing of the message does\n"
             "@    not necessarily mean the vanishing of the problem.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"));
        need_close = false;
      }

      if (need_close) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"));
      }

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
