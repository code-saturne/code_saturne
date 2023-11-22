/*============================================================================
 * Symmetry boundary condition management.
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

#include "cs_ale.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_set_coeffs_symmetry.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Boundary conditions for symmetry (icodcl = 4) for a scalar.
 *
 * \param[in]  f_sc  scalar field
 */
/*---------------------------------------------------------------------------*/

static void
_boundary_conditions_set_coeffs_symmetry_scalar(cs_field_t  *f_sc)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)fvq->b_face_u_normal;

  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const cs_real_t cp0 = fluid_props->cp0;
  const int icp = fluid_props->icp;

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kturt  = cs_field_key_id("turbulent_flux_model");
  const int kscacp = cs_field_key_id("is_temperature");

  const cs_real_t *viscl = CS_F_(mu)->val;
  const int ifcvsl = cs_field_get_key_int(f_sc, kivisl);

  /* Get the turbulent flux model for the scalar */

  const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
  const cs_real_t ctheta = cs_field_get_key_double(f_sc, kctheta);
  const int turb_flux_model = cs_field_get_key_int(f_sc, kturt);

  cs_field_t *f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
  cs_real_6_t *visten = (cs_real_6_t *)f_a_t_visc->val;
  cs_real_t *crom = crom = CS_F_(rho)->val;

  const cs_real_t *cpro_cp = NULL, *cpro_cv = NULL;
  if (icp >= 0)
    cpro_cp = CS_F_(cp)->val;

  cs_field_t *f_id_cv = cs_field_by_name_try("isobaric_heat_capacity");
  if (f_id_cv != NULL)
    cpro_cv = f_id_cv->val;

  /* Reference diffusivity of the primal scalar */
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  cs_real_t visls_0 = cs_field_get_key_double(f_sc, kvisl0);

  const cs_real_t *viscls = NULL;
  if (ifcvsl >= 0)
    viscls = cs_field_by_id(ifcvsl)->val;

  /* Does the scalar behave as a temperature ? */
  const int iscacp = cs_field_get_key_int(f_sc, kscacp);

  /* Turbulent diffusive flux of the scalar T
     (blending factor so that the component v'T' have only
     mu_T/(mu+mu_T)* Phi_T) */

  /* Turbulent flux */
  cs_field_t *f_tf
    = cs_field_by_composite_name(f_sc->name, "turbulent_flux");

  cs_real_3_t  *coefa_tf = (cs_real_3_t  *)f_tf->bc_coeffs->a;
  cs_real_33_t *coefb_tf = (cs_real_33_t *)f_tf->bc_coeffs->b;
  cs_real_3_t  *cofaf_tf = (cs_real_3_t  *)f_tf->bc_coeffs->af;
  cs_real_33_t *cofbf_tf = (cs_real_33_t *)f_tf->bc_coeffs->bf;
  cs_real_3_t  *cofar_tf = (cs_real_3_t  *)f_tf->bc_coeffs->ad;
  cs_real_33_t *cofbr_tf = (cs_real_33_t *)f_tf->bc_coeffs->bd;

  const int *icodcl_vel = CS_F_(vel)->bc_coeffs->icodcl;

  /* Loop on boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test on symmetry boundary condition: start */
    if (icodcl_vel[f_id] == 4) {

      /* Geometric quantities */
      const cs_lnum_t c_id = b_face_cells[f_id];
      const cs_real_t *nn = b_face_u_normal[f_id];
      const cs_real_t distbf = b_dist[f_id];

      /* Physical Properties */
      const cs_real_t visclc = viscl[c_id];

      cs_real_t cpp = 1.;
      if (iscacp == 1)
        cpp = (icp >= 0) ? cpro_cp[c_id] : cp0;
      else if (iscacp == 2)
        cpp = cpro_cv[c_id];

      const cs_real_t rkl = (ifcvsl < 0) ? visls_0/cpp : viscls[c_id]/cpp;

      /* FIXME for EB DFM */
      cs_real_t hintt[6] = {0., 0., 0., 0., 0., 0.};
      for (int isou = 0; isou < 6; isou++) {

        if (isou <= 3)
          hintt[isou] =   (0.5 * (visclc + rkl)
                        + ctheta*visten[c_id][isou]/cs_turb_csrij) / distbf;
        else
          hintt[isou] = ctheta * visten[c_id][isou] / cs_turb_csrij / distbf;

      }

      /* Gradient BCs */
      coefa_tf[f_id][0] = 0.0;
      coefa_tf[f_id][1] = 0.0;
      coefa_tf[f_id][2] = 0.0;

      coefb_tf[f_id][0][0] = 1.0 - nn[0]*nn[0];
      coefb_tf[f_id][1][1] = 1.0 - nn[1]*nn[1];
      coefb_tf[f_id][2][2] = 1.0 - nn[2]*nn[2];

      coefb_tf[f_id][0][1] = - nn[0] * nn[1];
      coefb_tf[f_id][0][2] = - nn[0] * nn[2];
      coefb_tf[f_id][1][0] = - nn[1] * nn[0];
      coefb_tf[f_id][1][2] = - nn[1] * nn[2];
      coefb_tf[f_id][2][0] = - nn[2] * nn[0];
      coefb_tf[f_id][2][1] = - nn[2] * nn[1];

      /* Flux BCs */
      cofaf_tf[f_id][0] = 0.0;
      cofaf_tf[f_id][1] = 0.0;
      cofaf_tf[f_id][2] = 0.0;

      /* Diagonal */

      cofbf_tf[f_id][0][0] =   hintt[0]*nn[0]*nn[0]
                             + hintt[3]*nn[0]*nn[1]
                             + hintt[5]*nn[0]*nn[2];

      cofbf_tf[f_id][1][1] =   hintt[3]*nn[0]*nn[1]
                             + hintt[1]*nn[1]*nn[1]
                             + hintt[4]*nn[1]*nn[2];

      cofbf_tf[f_id][2][2] =   hintt[5]*nn[0]*nn[2]
                             + hintt[4]*nn[1]*nn[2]
                             + hintt[2]*nn[2]*nn[2];

      /* Extra diagonal */

      cofbf_tf[f_id][1][0] =   hintt[0]*nn[0]*nn[1]
                             + hintt[3]*nn[1]*nn[1]
                             + hintt[5]*nn[1]*nn[2];

      cofbf_tf[f_id][0][1] =   hintt[3]*nn[0]*nn[0]
                             + hintt[1]*nn[1]*nn[0]
                             + hintt[4]*nn[0]*nn[2];

      cofbf_tf[f_id][2][0] =   hintt[0]*nn[0]*nn[2]
                             + hintt[3]*nn[1]*nn[2]
                             + hintt[5]*nn[2]*nn[2];

      cofbf_tf[f_id][0][2] =   hintt[5]*nn[0]*nn[0]
                             + hintt[4]*nn[1]*nn[0]
                             + hintt[2]*nn[2]*nn[0];

      cofbf_tf[f_id][2][1] =   hintt[3]*nn[0]*nn[2]
                             + hintt[1]*nn[1]*nn[2]
                             + hintt[4]*nn[2]*nn[2];

      cofbf_tf[f_id][1][2] =   hintt[5]*nn[0]*nn[1]
                             + hintt[4]*nn[1]*nn[1]
                             + hintt[2]*nn[2]*nn[1];

      /* Boundary conditions for thermal transport equation */
      for (int isou = 0; isou < 3; isou++) {
        cofar_tf[f_id][isou] = coefa_tf[f_id][isou];
        for (int jsou = 0; jsou < 3; jsou++)
          cofbr_tf[f_id][isou][jsou] = coefb_tf[f_id][isou][jsou];
      }

      /* EB-GGDH/AFM/DFM alpha boundary conditions */
      if (turb_flux_model == 11 || turb_flux_model == 21 || turb_flux_model == 31) {

        cs_field_t *f_al = cs_field_by_composite_name(f_sc->name, "alpha");
        cs_real_t *coefa_al = f_al->bc_coeffs->a;
        cs_real_t *coefb_al = f_al->bc_coeffs->b;
        cs_real_t *cofaf_al = f_al->bc_coeffs->af;
        cs_real_t *cofbf_al = f_al->bc_coeffs->bf;

        /* Dirichlet Boundary Condition
           ---------------------------- */

        const cs_real_t qimp = 0.0;

        const cs_real_t hint = 1.0 / distbf;

        cs_boundary_conditions_set_neumann_scalar(&coefa_al[f_id],
                                                  &cofaf_al[f_id],
                                                  &coefb_al[f_id],
                                                  &cofbf_al[f_id],
                                                  qimp,
                                                  hint);


      }

    } /* Test on velocity symmetry condition: end */

  } /* End loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*! \brief Boundary conditions for symmetry (icodcl = 4) for a vector.
 *
 * \param[in]  f_v  vector field
 */
/*----------------------------------------------------------------------------*/

static void
_boundary_conditions_set_coeffs_symmetry_vector(cs_field_t  *f_v)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_turb_model_type_t iturb = cs_glob_turb_model->iturb;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)fvq->b_face_u_normal;

  const int kivisl  = cs_field_key_id("diffusivity_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const int ifcvsl  = cs_field_get_key_int(f_v, kivisl);

  cs_equation_param_t *eqp_v = cs_field_get_equation_param(f_v);
  cs_real_3_t  *coefa_v = (cs_real_3_t  *)f_v->bc_coeffs->a;
  cs_real_33_t *coefb_v = (cs_real_33_t *)f_v->bc_coeffs->b;
  cs_real_3_t  *cofaf_v = (cs_real_3_t  *)f_v->bc_coeffs->af;
  cs_real_33_t *cofbf_v = (cs_real_33_t *)f_v->bc_coeffs->bf;

  cs_real_6_t *visten = NULL;
  if (eqp_v->idften & CS_ANISOTROPIC_DIFFUSION) {
    if (iturb != CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_field_t *f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
      visten = (cs_real_6_t *)f_a_t_visc->val;
    }
    else {/* EBRSM and (GGDH or AFM) */
      cs_field_t *f_vis = cs_field_by_name("anisotropic_turbulent_viscosity_scalar");
      visten = (cs_real_6_t *)f_vis->val;
    }
  }

  const cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *crom = crom = CS_F_(rho)->val;

  const cs_real_t *viscls = NULL;
  if (ifcvsl >= 0)
    viscls = cs_field_by_id(ifcvsl)->val;

  const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
  const cs_real_t ctheta = cs_field_get_key_double(f_v, kctheta);

  /* Retrieve turbulent Schmidt value for current scalar */
  const cs_real_t turb_schmidt = cs_field_get_key_double(f_v, ksigmas);

  /* Reference diffusivity */
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  cs_real_t visls_0 = cs_field_get_key_double(f_v, kvisl0);

  const int *icodcl_v = f_v->bc_coeffs->icodcl;

  /* Loop on boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test on symmetry boundary condition: start */
    if (icodcl_v[f_id] == 4) {

      /* Geometric quantities */
      const cs_lnum_t c_id = b_face_cells[f_id];
      const cs_real_t *nn = b_face_u_normal[f_id];
      const cs_real_t distbf = b_dist[f_id];

      const cs_real_t rkl = (ifcvsl < 0) ? visls_0 : viscls[c_id];

      /* Isotropic diffusivity */

      cs_real_t hintt[6] = {0., 0., 0., 0., 0., 0.};

      if (eqp_v->idften & CS_ISOTROPIC_DIFFUSION) { // a voir (zero)
        hintt[0] =   (eqp_v->idifft * cs_math_fmax(visct[c_id], 0)
                   / turb_schmidt + rkl) / distbf;

        hintt[1] = hintt[0];
        hintt[2] = hintt[0];
        hintt[3] = 0.0;
        hintt[4] = 0.0;
        hintt[5] = 0.0;
      }

      /* Symmetric tensor diffusivity */
      else if (eqp_v->idften & CS_ANISOTROPIC_DIFFUSION) {
        const cs_real_t temp = eqp_v->idifft * ctheta / cs_turb_csrij;

        hintt[0] = (temp * visten[c_id][0] + rkl) / distbf;
        hintt[1] = (temp * visten[c_id][1] + rkl) / distbf;
        hintt[2] = (temp * visten[c_id][2] + rkl) / distbf;
        hintt[3] =  temp * visten[c_id][3]        / distbf;
        hintt[4] =  temp * visten[c_id][4]        / distbf;
        hintt[5] =  temp * visten[c_id][5]        / distbf;
      }

      /* Gradient BCs */
      coefa_v[f_id][0] = 0.0;
      coefa_v[f_id][1] = 0.0;
      coefa_v[f_id][2] = 0.0;

      coefb_v[f_id][0][0] = 1.0 - nn[0] * nn[0];
      coefb_v[f_id][1][1] = 1.0 - nn[1] * nn[1];
      coefb_v[f_id][2][2] = 1.0 - nn[2] * nn[2];

      coefb_v[f_id][0][1] = - nn[0] * nn[1];
      coefb_v[f_id][0][2] = - nn[0] * nn[2];
      coefb_v[f_id][1][0] = - nn[1] * nn[0];
      coefb_v[f_id][1][2] = - nn[1] * nn[2];
      coefb_v[f_id][2][0] = - nn[2] * nn[0];
      coefb_v[f_id][2][1] = - nn[2] * nn[1];

      /* Flux BCs */
      cofaf_v[f_id][0] = 0.0;
      cofaf_v[f_id][1] = 0.0;
      cofaf_v[f_id][2] = 0.0;

      /* Diagonal */

      cofbf_v[f_id][0][0] =   hintt[0]*nn[0]*nn[0]
                            + hintt[3]*nn[0]*nn[1]
                            + hintt[5]*nn[0]*nn[2];

      cofbf_v[f_id][1][1] =   hintt[3]*nn[0]*nn[1]
                            + hintt[1]*nn[1]*nn[1]
                            + hintt[4]*nn[1]*nn[2];

      cofbf_v[f_id][2][2] =   hintt[5]*nn[0]*nn[2]
                            + hintt[4]*nn[1]*nn[2]
                            + hintt[2]*nn[2]*nn[2];

      /* Extra diagonal */

      cofbf_v[f_id][1][0] =   hintt[0]*nn[0]*nn[1]
                            + hintt[3]*nn[1]*nn[1]
                            + hintt[5]*nn[1]*nn[2];

      cofbf_v[f_id][0][1] =   hintt[0]*nn[0]*nn[1]
                            + hintt[3]*nn[1]*nn[1]
                            + hintt[5]*nn[1]*nn[2];

      cofbf_v[f_id][2][0] =   hintt[0]*nn[0]*nn[2]
                            + hintt[3]*nn[1]*nn[2]
                            + hintt[5]*nn[2]*nn[2];

      cofbf_v[f_id][0][2] =   hintt[0]*nn[0]*nn[2]
                            + hintt[3]*nn[1]*nn[2]
                            + hintt[5]*nn[2]*nn[2];

      cofbf_v[f_id][2][1] =   hintt[3]*nn[0]*nn[2]
                            + hintt[1]*nn[1]*nn[2]
                            + hintt[4]*nn[2]*nn[2];

      cofbf_v[f_id][1][2] =   hintt[3]*nn[0]*nn[2]
                            + hintt[1]*nn[1]*nn[2]
                            + hintt[4]*nn[2]*nn[2];

    } /* Test on velocity symmetry condition: end */

  } /* End loop on boundary faces */
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Symmetry boundary conditions for vectors and tensors.
 *
 * Corresponds to the code icodcl(ivar) = 4.
 *
 * Please refer to the
 * <a href="../../theory.pdf#clsyvt"><b>clsyvt</b></a> section of the
 * theory guide for more informations.
 *
 * \param[in]  velipb  value of the velocity at \f$ \centip \f$
 *                     of boundary cells
 * \param[in]  rijipb  value of \f$ R_{ij} \f$ at \f$ \centip \f$
 *                     of boundary cells
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs_symmetry(cs_real_t  velipb[][3],
                                           cs_real_t  rijipb[][6])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)mesh->b_face_cells;
  const cs_real_3_t *b_face_u_normal = (const cs_real_3_t *)fvq->b_face_u_normal;
  const cs_real_t   *b_dist         = fvq->b_dist;
  int               *isympa         = fvq->b_sym_flag;

  const int keysca  = cs_field_key_id("scalar_id");
  const int kturt   = cs_field_key_id("turbulent_flux_model");

  const int itytur = cs_glob_turb_model->itytur;
  const int idirsm = cs_glob_turb_rans_model->idirsm;
  const int n_fields = cs_field_n_fields();

  /* Initializations
     =============== */

  cs_real_6_t *visten = NULL;
  if (itytur == 3 && idirsm == 1) {
    cs_field_t *f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
    visten = (cs_real_6_t *)f_a_t_visc->val;
  }

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  /* Boundary Conditions */

  cs_field_t *vel = CS_F_(vel);
  cs_real_3_t  *coefa_vel = (cs_real_3_t  *)vel->bc_coeffs->a;
  cs_real_33_t *coefb_vel = (cs_real_33_t *)vel->bc_coeffs->b;
  cs_real_3_t  *cofaf_vel = (cs_real_3_t  *)vel->bc_coeffs->af;
  cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

  const int *icodcl_vel = vel->bc_coeffs->icodcl;
  cs_real_t *rcodcl1_vel = vel->bc_coeffs->rcodcl1;

  /* Loop over boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* Test for symmetry on velocity */
    if (icodcl_vel[f_id] != 4)
      continue;

    /* To cancel the mass flux */
    isympa[f_id] = 0;

    /* Geometric quantities */
    const cs_real_t *nn = b_face_u_normal[f_id];

    /* Local reference frame
       --------------------- */

    /* En ALE, we may have a displacement velocity of the face for which only
       the normal component counts (we continue to determine TX based on the
       absolute tangential velocity because the orientation of TX et T2X is
       not important for symmetry). */

    cs_real_t rcodcn = 0.0;
    if (cs_glob_ale > CS_ALE_NONE) {

      const cs_real_t rcodclxyz[3] = {rcodcl1_vel[n_b_faces*0 + f_id],
                                      rcodcl1_vel[n_b_faces*1 + f_id],
                                      rcodcl1_vel[n_b_faces*2 + f_id]};

      rcodcn = cs_math_3_dot_product(rcodclxyz, nn);

    }

    const cs_real_t upxyz[3] = {velipb[f_id][0],
                                velipb[f_id][1],
                                velipb[f_id][2]};

    cs_real_t eloglo[3][3], alpha[6][6];

    if (itytur == 3) {

      /* Relative tangential velocity */

      const cs_real_t usn = cs_math_3_dot_product(upxyz, nn);

      cs_real_t txyz[3] = {upxyz[0] - usn * nn[0],
                           upxyz[1] - usn * nn[1],
                           upxyz[2] - usn * nn[2]};

      /* Unit tangent
         If the velocity is zero,
         Tx, Ty, Tz is not used (we cancel the velocity), so we assign any
         value (zero for example) */
      cs_math_3_normalize(txyz, txyz);

      /* --> T2 = RN X T (where X is the cross product) */

      const cs_real_t t2xyz[3] = {nn[1]*txyz[2] - nn[2]*txyz[1],
                                  nn[2]*txyz[0] - nn[0]*txyz[2],
                                  nn[0]*txyz[1] - nn[1]*txyz[0]};

      /* Orthogonal matrix for change of reference frame ELOGLOij
         (from local to global reference frame)

                  | TX    TY    TZ |
         ELOGLO = |-RNX  -RNY  -RNZ|
                  | T2X   T2Y   T2Z|

         Its transpose ELOGLOt is its inverse. */

      eloglo[0][0] =  txyz[0];
      eloglo[1][0] = -nn[0];
      eloglo[2][0] =  t2xyz[0];

      eloglo[0][1] =  txyz[1];
      eloglo[1][1] = -nn[1];
      eloglo[2][1] =  t2xyz[1];

      eloglo[0][2] =  txyz[2];
      eloglo[1][2] = -nn[2];
      eloglo[2][2] =  t2xyz[2];

      /* Compute Reynolds stress transformation matrix */

      int clsyme = 1;
      cs_turbulence_bc_rij_transform(clsyme, eloglo, alpha);

    }

    /* Boundary conditions on the velocity
       (totaly implicit)
       The condition is a zero (except in ALE) Dirichlet on the normal component
       a homogenous Neumann on the other components
       ------------------------------------------------------------------------- */

    /* Coupled solving of the velocity components */

    const cs_lnum_t c_id = b_face_cells[f_id];

    /* Physical properties */
    const cs_real_t visclc = viscl[c_id];
    const cs_real_t visctc = visct[c_id];

    /* Geometrical quantity */
    const cs_real_t distbf = b_dist[f_id];

    const cs_real_t hint = (itytur == 3) ? visclc / distbf
                                         : (visclc + visctc) / distbf;

    /* Gradient BCs */
    coefa_vel[f_id][0] = rcodcn * nn[0];
    coefa_vel[f_id][1] = rcodcn * nn[1];
    coefa_vel[f_id][2] = rcodcn * nn[2];

    coefb_vel[f_id][0][0] = 1.0 - nn[0] * nn[0];
    coefb_vel[f_id][1][1] = 1.0 - nn[1] * nn[1];
    coefb_vel[f_id][2][2] = 1.0 - nn[2] * nn[2];

    coefb_vel[f_id][0][1] = - nn[0] * nn[1];
    coefb_vel[f_id][0][2] = - nn[0] * nn[2];
    coefb_vel[f_id][1][0] = - nn[1] * nn[0];
    coefb_vel[f_id][1][2] = - nn[1] * nn[2];
    coefb_vel[f_id][2][0] = - nn[2] * nn[0];;
    coefb_vel[f_id][2][1] = - nn[2] * nn[1];

    /* Flux BCs */
    cofaf_vel[f_id][0] = - hint * rcodcn * nn[0];
    cofaf_vel[f_id][1] = - hint * rcodcn * nn[1];
    cofaf_vel[f_id][2] = - hint * rcodcn * nn[2];

    cofbf_vel[f_id][0][0] = hint * nn[0] * nn[0];
    cofbf_vel[f_id][1][1] = hint * nn[1] * nn[1];
    cofbf_vel[f_id][2][2] = hint * nn[2] * nn[2];

    cofbf_vel[f_id][0][1] = hint * nn[0] * nn[1];
    cofbf_vel[f_id][0][2] = hint * nn[0] * nn[2];
    cofbf_vel[f_id][1][0] = hint * nn[1] * nn[0];
    cofbf_vel[f_id][1][2] = hint * nn[1] * nn[2];
    cofbf_vel[f_id][2][0] = hint * nn[2] * nn[0];
    cofbf_vel[f_id][2][1] = hint * nn[2] * nn[1];

    /* Boundary conditions on Rij (partially implicited)
       ================================================= */

    if (itytur == 3) {

      cs_field_t *rij = CS_F_(rij);
      cs_equation_param_t *eqp_rij = cs_field_get_equation_param(rij);

      cs_real_6_t  *coefa_rij = (cs_real_6_t  *)rij->bc_coeffs->a;
      cs_real_66_t *coefb_rij = (cs_real_66_t *)rij->bc_coeffs->b;
      cs_real_6_t  *cofaf_rij = (cs_real_6_t  *)rij->bc_coeffs->af;
      cs_real_66_t *cofbf_rij = (cs_real_66_t *)rij->bc_coeffs->bf;
      cs_real_6_t  *cofad_rij = (cs_real_6_t  *)rij->bc_coeffs->ad;
      cs_real_66_t *cofbd_rij = (cs_real_66_t *)rij->bc_coeffs->bd;

      cs_real_t visci[3][3], dist[3] = {0., 0., 0.}, hint_rij = 0.0;

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

        /* ||Ki.S||^2 */
        const cs_real_t viscis =   cs_math_pow2(  visci[0][0]*nn[0]
                                                + visci[1][0]*nn[1]
                                                + visci[2][0]*nn[2])
                                 + cs_math_pow2(  visci[0][1]*nn[0]
                                                + visci[1][1]*nn[1]
                                                + visci[2][1]*nn[2])
                                 + cs_math_pow2(  visci[0][2]*nn[0]
                                                + visci[1][2]*nn[1]
                                                + visci[2][2]*nn[2]);

        /* IF.Ki.S */
        cs_real_t fikis
          = (  cs_math_3_dot_product(dist, visci[0]) * nn[0]
             + cs_math_3_dot_product(dist, visci[1]) * nn[1]
             + cs_math_3_dot_product(dist, visci[2]) * nn[2]);

        /* Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
           NB: eps =1.d-1 must be consistent with
           cs_face_anisotropic_viscosity_scalar */
        fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distbf);

        hint_rij = viscis / fikis;
      }

      /* Scalar diffusivity */
      else
        hint_rij = (visclc + visctc * cs_turb_csrij / cs_turb_cmu) / distbf;

      /* Tensor Rij (Partially or totally implicited) */

      cs_real_t fcoefa[6], fcoefb[6], fcofad[6], fcofbd[6],  fcofaf[6], fcofbf[6];
      for (int isou = 0; isou < 6; isou++) {
        fcoefa[isou] = 0.;
        fcoefb[isou] = 0.;
        fcofad[isou] = 0.;
        fcofbd[isou] = 0.;
      }

      for (cs_lnum_t isou = 0; isou < 6; isou++) {

        /* Partial (or total if coupled) implicitation */

        if (cs_glob_turb_rans_model->irijco == 1) {
          coefa_rij[f_id][isou] = 0.0;
          cofaf_rij[f_id][isou] = 0.0;
          cofad_rij[f_id][isou] = 0.0;

          for (int ii = 0; ii < 6; ii++) {
            coefb_rij[f_id][isou][ii] = alpha[ii][isou];

            if (ii == isou)
              cofbf_rij[f_id][isou][ii]
                = hint_rij * (1.0 - coefb_rij[f_id][isou][ii]);
            else
              cofbf_rij[f_id][isou][ii] = - hint_rij*coefb_rij[f_id][isou][ii];

            cofbd_rij[f_id][isou][ii] = coefb_rij[f_id][isou][ii];
          }
        }
        else if (cs_glob_turb_rans_model->iclsyr == 1) {

          for (int ii = 0; ii < 6; ii++) {
            if (ii != isou)
              fcoefa[isou] = fcoefa[isou] + alpha[ii][isou] * rijipb[f_id][ii];
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

        /* Translate into Diffusive flux BCs */
        fcofaf[isou] = -hint_rij * fcoefa[isou];
        fcofbf[isou] =  hint_rij * (1.0 - fcoefb[isou]);

      }

      if (cs_glob_turb_rans_model->irijco != 1) {
        for (int isou = 0; isou < 6; isou++) {
          coefa_rij[f_id][isou] = fcoefa[isou];
          cofaf_rij[f_id][isou] = fcofaf[isou];
          cofad_rij[f_id][isou] = fcofad[isou];

          for (int ii = 0; ii < 6; ii++) {
            coefb_rij[f_id][isou][ii] = 0;
            cofbf_rij[f_id][isou][ii] = 0;
            cofbd_rij[f_id][isou][ii] = 0;
          }
          coefb_rij[f_id][isou][isou] = fcoefb[isou];
          cofbf_rij[f_id][isou][isou] = fcofbf[isou];
          cofbd_rij[f_id][isou][isou] = fcofbd[isou];
        }
      }

    } /* if (itytur == 3) */

  } /* End loop on boundary faces */

  /* Boundary conditions on transported vectors
     ========================================== */

  for (int ii = 0; ii < n_fields; ii++) {

    cs_field_t *f_scal = cs_field_by_id(ii);

    if (!(f_scal->type & CS_FIELD_VARIABLE))
      continue;
    if (cs_field_get_key_int(f_scal, keysca) <= 0)
      continue;

    /* Get the associated turbulent flux model */
    const int turb_flux_model = cs_field_get_key_int(f_scal, kturt);
    const int turb_flux_model_type = turb_flux_model / 10;

    /* u'T' */
    if (turb_flux_model_type == 3)
      _boundary_conditions_set_coeffs_symmetry_scalar(f_scal);

    /* additional transported vectors */
    if (f_scal->dim > 1)
      _boundary_conditions_set_coeffs_symmetry_vector(f_scal);

  }

  /* Symmetry boundary conditions for mesh velocity (ALE module)
     =========================================================== */

  if (cs_glob_ale == CS_ALE_LEGACY) {

    cs_field_t *m_vel = cs_field_by_name("mesh_velocity");

    cs_real_3_t  *claale = (cs_real_3_t  *)m_vel->bc_coeffs->a;
    cs_real_33_t *clbale = (cs_real_33_t *)m_vel->bc_coeffs->b;
    cs_real_3_t  *cfaale = (cs_real_3_t  *)m_vel->bc_coeffs->af;
    cs_real_33_t *cfbale = (cs_real_33_t *)m_vel->bc_coeffs->bf;

    int *icodcl_displ = m_vel->bc_coeffs->icodcl;
    cs_equation_param_t *eqp_displ = cs_field_get_equation_param(m_vel);

    const cs_real_t   *cpro_visma_s = NULL;
    const cs_real_6_t *cpro_visma_v = NULL;

    if (eqp_displ->idften & CS_ISOTROPIC_DIFFUSION)
      cpro_visma_s = CS_F_(vism)->val;
    else if (eqp_displ->idften & CS_ANISOTROPIC_DIFFUSION)
      cpro_visma_v = (const cs_real_6_t *)CS_F_(vism)->val;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      if (icodcl_displ[f_id] != 4)
        continue;

      const cs_lnum_t c_id = b_face_cells[f_id];

      /* For a sliding boundary, the normal velocity is enforced to zero
         whereas the other components have an Homogenous Neumann
         NB: no recontruction in I' here */

      /* Geometrical quantities */
      const cs_real_t distbf = b_dist[f_id];
      const cs_real_t *nn = b_face_u_normal[f_id];

      /* Physical properties */
      cs_real_6_t hintv = {0., 0., 0., 0., 0., 0.};

      if (eqp_displ->idften & CS_ISOTROPIC_DIFFUSION) {

        hintv[0] = cpro_visma_s[c_id]/distbf;
        hintv[1] = cpro_visma_s[c_id]/distbf;
        hintv[2] = cpro_visma_s[c_id]/distbf;
        hintv[3] = 0.0;
        hintv[4] = 0.0;
        hintv[5] = 0.0;

      }
      else if (eqp_displ->idften & CS_ANISOTROPIC_LEFT_DIFFUSION) {
        for (int ii = 0; ii < 6; ii++)
          hintv[ii] = cpro_visma_v[c_id][ii] / distbf;
      }

      /* Coupled solving of the velocity components */

      /* Gradient BCs */
      claale[f_id][0] = 0.0;
      claale[f_id][1] = 0.0;
      claale[f_id][2] = 0.0;

      clbale[f_id][0][0] = 1.0 - nn[0] * nn[0];
      clbale[f_id][1][1] = 1.0 - nn[1] * nn[1];
      clbale[f_id][2][2] = 1.0 - nn[2] * nn[2];

      clbale[f_id][0][1] = - nn[0] * nn[1];
      clbale[f_id][1][0] = - nn[1] * nn[0];
      clbale[f_id][0][2] = - nn[0] * nn[2];
      clbale[f_id][2][0] = - nn[2] * nn[0];
      clbale[f_id][1][2] = - nn[1] * nn[2];
      clbale[f_id][2][1] = - nn[2] * nn[1];

      /* Flux BCs */
      cfaale[f_id][0] = 0.0;
      cfaale[f_id][1] = 0.0;
      cfaale[f_id][2] = 0.0;

      cs_real_t rnn[6];
      rnn[0] = nn[0] * nn[0];
      rnn[1] = nn[1] * nn[1];
      rnn[2] = nn[2] * nn[2];
      rnn[3] = nn[0] * nn[1];
      rnn[4] = nn[1] * nn[2];
      rnn[5] = nn[0] * nn[2];

      cs_real_t htnn[6];
      cs_math_sym_33_product(hintv, rnn, htnn);
      cfbale[f_id][0][0] = htnn[0];
      cfbale[f_id][1][1] = htnn[1];
      cfbale[f_id][2][2] = htnn[2];
      cfbale[f_id][0][1] = htnn[3];
      cfbale[f_id][1][0] = htnn[3];
      cfbale[f_id][0][2] = htnn[5];
      cfbale[f_id][2][0] = htnn[5];
      cfbale[f_id][1][2] = htnn[4];
      cfbale[f_id][2][1] = htnn[4];

    } /* End loop on boundary faces */

  } /* End for ALE */
}

/*---------------------------------------------------------------------------- */

END_C_DECLS
