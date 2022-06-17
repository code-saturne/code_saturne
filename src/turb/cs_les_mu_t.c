/*============================================================================
 * Turbulent viscosity for LES models.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gradient.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_les_filter.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_les_mu_t.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the turbulent viscosity for
 *        a dynamic Smagorinsky LES model
 *
 * \f[ smago = \dfrac{L_{ij}M_{ij}}{M_{ij}M_{ij}} \f]
 *
 * \f[ \mu_T = \rho smago L^2  \sqrt{2 S_{ij}S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * Please refer to the
 * <a href="../../theory.pdf#dynsmago"><b>dynamic Smagorinsky model</b></a>
 * section of the theory guide for more informations.
 *
 * \param[out]  gradv the computed velocity gradients
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_dyn(cs_real_33_t  *gradv)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const int irovar = cs_glob_fluid_properties->irovar;

  /* Initialization
   * ============== */

  /*  Map field arrays */

  cs_real_t *visct =  CS_F_(mu_t)->val;
  const cs_real_t *crom  = CS_F_(rho)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  cs_real_t *cpro_smago
    = cs_field_by_name("smagorinsky_constant^2")->val;

  /* For the calculation of the viscosity of the sub-mesh */

  const cs_real_t xfil   = cs_turb_xlesfl;
  const cs_real_t xfil2  = cs_turb_xlesfd;
  const cs_real_t xa     = cs_turb_ales;
  const cs_real_t xb     = cs_turb_bles;
  const cs_real_t sqrt_2 = sqrt(2.0);
  const cs_real_t xsmgmx = cs_turb_csmago_max;
  const cs_real_t xsmgmn = cs_turb_csmago_min;

  /*  Allocate some work arrays */

  cs_real_t *w0, *w1, *xro, *xrof;
  cs_real_6_t *xmij;

  BFT_MALLOC(w0, n_cells_ext, cs_real_t);
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(xmij, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(xro, n_cells_ext, cs_real_t);
  BFT_MALLOC(xrof, n_cells_ext, cs_real_t);

  /* Take into account variable density case: Favre filtering
   * Constant density case: Reynolds filtering */

  if (irovar == 1)
    cs_array_copy_real(n_cells, 1, crom, xro);
  else
    cs_array_set_value_real(n_cells, 1, 1.0, xro);

  /* In case of constant density, xrof always 1.0 */

  cs_les_filter(1, xro, xrof);

  /* Calculation of velocity gradient and of
   * S11^2+S22^2+S33^2+2*(S12^2+S13^2+S23^2)
   *======================================== */

  /* Allocate temporary arrays for gradients calculation */

  cs_real_t *s_n, *sf_n;
  cs_real_6_t *w61, *w62;

  BFT_MALLOC(s_n, n_cells_ext, cs_real_t);
  BFT_MALLOC(sf_n, n_cells_ext, cs_real_t);
  BFT_MALLOC(w61, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(w62, n_cells_ext, cs_real_6_t);

  cs_field_gradient_vector(CS_F_(vel),
                           false, // no use_previous_t
                           1,     // inc
                           gradv);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    /* gradv[c_id][xyz][uvw] */

    const cs_real_t s11  = gradv[c_id][0][0];
    const cs_real_t s22  = gradv[c_id][1][1];
    const cs_real_t s33  = gradv[c_id][2][2];
    const cs_real_t dudy = gradv[c_id][0][1];
    const cs_real_t dudz = gradv[c_id][0][2];
    const cs_real_t dvdx = gradv[c_id][1][0];
    const cs_real_t dvdz = gradv[c_id][1][2];
    const cs_real_t dwdx = gradv[c_id][2][0];
    const cs_real_t dwdy = gradv[c_id][2][1];

    /* In the case of constant density, s11+s22+s33 is zero */

    xmij[c_id][0] = s11-irovar*1.0/3.0*(s11+s22+s33);
    xmij[c_id][1] = s22-irovar*1.0/3.0*(s11+s22+s33);
    xmij[c_id][2] = s33-irovar*1.0/3.0*(s11+s22+s33);
    xmij[c_id][3] = 0.5*(dudy+dvdx);
    xmij[c_id][4] = 0.5*(dudz+dwdx);
    xmij[c_id][5] = 0.5*(dvdz+dwdy);

    s_n[c_id] = sqrt_2 * sqrt(  cs_math_pow2(xmij[c_id][0])
                              + cs_math_pow2(xmij[c_id][1])
                              + cs_math_pow2(xmij[c_id][2])
                              + 2.0 * (  cs_math_pow2(xmij[c_id][3])
                                       + cs_math_pow2(xmij[c_id][4])
                                       + cs_math_pow2(xmij[c_id][5])));

    for (cs_lnum_t i = 0; i < 6; i++)
      w62[c_id][i] = xro[c_id] * xmij[c_id][i];
  }

  /* w62 temperarily contains rho*S */

  cs_les_filter(6, (cs_real_t*)w62, (cs_real_t*)w61);

  /* w61 <rho*S>/<rho>, sf_n is ||<rho*S>/<rho>|| */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t i = 0; i < 6; i++)
      w61[c_id][i] /= xrof[c_id];

    sf_n[c_id] = sqrt_2 * sqrt(  cs_math_pow2(w61[c_id][0])
                               + cs_math_pow2(w61[c_id][1])
                               + cs_math_pow2(w61[c_id][2])
                               +2.0 * (  cs_math_pow2(w61[c_id][3])
                                       + cs_math_pow2(w61[c_id][4])
                                       + cs_math_pow2(w61[c_id][5])));
  }

  /* Here XMIJ contains Sij
   *   S_n  contains ||S||
   *       sqrt(2)*sqrt(S11^2+S22^2+S33^2+2(S12^2+S13^2+S23^2))
   *   Sf_n contains ||SF||
   *       sqrt(2)*sqrt(S11F^2+S22F^2+S33F^2+2(S12F^2+S13F^2+S23F^2)) */

  /* Calculation of Mij
   *===================== */

  /* Reuse xmij as temporary array */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t delta = xfil * pow(xa*cell_vol[c_id], xb);
    w0[c_id] = delta;
    for (cs_lnum_t ii = 0; ii < 6; ii++)
      xmij[c_id][ii] *= -2.0 * xro[c_id] * cs_math_pow2(delta) * s_n[c_id];
  }

  /* w62 now contains <-2*rho*delta^2*||S||*S> */
  cs_les_filter(6, (cs_real_t*)xmij, (cs_real_t*)w62);

  /* Now compute final xmij value: M_ij = alpha_ij - beta_ij */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t delta = w0[c_id];
    const cs_real_t deltaf = xfil2*delta;
    for (cs_lnum_t ii = 0; ii < 6; ii++)
      xmij[c_id][ii] =    -2.0 * xro[c_id] * cs_math_pow2(deltaf)
                        * sf_n[c_id] * w61[c_id][ii]
                      - w62[c_id][ii];
  }

  BFT_FREE(w61);
  BFT_FREE(w62);

  /* Calculation of the dynamic Smagorinsky constant
   *================================================ */

  /* Allocate work arrays */
  cs_real_t *w2, *w3, *w4, *w5;
  cs_real_t *w6, *w7, *w8, *w9;

  BFT_MALLOC(w2, n_cells_ext, cs_real_t);
  BFT_MALLOC(w3, n_cells_ext, cs_real_t);
  BFT_MALLOC(w4, n_cells_ext, cs_real_t);
  BFT_MALLOC(w5, n_cells_ext, cs_real_t);
  BFT_MALLOC(w6, n_cells_ext, cs_real_t);
  BFT_MALLOC(w7, n_cells_ext, cs_real_t);
  BFT_MALLOC(w8, n_cells_ext, cs_real_t);
  BFT_MALLOC(w9, n_cells_ext, cs_real_t);

  /* Filtering the velocity and its square */

  // U^2
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][0]*vel[c_id][0];
  cs_les_filter(1, w0, w1);

  // V^2
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][1]*vel[c_id][1];
  cs_les_filter(1, w0, w2);

  // W^2
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][2]*vel[c_id][2];
  cs_les_filter(1, w0, w3);

  // UV
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][0]*vel[c_id][1];
  cs_les_filter(1, w0, w4);

  // UW
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][0]*vel[c_id][2];
  cs_les_filter(1, w0, w5);

  // VW
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][1]*vel[c_id][2];
  cs_les_filter(1, w0, w6);

  // U
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][0];
  cs_les_filter(1, w0, w7);
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w7[c_id] = w7[c_id]/xrof[c_id];

  // V
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][1];
  cs_les_filter(1, w0, w8);
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w8[c_id] = w8[c_id]/xrof[c_id];

  // W
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w0[c_id] = xro[c_id]*vel[c_id][2];
  cs_les_filter(1, w0, w9);
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++)
    w9[c_id] = w9[c_id]/xrof[c_id];

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    /* Calculation of Lij */
    const cs_real_t xl11 = w1[c_id]-xrof[c_id]*w7[c_id]*w7[c_id];
    const cs_real_t xl22 = w2[c_id]-xrof[c_id]*w8[c_id]*w8[c_id];
    const cs_real_t xl33 = w3[c_id]-xrof[c_id]*w9[c_id]*w9[c_id];
    const cs_real_t xl12 = w4[c_id]-xrof[c_id]*w7[c_id]*w8[c_id];
    const cs_real_t xl13 = w5[c_id]-xrof[c_id]*w7[c_id]*w9[c_id];
    const cs_real_t xl23 = w6[c_id]-xrof[c_id]*w8[c_id]*w9[c_id];

    const cs_real_t xm11 = xmij[c_id][0];
    const cs_real_t xm22 = xmij[c_id][1];
    const cs_real_t xm33 = xmij[c_id][2];
    const cs_real_t xm12 = xmij[c_id][3];
    const cs_real_t xm13 = xmij[c_id][4];
    const cs_real_t xm23 = xmij[c_id][5];

    /* Calculation of Mij :: Lij */
    w1[c_id] =   xm11*xl11 + 2.0*xm12*xl12 + 2.0*xm13*xl13
               +                 xm22*xl22 + 2.0*xm23*xl23
               +                                 xm33*xl33;

    /* Calculation of Mij :: Mij */
    w2[c_id] =   xm11*xm11 + 2.0*xm12*xm12 + 2.0*xm13*xm13
               +                 xm22*xm22 + 2.0*xm23*xm23
               +                                 xm33*xm33;

  }

  BFT_FREE(xmij);

  /* By default we compute a local average of the numerator and of the
     denominator, then only compute  the quotient. The user can override
     this  in cs_user_physical_properties_smagorinsky_c. */

  cs_les_filter(1, w1, w3);
  cs_les_filter(1, w2, w4);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (fabs(w4[c_id]) <= cs_math_epzero)
      cpro_smago[c_id] = xsmgmx;
    else
      cpro_smago[c_id] = w3[c_id]/w4[c_id];
  }

  cs_user_physical_properties_smagorinsky_c(cs_glob_domain, w1, w2);

  cs_gnum_t iclipc = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (cpro_smago[c_id] >= xsmgmx) {
      cpro_smago[c_id] = xsmgmx;
      iclipc += 1;
    }
    else if(cpro_smago[c_id] <= xsmgmn) {
      cpro_smago[c_id] = xsmgmn;
      iclipc += 1;
    }
  }

  /* Calculation of (dynamic) viscosity
   * ================================== */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t coef = cpro_smago[c_id];
    const cs_real_t delta = xfil * pow(xa*cell_vol[c_id], xb);
    visct[c_id] = crom[c_id]*coef*cs_math_pow2(delta)*s_n[c_id];
  }

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(CS_F_(vel));

  if (eqp->verbosity >= 1) {
    cs_real_t smagma = -1.0e12;
    cs_real_t smagmi =  1.0e12;
    cs_real_t smagmy =  0.0;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smagma = fmax(smagma, cpro_smago[c_id]);
      smagmi = fmin(smagmi, cpro_smago[c_id]);
      smagmy += cpro_smago[c_id]*cell_vol[c_id];
    }
    cs_parall_max(1, CS_REAL_TYPE, &smagma);
    cs_parall_min(1, CS_REAL_TYPE, &smagmi);
    cs_parall_sum(1, CS_REAL_TYPE, &smagmy);
    cs_parall_counter(&iclipc, 1);

    smagmy /= cs_glob_mesh_quantities->tot_vol;
    cs_log_printf(CS_LOG_DEFAULT,
                  _("N. clippings of the Smagorinsky constant %lu\n"
                    " --- Information on the squared Smagorinsky constant\n"
                    " --------------------------------\n"
                    " Mean value  Min value  Max value\n"
                    " --------------------------------\n"
                    " %e12.4, %e12.4, %e12.4\n"
                    " --------------------------------\n"),
                  (unsigned long)iclipc, smagmy, smagmi, smagma);
  }

  /* Clipping of the turbulent viscosity
   * =================================== */

  /* Clip turbulent viscosity so that it is always positive. */

  iclipc = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (visct[c_id] < 0.0) {
      visct[c_id] = 0.0;
      iclipc += 1;
    }
  }

  if (eqp->iwarni >= 1) {
    cs_parall_counter(&iclipc, 1);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("N. clippings for turbulent viscosity (mu_t>0): %lu\n"),
                  (unsigned long)iclipc);
  }

  /* Scalar turbulent model
   *======================= */

  /* In case of gas combustion, the SGS scalar flux constant and the turbulent
   * diffusivity are only evaluated with the mixture fraction, then applied
   * automatically to the other scalar equations */

  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
  const int key_sgs_sca_coef = cs_field_key_id("sgs_scalar_flux_coef_id");

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *fld = cs_field_by_id(f_id);

    int sc_id = cs_field_get_key_int(fld, keysca)-1;
    int iscavr = cs_field_get_key_int(fld, kscavr);
    /* For variance of a scalar, the turbulent diffusivity is not computed */
    if ((sc_id < 0) || (iscavr > -1))
      continue;

    /* For any scalar other than the mixture fraction in diffusion flames,
     * Dt is not computed either. TODO Soot may be an exception */
    if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] > -1) {
      if (fld != CS_F_(fm))
        continue;
    }

    int t_dif_id = cs_field_get_key_int(fld, key_turb_diff);
    int sca_dync_id = cs_field_get_key_int(fld, key_sgs_sca_coef);

    if ((t_dif_id >= 0) && (sca_dync_id >= 0)) {
      cs_real_t *cpro_turb_diff = cs_field_by_id(t_dif_id)->val;
      cs_real_t *cpro_sca_dync = cs_field_by_id(sca_dync_id)->val;

      /* Compute the Mi for scalar
       * ========================= */

      cs_real_3_t *grads, *gradsf;
      BFT_MALLOC(grads, n_cells_ext, cs_real_3_t);
      BFT_MALLOC(gradsf, n_cells_ext, cs_real_3_t);

      cs_real_t *coefas = fld->bc_coeffs->a;
      cs_real_t *coefbs = fld->bc_coeffs->b;

      cs_real_t *cvar_sca = fld->val;
      cs_field_gradient_scalar(fld,
                               false,     /* use previous t   */
                               1,         /* not on increment */
                               false,     /* recompute_cocg   */
                               grads);

      /* compute grad (<rho.Y>/<rho>) */
      cs_real_3_t *scami, *scamif;
      BFT_MALLOC(scami, n_cells_ext, cs_real_3_t);
      BFT_MALLOC(scamif, n_cells_ext, cs_real_3_t);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w0[c_id] = cvar_sca[c_id]*xro[c_id];
      cs_les_filter(1, w0, w4);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w4[c_id] = w4[c_id]/xrof[c_id];

      const cs_equation_param_t *eqp_fld
        = cs_field_get_equation_param_const(fld);

      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

      cs_gradient_type_by_imrgra(eqp_fld->imrgra,
                                 &gradient_type,
                                 &halo_type);

      cs_gradient_scalar("Work array",
                         gradient_type,
                         halo_type,
                         1,     /* inc */
                         true,  /* iccocg */
                         eqp_fld->nswrgr,
                         0,
                         0,
                         1,     /* w_stride */
                         eqp_fld->iwarni,
                         eqp_fld->imligr,
                         eqp_fld->epsrgr,
                         eqp_fld->climgr,
                         NULL,
                         coefas,
                         coefbs,
                         w4,
                         NULL,
                         NULL, /* internal coupling */
                         gradsf);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t delta  = xfil * pow(xa*cell_vol[c_id],xb);
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          scami[c_id][ii]
            = -xro[c_id]*cs_math_pow2(delta)*s_n[c_id]*grads[c_id][ii];
      }

      cs_les_filter(3, (cs_real_t *)scami, (cs_real_t *)scamif);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t deltaf  = xfil2* pow(xa*cell_vol[c_id],xb);
        for (int ii = 0; ii < 3; ii++)
          scami[c_id][ii] = - cs_math_pow2(deltaf)*xrof[c_id]
                              * sf_n[c_id]*gradsf[c_id][ii]
                            - scamif[c_id][ii];
      }

      /* Compute the Li for scalar
       * ========================= */

      /* rho*U*Y */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w0[c_id] = xro[c_id]*vel[c_id][0]*cvar_sca[c_id];
      }
      cs_les_filter(1, w0, w1);

      /* rho*V*Y */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w0[c_id] = xro[c_id]*vel[c_id][1]*cvar_sca[c_id];
      }
      cs_les_filter(1, w0, w2);

      /* rho*W*Y */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w0[c_id] = xro[c_id]*vel[c_id][2]*cvar_sca[c_id];
      }
      cs_les_filter(1, w0, w3);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t scal1 = w1[c_id] - xrof[c_id]*w7[c_id]*w4[c_id];
        const cs_real_t scal2 = w2[c_id] - xrof[c_id]*w8[c_id]*w4[c_id];
        const cs_real_t scal3 = w3[c_id] - xrof[c_id]*w9[c_id]*w4[c_id];

        w1[c_id] =   scal1*scami[c_id][0]
                   + scal2*scami[c_id][1]
                   + scal3*scami[c_id][2];
        w2[c_id] =   cs_math_pow2(scami[c_id][0])
                   + cs_math_pow2(scami[c_id][1])
                   + cs_math_pow2(scami[c_id][2]);
      }

      cs_les_filter(1, w1, w3);
      cs_les_filter(1, w2, w4);

      /*
       * Compute the SGS flux coefficient and SGS diffusivity
       * Cs >= 0, Dt >=0
       */

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if(fabs(w4[c_id]) <= cs_math_epzero)
          cpro_sca_dync[c_id] = 0.0;
        else
          cpro_sca_dync[c_id] = fmax(w3[c_id]/w4[c_id], 0.0);

        const cs_real_t delta  = xfil * pow(xa*cell_vol[c_id], xb);
        cpro_turb_diff[c_id] =   crom[c_id] * cpro_sca_dync[c_id]
                               * cs_math_pow2(delta) * s_n[c_id];
      }

      BFT_FREE(scami);
      BFT_FREE(scamif);
      BFT_FREE(grads);
      BFT_FREE(gradsf);
    }
  }

  /* Free memory */
  BFT_FREE(s_n);
  BFT_FREE(sf_n);

  BFT_FREE(w9);
  BFT_FREE(w8);
  BFT_FREE(w7);
  BFT_FREE(w6);
  BFT_FREE(w5);
  BFT_FREE(w4);
  BFT_FREE(w3);
  BFT_FREE(w2);
  BFT_FREE(w1);
  BFT_FREE(w0);

  BFT_FREE(xro);
  BFT_FREE(xrof);
}

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for a Smagorinsky LES model
 *
 * \f[ \mu_T = \rho (C_{S} l)^2  \sqrt{2 S_{ij}S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * \param[out]  gradv  computed velocity gradients (may be used by caller)
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_const(cs_real_33_t  *gradv)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  /* Initialization
   * ============== */

  cs_real_t *visct =  CS_F_(mu_t)->val;
  const cs_real_t *crom  = CS_F_(rho)->val;

  const cs_real_t xfil = cs_turb_xlesfl;
  const cs_real_t xa = cs_turb_ales;
  const cs_real_t xb = cs_turb_bles;

  const cs_real_t coef = cs_math_pow2(cs_turb_csmago)*sqrt(2.0);

  /* We need the velocity gradient */

  cs_field_gradient_vector(CS_F_(vel),
                           false,  /* no use_previous_t */
                           1,      /* inc */
                           gradv);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    /* S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2) */

    const cs_real_t s11  = gradv[c_id][0][0];
    const cs_real_t s22  = gradv[c_id][1][1];
    const cs_real_t s33  = gradv[c_id][2][2];
    const cs_real_t dudy = gradv[c_id][0][1];
    const cs_real_t dvdx = gradv[c_id][1][0];
    const cs_real_t dudz = gradv[c_id][0][2];
    const cs_real_t dwdx = gradv[c_id][2][0];
    const cs_real_t dvdz = gradv[c_id][1][2];
    const cs_real_t dwdy = gradv[c_id][2][1];

    visct[c_id] =   cs_math_pow2(s11)
                  + cs_math_pow2(s22)
                  + cs_math_pow2(s33)
                  + 0.5 * (  cs_math_pow2(dudy+dvdx)
                           + cs_math_pow2(dudz+dwdx)
                           + cs_math_pow2(dvdz+dwdy));

    /* Calculation of (dynamic) viscosity */

    cs_real_t delta  = xfil* (xa*pow(cell_vol[c_id],xb));
    delta = coef*cs_math_pow2(delta);
    visct[c_id] = crom[c_id] * delta * sqrt(visct[c_id]);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the turbulent viscosity for the WALE LES model.
 *
 * The turbulent viscosity is:
 * \f$ \mu_T = \rho (C_{wale} L)^2 * \dfrac{(\tens{S}:\tens{Sd})^{3/2}}
 *                                         {(\tens{S} :\tens{S})^(5/2)
 *                                         +(\tens{Sd}:\tens{Sd})^(5/4)} \f$
 * with \f$ \tens{S}  = \frac{1}{2}(\gradt \vect{u} + \transpose{\gradt \vect{u}})\f$
 * and  \f$ \tens{Sd} = \deviator{(\symmetric{(\tens{S}^2)})}\f$
 *
 * \param[out]  gradv the computed velocity gradients
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_wale(cs_real_33_t *restrict gradv)

{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  /* Initialization
   * ============== */

  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *crom = CS_F_(rho)->val;

  /* Computation of the velocity gradient
   * ==================================== */

  cs_field_gradient_vector(CS_F_(vel),
                           false,  /* no use_previous_t */
                           1,      /* inc */
                           gradv);

  cs_real_t dudx[3][3], kdelta[3][3], g2[3][3];

  /* Kronecker delta Dij */

  for (cs_lnum_t i = 0; i < 3; i++) {
    for (cs_lnum_t j = 0; j < 3; j++) {
      if (i == j)
        kdelta[i][j] = 1;
      else
        kdelta[i][j] = 0;
    }
  }

  const cs_real_t coef = sqrt(2.0)*cs_math_pow2(cs_turb_cwale);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    /* Dudx is interleaved, but not gradv...
     * gradv[c_id][xyz][uvw]
     * ====================================== */

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        dudx[i][j] = gradv[c_id][i][j];
    }

    cs_real_t s = 0, trace_g2 = 0;
    const cs_real_t third = 1.0/3.0;

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        /* s = 1/4 * (dUi/dXj + dUj/dXi) * (dUi/dXj + dUj/dXi) */
        s += 0.25*cs_math_pow2(dudx[i][j] + dudx[j][i]);

        /* g2 is the square tensor of the velocity gradient */
        g2[i][j] = 0.0;
        for (cs_lnum_t k = 0; k < 3; k++)
          g2[i][j] += dudx[i][k]*dudx[k][j];
      }

      trace_g2 += g2[i][i];
    }

    cs_real_t sd = 0;
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        /* Traceless symmetric part of the square of the velocity gradient tensor
         *   Sijd =   0.5*(dUi/dXk dUk/dXj + dUj/dXk dUk/dXi)
         *          - 1/3 Dij dUk/dXl dUl/dXk */
        cs_real_t sijd = 0.50*(g2[i][j]+g2[j][i])-third*kdelta[i][j]*trace_g2;
        sd += cs_math_pow2(sijd);
      }
    }

    /* Turbulent viscosity
     * =================== */

    /* Turbulent inverse time scale
     *   = (Sijd Sijd)^3/2 / [ (Sij Sij)^5/2 + (Sijd Sijd)^5/4 ] */

    cs_real_t sinv = pow(s, 2.5) + pow(sd, 1.25);
    cs_real_t con = 0;
    if (sinv > 0)
      con = pow(sd, 1.5)/sinv;

    cs_real_t delta = cs_turb_xlesfl * pow(cs_turb_ales*cell_vol[c_id],
                                           cs_turb_bles);
    delta = coef * cs_math_pow2(delta);

    visct[c_id] = crom[c_id] * delta * con;

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
