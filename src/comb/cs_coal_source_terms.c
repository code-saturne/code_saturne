/*============================================================================
 * Coal combustion model: source term computation
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_coal.h"
#include "cs_coal_ht_convert.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal_source_terms.c

  \brief Coal combustion model: source terms computation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute variance source terms for pulverized coal flame.
 *
 * Production and dissipation source term for variance
 * (explicit and implicit balance).
 *
 * \param[in]      fld_scal  pointer to scalar field
 * \param[in,out]  smbrs     explicit second member
 * \param[in,out]  rovsdt    implicit diagonal part
 */
/*----------------------------------------------------------------------------*/

static void
_coal_fp2st(const cs_field_t  *fld_scal,
            cs_real_t          smbrs[],
            cs_real_t          rovsdt[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  const cs_coal_model_t *cm = cs_glob_coal_model;

  /* Initialization
     -------------- */

  const int krvarfl = cs_field_key_id("variance_dissipation");
  const cs_real_t rvarfl = cs_field_get_key_double(fld_scal, krvarfl);

  cs_real_t *f1f2;
  BFT_MALLOC(f1f2, n_cells_ext, cs_real_t);
  cs_array_real_fill_zero(n_cells, f1f2);

  // The variance is not associated to a scalar but to f1+f2

  const cs_real_t *cvara_scal = fld_scal->val_pre;

  // Physical quantities

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  const cs_real_t *cpro_rom1 = cs_field_by_id(cm->irom1)->val;

  const cs_real_t *cvara_k = NULL;
  const cs_real_t *cvara_ep = NULL;
  const cs_real_t *cvara_omg = NULL;
  const cs_real_6_t *cvara_rij = NULL;

  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 5) {
    cvara_k = CS_F_(k)->val_pre;
    cvara_ep = CS_F_(eps)->val_pre;
  }
  else if (cs_glob_turb_model->itytur == 3) {
    cvara_rij = (const cs_real_6_t *)(CS_F_(k)->val_pre);
    cvara_ep = CS_F_(eps)->val_pre;
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {
    cvara_k = CS_F_(k)->val_pre;
    cvara_omg = CS_F_(omg)->val_pre;
  }

  /* Account for production and dissipation source terms by gradients
     ---------------------------------------------------------------- */

  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 3
      || cs_glob_turb_model->itytur == 5
      || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {

    // For lack of information on F1M+F2M, we take the same options as for F1M[0].

    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(cs_field_by_id(cm->if1m[0]));

    // Compute x1

    cs_real_t *x1;
    BFT_MALLOC(x1, n_cells_ext, cs_real_t);
    cs_array_real_set_scalar(n_cells_ext, 1.0, x1);

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t xmash = cm->xmash[class_id];
      const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
      const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
      const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[class_id])->val;
      const cs_real_t *cvar_xwtcl = NULL;
      if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
        cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        x1[c_id] -= (  cvar_xchcl[c_id]
                     + cvar_xckcl[c_id]
                     + cvar_xnpcl[c_id]*xmash);
        if (cvar_xwtcl != NULL)
          x1[c_id] -= cvar_xwtcl[c_id];
      }
    }

    // Compute F=F1+F2

    for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
      const cs_real_t *cvar_f1m_c = cs_field_by_id(cm->if1m[coal_id])->val;
      const cs_real_t *cvar_f2m_c = cs_field_by_id(cm->if2m[coal_id])->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        f1f2[c_id] += cvar_f1m_c[c_id] + cvar_f2m_c[c_id];
      }
    }

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      f1f2[c_id] /= x1[c_id];
    }

    // Compute gradient of f1f2

    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
    cs_gradient_type_by_imrgra(eqp->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_real_3_t *grad;
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    cs_gradient_scalar("f1f2",
                       gradient_type,
                       halo_type,
                       1.,            /* inc */
                       eqp->nswrgr,
                       0,             /* iphydp */
                       1,             /* w_stride */
                       eqp->verbosity,
                       (cs_gradient_limit_t)eqp->imligr,
                       eqp->epsrgr,
                       eqp->climgr,
                       NULL,          /* f_ext */
                       NULL,          /* bc_coeffs */
                       f1f2,
                       NULL,          /* c_weight */
                       NULL,          /* cpl */
                       grad);

    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    cs_real_t turb_schmidt = cs_field_get_key_double(fld_scal, ksigmas);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t xk = 0, xe = 0;
      if (cvara_rij != NULL) {
        xk = 0.5 * (  cvara_rij[c_id][0]
                    + cvara_rij[c_id][1]
                    + cvara_rij[c_id][2]);
        xe = cvara_ep[c_id];
      }
      else if (cvara_ep != NULL) {
        xk = cvara_k[c_id];
        xe = cvara_ep[c_id];
      }
      else if (cvara_omg != NULL) {
        xk = cvara_k[c_id];
        xe = cs_turb_cmu * xk * cvara_omg[c_id];
      }

      cs_real_t rhovst = cpro_rom1[c_id]*xe/(xk*rvarfl)*cell_f_vol[c_id];
      rovsdt[c_id] += fmax(0, rhovst);
      smbrs[c_id] += (  2.0 * visct[c_id] * cell_f_vol[c_id] / turb_schmidt
                      * cs_math_3_square_norm(grad[c_id]) * x1[c_id])
                     - rhovst * cvara_scal[c_id];

    }

    BFT_FREE(grad);
    BFT_FREE(x1);

  }

  /* Account for source terms relative to interfacial exchanges
     ---------------------------------------------------------- */

  // 2 versions available
  //   iold = 1 => old version
  //   iold = 2 => new version

  const int iold = 1;

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
    const cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;

    const cs_real_t diamdv = cm->diam20[class_id];

    if (iold == 1) {
      const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[class_id])->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t gdev  = -crom[c_id] * (cpro_cgd1[c_id] + cpro_cgd2[c_id])
                                      * cvar_xchcl[c_id];

        if (cvar_xnpcl[c_id] > cs_coal_epsilon) {
          cs_real_t  fsd =  1.0 - (1.0-f1f2[c_id])
                           * exp((cvar_xchcl[c_id]*(  cpro_cgd1[c_id]
                                                    + cpro_cgd2[c_id]))
                                 / (  2.0 * cs_math_pi *2.77e-4 * diamdv
                                    * cvar_xnpcl[c_id] * crom[c_id]));
          cs_real_t fdev = 1.0;

          /* Explicit ST */

          if (  (fsd - f1f2[c_id]) * (2.0*fdev - fsd - f1f2[c_id])
              > cs_coal_epsilon) {
            smbrs[c_id] +=   cell_f_vol[c_id] * gdev
                           * (fsd - f1f2[c_id]) * (2.0*fdev - fsd - f1f2[c_id]);
          }

        }
      }
    }
    else  { /* if iold == 2 */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t gdev  = -crom[c_id] * (cpro_cgd1[c_id] + cpro_cgd2[c_id])
                                      * cvar_xchcl[c_id];

        cs_real_t aux = gdev * cs_math_pow2(1.0 - f1f2[c_id]);

        /* Implicit ST: for now implicit in simple manner */

        cs_real_t rhovst = 0;
        if (fabs(f1f2[c_id] * (1.0-f1f2[c_id])) > cs_coal_epsilon) {
          rhovst =   aux*cvara_scal[c_id]
                   / cs_math_pow2(f1f2[c_id] * (1-f1f2[c_id]))
                   * cell_f_vol[c_id];
        }
        rovsdt[c_id] += fmax(0., rhovst);

        /* Explicit ST */

        smbrs[c_id] += aux*cell_f_vol[c_id] - rhovst*cvara_scal[c_id];
      }
    }
  } /* Loop on classes */

  // Free memory
  BFT_FREE(f1f2);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute scalar source terms for pulverized coal flame.
 *
 * \warning  the treatement is different from that of cs_user_source_terms.
 *
 * We solve: \f[ rovsdt D(var) = smbrs \f]
 *
 * rovsdt and smbrs already contain eventual user source terms.
 * So they have to be incremented and not erased.
 *
 * For stability reasons, only positive terms can be added in rovsdt.
 * There is no contraint for smbrs.
 *
 * In the case of a source term in \f$ cexp + cimp var \f$, it has to be written:
 *        - \f$ smbrs  = smbrs  + cexp + cimp var \f$
 *        - \f$ rovsdt = rovsdt + \max(-cimp,0) \f$
 *
 * Here are \f$ rovsdt \f$ and \f$ smbrs \f$ (they contain \f$ \rho volume\f$)
 *    smbrs in kg variable/s:
 *     \c i.e.: - for velocity            \f$ kg . m . s^{-2} \f$
 *              - for temperature         \f$ kg . [degres] . s^{-1} \f$
 *              - for enthalpy            \f$ J . s^{-1} \f$
 *              - rovsdt                  \f$ kg . s^{-1} \f$
 *
 * \param[in]      fld_id  scalar field id
 * \param[in,out]  smbrs   explicit second member
 * \param[in,out]  rovsdt  implicit diagonal part
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_source_terms_scalar(int        fld_id,
                            cs_real_t  smbrs[],
                            cs_real_t  rovsdt[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  const cs_coal_model_t *cm = cs_glob_coal_model;

  /* Initialization
   * -------------- */

  const cs_field_t *fld_scal = cs_field_by_id(fld_id);
  const cs_real_t *cvara_var = fld_scal->val_pre;
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(fld_scal);

  const int keyccl = cs_field_key_id("scalar_class");

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *cpro_cp = NULL;

  int icp = cs_glob_fluid_properties->icp;
  if (icp >= 0)
    cpro_cp = cs_field_by_id(icp)->val;

  const cs_real_t *cpro_temp = CS_F_(t)->val;
  const cs_real_t *cpro_rom1 = cs_field_by_id(cm->irom1)->val;

  const int ichx1 = cm->ichx1 -1;
  const int ichx2 = cm->ichx2 -1;
  const int ico = cm->ico -1;
  const int ih2s = cm->ih2s -1;
  const int ihy =cm->ihy -1;
  const int ihcn = cm->ihcn -1;
  const int inh3 = cm->inh3 -1;
  const int io2 = cm->io2 -1;
  const int ico2 = cm->ico2 -1;
  const int ih2o = cm->ih2o -1;
  const int in2 = cm->in2 -1;

  const cs_real_t *cpro_yo2 = cs_field_by_id(cm->iym1[io2])->val;
  const cs_real_t *cpro_yco2 = cs_field_by_id(cm->iym1[ico2])->val;
  const cs_real_t *cpro_yco = cs_field_by_id(cm->iym1[ico])->val;
  const cs_real_t *cpro_yh2o = cs_field_by_id(cm->iym1[ih2o])->val;
  const cs_real_t *cpro_mmel = cs_field_by_id(cm->immel)->val;

  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  /* Source term for gas enthalpy due to inter-phase fluxes */
  cs_real_t *smbrsh1 = cs_field_by_name("x_h_c_exp_st")->val;
  cs_real_t *rovsdth1 = cs_field_by_name("x_h_c_imp_st")->val;

  const cs_real_t *cvar_yno = NULL, *cvara_yno = NULL;
  const cs_real_t *cvara_yhcn = NULL, *cvara_ynh3 = NULL;

  if (cm->ieqnox == 1 && cs_glob_time_step->nt_cur > 1) {
    cvar_yno = cs_field_by_id(cm->iyno)->val;
    cvara_yno = cs_field_by_id(cm->iyno)->val_pre;
    cvara_yhcn = cs_field_by_id(cm->iyhcn)->val_pre;
    cvara_ynh3 = cs_field_by_id(cm->iynh3)->val_pre;
  }

  bool log_active = cs_log_default_is_active();

  const char log_st_fmt[] = N_("Source terms for variable %s\n\n");

  /* Aliases for simpler syntax */

  const cs_real_t *wmole = cm->wmole;
  const cs_real_t *wmolat = cm->wmolat;
  const int ngazem = CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  /* 2. Consideration of the source terms for variables
   *    relative to the classes of particles
   *--------------------------------------------------------------------------*/

  /* Source term for the mass fraction of reactive coal */

  if (fld_id >= cm->ixch[0] && fld_id <= cm->ixch[cm->nclacp -1]) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    int class_id = cs_field_get_key_int(fld_scal, keyccl) - 1;

    const cs_real_t *cpro_cgch = cs_field_by_id(cm->igmdch[class_id])->val;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      // W1 = - rho.GMDCH > 0
      cs_real_t xw1 = - crom[c_id]*cpro_cgch[c_id]*cell_f_vol[c_id];

      // Explicit and implicit parts of source terms
      rovsdt[c_id] += rovsdt[c_id] + cs_math_fmax(xw1, 0);
      smbrs[c_id]  -= xw1*cvara_var[c_id];
    }

  }

  /* 2.1 Source term for the mass fraction of coke
   *--------------------------------------------------------------------------*/

  if (fld_id >= cm->ixck[0] && fld_id <= cm->ixck[cm->nclacp -1]) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    int class_id = cs_field_get_key_int(fld_scal, keyccl) - 1;

    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cvara_xckcl = cs_field_by_id(cm->ixck[class_id])->val_pre;

    const cs_real_t *cpro_cgch = cs_field_by_id(cm->igmdch[class_id])->val;
    const cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
    const cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;
    const cs_real_t *cpro_cght = cs_field_by_id(cm->igmhet[class_id])->val;
    const cs_real_t *cpro_ghco2 = NULL, *cpro_ghh2o = NULL;

    if (cm->ihtco2 ==  1)
      cpro_ghco2 = cs_field_by_id(cm->ighco2[class_id])->val;
    if (cm->ihth2o ==  1)
      cpro_ghh2o = cs_field_by_id(cm->ighh2o[class_id])->val;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t exp_st = 0.;

      /* volatile formation minus coal consuming = char formation
       *  (Coke formation in French)
       *  NB: we take values at current and not previous time step
       *      to be conservative in mass */
      cs_real_t char_formation =   crom[c_id]*cvar_xchcl[c_id]*cell_f_vol[c_id]
                                 * (  cpro_cgd1[c_id] + cpro_cgd2[c_id]
                                    - cpro_cgch[c_id]);

      /* Compute the implict part of the Source term */
      if (cvara_xckcl[c_id] > cs_coal_epsilon) {

        // Reaction C(s) + O2 ---> 0.5CO
        exp_st = cpro_cght[c_id];

        // Reaction C(s) + CO2 ---> 2CO
        if (cm->ihtco2 == 1)
          exp_st += cpro_ghco2[c_id];

        // Reaction C(s) + H2O ---> CO + H2
        if (cm->ihth2o == 1)
          exp_st += cpro_ghh2o[c_id];

        exp_st = -2./3. * crom[c_id] * exp_st
                        * cell_f_vol[c_id] / pow(cvara_xckcl[c_id], (1./3.));
      }

      // Compute the explicit part of the Source term
      cs_real_t imp_st = -3./2. * exp_st * cvara_xckcl[c_id];

      rovsdt[c_id] += cs_math_fmax(exp_st, 0.);
      smbrs[c_id] += char_formation + imp_st;

    }

  }

  /* 2.2 Source term for the mass fraction of water
   *--------------------------------------------------------------------------*/

  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {

    if (fld_id >= cm->ixwt[0] && fld_id <= cm->ixwt[cm->nclacp -1]) {

      if (eqp->verbosity >= 1)
        bft_printf(_(log_st_fmt), fld_scal->name);

      int class_id = cs_field_get_key_int(fld_scal, keyccl) - 1;
      int coal_id = cm->ichcor[class_id] - 1;

      const cs_real_t xwatch = cm->xwatch[coal_id];
      const cs_real_t one_d_xwatch
        = (xwatch > cs_coal_epsilon) ? 1./cm->xwatch[coal_id] : 0;

      const cs_real_t *cpro_csec = cs_field_by_id(cm->igmsec[class_id])->val;
      const cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        // Calculation of explicit and implicit parts of source terms

        if (   cvara_var[c_id] > cs_coal_epsilon
            && xwatch > cs_coal_epsilon) {
          cs_real_t xw1 =  crom[c_id]*cpro_csec[c_id]*cell_f_vol[c_id]
                          *(1./cpro_x2[c_id])*one_d_xwatch;

          rovsdt[c_id] += cs_math_fmax(xw1, 0);
          smbrs[c_id]  -= xw1 * cvara_var[c_id];
        }
      }

    }

  }

  /* 2.3 Particle age source term
   *--------------------------------------------------------------------------*/

  if (cm->idrift >= 1) {

    if (strncmp(fld_scal->name, "n_p_age", 7) == 0) {

      // index of the coal particle class
      int class_id = cs_field_get_key_int(fld_scal, keyccl) - 1;

      const cs_real_t *cvar_xnpcl = cs_field_by_id(cm->inp[class_id])->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrs[c_id] += crom[c_id] * cell_f_vol[c_id]*cvar_xnpcl[c_id];
      }

    }

    else if (strncmp(fld_scal->name, "age", 3) == 0) {

      // Age of the bulk source term
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrs[c_id] += crom[c_id] * cell_f_vol[c_id];
      }

    }

  }

  /* 2.4 Particle velocity source terms
   *--------------------------------------------------------------------------*/

  if (cm->idrift == 1) {

    // index of the coal particle class
    int class_id = cs_field_get_key_int(fld_scal, keyccl) - 1;

    if (class_id >= 0) {
      const cs_real_t *cpro_cght = cs_field_by_id(cm->igmhet[class_id])->val;
      const cs_real_t *cpro_ghco2 = NULL, *cpro_ghh2o = NULL;

      if (cm->ihtco2 == 1)
        cpro_ghco2 = cs_field_by_id(cm->ighco2[class_id])->val;
      if (cm->ihth2o == 1) {
        cpro_ghh2o = cs_field_by_id(cm->ighh2o[class_id])->val;
      }

      const cs_real_t *cvara_coke = cs_field_by_id(cm->ixck[class_id])->val_pre;
      // Taup
      char name[25];
      snprintf(name, 24, "n_p_%02d", class_id+1);
      name[24] = '\0';

      cs_real_t *taup = cs_field_by_composite_name(name, "drift_tau")->val;

      const char *coo_prefix[3] = {"v_x_p_", "v_y_p_", "v_z_p_"};
      const char coo_c[3] = {'x', 'y', 'z'};

      for (cs_lnum_t coo_id = 0; coo_id < 3; coo_id++) {

        if (strncmp(fld_scal->name, coo_prefix[coo_id], 6) == 0) {

          snprintf(name, 24, "v_%c_p_%02d", coo_c[coo_id], class_id+1);
          name[24] = '\0';
          const cs_real_t *vp_xyz = cs_field_by_name(name)->val;

          snprintf(name, 24, "vg_lim_p_%02d", class_id+1);
          name[24] = '\0';
          const cs_real_3_t *vg_lim_pi
            = (const cs_real_3_t *)cs_field_by_name(name)->val;

          // Deviation velocity for the continuous phase
          const cs_real_3_t *vdc
            = (const cs_real_3_t *)cs_field_by_name("vd_c")->val;

          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

            cs_real_t smbrs1 = 0.;

            /* Drying and Devolatilization have no effect on Vp
             * (mass flux is only exiting the particle)
             * During Heterogeneous reactions:
             *  - Gas molecules, with the velocity Vc are are absorbed.
             *  - The product, with the velocity Vp, is released.
             * For Oxygen oxydation one O2 comes in and two CO go out
             * For CO2 gasification, one  CO2 comes in and two CO go out
             * For H2O gasification, one H2O comes in and one CO and one H2
             * goes out */
            if (cvara_coke[c_id] > cs_coal_epsilon) {
              smbrs1 += wmole[io2]/wmolat[cs_coal_atom_id_c]*cpro_cght[c_id];
              if (cm->ihtco2 == 1)
                smbrs1 +=   wmole[ico2]/cm->wmolat[cs_coal_atom_id_c]
                          * cpro_ghco2[c_id];
              if (cm->ihth2o == 1)
                smbrs1 += wmole[ih2o]/wmolat[cs_coal_atom_id_c]*cpro_ghh2o[c_id];
              smbrs1 *= pow(cvara_coke[c_id], 2./3.);
            }

            // relaxation to drop velocity
            smbrs1 =  crom[c_id]*cell_f_vol[c_id]*(1./taup[c_id]+smbrs1)
                     *(  vel[c_id][coo_id]+vdc[c_id][coo_id]
                       + vg_lim_pi[c_id][coo_id]-vp_xyz[c_id]);

            smbrs[c_id] += smbrs1;
            rovsdt[c_id] += crom[c_id]*cell_f_vol[c_id]/taup[c_id];

          } /* loop on cells */

        } /* name starts with one of "v_x_p_", "v_y_p_", "v_z_p_" */

      } /* loop on coordinate id */

    } /* test on class id */

  } /* cm->idrift == 1 */

  /* Source term for the enthalpies
   *--------------------------------------------------------------------------*/

  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells, cs_real_t);

  if (fld_id >= cm->ih2[0] && fld_id <= cm->ih2[cm->nclacp -1]) {

    // Initialization of the exchange terms for gas enthalpy
    if (fld_id == cm->ih2[0]) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrsh1[c_id] = 0.;
        rovsdth1[c_id] = 0.;
      }
    }

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    int class_id = cs_field_get_key_int(fld_scal, keyccl) - 1;
    int coal_id = cm->ichcor[class_id] - 1;

    const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;
    const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
    const cs_real_t *cvara_xckcl = cs_field_by_id(cm->ixck[class_id])->val_pre;

    const cs_real_t *cvar_xwtcl = NULL, *cvara_xwtcl = NULL;

    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
      cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;
      cvara_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val_pre;
    }

    const cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;
    const cs_real_t *cpro_rom2 = cs_field_by_id(cm->irom2[class_id])->val;
    const cs_real_t *cpro_diam2 = cs_field_by_id(cm->idiam2[class_id])->val;
    const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;
    const cs_real_t *cpro_cght = cs_field_by_id(cm->igmhet[class_id])->val;
    const cs_real_t *cpro_ghco2 = NULL, *cpro_ghh2o = NULL;

    if (cm->ihtco2 == 1)
      cpro_ghco2 = cs_field_by_id(cm->ighco2[class_id])->val;
    if (cm->ihth2o == 1)
      cpro_ghh2o = cs_field_by_id(cm->ighh2o[class_id])->val;

    const cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
    const cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;

    /* Contribution to the explicit and implicit balance
     * exchanges by molecular distribution
     * 6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2) */

    // Calculation of lambda in W1.

    const cs_real_t xnuss = 2.;

    const cs_field_t *fld_th = cs_thermal_model_field();
    const int kivisl  = cs_field_key_id("diffusivity_id");
    const int kvisl0 = cs_field_key_id("diffusivity_ref");
    const cs_real_t visls_0 = cs_field_get_key_double(fld_th, kvisl0);
    const cs_real_t cp0 = cs_glob_fluid_properties->cp0;

    const int ifcvsl = cs_field_get_key_int(fld_th, kivisl);
    const cs_real_t *cpro_viscls = NULL;
    if (ifcvsl >= 0)
      cpro_viscls = cs_field_by_id(ifcvsl)->val;

    if (cpro_viscls != NULL) {
      if (cpro_cp != NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          w1[c_id] = cpro_viscls[c_id] * cpro_cp[c_id];
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          w1[c_id] = cpro_viscls[c_id] * cp0;
      }
    }
    else {
      if (cpro_cp != NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          w1[c_id] = visls_0 * cpro_cp[c_id];
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          w1[c_id] = visls_0 * cp0;
      }
    }

    /* Contribution to the explicit and implicit balance
     * exchanges by molecular distribution
     * Remark: We use cpro_x2[c_id] because we want X2 at the iteration n */

    cs_real_t *cpro_rovsdt2 = cs_field_by_id(cm->igmtr[class_id])->val;

    const cs_real_t xashch = cm->xashch[coal_id];
    const cs_real_t diam20 = cm->diam20[class_id];
    const cs_real_t diam20_sq = cs_math_pow2(diam20);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      /* Calculation of diameter of the particles
       * d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5 */
      cs_real_t diam2 =   xashch*diam20_sq
                        + (1.-xashch)*cs_math_pow2(cpro_diam2[c_id]);

      cs_real_t aux =   6. * w1[c_id] * xnuss / diam2 / cpro_rom2[c_id]
                      * crom[c_id] * cell_f_vol[c_id];

      smbrs[c_id] += -aux*(cpro_temp2[c_id]-cpro_temp[c_id])*cpro_x2[c_id];
      smbrsh1[c_id] += aux*(cpro_temp2[c_id]-cpro_temp[c_id])*cpro_x2[c_id];

      /* Store the implicit part of the exchange so that we can compute a
       * conservative exhcange term when computing the gas enthalpy */
      cpro_rovsdt2[c_id] = aux / cm->cp2ch[coal_id];
      rovsdt[c_id] += cpro_rovsdt2[c_id];

    }

    /* Contribution to the explicit and implicit balances
     * of exchange term of energy between phases:
     * gama(dev1) H(mv1,T2)+gama(dev2) H(mv2,T2) */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
      cs_real_t f1mc[CS_COMBUSTION_MAX_COALS], f2mc[CS_COMBUSTION_MAX_COALS];

      /* Gama Dev1 et Gama Dev2 */

      cs_real_t gamdv1 = crom[c_id]*cvar_xchcl[c_id]*cpro_cgd1[c_id];
      cs_real_t gamdv2 = crom[c_id]*cvar_xchcl[c_id]*cpro_cgd2[c_id];

      /* H(mv1, T2) */

      for (int ige = 0; ige < ngazem; ige++)
        coefe[ige] = 0.;

      const cs_real_t t2 = cpro_temp2[c_id];
      cs_real_t den;

      int ichx1c_coal = cm->ichx1c[coal_id] - 1;

      den =   cm->a1[coal_id]*wmole[ichx1c_coal]
            + cm->b1[coal_id]*wmole[ico]
            + cm->c1[coal_id]*wmole[ih2o]
            + cm->d1[coal_id]*wmole[ih2s]
            + cm->e1[coal_id]*wmole[ihcn]
            + cm->f1[coal_id]*wmole[inh3];

      coefe[ichx1] = cm->a1[coal_id]*wmole[ichx1c_coal] / den;
      coefe[ico]   = cm->b1[coal_id]*wmole[ico]         / den;
      coefe[ih2o]  = cm->c1[coal_id]*wmole[ih2o]        / den;
      coefe[ih2s]  = cm->d1[coal_id]*wmole[ih2s]        / den;
      coefe[ihcn]  = cm->e1[coal_id]*wmole[ihcn]        / den;
      coefe[inh3]  = cm->f1[coal_id]*wmole[inh3]        / den;

      for (int icha = 0; icha < CS_COMBUSTION_MAX_COALS; icha++) {
        f1mc[icha] = 0;
        f2mc[icha] = 0;
      }
      f1mc[coal_id] = 1.;

      cs_real_t xhdev1
        = cs_coal_ht_convert_t_to_h_gas_by_yi_f1f2(t2, coefe, f1mc, f2mc);

      /* H(mv2, T2) */

      for (int ige = 0; ige < ngazem; ige++)
        coefe[ige] = 0.;

      int ichx2c_coal = cm->ichx2c[coal_id] - 1;

      den =   cm->a2[coal_id]*wmole[ichx2c_coal]
            + cm->b2[coal_id]*wmole[ico]
            + cm->c2[coal_id]*wmole[ih2o]
            + cm->d2[coal_id]*wmole[ih2s]
            + cm->e2[coal_id]*wmole[ihcn]
            + cm->f2[coal_id]*wmole[inh3];

      coefe[ichx2] = cm->a2[coal_id]*wmole[ichx2c_coal] / den;
      coefe[ico]   = cm->b2[coal_id]*wmole[ico]         / den;
      coefe[ih2o]  = cm->c2[coal_id]*wmole[ih2o]        / den;
      coefe[ih2s]  = cm->d2[coal_id]*wmole[ih2s]        / den;
      coefe[ihcn]  = cm->e2[coal_id]*wmole[ihcn]        / den;
      coefe[inh3]  = cm->f2[coal_id]*wmole[inh3]        / den;

      for (int icha = 0; icha < CS_COMBUSTION_MAX_COALS; icha++) {
        f1mc[icha] = 0;
        f2mc[icha] = 0;
      }
      f2mc[coal_id] = 1.;

      cs_real_t xhdev2
        = cs_coal_ht_convert_t_to_h_gas_by_yi_f1f2(t2, coefe, f1mc, f2mc);

      /*  Contribution to explicit and implicit balances */

      smbrs[c_id]   +=  (gamdv1*xhdev1+gamdv2*xhdev2)*cell_f_vol[c_id];
      smbrsh1[c_id] -=  (gamdv1*xhdev1+gamdv2*xhdev2)*cell_f_vol[c_id];
    }

    /* Heterogeneous combustion: C(s) + 02 ---> 0.5 C0
     * GamHET * (28/12 H(CO,T2)-16/12 H(O2,T1)) */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

      /* Calculation of HCO(T2) */

      for (int ige = 0; ige < ngazem; ige++)
        coefe[ige] = 0.;
      coefe[ico] = 1.;

      cs_real_t t2 = cpro_temp2[c_id];
      cs_real_t xhco
        = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

      /* Calculation of HO2(T1) */

      for (int ige = 0; ige < ngazem; ige++)
        coefe[ige] = 0.;
      coefe[io2] = 1.;

      cs_real_t t1 = cpro_temp[c_id];
      cs_real_t xho2
        = cs_coal_ht_convert_t_to_h_gas_by_yi(t1, coefe);

      /* Contribution to explicit and implicit balances */

      cs_real_t gamhet = 0;
      if (cvara_xckcl[c_id] > cs_coal_epsilon) {
        gamhet =   crom[c_id]*cpro_cght[c_id]
                 * (          pow(cvara_xckcl[c_id], 2./3.)
                    + 2./3. * (cvar_xckcl[c_id]-cvara_xckcl[c_id])
                            / pow(cvara_xckcl[c_id], 1./3.));
      }

      gamhet *=   (wmole[ico]*xhco - wmolat[cs_coal_atom_id_o]*xho2)
                / wmolat[cs_coal_atom_id_c]
                * cell_f_vol[c_id];

      smbrs[c_id] += gamhet;
      smbrsh1[c_id] -= gamhet;
    }

    /* Heterogeneous combustion: C(s) + C02 ---> 2 C0
     * GamHET * (56/12 H(CO,T2)-44/12 H(CO2,T1)) */

    if (cm->ihtco2 == 1) {

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        /* Calculation of HCO(T2) */

        cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ico] = 1.;

        cs_real_t t2 = cpro_temp2[c_id];
        cs_real_t xhco = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

        /* Calculation of HCO2(T1) */

        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ico2] = 1.;

        cs_real_t t1 = cpro_temp[c_id];
        cs_real_t xhco2 = cs_coal_ht_convert_t_to_h_gas_by_yi(t1, coefe);

        /* Contribution to explicit and implicit balances */

        cs_real_t gamhet = 0;
        if (cvara_xckcl[c_id] > cs_coal_epsilon) {
          gamhet =   crom[c_id]*cpro_ghco2[c_id]
                 * (          pow(cvara_xckcl[c_id], 2./3.)
                    + 2./3. * (cvar_xckcl[c_id]-cvara_xckcl[c_id])
                            / pow(cvara_xckcl[c_id], 1./3.));
        }

        gamhet *=   (2.*wmole[ico]*xhco-wmole[ico2]*xhco2)
                  / wmolat[cs_coal_atom_id_c]
                  * cell_f_vol[c_id];

        smbrs[c_id] += gamhet;
        smbrsh1[c_id] -= gamhet;
      }

    }

    /* Heterogeneous combustion: C(s) + H2O ---> CO + H2
     * GamHET * (28/12 H(CO,T2)+2/12 H(HY,T2) -18/12 H(H2O,T1)) */

    if (cm->ihth2o == 1) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        const cs_real_t t2 = cpro_temp2[c_id];

        /* Calculation of HCO(T2) */

        cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ico] = 1.;

        cs_real_t xhco = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

        /* Calculation of HH2(T2) */

        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ihy] = 1.;

        cs_real_t xhh2= cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

        /* Calculation of HH2O(T1) */

        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ih2o] = 1.;

        cs_real_t t1 = cpro_temp[c_id];
        cs_real_t xhh2o = cs_coal_ht_convert_t_to_h_gas_by_yi(t1, coefe);

        /* Contribution to explicit and implicit balances */

        cs_real_t gamhet = 0;
        if (cvara_xckcl[c_id] >cs_coal_epsilon) {
          gamhet =   crom[c_id]*cpro_ghh2o[c_id]
                 * (          pow(cvara_xckcl[c_id], 2./3.)
                    + 2./3. * (cvar_xckcl[c_id]-cvara_xckcl[c_id])
                            / pow(cvara_xckcl[c_id], 1./3.));
        }

        gamhet *=  (wmole[ico]*xhco+wmole[ihy]*xhh2 - wmole[ih2o]*xhh2o)
                  / wmolat[cs_coal_atom_id_c] * cell_f_vol[c_id];

        smbrs[c_id] += gamhet;
        smbrsh1[c_id] -= gamhet;
      }

    }

    /* Source term on H2 (coming from drying) */

    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {

      /* Contribution of interfacial source term interfacial to balances */

      const cs_real_t *cpro_csec = cs_field_by_id(cm->igmsec[class_id])->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        /* Calculation of H(H2O) at T2 */

        cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ih2o] = 1.;

        cs_real_t t2 = cpro_temp2[c_id];
        if (t2 > 100.+tkelvi) {
          t2 = 100.+tkelvi;
        }

        cs_real_t hh2ov
          = cs_coal_ht_convert_t_to_h_gas_by_yi_with_drying(t2, coefe);

        /* Contribution to explicit balance */

        cs_real_t aux = 0;
        if (   cvara_xwtcl[c_id]   > cs_coal_epsilon
            && cm->xwatch[coal_id] > cs_coal_epsilon) {
          aux = -   crom[c_id] * cpro_csec[c_id]
                  * (cvar_xwtcl[c_id] / cpro_x2[c_id])
                  * (1. / cm->xwatch[coal_id])*hh2ov;
        }

        smbrs[c_id] += aux * cell_f_vol[c_id];
        smbrsh1[c_id] -= aux * cell_f_vol[c_id];

      }

    }  /* drying */

  } /* Enthalpies of particles */

  /* 3. Taking into account source terms for relative variables in the mixture
   *--------------------------------------------------------------------------*/

  /* source terms from particles (convection and enthalpy
     driven by mass fluxes) */

  if (fld_id == cm->ihgas) {

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] += smbrsh1[c_id];
      rovsdt[c_id] += rovsdth1[c_id];
    }

    /* Explicit contribution due to implicit source term on particle
       * class enthalpy */
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cpro_rovsdt2 = cs_field_by_id(cm->igmtr[class_id])->val;
      const cs_real_t *cvar_x2h2 = cs_field_by_id(cm->ih2[class_id])->val;
      const cs_real_t *cvara_x2h2 = cs_field_by_id(cm->ih2[class_id])->val_pre;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrs[c_id] += cpro_rovsdt2[c_id]*(cvar_x2h2[c_id] - cvara_x2h2[c_id]);
      }
    }

  }

  /* Source term for light volatile materials */

  if (fld_id >= cm->if1m[0] && fld_id <= cm->if1m[cm->n_coals - 1]) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    /* Calculation of GMDEV1 = - Sum (rho.XCH.GMDV1) > 0  --> W1 */

    int numcha = fld_id - cm->if1m[0] + 1;

    cs_array_real_fill_zero(n_cells, w1);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cpro_cgd1 = cs_field_by_id(cm->igmdv1[class_id])->val;
      const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;

      if (cm->ichcor[class_id] == numcha) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          w1[c_id] -= crom[c_id]*cvar_xchcl[c_id] * cpro_cgd1[c_id];
        }
      }
    }

    /* Contribution of interfacial source term to balances */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] += cell_f_vol[c_id] * w1[c_id];
    }

  }

  /* Source terms for heavy volatile materials */

  if (fld_id >= cm->if2m[0] && fld_id <= cm->if2m[cm->n_coals - 1]) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    /* Calculation of GMDEV2 = - Sum (rho.XCH.GMDV2) > 0 --> W1 */

    int numcha = fld_id - cm->if2m[0] + 1;
    cs_array_real_fill_zero(n_cells, w1);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cpro_cgd2 = cs_field_by_id(cm->igmdv2[class_id])->val;
      const cs_real_t *cvar_xchcl = cs_field_by_id(cm->ixch[class_id])->val;

      if (cm->ichcor[class_id] == numcha) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          w1[c_id] -= crom[c_id]*cvar_xchcl[c_id] * cpro_cgd2[c_id];
        }
      }
    }

    /* Contribution of interfacial source term to explicit balance */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] += cell_f_vol[c_id] * w1[c_id];
    }

  }

  /* Source term for the tracer 7 (O2) (heterogeneous combustion by C) */

  if (fld_id == cm->if7m) {

    /* Remark: We take the same source term as for Xck to be conservative */

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    cs_array_real_fill_zero(n_cells, w1);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cpro_cght = cs_field_by_id(cm->igmhet[class_id])->val;
      const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
      const cs_real_t *cvara_xckcl = cs_field_by_id(cm->ixck[class_id])->val_pre;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (cvara_xckcl[c_id] > cs_coal_epsilon) {
          w1[c_id] -=  crom[c_id]*cpro_cght[c_id]
                       * (   pow(cvara_xckcl[c_id], 2./3.)
                          + 2./3.*(cvar_xckcl[c_id] - cvara_xckcl[c_id])
                                 / pow(cvara_xckcl[c_id], 1./3.));
        }
      }
    }

    /* Contribution of interfacial source term to balances */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] += cell_f_vol[c_id] * w1[c_id];
    }

  }

  /* Source term for the tracer 8 (CO2) (heterogeneous combustion by C) */

  if (cm->ihtco2 == 1 && fld_id == cm->if8m) {

    /* Remark: We take the same source term as for Xck to be conservative */

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    cs_array_real_fill_zero(n_cells, w1);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cpro_ghco2 = cs_field_by_id(cm->ighco2[class_id])->val;
      const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
      const cs_real_t *cvara_xckcl = cs_field_by_id(cm->ixck[class_id])->val_pre;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (cvara_xckcl[c_id] > cs_coal_epsilon) {
          w1[c_id] -=  crom[c_id]*cpro_ghco2[c_id]
                       * (   pow(cvara_xckcl[c_id], 2./3.)
                          + 2./3.*(cvar_xckcl[c_id] - cvara_xckcl[c_id])
                                 / pow(cvara_xckcl[c_id], 1./3.));
        }
      }
    }

    /* Contribution of interfacial source term to balances */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] += cell_f_vol[c_id] * w1[c_id];
    }

  }

  /* Source term for the tracer 9 (H2O) (heterogeneous combustion by H2O) */

  if (cm->ihth2o == 1 && fld_id == cm->if9m) {

    /* Remark: We take the same source term as for Xck to be conservative */

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    cs_array_real_fill_zero(n_cells, w1);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cpro_ghh2o = cs_field_by_id(cm->ighh2o[class_id])->val;
      const cs_real_t *cvar_xckcl = cs_field_by_id(cm->ixck[class_id])->val;
      const cs_real_t *cvara_xckcl = cs_field_by_id(cm->ixck[class_id])->val_pre;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (cvara_xckcl[c_id] > cs_coal_epsilon) {
          w1[c_id] -=  crom[c_id]*cpro_ghh2o[c_id]
                       * (   pow(cvara_xckcl[c_id], 2./3.)
                          + 2./3.*(cvar_xckcl[c_id] - cvara_xckcl[c_id])
                                 / pow(cvara_xckcl[c_id], 1./3.));
        }
      }
    }

    /* Contribution of interfacial source term to balances */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] += cell_f_vol[c_id] * w1[c_id];
    }

  }

  /* Source term for the fuel variance */

  if (fld_id == cm->ifvp2m) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    _coal_fp2st(fld_scal, smbrs, rovsdt);

  }

  /* Source term for the tracer 6 (Water coming from drying) */

  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING && fld_id == cm->if6m) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    /* Contribution of interfacial source term to balances */

    cs_array_real_fill_zero(n_cells, w1);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      const cs_real_t *cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;
      const cs_real_t *cvara_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val_pre;
      const cs_real_t *cpro_csec = cs_field_by_id(cm->igmsec[class_id])->val;
      const cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;
      const int coal_id = cm->ichcor[class_id] -1;
      const cs_real_t xwat_coal = cm->xwatch[coal_id];

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (   cvara_xwtcl[c_id] > cs_coal_epsilon
            && xwat_coal >  cs_coal_epsilon) {
          w1[c_id] +=  crom[c_id]*cpro_csec[c_id]
                      * (cvar_xwtcl[c_id] / cpro_x2[c_id])
                      * (1. / xwat_coal);
        }
      }

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrs[c_id] += cell_f_vol[c_id] * w1[c_id];
      }

    }

  }

  /* Source term for CO2 */

  if (cm->ieqco2 == 1 && fld_id == cm->iyco2) {

    if (eqp->verbosity >= 1)
      bft_printf(_(log_st_fmt), fld_scal->name);

    // Use "field_by_name" here rather than field pointers so
    // as to have a readable error message if used with an incompatible
    // turbulence model.
    const cs_real_t *cvara_k = cs_field_by_name("k")->val_pre;
    const cs_real_t *cvara_ep = cs_field_by_name("epsilon")->val_pre;

    cs_real_t **cpro_x2c;
    BFT_MALLOC(cpro_x2c, cm->nclacp, cs_real_t *);
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cpro_x2c[class_id] = cs_field_by_id(cm->ix2[class_id])->val;
    }

    /* Contribution of interfacial source term to balances */

    // Oxydation of CO
    // ---------------

    //  Dryer Glassman : XK0P in (mol/m3)**(-0.75) s-1
    //          xk0p = 1.26e10
    //          xk0p = 1.26e7 * (1.1)**(nt_cur)
    //          if (xk0p > 1.26e10) xk0p=1.26e10
    //          t0p  = 4807.
    //  Howard : xk0p in (mol/m3)**(-0.75) s-1
    //           xk0p = 4.11e9
    //           t0p  = 15090.
    //  Westbrook & Dryer

    const cs_real_t lnk0p = 23.256e0;
    const cs_real_t t0p  = 20096.;

    //  Hawkin and Smith Purdue University Engeneering Bulletin, i
    //  Research series 108 vol 33, n 3n 1949
    //  Kp = 10**(4.6-14833/T)
    //  Equilibrum constant in partial pressure [atm]
    //  XKOE is the decimal log of the pre-exponential constant
    //  TOE is NOT an activation temperature ... there is a lg(e)
    //  to return to Kc and to use concentrations (in mol/m3)
    //  Kc = (1/RT)**variation nb moles * Kp
    //  here Kc = sqrt(0.082*T)*Kp

    const cs_real_t l10k0e = 4.6e0;
    const cs_real_t t0e  = 14833.;

    // Dissociation of CO2 (Trinh Minh Chinh)
    // -------------------
    //          XK0M = 5.D8
    //          T0M  = 4807.
    //          XK0M = 0.
    //  Westbrook & Dryer

    const cs_real_t lnk0m = 20.03e0;
    const cs_real_t t0m  = 20096.;

    cs_real_t err1mx = 0.;

    // Number of iterations
    const int itermx = 500;
    // Number of convergent points

    int nberic = 0, nbpass = 0, nbarre = 0, nbimax = 0;

    // Precision on the convergence
    const cs_real_t errch = 1.e-8;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t xxco  = cpro_yco[c_id]/wmole[ico] * cpro_rom1[c_id];
      cs_real_t xxo2  = cpro_yo2[c_id]/wmole[io2] * cpro_rom1[c_id];
      cs_real_t xxco2 = cpro_yco2[c_id]/wmole[ico2] * cpro_rom1[c_id];
      cs_real_t xxh2o = cpro_yh2o[c_id]/wmole[ih2o] * cpro_rom1[c_id];

      xxco  = cs_math_fmax(xxco, 0);
      xxo2  = cs_math_fmax(xxo2, 0);
      xxco2 = cs_math_fmax(xxco2, 0);
      xxh2o = cs_math_fmax(xxh2o, 0);

      cs_real_t sqh2o = sqrt(xxh2o);

      cs_real_t xkp = exp(lnk0p - t0p/cpro_temp[c_id]);
      cs_real_t xkm = exp(lnk0m - t0m/cpro_temp[c_id]);

      cs_real_t xkpequ = pow(10., l10k0e - t0e/cpro_temp[c_id]);
      cs_real_t xkcequ = xkpequ / sqrt(8.32 * cpro_temp[c_id] / 1.015e5);

      // initialization by the transported state

      // cs_real_t anmr  = xxco2;
      cs_real_t xcom  = xxco + xxco2;
      cs_real_t xo2m  = xxo2 + 0.5*xxco2;

      // cs_real_t xo2eq  = 0.;
      // cs_real_t xcoeq  = 0;
      cs_real_t xco2eq = 0.;

      if (cpro_temp[c_id] > 1200.) {

        // Search for the equilibrum state.
        // Iterative search without control of convergence
        //  (to keep the parallelisation on meshes).
        // On the numver of moles of separating reaction
        //  the state before reaction (such as calculated by Cpcym)
        //  of the equilibrum state.
        //  anmr has to be the boundary between 0 and Min(XCOM,2.*XO2M)
        //  We look for the solution by bisection.

        cs_real_t anmr0 = 0.;
        cs_real_t anmr1 = cs_math_fmin(xcom, 2.*xo2m);
        cs_real_t anmr2 = 0.5*(anmr0+anmr1);
        int iterch = 0;
        cs_real_t fn2 = 1.;
        cs_real_t fn0 = -0.5                                     * pow(anmr0, 3)
                       + (     xcom     + xo2m - pow(xkcequ, 2)) * pow(anmr0, 2)
                       - (  .5*xcom  + 2.*xo2m)*xcom             * anmr0
                       +   pow(xcom, 2) * xo2m;
        cs_real_t fn1 = -0.5                                     * pow(anmr1, 3)
                       + (     xcom     + xo2m - pow(xkcequ, 2)) * pow(anmr1, 2)
                       - (  .5*xcom  + 2.*xo2m)*xcom             * anmr1
                       +   pow(xcom, 2) * xo2m;

        if (xo2m > 1.e-6) {
          while (iterch < itermx && fn2 > errch) {
            anmr2 = 0.5*(anmr0+anmr1);
            fn2 = -0.5                                   * pow(anmr2, 3)
                   + (   xcom   + xo2m - pow(xkcequ, 2)) * pow(anmr2, 2)
                   - (.5*xcom +2.*xo2m)*xcom             * anmr2
                   + pow(xcom, 2) * xo2m;
            if (fn0*fn2 > 0.) {
              anmr0 = anmr2;
              fn0 = fn2;
            }
            else if (fn1*fn2 > 0.) {
              anmr1 = anmr2;
              fn1 = fn2;
            }
            else if (fn0*fn1 > 0.) {
              iterch = itermx;
              anmr2 = cs_math_fmin(xcom, 2.*xo2m);
              nbarre = nbarre + 1;
            }
            iterch = iterch + 1;
          } /* end while */

          if (iterch >= itermx)
            nberic = nberic + 1;
          else
            nbimax = CS_MAX(nbimax, iterch);
          err1mx = cs_math_fmax(err1mx, fn2);

          xco2eq = anmr2;
          // xcoeq  = xcom - anmr2;
          // xo2eq  = xo2m - 0.5 * anmr2;
        }
        else {
          // xo2eq  = 0.;
          // xcoeq  = xxco;
          xco2eq = 0.;
        }
      }

      else {
        xco2eq = cs_math_fmin(xcom, 2.*xo2m);
        // xo2eq  = xo2m - 0.5*xco2eq;
        // xcoeq  = xcom - xco2eq;
      }

      cs_real_t xden = 0;

      if (xco2eq > xxco2) {
        // oxydation
        xden = xkp * sqh2o * pow(xxo2, 0.25);
      }
      else {
        // dissociation
        xden = xkm;
      }
      if (cs_math_fabs(xden) > 0.) {

        cs_real_t tauchi = 1./xden;
        cs_real_t tautur = cvara_k[c_id]/cvara_ep[c_id];

        cs_real_t x2 = 0.;
        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          x2 += cpro_x2c[class_id][c_id];
        }

        if (cm->ieqco2 == 1) {
          // We transport CO2
          smbrs[c_id] +=   wmole[ico2] / cpro_rom1[c_id]
                         * (xco2eq-xxco2)/(tauchi+tautur) * (1.-x2)
                         * cell_f_vol[c_id] * crom[c_id];
        }
        else if (cm->ieqco2 == 2) {
          // We transport CO
          smbrs[c_id] +=   wmole[ico] / cpro_rom1[c_id]
                         * (xco2eq-xxco)/(tauchi+tautur) * (1.-x2)
                         * cell_f_vol[c_id] * crom[c_id];
        }

        w1[c_id] = cell_f_vol[c_id]*crom[c_id]/(tauchi+tautur);
        rovsdt[c_id] += cs_math_fmax(w1[c_id], 0);

      }

    } /* loop on cells */

    if (log_active) {
      cs_gnum_t cpt[] = {nberic, nbarre, nbpass};
      cs_parall_counter(cpt, 3);

      cs_parall_max(1, CS_INT_TYPE, &nbimax);
      cs_parall_max(1, CS_REAL_TYPE, &err1mx);

      cs_log_printf(CS_LOG_DEFAULT,
                    _(" Max Error = %g\n"
                      " no Points   %llu %llu %llu\n"
                      " Iter max number %d\n"),
                    err1mx, (unsigned long long)cpt[0],
                    (unsigned long long)cpt[1], (unsigned long long)cpt[2],
                    nbimax);
    }

    BFT_FREE(cpro_x2c);

    /* Source term: heterogeneous combustion by CO2 */

    if (cm->ihtco2 == 1) {

      // Arrays of pointers containing the fields values for each class
      // (loop on cells outside loop on classes)

      cs_real_t **cvara_xck, **cpro_ghco2a;
      BFT_MALLOC(cvara_xck, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cpro_ghco2a, cm->nclacp, cs_real_t *);
      for (int class_id = 0; class_id < cm->nclacp; class_id++) {
        cvara_xck[class_id] = cs_field_by_id(cm->ixck[class_id])->val_pre;
        cpro_ghco2a[class_id] = cs_field_by_id(cm->ighco2[class_id])->val;
      }

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t aux = 0.;
        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          aux +=  crom[c_id] * cpro_ghco2a[class_id][c_id]
                * pow(cvara_xck[class_id][c_id], 2./3.) * cell_f_vol[c_id];

        }

        rovsdt[c_id] -= aux*(wmole[ico2]/0.012);
      }

      BFT_FREE(cvara_xck);
      BFT_FREE(cpro_ghco2a);

    }

  }

  /* Source term for Enth_Ox
   *                 HCN and NO: only from the second iteration */

  if (cm->ieqnox== 1 && cs_glob_time_step->nt_cur >  1) {

    /* Terms on Oxydant enthalpy */

    if (fld_id == cm->ihox) {

      /* Arrays of pointers to field value arrays for each class */

      cs_real_t **cvar_xck = NULL, **cvara_xck = NULL, **cvar_xch = NULL;
      cs_real_t **cvar_xnp = NULL, **cpro_t2 = NULL, **cpro_gmhet = NULL;
      cs_real_t **cpro_ghco2a = NULL, **cpro_ghh2oa = NULL, **cvar_xwt = NULL;

      BFT_MALLOC(cvar_xck, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cvara_xck, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cvar_xch, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cvar_xnp, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cpro_t2, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cpro_gmhet, cm->nclacp, cs_real_t *);

      if (cm->ihtco2 == 1)
        BFT_MALLOC(cpro_ghco2a, cm->nclacp, cs_real_t *);

      if (cm->ihth2o == 1)
        BFT_MALLOC(cpro_ghh2oa, cm->nclacp, cs_real_t *);

      if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
        BFT_MALLOC(cvar_xwt, cm->nclacp, cs_real_t *);

      for (int class_id = 0; class_id < cm->nclacp; class_id++) {
        cvar_xck[class_id] = cs_field_by_id(cm->ixck[class_id])->val;
        cvara_xck[class_id] = cs_field_by_id(cm->ixck[class_id])->val_pre;
        cvar_xch[class_id] = cs_field_by_id(cm->ixch[class_id])->val;
        cvar_xnp[class_id] = cs_field_by_id(cm->inp[class_id])->val;
        cpro_t2[class_id] = cs_field_by_id(cm->itemp2[class_id])->val;
        cpro_gmhet[class_id] = cs_field_by_id(cm->igmhet[class_id])->val;
        if (cpro_ghco2a != NULL)
          cpro_ghco2a[class_id] = cs_field_by_id(cm->ighco2[class_id])->val;
        if (cpro_ghh2oa != NULL)
          cpro_ghh2oa[class_id] = cs_field_by_id(cm->ighh2o[class_id])->val;
        if (cvar_xwt != NULL)
          cvar_xwt[class_id] = cs_field_by_id(cm->ixwt[class_id])->val;
      }

      /* Calculation of T2 average on particles */

      cs_real_t tfuelmin =  HUGE_VAL;
      cs_real_t tfuelmax = -HUGE_VAL;

      cs_real_t *tfuel;
      BFT_MALLOC(tfuel, n_cells, cs_real_t);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        cs_real_t xmx2 = 0.;
        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          xmx2 +=   cvar_xck[class_id][c_id]
                  + cvar_xch[class_id][c_id]
                  + cvar_xnp[class_id][c_id]*cm->xmash[class_id];
          if (cvar_xwt != NULL)
            xmx2 += cvar_xwt[class_id][c_id];
        }

        if (xmx2 > 0.) {
          tfuel[c_id] = 0.;
          for (int class_id = 0; class_id < cm->nclacp; class_id++) {
            tfuel[c_id] += (  cvar_xck[class_id][c_id]
                            + cvar_xch[class_id][c_id]
                            + cvar_xnp[class_id][c_id]*cm->xmash[class_id])
                          * cpro_t2[class_id][c_id];
            if (cvar_xwt != NULL)
              tfuel[c_id] += cvar_xwt[class_id][c_id] * cpro_t2[class_id][c_id];
          }

          tfuel[c_id] /= xmx2;
        }
        else {
          tfuel[c_id] = cpro_temp[c_id];
        }

        tfuelmin = cs_math_fmin(tfuel[c_id], tfuelmin);
        tfuelmax = cs_math_fmax(tfuel[c_id], tfuelmax);
      }

      if (log_active) {
        cs_real_t buf[2] = {-tfuelmin, tfuelmax};
        cs_parall_max(2, CS_REAL_TYPE, buf);
        tfuelmin = -buf[0];
        tfuelmax = buf[1];
        cs_log_printf(CS_LOG_DEFAULT,
                    _(" Min max of Tfuel for Hoxy %g %g\n"),
                      tfuelmin, tfuelmax);
      }

      /* Heterogeneous combustion: C + O2 ---> 0.5 CO */

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        // Calculation of HCO(T2)

        cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[ico] = 1.;

        cs_real_t t2 = tfuel[c_id];
        cs_real_t xhco = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

        //  Calculation of HO2(T1)

        for (int ige = 0; ige < ngazem; ige++)
          coefe[ige] = 0.;
        coefe[io2] = 1.;
        cs_real_t t1 = cpro_temp[c_id];
        cs_real_t xho2 = cs_coal_ht_convert_t_to_h_gas_by_yi(t1, coefe);

        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          cs_real_t gamhet;
          if (cvara_xck[class_id][c_id] > cs_coal_epsilon) {
            gamhet =   crom[c_id]*cpro_gmhet[class_id][c_id]
                     * (  pow(cvara_xck[class_id][c_id], 2./3.)
                        + 2./3. * (  cvar_xck[class_id][c_id]
                                   - cvara_xck[class_id][c_id])
                                / pow(cvara_xck[class_id][c_id], 1./3.));
            smbrs[c_id] -= gamhet *(  wmole[ico]*xhco
                                    - wmolat[cs_coal_atom_id_o]*xho2)
                                  / wmolat[cs_coal_atom_id_c]
                                  * cell_f_vol[c_id];
          }
        }

      }

      /*  Heterogeneous combustion: C + CO2 ---> 2 CO */

      //  Calculation of HO2(T1)

      if (cm->ihtco2 == 1) {

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

          cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

          // Calculation of HCO(T2)

          for (int ige = 0; ige < ngazem; ige++)
            coefe[ige] = 0.;
          coefe[ico] = 1.;

          cs_real_t t2   = tfuel[c_id];
          cs_real_t xhco = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

          //  Calculation of HCO2(T1)

          for (int ige = 0; ige < ngazem; ige++)
            coefe[ige] = 0.;
          coefe[ico2] = 1.;
          cs_real_t t1    = cpro_temp[c_id];
          cs_real_t xhco2 = cs_coal_ht_convert_t_to_h_gas_by_yi(t1, coefe);

          for (int class_id = 0; class_id < cm->nclacp; class_id++) {
            cs_real_t gamhet;
            if (cvara_xck[class_id][c_id] > cs_coal_epsilon) {
              gamhet =  crom[c_id]*cpro_ghco2a[class_id][c_id]
                       * (  pow(cvara_xck[class_id][c_id], 2./3.)
                          + 2./3. * (  cvar_xck[class_id][c_id]
                                     - cvara_xck[class_id][c_id])
                                  / pow(cvara_xck[class_id][c_id], 1./3.));
              smbrs[c_id] -=  gamhet * (2.*wmole[ico]*xhco-wmole[ico2]*xhco2)
                                     / wmolat[cs_coal_atom_id_c]
                                     * cell_f_vol[c_id];
            }
          }

        }

      } /* ihtco2 */

      /* Heterogeneous combustion: C + H2O ---> CO + H2 */

      // Calculation of HO2(T1)

      if (cm->ihth2o == 1) {

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

          const cs_real_t t2 = tfuel[c_id];

          cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

          //  Calculation of HCO(T2)

          for (int ige = 0; ige < ngazem; ige++)
            coefe[ige] = 0.;
          coefe[ico] = 1.;

          cs_real_t xhco = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

          //  Calculation of HH2(T2)

          for (int ige = 0; ige < ngazem; ige++)
            coefe[ige] = 0.;
          coefe[ihy] = 1.;

          cs_real_t xhh2 = cs_coal_ht_convert_t_to_h_gas_by_yi(t2, coefe);

          // Calculation of HH2O(T1)

          for (int ige = 0; ige < ngazem; ige++)
            coefe[ige] = 0.;
          coefe[ih2o] = 1.;
          cs_real_t t1    = cpro_temp[c_id];
          cs_real_t xhh2o = cs_coal_ht_convert_t_to_h_gas_by_yi(t1, coefe);

          for (int class_id = 0; class_id < cm->nclacp; class_id++) {
            cs_real_t gamhet;
            if (cvara_xck[class_id][c_id] > cs_coal_epsilon) {
              gamhet =   crom[c_id]*cpro_ghh2oa[class_id][c_id]
                       * (  pow(cvara_xck[class_id][c_id], 2./3.)
                          + 2./3. * (  cvar_xck[class_id][c_id]
                                     - cvara_xck[class_id][c_id])
                                  / pow(cvara_xck[class_id][c_id], 1./3.));
              smbrs[c_id] -= gamhet  * (  wmole[ico]*xhco
                                        + wmole[ihy]*xhh2
                                        - wmole[ih2o]*xhh2o)
                                     / wmolat[cs_coal_atom_id_c]
                                     * cell_f_vol[c_id];
            }
          }

        }

      } /* ihth2o */

      BFT_FREE(tfuel);

      /* Drying */

      if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {

        for (int class_id = 0; class_id < cm->nclacp; class_id++) {

          int coal_id = cm->ichcor[class_id] -1;

          const cs_real_t *cvar_xwtcl = cs_field_by_id(cm->ixwt[class_id])->val;
          const cs_real_t *cvara_xwtcl
            = cs_field_by_id(cm->ixwt[class_id])->val_pre;
          const cs_real_t *cpro_csec = cs_field_by_id(cm->igmsec[class_id])->val;
          const cs_real_t *cpro_temp2 = cs_field_by_id(cm->itemp2[class_id])->val;
          const cs_real_t *cpro_x2 = cs_field_by_id(cm->ix2[class_id])->val;

          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

            cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

            // Calculation of H(H2O) at T2

            for (int ige = 0; ige < ngazem; ige++)
              coefe[ige] = 0.;
            coefe[ih2o] = 1.;

            cs_real_t t2 = cpro_temp2[c_id];
            cs_real_t hh2ov
              = cs_coal_ht_convert_t_to_h_gas_by_yi_with_drying(t2, coefe);

            // Contribution to explicit balance

            if (   cvara_xwtcl[c_id] > cs_coal_epsilon
                && cm->xwatch[coal_id] > cs_coal_epsilon) {
              cs_real_t aux =   crom[c_id] * cpro_csec[c_id]
                              * (cvar_xwtcl[c_id] / cpro_x2[c_id])
                              * (1. / cm->xwatch[coal_id]) *hh2ov;

              smbrs[c_id] -= aux * cell_f_vol[c_id];
            }

          } /* loop on cells */

        } /* loop on classes */

      } /* drying */

      BFT_FREE(cvar_xck);
      BFT_FREE(cvara_xck);
      BFT_FREE(cvar_xch);
      BFT_FREE(cvar_xnp);
      BFT_FREE(cpro_t2);
      BFT_FREE(cpro_gmhet);
      BFT_FREE(cpro_ghco2a);
      BFT_FREE(cpro_ghh2oa);
      BFT_FREE(cvar_xwt);
    }
  }

  BFT_FREE(w1);  // Not needed after this point

  /* Source terms on Y_HCN and Y_NO */

  if (   (cm->ieqnox== 1 && cm->imdnox == 0 && cs_glob_time_step->nt_cur > 1)
      && (fld_id == cm->iyhcn || fld_id == cm->iyno)) {

    // Pointers to source terms
    const cs_real_t *cpro_exp1 = cs_field_by_id(cm->ighcn1)->val;
    const cs_real_t *cpro_exp2 = cs_field_by_id(cm->ighcn2)->val;
    const cs_real_t *cpro_exp3 = cs_field_by_id(cm->ignoth)->val;

    // Molar Mass
    const cs_real_t wmhcn = wmole[ihcn];
    const cs_real_t wmno  = 0.030;
    const cs_real_t wmo2  = wmole[io2];

    if (fld_id == cm->iyhcn) {

      // HCN source term

      if (eqp->verbosity >= 1)
        bft_printf(_(log_st_fmt), fld_scal->name);

      // Arrays of pointers containing the fields values for each class
      // (loop on cells outside loop on classes)

      cs_real_t **cvara_xck, **cvara_xch, **cpro_gmhet;
      cs_real_t **cpro_gmdv1, **cpro_gmdv2;
      BFT_MALLOC(cvara_xck, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cvara_xch, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cpro_gmhet, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cpro_gmdv1, cm->nclacp, cs_real_t *);
      BFT_MALLOC(cpro_gmdv2, cm->nclacp, cs_real_t *);

      for (int class_id = 0; class_id < cm->nclacp; class_id++) {
        cvara_xck[class_id] = cs_field_by_id(cm->ixck[class_id])->val_pre;
        cvara_xch[class_id] = cs_field_by_id(cm->ixch[class_id])->val_pre;
        cpro_gmdv1[class_id] = cs_field_by_id(cm->igmdv1[class_id])->val;
        cpro_gmdv2[class_id] = cs_field_by_id(cm->igmdv2[class_id])->val;
        cpro_gmhet[class_id] = cs_field_by_id(cm->igmhet[class_id])->val;
      }

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        cs_real_t gmdev1[CS_COMBUSTION_MAX_COALS];
        cs_real_t gmdev2[CS_COMBUSTION_MAX_COALS];
        cs_real_t gmhet[CS_COMBUSTION_MAX_COALS];

        cs_real_t wmel = cpro_mmel[c_id];
        cs_real_t xo2 = cpro_yo2[c_id] * wmel / wmo2;

        cs_real_t aux
          =   cell_f_vol[c_id]*crom[c_id]
            * (  cpro_exp2[c_id]
               + cpro_exp1[c_id]*cvara_yno[c_id]*cpro_mmel[c_id]/wmno);

        smbrs[c_id] -= aux * cvara_var[c_id];
        rovsdt[c_id] += aux;

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
          gmdev1[coal_id] = 0.;
          gmdev2[coal_id] = 0.;
          gmhet [coal_id] = 0.;
        }

        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          int coal_id = cm->ichcor[class_id] -1;

          gmdev1[coal_id] +=   cpro_gmdv1[class_id][c_id] * crom[c_id]
                             * cvara_xch[class_id][c_id];
          gmdev2[coal_id] +=   cpro_gmdv2[class_id][c_id] * crom[c_id]
                             * cvara_xch[class_id][c_id];
          gmhet[coal_id] +=    cpro_gmhet[class_id][c_id]  * crom[c_id]
                             * pow(cvara_xck[class_id][c_id], 2./3.);
        }

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
          // % of pure nitrogen in the coal

          aux = -   cell_f_vol[c_id] * cm->fn[coal_id]
                   * wmhcn / (wmole[in2]/2.)
                   * (cm->qpr[coal_id] * (gmdev1[coal_id] + gmdev2[coal_id]));
          if (xo2 > 0.03) {
            aux -=   cell_f_vol[c_id] * cm->fn[coal_id]
                   * wmhcn / (wmole[in2]/2.)
                   * (1. - cm->qpr[coal_id] * cm->y2ch[coal_id])
                   / (1 - cm->y2ch[coal_id]) * gmhet[coal_id]
                                             * (1. - cm->xashch[coal_id]);
          }
          smbrs[c_id] += aux;
        }

      } /* loop on cells */

      BFT_FREE(cvara_xck);
      BFT_FREE(cvara_xch);
      BFT_FREE(cpro_gmhet);
      BFT_FREE(cpro_gmdv1);
      BFT_FREE(cpro_gmdv2);

    } /* fld_id == iyhcn */

    if (fld_id == cm->iyno) {

      // NO source term

      const cs_real_t *cpro_yn2 = cs_field_by_id(cm->iym1[in2])->val;

      if (eqp->verbosity >= 1)
        bft_printf(_(log_st_fmt), fld_scal->name);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t aux1 =   cell_f_vol[c_id] * crom[c_id]
                         * cpro_exp1[c_id] * cvara_yhcn[c_id]
                         * cpro_mmel[c_id] / wmhcn;
        cs_real_t aux2 =   cell_f_vol[c_id] * crom[c_id]
                         * cpro_exp2[c_id] * cvara_yhcn[c_id]
                         * wmno / wmhcn;
        cs_real_t aux3 =   cell_f_vol[c_id] * pow(crom[c_id], 1.5)
                         * cpro_exp3[c_id]  * cpro_yn2[c_id];

        smbrs[c_id] += - aux1*cvara_var[c_id] + aux2 + aux3;
        rovsdt[c_id] += aux1;
      }

    } /* fld_id == iyno */

  }  /* ieqnox == 1, imdnox == 0, nt_cur > 1, and fld_id in (iyhcn,  iyno) */

  if (   (cm->ieqnox== 1 && cm->imdnox == 1 && cs_glob_time_step->nt_cur > 1)
      && (fld_id == cm->iyhcn || fld_id == cm->iyno || fld_id == cm->iynh3)) {

    /* Source terms on Y_HCN and Y_NO */

    // Arrays of pointers containing the fields values for each class
    // (loop on cells outside loop on classes)

    cs_real_t **cvara_xck, **cvara_xch, **cpro_t2, **cpro_gmhet;
    cs_real_t **cpro_gmdv1, **cpro_gmdv2;
    BFT_MALLOC(cvara_xck, cm->nclacp, cs_real_t *);
    BFT_MALLOC(cvara_xch, cm->nclacp, cs_real_t *);
    BFT_MALLOC(cpro_t2, cm->nclacp, cs_real_t *);
    BFT_MALLOC(cpro_gmhet, cm->nclacp, cs_real_t *);
    BFT_MALLOC(cpro_gmdv1, cm->nclacp, cs_real_t *);
    BFT_MALLOC(cpro_gmdv2, cm->nclacp, cs_real_t *);

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cvara_xck[class_id] = cs_field_by_id(cm->ixck[class_id])->val_pre;
      cvara_xch[class_id] = cs_field_by_id(cm->ixch[class_id])->val_pre;
      cpro_t2[class_id] = cs_field_by_id(cm->itemp2[class_id])->val;
      cpro_gmhet[class_id] = cs_field_by_id(cm->igmhet[class_id])->val;
      cpro_gmdv1[class_id] = cs_field_by_id(cm->igmdv1[class_id])->val;
      cpro_gmdv2[class_id] = cs_field_by_id(cm->igmdv2[class_id])->val;
    }

    // Pointers to NO gas phase source terms
    const cs_real_t *cpro_exp1 = cs_field_by_id(cm->ighcn1)->val;
    const cs_real_t *cpro_exp2 = cs_field_by_id(cm->ighcn2)->val;
    const cs_real_t *cpro_exp3 = cs_field_by_id(cm->ignoth)->val;
    const cs_real_t *cpro_exp4 = cs_field_by_id(cm->ignh31)->val;
    const cs_real_t *cpro_exp5 = cs_field_by_id(cm->ignh32)->val;
    const cs_real_t *cpro_exprb = cs_field_by_id(cm->igrb)->val;

    cs_real_t *cpro_cnorb = cs_field_by_id(cm->ifnh3d)->val;
    cs_real_t *cpro_fnoch = cs_field_by_id(cm->ifnh3c)->val;
    cs_real_t *cpro_cnohc = cs_field_by_id(cm->icnohc)->val;
    cs_real_t *cpro_fnohc = cs_field_by_id(cm->ifnohc)->val;
    cs_real_t *cpro_fnonh = cs_field_by_id(cm->ifnonh)->val;
    cs_real_t *cpro_fnoth = cs_field_by_id(cm->ifnoth)->val;
    cs_real_t *cpro_cnonh = cs_field_by_id(cm->icnonh)->val;
    cs_real_t *cpro_fnh3d = cs_field_by_id(cm->ifnh3d)->val;
    cs_real_t *cpro_fnh3c = cs_field_by_id(cm->ifnh3c)->val;

    // Pointer to CHx1 and CHx2
    const cs_real_t *cpro_cyf1 = cs_field_by_id(cm->iym1[0])->val;
    const cs_real_t *cpro_cyf2 = cs_field_by_id(cm->iym1[1])->val;

    cs_real_t *cpro_fhcnr = cs_field_by_id(cm->ifhcnr)->val;
    cs_real_t *cpro_fhcnd = cs_field_by_id(cm->ifhcnd)->val;
    cs_real_t *cpro_fhcnc = cs_field_by_id(cm->ifhcnc)->val;

    // Molar mass

    const cs_real_t wmhcn = wmole[ihcn];
    const cs_real_t wmno  = 0.030;
    const cs_real_t wmnh3 = wmole[inh3];

    // Aliases to model arrays

    const cs_real_t *teno = cm->teno;
    const cs_real_t *chi2 = cm->chi2;

    if (fld_id == cm->iyhcn) {

      cs_array_real_fill_zero(n_cells, cpro_fhcnr);
      cs_array_real_fill_zero(n_cells, cpro_fhcnd);
      cs_array_real_fill_zero(n_cells, cpro_fhcnc);

      // HCN source term

      if (eqp->verbosity >= 1)
        bft_printf(_(log_st_fmt), fld_scal->name);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        // Molar mass of the gas mixture
        cs_real_t wmel = cpro_mmel[c_id];

        //  Coefficient of reactions HCN + O2 et HCN + NO
        cs_real_t aux =   cell_f_vol[c_id] * crom[c_id]
                        * (  cpro_exp2[c_id]
                           + cpro_exp1[c_id]*cvara_yno[c_id]*wmel/wmno);

        smbrs[c_id]  -= aux * cvara_var[c_id];
        rovsdt[c_id] += aux;

        // Reburning ?

        if (cm->irb == 1) {  // Chen's model

          // Remark: no expression depending on coal id in this loop
          for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

            cs_real_t ychx =   (cpro_cyf1[c_id] * wmel/cm->wmchx1)
                             + (cpro_cyf2[c_id] * wmel/cm->wmchx2);

            aux =   cell_f_vol[c_id] * wmhcn * cpro_exprb[c_id]
                  * cvara_yno[c_id] * wmel / wmno * ychx;

            smbrs[c_id] += aux;
            cpro_fhcnr[c_id] += aux;

          }

        }
        else if (cm->irb == 2) {  //  Dimitriou's model

          for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

            const cs_real_t *chxi = NULL;
            const cs_real_t *cpro_cyfi = NULL;
            cs_real_t wmchxi = 0;

            // Reburning by CHx1 or CHx2
            if (cpro_cyf1[c_id] > 0.) {
              chxi = cm->chx1;
              cpro_cyfi = cpro_cyf1;
              wmchxi = cm->wmchx1;
            }
            else if (cpro_cyf2[c_id] > 0.) {
              chxi = cm->chx2;
              cpro_cyfi = cpro_cyf2;
              wmchxi = cm->wmchx2;
            }

            if (chxi != NULL) {

              cs_real_t core1 = 0, core2 = 0, para2 = 0;

              // Number of secgments of the temperature discretization
              for (int ii = 0; ii < 7; ii++) {

                // We look for the interval teno[ii] < Tgas < teno[ii+1]
                if (   cpro_temp[c_id] >= teno[ii]
                    && cpro_temp[c_id] <  teno[ii+1]) {

                  // JJ indicates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  for (int jj = 0; jj < 4; jj++) {

                    cs_real_t s =   (cpro_temp[c_id] - teno[ii])
                                  / (teno[ii+1] - teno[ii]);

                    // We look for the interval jj < chx1[coal_id] < jj + 1
                    if (chxi[coal_id] >= 4.) {
                      core1 =   cm->ka[ii][3]
                              + (cm->ka[ii+1][3] - cm->ka[ii][3]) * s;
                      core2 = cm->kb[ii][3]
                              + (cm->kb[ii+1][3] - cm->kb[ii][3]) * s;
                      para2 = chi2[ii] + (chi2[ii+1] - chi2[ii]) * s;
                    }
                    else if (chxi[coal_id] <= 1.) {
                      core1 =   cm->ka[ii][0]
                              + (cm->ka[ii+1][0] - cm->ka[ii][0]) * s;
                      core2 = cm->kb[ii][0]
                              + (cm->kb[ii+1][0] - cm->kb[ii][0]) * s;
                      para2 = chi2[ii] + (chi2[ii+1] - chi2[ii]) * s;
                    }
                    else if (   chxi[coal_id] >= jj+1
                             && chxi[coal_id] < jj+2) {
                      core1 =   cm->ka[ii][jj]
                              + (cm->ka[ii+1][jj+1] - cm->ka[ii][jj]) * s;
                      core2 = cm->kb[ii][jj]
                              + (cm->kb[ii+1][jj+1] - cm->kb[ii][jj]) * s;
                      para2 = chi2[ii] + (chi2[ii+1] - chi2[ii]) * s;
                    }

                  } // loop on jj

                } // in teno interval

              } // segments of temperature discretization

              // Reburning by CHx1 or CHx2
              if (chxi[coal_id] >= 3.) {
                aux =   (cell_f_vol[c_id] * wmhcn)
                      * ((core1 + core2) * para2)
                      * (cvar_yno[c_id] * crom[c_id] / wmno)
                      * (cpro_cyfi[c_id] * cpro_rom1[c_id] / wmchxi);
              }
              else {
                aux =   (cell_f_vol[c_id] * wmhcn)
                      * (core1 + core2)
                      * (cvar_yno[c_id] * crom[c_id] / wmno)
                      * (cpro_cyfi[c_id] * cpro_rom1[c_id] / wmchxi);
              }

              smbrs[c_id] += aux;
              cpro_fhcnr[c_id] += aux;

            } // Reburning by CHx1 of CHx2

          } /* loop on coals */

        } // cm->irb == 2, Dimitiou's model

        //  Initialization of variables

        cs_real_t gmdev1[CS_COMBUSTION_MAX_COALS];
        cs_real_t gmdev2[CS_COMBUSTION_MAX_COALS];
        cs_real_t gmhet[CS_COMBUSTION_MAX_COALS];

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
          gmdev1[coal_id] = 0.;
          gmdev2[coal_id] = 0.;
          gmhet [coal_id] = 0.;
        }

        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          const int coal_id = cm->ichcor[class_id] -1;

          cs_real_t mckcl1 =   (1.-cm->y1ch[coal_id]) * cm->a1ch[coal_id]
                             * exp(  -cm->e1ch[coal_id]
                                   / (  cs_physical_constants_r
                                      * cpro_t2[class_id][c_id]));

          cs_real_t mckcl2 =   (1.-cm->y2ch[coal_id]) * cm->a2ch[coal_id]
                             * exp(  -cm->e2ch[coal_id]
                                   / (  cs_physical_constants_r
                                      * cpro_t2[class_id][c_id]));

          // Forming rate of the first pyrolisis reaction
          gmdev1[coal_id] +=   cpro_gmdv1[class_id][c_id] * crom[c_id]
                             * cvara_xch[class_id][c_id];

          // Forming rate of the second pyrolisis reaction
          gmdev2[coal_id] +=   cpro_gmdv2[class_id][c_id] * crom[c_id]
                             * cvara_xch[class_id][c_id];

          if (cvara_xck[class_id][c_id] > cs_coal_epsilon) {
            // Reaction rate of the heterogeneous combustion
            gmhet[coal_id] +=     cpro_gmhet[class_id][c_id] * crom[c_id]
                                * pow(  cvara_xck[class_id][c_id]
                                      * (  (1./(mckcl2/mckcl1 + 1.))
                                           * cm->yhcnc1[coal_id]
                                         + (1./(mckcl1/mckcl2 + 1.))
                                           * cm->yhcnc2[coal_id]), 2./3.);
          }
        }

        //  Modified source term (new model of NOx)

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

          // Release of HCN during devolatilization
          aux =   -cell_f_vol[c_id]
                * (  gmdev1[coal_id] * cm->yhcnle[coal_id]
                   + gmdev2[coal_id] * cm->yhcnlo[coal_id]);

          // Release of HCN during the heterogeneous combustion according
          // to the value repnck[coal_id]

          aux -= cell_f_vol[c_id] * gmhet[coal_id];
          smbrs[c_id] += aux;

          // Source terms displaying
          cpro_fhcnd[c_id] -=  cell_f_vol[c_id]
                              * (  gmdev1[coal_id] * cm->yhcnle[coal_id]
                                 + gmdev2[coal_id] * cm->yhcnlo[coal_id]);

          cpro_fhcnc[c_id] -= cell_f_vol[c_id] * gmhet[coal_id];

        } // loop on coals

      } // loop on cells

    } /* fld_id == iyhcn */

    else if (fld_id == cm->iynh3) {

      cs_array_real_fill_zero(n_cells, cpro_fnh3d);

      // NH3 source term

      if (eqp->verbosity >= 1)
        bft_printf(_(log_st_fmt), fld_scal->name);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        // Molar mass of the gaseous mixture
        cs_real_t wmel = cpro_mmel[c_id];

        //  Coefficient of reactions NH3 + O2 and NH3 + NO
        cs_real_t aux  =   cell_f_vol[c_id] * crom[c_id]
                         * (  cpro_exp4[c_id]
                            + cpro_exp5[c_id]*cvara_yno[c_id]*wmel/wmno);

        smbrs[c_id]  -= aux * cvara_var[c_id];
        rovsdt[c_id] += aux;

        //  Initialization of variables

        cs_real_t gmdev1[CS_COMBUSTION_MAX_COALS];
        cs_real_t gmdev2[CS_COMBUSTION_MAX_COALS];

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
          gmdev1[coal_id] = 0.;
          gmdev2[coal_id] = 0.;
        }

        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          const int coal_id = cm->ichcor[class_id] -1;

          // Forming rate of the first pyrolisis reaction
          gmdev1[coal_id] +=   cpro_gmdv1[class_id][c_id] * crom[c_id]
                             * cvara_xch[class_id][c_id];

          // Forming rate of the second pyrolisis reaction
          gmdev2[coal_id] +=   cpro_gmdv2[class_id][c_id] * crom[c_id]
                             * cvara_xch[class_id][c_id];
        }

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

          // Release of NH3 during the devolatization.
          aux =   - cell_f_vol[c_id]
                * (  gmdev1[coal_id]*cm->ynh3le[coal_id]
                   + gmdev2[coal_id]*cm->ynh3lo[coal_id]);

          smbrs[c_id] += aux;

          // Source terms display (FIXME: an intensive field would be better)
          cpro_fnh3d[c_id] -=   cell_f_vol[c_id]
                              * (  gmdev1[coal_id]*cm->ynh3le[coal_id]
                                 + gmdev2[coal_id]*cm->ynh3lo[coal_id]);

          cpro_fnh3c[c_id] = 0;

        }

      }

    } /* fld_id == iynh3 */

    if (fld_id == cm->iyno) {

      // NO source term

      const cs_real_t *cpro_yn2 = cs_field_by_id(cm->iym1[in2])->val;

      cs_array_real_fill_zero(n_cells, cpro_cnorb);
      cs_array_real_fill_zero(n_cells, cpro_fnoch);

      if (eqp->verbosity >= 1)
        bft_printf(_(log_st_fmt), fld_scal->name);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        // Molar mass of the gaseous mixture
        const cs_real_t wmel = cpro_mmel[c_id];

        // Coefficient of reaction HCN + NO
        cs_real_t aux1 =   cell_f_vol[c_id] * crom[c_id]
                         * cpro_exp1[c_id] * cvara_yhcn[c_id]
                         * wmel / wmhcn;

        cpro_cnohc[c_id] = aux1 * cvara_var[c_id];

        // Coefficient of reaction HCN + O2
        cs_real_t aux2 =   cell_f_vol[c_id] * crom[c_id]
                         * cpro_exp2[c_id] * cvara_yhcn[c_id]
                         * wmno / wmhcn;

        cpro_fnohc[c_id] = aux2;

        // Coefficient of thermal NO
        //FIXME why is the N mass fraction not converted to a molar fraction ?
        cs_real_t aux3 =   cell_f_vol[c_id] * pow(crom[c_id], 1.5)
                         * cpro_exp3[c_id] * cpro_yn2[c_id];

        cpro_fnoth[c_id] = aux3;

        // Coefficient of reaction NH3 + O2 --> NO + ...
        cs_real_t aux4 =   cell_f_vol[c_id] * crom[c_id]
                         * cpro_exp4[c_id] * cvara_ynh3[c_id]
                         * wmno / wmnh3;

        cpro_fnonh[c_id] = aux4;

        //  Coefficient of reaction NH3 + NO --> N2 + ...
        cs_real_t aux5 =   cell_f_vol[c_id] * crom[c_id]
                         * cpro_exp5[c_id]*cvara_ynh3[c_id]
                         * wmel / wmnh3;

        cpro_cnonh[c_id] = aux5*cvara_var[c_id];

        // Reburning ?
        if (cm->irb == 1) { //  Chen's model

          // Remark: no expression depending on coal id in this loop
          for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

             cs_real_t ychx =   (cpro_cyf1[c_id] * wmel / cm->wmchx1)
                              + (cpro_cyf2[c_id] * wmel / cm->wmchx2);

             cs_real_t aux =   cell_f_vol[c_id] * wmhcn*cpro_exprb[c_id]
                             * cvara_yno[c_id] * wmel / wmno  * ychx;

             smbrs[c_id] -= aux;

             cpro_cnorb[c_id] = cpro_cnorb[c_id] + aux;

          }

        }
        else if (cm->irb == 2) {  //  Dimitiou's model

          for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {

            const cs_real_t *chxi = NULL;
            const cs_real_t *cpro_cyfi = NULL;
            cs_real_t wmchxi = 0;

            // Reburning by CHx1 or CHx2
            if (cpro_cyf1[c_id] > 0.) {
              chxi = cm->chx1;
              cpro_cyfi = cpro_cyf1;
              wmchxi = cm->wmchx1;
            }
            else if (cpro_cyf2[c_id] > 0.) {
              chxi = cm->chx2;
              cpro_cyfi = cpro_cyf2;
              wmchxi = cm->wmchx2;
            }

            if (chxi != NULL) {

              cs_real_t core1 = 0, core2 = 0, core3 = 0, para2 = 0;

              // Number of secgments of the temperature discretization
              for (int ii = 0; ii < 7; ii++) {

                // We look for the interval teno[ii] < Tgas < teno[ii+1]
                if (   cpro_temp[c_id] >= teno[ii]
                    && cpro_temp[c_id] <  teno[ii+1]) {

                  // JJ indicates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  for (int jj = 0; jj < 4; jj++) {

                    cs_real_t s =   (cpro_temp[c_id] - teno[ii])
                                  / (teno[ii+1] - teno[ii]);

                    // We look for the interval jj < chx1[coal_id] < jj + 1
                    if (chxi[coal_id] >= 4.) {
                      core1 =   cm->ka[ii][3]
                              + (cm->ka[ii+1][3] - cm->ka[ii][3]) * s;
                      core2 = cm->kb[ii][3]
                              + (cm->kb[ii+1][3] - cm->kb[ii][3]) * s;
                      core3 = cm->kc[ii][3]
                              + (cm->kc[ii+1][3] - cm->kc[ii][3]) * s;
                      para2 = chi2[ii] + (chi2[ii+1] - chi2[ii]) * s;
                    }
                    else if (chxi[coal_id] <= 1.) {
                      core1 =   cm->ka[ii][0]
                              + (cm->ka[ii+1][0] - cm->ka[ii][0]) * s;
                      core2 = cm->kb[ii][0]
                              + (cm->kb[ii+1][0] - cm->kb[ii][0]) * s;
                      core3 = cm->kc[ii][0]
                              + (cm->kc[ii+1][0] - cm->kc[ii][0]) * s;
                      para2 = chi2[ii] + (chi2[ii+1] - chi2[ii]) * s;
                    }
                    else if (   chxi[coal_id] >= jj+1
                             && chxi[coal_id] < jj+2) {
                      core1 =   cm->ka[ii][jj]
                              + (cm->ka[ii+1][jj+1] - cm->ka[ii][jj]) * s;
                      core2 = cm->kb[ii][jj]
                              + (cm->kb[ii+1][jj+1] - cm->kb[ii][jj]) * s;
                      core3 = cm->kc[ii][jj]
                              + (cm->kc[ii+1][jj+1] - cm->kc[ii][jj]) * s;
                      para2 = chi2[ii] + (chi2[ii+1] - chi2[ii]) * s;
                    }

                  } // loop on jj

                } // in segment

              } // loop on segments

              //  Reburning by CHx1 or CHx2
              cs_real_t auxrbi = 0;

              if (chxi[coal_id] >= 3.) {
                auxrbi =   (cell_f_vol[c_id] * wmno)
                         * ((core1 + core2 + core3) * para2)
                         * (cpro_cyfi[c_id] * cpro_rom1[c_id] / wmchxi)
                         * (cvar_yno[c_id] * crom[c_id] / wmno);
              }
              else {
                auxrbi =   (cell_f_vol[c_id] * wmno)
                         * (core1 + core2 + core3)
                         * (cpro_cyfi[c_id]*cpro_rom1[c_id] / wmchxi)
                         * (cvar_yno[c_id] * crom[c_id] / wmno);
              }

              smbrs[c_id] -= auxrbi;
              cpro_cnorb[c_id] += auxrbi;

            } /* cpro_cyf1[c_id] > 0. || cpro_cyf2[c_id] > 0. */

          } /* loop on coals */

        } /* irb == 2; Dimitious's model */

        //  Initialization

        cs_real_t gmhet[CS_COMBUSTION_MAX_COALS];
        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
          gmhet [coal_id] = 0.;
        }

        for (int class_id = 0; class_id < cm->nclacp; class_id++) {
          int coal_id = cm->ichcor[class_id] -1;

          if (cvara_xck[class_id][c_id] > cs_coal_epsilon) {

            cs_real_t mckcl1 =   (1.-cm->y1ch[coal_id]) * cm->a1ch[coal_id]
                               * exp(  - cm->e1ch[coal_id]
                                     / (  cs_physical_constants_r
                                        * cpro_t2[class_id][c_id]));

            cs_real_t mckcl2 =   (1.-cm->y2ch[coal_id]) * cm->a2ch[coal_id]
                               * exp(  - cm->e2ch[coal_id]
                                     / (  cs_physical_constants_r
                                        * cpro_t2[class_id][c_id]));

            // Reaction rate of the heterogeneous combustion
            gmhet[coal_id]
              +=  cpro_gmhet[class_id][c_id] * crom[c_id]
                * pow(  cvara_xck[class_id][c_id]
                      * (  (1./(mckcl2/mckcl1+1.))*cm->ynoch1[coal_id]
                         + (1./(mckcl1/mckcl2+1.))*cm->ynoch2[coal_id]),
                        2./3.);
          }

        } // loop on coals

        //  Coefficient of released NO during the heterogeneous combustion

        for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
          cs_real_t auxhet = -cell_f_vol[c_id] * gmhet[coal_id];
          cpro_fnoch[c_id] = cpro_fnoch[c_id] + auxhet;

          smbrs[c_id] +=  - aux1 * cvara_var[c_id]
                          - aux5 * cvara_var[c_id]
                          + aux2 + aux3 + aux4 + auxhet;

          rovsdt[c_id] += aux1 + aux5;
        }

      } /* loop on cells */

    } /* fld_id == iyno */

    BFT_FREE(cvara_xck);
    BFT_FREE(cvara_xch);
    BFT_FREE(cpro_t2);
    BFT_FREE(cpro_gmhet);
    BFT_FREE(cpro_gmdv1);
    BFT_FREE(cpro_gmdv2);

  } /* (ieqnox == 1, imdnox == 1, nt_cur > 1)
       and (fld_id in (yhcn, iyno, iynh3)) */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
