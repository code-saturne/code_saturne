/*============================================================================
 * Defines the source terms for the soot mass fraction and the precursor
 * number for soot model of Moss et al for one time step.
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

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_boundary_conditions.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_steady_laminar_flamelet_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_combustion_reconstruct_variance(int iprev);

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_steady_laminar_flamelet_source_terms.c
        Specific physic routine: STE/VTE and progress variable equations.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Defines the source terms for the soot mass fraction and the precursor
 *        number for soot model of Moss et al for one time step.
 *
 *  The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 *  \f$ rovsdt \f$ et \f$ smbrs \f$ could already contain source term
 *  and don't have to be erased but incremented.
 *
 *  For stability sake, only positive terms should be add in \f$ rovsdt \f$.
 *  There is no constrain for \f$ smbrs \f$.
 *
 *  For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
 *           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
 *           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
 *
 *  Here are set \f$ rovsdt \f$ and \f$ smbrs \f$ containning \f$ \rho \Omega \f$
 *   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
 *     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.C.s^{-1} \f$,
 *     for enthalpy: \f$ J.s^{-1} \f$)
 *   - \f$ rovsdt \f$ en \f$ kg.s^{-1} \f$
 *
 * \param[in]      fld_id        field id
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_steady_laminar_flamelet_source_terms(int        fld_id,
                                        cs_real_t  smbrs[],
                                        cs_real_t  rovsdt[])
{

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_real_t *cell_vol = fvq->cell_vol;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  const int *bc_type = cs_glob_bc_type;

  const int kivisl  = cs_field_key_id("diffusivity_id");
  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");

  /* Initialization
     -------------- */

  /* Coef. of SGS kinetic energy used for the variance dissipation computation */
  const cs_real_t coef_k = 7.0e-2;

  cs_field_t *f_sc = cs_field_by_id(fld_id);
  cs_field_t *f_fm = CS_F_(fm);

  cs_real_t *scal_pre = f_sc->val_pre;
  cs_equation_param_t *eqp_sc = cs_field_get_equation_param(f_sc);

  const int ifcvsl = cs_field_get_key_int(f_sc, kivisl);
  const cs_real_t *viscls = NULL;
  if (ifcvsl >= 0)
    viscls = cs_field_by_id(ifcvsl)->val;

  /* Writings
     -------- */

  if (eqp_sc->verbosity >= 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  _("Specific physics source terms for the variable %s\n"),
                  f_sc->name);

  cs_field_t *f_recvr = cs_field_by_name_try("reconstructed_fp2m");
  cs_field_t *f_fp2m = CS_F_(fp2m);

  cs_real_t *turb_diff = NULL;
  if (cs_glob_turb_model->model == 41) {
    /* Retrieve turbulent diffusivity value for the mixture fraction */
    const int t_dif_id = cs_field_get_key_int(f_fm, key_turb_diff);
    if (t_dif_id > -1)
      turb_diff = cs_field_by_id(t_dif_id)->val;

  /* --- Cuenot et al.:
   * STE: Prod := 0
   *      Disp := - (D + Dtur)/(C_k * Delta_les**2)*fp2m
   * VTE: Prod := 2*rho*(D + Dtur)*|grad(Z)|**2
   *      Disp := - (D + Dtur)/(C_k * Delta_les**2)*fp2m
   *
   * --- Pierce:
   * Progress variable equation:
   *      Prod := flamelet_lib(fm, fp2m, ki, progvar)
   */

    /* For the moment, this model for source computation
       is only available in LES */

    if (f_sc == f_fp2m) {

      cs_real_t *fp2m = f_fp2m->val_pre;

      /* Allocate a temporary array for the gradient reconstruction */
      cs_real_3_t *grad;
      BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

      cs_real_t *coefa_p, *coefb_p;
      BFT_MALLOC(coefa_p, n_b_faces, cs_real_t);
      BFT_MALLOC(coefb_p, n_b_faces, cs_real_t);

      /* Homogeneous Neumann on convective inlet on the
         production term for the variance */

      cs_real_t *coefap = f_fm->bc_coeffs->a;
      cs_real_t *coefbp = f_fm->bc_coeffs->b;

      /* Overwrite diffusion at inlets */
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        coefa_p[f_id] = coefap[f_id];
        coefb_p[f_id] = coefbp[f_id];

        if (bc_type[f_id] == CS_CONVECTIVE_INLET) {
          coefap[f_id] = 0.;
          coefbp[f_id] = 1.;
        }

      }

      cs_field_gradient_scalar(f_fm,
                               true, /* use_previous_t */
                               1,    /* inc */
                               grad);

      /* Put back the value */
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        coefap[f_id] = coefa_p[f_id];
        coefbp[f_id] = coefb_p[f_id];
      }

      BFT_FREE(coefa_p);
      BFT_FREE(coefb_p);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        const cs_real_t delta_les
          = cs_turb_xlesfl * pow(cs_turb_ales*cell_vol[c_id], cs_turb_bles);

        const cs_real_t cexp
          =   2.0 * (turb_diff[c_id] + viscls[c_id]) * cell_f_vol[c_id]
            * cs_math_3_dot_product(grad[c_id], grad[c_id])
            - ((turb_diff[c_id] + viscls[c_id])/ (coef_k
            * cs_math_pow2(delta_les)) * fp2m[c_id]) * cell_f_vol[c_id];

        const cs_real_t cimp = 0.;
        smbrs[c_id]  += cexp + cimp * scal_pre[c_id];
        rovsdt[c_id] += cs_math_fmax(-cimp, 0.);
      }

      BFT_FREE(grad);
    }
    else if (f_sc == cs_field_by_name_try("mixture_fraction_2nd_moment")) {

      cs_f_combustion_reconstruct_variance(1);

      cs_real_t *fp2m = f_recvr->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        const cs_real_t delta_les
          = cs_turb_xlesfl *pow(cs_turb_ales*cell_vol[c_id], cs_turb_bles);

        const cs_real_t cexp = - (  (turb_diff[c_id] + viscls[c_id])
                                  / (coef_k*cs_math_pow2(delta_les))*fp2m[c_id])
                               * cell_f_vol[c_id];

        const cs_real_t cimp = 0.;
        smbrs[c_id]  += cexp + cimp * scal_pre[c_id];
        rovsdt[c_id] += cs_math_fmax(-cimp, 0.);
      }
    }

  } /* End test on model = 41 */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 2) {

    if (f_sc == cs_field_by_name_try("progress_variable")) {
      cs_real_t *omega_c = cs_field_by_name("omega_c")->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t cexp = omega_c[c_id];
        const cs_real_t cimp = 0.;

        smbrs[c_id]  += cexp + cimp * scal_pre[c_id];
        rovsdt[c_id] += cs_math_fmax(-cimp, 0.);
      }

    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
