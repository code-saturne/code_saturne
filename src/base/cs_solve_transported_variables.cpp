/*============================================================================
 * Resolution of source term convection diffusion equations
 * for scalars in a time step.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"
#include "cs_base.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_cf_energy.h"
#include "cs_cf_model.h"
#include "cs_elec_model.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_hgn_source_terms_step.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_restart.h"
#include "cs_solve_equation.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solve_transported_variables.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_solve_transported.c
  * \brief  Resolution of source term convection diffusion equations
  *         for scalars in a time step.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_kinetics_rates_compute(void);

void
cs_f_specific_physic_init(void);

void
cs_f_coal_masstransfer(void);

void
cs_f_max_mid_min_progvar(const cs_real_t  *zmo,
                         cs_real_t        *cmax,
                         cs_real_t        *cmid,
                         cs_real_t        *cmin);

void
cs_f_compute_gaseous_chemistry(cs_real_t dt[]);

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool _initialized = false;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *
 * \brief Resolution of source term convection diffusion equations
 *        for scalars in a time step.
 *
 * \param[in]     iterns        Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_transported_variables(int iterns)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;

  const cs_atmo_chemistry_t *atmo_chem = cs_glob_atmo_chemistry;
  const int nespg = atmo_chem->n_species;

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  cs_real_t *dt = CS_F_(dt)->val;
  const int n_fields = cs_field_n_fields();

  /* Initialization
     -------------- */

  /* Allocate temporary arrays for the species resolution */
  cs_real_t *i_visc, *b_visc;
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);

  /* Atmospheric chemistry => all chemical fields are not buoyant */

  if (   atmo_chem->model >= 1
      && atmo_chem->aerosol_model == CS_ATMO_AEROSOL_OFF
      && nespg > 0 && iterns == -1) {
    /* Computation of kinetics rates */
    cs_f_kinetics_rates_compute();
  }

  /* Handle model or specific physics scalars
     ---------------------------------------- */

  int nscapp = 0;

  for (int ii = 0; ii < n_fields; ii++) {

    cs_field_t *f_scal = cs_field_by_id(ii);

    if (!(f_scal->type & CS_FIELD_VARIABLE))
      continue;
    if (cs_field_get_key_int(f_scal, keysca) <= 0)
      continue;
    if (f_scal->type & CS_FIELD_USER)
      continue;

    nscapp += 1;
  }

  if (nscapp > 0) {

    if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] >= 1) {

      cs_f_specific_physic_init();

      if (_initialized == false && !cs_restart_present()) {
        _initialized = true;

        for (int ii = 0; ii < n_fields; ii++) {
          cs_field_t *f_scal = cs_field_by_id(ii);

          if (!(f_scal->type & CS_FIELD_VARIABLE))
            continue;
          if (cs_field_get_key_int(f_scal, keysca) <= 0)
            continue;
          if (f_scal->type & CS_FIELD_USER)
            continue;

          cs_real_t *cvar_var = f_scal->val;
          cs_real_t *cvara_var = f_scal->val_pre;
          int scal_dim = f_scal->dim;

          for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
            for (cs_lnum_t k = 0; k < scal_dim; k++)
              cvara_var[c_id*scal_dim + k] = cvar_var[c_id*scal_dim + k];

        }
      }
    }

    /* TS computations related to coal physics
       GMDEV1, GMDEV2, GMHET, GMDCH */

    if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] != -1)
      cs_f_coal_masstransfer();

    /* WARNING : For the clipping with ICLP = 1, scalars must be
       solved before their associated variance

       Loop over specific physical scalars.
       Instead, we can imagine coupled resolutions.
       Here, we give just one example.
    */

    for (int ii = 0; ii < n_fields; ii++) {
      cs_field_t *f_scal = cs_field_by_id(ii);

      if (!(f_scal->type & CS_FIELD_VARIABLE))
        continue;
      if (cs_field_get_key_int(f_scal, keysca) <= 0)
        continue;
      if (f_scal->type & CS_FIELD_USER)
        continue;

      /* Compressible scheme without shock:
         ---> Special processing for density, temperature and energy

         ISPECF indicator will be non-zero if the lower scalar is not
         solved via solve_equation_scalar.
      */

      int ispecf = 0;

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0 && iterns == -1) {

        if (f_scal == CS_F_(t_kelvin))
          ispecf = 1;
        else if (f_scal == CS_F_(e_tot)) {
          ispecf = 2;
          cs_cf_energy(f_scal->id);
        }
      }

      /* For compressible physics, already solved variables
         (or that should not be solved) are not solved.
         For other physics, we solve the scalars in order (no tests required).
      */

      if (ispecf == 0) {

        /* Variances and scalars

           iscavr: scalar/variance
           itspdv: computation of additional production and dissipation
                   terms, or not (itspdv = 1 : yes, 0 : no)

         if iscavr = -1
           scalar
           itspdv = 0
         else
           variance
           if iscavr > 0 and iscavr < nscal+1
             itspdv = 1
           else
             for the moment, we're stopping. Ultimately, combustionists
             will be able to give their own scalar associated with the
             variance and possibly reconstructed outside solve_equation_scalar.
             Additional work arrays could be constructed.
           end if
         end if
        */

        const int iscavr = cs_field_get_key_int(f_scal, kscavr);

        int itspdv = -1;
        if (iscavr == -1)
          itspdv = 0;
        else if (iscavr > 0)
          itspdv = 1;
        else {

          bft_error(__FILE__, __LINE__, 0,
          _("%s:Abort while solving scalars equations\n"
            "    ========\n"
            "    Scalar name = %s, id = %d\n"
            "    iscavr must be stricly positive or -1 integer\n"
            "    its value is %d\n"
            "\n\n"
            "  If iscavr(I) = -1, the scalar I is not a variance\n"
            "  If iscavr(I) > 0, the scalar I is a variance:\n"
            "    it is the variance of the fluctuations of the scalar J.\n"
            "    whose number is iscavr(I)\n"
            "\n"
            "  Check parameters"), __func__,
                    f_scal->name, f_scal->id, iscavr);
        }

        /* Specific process BC for gas combustion: steady laminar flamelet */

        if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0) {
          cs_field_t *f_fm = CS_F_(fm);

          cs_real_t *coefa_fm = f_fm->bc_coeffs->a;
          cs_real_t *coefb_fm = f_fm->bc_coeffs->b;
          cs_real_t *cvar_fm = f_fm->val;

          const int ifcvsl = cs_field_get_key_int(f_scal, kivisl);
          const cs_real_t *viscls = nullptr;
          if (ifcvsl >= 0)
            viscls = cs_field_by_id(ifcvsl)->val;

          /* Second moment of the mixing rate (only if mode_fp2m == 1) */
          if (f_scal == cs_field_by_name_try("mixture_fraction_2nd_moment")) {
            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              const cs_lnum_t c_id = b_face_cells[f_id];

              cs_real_t fmb = coefa_fm[f_id] + coefb_fm[f_id]*cvar_fm[c_id];
              cs_real_t hint = viscls[c_id] / b_dist[f_id];
              cs_real_t pimp = cs_math_pow2(fmb);
              cs_real_t hext = cs_math_infinite_r;

              cs_boundary_conditions_set_dirichlet_scalar(f_id,
                                                          f_scal->bc_coeffs,
                                                          pimp,
                                                          hint,
                                                          hext);
            }
          }

          if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 2) {

            if (f_scal == cs_field_by_name("progress_variable")) {
              for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
                const cs_lnum_t c_id = b_face_cells[f_id];

                cs_real_t fmb = coefa_fm[f_id] + coefb_fm[f_id]*cvar_fm[c_id];
                cs_real_t hint = viscls[c_id] / b_dist[f_id];

                cs_real_t cmax, cmid, cmin;
                cs_f_max_mid_min_progvar(&fmb, &cmax, &cmid, &cmin);
                cs_real_t pimp = cmid;
                cs_real_t hext = cs_math_infinite_r;

                cs_boundary_conditions_set_dirichlet_scalar(f_id,
                                                            f_scal->bc_coeffs,
                                                            pimp,
                                                            hint,
                                                            hext);
              }
            }
          } /* End CS_COMBUSTION_SLFM >= 2 */

        } /* End CS_COMBUSTION_SLFM >= 0 */

        if (f_scal->dim == 1) {
          cs_solve_equation_scalar(f_scal,
                                   iterns,
                                   itspdv,
                                   i_visc,
                                   b_visc);
        }
        else
          cs_solve_equation_vector(f_scal,
                                   iterns,
                                   i_visc,
                                   b_visc);

        /* --> Electrical versions
           Joule effect
           Electric arc
           Ion Conduction

           Real and imaginary j, E, j.E are computed
        */

        if (   (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 1
                || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1)
             && iterns == -1) {

          /*  We use the fact that the scalars are in the order
              H, PotR, [PotI], [A] to compute j, E and j.E after
              determining PotR [and PotI].
          */
          bool icalc = false;

          /* We can go after PotR if we are in arc or Joule without PotI */

          if (   cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1
              || cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 1
              || cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 3) {

            cs_field_t *f_elec_pot_r = cs_field_by_name("elec_pot_r");
            if (f_scal == f_elec_pot_r)
              icalc = true;
          }

          /*  We go after PotI if we are on Joule with PotI */
          if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 2
              || cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 4) {

            cs_field_t *f_elec_pot_i = cs_field_by_name("elec_pot_i");
            if (f_scal == f_elec_pot_i)
              icalc = true;
          }

          if (icalc) {
            /* Compute j, E and j.E */
            cs_elec_compute_fields(m, 1);

            /* Readjust electric variables j, j.E (and Pot, E) */
            if (   cs_glob_elec_option->ielcor == 1
                && cs_get_glob_time_step()->nt_cur > 1)
              cs_elec_scaling_function(m, fvq, dt);

          }

        }/* End CS_JOULE_EFFECT and CS_ELECTRIC_ARCS */

      } /* End test on specific physical scalar */

    } /* End of loop on specific physical models */

  } /* End if nscapp > 0 */

  /* Electric arcs:
     computation of magnetic field B and Laplace effect jxB */

  if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1 && iterns == -1)
    /* We use the fact that the scalars are in the order
       H, PotR, [PotI], [A] to compute A, B, jxB after
       determining and recalibrating j */
    cs_elec_compute_fields(m, 2);

  /* Compressible homogeneous two-phase model:
     return to equilibrium source term step
     for volume, mass, energy fractions */
  if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 1
       && cs_glob_cf_model->hgn_relax_eq_st >= 0)
    cs_hgn_source_terms_step(m);

  /* Standard user-defined scalars processing
     (Even if they're numbered first, they can be processed
     last if you want. We can also imagine by-passing this step
     if the particular physics model requires it.
  */

  /* Loop on user-defined scalars */

  for (int ii = 0; ii < n_fields; ii++) {

    cs_field_t *f_scal = cs_field_by_id(ii);

    if (!(f_scal->type & CS_FIELD_VARIABLE))
      continue;
    if (cs_field_get_key_int(f_scal, keysca) <= 0)
      continue;
    if (!(f_scal->type & CS_FIELD_USER))
      continue;

    /* Variances and scalars

    iscavr: scalaire/variance
    itspdv: computation of additional production and dissipation terms, or not
           (itspdv = 1 : yes, 0 : no)

         if iscavr = 0
           scalar
           itspdv = 0
         else
           variance
           if iscavr > 0 and iscavr < nscal+1
             itspdv = 1
           else
             for the moment, we're stopping. Ultimately, combustionists
             will be able to give their own scalar associated with the
             variance and possibly reconstructed outside solve_equation_scalar.
             Additional work arrays could be constructed.
           end if
         end if
    */

    const int iscavr = cs_field_get_key_int(f_scal, kscavr);
    int itspdv = -1;

    if (iscavr == -1)
      itspdv = 0;
    else if (iscavr > 0)
      itspdv = 1;
    else {
      bft_error(__FILE__, __LINE__, 0,
          _("%s:Abort while solving user-defined scalars equations\n"
            "    ========\n"
            "    Scalar name = %s, id = %d\n"
            "    iscavr must be stricly positive or -1 integer\n"
            "    its value is %d\n"
            "\n\n"
            "  If iscavr(I) = -1, the scalar I is not a variance\n"
            "  If iscavr(I) > 0, the scalar I is a variance:\n"
            "    it is the variance of the fluctuations of the scalar J.\n"
            "    whose number is iscavr(I)\n"
            "\n"
            "  Check parameters"), __func__,
            f_scal->name, f_scal->id, iscavr);
    }

    if (f_scal->dim == 1)
      cs_solve_equation_scalar(f_scal,
                               iterns,
                               itspdv,
                               i_visc,
                               b_visc);
    else
      cs_solve_equation_vector(f_scal,
                               iterns,
                               i_visc,
                               b_visc);

  } /* End of loop on user-defined scalars */

  /* Atmospheric gaseous chemistry
     Resolution of chemical evolution of species */

  if (   atmo_chem->model >= 1
      && atmo_chem->aerosol_model == CS_ATMO_AEROSOL_OFF
      && nespg > 0 && iterns == -1)
    cs_f_compute_gaseous_chemistry(dt);

  /* Atmospheric gas + aerosol chemistry */
  if (   atmo_chem->model >= 1
      && atmo_chem->aerosol_model != CS_ATMO_AEROSOL_OFF && iterns == -1)
   cs_atmo_aerosol_time_advance();

  /* Free memory */
  CS_FREE_HD(i_visc);
  CS_FREE_HD(b_visc);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
