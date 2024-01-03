/*============================================================================
 * Boundary condition check.
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
#include "cs_boundary_conditions.h"
#include "cs_field_pointer.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_check.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_boundary_conditions_check.c
        Check boundary condition code.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief synchronize boundary condition error logging across MPI ranks.
 *
 * \param[in, out] nerloc       number of errors (local rank in, global out)
 * \param[in]      nerrcd       number of codes saved at error faces
 * \param[in, out] errcod       codes saved at one error face (local in,
 *                              broadcast out)
 */
/*----------------------------------------------------------------------------*/

static void
_synchronize_boundary_conditions_error(cs_gnum_t  nerloc,
                                       int        nerrcd,
                                       int        errcod[])
{
  if (cs_glob_n_ranks > 1) {
    int irkerr = -1;
    if (nerloc > 0)
      irkerr = cs_glob_rank_id;

    cs_parall_counter(&nerloc, 1);
    if (nerloc != 0) {
      cs_parall_max(1, CS_INT_TYPE, &irkerr);
      cs_parall_bcast(irkerr, nerrcd, CS_INT_TYPE, errcod);
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check boundary condition code.
 *
 * \param[in,out]  bc_type      face boundary condition type
 * \param[in,out]  ale_bc_type  ale boundary condition type
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_check(int  bc_type[],
                             int  ale_bc_type[])

{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const int n_fields = cs_field_n_fields();

  cs_turb_model_type_t iturb = cs_glob_turb_model->iturb;
  int itytur = cs_glob_turb_model->itytur;

  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;

  /* Type ids and names */

  const int type_id_vel[10] = {1, 2, 3, 4, 5, 6, 9, 11, 13, 14};
  const int type_id_p[7]    = {1, 2, 3, 11, 12, 13, 15};
  const int type_id_turb[7] = {1, 2, 3, 5, 6, 13, 4};
  const int type_id_sc[11]  = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15};

  /* Initializations
     =============== */

  cs_real_t *bpro_rough = NULL, *bpro_rough_t = NULL;

  cs_field_t *f_rough = cs_field_by_name_try("boundary_roughness");
  if (f_rough != NULL) {
    bpro_rough = f_rough->val;
    bpro_rough_t = f_rough->val;
  }

  cs_field_t *f_rough_t = cs_field_by_name_try("boundary_thermal_roughness");
  if (f_rough_t != NULL)
    bpro_rough_t = f_rough_t->val;

  int *icodcl_vel = CS_F_(vel)->bc_coeffs->icodcl;

  /* Boundary conditions check
     ========================= */

  /*  Initialization
      -------------- */

  /* In cs_user_boundary_conditions, we remain flexible on
     boundary conditions specifications for variables.
     Nevertheless, to limit the range of tests, we give
     the following constraints:

     - no friction conditions on pressure
     - consistency between velocity and pressure boundary conditions
     - consistency between velocity and turbulence boundary conditions
  */

  /* Check that all BC are initialized
     --------------------------------- */

  char name_init[32];
  cs_gnum_t n_init_error = 0;
  int icodcl_init = -1;

  /* First fast loop */

  bool indef = false;

  for (int ii = 0; ii < n_fields; ii++) {

    cs_field_t *f_var = cs_field_by_id(ii);

    if (!(f_var->type & CS_FIELD_VARIABLE))
      continue;

    if (f_var->type & CS_FIELD_CDO)
      continue;

    const int *icodcl_var = (const int *)f_var->bc_coeffs->icodcl;
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (icodcl_var[f_id] == 0)
        indef = true;
    }

  }

  /* Second (more costly) loop if problem detected above */

  if (indef) {

    for (int ii = 0; ii < n_fields; ii++) {

      cs_field_t *f_var = cs_field_by_id(ii);

      if (!(f_var->type & CS_FIELD_VARIABLE))
        continue;

      if (f_var->type & CS_FIELD_CDO)
        continue;

      const int *icodcl_var = (const int *)f_var->bc_coeffs->icodcl;
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (icodcl_var[f_id] == 0) {
          if (bc_type[f_id] > 0)
            bc_type[f_id] = - bc_type[f_id];

          strncpy(name_init, f_var->name, 31); name_init[31] = '\0';
          icodcl_init = f_id;
          n_init_error = n_init_error + 1;
        }
      }
    }

  }

  /* Checking admissibility of boundary conditions
     --------------------------------------------- */

  cs_gnum_t n_vel_error = 0;
  int icodcl_v = -1;

  bool iok_rough = false;

  /* Eligible conditions for velocity components */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    bool eligible = false;
    for (int i = 0; i < 10; i++) {
      if (icodcl_vel[f_id] == type_id_vel[i])
        eligible = true;
    }

    if (!(eligible)) {

      if (bc_type[f_id] > 0)
        bc_type[f_id] = - bc_type[f_id];

      icodcl_v = icodcl_vel[f_id];
      n_vel_error = n_vel_error + 1;
    }

    /* Check roughness if rough wall function */
    if (icodcl_vel[f_id] == 6) {

      iok_rough = false;
      if (f_rough == NULL)
        iok_rough = true;
      else if (bpro_rough[f_id] <= 0.0)
        iok_rough = true;

      if (iok_rough) {
        if (bc_type[f_id] > 0)
          bc_type[f_id] = - bc_type[f_id];

        icodcl_v = icodcl_vel[f_id];
        n_vel_error = n_vel_error + 1;
      }
    }
  }

  /* Eligible conditions for pressure */

  cs_gnum_t n_p_error = 0;
  int icodcl_pr = -1;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    int *icodcl_p = CS_F_(p)->bc_coeffs->icodcl;

    bool eligible = false;
    for (int i = 0; i < 7; i++) {
      if (icodcl_p[f_id] == type_id_p[i])
        eligible = true;
    }

    if (!(eligible)) {

      if (bc_type[f_id] > 0)
        bc_type[f_id] = - bc_type[f_id];

      icodcl_pr = icodcl_p[f_id];
      n_p_error = n_p_error + 1;

    }
  }

  /* Eligible conditions for turbulence model */

  cs_field_t *f_turb_list[8] = {CS_F_(omg), CS_F_(phi), CS_F_(eps),
                                CS_F_(k), CS_F_(f_bar), CS_F_(alp_bl),
                                CS_F_(nusa), CS_F_(rij)};
  cs_field_t *f_turb_fields[8];
  int        *turb_icodcl[8];

  int cpt_turb = 0;

  for (int ii = 0; ii < 8; ii++) {
    const cs_field_t *f_turb = f_turb_list[ii];

    /* Map pointers for later access */
    if (f_turb != NULL) {
      if (!(f_turb->type & CS_FIELD_VARIABLE))
        continue;
      if (f_turb->type & CS_FIELD_CDO)
        continue;
      f_turb_fields[cpt_turb] = f_turb_list[ii];
      turb_icodcl[cpt_turb] = f_turb->bc_coeffs->icodcl;
      cpt_turb += 1;
    }
  }

  int id_turb = -1;
  cs_gnum_t n_turb_error = 0;
  int icodcl_turb = -1;

  for (int ii = 0; ii < cpt_turb; ii++) {

    cs_field_t *f_turb = f_turb_fields[ii];
    int *icodcl = f_turb->bc_coeffs->icodcl;

    cs_gnum_t n_field_errors = 0;
    int n_allowed_codes = (f_turb == CS_F_(rij)) ? 7 : 6;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      bool bc_eligible = false;

      for (int i = 0; i < n_allowed_codes; i++) {
        if (icodcl[f_id] == type_id_turb[i])
          bc_eligible = true;
      }

      if (!(bc_eligible)) {
        if (bc_type[f_id] > 0)
          bc_type[f_id] = - bc_type[f_id];
        n_field_errors += 1;
      }

    } /* End loop en boundary faces */

    if (n_field_errors > 0)
      id_turb = f_turb->id;

    n_turb_error += n_field_errors;

  } /* End loop on turbulent variables */

  /* No rough wall with EBRSM */

  cs_gnum_t n_ebrsm_rough_errors = 0;

  if (iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

    for (int jj = 0; jj < n_fields; jj++) {

      cs_field_t *f_var = cs_field_by_id(jj);

      if (!(f_var->type & CS_FIELD_VARIABLE))
        continue;

      int *icodcl_var = f_var->bc_coeffs->icodcl;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

        if (icodcl_var[f_id] == 6) {
          if (bc_type[f_id] > 0)
            bc_type[f_id] = - bc_type[f_id];

          n_ebrsm_rough_errors += 1;

        }
      }

      if (n_ebrsm_rough_errors > 0)
        id_turb = f_var->id;

      n_turb_error += n_ebrsm_rough_errors;

    }
  } /* End for EBRSM */

  /* Admissible conditions for additional transported variables
     (scalars or vectors) */

  const char *sc_name = NULL, *sc_vf_name = NULL;
  int icodcl_scal = -1, icodcl_vf_sc = -1;
  cs_gnum_t n_scal_error = 0, n_scal_vf_error = 0;
  bool iok_rough_sc = false;

  for (int ii = 0; ii < n_fields; ii++) {

    cs_field_t *f_sc = cs_field_by_id(ii);

    if (!(f_sc->type & CS_FIELD_VARIABLE))
      continue;
    if (f_sc->bc_coeffs == NULL)
      continue;
    if (cs_field_get_key_int(f_sc, keysca) <= 0)
      continue;

    const int iscavr = cs_field_get_key_int(f_sc, kscavr);
    int *icodcl_sc = f_sc->bc_coeffs->icodcl;

    cs_gnum_t n_field_errors = 0;
    cs_gnum_t n_field_vf_errors = 0;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      bool eligible_sc = false;
      for (int i = 0; i < 11; i++) {
        if (icodcl_sc[f_id] == type_id_sc[i])
          eligible_sc = true;
      }

      if (   !(eligible_sc)
          || (icodcl_sc[f_id] == 12 && f_sc->dim > 1)
          || (icodcl_sc[f_id] == 4  && f_sc->dim == 1)
          || (icodcl_sc[f_id] == 11 && f_sc->dim == 1)
          || (icodcl_sc[f_id] == 14 && f_sc->dim == 1)) {

        if (bc_type[f_id] > 0)
          bc_type[f_id] = - bc_type[f_id];

        sc_name = f_sc->name;
        icodcl_scal = icodcl_sc[f_id];
        n_field_errors += 1;
      }

      if (iscavr >= 0) {
        if (   icodcl_sc[f_id] == 5
            || icodcl_sc[f_id] == 6
            || icodcl_sc[f_id] == 15) {
          if (bc_type[f_id] > 0)
            bc_type[f_id] = - bc_type[f_id];

          sc_vf_name = f_sc->name;
          icodcl_vf_sc = icodcl_sc[f_id];
          n_field_vf_errors += 1;
        }
      }

      /* Check roughness if rough wall function */
      if (icodcl_sc[f_id] == 6) {
        bool ierr_rough_sc = false;

        if (f_rough_t == NULL)
          ierr_rough_sc = true;
        else if (bpro_rough_t[f_id] <= 0.0)
          ierr_rough_sc = true;

        if (ierr_rough_sc) {
          if (bc_type[f_id] > 0)
            bc_type[f_id] = - bc_type[f_id];

          icodcl_scal = icodcl_sc[f_id];
          n_field_errors += 1;
        }
      }
    }

    n_scal_error += n_field_errors;
    n_scal_vf_error += n_field_vf_errors;

  }

  /* Check if the gravity is non zero in case of free-surface */
  int iok_ale = 0;
  if (cs_glob_ale >= CS_ALE_LEGACY) {

    const cs_real_t grav2 = cs_math_3_dot_product(gxyz, gxyz);
    if (grav2 <= cs_math_pow2(cs_math_epzero)) {

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (ale_bc_type[f_id] == CS_BOUNDARY_ALE_FREE_SURFACE)
        iok_ale = 1;
      }

      cs_parall_max(1, CS_INT_TYPE, &iok_ale);
      if (iok_ale != 0) {
        bft_error(__FILE__, __LINE__, 0,
                  _("Free surface boundary conditions in ALE\n"
                    "must be combined with a non zero gravity term!\n"
                    "\n"
                    "Check parameters or modify the gravity."));
      }
    }

  }

  /* Inter-variable consistency check
     -------------------------------- */

  /* Velocity-pressure consistency */

  /* Remarks :
       No strict speed/pressure consistency rule.
       In the past, a Dirichlet was imposed on the pressure to
       inlet/outlet, but this doesn't seem imperative. The test
       is left as a comment to be adapted if necessary.

       if (icodcl_vel == 9) then
         if (icodcl(ipriph) != 1) then
           nstoup = nstoup + 1
         endif
       endif
  */

  /* Velocity-turbulence consistency */

  cs_gnum_t n_turb_consis_error = 0;
  int icodcl_consis_turb[2] = {-1, -1};

  int id_turb_consis = -1;
  int _icodcl[2] = {5, 6};

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    for (int icod = 0; icod < 2; icod++) {

      int n_icod_turb_wall = 0;

      for (int ii = 0; ii < cpt_turb; ii++) {
        int *icodcl = turb_icodcl[ii];

        if (icodcl[f_id] == _icodcl[icod]) {
          n_icod_turb_wall += 1;
        }

      }

      if (   icodcl_vel[f_id] == _icodcl[icod]
          || n_icod_turb_wall > 0) {

        /* Verify if icodcl[f_id] == _icodcl[icod] for
           velocity and all turbulent variables */

        if (n_icod_turb_wall != cpt_turb) {
          for (int ii = 0; ii < cpt_turb; ii++) {
            cs_field_t *f_turb = f_turb_fields[ii];
            int *icodcl = turb_icodcl[ii];

            if (icodcl[f_id] != _icodcl[icod]) {
              id_turb_consis = f_turb->id;
              icodcl_consis_turb[0] = icodcl[f_id];
            }
          }
        }

        if (   icodcl_vel[f_id] != _icodcl[icod]
            || n_icod_turb_wall != cpt_turb) {

          if (bc_type[f_id] > 0) {
            bc_type[f_id] = - bc_type[f_id];
          }
          icodcl_consis_turb[1] = icodcl_vel[f_id];
          n_turb_consis_error += 1;
        }

      }

    } /* End rough and smooth wall process */

    if (CS_F_(rij) != NULL) {

      cs_field_t *f_rij = CS_F_(rij);
      if (    (f_rij->type & CS_FIELD_VARIABLE)
          && !(f_rij->type & CS_FIELD_CDO)) {

        int *icodcl_rij = CS_F_(rij)->bc_coeffs->icodcl;

        if (icodcl_vel[f_id] == 4 || icodcl_rij[f_id] == 4) {

          bool icod_turb_outlet = true;

          for (int ii = 0; ii < cpt_turb; ii++) {

            cs_field_t *f_turb = f_turb_fields[ii];
            if (f_turb == CS_F_(rij))
              continue;
            int *icodcl = turb_icodcl[ii];

            if (icodcl[f_id] != 3) {
              icod_turb_outlet = false;
              id_turb_consis = ii;
              icodcl_consis_turb[0] = icodcl[f_id];
            }

          }

          if (   icodcl_vel[f_id] != 4
              || icodcl_rij[f_id] != 4
              || icod_turb_outlet == false) {

            if (bc_type[f_id] > 0) {
              bc_type[f_id] = - bc_type[f_id];
            }

            if (icodcl_rij[f_id] != 4)
              icodcl_consis_turb[0] = icodcl_rij[f_id];

            icodcl_consis_turb[1] = icodcl_vel[f_id];
            n_turb_consis_error += 1;
          }
        }
      }

    }

  } /* End loop on boundary faces */

  /* Velocity-scalars consistency */

  cs_gnum_t n_sc_consis_error = 0;
  int icodcl_consis_sc[2] = {-1, -1};
  int id_sc_consis = -1;

  if (itytur == 2 || itytur == 3) {

    for (int ii = 0; ii < n_fields; ii++) {

      cs_field_t *f_sc = cs_field_by_id(ii);

      if (!(f_sc->type & CS_FIELD_VARIABLE))
        continue;
      if (cs_field_get_key_int(f_sc, keysca) <= 0)
        continue;

      int *icodcl_sc = f_sc->bc_coeffs->icodcl;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

        if (   icodcl_sc[f_id] == 5
            && icodcl_vel[f_id] != 5) {

          if (bc_type[f_id] > 0)
            bc_type[f_id] = - bc_type[f_id];

          icodcl_consis_sc[0] = icodcl_sc[f_id];
          icodcl_consis_sc[1] = icodcl_vel[f_id];
          id_sc_consis = f_id;
          n_sc_consis_error += 1;
        }

      }
    }
  }

  /* Summary
     ------- */

  int error = 0;

  if (n_init_error > 0 || n_vel_error > 0 || n_p_error > 0 ||
      n_turb_error > 0 || n_scal_error > 0 || n_scal_vf_error > 0 ||
      n_turb_consis_error > 0 || n_sc_consis_error > 0)
    error = 1;

  cs_parall_max(1, CS_INT_TYPE, &error);

  if (error != 0) {

    _synchronize_boundary_conditions_error(n_init_error, 1, &icodcl_init);
    if (n_init_error != 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Uninitialized boundary conditions\n"
           "@   Number of boundary faces: %llu\n"
           "@     variable: %s\n"
           "@     icodcl last face: %d\n"
           "@\n"), (unsigned long long)n_init_error, name_init, icodcl_init);

    }

    _synchronize_boundary_conditions_error(n_vel_error, 1, &icodcl_v);

    if (n_vel_error != 0) {
      char string[10];
      if (iok_rough)
        strncpy(string, "roughness", 10);
      else
        strncpy(string, "velocity", 9);

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Unexpected boundary conditions\n"
           "@   Number of boundary faces: %llu\n"
           "@     variable: %s\n"
           "@     icodcl for last face: %d\n"
           "@\n"), (unsigned long long)n_init_error, string, icodcl_v);
    }

    _synchronize_boundary_conditions_error(n_p_error, 1, &icodcl_pr);

    if (n_p_error != 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Unexpected boundary conditions\n"
           "@   Number of boundary faces: %llu\n"
           "@     variable: pressure\n"
           "@     icodcl for last face: %d\n"
           "@\n"), (unsigned long long)n_p_error, icodcl_pr);
    }

    _synchronize_boundary_conditions_error(n_turb_error, 1, &icodcl_turb);

    if (n_turb_error != 0) {

      const cs_field_t *f_turb = cs_field_by_id(id_turb);
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Unexpected boundary conditions\n"
           "@   Number of boundary faces: %llu\n"
           "@     variable: %s\n"
           "@     icodcl for last face: %d\n"
           "@\n"),
         (unsigned long long)n_turb_error, f_turb->name, icodcl_turb);

      cs_parall_counter(&n_ebrsm_rough_errors, 1);

      if (n_ebrsm_rough_errors > 0) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@\n"
             "@ @@ Error: abort after boundary conditions check\n"
             "@    ======\n"
             "@\n"
             "@    Rough wall boundary conditions incompatible\n"
             "@      with EBRSM model.\n"
             "@\n"
             "@    Verify the parameters and boundary conditions.\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"));
      }

    }

    /* Admissible conditions for additional transporter variables */

    _synchronize_boundary_conditions_error(n_scal_error, 1, &icodcl_scal);

    if (n_scal_error != 0) {
      char string[20];

      if (iok_rough_sc)
        strncpy(string, "roughness", 10);
      else
        strncpy(string, sc_name, 20);

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Unexpected boundary conditions\n"
           "@   Number of boundary faces: %llu\n"
           "@     variable: %s\n"
           "@     icodcl for last face: %d\n"
           "@\n"), (unsigned long long)n_scal_error, string, icodcl_scal);
    }

    _synchronize_boundary_conditions_error(n_scal_vf_error, 1, &icodcl_vf_sc);
    if (n_scal_vf_error != 0) {
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Unexpected boundary conditions\n"
           "@   Number of boundary faces: %llu\n"
           "@     variable: %s\n"
           "@     icodcl for last face: %d\n"
           "@\n"),
         (unsigned long long)n_scal_vf_error, sc_vf_name, icodcl_vf_sc);
    }

    /* Velocity-turbulence consistency */

    _synchronize_boundary_conditions_error(n_turb_consis_error, 2, icodcl_consis_turb);
    if (n_turb_consis_error != 0) {
      const cs_field_t *f_turb = cs_field_by_id(id_turb_consis);

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Incoherency boundary conditions velocity-turbulent variables\n"
           "@   Total consistency error number: %llu\n"
           "@\n"
           "@ Turbulent variable (last): %s\n"
           "@   icodcl for last face: %d\n"
           "@   velocity icodcl: %d\n"
           "@\n"),
         (unsigned long long)n_turb_consis_error, f_turb->name,
         icodcl_consis_turb[0], icodcl_consis_turb[1]);
    }

     /* Velocity-scalars consistency */

    _synchronize_boundary_conditions_error(n_sc_consis_error, 2, icodcl_consis_sc);
    if (n_sc_consis_error != 0) {
      const cs_field_t *f_sc = cs_field_by_id(id_sc_consis);

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ Incoherent boundary conditions velocity-scalars\n"
           "@   Total consistency error number: %llu\n"
           "@\n"
           "@ Scalar example: %s\n"
           "@   icodcl scalar last face: %d\n"
           "@   icodcl velocity: %d\n"
           "@\n"),
         (unsigned long long)n_sc_consis_error, f_sc->name,
         icodcl_consis_sc[0], icodcl_consis_sc[1]);

    }

    if (   n_init_error > 0
        || n_scal_error > 0
        || n_scal_vf_error > 0
        || n_sc_consis_error > 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n"
           "@ @@ Error: abort after boundary conditions check\n"
           "@    ======\n"
           "@\n"
           "@    Uninitialized boundary conditions        : %llu\n"
           "@    Unexpected  boundary conditions:\n"
           "@      on the scalars                       : %llu\n"
           "@      on the scalars representing\n"
           "@                                      a variance  : %llu\n"
           "@    Incoherencies:\n"
           "@      between velocity and scalars         : %llu\n"
           "@\n"
           "@    The calculation will not be run.\n"
           "@\n"
           "@    Verify the parameters given via the interface or\n"
           "@      cs_user_boundary_conditions.\n"
           "@\n\n"),
         (unsigned long long)n_init_error, (unsigned long long)n_scal_error,
         (unsigned long long)n_scal_vf_error,
         (unsigned long long)n_sc_consis_error);

    }

    if (   n_vel_error > 0
        || n_p_error > 0
        || n_turb_error > 0
        || n_turb_consis_error > 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n"
           "@ @@ Error: abort after boundary conditions check\n"
           "@    ======\n"
           "@\n"
           "@    Unexpected boundary conditions:\n"
           "@      on the velocity                      : %llu\n"
           "@      on the pressure                      : %llu\n"
           "@      on turbulent variable                : %llu\n"),
         (unsigned long long)n_vel_error, (unsigned long long)n_p_error,
         (unsigned long long)n_turb_error);

      if (n_turb_error > 0) {
        const cs_field_t *f_turb = cs_field_by_id(id_turb);
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@  --> turbulent variable: %s\n"
             "@\n"), f_turb->name);
      }
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@     Incoherencies between velocity and\n"
           "@       turbulent variable                  : %llu\n"),
         (unsigned long long)n_turb_consis_error);

      if (n_turb_consis_error > 0) {
        const cs_field_t *f_turb_consis = cs_field_by_id(id_turb);
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("@ --> turbulent variable: %s\n"), f_turb_consis->name);
      }
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@    Verify the parameters and boundary conditions.\n"
           "@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n"));

    }

    cs_boundary_conditions_error(bc_type, NULL);
  }
}

/*---------------------------------------------------------------------------- */

END_C_DECLS
