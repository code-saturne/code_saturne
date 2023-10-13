/*============================================================================
 * Radiation solver boundary conditions treatment.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

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

#include "cs_array.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_boundary_zone.h"
#include "cs_log.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_ht_convert.h"
#include "cs_internal_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_physical_constants.h"

#include "cs_gui_radiative_transfer.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_wall_flux.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_bcs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_rad_transfer_bcs.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

int ipacli = 0;

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for fortran API
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Adjust radiative BC's for internal coupling.
 *
 * If emissivity has not been defined on solid faces, it is set to 0
 * on those faces.
 *
 * \param[in, out]  cpl     internal coupling structure
 * \param[in, out]  isothm  internal coupling BC type
 * \param[in, out]  xlamp   conductivity (W/m/K)
 * \param[in, out]  beps    emissivity
 * \param[in, out]  textp   outside temperature (K)
 */
/*----------------------------------------------------------------------------*/

static void
_set_internal_coupling_bcs(cs_internal_coupling_t  *cpl,
                           int                      isothm[],
                           cs_real_t                xlamp[],
                           cs_real_t                beps[],
                           cs_real_t                textp[])
{
  const cs_mesh_t *m = cs_glob_mesh;

  cs_lnum_t  n_local = 0, n_distant = 0;
  bool have_unset = false;
  const cs_lnum_t *faces_local = NULL, *faces_distant = NULL;

  cs_internal_coupling_coupled_faces(cpl,
                                     &n_local,
                                     &faces_local,
                                     &n_distant,
                                     &faces_distant);

  for (cs_lnum_t i = 0; i < n_distant; i++) {
    cs_lnum_t face_id = faces_local[i];
    isothm[face_id] = CS_BOUNDARY_RAD_WALL_GRAY;
    xlamp[face_id] = -cs_math_big_r;
    textp[face_id] = -cs_math_big_r;
    if (beps[face_id] < 0)
      have_unset = true;
  }

  if (have_unset) {
    int *is_solid = NULL;

    if (n_distant > 0) {
      cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
      BFT_MALLOC(is_solid, n_cells_ext, int);
      for (cs_lnum_t i = 0; i < n_cells_ext; i++)
        is_solid[i] = 0;
      cs_volume_zone_tag_cell_type(CS_VOLUME_ZONE_SOLID, 1, is_solid);
    }

    for (cs_lnum_t i = 0; i < n_distant; i++) {
      cs_lnum_t face_id = faces_local[i];
      if (beps[face_id] < 0) {
        cs_lnum_t cell_id = m->b_face_cells[face_id];
        if (is_solid[cell_id])
          beps[face_id] = 0;
      }
    }

    BFT_FREE(is_solid);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief synchronize radiative boundary condition error logging across
 *        MPI ranks.
 *
 * \param[in, out] nerloc       number of errors (local rank in, global out)
 * \param[in]      nerrcd       number of codes saved at error faces
 * \param[in, out] znferr       zone number for one error face (local in,
 *                              broadcast out)
 * \param[in, out] rvferr       values saved at one error face (local in,
 *                              broadcast out)
 */
/*----------------------------------------------------------------------------*/

inline static void
_sync_rad_bc_err(cs_gnum_t  nerloc[],
                 int        nerrcd,
                 int       *znferr,
                 cs_real_t  rvferr[])
{
  if (cs_glob_rank_id >= 0) {

    int irkerr = -1;

    if (*nerloc > 0)
      irkerr  = cs_glob_rank_id;

    cs_parall_sum(1, CS_GNUM_TYPE, nerloc);

    if (*nerloc != 0) {
      cs_parall_max(1, CS_INT_TYPE, &irkerr);


      if (rvferr != NULL) {
        cs_parall_bcast(irkerr, nerrcd, CS_DOUBLE, rvferr);
        cs_parall_bcast(irkerr, 1, CS_INT_TYPE, znferr);
      }
      else
        cs_parall_bcast(irkerr, nerrcd, CS_INT_TYPE, znferr);

    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute wall temperature for radiative transfer, and update BCs.
 *
 *   1) Compute wall temperature for radiative transfer
 *
 *   2) Update BCs for the energy computation
 *
 *   \param[in]     bc_type       face boundary condition type
 *----------------------------------------------------------------------------*/

void
cs_rad_transfer_bcs(int bc_type[])
{
  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;
  const cs_time_step_t *ts = cs_glob_time_step;

  /* By pass BCs if time step is not active */
  bool is_active = cs_time_control_is_active(&(rt_params->time_control),ts);

  cs_real_t *dt = CS_F_(dt)->val;

  cs_real_t stephn = cs_physical_constants_stephan;

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  /* Allocate temporary arrays */

  int  *isothm;
  cs_real_t *tempk, *text, *twall;
  BFT_MALLOC(isothm, n_b_faces, int);
  BFT_MALLOC(tempk, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(text, n_b_faces, cs_real_t);
  BFT_MALLOC(twall, n_b_faces, cs_real_t);

  /* Map field arrays */
  cs_field_t *f_b_temp = cs_field_by_name_try("boundary_temperature");
  cs_field_t *f_bqinci = cs_field_by_name_try("rad_incident_flux");
  cs_field_t *f_bxlam  = cs_field_by_name_try("wall_thermal_conductivity");
  cs_field_t *f_bepa   = cs_field_by_name_try("wall_thickness");
  cs_field_t *f_beps   = cs_field_by_name_try("emissivity");

  /* Call counter  */
  ipacli++;

  /* Indicator: if not restart and first time step */
  int ideb = 0;

  /* --> Min and Max values of the temperature (in Kelvin)   */

  cs_real_t tmin = 0.0;
  cs_real_t tmax = cs_math_big_r + tkelvi;

  /* Relaxation coefficient: 0 < tx <= 1

   * To compute the wall temperature, compute a temperature increment
   * DeltaT between the current step n et previous step n-1, then compute:
   *    n    n-1                                 n-1
   *   T  = T    + DeltaT if the ratio  DeltaT/T    =< tx,
   *   otherwise:
   *    n    n-1                      n-1             n-1
   *   T  = T    * (1 + tx *((DeltaT/T   ) / |DeltaT/T   |))
   */

  cs_real_t tx = 0.1;

  /* Temperature scale */
  cs_real_t xmtk = 0.;
  if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
    xmtk = tkelvi;

  /* Wall temperature */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (   bc_type[face_id] == CS_SMOOTHWALL
        || bc_type[face_id] == CS_ROUGHWALL)
      twall[face_id] = f_b_temp->val[face_id] + xmtk;
    else
      twall[face_id] = 0.0;
  }

  /* Default initialization */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    isothm[face_id] = -1;
    f_bxlam->val[face_id] = -cs_math_big_r;
    f_bepa->val[face_id]  = -cs_math_big_r;
    f_beps->val[face_id]  = -cs_math_big_r;
    text[face_id] = -cs_math_big_r;
  }

  /* Index of the thermal variable */

  cs_field_t *fth = cs_thermal_model_field();
  int *th_icodcl = fth->bc_coeffs->icodcl;
  cs_real_t *th_rcodcl1 = fth->bc_coeffs->rcodcl1;
  cs_real_t *th_rcodcl2 = fth->bc_coeffs->rcodcl2;
  cs_real_t *th_rcodcl3 = fth->bc_coeffs->rcodcl3;

  cs_field_t *f_hgas = cs_field_by_name_try("x_c_h");

  int *hg_icodcl = NULL;
  cs_real_t *hg_rcodcl1 = NULL;
  cs_real_t *hg_rcodcl2 = NULL;
  cs_real_t *hg_rcodcl3 = NULL;

  if (f_hgas != NULL) {
    hg_icodcl = f_hgas->bc_coeffs->icodcl;
    hg_rcodcl1 = f_hgas->bc_coeffs->rcodcl1;
    hg_rcodcl2 = f_hgas->bc_coeffs->rcodcl2;
    hg_rcodcl3 = f_hgas->bc_coeffs->rcodcl3;
  }

  /* Pointers to specific fields    */
  cs_field_t *f_bfconv = cs_field_by_name_try("rad_convective_flux");
  cs_field_t *f_bhconv = cs_field_by_name_try("rad_exchange_coefficient");

  /* If no restart info is available, then initialization at first pass,
   * for qincid:
   *   qincid <= stephn*tint**4
   * (if qincid is set to zero, there will be a deficit on the radiance
   *  BC at the first time step with DOM). */

  if (   ipacli == 1 && is_active
      && cs_glob_rad_transfer_params->restart == false) {

    /* If not a restart and first time step. */
    ideb = 1;

    cs_real_t *rad_st_impl = CS_FI_(rad_ist, 0)->val;
    cs_real_t *rad_st_expl = CS_FI_(rad_est, 0)->val;

    cs_array_real_fill_zero(cs_glob_mesh->n_cells_with_ghosts, rad_st_impl);
    cs_array_real_fill_zero(cs_glob_mesh->n_cells_with_ghosts, rad_st_expl);

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      f_bhconv->val[face_id] = 0.0;
      f_bfconv->val[face_id] = 0.0;
    }

    /* User definitions */

    cs_gui_radiative_transfer_bcs(bc_type,
                                  isothm,
                                  f_beps->val,
                                  f_bepa->val,
                                  text,
                                  f_bxlam->val);

    cs_user_radiative_transfer_bcs(cs_glob_domain,
                                   bc_type,
                                   isothm,
                                   &tmin,
                                   &tmax,
                                   &tx,
                                   dt,
                                   twall,
                                   f_bqinci->val,
                                   f_bhconv->val,
                                   f_bfconv->val,
                                   f_bxlam->val,
                                   f_bepa->val,
                                   f_beps->val,
                                   text);

    //FIXME we do the contrary of the message...
    cs_log_printf(CS_LOG_DEFAULT,
                  _("\n"
                    "   ** Information on the radiative module\n"
                    "      -----------------------------------\n"
                    "    Initialization of the wall temperature\n"
                    "    with incident flux at walls (qincid).\n"));

    /* twall in Kelvin and qincid in W/m2  */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   bc_type[face_id] == CS_SMOOTHWALL
          || bc_type[face_id] == CS_ROUGHWALL) {
        f_bqinci->val[face_id] = stephn * cs_math_pow4(twall[face_id]);
      }
      else {
        twall[face_id]         = 0.0;
        f_bqinci->val[face_id] = 0.0;
      }
    }
  }

  /* Values for boundary faces */

  cs_gui_radiative_transfer_bcs(bc_type,
                                isothm,
                                f_beps->val,
                                f_bepa->val,
                                text,
                                f_bxlam->val);

  cs_user_radiative_transfer_bcs(cs_glob_domain,
                                 bc_type,
                                 isothm,
                                 &tmin,
                                 &tmax,
                                 &tx,
                                 dt,
                                 twall,
                                 f_bqinci->val,
                                 f_bhconv->val,
                                 f_bfconv->val,
                                 f_bxlam->val,
                                 f_bepa->val,
                                 f_beps->val,
                                 text);

  /* Internal coupling settings */
  if (cs_internal_coupling_n_couplings() > 0) {
    cs_internal_coupling_t *cpl = NULL;

    cs_field_t *tf = cs_thermal_model_field();
    if (tf != NULL) {
      const int coupling_key_id = cs_field_key_id("coupling_entity");
      int coupling_id = cs_field_get_key_int(tf, coupling_key_id);
      if (coupling_id >= 0)
        cpl = cs_internal_coupling_by_id(coupling_id);
    }

    if (cpl != NULL)
      _set_internal_coupling_bcs(cpl,
                                 isothm,
                                 f_bxlam->val,
                                 f_beps->val,
                                 text);
  }

  /* Check user BC definitions */

  const int *face_zone_id = cs_boundary_zone_face_zone_id();

  /* Error counter */
  cs_gnum_t nrferr[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int       icoerr[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  cs_real_t rvferr[27];

  {
    /* Error if isothm not defined on wall, or defined on non-wall */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   bc_type[face_id] == CS_SMOOTHWALL
          || bc_type[face_id] == CS_ROUGHWALL) {
        if (isothm[face_id] == -1) {
          nrferr[2]++;
          icoerr[2]    = face_zone_id[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }
      else {
        if (   isothm[face_id] != -1
            && isothm[face_id] != cs_glob_rad_transfer_params->ifinfe) {
          nrferr[3]++;
          icoerr[3]    =  face_zone_id[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }
    }

    /* If incorrect physical value given, error */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      int rad_bc_code = isothm[face_id];

      if (rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0) {
          nrferr[4]++;
          icoerr[4]    = face_zone_id[face_id];
          rvferr[0]    = f_beps->val[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0
            || f_bxlam->val[face_id] <= 0.0
            || f_bepa->val[face_id] <= 0.0
            || text[face_id] <= 0.0) {
          nrferr[5]++;
          icoerr[5]    = face_zone_id[face_id];
          rvferr[2]    = f_beps->val[face_id];
          rvferr[3]    = f_bxlam->val[face_id];
          rvferr[4]    = f_bepa->val[face_id];
          rvferr[5]    = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T) {
        if (   f_bxlam->val[face_id] <= 0.0
            || f_bepa->val[face_id] <= 0.0
            || text[face_id] <= 0.0) {
          nrferr[6]++;
          icoerr[6]    = face_zone_id[face_id];
          rvferr[7]    = f_bxlam->val[face_id];
          rvferr[8]    = f_bepa->val[face_id];
          rvferr[9]    = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0) {
          nrferr[7]++;
          icoerr[7]    = face_zone_id[face_id];
          rvferr[11]   = f_beps->val[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY_1D_T) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0) {
          nrferr[13]++;
          icoerr[15]   = face_zone_id[face_id];
          rvferr[25]   = f_beps->val[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (   rad_bc_code !=  -1
               && rad_bc_code != CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX
               && rad_bc_code != cs_glob_rad_transfer_params->ifinfe) {
        nrferr[8]++;
        icoerr[9]    = face_zone_id[face_id];
        icoerr[10]   = rad_bc_code;
        bc_type[face_id] = -CS_ABS(bc_type[face_id]);
      }
    }

    /* If value defined for no reason, error */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY) {
        if (   f_bxlam->val[face_id] > 0.0
            || f_bepa->val[face_id] > 0.0
            || text[face_id] > 0.0) {
          nrferr[9]++;
          icoerr[11]   = face_zone_id[face_id];
          rvferr[14]   = f_bxlam->val[face_id];
          rvferr[15]   = f_bepa->val[face_id];
          rvferr[16]   = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T) {
        if (f_beps->val[face_id] > 0.0) {
          nrferr[10]++;
          icoerr[12]   = face_zone_id[face_id];
          rvferr[17]   = f_beps->val[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
        f_beps->val[face_id] = 0.;
      }

      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
        if (   f_bxlam->val[face_id] > 0.0
            || f_bepa->val[face_id] > 0.0
            || text[face_id] > 0.0) {
          nrferr[11]++;
          icoerr[13]   = face_zone_id[face_id];
          rvferr[18]   = f_bxlam->val[face_id];
          rvferr[19]   = f_bepa->val[face_id];
          rvferr[20]   = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX) {
        if (   f_beps->val[face_id] > 0.0
            || f_bxlam->val[face_id] > 0.0
            || f_bepa->val[face_id] > 0.0
            || text[face_id] > 0.0) {
          nrferr[12]++;
          icoerr[14]   = face_zone_id[face_id];
          rvferr[21]   = f_beps->val[face_id];
          rvferr[22]   = f_bxlam->val[face_id];
          rvferr[23]   = f_bepa->val[face_id];
          rvferr[24]   = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
        f_beps->val[face_id] = 0.;
      }

    }

  }

  /* Error logging */

  int iok = 0;
  for (int ii = 0; ii < 14; ii++) {
    if (nrferr[ii] > 0){
      iok = 1;
      break;
    }
  }

  if (cs_glob_rank_id >= 0)
    cs_parall_max(1, CS_INT_TYPE, &iok);

  if (iok != 0) {

    _sync_rad_bc_err(&nrferr[2], 1, &icoerr[2], NULL);
    if (nrferr[2] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "isothp must be defined on all wall faces.\n"
                      "It has not been defined for %llu faces.\n\n"
                      "  last undefined face value:  %d\n"),
                    _("Radiative boundary conditions errors:\n"),
                    (unsigned long long)nrferr[2],
                    icoerr[2]);
    }

    _sync_rad_bc_err(&nrferr[3], 1, &icoerr[3], NULL);
    if (nrferr[3] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "isothp must be defined only on wall faces.\n"
                      "It has been defined for %llu non-wall faces.\n\n"
                      "  last such face value:  %d\n"),
                    _("Radiative boundary conditions errors:\n"),
                    (unsigned long long)nrferr[3],
                    icoerr[3]);
    }

    _sync_rad_bc_err(&nrferr[4], 2, &icoerr[4], &rvferr[0]);
    if (nrferr[4] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "CS_BOUNDARY_RAD_WALL_GRAY",
                    (unsigned long long)nrferr[4],
                    icoerr[4],
                    rvferr[0]);
    }

    _sync_rad_bc_err(&nrferr[5], 5, &icoerr[5], &rvferr[2]);
    if (nrferr[5] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "  XLAMP, EPAP, and TEXTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T",
                    (unsigned long long)nrferr[5],
                    icoerr[5],
                    rvferr[2],
                    rvferr[3],
                    rvferr[4],
                    rvferr[5]);
    }

    _sync_rad_bc_err (&nrferr[6], 4, &icoerr[6], &rvferr[7]);
    if (nrferr[6] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  XLAMP, EPAP and TEXTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T",
                    (unsigned long long)nrferr[6],
                    icoerr[6],
                    rvferr[7],
                    rvferr[8],
                    rvferr[9]);
    }

    _sync_rad_bc_err(&nrferr[7], 2, &icoerr[7], &rvferr[11]);
    if (nrferr[7] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX",
                    (unsigned long long)nrferr[7],
                    icoerr[7],
                    rvferr[11]);
    }

    _sync_rad_bc_err(&nrferr[13], 2, &icoerr[14], &rvferr[25]);
    if (nrferr[13] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"),
                    _("Radiative boundary conditions errors"),
                    "CS_BOUNDARY_RAD_1D_WALL_T",
                    (unsigned long long)nrferr[13],
                    icoerr[15],
                    rvferr[25]);
    }

    _sync_rad_bc_err(&nrferr[8], 2, &icoerr[9], NULL);
    if (nrferr[8] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "Forbidden value of ISOTHM for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    ISOTHM = %d\n"),
                    _("Radiative boundary conditions errors"),
                    (unsigned long long)nrferr[8],
                    icoerr[9],
                    icoerr[10]);
    }

    _sync_rad_bc_err(&nrferr[9], 3, &icoerr[11], &rvferr[14]);
    if (nrferr[9] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  XLAMP, EPAP and TEXTP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors"),
                    "CS_BOUNDARY_RAD_WALL_GRAY",
                    (unsigned long long)nrferr[9],
                    icoerr[11],
                    rvferr[14],
                    rvferr[15],
                    rvferr[16]);
    }

    _sync_rad_bc_err(&nrferr[10], 1, &icoerr[12], &rvferr[17]);
    if (nrferr[10] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  EPSP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP = %12.4e\n"),
                    _("Radiative boundary conditions errors"),
                    "CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T",
                    (unsigned long long)nrferr[10],
                    icoerr[12],
                    rvferr[17]);
    }

    _sync_rad_bc_err(&nrferr[11], 3, &icoerr[13], &rvferr[18]);
    if (nrferr[11] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  XLAMP, EPAP and TEXTP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors"),
                    "CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX",
                    (unsigned long long)nrferr[11],
                    icoerr[13],
                    rvferr[18],
                    rvferr[19],
                    rvferr[20]);
    }

    _sync_rad_bc_err(&nrferr[12], 4, &icoerr[14], &rvferr[21]);
    if (nrferr[12] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n%s\n\n"
                      "With isothp = %s\n"
                      "  EPSP, XLAMP, EPAP and TEXTP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors"),
                    "CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX",
                    (unsigned long long)nrferr[12],
                    icoerr[14],
                    rvferr[21],
                    rvferr[22],
                    rvferr[23],
                    rvferr[24]);
    }

    cs_boundary_conditions_error(bc_type, NULL);

  }

  /* Complete user definitions: icodcl and EPS (when zero) */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    int icodw;
    if (bc_type[face_id] == CS_ROUGHWALL)
      icodw = 6;
    else
      icodw = 5;

    int rad_bc_code = isothm[face_id];

    if (   rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY
        || rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY_1D_T) {
      int t_bc_code = th_icodcl[face_id];
      /* Negativ icodcl if no conversion is needed
       * (BCs directly expressed in term of solved variable) */
      int sgn = (t_bc_code < 0) ? - 1 : 1;
      if (sgn*t_bc_code != icodw) {
        if (th_icodcl[face_id] == 15)
          th_icodcl[face_id] = 15*sgn;
        else
          th_icodcl[face_id] = icodw*sgn;
      }
      if (f_hgas != NULL)
        hg_icodcl[face_id] = icodw;
    }

    else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T) {
      th_icodcl[face_id]   = icodw;
      if (f_hgas != NULL)
        hg_icodcl[face_id] = icodw;
    }

    else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T) {
      th_icodcl[face_id]   = icodw;
      if (f_hgas != NULL)
        hg_icodcl[face_id] = icodw;
      f_beps->val[face_id]     = 0.0;
    }

    else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
      th_icodcl[face_id]   = icodw;
      if (f_hgas != NULL)
        hg_icodcl[face_id] = icodw;
    }

    else if (rad_bc_code == CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX) {
      th_icodcl[face_id]   = 3;
      if (f_hgas != NULL)
        hg_icodcl[face_id] = 3;
      f_beps->val[face_id]     = 0.0;
    }

  }

  /* Save temperature (in Kelvin) in tempk */

  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE) {

    cs_field_t *f_temp = CS_F_(t);
    if (f_temp == NULL)
      f_temp = CS_FI_(t, 0);

    cs_real_t *cval_t = f_temp->val;
    if (f_temp->n_time_vals > 1)
      cval_t = f_temp->vals[1];

    if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS) {

      for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
        tempk[iel] = cval_t[iel] + tkelvi;

    }
    else if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_KELVIN) {

      /* val index to access, necessary for compatibility with neptune_cfd */
      int tval_id = f_temp->n_time_vals - 1;
      for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
        tempk[iel] = f_temp->vals[tval_id][iel];

    }

  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    cs_field_t *f_enthalpy = CS_F_(h);

    /* Results: T to K */
    cs_ht_convert_h_to_t_cells(f_enthalpy->val, tempk);

  }

  /* Compute wall temperature
     ------------------------
     In all cases hfconv contains Lambda * Hturb / distance
     (hfconv: W/(m2 K); Hturb is dimensionless)
     (at first pass, it is zero)
     -> Compute convective flux, by which we mean:
     convective flux parallel to wall; we assume the wall is watertight
     The flux is computed in cs_boundary_condition_set_coeffs/clptur, except at
     the first pass with no restart, as cs_rad_transfer_bcs is called first.
   */

  if (ideb == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (isothm[face_id] !=  -1)
        f_bfconv->val[face_id] *=   tempk[cs_glob_mesh->b_face_cells[face_id]]
                                  - twall[face_id];
    }
  }

  /* The cases where twall must be computed first, are at first pass without
     restart, cases with prescribed temperature twall */

  if (ideb == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T
          || isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T
          || isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX
          || isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_1D_T) {
        isothm[face_id] = CS_BOUNDARY_RAD_WALL_GRAY;
        th_rcodcl1[face_id] = twall[face_id] - xmtk;
      }
    }
  }

  if (ideb == 0)
    cs_rad_transfer_compute_wall_t(isothm,
                                   tmin,
                                   tmax,
                                   tx,
                                   f_bqinci->val,
                                   text,
                                   f_bxlam->val,
                                   f_bepa->val,
                                   f_beps->val,
                                   f_bhconv->val,
                                   f_bfconv->val,
                                   tempk,
                                   twall);

  /* Change user boundary conditions */

  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE) {

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T) {
        th_rcodcl1[face_id] = text[face_id] - xmtk;
        th_rcodcl2[face_id] = f_bxlam->val[face_id] / f_bepa->val[face_id];
        th_rcodcl3[face_id] = 0.0;
      }

      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX) {
        th_icodcl[face_id] = 3;
        th_rcodcl1[face_id] = 0.0;
        th_rcodcl2[face_id] = cs_math_infinite_r;
      }
      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
        /* Update wall temperature to be imposed */
        th_rcodcl1[face_id] = twall[face_id] - xmtk;
      }

    }

  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    /* Read user data;
     * convert twall to enthalpy at boundary */

    cs_lnum_t *lstfac;
    BFT_MALLOC(lstfac, n_b_faces, cs_lnum_t);

    cs_real_t *wall_enth = NULL, *ext_enth = NULL;
    BFT_MALLOC(wall_enth, n_b_faces, cs_real_t);
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      wall_enth[face_id] = 0.;

    cs_lnum_t nlst = 0;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY
          || isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T
          || isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
        if (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL) {
          lstfac[nlst] = face_id;
          nlst++;
        }
      }
    }
    if (nlst > 0)
      cs_ht_convert_t_to_h_faces_l(nlst, lstfac, twall, wall_enth);

    nlst = 0;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T) {
        if (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL) {
          lstfac[nlst] = face_id;
          nlst++;
        }
      }
    }

    if (nlst > 0) {
      BFT_MALLOC(ext_enth, n_b_faces, cs_real_t);
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
        ext_enth[face_id] = 0.;

      cs_ht_convert_t_to_h_faces_l(nlst, lstfac, text, ext_enth);
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY) {
        if (f_hgas != NULL) {
          hg_rcodcl1[face_id] = wall_enth[face_id];
          hg_rcodcl2[face_id] = cs_math_infinite_r;
          hg_rcodcl3[face_id] = 0.0;
        }
      }
      else if (  isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T
              || isothm[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
        /* use twall instead of wall_enth for thermal scalar
         * to avoid extra conversions */
        th_icodcl[face_id] *= -1;
        th_rcodcl1[face_id] = twall[face_id] - xmtk;
        th_rcodcl2[face_id] = cs_math_infinite_r;
        th_rcodcl3[face_id] = 0.0;

        if (f_hgas != NULL) {
          hg_rcodcl1[face_id] = wall_enth[face_id];
          hg_rcodcl2[face_id] = cs_math_infinite_r;
          hg_rcodcl3[face_id] = 0.0;
        }
      }
      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T) {
        th_rcodcl1[face_id] = ext_enth[face_id];
        /* hext  */
        th_rcodcl2[face_id] = f_bxlam->val[face_id] / f_bepa->val[face_id];
        th_rcodcl3[face_id] = 0.0;
        if (f_hgas != NULL) {
          hg_rcodcl1[face_id] = ext_enth[face_id];
          hg_rcodcl2[face_id] = f_bxlam->val[face_id] / f_bepa->val[face_id];
          hg_rcodcl3[face_id] = 0.0;
        }
      }
      else if (isothm[face_id] == CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX) {
        th_icodcl[face_id] = 3;
        th_rcodcl1[face_id] = 0.0;
        th_rcodcl2[face_id] = cs_math_infinite_r;
        if (f_hgas != NULL) {
          hg_icodcl[face_id] = 3;
          hg_rcodcl1[face_id] = 0.0;
          hg_rcodcl2[face_id] = cs_math_infinite_r;
        }
      }
    }

    BFT_FREE(lstfac);
    BFT_FREE(ext_enth);
    BFT_FREE(wall_enth);
  }

  /* Update boundary temperature field   */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (   bc_type[face_id] == CS_SMOOTHWALL
        || bc_type[face_id] == CS_ROUGHWALL) {
      //FIXME check if this is useful; or simply updating rcodcl instead ?
      f_b_temp->val[face_id] = twall[face_id] - xmtk;
    }

    else {
      /* We use f_beps values as an indicator for open boundary types
       * (as this value si otherwise onlyused for wall conditons and we do not
       * yet have "revised" radiative boundary condition types */
      if (   cs_glob_rad_transfer_params->atmo_model & CS_RAD_ATMO_3D_INFRARED
          || isothm[face_id] == cs_glob_rad_transfer_params->ifinfe)
        f_beps->val[face_id] = 1.0;
    }
  }

  /* Free memory */

  BFT_FREE(isothm);
  BFT_FREE(tempk);
  BFT_FREE(text);
  BFT_FREE(twall);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Boundary conditions for DO and P-1 models
 *
 *  The coefap array stores the intensity for each boundary face,
 *  depending of the nature of the boundary (Dirichlet condition).
 *  The intensity of radiation is defined as the rate of emitted
 *  energy from unit surface area through a unit solid angle.
 *
 *   1/ Gray wall: isotropic radiation field.
 *
 *   \f$ coefap =  \epsilon.\sigma.twall^4 / \pi + (1-\epsilon).qincid / \pi \f$
 *
 *   which is the sum of the wall emission and reflecting flux
 *   (eps=1: black wall; eps=0: reflecting wall).
 *
 *   2/ Free boundary: condition to mimic infinite domain
 *
 * \param[in]  bc_type         boundary face types
 * \param[in]  vect_s          direction vector or NULL
 * \param[in]  ckmel           Absoprtion coefficcient of the mixture
 *                             gas-particules of coal or NULL
 * \param[in]  bpro_eps        Boundary emissivity, or NULL for solar radiation
 * \param[in]  w_gg            Weights of the i-th gray gas at boundaries
 * \param[in]  gg_id           number of the i-th grey gas
 * \param[out] coefap          boundary conditions for intensity or P-1 model
 * \param[out] coefbp          boundary conditions for intensity or P-1 model
 * \param[out] cofafp          boundary conditions for intensity or P-1 model
 * \param[out] cofbfp          boundary conditions for intensity or P-1 model
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_bc_coeffs(int        bc_type[],
                          cs_real_t  vect_s[3],
                          cs_real_t  ckmel[],
                          cs_real_t  bpro_eps[],
                          cs_real_t  w_gg[],
                          int        gg_id,
                          cs_real_t  coefap[],
                          cs_real_t  coefbp[],
                          cs_real_t  cofafp[],
                          cs_real_t  cofbfp[])
{
  cs_real_t stephn = cs_physical_constants_stephan;
  cs_real_t onedpi  = 1.0 / cs_math_pi;

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  cs_real_3_t *b_face_normal
    = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;

  /* Initialization */

  /* Pointer to the spectral flux density field */
  cs_field_t *f_qinspe = cs_field_by_name_try("spectral_rad_incident_flux");

  cs_real_t *q_incid = NULL;
  cs_lnum_t stride = 1;
  if (f_qinspe != NULL) {
    q_incid = f_qinspe->val;
    stride = f_qinspe->dim;
  }
  else
    q_incid = cs_field_by_name("rad_incident_flux")->val;


  /* Pointer to the wall emissivity field */
  cs_field_t *f_eps = cs_field_by_name("emissivity");

  /* Wall temperature */
  cs_field_t *f_tempb = CS_F_(t_b);
  cs_real_t xmtk = 0;
  if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
    xmtk = cs_physical_constants_celsius_to_kelvin;

  /* -> Initialization to a non-admissible value for testing after raycll   */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
    coefap[face_id] = -cs_math_big_r;

  /* Boundary conditions for DO model
   * coefap must be filled with the intensity */

  cs_real_t qpatmp;

  if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_DOM) {

    const cs_real_t *grav = cs_glob_physical_constants->gravity;
    const cs_real_t g_norm = cs_math_3_norm(grav);
    const cs_real_t d_g = (g_norm > 0.) ? 1./g_norm : 0.;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      /* Copy the appropriate flux density to the local variable qpatmp*/

      /* Value of the flux density at the boundary face */
      qpatmp = q_incid[face_id * stride + gg_id];

      /* Dirichlet Boundary Conditions  */

      cs_real_t hint = 0.0;

      /* Open boundary conditions */

      if (   bc_type[face_id] == CS_INLET
          || bc_type[face_id] == CS_CONVECTIVE_INLET
          || bc_type[face_id] == CS_OUTLET
          || bc_type[face_id] == CS_FREE_INLET
          || bc_type[face_id] == CS_SYMMETRY) {

        /* Legacy open boundary conditions if (eps < 0)
           TODO use boundary definitions from cs_boundary.h instead
           of this temporary hack. */

        if (f_eps->val[face_id] < 0.5) {

          cs_real_t pimp;
          /* Symmetry: pseudo-reflecting boundary conditions
           * true symmetry would require coupling directions */
          if (bc_type[face_id] == CS_SYMMETRY)
            pimp  = qpatmp * onedpi;
          /* Inlet/Outlet face: entering intensity fixed to zero
           * (warning: the treatment is different from than of P-1 model) */
          else
            pimp  = cs_math_epzero;

          cs_boundary_conditions_set_dirichlet_scalar(&coefap[face_id],
                                                      &cofafp[face_id],
                                                      &coefbp[face_id],
                                                      &cofbfp[face_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);
          continue;
        }

        /* Inlet/Outlet face: entering intensity fixed to zero flux
         * which mimics an infinite extrusion
         * (warning: the treatment is different from than of P-1 model)
         * Symmetry: the reflecting boundary conditions (eps=0)
         * is not feasible (because that would required to couple directions)
         * So we impose a homogeneous Neumann that models an infinite
         * extrusion */
        bool neumann = true;

        /* Entering intensity fixed to zero
         * if the normal to the face is colinear to the direction
         * or if coupled to a 1D module (Atmospheric module)
         * (warning: the treatment is different from than of P-1 model) */
        if (vect_s != NULL) {
          cs_real_t normal[3];
          cs_math_3_normalize(b_face_normal[face_id], normal);
          cs_real_t vs_dot_n = cs_math_3_dot_product(vect_s, normal);
          if (CS_ABS(vs_dot_n) < cs_math_epzero) /* entering */
            neumann = false;

          /* Top of the domain */
          if (cs_glob_rad_transfer_params->atmo_model != CS_RAD_ATMO_3D_NONE) {
            cs_real_t g_dot_n_norm = cs_math_3_dot_product(grav, normal) * d_g;
            if (g_dot_n_norm < -0.5)
              neumann = false;
          }
        }

        if (neumann) {
          cs_real_t qimp  = 0.;
          hint = 1.;

          cs_boundary_conditions_set_neumann_scalar(&coefap[face_id],
                                                    &cofafp[face_id],
                                                    &coefbp[face_id],
                                                    &cofbfp[face_id],
                                                    qimp,
                                                    hint);
        }
        else {
          /* Get boundary values:
           * for atmospheric, get the one computed by the 1D model */
          cs_real_t pimp = qpatmp * onedpi;

          cs_boundary_conditions_set_dirichlet_scalar(&coefap[face_id],
                                                      &cofafp[face_id],
                                                      &coefbp[face_id],
                                                      &cofbfp[face_id],
                                                      pimp,
                                                      hint,
                                                      cs_math_infinite_r);

        }

      }

      /* Wall boundary face: calculated intensity */
      else if (   bc_type[face_id] == CS_SMOOTHWALL
               || bc_type[face_id] == CS_ROUGHWALL) {
        cs_real_t twall = f_tempb->val[face_id] + xmtk;
        /* Remember: In case of the usage of the standard radiation
           models of code_saturne, w_gg=1  */
        cs_real_t pimp;
        if (bpro_eps != NULL)
          pimp = bpro_eps[face_id] * stephn * cs_math_pow4(twall)
                 * onedpi * w_gg[face_id + gg_id * n_b_faces]
               + (1.0 - bpro_eps[face_id]) * qpatmp * onedpi;
        /* For solar radiation,
         * Note that albedo is already taken into account
         * in the incident flux */
        else
          pimp = qpatmp * onedpi;

        cs_boundary_conditions_set_dirichlet_scalar(&coefap[face_id],
                                                    &cofafp[face_id],
                                                    &coefbp[face_id],
                                                    &cofbfp[face_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

      }
      /* Error if there are forgotten faces */
      else
        bc_type[face_id] = -CS_ABS(bc_type[face_id]);

    }

  }

  /* Boundary conditions for P-1 model:
   * coefap and coefbp must be filled */

  else if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_P1) {

    const cs_real_t *b_dist =  cs_glob_mesh_quantities->b_dist;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_lnum_t iel = cs_glob_mesh->b_face_cells[face_id];

      cs_real_t hint = 1.0 / (ckmel[iel] * b_dist[face_id]);

      /* Symmetry or reflecting wall (EPS = 0) : zero flux */
      if (   bc_type[face_id] == CS_SYMMETRY
          || (   (   bc_type[face_id] == CS_SMOOTHWALL
                  || bc_type[face_id] == CS_ROUGHWALL)
              && bpro_eps[face_id] <= 0.0) ) {
        cs_real_t qimp = 0.;
        cs_boundary_conditions_set_neumann_scalar(&coefap[face_id],
                                                  &cofafp[face_id],
                                                  &coefbp[face_id],
                                                  &cofbfp[face_id],
                                                  qimp,
                                                  hint);
      }

      /* Inlet/Outlet faces: zero flux
       * (warning: the treatment is different from than of DO model) */
      else if (   bc_type[face_id] == CS_INLET
               || bc_type[face_id] == CS_CONVECTIVE_INLET
               || bc_type[face_id] == CS_OUTLET
               || bc_type[face_id] == CS_FREE_INLET) {
        cs_real_t qimp = 0.;
        cs_boundary_conditions_set_neumann_scalar(&coefap[face_id],
                                                  &cofafp[face_id],
                                                  &coefbp[face_id],
                                                  &cofbfp[face_id],
                                                  qimp,
                                                  hint);
      }

      /*  Wall boundary faces */
      else if (   bc_type[face_id] == CS_SMOOTHWALL
               || bc_type[face_id] == CS_ROUGHWALL) {
        cs_real_t twall = f_tempb->val[face_id] + xmtk;
        cs_real_t distbf  = cs_glob_mesh_quantities->b_dist[face_id];
        cs_real_t xit = 1.5 * distbf * ckmel[iel]
                            * (2.0 / (2.0 - bpro_eps[face_id]) - 1.0);
        cs_real_t cfl = 1.0 / xit;
        cs_real_t pimp =   cs_math_pow4(twall)
                         * w_gg[face_id + gg_id * n_b_faces];
        cs_boundary_conditions_set_convective_outlet_scalar(&coefap[face_id],
                                                            &cofafp[face_id],
                                                            &coefbp[face_id],
                                                            &cofbfp[face_id],
                                                            pimp,
                                                            cfl,
                                                            hint);
      }

      /* 2.4 - Stop if there are forgotten faces
       *       ---------------------------------  */
      else
        bc_type[face_id] = -CS_ABS(bc_type[face_id]);
    }
  }

  cs_boundary_conditions_error(bc_type, NULL);

  /* --> Check radiance boundary conditions arrays
   *     In the case of the P-1 approximation, the value of coefap can be
   *     large (of the order of twall**4), hence the value of coefmn */

  cs_real_t xlimit = -cs_math_big_r * 0.1;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (coefap[face_id] <= xlimit)
      bc_type[face_id] = -CS_ABS(bc_type[face_id]);
  }

  cs_boundary_conditions_error(bc_type, "Radiance BC values");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
