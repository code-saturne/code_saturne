/*============================================================================
 * Radiation solver boundary conditions treatment.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_boundary_zone.h"
#include "cs_log.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_boundary_conditions.h"
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
 *   \param[in]     nvar          total number of variable BC's
 *   \param[in,out] icodcl        face boundary condition code:
 *                                 - 1 Dirichlet
 *                                 - 2 Radiative outlet
 *                                 - 3 Neumann
 *                                 - 4 sliding and
 *                                   \f$ \vect{u} \cdot \vect{n} = 0 \f$
 *                                 - 5 smooth wall and
 *                                   \f$ \vect{u} \cdot \vect{n} = 0 \f$
 *                                 - 6 rough wall and
 *                                   \f$ \vect{u} \cdot \vect{n} = 0 \f$
 *                                 - 9 free inlet/outlet
 *                                   (input mass flux blocked to 0)
 *                                 - 13 Dirichlet for the advection operator and
 *                                      Neumann for the diffusion operator
 *   \param[in]     bc_type       face boundary condition type
 *   \param[in]     dt            time step (per cell)
 *   \param[in,out] rcodcl        boundary condition values:
 *                                 - rcodcl(1) value of the dirichlet
 *                                 - rcodcl(2) value of the exterior exchange
 *                                   coefficient (infinite if no exchange)
 *                                 - rcodcl(3) value flux density
 *                                   (negative if gain) in w/m2 or roughness
 *                                   in m if icodcl=6
 *                                   -# for the velocity \f$ (\mu+\mu_T)
 *                                      \gradv \vect{u} \cdot \vect{n}  \f$
 *                                   -# for the pressure \f$ \Delta t
 *                                      \grad P \cdot \vect{n}  \f$
 *                                   -# for a scalar \f$ cp \left( K +
 *                                       \dfrac{K_T}{\sigma_T} \right)
 *                                       \grad T \cdot \vect{n} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_bcs(int         nvar,
                    int         bc_type[],
                    int         icodcl[],
                    cs_real_t   dt[],
                    cs_real_t   rcodcl[])
{
  cs_real_t stephn = cs_physical_constants_stephan;

  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  /* Allocate temporary arrays */

  cs_int_t  *isothm, *lstfac;
  cs_real_t *tempk, *thwall, *text, *tint, *twall;
  BFT_MALLOC(isothm, n_b_faces, cs_int_t);
  BFT_MALLOC(lstfac, n_b_faces, cs_int_t);
  BFT_MALLOC(tempk, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(thwall, n_b_faces, cs_real_t);
  BFT_MALLOC(text, n_b_faces, cs_real_t);
  BFT_MALLOC(tint, n_b_faces, cs_real_t);
  BFT_MALLOC(twall, n_b_faces, cs_real_t);

  /* Map field arrays     */
  cs_field_t *f_b_temp = cs_field_by_name_try("boundary_temperature");
  cs_field_t *f_bqinci = cs_field_by_name_try("rad_incident_flux");
  cs_field_t *f_bxlam  = cs_field_by_name_try("wall_thermal_conductivity");
  cs_field_t *f_bepa   = cs_field_by_name_try("wall_thickness");
  cs_field_t *f_beps   = cs_field_by_name_try("emissivity");
  cs_field_t *f_bfnet  = cs_field_by_name_try("rad_net_flux");

  /* Call counter  */
  ipacli++;

  /* Indicator: if not restart and first time step */
  int ideb = 0;

  /* --> Min and Max values of the temperature (in Kelvin)   */

  cs_real_t tmin = 0.0;
  cs_real_t tmax = cs_math_big_r + tkelvi;

  /* Relaxation coefficient: 0 < tx <= 1

   * To compute the wall temperature, compute a tmperature increment
   * DeltaT between the current step n et previous step n-1, then compute:
   *    n    n-1                                 n-1
   *   T  = T    + DeltaT si le rapport DeltaT/T    =< tx, sinon
   *    n    n-1                      n-1             n-1
   *   T  = T    * (1 + tx *((DeltaT/T   ) / |DeltaT/T   |))
   */

  cs_real_t tx = 0.1;

  /* Temperature scale */
  cs_real_t xmtk;
  if (cs_glob_thermal_model->itpscl == 2)
    xmtk = -tkelvi;
  else if (cs_glob_thermal_model->itpscl == 1)
    xmtk = 0.0;

  /* Wall temperature */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (   bc_type[face_id] == CS_SMOOTHWALL
        || bc_type[face_id] == CS_ROUGHWALL)
      twall[face_id] = f_b_temp->val[face_id] - xmtk;
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
    tint[face_id] = -cs_math_big_r;
  }

  /* Index of the thermal variable */

  cs_field_t *fth = cs_thermal_model_field();
  const cs_lnum_t ivart
    = cs_field_get_key_int(fth, cs_field_key_id("variable_id")) - 1;

  cs_field_t *f_hgas = cs_field_by_name_try("x_c_h");
  int ivahg = -1;
  if (f_hgas != NULL) {
    const int var_key_id = cs_field_key_id("variable_id");
    ivahg = cs_field_get_key_int(f_hgas, var_key_id);
  }

  /* Pointers to specific fields    */
  cs_field_t *f_bfconv = cs_field_by_name_try("rad_convective_flux");
  cs_field_t *f_bhconv = cs_field_by_name_try("rad_exchange_coefficient");

  /* If no restart info is available, then initialization at first pass,
   * for twall and qincid:
   *   twall <= tint
   *   qincid <= stephn*tint**4
   * (if qincid is et to zero, there will be a deficit on the luminance
   *  BC at the first time step with DOM). */

  if (   ipacli == 1
      && cs_glob_rad_transfer_params->restart == false) {

    /* If not a restart and first time step. */
    ideb = 1;

    cs_real_t *rad_st_impl = CS_FI_(rad_ist, 0)->val;
    cs_real_t *rad_st_expl = CS_FI_(rad_est, 0)->val;

    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells_with_ghosts; iel++) {
      rad_st_impl[iel] = 0.0;
      rad_st_expl[iel] = 0.0;
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      f_bhconv->val[face_id] = 0.0;
      f_bfconv->val[face_id] = 0.0;
    }

    /*        On utilise TBORD comme auxiliaire pour l'appel a USRAY2
     *          pour etre sur que TPAROI ne sera pas modifie
     *          (puisqu'on a TBORD libre)
     *        On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
     *          pour etre sur que QINCID ne sera pas modifie
     *          (puisqu'on a FLUNET libre) */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      thwall[face_id]        = 0.0;
      f_bfnet->val[face_id]  = 0.0;
    }

    /* User definitions */

    cs_gui_radiative_transfer_bcs(bc_type,
                                  nvar,
                                  ivart,
                                  isothm,
                                  f_beps->val,
                                  f_bepa->val,
                                  tint,
                                  text,
                                  f_bxlam->val,
                                  rcodcl);

    cs_user_radiative_transfer_bcs(nvar,
                                   bc_type,
                                   icodcl,
                                   isothm,
                                   &tmin,
                                   &tmax,
                                   &tx,
                                   dt,
                                   rcodcl,
                                   thwall,
                                   f_bfnet->val,
                                   f_bhconv->val,
                                   f_bfconv->val,
                                   f_bxlam->val,
                                   f_bepa->val,
                                   f_beps->val,
                                   text,
                                   tint);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("\n"
                    "   ** Information on the radiative module\n"
                    "      -----------------------------------\n"
                    "    Initialization of the wall temperature\n"
                    "    with user profile (tintp)\n"
                    "    and incident flux at walls (qincid).\n"));

    /* Tparoi en Kelvin et QINCID en W/m2  */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   bc_type[face_id] == CS_SMOOTHWALL
          || bc_type[face_id] == CS_ROUGHWALL) {
        twall[face_id]         = tint[face_id];
        f_bqinci->val[face_id] = stephn * cs_math_pow4(tint[face_id]);
      }
      else {
        twall[face_id]         = 0.0;
        f_bqinci->val[face_id] = 0.0;
      }
    }

  }

  /* Values for boundary faces */

  /* We use bfnet as an auxiliary fo the call to
     cs_user_radiative_transfer_bcs to make sure qincid is not modified
     (as bfnet is free) */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    thwall[face_id]       = twall[face_id];
    f_bfnet->val[face_id] = f_bqinci->val[face_id];
  }

  cs_gui_radiative_transfer_bcs(bc_type,
                                nvar,
                                ivart,
                                isothm,
                                f_beps->val,
                                f_bepa->val,
                                tint,
                                text,
                                f_bxlam->val,
                                rcodcl);

  cs_user_radiative_transfer_bcs(nvar,
                                 bc_type,
                                 icodcl,
                                 isothm,
                                 &tmin,
                                 &tmax,
                                 &tx,
                                 dt,
                                 rcodcl,
                                 thwall,
                                 f_bfnet->val,
                                 f_bhconv->val,
                                 f_bfconv->val,
                                 f_bxlam->val,
                                 f_bepa->val,
                                 f_beps->val,
                                 text,
                                 tint);

  /* Check user BC definitions */

  const int *face_zone_id = cs_boundary_zone_face_zone_id();

  /* Error counter */
  cs_gnum_t nrferr[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int       icoerr[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  cs_real_t rvferr[27];

  {
    /* Error if isothm not defined on wall, or defined on non-wall */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (  (   bc_type[face_id] == CS_SMOOTHWALL
             || bc_type[face_id] == CS_ROUGHWALL)
          && isothm[face_id] ==  -1) {
        nrferr[2]++;
        icoerr[2]    = face_zone_id[face_id];
        bc_type[face_id] = -CS_ABS(bc_type[face_id]);
      }
      else if (   bc_type[face_id] != CS_SMOOTHWALL
               && bc_type[face_id] != CS_ROUGHWALL
               && isothm[face_id] != -1) {
        nrferr[3]++;
        icoerr[3]    =  face_zone_id[face_id];
        bc_type[face_id] = -CS_ABS(bc_type[face_id]);
      }
    }

    /* If incorrect physical value given, error */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      if ( isothm[face_id] == cs_glob_rad_transfer_params->itpimp){
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0
            || tint[face_id] <= 0.0) {
          nrferr[4]++;
          icoerr[4]    = face_zone_id[face_id];
          rvferr[0]    = f_beps->val[face_id];
          rvferr[1]    = tint[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->ipgrno) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0
            || f_bxlam->val[face_id] <= 0.0
            || f_bepa->val[face_id] <= 0.0
            || text[face_id] <= 0.0
            || tint[face_id] <= 0.0) {
          nrferr[5]++;
          icoerr[5]    = face_zone_id[face_id];
          rvferr[2]    = f_beps->val[face_id];
          rvferr[3]    = f_bxlam->val[face_id];
          rvferr[4]    = f_bepa->val[face_id];
          rvferr[5]    = text[face_id];
          rvferr[6]    = tint[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->iprefl) {
        if (   f_bxlam->val[face_id] <= 0.0
            || f_bepa->val[face_id] <= 0.0
            || text[face_id] <= 0.0
            || tint[face_id] <= 0.0) {
          nrferr[6]++;
          icoerr[6]    = face_zone_id[face_id];
          rvferr[7]    = f_bxlam->val[face_id];
          rvferr[8]    = f_bepa->val[face_id];
          rvferr[9]    = text[face_id];
          rvferr[10]   = tint[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->ifgrno) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0
            || tint[face_id] <= 0.0) {
          nrferr[7]++;
          icoerr[7]    = face_zone_id[face_id];
          rvferr[11]   = f_beps->val[face_id];
          rvferr[12]   = tint[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->ifrefl) {
        if (tint[face_id] <= 0.0) {
          nrferr[8]++;
          icoerr[8]    = face_zone_id[face_id];
          rvferr[13]   = tint[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->itpt1d) {
        if (   f_beps->val[face_id] < 0.0
            || f_beps->val[face_id] > 1.0
            || tint[face_id] <= 0.0) {
          nrferr[14]++;
          icoerr[15]   = face_zone_id[face_id];
          rvferr[25]   = f_beps->val[face_id];
          rvferr[26]   = tint[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] !=  -1) {
        nrferr[9]++;
        icoerr[9]    = face_zone_id[face_id];
        icoerr[10]   = isothm[face_id];
        bc_type[face_id] = -CS_ABS(bc_type[face_id]);
      }
    }

    /* If value defined for no reason, error */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      if (isothm[face_id] == cs_glob_rad_transfer_params->itpimp) {
        if (   f_bxlam->val[face_id] > 0.0
            || f_bepa->val[face_id] > 0.0
            || text[face_id] > 0.0) {
          nrferr[10]++;
          icoerr[11]   = face_zone_id[face_id];
          rvferr[14]   = f_bxlam->val[face_id];
          rvferr[15]   = f_bepa->val[face_id];
          rvferr[16]   = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->iprefl) {
        if (f_beps->val[face_id] >= 0.0) {
          nrferr[11]++;
          icoerr[12]   = face_zone_id[face_id];
          rvferr[17]   = f_beps->val[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->ifgrno) {
        if (   f_bxlam->val[face_id] > 0.0
            || f_bepa->val[face_id] > 0.0
            || text[face_id] > 0.0) {
          nrferr[12]++;
          icoerr[13]   = face_zone_id[face_id];
          rvferr[18]   = f_bxlam->val[face_id];
          rvferr[19]   = f_bepa->val[face_id];
          rvferr[20]   = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->ifrefl) {
        if (   f_beps->val[face_id] >= 0.0
            || f_bxlam->val[face_id] > 0.0
            || f_bepa->val[face_id] > 0.0
            || text[face_id] > 0.0) {
          nrferr[13]++;
          icoerr[14]   = face_zone_id[face_id];
          rvferr[21]   = f_beps->val[face_id];
          rvferr[22]   = f_bxlam->val[face_id];
          rvferr[23]   = f_bepa->val[face_id];
          rvferr[24]   = text[face_id];
          bc_type[face_id] = -CS_ABS(bc_type[face_id]);
        }
      }

    }

  }

  /* Error logging */

  int iok = 0;
  for (int ii = 0; ii < 15; ii++) {
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
                    _("%s:\n"
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
                    _("%s:\n"
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
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "  TINTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    TINTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "ITPIMP",
                    (unsigned long long)nrferr[4],
                    icoerr[4],
                    rvferr[0],
                    rvferr[1]);
    }

    _sync_rad_bc_err(&nrferr[5], 5, &icoerr[5], &rvferr[2]);
    if (nrferr[5] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "  XLAMP, EPAP, TEXTP, and TINTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"
                      "    TINTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IPGRNO",
                    (unsigned long long)nrferr[5],
                    icoerr[5],
                    rvferr[2],
                    rvferr[3],
                    rvferr[4],
                    rvferr[5],
                    rvferr[6]);
    }

    _sync_rad_bc_err (&nrferr[6], 4, &icoerr[6], &rvferr[7]);
    if (nrferr[6] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  XLAMP, EPAP, TEXTP, and TINTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"
                      "    TINTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IPREFL",
                    (unsigned long long)nrferr[6],
                    icoerr[6],
                    rvferr[7],
                    rvferr[8],
                    rvferr[9],
                    rvferr[10]);
    }

    _sync_rad_bc_err(&nrferr[7], 2, &icoerr[7], &rvferr[11]);
    if (nrferr[7] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "  TINTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    TINTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IFGRNO",
                    (unsigned long long)nrferr[7],
                    icoerr[7],
                    rvferr[11],
                    rvferr[12]);
    }

    _sync_rad_bc_err(&nrferr[8], 1, &icoerr[8], &rvferr[13]);
    if (nrferr[8] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  TINTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    TINTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IFREFL",
                    (unsigned long long)nrferr[8],
                    icoerr[8],
                    rvferr[13]);
    }

    _sync_rad_bc_err(&nrferr[14], 2, &icoerr[14], &rvferr[25]);
    if (nrferr[14] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  EPSP  must be in range [0.; 1.]\n"
                      "  TINTP must be > 0\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    TINTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "ITPT1D",
                    (unsigned long long)nrferr[14],
                    icoerr[15],
                    rvferr[25],
                    rvferr[26]);
    }

    _sync_rad_bc_err(&nrferr[9], 2, &icoerr[9], NULL);
    if (nrferr[9] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "Forbidden value of ISOTHM for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    ISOTHM = %d\n"),
                    _("Radiative boundary conditions errors:\n"),
                    (unsigned long long)nrferr[9],
                    icoerr[9],
                    icoerr[10]);
    }

    _sync_rad_bc_err(&nrferr[10], 3, &icoerr[11], &rvferr[14]);
    if (nrferr[10] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  XLAMP, EPAP and TEXTP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "ITPIMP",
                    (unsigned long long)nrferr[10],
                    icoerr[11],
                    rvferr[14],
                    rvferr[15],
                    rvferr[16]);
    }

    _sync_rad_bc_err(&nrferr[11], 1, &icoerr[12], &rvferr[17]);
    if (nrferr[11] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  EPSP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IPREFL",
                    (unsigned long long)nrferr[11],
                    icoerr[12],
                    rvferr[17]);
    }

    _sync_rad_bc_err(&nrferr[12], 3, &icoerr[13], &rvferr[18]);
    if (nrferr[12] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  XLAMP, EPAP and TEXTP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IFGRNO",
                    (unsigned long long)nrferr[12],
                    icoerr[13],
                    rvferr[18],
                    rvferr[19],
                    rvferr[20]);
    }

    _sync_rad_bc_err(&nrferr[13], 4, &icoerr[14], &rvferr[21]);
    if (nrferr[13] > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%s:\n"
                      "With isothp = %s\n"
                      "  EPSP, XLAMP, EPAP and TEXTP must not be defined\n\n"
                      "This is not the case for %llu faces.\n\n"
                      "  last such face (zone  %d):\n"
                      "    EPSP  = %12.4e\n"
                      "    XLAMP = %12.4e\n"
                      "    EPAP  = %12.4e\n"
                      "    TEXTP = %12.4e\n"),
                    _("Radiative boundary conditions errors:\n"),
                    "IFREFL",
                    (unsigned long long)nrferr[13],
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

    if (isothm[face_id] == cs_glob_rad_transfer_params->itpimp){
      icodcl[ivart*n_b_faces + face_id]   = icodw;
      if (ivahg >= 0)
        icodcl[(ivahg - 1)*n_b_faces + face_id] = icodw;
    }

    else if (isothm[face_id] == cs_glob_rad_transfer_params->ipgrno){
      icodcl[ivart*n_b_faces + face_id]   = icodw;
      if (ivahg >= 0)
        icodcl[(ivahg - 1)*n_b_faces + face_id] = icodw;
    }

    else if(isothm[face_id] == cs_glob_rad_transfer_params->iprefl){
      icodcl[ivart*n_b_faces + face_id]   = icodw;
      if (ivahg >= 0)
        icodcl[(ivahg - 1)*n_b_faces + face_id] = icodw;
      f_beps->val[face_id]     = 0.0;
    }

    else if(isothm[face_id] == cs_glob_rad_transfer_params->ifgrno){
      icodcl[ivart*n_b_faces + face_id]   = icodw;
      if (ivahg >= 0)
        icodcl[(ivahg - 1)*n_b_faces + face_id] = icodw;
    }

    else if(isothm[face_id] == cs_glob_rad_transfer_params->ifrefl){
      icodcl[ivart*n_b_faces + face_id]   = 3;
      if (ivahg >= 0)
        icodcl[(ivahg - 1)*n_b_faces + face_id] = 3;
      f_beps->val[face_id]     = 0.0;
    }

    else if(isothm[face_id] == cs_glob_rad_transfer_params->itpt1d){
      icodcl[ivart*n_b_faces + face_id]   = icodw;
      if (ivahg >= 0)
        icodcl[(ivahg - 1)*n_b_faces + face_id] = icodw;
    }
  }

  /* Save temperature (in Kelvin) in tempk */

  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE) {

    cs_field_t *f_temp = CS_F_(t);
    if (f_temp == NULL)
      f_temp = CS_FI_(t,0);

    if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS) {

      for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
        tempk[iel] = f_temp->vals[1][iel] + tkelvi;

    }
    else if (cs_glob_thermal_model->itpscl == 1) {

      /* val index to access, necessary for compatibility with neptune */
      int tval_id = f_temp->n_time_vals - 1;
      for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
        tempk[iel] = f_temp->vals[tval_id][iel];

    }

  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    cs_field_t *f_enthalpy = CS_F_(h);

    /* Resultat : T en K    */
    CS_PROCF(c_h_to_t, C_H_TO_T)(f_enthalpy->val, tempk);

  }

  /* Compute wall temperature
     ------------------------
     In all cases hfconv contains Lambda * Hturb / distance
     (hfconv: W/(m2 K); Hturb is dimensionless)
     (at first pass, it is zero)
     -> Compute convective flux, by which we mean:
     convective flux parallel to wall; we assume the wall is watertight
     The flux is computed in dans condli/clptur, except at the first
     pass with no restart, as cs_rad_transfer_bcs is called first.
   */

  if (ideb == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (isothm[face_id] !=  -1)
        f_bfconv->val[face_id] *= tempk[cs_glob_mesh->b_face_cells[face_id]] - twall[face_id];
    }
  }

  /* The cases where twall must be computed first, are at first pass without
     restart, cases with prescribed temperature twall = tint */

  if (ideb == 1) {
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   isothm[face_id] == cs_glob_rad_transfer_params->ipgrno
          || isothm[face_id] == cs_glob_rad_transfer_params->iprefl
          || isothm[face_id] == cs_glob_rad_transfer_params->ifgrno
          || isothm[face_id] == cs_glob_rad_transfer_params->itpt1d)
        isothm[face_id] = cs_glob_rad_transfer_params->itpimp;
    }
  }

  if (ideb == 0)
    cs_rad_transfer_wall_flux(nvar,
                              ivart,
                              isothm,
                              &tmin,
                              &tmax,
                              &tx,
                              rcodcl,
                              twall,
                              f_bqinci->val,
                              text,
                              tint,
                              f_bxlam->val,
                              f_bepa->val,
                              f_beps->val,
                              f_bhconv->val,
                              f_bfconv->val,
                              tempk);

  /* Change user boundary conditions */

  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE) {

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      if (   isothm[face_id] == cs_glob_rad_transfer_params->itpimp
          || isothm[face_id] == cs_glob_rad_transfer_params->ipgrno
          || isothm[face_id] == cs_glob_rad_transfer_params->ifgrno) {
        rcodcl[0*n_b_faces*nvar + ivart*n_b_faces + face_id] = twall[face_id] + xmtk;
        rcodcl[1*n_b_faces*nvar + ivart*n_b_faces + face_id] = cs_math_infinite_r;
        rcodcl[2*n_b_faces*nvar + ivart*n_b_faces + face_id] = 0.0;
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->iprefl) {
        rcodcl[0*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = text[face_id] + xmtk;
        rcodcl[1*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = f_bxlam->val[face_id] / f_bepa->val[face_id];
        rcodcl[2*n_b_faces*nvar + ivart*n_b_faces + face_id] = 0.0;
      }

      else if (isothm[face_id] == cs_glob_rad_transfer_params->ifrefl) {
        icodcl[ivart*n_b_faces + face_id] = 3;
        rcodcl[0*n_b_faces*nvar + ivart*n_b_faces + face_id] = 0.0;
        rcodcl[1*n_b_faces*nvar + ivart*n_b_faces + face_id] = cs_math_infinite_r;
      }

    }

  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    /* Read user data;
     * convert twall to enthalpy at boundary, saved in flunet,
     * which is used as an auxiliary */

    int mode = 0;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   isothm[face_id] == cs_glob_rad_transfer_params->itpimp
          || isothm[face_id] == cs_glob_rad_transfer_params->ipgrno
          || isothm[face_id] == cs_glob_rad_transfer_params->ifgrno)
        mode  =  -1;
    }

    if (mode ==  -1) {
      cs_lnum_t nlst = 0;
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        if (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL) {
          lstfac[nlst] = face_id + 1; // for compatibility purpose with b_t_to_h
          nlst++;
        }
      }
      CS_PROCF(b_t_to_h, B_T_TO_H)(&nlst, lstfac, twall, f_bfnet->val);
    }

    mode = 0;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (isothm[face_id] == cs_glob_rad_transfer_params->iprefl) {
        mode = -1;
      }
    }

    if (mode == -1) {
      cs_lnum_t nlst = 0;
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        if (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL) {
          lstfac[nlst] = face_id + 1;// for compatibility purpose with b_t_to_h
          nlst++;
        }
      }
      CS_PROCF(b_t_to_h, B_T_TO_H) (&nlst, lstfac, text, thwall);
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   isothm[face_id] == cs_glob_rad_transfer_params->itpimp
          || isothm[face_id] == cs_glob_rad_transfer_params->ipgrno
          || isothm[face_id] == cs_glob_rad_transfer_params->ifgrno) {
        rcodcl[0*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = f_bfnet->val[face_id];
        rcodcl[1*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = cs_math_infinite_r;
        rcodcl[2*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = 0.0;
        if (ivahg >= 0) {
          rcodcl[0*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = f_bfnet->val[face_id];
          rcodcl[1*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = cs_math_infinite_r;
          rcodcl[2*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = 0.0;
        }
      }
      else if (isothm[face_id] == cs_glob_rad_transfer_params->iprefl) {
        rcodcl[0*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = thwall[face_id];
        /* hext  */
        rcodcl[1*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = f_bxlam->val[face_id] / f_bepa->val[face_id];
        rcodcl[2*n_b_faces*nvar + ivart*n_b_faces + face_id] = 0.0;
        if (ivahg >= 0) {
          rcodcl[0*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = thwall[face_id];
          rcodcl[1*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = f_bxlam->val[face_id] / f_bepa->val[face_id];
          rcodcl[2*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = 0.0;
        }
      }
      else if (isothm[face_id] == cs_glob_rad_transfer_params->ifrefl) {
        icodcl[ivart*n_b_faces + face_id] = 3;
        rcodcl[0*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = 0.0;
        rcodcl[1*n_b_faces*nvar + ivart*n_b_faces + face_id]
          = cs_math_infinite_r;
        if (ivahg >= 0) {
          icodcl[(ivahg - 1)*n_b_faces + face_id] = 3;
          rcodcl[0*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = 0.0;
          rcodcl[1*n_b_faces*nvar + (ivahg - 1)*n_b_faces + face_id]
            = cs_math_infinite_r;
        }
      }
    }
  }

  /* Update boundary temperature field   */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (   bc_type[face_id] == CS_SMOOTHWALL
        || bc_type[face_id] == CS_ROUGHWALL) {

      if (cs_glob_thermal_model->itpscl == 2)
        f_b_temp->val[face_id] = twall[face_id] - tkelvi;
      else
        f_b_temp->val[face_id] = twall[face_id];

    }

  }

  /* Free memory     */
  BFT_FREE(isothm);
  BFT_FREE(lstfac);
  BFT_FREE(tempk);
  BFT_FREE(thwall);
  BFT_FREE(text);
  BFT_FREE(tint);
  BFT_FREE(twall);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Boundary conditions for DO and P-1 models
 *
 * 1. Boundary conditions for the radiative intensity (DO model)
 * --------------------------------------------------------------
 *      The array coefap stores the intensity for each boundary faces,
 *        depending of the natur of the boundary (Dirichlet condition).
 *      The intensity of radiation is defined as the rate of emitted
 *        energy from unit surface area through unit solid angle.
 *
 *   1/ Gray wall: isotropic radiation field.
 *                                   4
 *                     eps.sig.twall         (1-eps).qincid
 *       coefap   =    --------------    +    --------------
 *                           pi                     pi
 *   wall intensity     wall emission           reflecting flux.
 *      (eps=1: black wall; eps=0: reflecting wall)
 *   2/ Free boundary: condition to mimic infinite domain
 *
 * 2. Boundary conditions for the P-1 model
 * ----------------------------------------
 *
 * \param[in]  bc_type         boundary face types
 * \param[in]  vect_s          direction vector or NULL
 * \param[out] coefap          boundary conditions for intensity or P-1 model
 * \param[out] coefbp          boundary conditions for intensity or P-1 model
 * \param[out] cofafp          boundary conditions for intensity or P-1 model
 * \param[out] cofbfp          boundary conditions for intensity or P-1 model
 * \param[in]  ckmel           coeff d'absorption du melange
 *                               gaz-particules de charbon
 * \param[in]  w_gg            Weights of the i-th gray gas at boundaries
 * \param[in]  gg_id           number of the i-th grey gas
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_bc_coeffs(int        bc_type[],
                          cs_real_t  vect_s[3],
                          cs_real_t  coefap[],
                          cs_real_t  coefbp[],
                          cs_real_t  cofafp[],
                          cs_real_t  cofbfp[],
                          cs_real_t  ckmel[],
                          cs_real_t  w_gg[],
                          int        gg_id)
{
  cs_real_t stephn = cs_physical_constants_stephan;
  cs_real_t onedpi  = 1.0 / cs_math_pi;

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  cs_real_3_t *b_face_normal
    = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;

  /* Initialization */

  /* Pointer to the spectral flux density field */
  cs_field_t *f_qinspe = NULL;
  if (   cs_glob_rad_transfer_params->imoadf >= 1
      || cs_glob_rad_transfer_params->imfsck == 1)
    f_qinspe = cs_field_by_name_try("spectral_rad_incident_flux");

  /* Pointer to the radiative incident flux field */
  cs_field_t *f_qincid = cs_field_by_name("rad_incident_flux");

  /* Pointer to the wall emissivity field */
  cs_field_t *f_eps = cs_field_by_name("emissivity");

  /* Wall temperature */
  cs_field_t *f_tempb = CS_F_(t_b);
  cs_real_t xptk;
  if (cs_glob_thermal_model->itpscl == 2)
    xptk = cs_physical_constants_celsius_to_kelvin;
  else
    xptk = 0.0;

  /* -> Initialization to a non-admissible value for testing after raycll   */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
    coefap[face_id] = -cs_math_big_r;

  /* Boundary conditions for DO model
   * coefap must be filled with the intensity */

  cs_real_t qpatmp;

  if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_DOM) {

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      /* Copy the appropriate flux density to the local variable qpatmp*/

      /* Value of the flux density at the boundary face     */
      if (   cs_glob_rad_transfer_params->imoadf >= 1
          || cs_glob_rad_transfer_params->imfsck == 1)
        qpatmp = f_qinspe->val[face_id * f_qinspe->dim + gg_id];
      else
        qpatmp = f_qincid->val[face_id];

      /* Dirichlet Boundary Conditions  */

      cs_real_t hint = 0.0;

      if (vect_s != NULL) {
        cs_real_3_t normal;
        cs_math_3_normalise(b_face_normal[face_id], normal);
        if ((bc_type[face_id] == CS_INLET
          || bc_type[face_id] == CS_CONVECTIVE_INLET
          || bc_type[face_id] == CS_OUTLET
          || bc_type[face_id] == CS_FREE_INLET)
            && (CS_ABS(cs_math_3_dot_product(vect_s,  normal)) <cs_math_epzero)) {

        }
      }
      /* Inlet/Outlet face: entering intensity fixed to zero flux
       * which mimic an infinite extrusion
       * (warning: the treatment is different from than of P-1 model)
       * Symmetry: the reflecting boundary conditions (eps=0)
       * is not feasible (because that would required to couple directions)
       * So we impose an homogeneous Neumann that model an infinite extrusion */
      if ( bc_type[face_id] == CS_INLET
        || bc_type[face_id] == CS_CONVECTIVE_INLET
        || bc_type[face_id] == CS_OUTLET
        || bc_type[face_id] == CS_FREE_INLET
        || bc_type[face_id] == CS_SYMMETRY) {


        bool neumann = true;

        /* Entering intensity fixed to zero
         * if the normal to the face is colinear to the direction
         * or if coupled to a 1D module (Atmospheric module)
         * (warning: the treatment is different from than of P-1 model) */
        if (vect_s != NULL) {
          cs_real_3_t normal;
          cs_math_3_normalise(b_face_normal[face_id], normal);
          if (CS_ABS(cs_math_3_dot_product(vect_s,  normal)) <cs_math_epzero)
            neumann = false;
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
          //TODO imposed value at the top for atmospheric
          cs_real_t pimp  = 0.;

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
        cs_real_t twall = f_tempb->val[face_id] + xptk;
        /* Remember: In case of the usage of the standard radiation
           models of Code_Saturne, w_gg=1  */
        cs_real_t pimp = f_eps->val[face_id] * stephn * cs_math_pow4(twall)
                       * onedpi * w_gg[face_id + gg_id * n_b_faces]
                       + (1.0 - f_eps->val[face_id]) * qpatmp * onedpi;

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
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_lnum_t iel = cs_glob_mesh->b_face_cells[face_id];

      cs_real_t hint = 1.0 / (ckmel[iel] * cs_glob_mesh_quantities->b_dist[face_id]);

      /* Symmetry or reflecting wall (EPS = 0) : zero flux */
      if (   bc_type[face_id] == CS_SYMMETRY
          || (   (   bc_type[face_id] == CS_SMOOTHWALL
                  || bc_type[face_id] == CS_ROUGHWALL)
              && f_eps->val[face_id] <= 0.0) ) {
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
        cs_real_t twall = f_tempb->val[face_id] + xptk;
        cs_real_t distbf  = cs_glob_mesh_quantities->b_dist[face_id];
        cs_real_t xit = 1.5 * distbf * ckmel[iel]
                            * (2.0 / (2.0 - f_eps->val[face_id]) - 1.0);
        cs_real_t cfl = 1.0 / xit;
        cs_real_t pimp = cs_math_pow4(twall) * w_gg[face_id + gg_id * n_b_faces];
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

  /* --> Check luminance boundary conditions arrays
   *     In the case of the P-1 approximation, the value of coefap can be
   *     large (of the order of twall**4), hence the value of coefmn */

  cs_real_t xlimit = -cs_math_big_r * 0.1;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (coefap[face_id] <= xlimit)
      bc_type[face_id] = -CS_ABS(bc_type[face_id]);
  }

  cs_boundary_conditions_error(bc_type, "Luminance BC values");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
