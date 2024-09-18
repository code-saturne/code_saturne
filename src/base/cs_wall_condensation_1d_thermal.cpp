/*============================================================================
 * Base wall condensation model data.
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_defs.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "cs_base.h"
#include "cs_log_iteration.h"
#include "cs_math.h"
#include "cs_wall_condensation.h" // not great...
#include "cs_wall_condensation_1d_thermal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

// TODO : to remove when the general 1D thermal model replaces
// the condensation-specific 1D thermal model

static cs_wall_cond_1d_thermal_t _wall_cond_1d_thermal = { .nzones = 0,
                                                           .znmurx = 0,
                                                           .ztheta = nullptr,
                                                           .zdxmin = nullptr,
                                                           .znmur  = nullptr,
                                                           .zepais = nullptr,
                                                           .ztpar0 = nullptr,

                                                           .zhext  = nullptr,
                                                           .ztext  = nullptr,
                                                           .zrob   = nullptr,
                                                           .zcondb = nullptr,
                                                           .zcpb   = nullptr,

                                                           .zdxp  = nullptr,
                                                           .ztmur = nullptr };

static cs_wall_cond_0d_thermal_t _wall_cond_0d_thermal = {
  .nvolumes         = 0,
  .volume_thickness = nullptr,
  .volume_t         = nullptr,
  .volume_rho       = nullptr,
  .volume_cp        = nullptr,
  .volume_lambda    = nullptr,
  .volume_mass      = nullptr,
  .volume_surf      = nullptr,
  .volume_t0        = nullptr,
  .volume_measure   = nullptr
};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

// TODO : to remove when the general 1D thermal model replaces
// the condensation-specific 1D thermal model
const cs_wall_cond_1d_thermal_t *cs_glob_wall_cond_1d_thermal =
  &_wall_cond_1d_thermal;

const cs_wall_cond_0d_thermal_t *cs_glob_wall_cond_0d_thermal =
  &_wall_cond_0d_thermal;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void
_log_debug(void)
{
  const cs_lnum_t  nfbpcd  = cs_glob_wall_condensation->nfbpcd;
  const cs_lnum_t *izzftcd = cs_glob_wall_condensation->izzftcd;

  cs_real_t *zdxp = _wall_cond_1d_thermal.zdxp;

  const int        nzones = _wall_cond_1d_thermal.nzones;
  const cs_lnum_t *znmur  = _wall_cond_1d_thermal.znmur;
  const cs_real_t *zdxmin = _wall_cond_1d_thermal.zdxmin;
  const cs_real_t *zepais = _wall_cond_1d_thermal.zepais;

  for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
    const cs_lnum_t iz      = izzftcd[ii];
    const cs_real_t _zdxmin = zdxmin[iz];
    const cs_real_t _zepais = zepais[iz];
    const cs_lnum_t _znmur  = znmur[iz] - 1;

    int             iter  = 0;
    cs_real_t       r1    = 2.0;
    cs_real_t       delta = 1.e5;
    const cs_real_t epsi  = 0.0001;

    while (delta > epsi && iter < 100) {
      iter++;
      const cs_real_t r0 = r1;
      r1 = pow(1.0 + (_zepais * (r0 - 1.0)) / _zdxmin, 1.0 / (cs_real_t)_znmur);
      const cs_real_t epai1 = _zdxmin * (pow(r1, _znmur) - 1.0) / (r1 - 1.0);
      delta                 = fabs(epai1 - _zepais) / _zepais;
    }

    bft_printf("------------------------------------------------\n"
               "1-D mesh generation of the thermal model\n"
               "this one is coupled with the condensation model.\n"
               "------------------------------------------------\n"
               "     geometric ratio : %10.07le\n",
               r1);

    cs_real_t r0 = 0.0;
    for (cs_lnum_t kk = 0; kk < _znmur; kk++) {
      r0 += zdxp[iz + kk * nzones];
      if (kk == 0)
        bft_printf("     cell id     cell size      distance to the wall\n");
      bft_printf("           %d          %10.07le        %10.07le\n",
                 kk,
                 zdxp[iz + kk * nzones],
                 r0);
    }
  }
}

/*! (DOXYGEN_SHOULD_SK.IP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for wall condensation models.
 *
 * \param[in] nfbpcd   number of faces with wall condensation
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_create(int nzones)
{
  _wall_cond_1d_thermal.nzones = nzones;

  BFT_MALLOC(_wall_cond_1d_thermal.znmur, nzones, cs_lnum_t);
  BFT_MALLOC(_wall_cond_1d_thermal.ztheta, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.zdxmin, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.zepais, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.zrob, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.zcondb, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.zcpb, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.zhext, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.ztext, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.ztpar0, nzones, cs_real_t);

  for (cs_lnum_t iz = 0; iz < _wall_cond_1d_thermal.nzones; iz++) {
    _wall_cond_1d_thermal.znmur[iz]  = 0;
    _wall_cond_1d_thermal.ztheta[iz] = 0.;
    _wall_cond_1d_thermal.zdxmin[iz] = 0.;
    _wall_cond_1d_thermal.zepais[iz] = 0.;
    _wall_cond_1d_thermal.zrob[iz]   = 0.;
    _wall_cond_1d_thermal.zcondb[iz] = 0.;
    _wall_cond_1d_thermal.zcpb[iz]   = 0.;
    _wall_cond_1d_thermal.zhext[iz]  = 0.;
    _wall_cond_1d_thermal.ztext[iz]  = 0.;
    _wall_cond_1d_thermal.ztpar0[iz] = 0.;
  }

  cs_base_at_finalize(cs_wall_condensation_1d_thermal_free);
}

void
cs_wall_condensation_1d_thermal_mesh_create(int znmurx, int nfbpcd, int nzones)
{
  _wall_cond_1d_thermal.znmurx = znmurx;

  BFT_MALLOC(_wall_cond_1d_thermal.zdxp, nzones * znmurx, cs_real_t);
  BFT_MALLOC(_wall_cond_1d_thermal.ztmur, nfbpcd * znmurx, cs_real_t);

  for (int im = 0; im < znmurx; im++) {
    for (cs_lnum_t ieltcd = 0; ieltcd < nfbpcd; ieltcd++) {
      _wall_cond_1d_thermal.ztmur[ieltcd * znmurx + im] = 0.0;
    }
    for (int iz = 0; iz < nzones; iz++) {
      _wall_cond_1d_thermal.zdxp[iz * znmurx + im] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation 1d thermal models.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_free(void)
{
  BFT_FREE(_wall_cond_1d_thermal.znmur);
  BFT_FREE(_wall_cond_1d_thermal.ztheta);
  BFT_FREE(_wall_cond_1d_thermal.zdxmin);
  BFT_FREE(_wall_cond_1d_thermal.zepais);
  BFT_FREE(_wall_cond_1d_thermal.zrob);
  BFT_FREE(_wall_cond_1d_thermal.zcondb);
  BFT_FREE(_wall_cond_1d_thermal.zcpb);
  BFT_FREE(_wall_cond_1d_thermal.zhext);
  BFT_FREE(_wall_cond_1d_thermal.ztext);
  BFT_FREE(_wall_cond_1d_thermal.ztpar0);
  BFT_FREE(_wall_cond_1d_thermal.zdxp);
  BFT_FREE(_wall_cond_1d_thermal.ztmur);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for 0d wall condensation thermal model.
 *
 * \param[in] nvolumes  number of volumes with metal mass condensation
 * \param[in] ncmast    number of cells with metal mass condensation

 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_0d_thermal_create(cs_lnum_t nvolumes, cs_lnum_t ncmast)
{
  _wall_cond_0d_thermal.nvolumes = nvolumes;
  BFT_MALLOC(_wall_cond_0d_thermal.volume_t, ncmast, cs_real_2_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_thickness, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_rho, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_cp, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_lambda, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_mass, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_surf, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_t0, nvolumes, cs_real_t);
  BFT_MALLOC(_wall_cond_0d_thermal.volume_measure, nvolumes, cs_real_t);

  memset(_wall_cond_0d_thermal.volume_t, 0, ncmast * sizeof(cs_real_2_t));
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_thickness);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_rho);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_cp);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_lambda);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_mass);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_surf);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_t0);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_measure);

  cs_base_at_finalize(cs_wall_condensation_0d_thermal_free);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free structure related to wall condensation 0d thermal model.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_0d_thermal_free(void)
{
  BFT_FREE(_wall_cond_0d_thermal.volume_t);
  BFT_FREE(_wall_cond_0d_thermal.volume_thickness);
  BFT_FREE(_wall_cond_0d_thermal.volume_rho);
  BFT_FREE(_wall_cond_0d_thermal.volume_cp);
  BFT_FREE(_wall_cond_0d_thermal.volume_lambda);
  BFT_FREE(_wall_cond_0d_thermal.volume_mass);
  BFT_FREE(_wall_cond_0d_thermal.volume_surf);
  BFT_FREE(_wall_cond_0d_thermal.volume_t0);
  BFT_FREE(_wall_cond_0d_thermal.volume_measure);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writeable access to _wall_cond_1d_thermal structure.
 *
 * \return pointer to global wall_cond_1d_thermal structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_1d_thermal_t *
cs_get_glob_wall_cond_1d_thermal(void)
{
  return &_wall_cond_1d_thermal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writeable access to _wall_cond_0d_thermal structure.
 *
 * \return pointer to global wall_cond_0d_thermal structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_0d_thermal_t *
cs_get_glob_wall_cond_0d_thermal(void)
{
  return &_wall_cond_0d_thermal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Use 0-D thermal model to solve the temperature and themal flux
 *        at the volume structure walls
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_0d_thermal_solve()
{
  const cs_real_t *restrict dt = CS_F_(dt)->val;

  const cs_wall_condensation_t *wall_cond = cs_glob_wall_condensation;

  cs_field_t *f          = cs_field_by_name("pressure");
  const int   var_id_key = cs_field_key_id("variable_id");
  const int   ipr        = cs_field_get_key_int(f, var_id_key) - 1;

  const cs_real_t xlcond = 2278.0e3;

  cs_real_t   *flxmst    = wall_cond->flxmst;
  cs_lnum_t   *itagms    = wall_cond->itagms;
  cs_real_t   *mass      = _wall_cond_0d_thermal.volume_mass;
  cs_real_t   *thickness = _wall_cond_0d_thermal.volume_thickness;
  cs_real_t   *cp        = _wall_cond_0d_thermal.volume_cp;
  cs_real_t   *lambda    = _wall_cond_0d_thermal.volume_lambda;
  cs_real_t   *surf      = _wall_cond_0d_thermal.volume_surf;
  cs_real_2_t *t         = _wall_cond_0d_thermal.volume_t;

  cs_lnum_t temperature_is_var = 0;

  cs_real_t tau_min = cs_math_infinite_r, tau_max = -cs_math_infinite_r;
  cs_real_t tmin[2] = { cs_math_infinite_r, cs_math_infinite_r };
  cs_real_t tmax[2] = { -cs_math_infinite_r, -cs_math_infinite_r };

  for (cs_lnum_t ii = 0; ii < wall_cond->ncmast; ii++) {
    cs_lnum_t c_id   = wall_cond->ltmast[ii];
    cs_lnum_t vol_id = wall_cond->izmast[ii];

    temperature_is_var = itagms[vol_id];

    if (temperature_is_var == 1) {
      /* Explicit flux recovered at the fluid interfaces */

      cs_real_t flux =
        flxmst[ii] - wall_cond->svcond[ipr * wall_cond->ncmast + ii] * xlcond;

      /* Temperature at the volume structures wall and the symmetry
         on the middle of the volume structure cell */

      cs_real_t t_wall = t[ii][0], t_sym = t[ii][1];

      /* Characteristic time of heat propagation on a half
       * thick structure */
      cs_real_t inv_tau =
        lambda[vol_id] * surf[vol_id] /
        (thickness[vol_id] / 2.0 * mass[vol_id] * cp[vol_id] / 2.0);

      tau_min = CS_MIN(tau_min, 1.0 / inv_tau);
      tau_max = CS_MAX(tau_max, 1.0 / inv_tau);

      /* Solve a 0-D unsteady conduction problem at the fluid and
       * the symmetry frontiers with both equations given t_0 and t_1
       * respectively with t_0 the temperature past to the condensation
       * correlations for the volume structures modelling.*/

      /* Compute t_0(n+1) near the mass wall */
      t[ii][0] = t_wall + dt[c_id] * inv_tau *
                            (flux * thickness[vol_id] / (2.0 * lambda[vol_id]) +
                             t_sym - t_wall);

      /* Compute t_1(n+1) near the symmetry */
      t[ii][1] = t_sym + dt[c_id] * inv_tau * (t_wall - t_sym);

      tmin[0] = CS_MIN(tmin[0], t[ii][0]);
      tmin[1] = CS_MIN(tmin[1], t[ii][1]);

      tmax[0] = CS_MAX(tmax[0], t[ii][0]);
      tmax[1] = CS_MAX(tmax[1], t[ii][1]);
    }
  }

  cs_parall_min(1, CS_REAL_TYPE, &tau_min);
  cs_parall_max(1, CS_REAL_TYPE, &tau_max);
  cs_parall_min(2, CS_REAL_TYPE, tmin);
  cs_parall_max(2, CS_REAL_TYPE, tmax);

  cs_parall_counter_max(&temperature_is_var, 1);

  if (cs_log_default_is_active() && temperature_is_var == 1) {
    bft_printf(" ======================================== \n"
               "  Resolution of the 0-D thermal problem   \n"
               "  coupled with condensation on volume     \n"
               "  structures                              \n"
               " ======================================== \n");

    bft_printf(" Min/Max temperature at fluid interfaces:"
               " %15.7e   %15.7e\n",
               tmin[0],
               tmax[0]);
    bft_printf(" Min/Max temperature at volume symetry planes:"
               " %15.7e   %15.7e\n",
               tmin[1],
               tmax[1]);
    bft_printf(" Min/Max heat propagation characteristic time:"
               " %15.7e   %15.7e\n",
               tau_min,
               tau_max);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Used to generate the 1-D mesh and initialize
 *        the temperature field of the thermal model
 *        coupled with condensation model.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_mesh_initialize(void)
{
  const cs_lnum_t  nfbpcd  = cs_glob_wall_condensation->nfbpcd;
  const cs_lnum_t *izzftcd = cs_glob_wall_condensation->izzftcd;

  cs_real_t *zdxp = _wall_cond_1d_thermal.zdxp;

  const int        nzones = _wall_cond_1d_thermal.nzones;
  const cs_lnum_t *znmur  = _wall_cond_1d_thermal.znmur;
  const cs_real_t *zdxmin = _wall_cond_1d_thermal.zdxmin;
  const cs_real_t *zepais = _wall_cond_1d_thermal.zepais;

  int  iter      = 0;
  bool log_debug = false;

#pragma omp parallel for reduction(+ : iter) if (nfbpcd > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
    const cs_lnum_t iz     = izzftcd[ii];
    const cs_lnum_t _znmur = znmur[iz] - 1;

    /* Generate a homegeneous 1-D mesh with constant space step
         -------------------------------------------------------- */
    if ((zdxmin[iz] <= 0.0) ||
        (zdxmin[iz] > zepais[iz] / (cs_real_t)(_znmur))) {
      for (cs_lnum_t kk = 0; kk < _znmur; kk++)
        zdxp[iz + kk * nzones] = zepais[iz] / (cs_real_t)_znmur;
    }
    /* Generate a heterogeneous 1-D mesh with variable space step
       ---------------------------------------------------------- */
    else {
      /* Compute the geometric ratio with a iterative method */
      iter                  = 0;
      cs_real_t       r1    = 2.0;
      cs_real_t       delta = 1.e5;
      const cs_real_t epsi  = 0.0001;

      while (delta > epsi && iter < 100) {
        iter++;
        const cs_real_t r0 = r1;
        r1                 = pow(1.0 + (zepais[iz] * (r0 - 1.0)) / zdxmin[iz],
                 1.0 / (cs_real_t)_znmur);
        const cs_real_t epai1 =
          zdxmin[iz] * (pow(r1, _znmur) - 1.0) / (r1 - 1.0);
        delta = fabs(epai1 - zepais[iz]) / zepais[iz];
      }

      /* Compute the final 1-D mesh of the thermal model  */
      zdxp[iz] = zdxmin[iz];
      for (cs_lnum_t kk = 1; kk < _znmur; kk++)
        zdxp[iz + kk * nzones] = zdxp[iz + (kk - 1) * nzones] * r1;
    }

  } // end loop on faces with condensation source terms

  if (iter > 99)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error with the 1-D mesh Generation.\n"));

  /* FIXME add verbosity control, or log only total/mean values. */
  if (log_debug)
    _log_debug();

  /* If this is restarted computation, do not reinitialize values */
  if (cs_glob_time_step->nt_prev > 0)
    return;

  /* Initialization of the 1D temperature field for the
     thermal model which is coupled with the condensation model */

  cs_real_t       *ztmur  = _wall_cond_1d_thermal.ztmur;
  const cs_real_t *ztpar0 = _wall_cond_1d_thermal.ztpar0;

#pragma omp parallel for if (nfbpcd > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
    const cs_lnum_t iz = izzftcd[ii];
    for (cs_lnum_t kk = 0; kk < znmur[iz]; kk++)
      ztmur[ii + kk * nfbpcd] = ztpar0[iz];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief The 1D thermal model to compute the temperature to impose
 *        at the cold wall. This one is used by the COPAIN model to estimate
 *        he heat flux at the wall where the condensation occurs.
 *
 *        Is used to compute at each face the
 *        \f$T^{fb}_{\mbox{mur}} \f$ at cold wall.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_compute_temperature(void)
{
  const cs_lnum_t *restrict b_face_cells =
    (const cs_lnum_t *restrict)cs_glob_mesh->b_face_cells;

  const cs_real_t *dt = CS_F_(dt)->val;

  const cs_lnum_t nfbpcd = cs_glob_wall_condensation->nfbpcd;

  const cs_lnum_t *ifbpcd  = cs_glob_wall_condensation->ifbpcd;
  const cs_lnum_t *izzftcd = cs_glob_wall_condensation->izzftcd;
  const cs_lnum_t *iztag1d = cs_glob_wall_condensation->iztag1d;

  const cs_real_t *flthr  = cs_glob_wall_condensation->flthr;
  const cs_real_t *dflthr = cs_glob_wall_condensation->dflthr;

  cs_real_t *ztmur = _wall_cond_1d_thermal.ztmur;

  const int        znmurx = _wall_cond_1d_thermal.znmurx;
  const int        nzones = _wall_cond_1d_thermal.nzones;
  const cs_lnum_t *znmur  = _wall_cond_1d_thermal.znmur;

  const cs_real_t *zcpb   = _wall_cond_1d_thermal.zcpb;
  const cs_real_t *zrob   = _wall_cond_1d_thermal.zrob;
  const cs_real_t *zdxp   = _wall_cond_1d_thermal.zdxp;
  const cs_real_t *zhext  = _wall_cond_1d_thermal.zhext;
  const cs_real_t *ztext  = _wall_cond_1d_thermal.ztext;
  const cs_real_t *ztheta = _wall_cond_1d_thermal.ztheta;
  const cs_real_t *zcondb = _wall_cond_1d_thermal.zcondb;

  /* Resolution of the 1-D thermal problem coupled with condensation
     ---------------------------------------------------------------*/

  cs_real_t dtmur[znmurx];
  cs_real_t da[znmurx], xsm[znmurx], xa[znmurx][2];

  for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
    const cs_lnum_t iz = izzftcd[ii];

    if (iztag1d[iz] != 1)
      continue;

    const cs_lnum_t face_id = ifbpcd[ii];
    const cs_lnum_t c_id    = b_face_cells[face_id];

    const cs_real_t rocp = zrob[iz] * zcpb[iz];
    /* we recover the flow and the derivative of the flow (implicity) */
    const cs_real_t phi  = flthr[ii];
    const cs_real_t dphi = dflthr[ii];

    for (cs_lnum_t kk = 1; kk < znmur[iz] - 1; kk++) {
      const cs_real_t dxv =
        0.5 * (zdxp[iz + (kk - 1) * nzones] + zdxp[iz + kk * nzones]);

      da[kk] = rocp / dt[c_id] +
               ztheta[iz] * zcondb[iz] / (zdxp[iz + (kk - 1) * nzones] * dxv) +
               ztheta[iz] * zcondb[iz] / (zdxp[iz + kk * nzones] * dxv);
      xa[kk][0] =
        -ztheta[iz] * zcondb[iz] / (zdxp[iz + (kk - 1) * nzones] * dxv);
      xa[kk][1] = -ztheta[iz] * zcondb[iz] / (zdxp[iz + kk * nzones] * dxv);
      xsm[kk] =
        zcondb[iz] *
        (ztmur[ii + (kk + 1) * nfbpcd] / (zdxp[iz + kk * nzones] * dxv) -
         ztmur[ii + kk * nfbpcd] / (zdxp[iz + kk * nzones] * dxv) -
         ztmur[ii + kk * nfbpcd] / (zdxp[iz + (kk - 1) * nzones] * dxv) +
         ztmur[ii + (kk - 1) * nfbpcd] / (zdxp[iz + (kk - 1) * nzones] * dxv));
    }

    /* fluide side */
    cs_real_t dx  = zdxp[iz];
    cs_real_t dx2 = zdxp[iz] * zdxp[iz];
    da[0] =
      rocp / dt[c_id] + ztheta[iz] * 2.0 * zcondb[iz] / dx2 + 2.0 * dphi / dx;
    xa[0][0] = 0.0;
    xa[0][1] = -ztheta[iz] * 2.0 * zcondb[iz] / dx2;
    xsm[0]   = 2.0 * zcondb[iz] / dx2 * (ztmur[ii + nfbpcd] - ztmur[ii]) +
             (2.0 / dx) * phi;

    /* extern side */
    const cs_lnum_t k = znmur[iz] - 1;
    dx                = zdxp[iz + (k - 1) * nzones];
    dx2   = zdxp[iz + (k - 1) * nzones] * zdxp[iz + (k - 1) * nzones];
    da[k] = rocp / dt[c_id] + ztheta[iz] * 2.0 * zcondb[iz] / dx2 +
            2.0 * zhext[iz] / dx;
    xa[k][0] = -ztheta[iz] * 2.0 * zcondb[iz] / dx2;
    xa[k][1] = 0.0;
    xsm[k]   = 2.0 * zcondb[iz] / dx2 *
               (ztmur[ii + (k - 1) * nfbpcd] - ztmur[ii + k * nfbpcd]) -
             (2.0 / dx) * zhext[iz] * (ztmur[ii + k * nfbpcd] - ztext[iz]);

    /* Resolution on increment */
    for (cs_lnum_t kk = 0; kk < znmur[iz]; kk++)
      dtmur[kk] = 0.0;

    dtmur[0] = (xsm[0] + xa[0][1] * dtmur[1]) / da[0];
    for (cs_lnum_t kk = 1; kk < znmur[iz] - 1; kk++)
      // Not a real theta scheme:
      // dtmur[kk+1] corresponds to the previous iteration
      dtmur[kk] =
        (xsm[kk] + xa[kk][0] * dtmur[kk - 1] + xa[kk][1] * dtmur[kk + 1]) /
        da[kk];
    dtmur[k] = (xsm[k] + xa[k][0] * dtmur[k - 1]) / da[k];

    // Temperature update
    for (cs_lnum_t kk = 0; kk < znmur[iz]; kk++)
      ztmur[ii + kk * nfbpcd] += dtmur[kk];

  } // end loop on faces with condensation source terms

  if (!(cs_log_default_is_active()))
    return;

  cs_real_t tpminf[nzones], tpmaxf[nzones], tpminp[nzones], tpmaxp[nzones];

  for (int z_id = 0; z_id < nzones; z_id++) {
    tpminf[z_id] = +1.e20;
    tpmaxf[z_id] = -1.e20;
    tpminp[z_id] = +1.e20;
    tpmaxp[z_id] = -1.e20;
  }

  for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
    const cs_lnum_t iz = izzftcd[ii];

    if (iztag1d[iz] != 1)
      continue;

    const cs_lnum_t kk = znmur[iz] - 1;
    tpminf[iz]         = cs_math_fmin(tpminf[iz], ztmur[ii]);
    tpmaxf[iz]         = cs_math_fmax(tpmaxf[iz], ztmur[ii]);
    tpminp[iz]         = cs_math_fmin(tpminp[iz], ztmur[ii + kk * nfbpcd]);
    tpmaxp[iz]         = cs_math_fmax(tpmaxp[iz], ztmur[ii + kk * nfbpcd]);
  }

  cs_parall_min(nzones, CS_REAL_TYPE, tpminf);
  cs_parall_min(nzones, CS_REAL_TYPE, tpminp);
  cs_parall_max(nzones, CS_REAL_TYPE, tpmaxf);
  cs_parall_max(nzones, CS_REAL_TYPE, tpmaxp);

  bft_printf("=====================================\n"
             "Resolution of the 1-D thermal problem\n"
             " coupled with the condensation model\n"
             "=====================================\n"
             "------------------------------------------"
             "----------------------------------------------\n"
             "time         izones     Tp_f  (min)      Tp_f"
             "   (max)      Tp_ext(min)      Tp_ext (max)\n"
             "(s)                      (C)             (C)"
             "               (C)             (C)        \n"
             "------------------------------------------"
             "----------------------------------------------\n");

  for (int z_id = 0; z_id < nzones; z_id++)
    bft_printf("  %10.06f        %d      %10.06f        %10.06f        %10.06f"
               "        %10.06f\n",
               cs_glob_time_step->t_cur,
               z_id,
               tpminf[z_id],
               tpmaxf[z_id],
               tpminp[z_id],
               tpmaxp[z_id]);
  bft_printf("------------------------------------------"
             "----------------------------------------------\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
