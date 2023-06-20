/*============================================================================
 * Base wall condensation model data.
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

#include "cs_log_iteration.h"
#include "cs_math.h"
#include "cs_array.h"
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

static cs_wall_cond_1d_thermal_t _wall_cond_1d_thermal = {.nzones = 0,
                                                          .znmurx = 0,
                                                          .ztheta = NULL,
                                                          .zdxmin = NULL,
                                                          .znmur  = NULL,
                                                          .zepais = NULL,
                                                          .ztpar0 = NULL,

                                                          .zhext  = NULL,
                                                          .ztext  = NULL,
                                                          .zrob   = NULL,
                                                          .zcondb = NULL,
                                                          .zcpb   = NULL,

                                                          .zdxp  = NULL,
                                                          .ztmur = NULL };

static cs_wall_cond_0d_thermal_t
  _wall_cond_0d_thermal = {.nvolumes         = 0,
                           .volume_thickness = NULL,
                           .volume_t         = NULL,
                           .volume_rho       = NULL,
                           .volume_cp        = NULL,
                           .volume_lambda    = NULL,
                           .volume_mass      = NULL,
                           .volume_surf      = NULL,
                           .volume_t0        = NULL,
                           .volume_measure   = NULL};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

// TODO : to remove when the general 1D thermal model replaces
// the condensation-specific 1D thermal model
const cs_wall_cond_1d_thermal_t *cs_glob_wall_cond_1d_thermal
  = &_wall_cond_1d_thermal;

const cs_wall_cond_0d_thermal_t *cs_glob_wall_cond_0d_thermal
  = &_wall_cond_0d_thermal;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void cs_f_wall_condensation_1d_thermal_get_pointers(cs_lnum_t **znmur,
                                                    cs_real_t **ztheta,
                                                    cs_real_t **zdxmin,
                                                    cs_real_t **zepais,
                                                    cs_real_t **zrob,
                                                    cs_real_t **zcondb,
                                                    cs_real_t **zcpb,
                                                    cs_real_t **zhext,
                                                    cs_real_t **ztext,
                                                    cs_real_t **ztpar0);

void cs_f_wall_condensation_1d_thermal_get_mesh_pointers(int **      znmurx,
                                                         cs_real_t **zdxp,
                                                         cs_real_t **ztmur);


void cs_f_wall_condensation_0d_thermal_get_pointers(cs_real_t   **volume_thickness,
                                                    cs_real_2_t **volume_t,
                                                    cs_real_t   **volume_rho,
                                                    cs_real_t   **volume_cp,
                                                    cs_real_t   **volume_lambda,
                                                    cs_real_t   **volume_mass,
                                                    cs_real_t   **volume_surf,
                                                    cs_real_t   **volume_measure,
                                                    cs_real_t   **volume_t0);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_wall_condensation_1d_thermal_get_pointers(cs_lnum_t **znmur,
                                               cs_real_t **ztheta,
                                               cs_real_t **zdxmin,
                                               cs_real_t **zepais,
                                               cs_real_t **zrob,
                                               cs_real_t **zcondb,
                                               cs_real_t **zcpb,
                                               cs_real_t **zhext,
                                               cs_real_t **ztext,
                                               cs_real_t **ztpar0)
{
  *znmur  = _wall_cond_1d_thermal.znmur;
  *ztheta = _wall_cond_1d_thermal.ztheta;
  *zdxmin = _wall_cond_1d_thermal.zdxmin;
  *zepais = _wall_cond_1d_thermal.zepais;
  *zrob   = _wall_cond_1d_thermal.zrob;
  *zcondb = _wall_cond_1d_thermal.zcondb;
  *zcpb   = _wall_cond_1d_thermal.zcpb;
  *zhext  = _wall_cond_1d_thermal.zhext;
  *ztext  = _wall_cond_1d_thermal.ztext;
  *ztpar0 = _wall_cond_1d_thermal.ztpar0;
}

void
cs_f_wall_condensation_1d_thermal_get_mesh_pointers(int **      znmurx,
                                                    cs_real_t **zdxp,
                                                    cs_real_t **ztmur)
{
  *znmurx = &(_wall_cond_1d_thermal.znmurx);
  *zdxp   = _wall_cond_1d_thermal.zdxp;
  *ztmur  = _wall_cond_1d_thermal.ztmur;
}

void cs_f_wall_condensation_0d_thermal_get_pointers(cs_real_t   **volume_thickness,
                                                    cs_real_2_t **volume_t,
                                                    cs_real_t   **volume_rho,
                                                    cs_real_t   **volume_cp,
                                                    cs_real_t   **volume_lambda,
                                                    cs_real_t   **volume_mass,
                                                    cs_real_t   **volume_surf,
                                                    cs_real_t   **volume_measure,
                                                    cs_real_t   **volume_t0)
{
  *volume_thickness    = _wall_cond_0d_thermal.volume_thickness;
  *volume_t            = _wall_cond_0d_thermal.volume_t;
  *volume_rho          = _wall_cond_0d_thermal.volume_rho;
  *volume_cp           = _wall_cond_0d_thermal.volume_cp;
  *volume_lambda       = _wall_cond_0d_thermal.volume_lambda;
  *volume_mass         = _wall_cond_0d_thermal.volume_mass;
  *volume_surf         = _wall_cond_0d_thermal.volume_surf;
  *volume_measure      = _wall_cond_0d_thermal.volume_measure;
  *volume_t0           = _wall_cond_0d_thermal.volume_t0;
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
 * \param[in] nvar     number of variables (?)
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
cs_wall_condensation_0d_thermal_create(cs_lnum_t nvolumes,
                                       cs_lnum_t ncmast)

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

  memset(_wall_cond_0d_thermal.volume_t, 0, ncmast*sizeof(cs_real_2_t));
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_thickness);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_rho);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_cp);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_lambda);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_mass);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_surf);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_t0);
  cs_array_real_fill_zero(nvolumes, _wall_cond_0d_thermal.volume_measure);
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

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_wall_cond_t *wall_cond = cs_glob_wall_cond;

  cs_field_t *f          = cs_field_by_name("pressure");
  const int   var_id_key = cs_field_key_id("variable_id");
  const int   ipr        = cs_field_get_key_int(f, var_id_key) - 1;

  const cs_real_t xlcond  = 2278.0e3;

  cs_real_t *flxmst = wall_cond->flxmst;
  cs_lnum_t *itagms = wall_cond->itagms;
  cs_real_t *mass = _wall_cond_0d_thermal.volume_mass;
  cs_real_t *measure = _wall_cond_0d_thermal.volume_measure;
  cs_real_t *rho = _wall_cond_0d_thermal.volume_rho;
  cs_real_t *thickness = _wall_cond_0d_thermal.volume_thickness;
  cs_real_t *cp = _wall_cond_0d_thermal.volume_cp;
  cs_real_t *lambda = _wall_cond_0d_thermal.volume_lambda;
  cs_real_t *surf = _wall_cond_0d_thermal.volume_surf;
  cs_real_t *t0 = _wall_cond_0d_thermal.volume_t0;
  cs_real_2_t *t = _wall_cond_0d_thermal.volume_t;

  cs_lnum_t temperature_is_var = 0;

  cs_real_t tau_min = cs_math_infinite_r, tau_max = -cs_math_infinite_r;
  cs_real_t tmin[2] = {cs_math_infinite_r, cs_math_infinite_r};
  cs_real_t tmax[2] = {-cs_math_infinite_r, -cs_math_infinite_r};

  for (cs_lnum_t ii = 0; ii < wall_cond->ncmast; ii++) {

    cs_lnum_t c_id = wall_cond->ltmast[ii];
    cs_lnum_t vol_id = wall_cond->izmast[ii];

    temperature_is_var = itagms[vol_id];

    if (temperature_is_var == 1) {

      /* Explicit flux recovered at the fluid interfaces */

      cs_real_t flux = flxmst[ii]
                     - wall_cond->svcond[ipr*wall_cond->ncmast + ii]*xlcond;

      /* Temperature at the volume structures wall and the symmetry
         on the middle of the volume structure cell */

      cs_real_t t_wall = t[ii][0], t_sym = t[ii][1];

      /* Characteristic time of heat propagation on a half
       * thick structure */
      cs_real_t inv_tau = lambda[vol_id]*surf[vol_id]
                        / (thickness[vol_id]/2.0*mass[vol_id]*cp[vol_id]/2.0);

      tau_min = CS_MIN(tau_min, 1.0/inv_tau);
      tau_max = CS_MAX(tau_max, 1.0/inv_tau);

      /* Solve a 0-D unsteady conduction problem at the fluid and
       * the symmetry frontiers with both equations given t_0 and t_1
       * respectively with t_0 the temperature past to the condensation
       * correlations for the volume structures modelling.*/

      /* Compute t_0(n+1) near the mass wall */
      t[ii][0] = t_wall + dt[c_id]*inv_tau
               * (flux*thickness[vol_id]/(2.0*lambda[vol_id]) + t_sym - t_wall);

      /* Compute t_1(n+1) near the symmetry */
      t[ii][1] = t_sym + dt[c_id]*inv_tau*(t_wall-t_sym);

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
               " %15.7e   %15.7e\n", tmin[0], tmax[0]);
    bft_printf(" Min/Max temperature at volume symetry planes:"
               " %15.7e   %15.7e\n", tmin[1], tmax[1]);
    bft_printf(" Min/Max heat propagation characteristic time:"
               " %15.7e   %15.7e\n", tau_min, tau_max);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
