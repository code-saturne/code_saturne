/*============================================================================
 * Base wall condensation model data.
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
#include "cs_wall_condensation.h" // not great...
#include "cs_wall_condensation_1d_thermal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*!
 * Karman constant. (= 0.42)
 *
 * Useful if and only if \ref iturb >= 10.
 *  (mixing length, \f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$,
 * LES, v2f or \f$k-\omega\f$).
 */
// const double cs_turb_xkappa = 0.42;

//
// TODO : to remove when the general 1D thermal model replaces
// the condensation-specific 1D thermal model

static cs_wall_cond_1d_thermal_t _wall_cond_thermal = { .nzones = 0,
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

// TODO : to remove when the general 1D thermal model replaces
// the condensation-specific 1D thermal model
const cs_wall_cond_1d_thermal_t *cs_glob_wall_cond_1d_thermal
  = &_wall_cond_thermal;

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
  *znmur  = _wall_cond_thermal.znmur;
  *ztheta = _wall_cond_thermal.ztheta;
  *zdxmin = _wall_cond_thermal.zdxmin;
  *zepais = _wall_cond_thermal.zepais;
  *zrob   = _wall_cond_thermal.zrob;
  *zcondb = _wall_cond_thermal.zcondb;
  *zcpb   = _wall_cond_thermal.zcpb;
  *zhext  = _wall_cond_thermal.zhext;
  *ztext  = _wall_cond_thermal.ztext;
  *ztpar0 = _wall_cond_thermal.ztpar0;
}

void
cs_f_wall_condensation_1d_thermal_get_mesh_pointers(int **      znmurx,
                                                    cs_real_t **zdxp,
                                                    cs_real_t **ztmur)
{
  *znmurx = &(_wall_cond_thermal.znmurx);
  *zdxp   = _wall_cond_thermal.zdxp;
  *ztmur  = _wall_cond_thermal.ztmur;
}

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
  _wall_cond_thermal.nzones = nzones;

  BFT_MALLOC(_wall_cond_thermal.znmur, nzones, cs_lnum_t);
  BFT_MALLOC(_wall_cond_thermal.ztheta, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.zdxmin, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.zepais, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.zrob, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.zcondb, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.zcpb, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.zhext, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.ztext, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.ztpar0, nzones, cs_real_t);

  for (cs_lnum_t iz = 0; iz < _wall_cond_thermal.nzones; iz++) {
    _wall_cond_thermal.znmur[iz]  = 0;
    _wall_cond_thermal.ztheta[iz] = 0.;
    _wall_cond_thermal.zdxmin[iz] = 0.;
    _wall_cond_thermal.zepais[iz] = 0.;
    _wall_cond_thermal.zrob[iz]   = 0.;
    _wall_cond_thermal.zcondb[iz] = 0.;
    _wall_cond_thermal.zcpb[iz]   = 0.;
    _wall_cond_thermal.zhext[iz]  = 0.;
    _wall_cond_thermal.ztext[iz]  = 0.;
    _wall_cond_thermal.ztpar0[iz] = 0.;
  }
}

void
cs_wall_condensation_1d_thermal_mesh_create(int znmurx, int nfbpcd, int nzones)
{
  _wall_cond_thermal.znmurx = znmurx;

  BFT_MALLOC(_wall_cond_thermal.zdxp, nzones * znmurx, cs_real_t);
  BFT_MALLOC(_wall_cond_thermal.ztmur, nfbpcd * znmurx, cs_real_t);

  for (int im = 0; im < znmurx; im++) {
    for (cs_lnum_t ieltcd = 0; ieltcd < nfbpcd; ieltcd++) {
      _wall_cond_thermal.ztmur[ieltcd * znmurx + im] = 0.0;
    }
    for (int iz = 0; iz < nzones; iz++) {
      _wall_cond_thermal.zdxp[iz * znmurx + im] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation models.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_free(void)
{
  BFT_FREE(_wall_cond_thermal.znmur);
  BFT_FREE(_wall_cond_thermal.ztheta);
  BFT_FREE(_wall_cond_thermal.zdxmin);
  BFT_FREE(_wall_cond_thermal.zepais);
  BFT_FREE(_wall_cond_thermal.zrob);
  BFT_FREE(_wall_cond_thermal.zcondb);
  BFT_FREE(_wall_cond_thermal.zcpb);
  BFT_FREE(_wall_cond_thermal.zhext);
  BFT_FREE(_wall_cond_thermal.ztext);
  BFT_FREE(_wall_cond_thermal.ztpar0);
  BFT_FREE(_wall_cond_thermal.zdxp);
  BFT_FREE(_wall_cond_thermal.ztmur);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writeable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_1d_thermal_t *
cs_get_glob_wall_cond_1d_thermal(void)
{
  return &_wall_cond_thermal;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
