/*============================================================================
 * Base wall condensation model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_defs.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_math.h"
#include "cs_log_iteration.h"
#include "cs_wall_condensation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*!
 * Karman constant. (= 0.42)
 *
 * Useful if and only if \ref iturb >= 10.
 *  (mixing length, \f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$,
 * LES, v2f or \f$k-\omega\f$).
 */

//const double cs_turb_xkappa = 0.42;

static cs_wall_cond_t _wall_cond =
{
  .icondb     = -1,
  .model      = CS_WALL_COND_MODEL_COPAIN,
  .regime     = CS_WALL_COND_REGIME_NATURAL_CONVECTION,

   // Mesh related quantities
   // TODO : clean unnecessary quantities
  .nfbpcd     = 0,
  .ifbpcd     = NULL,
  .itypcd     = NULL,
  .izzftcd    = NULL,
  .spcond     = NULL,
  .hpcond     = NULL,
  .twall_cond = NULL,
  .thermal_condensation_flux = NULL,
  .flthr      = NULL,
  .dflthr     = NULL,

  // Zone related quantities
  .nzones     = -1,
  .izcophc    = NULL,
  .izcophg    = NULL,
  .iztag1d    = NULL,
  .ztpar      = NULL,
  .zxrefcond  = NULL,
  .zprojcond  = NULL
};

const cs_wall_cond_t *cs_glob_wall_cond = &_wall_cond;

/*============================================================================
 * Fortran function prototypes
 *============================================================================*/

void condensation_copain_model
(
  const int nvar,
  const int nfbpcd,
  const int ifbpcd[],
  const int izzftcd[],
  cs_real_t spcond[],
  cs_real_t hpcond[],
  const int icondb_regime
);

void condensation_copain_benteboula_dabbene_model
(
  const int nvar,
  const int nfbpcd,
  const int ifbpcd[],
  const int izzftcd[],
  cs_real_t spcond[],
  cs_real_t hpcond[],
  const int icondb_regime
);

void condensation_uchida_model
(
  const int nvar,
  const int nfbpcd,
  const int ifbpcd[],
  const int izzftcd[],
  cs_real_t spcond[],
  cs_real_t hpcond[],
  const int icondb_regime
);

void condensation_dehbi_model
(
  const int nvar,
  const int nfbpcd,
  const int ifbpcd[],
  const int izzftcd[],
  cs_real_t spcond[],
  cs_real_t hpcond[],
  const int icondb_regime
);

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_wall_condensation_get_model_pointers(int **icondb,
                                          cs_wall_cond_model_t **icondb_model,
                                          cs_wall_cond_regime_t **icondb_regime);

void
cs_f_wall_condensation_get_size_pointers(cs_lnum_t **nfbpcd, cs_lnum_t **nzones);

void
cs_f_wall_condensation_get_pointers(cs_lnum_t **ifbpcd, cs_lnum_t **itypcd,
                                    cs_lnum_t **izzftcd, cs_real_t **spcond,
                                    cs_real_t **hpcond, cs_real_t **twall_cond,
                                    cs_real_t **thermflux, cs_real_t **flthr,
                                    cs_real_t **dflthr, cs_lnum_t **izcophc,
                                    cs_lnum_t **izcophg, cs_lnum_t **iztag1d,
                                    cs_real_t **ztpar, cs_real_t **xrefcond,
                                    cs_real_t **projcond);


void
cs_wall_condensation_set_model(cs_wall_cond_model_t model);

void
cs_wall_condensation_set_regime(cs_wall_cond_regime_t regime);

void
cs_wall_condensation_set_onoff_state(int icondb);

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void compute_mix_properties(void);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_wall_condensation_get_model_pointers(int **icondb,
                                          cs_wall_cond_model_t **icondb_model,
                                          cs_wall_cond_regime_t **icondb_regime)
{
  *icondb        = &(_wall_cond.icondb);
  *icondb_model  = &(_wall_cond.model);
  *icondb_regime = &(_wall_cond.regime);
}

void
cs_f_wall_condensation_get_size_pointers(cs_lnum_t **nfbpcd, cs_lnum_t **nzones)
{
  *nfbpcd = &(_wall_cond.nfbpcd);
  *nzones = &(_wall_cond.nzones);
}

void
cs_f_wall_condensation_get_pointers(cs_lnum_t **ifbpcd, cs_lnum_t **itypcd,
                                    cs_lnum_t **izzftcd, cs_real_t **spcond,
                                    cs_real_t **hpcond, cs_real_t **twall_cond,
                                    cs_real_t **thermflux, cs_real_t **flthr,
                                    cs_real_t **dflthr, cs_lnum_t **izcophc,
                                    cs_lnum_t **izcophg, cs_lnum_t **iztag1d,
                                    cs_real_t **ztpar, cs_real_t **zxrefcond,
                                    cs_real_t **zprojcond)
{
  *ifbpcd = _wall_cond.ifbpcd ;
  *itypcd = _wall_cond.itypcd;
  *izzftcd = _wall_cond.izzftcd;
  *spcond = _wall_cond.spcond ;
  *hpcond = _wall_cond.hpcond ;
  *twall_cond = _wall_cond.twall_cond;
  *thermflux = _wall_cond.thermal_condensation_flux;
  *flthr = _wall_cond.flthr;
  *dflthr = _wall_cond.dflthr;
  *izcophc = _wall_cond.izcophc;
  *izcophg = _wall_cond.izcophg;
  *iztag1d = _wall_cond.iztag1d;
  *ztpar = _wall_cond.ztpar;
  *zxrefcond = _wall_cond.zxrefcond;
  *zprojcond = _wall_cond.zprojcond;
}

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the wall condensation model
 *
 * \param[in] model    integer corresponding to the desired model
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_model(cs_wall_cond_model_t  model)
{
  _wall_cond.model = model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the wall condensation regime
 *
 * \param[in] model    integer corresponding to the desired model
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_regime(cs_wall_cond_regime_t  regime)
{
  _wall_cond.regime = regime;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the onoff state of wall condensation modeling
 *
 * \param[in] icondb integer corresponding to the onoff state (-1 : off, 0: on)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_onoff_state(int  icondb)
{
  _wall_cond.icondb = icondb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for wall condensation models.
 *
 * \param[in] nfbpcd   number of faces with wall condensation
 * \param[in] nzones   number of zones with wall condensation
 * \param[in] nvar     number of variables (?)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_create(cs_lnum_t  nfbpcd,
                            cs_lnum_t  nzones,
                            cs_lnum_t  nvar)
{
  _wall_cond.nfbpcd = nfbpcd;
  if (nzones < 1) {
    _wall_cond.nzones = 1;
  }
  else {
    _wall_cond.nzones = nzones;
  }

  // Mesh related quantities
  BFT_MALLOC(_wall_cond.ifbpcd, nfbpcd, cs_lnum_t);
  BFT_MALLOC(_wall_cond.itypcd, nfbpcd*nvar, cs_lnum_t);
  BFT_MALLOC(_wall_cond.izzftcd, nfbpcd, cs_lnum_t);
  BFT_MALLOC(_wall_cond.spcond, nfbpcd*nvar, cs_real_t);
  BFT_MALLOC(_wall_cond.hpcond, nfbpcd, cs_real_t);
  BFT_MALLOC(_wall_cond.twall_cond, nfbpcd, cs_real_t);
  BFT_MALLOC(_wall_cond.thermal_condensation_flux, nfbpcd, cs_real_t);
  BFT_MALLOC(_wall_cond.flthr, nfbpcd, cs_real_t);
  BFT_MALLOC(_wall_cond.dflthr, nfbpcd, cs_real_t);

  // Zone related quantities
  BFT_MALLOC(_wall_cond.izcophc, nzones, cs_lnum_t);
  BFT_MALLOC(_wall_cond.izcophg, nzones, cs_lnum_t);
  BFT_MALLOC(_wall_cond.iztag1d, nzones, cs_lnum_t);
  BFT_MALLOC(_wall_cond.ztpar, nzones, cs_real_t);
  BFT_MALLOC(_wall_cond.zxrefcond, 3*nzones, cs_real_t);
  BFT_MALLOC(_wall_cond.zprojcond, 3*nzones, cs_real_t);

  for (cs_lnum_t i= 0 ; i< nfbpcd; i++) {
    _wall_cond.ifbpcd[i] = 0;
    _wall_cond.hpcond[i] = 0.0;
    _wall_cond.twall_cond[i] = 0.0;
    _wall_cond.thermal_condensation_flux[i] = 0.0;
    _wall_cond.flthr[i] = 0.0;
    _wall_cond.dflthr[i] = 0.0;
    if (_wall_cond.nzones <= 1) {
      _wall_cond.izzftcd[i] = 1;
    }
    else {
      _wall_cond.izzftcd[i] = 0;
    }
    for (cs_lnum_t j=0; j<nvar; j++) {
      _wall_cond.itypcd[nvar*i+j] = 0;
      _wall_cond.spcond[nvar*i+j] = 0.0;
    }
  }

  for (cs_lnum_t i=0; i<nzones; i++) {
    _wall_cond.izcophc[i] = 0;
    _wall_cond.izcophg[i] = 0;
    _wall_cond.iztag1d[i] = 0;
    _wall_cond.ztpar[i] = -1.0;
    _wall_cond.zxrefcond[3*i] = 0.0;
    _wall_cond.zxrefcond[3*i+1] = 0.0;
    _wall_cond.zxrefcond[3*i+2] = 0.0;
    _wall_cond.zprojcond[3*i] = 0.0;
    _wall_cond.zprojcond[3*i+1] = 0.0;
    _wall_cond.zprojcond[3*i+2] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation models
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_free(void)
{
  BFT_FREE(_wall_cond.ifbpcd);
  BFT_FREE(_wall_cond.itypcd);
  BFT_FREE(_wall_cond.izzftcd);
  BFT_FREE(_wall_cond.spcond);
  BFT_FREE(_wall_cond.hpcond);
  BFT_FREE(_wall_cond.twall_cond);
  BFT_FREE(_wall_cond.thermal_condensation_flux);
  BFT_FREE(_wall_cond.flthr);
  BFT_FREE(_wall_cond.dflthr);

  BFT_FREE(_wall_cond.izcophc);
  BFT_FREE(_wall_cond.izcophg);
  BFT_FREE(_wall_cond.iztag1d);
  BFT_FREE(_wall_cond.ztpar);
  BFT_FREE(_wall_cond.zxrefcond);
  BFT_FREE(_wall_cond.zprojcond);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the wall condensation source terms.
 *
 * \param[in] nvar     number of variables (?)
 * \param[in] izzftcd  pointer to the table connecting faces to their
 *                     condensation zone
 *
 * \return
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_compute(int        nvar,
                             cs_lnum_t  nfbpcd,
                             cs_lnum_t  ifbpcd[],
                             int        izzftcd[],
                             cs_real_t  spcond[],
                             cs_real_t  hpcond[])
{
  int icondb_regime = _wall_cond.regime;

  switch (_wall_cond.model) {
    case CS_WALL_COND_MODEL_NONE:
      return;
    case CS_WALL_COND_MODEL_COPAIN:
      condensation_copain_model(nvar, nfbpcd, ifbpcd, izzftcd, spcond,
		      hpcond, icondb_regime);
      return;
    case CS_WALL_COND_MODEL_COPAIN_BD:
      condensation_copain_benteboula_dabbene_model(nvar, nfbpcd, ifbpcd,
		      izzftcd, spcond, hpcond, icondb_regime);
      return;
    case CS_WALL_COND_MODEL_UCHIDA:
      condensation_uchida_model(nvar, nfbpcd, ifbpcd, izzftcd, spcond,
		      hpcond, icondb_regime);
      return;
    case CS_WALL_COND_MODEL_DEHBI:
      condensation_dehbi_model(nvar, nfbpcd, ifbpcd, izzftcd, spcond,
		      hpcond, icondb_regime);
      return;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_t *
cs_get_glob_wall_cond(void)
{
  return &_wall_cond;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
