#ifndef __CS_WALLCONDENSATION_MODEL_H__
#define __CS_WALLCONDENSATION_MODEL_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {
  CS_WALL_COND_MODEL_COPAIN    = 0,
  CS_WALL_COND_MODEL_COPAIN_BD = 1,
  CS_WALL_COND_MODEL_UCHIDA    = 2,
  CS_WALL_COND_MODEL_DEHBI     = 3,
} cs_wall_cond_natural_conv_model_t;

typedef enum {
  CS_WALL_COND_MODEL_WALL_LAW    = 0,
  CS_WALL_COND_MODEL_SCHLICHTING = 1
} cs_wall_cond_forced_conv_model_t;

typedef enum {
  CS_WALL_COND_MIXED_MAX       = 0,
  CS_WALL_COND_MIXED_INCROPERA = 1
} cs_wall_cond_mixed_conv_model_t;

typedef struct {
  int icondb; // Switch used to activate wall condensation (0 : activated)
  // Model type
  cs_wall_cond_natural_conv_model_t natural_conv_model;
  cs_wall_cond_forced_conv_model_t  forced_conv_model;
  cs_wall_cond_mixed_conv_model_t   mixed_conv_model;

  // Mesh-related information
  cs_lnum_t  nfbpcd;
  cs_lnum_t *ifbpcd;
  cs_lnum_t *itypcd;
  cs_lnum_t *izzftcd;
  cs_real_t *spcond;
  cs_real_t *hpcond;
  cs_real_t *twall_cond;
  cs_real_t *thermal_condensation_flux;
  cs_real_t *convective_htc;
  cs_real_t *condensation_htc;
  cs_real_t *total_htc;
  cs_real_t *flthr;
  cs_real_t *dflthr;

  // Zone related quantities
  cs_lnum_t  nzones;
  cs_lnum_t *izcophc;
  cs_lnum_t *izcophg;
  cs_lnum_t *iztag1d;
  cs_real_t *ztpar;
  cs_real_t *zxrefcond;
  cs_real_t *zprojcond;
} cs_wall_cond_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to wall condensation descriptor structure */
extern const cs_wall_cond_t *cs_glob_wall_cond;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the wall condensation model
 *
 * \param[in] model    integer corresponding to the desired model
 */
/*----------------------------------------------------------------------------*/
void cs_wall_condensation_set_model(cs_wall_cond_natural_conv_model_t model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the onoff state of wall condensation modeling
 *
 * \param[in] icondb integer corresponding to the onoff state (-1 : off, 0: on)
 */
/*----------------------------------------------------------------------------*/
void cs_wall_condensation_set_onoff_state(int icondb);

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
cs_wall_condensation_create(cs_lnum_t nfbpcd, cs_lnum_t nzones, cs_lnum_t nvar);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation models
 */
/*----------------------------------------------------------------------------*/
void cs_wall_condensation_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the wall condensation source terms.
 *
 * \return
 */
/*----------------------------------------------------------------------------*/
void cs_wall_condensation_compute(cs_real_t total_htc[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output statistics about wall condensation source terms (for user log)
 *
 * \return
 */
/*----------------------------------------------------------------------------*/
void cs_wall_condensation_log(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_t *cs_get_glob_wall_cond(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALLCONDESATION_MODEL_H__ */
