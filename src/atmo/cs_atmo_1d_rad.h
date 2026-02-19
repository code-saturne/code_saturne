#ifndef CS_ATMO_1D_RAD_H
#define CS_ATMO_1D_RAD_H

/*============================================================================
 * Atmospheric radiative fluxes for 1D scheme.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_base.h"

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*! \file cs_atmo_1d_rad.h */

/*----------------------------------------------------------------------------
 * 1-D atmospheric radiative model  option
 *----------------------------------------------------------------------------*/

typedef struct {

  /*! 1-D radiative model (0 off, 1 on) */
  int radiative_model_1d;
  /*! 1-D radiative model: number of vertical arrays */
  int nvert;
  /*! 1-D radiative model: number of levels (up to the top of the domain) */
  int nlevels;
  /*! 1-D radiative model: number of levels (up to 11000 m)
    (automatically computed) */
  int nlevels_max;
  /*! 1D radiative model pass frequency (1 valu bu default)*/
  int frequency;

  /*! horizontal coordinates of the vertical grid */
  cs_real_t *xy;

  /*! vertical grid for 1D radiative scheme */
  cs_real_t *z;
  /*! absorption for CO2 + 03 */
  cs_real_t *acinfe;
  /*! differential absorption for CO2 + 03 */
  cs_real_t *dacinfe;
  /*! absorption for CO2 only */
  cs_real_t *aco2;
  cs_real_t *aco2s;
  /*! differential absorption for CO2 only */
  cs_real_t *daco2;
  cs_real_t *daco2s;
  /*! as acinfe, downwing flux */
  cs_real_t *acsup;
  cs_real_t *acsups;
  cs_real_t *dacsup;
  cs_real_t *dacsups;
  /*! internal variable for 1D radiative model */
  cs_real_t *tauzq;
  /*! internal variable for 1D radiative model */
  cs_real_t *tauz;
  /*! vertical grid for 1D radiative scheme
   * (staggered grid associated to faces) */
  cs_real_t *zq;
  /*! flux divergence of IR radiation */
  cs_real_t *ir_div;
  /*! flux divergence of solar radiation */
  cs_real_t *sol_div;
  /*! Upward and downward radiative fluxes (infrared, solar)
    along each vertical */
  cs_real_t *iru;
  cs_real_t *ird;
  cs_real_t *solu;
  cs_real_t *sold;

  /*! 1D profiles of total water mass fraction along each vertical */
  cs_real_t *qw;
  /*! 1D profiles of liquid water mass fraction along each vertical */
  cs_real_t *ql;
  /*! 1D profiles of vapor water mass fraction along each vertical */
  cs_real_t *qv;
  /*! 1D profiles of number of droplets along each vertical */
  cs_real_t *nc;
  /*! 1D profiles of nebulosity along each vertical */
  cs_real_t *fn;
  /*! 1D profiles of aerosols along each vertical */
  cs_real_t *aerosols;

  /*! Value of ground albedo for each vertical */
  cs_real_t *albedo0;
  /*! Value of ground emissivity for each vertical */
  cs_real_t *emissi0;
  /*! Value of ground temperature for each vertical */
  cs_real_t *temp0;
  /*! Value of ground potential temperature for each vertical */
  cs_real_t *theta0;
  /*! Value of ground total water mass fraction for each vertical */
  cs_real_t *qw0;
  /*! Value of ground pressure for each vertical */
  cs_real_t *p0;
  /*! Value of ground density for each vertical */
  cs_real_t *rho0;

} cs_atmo_1d_rad_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to 1-D atmospheric radiative options structure */
extern cs_atmo_1d_rad_t *cs_glob_atmo_1d_rad;

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief free arrays for atmo 1-D radiative module
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_ATMO_1D_RAD */
