#ifndef __CS_WALLCONDENSATION_1DTHERMAL_H__
#define __CS_WALLCONDENSATION_1DTHERMAL_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/


BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

typedef struct {

  /*! number of the zones with a specific condensation source terms
      depending on the wall temperature and material properties */
  int        nzones;
  /*! Maximal number of discretized points */
  int        znmurx;
  /*! ztheta-scheme of the 1-D thermal model
      - 0 : explicit scheme
      - 1 : implicit scheme */
  cs_real_t *ztheta;
  /*! the minimal space step of 1-D thermal model
      by default equal to 0 with a homogeneus space step.
      this numerical parameter is used to impose a
      geometric progression ratio of the mesh refinement */
  cs_real_t *zdxmin;
  /*! number of discretized points */
  cs_lnum_t *znmur;
  /*! the wall thickness */
  cs_real_t *zepais;
  /*! initial temperature */
  cs_real_t *ztpar0;
  /*! exterior exchange coefficient */
  cs_real_t *zhext;
  /*! exterior temperature */
  cs_real_t *ztext;
  /*! concrete density */
  cs_real_t *zrob;
  /*! concrete conductivity coefficient */
  cs_real_t *zcondb;
  /*!  concrete specific heat coefficient */
  cs_real_t *zcpb;
  /*! initial temperature */
  cs_real_t *ztpar;
  /*! space step */
  cs_real_t *zdxp;
  /*! wall temperature */
  cs_real_t *ztmur;

} cs_wall_cond_1d_thermal_t;

typedef struct {

  /*! number of the volume strutures with a specific condensation source terms */
  cs_lnum_t    nvolumes;
  /*! thickness */
  cs_real_t   *volume_thickness;
  /*! wall temperature */
  cs_real_2_t *volume_t;
  /*! the density (kg.m-3) */
  cs_real_t   *volume_rho;
  /*! specific heat coefficient (J.kg-1.C-1) */
  cs_real_t   *volume_cp;
  /*! conductivity coefficient (W.m-1.C-1) */
  cs_real_t   *volume_lambda;
  /*! metal mass (kg) */
  cs_real_t   *volume_mass;
  /*! exchange surface (m2) */
  cs_real_t   *volume_surf;
  /*! wall temperature */
  cs_real_t   *volume_t0;
  /*! volume mesaure (m3) */
  cs_real_t   *volume_measure;

} cs_wall_cond_0d_thermal_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to wall condensation descriptor structure */

extern const cs_wall_cond_1d_thermal_t *cs_glob_wall_cond_1d_thermal;
extern const cs_wall_cond_0d_thermal_t *cs_glob_wall_cond_0d_thermal;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for wall condensation thermal models.
 *
 * \param[in] nzones number of zones
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_create(int nzones);

void
cs_wall_condensation_1d_thermal_mesh_create(int  znmurx,
                                            int  nfbpcd,
                                            int  nzones);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation models.
 */
/*----------------------------------------------------------------------------*/

void cs_wall_condensation_1d_thermal_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writeable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_1d_thermal_t *
cs_get_glob_wall_cond_1d_thermal(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for 0D wall condensation thermal models.
 *
 * \param[in] nvolumes number of volumes
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_0d_thermal_create(cs_lnum_t  nvolumes,
                                       cs_lnum_t  ncmast);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free structures related to 0D wall condensation models.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_0d_thermal_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writeable access to cs_wall_cond_0d_thermal_t structure.
 *
 * \return pointer to global cs_glob_wall_cond_0d_thermal structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_0d_thermal_t *
cs_get_glob_wall_cond_0d_thermal(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Use 0-D thermal model to solve the temperature and themal flux
 *        at the volume structure walls
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_0d_thermal_solve(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Used to generate the 1-D mesh and initialize
 *        the temperature field of the thermal model
 *        coupled with condensation model.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_1d_thermal_mesh_initialize(void);

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
cs_wall_condensation_1d_thermal_compute_temperature(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALLCONDENSATION_1DTHERMAL_H__ */
