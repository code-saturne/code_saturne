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
  /*! aerosols active (0 off, 1 on) */
  int has_aerosol;
  /*! output verbosity (0 off, 1 on, etc.) */
  int verbosity;

  /*! internal variable for 1D radiative model */
  cs_real_t tausup;

  /*! adimensional :  aod_o3_tot=0.2 other referenced values are  0.10, 0.16 */
  cs_real_t aod_o3_tot;
  /*! adimensional :  aod_h2o_tot=0.10 other referenced values are  0.06, 0.08 */
  cs_real_t aod_h2o_tot;
  /*! total aerosol optical depth in the IR domain for thermal radiation
      deduced from aeronet data */
  cs_real_t aod_ir;
  /*!  CO2 concentration in cm NTP with correction to the ratio of
    molar masses for Mco2 and Mair
    default is 350ppm */
  cs_real_t conco2;
  /*! Asymmetry factor for O3 (non-dimensional)
      climatic value gaero_o3=0.66 */
  cs_real_t gaero_o3;
  /*! Asymmetry factor for H2O (non-dimensional)
      climatic value gaero_h2o=0.64 */
  cs_real_t gaero_h2o;
  /*! Single scattering albedo for O3 (non-dimensional)
     climatic value piaero_o3=0.84, other referenced values are 0.963 */
  cs_real_t piaero_o3;
  /*! Single scattering albedo for H2O (non-dimensional)
     climatic value piaero_h2o=0.84, other referenced values are 0.964 */
  cs_real_t piaero_h2o;
  /*! Fraction of Black carbon (non-dimensional): black_carbon_frac=1.d-8 for no BC */
  cs_real_t black_carbon_frac;
  /*! Maximal height for aerosol distribution on the vertical
    important should be <= zqq(kmray-1);
    in meters : referenced value: zaero=6000 */
  cs_real_t zaero;

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
 * \brief Allocate arrays for atmo 1-D radiative module
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_atmo_1d_rad_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief free arrays for atmo 1-D radiative module
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_atmo_1d_rad_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute radiative source term for the atmospheric model (1D scheme).
 *
 *         Computes the source term for scalar equations from radiative forcing
 *         (UV and IR radiative fluxes) using a 1D formulation.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_compute_radiative_source_term_1d(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute solar fluxes for both clear and cloudy atmosphere following
 * Lacis and Hansen (1974). The multiple diffusion is taken into account by an
 * addition method and overlapping between water vapor and liquid water with k
 * distribution method.
 * Some improvements from original version concerns:
 * - introduction of cloud fraction with hazardous recovering
 * - introduction of aerosol diffusion in the same way as for cloud droplets
 *   but with specific optical properties for aerosols.
 *
 * \param[in]   ivertc      index of vertical profile
 * \param[in]   k1          index for ground level
 * \param[in]   kmray       vertical levels number for radiation
 * \param[in]   heuray      Universal time (Hour)
 * \param[in]   imer1       sea index
 * \param[in]   albe        albedo
 * \param[in]   qqv         optical depth for water vapor (0,zqq)
 * \param[in]   qqvinf      idem qqv but for altitude above 11000m
 * \param[in]   zqq         vertical levels (interfaces)
 * \param[in]   zray        vertical levels (volumes)
 * \param[in]   qvray       specific humidity for water vapor
 * \param[in]   qlray       specific humidity for liquid water
 * \param[in]   fneray      cloud fraction
 * \param[in]   romray      air density
 * \param[in]   preray      pressure
 * \param[in]   temray      temperature
 * \param[out]  fos         global downward solar flux at the ground
 * \param[out]  rayst       flux divergence of solar radiation
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_atmo_1d_rad_compute_solar(const int       ivertc,
                             const int       k1,
                             const int       kmray,
                             const cs_real_t heuray,
                             const int       imer1,
                             cs_real_t       *albe,
                             cs_real_t       qqv[],
                             const cs_real_t qqvinf,
                             const cs_real_t zqq[],
                             const cs_real_t zray[],
                             const cs_real_t qvray[],
                             cs_real_t       qlray[],
                             cs_real_t       fneray[],
                             const cs_real_t romray[],
                             const cs_real_t preray[],
                             const cs_real_t temray[],
                             cs_real_t       *fos,
                             cs_real_t       rayst[],
                             const cs_real_t ncray[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute infrared flux divergence profile and downward flux at
 *        ground level relying on a 1D radiative scheme.
 *
 * \param[in]  ivertc      index of vertical profile
 * \param[in]  k1          index for ground level
 * \param[in]  kmray       vertical levels number for radiation
 * \param[in]  emis        ground surface emissivity
 * \param[out] qqv         optical depth for water vapor (0,zqq)
 * \param[out] qqqv        idem qqv but for intermediates vertical levels (zray)
 * \param[in]  qqvinf      idem qqv but for altitude above 11000m
 * \param[in]  zqq         vertical levels (interfaces)
 * \param[in]  zray        vertical levels (volumes)
 * \param[in]  temray      temperature in Celsius
 * \param[in]  qvray       specific humidity for water vapor
 * \param[in]  qlray       specific humidity for liquid water
 * \param[in]  fnerir      cloud fraction
 * \param[in]  romray      air density
 * \param[in]  preray      pressure
 * \param[in]  aeroso      aerosol concentration in micro-g/m3
 * \param[in]  t_surf      surface temperature
 * \param[in]  p_surf      surface pressure
 * \param[out] foir        downward IR flux at the ground
 * \param[out] rayi        IR flux divergence
 * \param[in]  ncray       Number of droplets interpolated on vertical grid
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_atmo_1d_rad_compute_infrared(const int        ivertc,
                                const int        k1,
                                const int        kmray,
                                const cs_real_t  emis,
                                cs_real_t        qqv[],
                                cs_real_t        qqqv[],
                                cs_real_t        *qqvinf,
                                cs_real_t        zqq[],
                                cs_real_t        zray[],
                                const cs_real_t  temray[],
                                const cs_real_t  qvray[],
                                const cs_real_t  qlray[],
                                cs_real_t        fnerir[],
                                const cs_real_t  romray[],
                                const cs_real_t  preray[],
                                const cs_real_t  aeroso[],
                                const cs_real_t  t_surf,
                                const cs_real_t  p_surf,
                                cs_real_t        *foir,
                                cs_real_t        rayi[],
                                const cs_real_t  ncray[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute radiative fluxes for the atmospheric model.
 *
 * Computes the source term for scalar equations from radiative forcing
 * (UV and IR radiative fluxes) with a 1D scheme.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_atmo_1d_rad_source_term(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_ATMO_1D_RAD */
