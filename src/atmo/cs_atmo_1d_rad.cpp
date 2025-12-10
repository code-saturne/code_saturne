/*============================================================================
 * Atmospheric radiative fluxes for 1D scheme
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_array.h"
#include "base/cs_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_1d_rad.h"

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* local atmo 1d radiative structure */
static cs_atmo_1d_rad_t _atmo_1d_rad = {
  .radiative_model_1d = 0,
  .nvert = 1,
  .nlevels = 20,
  .nlevels_max = 0,
  .frequency = 1,
  .xy = nullptr,
  .z = nullptr,
  .acinfe = nullptr,
  .dacinfe = nullptr,
  .aco2 = nullptr,
  .aco2s = nullptr,
  .daco2 = nullptr,
  .daco2s = nullptr,
  .acsup = nullptr,
  .acsups = nullptr,
  .dacsup = nullptr,
  .dacsups = nullptr,
  .tauzq = nullptr,
  .tauz = nullptr,
  .zq = nullptr,
  .ir_div = nullptr,
  .sol_div = nullptr,
  .iru = nullptr,
  .ird = nullptr,
  .solu = nullptr,
  .sold = nullptr,
  .qw = nullptr,
  .ql = nullptr,
  .qv = nullptr,
  .nc = nullptr,
  .fn = nullptr,
  .aerosols = nullptr,
  .albedo0 = nullptr,
  .emissi0 = nullptr,
  .temp0   = nullptr,
  .theta0  = nullptr,
  .qw0     = nullptr,
  .p0      = nullptr,
  .rho0    = nullptr
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_atmo_1d_rad_t *cs_glob_atmo_1d_rad = &_atmo_1d_rad;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_atmo_get_radiative_model_1d(int  **radiative_model_1d,
                                 int  **nvertint,
                                 int  **kvert,
                                 int  **kmx,
                                 int  **nfatr1);

void
cs_f_atmo_arrays_get_radiative_model_1d(cs_real_t **xyvert,
                                        cs_real_t **zray,
                                        cs_real_t **acinfe,
                                        cs_real_t **dacinfe,
                                        cs_real_t **aco2,
                                        cs_real_t **aco2s,
                                        cs_real_t **daco2,
                                        cs_real_t **daco2s,
                                        cs_real_t **acsup,
                                        cs_real_t **acsups,
                                        cs_real_t **dacsup,
                                        cs_real_t **dacsups,
                                        cs_real_t **tauzq,
                                        cs_real_t **tauz,
                                        cs_real_t **zq,
                                        cs_real_t **rayi,
                                        cs_real_t **rayst,
                                        cs_real_t **iru,
                                        cs_real_t **ird,
                                        cs_real_t **solu,
                                        cs_real_t **sold,
                                        cs_real_t **ground_albedo,
                                        cs_real_t **ground_emissi,
                                        cs_real_t **ground_ttground,
                                        cs_real_t **ground_tpground,
                                        cs_real_t **ground_totwat,
                                        cs_real_t **ground_pressure,
                                        cs_real_t **ground_density,
                                        int         dim_xyvert[2],
                                        int         dim_kmx2[2],
                                        int         dim_kmx_nvert[2]);

void
cs_f_atmo_rad_1d_arrays_get_pointers(cs_real_t **qwvert,
                                     cs_real_t **qlvert,
                                     cs_real_t **qvvert,
                                     cs_real_t **ncvert,
                                     cs_real_t **fnvert,
                                     cs_real_t **aevert);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Access pointers for Fortran mapping.
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_radiative_model_1d(int  **radiative_model_1d,
                                 int  **nvert,
                                 int  **kvert,
                                 int  **kmx,
                                 int  **nfatr1)
{
  *radiative_model_1d = &(_atmo_1d_rad.radiative_model_1d);
  *nvert = &(_atmo_1d_rad.nvert);
  *kvert = &(_atmo_1d_rad.nlevels);
  *kmx = &(_atmo_1d_rad.nlevels_max);
  *nfatr1 = &(_atmo_1d_rad.frequency);
}

void
cs_f_atmo_rad_1d_arrays_get_pointers(cs_real_t **qwvert,
                                     cs_real_t **qlvert,
                                     cs_real_t **qvvert,
                                     cs_real_t **ncvert,
                                     cs_real_t **fnvert,
                                     cs_real_t **aevert)
{
  *qwvert = _atmo_1d_rad.qw;
  *qlvert = _atmo_1d_rad.ql;
  *qvvert = _atmo_1d_rad.qv;
  *ncvert = _atmo_1d_rad.nc;
  *fnvert = _atmo_1d_rad.fn;
  *aevert = _atmo_1d_rad.aerosols;
}

void
cs_f_atmo_arrays_get_radiative_model_1d(cs_real_t **xyvert,
                                        cs_real_t **zray,
                                        cs_real_t **acinfe,
                                        cs_real_t **dacinfe,
                                        cs_real_t **aco2,
                                        cs_real_t **aco2s,
                                        cs_real_t **daco2,
                                        cs_real_t **daco2s,
                                        cs_real_t **acsup,
                                        cs_real_t **acsups,
                                        cs_real_t **dacsup,
                                        cs_real_t **dacsups,
                                        cs_real_t **tauzq,
                                        cs_real_t **tauz,
                                        cs_real_t **zq,
                                        cs_real_t **rayi,
                                        cs_real_t **rayst,
                                        cs_real_t **iru,
                                        cs_real_t **ird,
                                        cs_real_t **solu,
                                        cs_real_t **sold,
                                        cs_real_t **ground_albedo,
                                        cs_real_t **ground_emissi,
                                        cs_real_t **ground_ttground,
                                        cs_real_t **ground_tpground,
                                        cs_real_t **ground_totwat,
                                        cs_real_t **ground_pressure,
                                        cs_real_t **ground_density,
                                        int         dim_xyvert[2],
                                        int         dim_kmx2[2],
                                        int         dim_kmx_nvert[2])
{
  int n_vert = 0;
  int n_level = 0;

  if (_atmo_1d_rad.radiative_model_1d != 0) {
    n_level = cs::max(1, _atmo_1d_rad.nlevels_max);
    n_vert = cs::max(1, _atmo_1d_rad.nvert);
  }

  if (         _atmo_1d_rad.xy == nullptr)
    CS_MALLOC(_atmo_1d_rad.xy , 3*n_vert, cs_real_t);
  if (         _atmo_1d_rad.z == nullptr)
    CS_MALLOC(_atmo_1d_rad.z , n_level, cs_real_t);
  if (         _atmo_1d_rad.acinfe == nullptr)
    CS_MALLOC(_atmo_1d_rad.acinfe , n_level, cs_real_t);
  if (         _atmo_1d_rad.dacinfe == nullptr)
    CS_MALLOC(_atmo_1d_rad.dacinfe , n_level, cs_real_t);
  if (         _atmo_1d_rad.aco2 == nullptr)
    CS_MALLOC(_atmo_1d_rad.aco2, n_level*n_level, cs_real_t);
  if (         _atmo_1d_rad.aco2s == nullptr)
    CS_MALLOC(_atmo_1d_rad.aco2s, n_level*n_level, cs_real_t);
  if (         _atmo_1d_rad.daco2 == nullptr)
    CS_MALLOC(_atmo_1d_rad.daco2, n_level*n_level, cs_real_t);
  if (         _atmo_1d_rad.daco2s == nullptr)
    CS_MALLOC(_atmo_1d_rad.daco2s, n_level*n_level, cs_real_t);
  if (         _atmo_1d_rad.acsup == nullptr)
    CS_MALLOC(_atmo_1d_rad.acsup, n_level, cs_real_t);
  if (         _atmo_1d_rad.acsups == nullptr)
    CS_MALLOC(_atmo_1d_rad.acsups, n_level, cs_real_t);
  if (         _atmo_1d_rad.dacsup == nullptr)
    CS_MALLOC(_atmo_1d_rad.dacsup, n_level, cs_real_t);
  if (         _atmo_1d_rad.dacsups == nullptr)
    CS_MALLOC(_atmo_1d_rad.dacsups, n_level, cs_real_t);
  if (         _atmo_1d_rad.tauzq == nullptr)
    CS_MALLOC(_atmo_1d_rad.tauzq, n_level+1, cs_real_t);
  if (         _atmo_1d_rad.tauz == nullptr)
    CS_MALLOC(_atmo_1d_rad.tauz, n_level+1, cs_real_t);
  if (         _atmo_1d_rad.zq == nullptr)
    CS_MALLOC(_atmo_1d_rad.zq, n_level+1, cs_real_t);
  if (         _atmo_1d_rad.ir_div == nullptr)
    CS_MALLOC(_atmo_1d_rad.ir_div, n_level * n_vert, cs_real_t);
  if (         _atmo_1d_rad.sol_div == nullptr)
    CS_MALLOC(_atmo_1d_rad.sol_div, n_level * n_vert, cs_real_t);
  if (         _atmo_1d_rad.iru == nullptr)
    CS_MALLOC(_atmo_1d_rad.iru, n_level * n_vert, cs_real_t);
  if (         _atmo_1d_rad.ird == nullptr) {
    CS_MALLOC(_atmo_1d_rad.ird, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.ird);
  }
  if (         _atmo_1d_rad.solu == nullptr)
    CS_MALLOC(_atmo_1d_rad.solu, n_level * n_vert, cs_real_t);
  if (         _atmo_1d_rad.sold == nullptr) {
    CS_MALLOC(_atmo_1d_rad.sold, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.sold);
  }
  if (         _atmo_1d_rad.qw == nullptr) {
    CS_MALLOC(_atmo_1d_rad.qw, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.qw);
  }
  if (         _atmo_1d_rad.ql == nullptr) {
    CS_MALLOC(_atmo_1d_rad.ql, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.ql);
  }
  if (         _atmo_1d_rad.qv == nullptr) {
    CS_MALLOC(_atmo_1d_rad.qv, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.qv);
  }
  if (         _atmo_1d_rad.nc == nullptr) {
    CS_MALLOC(_atmo_1d_rad.nc, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.nc);
  }
  if (         _atmo_1d_rad.fn == nullptr) {
    CS_MALLOC(_atmo_1d_rad.fn, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.fn);
  }
  if (         _atmo_1d_rad.aerosols == nullptr) {
    CS_MALLOC(_atmo_1d_rad.aerosols, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_1d_rad.aerosols);
  }
  if (         _atmo_1d_rad.albedo0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.albedo0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.albedo0);
  }
  if (         _atmo_1d_rad.emissi0== nullptr) {
    CS_MALLOC(_atmo_1d_rad.emissi0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.emissi0);
  }
  if (         _atmo_1d_rad.temp0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.temp0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.temp0);
  }
  if (         _atmo_1d_rad.theta0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.theta0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.theta0);
  }
  if (         _atmo_1d_rad.qw0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.qw0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.qw0);
  }
  if (         _atmo_1d_rad.p0  == nullptr) {
    CS_MALLOC(_atmo_1d_rad.p0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.p0);
  }
  if (         _atmo_1d_rad.rho0 == nullptr) {
    CS_MALLOC(_atmo_1d_rad.rho0, n_vert, cs_real_t);
    cs_array_real_fill_zero(n_vert, _atmo_1d_rad.rho0);
  }

  *xyvert = _atmo_1d_rad.xy;
  *zray   = _atmo_1d_rad.z;
  *acinfe = _atmo_1d_rad.acinfe;
  *dacinfe = _atmo_1d_rad.dacinfe;
  *aco2   = _atmo_1d_rad.aco2  ;
  *aco2s  = _atmo_1d_rad.aco2s ;
  *daco2  = _atmo_1d_rad.daco2 ;
  *daco2s = _atmo_1d_rad.daco2s;
  *acsup  = _atmo_1d_rad.acsup ;
  *acsups = _atmo_1d_rad.acsups;
  *dacsup = _atmo_1d_rad.dacsup;
  *dacsups = _atmo_1d_rad.dacsups;
  *tauzq  = _atmo_1d_rad.tauzq ;
  *tauz   = _atmo_1d_rad.tauz  ;
  *zq     = _atmo_1d_rad.zq    ;
  *rayi   = _atmo_1d_rad.ir_div  ;
  *rayst  = _atmo_1d_rad.sol_div ;
  *iru    = _atmo_1d_rad.iru   ;
  *ird    = _atmo_1d_rad.ird   ;
  *solu   = _atmo_1d_rad.solu  ;
  *sold   = _atmo_1d_rad.sold  ;

  /* ground level arrays, of size n_vert */
  *ground_albedo   = _atmo_1d_rad.albedo0;
  *ground_emissi   = _atmo_1d_rad.emissi0;
  *ground_ttground   = _atmo_1d_rad.temp0;
  *ground_tpground   = _atmo_1d_rad.theta0;
  *ground_totwat   = _atmo_1d_rad.qw0;
  *ground_pressure = _atmo_1d_rad.p0;
  *ground_density  = _atmo_1d_rad.rho0;

  dim_kmx2[0]      = _atmo_1d_rad.nlevels_max;
  dim_kmx2[1]      = _atmo_1d_rad.nlevels_max;
  dim_xyvert[0]    = _atmo_1d_rad.nvert;
  dim_xyvert[1]    = 3;
  dim_kmx_nvert[0] = _atmo_1d_rad.nlevels_max;
  dim_kmx_nvert[1] = _atmo_1d_rad.nvert;

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief free arrays for atmo 1-D radiative module
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_1d_rad_finalize(void)
{
  CS_FREE(_atmo_1d_rad.xy);
  CS_FREE(_atmo_1d_rad.z);
  CS_FREE(_atmo_1d_rad.acinfe);
  CS_FREE(_atmo_1d_rad.dacinfe);
  CS_FREE(_atmo_1d_rad.aco2);
  CS_FREE(_atmo_1d_rad.aco2s);
  CS_FREE(_atmo_1d_rad.daco2);
  CS_FREE(_atmo_1d_rad.daco2s);
  CS_FREE(_atmo_1d_rad.acsup);
  CS_FREE(_atmo_1d_rad.acsups);
  CS_FREE(_atmo_1d_rad.dacsup);
  CS_FREE(_atmo_1d_rad.dacsups);
  CS_FREE(_atmo_1d_rad.tauzq);
  CS_FREE(_atmo_1d_rad.tauz);
  CS_FREE(_atmo_1d_rad.zq);
  CS_FREE(_atmo_1d_rad.ir_div);
  CS_FREE(_atmo_1d_rad.sol_div);
  CS_FREE(_atmo_1d_rad.iru);
  CS_FREE(_atmo_1d_rad.ird);
  CS_FREE(_atmo_1d_rad.solu);
  CS_FREE(_atmo_1d_rad.sold);
  CS_FREE(_atmo_1d_rad.albedo0);
  CS_FREE(_atmo_1d_rad.emissi0);
  CS_FREE(_atmo_1d_rad.temp0);
  CS_FREE(_atmo_1d_rad.theta0);
  CS_FREE(_atmo_1d_rad.qw0);
  CS_FREE(_atmo_1d_rad.p0);
  CS_FREE(_atmo_1d_rad.rho0);
  CS_FREE(_atmo_1d_rad.aerosols);
  CS_FREE(_atmo_1d_rad.fn);
  CS_FREE(_atmo_1d_rad.nc);
  CS_FREE(_atmo_1d_rad.qv);
  CS_FREE(_atmo_1d_rad.ql);
  CS_FREE(_atmo_1d_rad.qw);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
