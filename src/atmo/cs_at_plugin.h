#ifndef __CS_ATPLUGIN_H__
#define __CS_ATPLUGIN_H__

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_at_plugin.h

  \brief Plugin to get aerosol and compute_coagulation_coefficient functions
    from SIREAM library (ENPC - INRIA - EDF R&D)
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/* Plug-in to get aerosol function
   from SIREAM library (ENPC - INRIA - EDF R&D) */

void CS_PROCF(plug_aerosol, PLUG_AEROSOL)
(
 cs_int_t   *nx,
 cs_int_t   *ny,
 cs_int_t   *nz,
 cs_int_t   *ns,
 cs_real_t  *ts,
 cs_real_t  *dlhumid,
 cs_real_t  *dltemp,
 cs_real_t  *dlpress,
 cs_real_t  *delta_t,
 cs_real_t  *dlconc,
 cs_int_t   *noptions_aer,
 cs_int_t   *options_aer,
 cs_int_t   *ns_aer,
 cs_int_t   *nbin_aer,
 cs_int_t   *ncycle_aer,
 cs_real_t  *bin_bound_aer,
 cs_real_t  *fixed_density_aer,
 cs_real_t  *density_aer,
 cs_int_t   *couples_coag,
 cs_int_t   *first_index_coag,
 cs_int_t   *second_index_coag,
 cs_real_t  *coefficient_coag,
 cs_real_t  *dlconc_aer,
 cs_real_t  *dlnum_aer
);

/* Plug-in to get compute_coagulation_coefficient function
   from SIREAM library (ENPC - INRIA - EDF R&D) */

void CS_PROCF(plug_compute_coagulation_coefficient,
              PLUG_COMPUTE_COAGULATION_COEFFICIENT)
(
 cs_int_t   *nbin_aer,
 cs_real_t  *bin_bound,
 cs_int_t   *couple,
 cs_int_t   *first_index,
 cs_int_t   *second_index,
 cs_real_t  *partition_coefficient
 );

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATPLUGIN_H__ */
