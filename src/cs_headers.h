#ifndef __CS_HEADERS_H__
#define __CS_HEADERS_H__

/*============================================================================
 * Global code_saturne headers file for easier include
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/* Disable some compiler warnings */

#if defined __INTEL_COMPILER
#pragma warning disable 2259
#endif

/* Include headers by groups */

#include "bft/cs_bft_headers.h"
#include "base/cs_base_headers.h"
#include "fvm/fvm_headers.h"
#include "mesh/cs_mesh_headers.h"
#include "alge/cs_alge_headers.h"
#include "atmo/cs_atmo_headers.h"
#include "cfbl/cs_cfbl_headers.h"
#include "cdo/cs_cdo_headers.h"
#include "cogz/cs_cogz_headers.h"
#include "comb/cs_comb_headers.h"
#include "ctwr/cs_ctwr_headers.h"
#include "elec/cs_elec_headers.h"
#include "gui/cs_gui_headers.h"
#include "gwf/cs_gwf_headers.h"
#include "lagr/cs_lagr_headers.h"
#include "meg/cs_meg_headers.h"
#include "pprt/cs_pprt_headers.h"
#include "rayt/cs_rad_headers.h"
#include "turb/cs_turbulence_headers.h"

/*----------------------------------------------------------------------------*/

#if defined(USE_NEPTUNE_CFD)
#include "nc_c_headers.h"
#endif

/*----------------------------------------------------------------------------*/

#endif /* __CS_HEADERS_H__ */
