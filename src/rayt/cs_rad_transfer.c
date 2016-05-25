/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_timer.h"

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_rad_transfer.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

cs_rad_transfer_params_t _rt_params = {.iirayo = 0,
                                       .nrphas = 0,
                                       .iimpar = 0,
                                       .iimlum = 0,
                                       .imodak = 0,
                                       .imoadf = 0,
                                       .iwrp1t = 0,
                                       .imfsck = 0,
                                       .xnp1mx = 0.,
                                       .idiver = 0,
                                       .i_quadrature = 0,
                                       .ndirec = 0,
                                       .ndirs = 0,
                                       .sxyz = NULL,
                                       .angsol = NULL,
                                       .restart = 0,
                                       .nfreqr = 0,
                                       .nwsgg = 0,
                                       .wq = NULL,
                                       .nbzrdm = 2000,
                                       .nozrdm = 2000,
                                       .nozarm = 0,
                                       .nzfrad = 0,
                                       .ilzrad = NULL,
                                       .itpimp = 1,
                                       .ipgrno = 21,
                                       .iprefl = 22,
                                       .ifgrno = 31,
                                       .ifrefl = 32,
                                       .itpt1d = 4 };

cs_rad_transfer_params_t *cs_glob_rad_transfer_params = &_rt_params;

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_rad_transfer_get_pointers(int  **p_iirayo,
                             int  **p_nfreqr);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_rad_transfer_get_pointers(int  **p_iirayo,
                             int  **p_nfreqr)
{
  *p_iirayo = &_rt_params.iirayo;
  *p_nfreqr = &_rt_params.nfreqr;
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize radiative transfer module.
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_finalize(void)
{
  BFT_FREE(_rt_params.sxyz);
  BFT_FREE(_rt_params.angsol);
  BFT_FREE(_rt_params.wq);
  BFT_FREE(_rt_params.ilzrad);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
