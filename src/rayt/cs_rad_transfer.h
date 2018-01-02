#ifndef __CS_RAD_TRANSFER_H__
#define __CS_RAD_TRANSFER_H__

/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

typedef struct {

  int           iirayo;
  int           nrphas;
  int           iimpar;
  int           iimlum;
  int           imodak;
  int           imoadf;
  int           iwrp1t;
  int           imfsck;
  double        xnp1mx;
  int           idiver;
  int           i_quadrature;
  int           ndirec;
  int           ndirs;
  cs_real_3_t  *sxyz;
  cs_real_t    *angsol;
  int           restart;
  int           nfreqr;
  int           nwsgg;
  cs_real_t    *wq;
  int           nzfrad;
  int           itpimp;
  int           ipgrno;
  int           iprefl;
  int           ifgrno;
  int           ifrefl;
  int           itpt1d;

} cs_rad_transfer_params_t;

extern cs_rad_transfer_params_t *cs_glob_rad_transfer_params;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize radiative transfer module.
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_H__ */
