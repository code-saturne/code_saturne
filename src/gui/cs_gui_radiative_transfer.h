#ifndef __CS_GUI_RADIATIVE_TRANSFER_H__
#define __CS_GUI_RADIATIVE_TRANSFER_H__

/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Free GUI-defined radiative transfer parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfers_finalize(void);

/*----------------------------------------------------------------------------
 * Read GUI-defined radiative transfer parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_parameters(void);

/*----------------------------------------------------------------------------
 * Set the radiative transfer absorption coefficient
 *
 * parameters:
 *   ck --> absorption coefficient at cells
 *----------------------------------------------------------------------------*/

void
cs_gui_rad_transfer_absorption(cs_real_t  ck[]);

/*----------------------------------------------------------------------------
 * Postprocessing settings for radiative transfer
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_postprocess(void);

/*----------------------------------------------------------------------------
 * Radiative transfer boundary conditions
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_bcs(const    int   itypfb[],
                              int            nvar,
                              int            ivart,
                              int           *isothp,
                              double        *epsp,
                              double        *epap,
                              double        *tintp,
                              double        *textp,
                              double        *xlamp,
                              double        *rcodcl);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_RADIATIVE_TRANSFER_H__ */
