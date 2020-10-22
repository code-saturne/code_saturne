#ifndef __CS_CLIP_KE_H__
#define __CS_CLIP_KE_H__

/*============================================================================
 * Clipping of the turbulent kinetic energy and the turbulent
 * dissipation.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of the turbulent kinetic energy and turbulent dissipation.
 *
 * \param[in]     n_cells  number of cells
 * \param[in]     iclip    indicator = 0 if viscl0 is used
 *                         otherwise viscl is used.
 */
/*----------------------------------------------------------------------------*/

void
cs_clip_ke(cs_lnum_t  n_cells,
           int        iclip);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CLIP_KE_H__ */
