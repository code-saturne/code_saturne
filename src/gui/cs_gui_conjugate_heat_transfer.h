#ifndef __CS_GUI_CONJUGATE_HEAT_TRANSFER_H__
#define __CS_GUI_CONJUGATE_HEAT_TRANSFER_H__

/*============================================================================
 * Management of the GUI parameters file: conjugate heat transfer
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

BEGIN_C_DECLS

/*============================================================================
 * Public Function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define SYRTHES couplings based on the GUI-generated setup.
 *
 * In the case of a single code_saturne and single SYRTHES instance, the
 * syrthes_name argument is ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument.
 */
/*----------------------------------------------------------------------------*/

void cs_gui_syrthes_coupling(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_CONJUGATE_HEAT_TRANSFER_H__ */
