/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_GUI_CONJUGATE_HEAT_TRANSFER_H__
#define __CS_GUI_CONJUGATE_HEAT_TRANSFER_H__

/*============================================================================
 * Management of the GUI parameters file: conjugate heat transfer
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

BEGIN_C_DECLS

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * syrthes_app_num and syrthes_name arguments are ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument, then
 * on the syrthes_app_num argument. If syrthes_name is empty, matching will
 * be based on syrthes_app_num only.
 *
 * subroutine uisyrc
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uisyrc, UISYRC) (void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_CONJUGATE_HEAT_TRANSFER_H__ */
