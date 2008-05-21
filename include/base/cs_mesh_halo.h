/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_MESH_HALO_H__
#define __CS_MESH_HALO_H__

/*============================================================================
 * Structure and function headers handling with ghost cells
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_interface.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define halo structures for internal and distant ghost cells.
 *
 * parameters:
 *   mesh             -->  pointer to cs_mesh_t structure
 *   interface_set    -->  pointer to fvm_interface_set_t structure.
 *   p_gcell_vtx_idx  <--  pointer to the connectivity index
 *   p_gcell_vtx_lst  <--  pointer to the connectivity list
 *---------------------------------------------------------------------------*/

void
cs_mesh_halo_define(cs_mesh_t            *mesh,
                    fvm_interface_set_t  *interface_set,
                    cs_int_t             *p_gcell_vtx_idx[],
                    cs_int_t             *p_gcell_vtx_lst[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_MESH_HALO_H__ */
