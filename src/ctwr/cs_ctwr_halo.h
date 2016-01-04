#ifndef __CS_CTWR_HALO_H__
#define __CS_CTWR_HALO_H__

/*============================================================================
 * Structure and function headers handling with ghost cells
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_interface.h"

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*---------------------------------------------------------------------------
 * Reverse "ghost cells -> vertex" connectivity into "vertex -> ghost cells"
 * connectivity for out_halo elements.
 * Build the connectivity list.
 *
 * parameters:
 *---------------------------------------------------------------------------*/

void
cs_reverse_vtx_faces_connect(const fvm_nodal_t   *this_nodal,
                             cs_int_t   *faces_vtx_idx[],
                             cs_int_t   *faces_vtx_lst[]);

/*----------------------------------------------------------------------------
 * Define halo structures for internal and distant ghost cells.
 *
 * parameters:
 *   mesh                 -->  pointer to cs_mesh_t structure
 *   interface_set        -->  pointer to cs_interface_set_t structure.
 *---------------------------------------------------------------------------*/

void
cs_ctwr_halo_define(cs_ctwr_zone_t      *ct,
                    cs_interface_set_t  *interface_set);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_HALO_H__ */
