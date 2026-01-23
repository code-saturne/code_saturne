#ifndef CS_REDISTRIBUTE_H
#define CS_REDISTRIBUTE_H

/*============================================================================
 * Repartitioning and redistribution mesh data and fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Redistribute mesh and fields based on a cell destination rank map.
 *
 * If no cell map is given, a random one is created.
 *
 * \param[in]  cell_dest_rank  destination rank for each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_redistribute(int  cell_dest_rank[]);

/*----------------------------------------------------------------------------*/

#endif /* CS_REDISTRIBUTE_H */
