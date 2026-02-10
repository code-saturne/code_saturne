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

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {
  int *c_r_level;
  int *c_r_flag;
} cs_redistribute_data_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

extern cs_redistribute_data_t  *cs_glob_redistribute_data;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Redistribute mesh, fields and data based on a cell destination rank
 * map.
 *
 * If no cell map is given, a random one is created.
 * If pointer to data is not null, this function will try to redistribute all
 * its internal buffers.
 *
 * \param[in]     cell_dest_rank  destination rank for each cell
 * \param[inout]  data            pointer to redistribution data structure
 */
/*----------------------------------------------------------------------------*/

void
cs_redistribute(const int                cell_dest_rank[],
                cs_redistribute_data_t  *data);

/*----------------------------------------------------------------------------*/

#endif /* CS_REDISTRIBUTE_H */
