#ifndef __CS_SEARCH_H__
#define __CS_SEARCH_H__

/*============================================================================
 * Search elements in arrays
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search. Array or list must be ordered.
 *
 * parameters:
 *   size   <--  size of list
 *   gnum   <--  find index for this number
 *   lst    <--  list of ordered global numbers to scan
 *
 * returns:
 *   id associated to the current number. If not found, returned -1.
 *---------------------------------------------------------------------------*/

int
cs_search_g_binary(size_t             size,
                   cs_gnum_t          gnum,
                   const cs_gnum_t    lst[]);

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search (binary search). Array or list must be ordered.
 *
 * parameters:
 *   size   <--  size of list
 *   num    <--  find index for this number
 *   lst    <--  list of ordered numbers to scan
 *
 * returns:
 *   id associated to the current number. If not found, return -1.
 *---------------------------------------------------------------------------*/

int
cs_search_binary(size_t           size,
                 cs_lnum_t        num,
                 const cs_lnum_t  lst[]);

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search (binary search). Index must be ordered and without null range.
 *
 * parameters:
 *   size   <--  size of index -1
 *   gnum   <--  number for which we want the position in index
 *   index  <--  index array
 *
 * returns:
 *   id in index of gnum.  If not found, returned -1.
 *---------------------------------------------------------------------------*/

int
cs_search_gindex_binary(size_t             size,
                        cs_gnum_t          gnum,
                        const cs_gnum_t    index[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SEARCH_H__ */
