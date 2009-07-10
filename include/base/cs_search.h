/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2009 EDF S.A., France
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

#ifndef __CS_SEARCH_H__
#define __CS_SEARCH_H__

/*============================================================================
 * Search elements in arrays
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

#include <fvm_defs.h>

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
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search. Array or list must be ordered.
 *
 * parameters:
 *   start  <--  start search from this index
 *   end    <--  end search to this index
 *   gnum   <--  find index for this number
 *   lst    <--  list of ordered global numbers to scan
 *
 * returns:
 *   id associated to the current number. If not found, returned -1.
 *---------------------------------------------------------------------------*/

int
cs_search_g_binary(int                start,
                   int                end,
                   fvm_gnum_t         gnum,
                   const fvm_gnum_t   lst[]);

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search. Array or list must be ordered.
 *
 * parameters:
 *   start  <--  start search from this index
 *   end    <--  end search to this index
 *   num    <--  find index for this number
 *   lst    <--  list of ordered numbers to scan
 *
 * returns:
 *   id associated to the current number. If not found, returned -1.
 *---------------------------------------------------------------------------*/

int
cs_search_binary(int              start,
                 int              end,
                 cs_int_t         num,
                 const cs_int_t   lst[]);

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search. Index must be ordered and without null range.
 *
 * parameters:
 *   start    <--  start search from this index
 *   end      <--  end search to this index
 *   gnum     <--  number for which we want the position in index
 *   index    <--  index array
 *
 * returns:
 *   id in index of gnum.  If not found, returned -1.
 *---------------------------------------------------------------------------*/

int
cs_search_gindex_binary(int                start,
                        int                end,
                        fvm_gnum_t         gnum,
                        const fvm_gnum_t   index[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SEARCH_H__ */
