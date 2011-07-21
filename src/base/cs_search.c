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

/*============================================================================
 * Subroutines to search elements in arrays
 *===========================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *---------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_search.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search (binary search). Array or list must be ordered.
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

inline static int
_binary_search(int              start,
               int              end,
               cs_int_t         num,
               const cs_int_t   lst[])
{
  if (lst[start] == num)
    return start;

  else if (lst[end] == num)
    return end;

  else {

    int  range = (end - start)/2;
    int  middle = start + range;

    assert(middle <= end);

    if (range == 0)
      return -1;

    if ( lst[middle] > num )
      return _binary_search(start, middle, num, lst);
    else
      return _binary_search(middle, end, num, lst);

  }

}

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search (binary search). Array or list must be ordered.
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

inline static int
_binary_gsearch(int                start,
                int                end,
                cs_gnum_t          gnum,
                const cs_gnum_t    lst[])
{
  if (lst[start] == gnum)
    return start;

  else if (lst[end] == gnum)
    return end;

  else {

    int  range = (end - start)/2;
    int  middle = start + range;

    if (range == 0)
      return -1;

    if ( lst[middle] > gnum )
      return _binary_gsearch(start, middle, gnum, lst);
    else
      return _binary_gsearch(middle, end, gnum, lst);

  }

}

/*----------------------------------------------------------------------------
 * Get the position inside an array related to a value thanks to a binary
 * search (binary search). Index must be ordered.
 *
 * parameters:
 *   start    <--  start search from this index
 *   end      <--  end search to this index
 *   gnum     <--  number for which we want the position in index
 *   index    <--  index array
 *
 * returns:
 *   id in index of gnum.  If not found, returned -1.
 *----------------------------------------------------------------------------*/

inline static int
_binary_index_gsearch(int                 start,
                      int                 end,
                      cs_gnum_t           gnum,
                      const cs_gnum_t     index[])
{
  if (end - start < 2)
    return start;

  else {

    int  middle_id = (end - start)/2 + start;

    if (gnum < index[middle_id])
      return _binary_index_gsearch(start, middle_id, gnum,index);
    else
      return _binary_index_gsearch(middle_id, end, gnum, index);

  }

}

/*============================================================================
 * Public function definitions
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
                   const cs_gnum_t    lst[])
{
  return _binary_gsearch(0, size - 1, gnum, lst);
}

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
                 cs_int_t         num,
                 const cs_int_t   lst[])
{
  return _binary_search(0, size - 1, num, lst);
}

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
                        const cs_gnum_t    index[])
{
  return _binary_index_gsearch(0, size, gnum, index);
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
