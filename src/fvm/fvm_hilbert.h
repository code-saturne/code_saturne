#ifndef __FVM_HILBERT_H__
#define __FVM_HILBERT_H__

/*============================================================================
 * Hilbert space-filling curve construction for coordinates.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/* Hilbert code
   (could be switched from double to long double for extended range) */

typedef double  fvm_hilbert_code_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine the global extents associated with a set of coordinates
 *
 * parameters:
 *   dim       <-- spatial dimension
 *   n_coords  <-- local number of coordinates
 *   coords    <-- entity coordinates; size: n_entities*dim (interlaced)
 *   g_extents --> global extents (size: dim*2)
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void
fvm_hilbert_get_coord_extents(int               dim,
                              size_t            n_coords,
                              const cs_coord_t  coords[],
                              cs_coord_t        g_extents[],
                              MPI_Comm          comm);
#else
void
fvm_hilbert_get_coord_extents(int               dim,
                              size_t            n_coords,
                              const cs_coord_t  coords[],
                              cs_coord_t        g_extents[]);
#endif

/*----------------------------------------------------------------------------
 * Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   h_code   --> array of corresponding Hilbert codes (size: n_coords)
 *----------------------------------------------------------------------------*/

void
fvm_hilbert_encode_coords(int                 dim,
                          const cs_coord_t    extents[],
                          cs_lnum_t           n_coords,
                          const cs_coord_t    coords[],
                          fvm_hilbert_code_t  h_code[]);

/*----------------------------------------------------------------------------
 * Locally order a list of Hilbert ids.
 *
 * This variant uses an encoding into floating-point numbers. In 3D, this
 * limits us to 19 levels for a double, though using a long double could
 * increase this range.
 *
 * parameters:
 *   n_codes       <-- number of Hilbert ids to order
 *   hilbert_codes <-- array of Hilbert ids to order
 *   order         --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
fvm_hilbert_local_order(cs_lnum_t                 n_codes,
                        const fvm_hilbert_code_t  hilbert_codes[],
                        cs_lnum_t                 order[]);

/*----------------------------------------------------------------------------
 * Locally order a list of coordinates based on their Hilbert code.
 *
 * This variant may use a maximum depth of 32 levels, and switches
 * to lexicographical ordering if this is not enough.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   order    --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
fvm_hilbert_local_order_coords(int                dim,
                               const cs_coord_t   extents[],
                               cs_lnum_t          n_coords,
                               const cs_coord_t   coords[],
                               cs_lnum_t          order[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for conversion of a double precision value in
 * range [0, 1] to a given Hilbert code.
 *
 * \param[in]   s      coordinate between 0 and 1
 * \param[out]  elt    pointer to element
 * \param[in]   input  pointer to optional (untyped) value or structure.
 */
/*----------------------------------------------------------------------------*/

void
fvm_hilbert_s_to_code(double       s,
                      void        *elt,
                      const void  *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for comparison of 2 Hilbert codes.
 *
 * This function is the same type as that used by qsort_r.
 *
 * \param[in]  elt1   coordinate between 0 and 1
 * \param[in]  elt2   pointer to optional (untyped) value or structure.
 * \param[in]  input  pointer to optional (untyped) value or structure.
 *
 * \return < 0 if elt1 < elt2, 0 if elt1 == elt2, > 0 if elt1 > elt2
 */
/*----------------------------------------------------------------------------*/

int
fvm_hilbert_compare(const void  *elt1,
                    const void  *elt2,
                    const void  *input);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_HILBERT_H__ */
