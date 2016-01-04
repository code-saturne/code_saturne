#ifndef __FVM_HILBERT_H__
#define __FVM_HILBERT_H__

/*============================================================================
 * Hilbert space-filling curve construction for coordinates.
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
   (could be switched from double to long double for extened range) */

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

/*----------------------------------------------------------------------------
 * Get the quantile associated to a Hilbert code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * parameters:
 *   n_quantiles    <-- number of quantiles
 *   code           <-- code we are searching for
 *   quantile_start <-- first Hilbert code in each quantile (size: n_quantiles)
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

size_t
fvm_hilbert_quantile_search(size_t              n_quantiles,
                            fvm_hilbert_code_t  code,
                            fvm_hilbert_code_t  quantile_start[]);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build a global Hilbert encoding rank index.
 *
 * The rank_index[i] contains the first Hilbert code assigned to rank [i].
 *
 * parameters:
 *   dim          <-- 1D, 2D or 3D
 *   n_codes      <-- number of Hilbert codes to be indexed
 *   hilbert_code <-- array of Hilbert codes to be indexed
 *   weight       <-- weighting related to each code
 *   order        <-- ordering array
 *   rank_index   <-> pointer to the global Hilbert encoding rank index
 *   comm         <-- MPI communicator on which we build the global index
 *
 * returns:
 *  the fit related to the Hilbert encoding distribution (lower is better).
 *----------------------------------------------------------------------------*/

double
fvm_hilbert_build_rank_index(int                       dim,
                             cs_lnum_t                 n_codes,
                             const fvm_hilbert_code_t  hilbert_code[],
                             const cs_lnum_t           weight[],
                             const cs_lnum_t           order[],
                             fvm_hilbert_code_t        rank_index[],
                             MPI_Comm                  comm);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_HILBERT_H__ */
