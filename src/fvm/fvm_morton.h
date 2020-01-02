#ifndef __FVM_MORTON_H__
#define __FVM_MORTON_H__

/*============================================================================
 * Morton encoding for 2D or 3D coordinates.
 *============================================================================*/

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef enum {

  FVM_MORTON_EQUAL_ID,
  FVM_MORTON_SAME_ANCHOR,
  FVM_MORTON_DIFFERENT_ID

} fvm_morton_compare_t;

typedef unsigned int   fvm_morton_int_t;

typedef struct {

  fvm_morton_int_t   L;     /* Level in the tree structure */
  fvm_morton_int_t   X[3];  /* X, Y, Z coordinates in Cartesian grid */

} fvm_morton_code_t;

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
fvm_morton_get_coord_extents(int               dim,
                             size_t            n_coords,
                             const cs_coord_t  coords[],
                             cs_coord_t        g_extents[],
                             MPI_Comm          comm);

#else

void
fvm_morton_get_coord_extents(int               dim,
                             size_t            n_coords,
                             const cs_coord_t  coords[],
                             cs_coord_t        g_extents[]);

#endif

/*----------------------------------------------------------------------------
 * Determine the global extents associated with a set of local extents
 *
 * parameters:
 *   dim       <-- spatial dimension
 *   n_extents <-- local number of coordinates
 *   extents   <-- entity coordinates; size: n_entities*dim*2 (interlaced)
 *   g_extents --> global extents (size: dim*2)
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
fvm_morton_get_global_extents(int               dim,
                              size_t            n_extents,
                              const cs_coord_t  extents[],
                              cs_coord_t        g_extents[],
                              MPI_Comm          comm);

#else

void
fvm_morton_get_global_extents(int               dim,
                              size_t            n_extents,
                              const cs_coord_t  extents[],
                              cs_coord_t        g_extents[]);

#endif

/*----------------------------------------------------------------------------
 * Build a Morton code according to the level in an octree grid and its
 * coordinates in the grid.
 *
 * parameters:
 *   dim    <-- 1D, 2D or 3D
 *   level  <-- level in the grid
 *   coords <-- coordinates in the grid (normalized)
 *
 * returns:
 *  a Morton code
 *----------------------------------------------------------------------------*/

fvm_morton_code_t
fvm_morton_encode(int                dim,
                  fvm_morton_int_t   level,
                  const cs_coord_t   coords[]);

/*----------------------------------------------------------------------------
 * Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   level    <-- level in the grid
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   m_code   --> array of corresponding Morton codes
 *----------------------------------------------------------------------------*/

void
fvm_morton_encode_coords(int                dim,
                         fvm_morton_int_t   level,
                         const cs_coord_t   extents[],
                         size_t             n_coords,
                         const cs_coord_t   coords[],
                         fvm_morton_code_t  m_code[]);

/*----------------------------------------------------------------------------
 * Given a Morton code in the grid, compute the Morton codes of its
 * children when refining the grid by one level.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   parent   <-- Morton code associated with parent
 *   children --> array of children Morton codes
 *                (size: 8 in 3D, 4 in 2D, 2 in 1D)
 *----------------------------------------------------------------------------*/

void
fvm_morton_get_children(int                dim,
                        fvm_morton_code_t  parent,
                        fvm_morton_code_t  children[]);

/*----------------------------------------------------------------------------
 * Compare two Morton encoding and check if these two codes are equal,
 * different or shared the same anchor.
 *
 * parameters:
 *   dim    <-- 2D or 3D
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  a type on the kind of relation between the two Morton encodings.
 *----------------------------------------------------------------------------*/

fvm_morton_compare_t
fvm_morton_compare(int                dim,
                   fvm_morton_code_t  code_a,
                   fvm_morton_code_t  code_b);

/*----------------------------------------------------------------------------
 * Locally order a list of Morton ids.
 *
 * parameters:
 *   n_codes      <-- number of Morton ids to order
 *   morton_codes <-- array of Morton ids to order
 *   order        --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
fvm_morton_local_order(cs_lnum_t                n_codes,
                       const fvm_morton_code_t  morton_codes[],
                       cs_lnum_t                order[]);

/*----------------------------------------------------------------------------
 * Locally sort a list of Morton ids.
 *
 * parameters:
 *   n_codes      <-- number of Morton ids to order
 *   morton_codes <-> array of Morton ids to sort
 *----------------------------------------------------------------------------*/

void
fvm_morton_local_sort(cs_lnum_t          n_codes,
                      fvm_morton_code_t  morton_codes[]);

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater than Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

bool
fvm_morton_a_gt_b(fvm_morton_code_t  a,
                  fvm_morton_code_t  b);

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater or equal to Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

bool
fvm_morton_a_ge_b(fvm_morton_code_t  a,
                  fvm_morton_code_t  b);

/*----------------------------------------------------------------------------
 * Get the index associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the array.
 *
 * parameters:
 *   size  <-- size of the array
 *   code  <-- code we are searching for
 *   codes <-- array of Morton codes
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

int
fvm_morton_binary_search(cs_lnum_t           size,
                         fvm_morton_code_t   code,
                         fvm_morton_code_t  *codes);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build a global Morton encoding rank index.
 *
 * The rank_index[i] contains the first Morton code assigned to rank [i].
 *
 * parameters:
 *   dim         <-- 1D, 2D or 3D
 *   gmax_level  <-- level in octree used to build the Morton encoding
 *   n_codes     <-- number of Morton codes to be indexed
 *   morton_code <-- array of Morton codes to be indexed
 *   weight      <-- weighting related to each code
 *   order       <-- ordering array
 *   rank_index  <-> pointer to the global Morton encoding rank index
 *   comm        <-- MPI communicator on which we build the global index
 *
 * returns:
 *  the fit related to the Morton encoding distribution (lower is better).
 *----------------------------------------------------------------------------*/

double
fvm_morton_build_rank_index(int                      dim,
                            int                      gmax_level,
                            cs_gnum_t                n_codes,
                            const fvm_morton_code_t  code[],
                            const cs_lnum_t          weight[],
                            const cs_lnum_t          order[],
                            fvm_morton_code_t        rank_index[],
                            MPI_Comm                 comm);

#endif /* if HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for conversion of a double precision value in
 * range [0, 1] to a given Morton code.
 *
 * \param[in]   s      coordinate between 0 and 1
 * \param[out]  elt    pointer to element
 * \param[in]   input  pointer to optional (untyped) value or structure;
 *                     here, this is an interger representing the spatial
 *                     dimension.
 */
/*----------------------------------------------------------------------------*/

void
fvm_morton_s_to_code(double       s,
                     void        *elt,
                     const void  *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for comparison of 2 Morton codes.
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
fvm_morton_compare_o(const void  *elt1,
                     const void  *elt2,
                     const void  *input);

/*----------------------------------------------------------------------------
 * Dump a Morton to standard output or to a file.
 *
 * parameters:
 *   dim  <-- 2D or 3D
 *   code <-- Morton code to dump
 *----------------------------------------------------------------------------*/

void
fvm_morton_dump(int                 dim,
                fvm_morton_code_t   code);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_MORTON_H__ */
