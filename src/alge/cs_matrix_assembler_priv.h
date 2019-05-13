#ifndef __CS_MATRIX_ASSEMBLER_PRIV_H__
#define __CS_MATRIX_ASSEMBLER_PRIV_H__

/*============================================================================
 * Incremental or general construction of matrix.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_rank_neighbors.h"
#include "cs_timer.h"

#include "cs_matrix.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure used to pre-build a matrix
 *----------------------------------------------------------------------------*/

struct _cs_matrix_assembler_t {

  bool        separate_diag;     /* is diagonal handled separately ? */

  int         flags;             /* sum (bitwise or) of option constants */

  cs_gnum_t   l_range[2];        /* local global row range */
  cs_gnum_t   n_g_rows;          /* global number of rows */
  cs_lnum_t   n_rows;            /* local number of rows */

  cs_lnum_t   size;              /* current insertion array size */
  cs_lnum_t   max_size;          /* maximum insertion array size */

  const cs_lnum_t  *r_idx;       /* shared main row index (0 to n-1) */
  const cs_lnum_t  *c_id;        /* shared main column ids (0 to n-1) */

  cs_lnum_t   *_r_idx;           /* private main row index (0 to n-1) */
  cs_lnum_t   *_c_id;            /* private main column ids (0 to n-1) */

  cs_lnum_t   *d_r_idx;          /* distant row index (0 to n-1) */
  cs_gnum_t   *d_g_c_id;         /* distant global column ids (0 to n-1) */

  cs_gnum_t   *g_rc_id;          /* global row and column ids
                                    (local and distant) */

#if defined(HAVE_MPI)

  /* Distant columns associated with local rows */

  /* Metadata for exchange of matrix coefficient values with other ranks */

  int           n_coeff_ranks;           /* number of MPI ranks with which
                                            coefficients are exchanged */
  int          *coeff_rank;              /* ranks with which coefficients are
                                            exchanged */

  cs_lnum_t     coeff_send_size;         /* number of coefficients to send */
  cs_lnum_t     coeff_recv_size;         /* number of coefficients to receive */

  cs_lnum_t     coeff_send_n_rows;       /* number of matching rows */
  cs_lnum_t    *coeff_send_index;        /* index of sent coefficient rows */
  cs_gnum_t    *coeff_send_row_g_id;     /* global ids matching rows (ordered) */
  cs_gnum_t    *coeff_send_col_g_id;     /* global ids matching columns
                                            (ordered) */

  cs_lnum_t    *coeff_rank_send_index;   /* index of data to send */
  cs_lnum_t    *coeff_rank_recv_index;   /* index of data to receive */

  cs_lnum_t    *coeff_recv_row_id;       /* local row ids associated with
                                            received data; */
  cs_lnum_t    *coeff_recv_col_idx;      /* local column index associated with
                                            received data; the column numbering
                                            implicitely assumes local terms
                                            first, distant ones second */
  cs_gnum_t    *coeff_recv_col_g_id;     /* global column id couples
                                            associated with received data */

  /* Associated communicator */

  MPI_Comm     comm;             /* associated MPI communicator */

  /* Statistics */

  int          n_ranks_init[2];  /* Number of ranks for initial exchange
                                    for distant rows then columns */
#endif /* HAVE_MPI */

  /* Associated vector ghost element info */

  const cs_halo_t  *halo;                /* shared halo for associated vectors */
  cs_halo_t        *_halo;               /* private halo for associated vectors */

  cs_lnum_t         n_e_g_ids;           /* number of external global ids */
  cs_gnum_t        *e_g_id;              /* global ids associated with halo
                                            elements (size: n_e_g_ids */

};

/*----------------------------------------------------------------------------
 * Structure managing matrix coefficient contributions.
 *----------------------------------------------------------------------------*/

struct _cs_matrix_assembler_values_t {

  const  cs_matrix_assembler_t  *ma;  /* associated matrix assembler */

  bool        separate_diag;          /* is diagonal handled separately ? */
  bool        final_assembly;         /* are we ready for final assembly ? */

  cs_lnum_t   db_size[4];             /* Diag Block size including padding:
                                         0: useful block size
                                         1: vector block extents
                                         2: matrix line extents
                                         3: matrix line*column extents */

  cs_lnum_t   eb_size[4];             /* Extradiag block size including padding:
                                         0: useful block size
                                         1: vector block extents
                                         2: matrix line extents
                                         3: matrix line*column extents */

  cs_lnum_t  *diag_idx;               /* Local index of diagonal in each row
                                         when conversion beween separate
                                         diagonal and included diagonal is
                                         required */

  /* Accumulated contributions to distant rows, indexed as per
     coeff_send_index of the matching assembler structure */

#if defined(HAVE_MPI)

  cs_real_t  *coeff_send;

#endif

  /* Matching strcuture and function pointers; some function type may not be
     useful for certain matrix structures or libraries. */

  void                                 *matrix;          /* pointer to
                                                            matrix structure */

  cs_matrix_assembler_values_init_t    *init;
  cs_matrix_assembler_values_add_t     *add_values;
  cs_matrix_assembler_values_add_g_t   *add_values_g;
  cs_matrix_assembler_values_begin_t   *assembly_begin;  /* optional */
  cs_matrix_assembler_values_end_t     *assembly_end;    /* optional */

};

/*============================================================================
 * Public inline function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given local id in a given array of
 *        ordered local ids, when the id might not be present
 *
 * We assume the id is present in the array.
 *
 * \param[in]  l_id_array size  array_size
 * \param[in]  l_id             local id to search for
 * \param[in]  l_id_array       ordered unique local ids array
 *
 * \return  index of l_id in l_id_array, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_l_id_binary_search(cs_lnum_t        l_id_array_size,
                    cs_lnum_t        l_id,
                    const cs_lnum_t  l_id_array[])
{
  if (l_id_array_size < 1)
    return -1;

  cs_lnum_t start_id = 0;
  cs_lnum_t end_id = l_id_array_size - 1;
  cs_lnum_t mid_id = end_id/2;
  while (start_id < end_id) {
    if (l_id_array[mid_id] < l_id)
      start_id = mid_id + 1;
    else if (l_id_array[mid_id] > l_id)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + (end_id - start_id)/2;
  }
  if (l_id_array[mid_id] != l_id)
    mid_id = -1;

  return mid_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given global id in a given array of
 *        ordered global ids, when we know the id is present.
 *
 * We assume the id is present in the array.
 *
 * \param[in]  g_id_array size  array_size
 * \param[in]  g_id             global id to search for
 * \param[in]  g_id_array       ordered unique global ids array
 *
 * \return  index of g_id in g_id_array.
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_g_id_binary_find(cs_lnum_t        g_id_array_size,
                  cs_gnum_t        g_id,
                  const cs_gnum_t  g_id_array[])
{
  cs_lnum_t start_id = 0;
  cs_lnum_t end_id = g_id_array_size - 1;
  cs_lnum_t mid_id = (end_id -start_id) / 2;
  while (start_id < end_id) {
    if (g_id_array[mid_id] < g_id)
      start_id = mid_id + 1;
    else if (g_id_array[mid_id] > g_id)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }
  assert(g_id_array[mid_id] == g_id);

  return mid_id;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_ASSEMBLER_PRIV_H__ */
