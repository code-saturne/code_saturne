#ifndef __CS_MATRIX_ASSEMBLER_H__
#define __CS_MATRIX_ASSEMBLER_H__

/*============================================================================
 * Incremental or general construction of matrix.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * Matrix assembler option flags
 */

#define CS_MATRIX_DISTANT_ROW_USE_COL_IDX     (1 << 0)
#define CS_MATRIX_DISTANT_ROW_USE_COL_G_ID    (1 << 1)

#define CS_MATRIX_EXTERNAL_HALO               (1 << 2)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Structure used to pre-build a matrix structure */

typedef struct _cs_matrix_assembler_t  cs_matrix_assembler_t;

/*! Structure used to set matrix coefficients */

typedef struct _cs_matrix_assembler_values_t  cs_matrix_assembler_values_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for initialization of matrix coefficients using
 *        local row ids and column indexes.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix   untyped pointer to matrix description structure
 * \param[in]       db size  optional diagonal block sizes
 * \param[in]       eb size  optional extra-diagonal block sizes
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_matrix_assembler_values_init_t) (void        *matrix,
                                     const int   *db_size,
                                     const int   *eb_size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for addition to matrix coefficients using
 *        local row ids and column indexes.
 *
 * Values whose associated row index is negative should be ignored;
 * Values whose column index is -1 are assumed to be assigned to a
 * separately stored diagonal. Other indexes should be valid.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in, out]  matrix   untyped pointer to matrix description structure
 * \param[in]       n        number of values to add
 * \param[in]       stride   associated data block size
 * \param[in]       row_id   associated local row ids
 * \param[in]       col_idx  associated local column indexes
 * \param[in]       val      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_matrix_assembler_values_add_t) (void             *matrix,
                                    cs_lnum_t         n,
                                    cs_lnum_t         stride,
                                    const cs_lnum_t   row_id[],
                                    const cs_lnum_t   col_idx[],
                                    const cs_real_t   vals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for addition to matrix coefficients using
 *        global row ids and column indexes.
 *
 * Values whose global row id is outside the local range should be ignored
 * (this allows simply filtering values before calling this function,
 * avoiding an extra copy (a copy will probably be done by the function
 * anyways).
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix    untyped pointer to matrix description structure
 * \param[in]       n         number of values to add
 * \param[in]       stride    associated data block size
 * \param[in]       g_row_id  associated global row ids
 * \param[in]       g_col_id  associated global column ids
 * \param[in]       val       pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_matrix_assembler_values_add_g_t) (void            *matrix,
                                      cs_lnum_t        n,
                                      cs_lnum_t        stride,
                                      const cs_gnum_t  row_g_id[],
                                      const cs_gnum_t  col_g_id[],
                                      const cs_real_t  vals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer to start the final assembly of matrix coefficients.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix  untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_matrix_assembler_values_begin_t) (void  *matrix);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer to complete the final assembly of matrix
 *        coefficients.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix  untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_matrix_assembler_values_end_t) (void  *matrix);

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix assembler structure.
 *
 * The associated matrix structure data will initially be empty, though
 * the range of rows associated with the current rank must be defined
 * immediately.
 *
 * \param[in]  l_range        global id range [min, max[ for local rank
 * \param[in]  separate_diag  if true, diagonal terms are handled separately
 *
 * \return  pointer to created matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_t *
cs_matrix_assembler_create(const cs_gnum_t  l_range[2],
                           bool             separate_diag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix assembler structure based on a given connectivity
 *        and associated halo structure.
 *
 * The assembler may not be later modified when built this way.
 *
 * The connectivity arrays and halo will be shared by the caller, so their
 * life cycle must be at least as long as the matrix assembler structure.
 *
 * \param[in]  n_rows         number fo local rows
 * \param[in]  separate_diag  if true, diagonal terms are handled separately
 * \param[in]  row_idx        matrix row index
 * \param[in]  col_id         matrix column indexes
 * \param[in]  halo           associated halo structure
 *
 * \return  pointer to created matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_t *
cs_matrix_assembler_create_from_shared(cs_lnum_t         n_rows,
                                       bool              separate_diag,
                                       const cs_lnum_t   row_idx[],
                                       const cs_lnum_t   col_id[],
                                       const cs_halo_t  *halo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a matrix assembler structure.
 *
 * \param[in, out]  ma  pointer to matrix assembler structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_destroy(cs_matrix_assembler_t  **ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add entries to a matrix assembler structure.
 *
 * This function should be called by a single thread for a given assembler.
 *
 * \param[in, out]  ma        pointer to matrix assembler structure
 * \param[in]       n         number of entries
 * \param[in]       g_col_id  global column ids associated with entries
 * \param[in]       g_row_id  global row ids associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_add_g_ids(cs_matrix_assembler_t  *ma,
                              cs_lnum_t               n,
                              cs_gnum_t               g_row_id[],
                              cs_gnum_t               g_col_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute internal structures required by a matrix assembler.
 *
 * The associated vector halo is also computed at this stage.
 *
 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_compute(cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set option flags for a given matrix assembler structure.
 *
 * When used, this function should be called before defining entries
 * (in case of future option additions) and in any case before computing
 * the final structure.
 *
 * Flags are defined as a sum (bitwise or) of constants, described below:
 * - \ref CS_MATRIX_DISTANT_ROW_USE_COL_IDX indicates that column indexes
 *   matching distant row data are computed and maintained, so the
 *   associated coefficient contributions can be added more efficiently.
 * - \ref CS_MATRIX_DISTANT_ROW_USE_COL_G_ID  indicates that columns
 *   global ids matching distant row data are maintained directly,
 *   so no lookup to a global id description of matrix columns is needed;
 *   this option is useful only when using an external library allowing
 *   incremental setting or accumulation of coefficients using global
 *   column ids, such as with PETSc's MatSetValues. When building matrices
 *   locally first (which includes most cases, whether internal matrices,
 *   or using library conversion or import functions such as PETSc's
 *   MatCreateMPIAIJWithArrays) this is slightly less efficient, as
 *   column ids will need to be matched to each row's columns a second
 *   time for those (exchanged) coefficients.
 * - \ref CS_MATRIX_EXTERNAL_HALO indicates that we do not need to
 *   build an associated halo structure, as it will be built externally
 *   (i.e. by an external library such as PETSc, HYPRE, ...) using its
 *   own equivalent structures.
 *
 * \param[in, out]  ma     pointer to matrix assembler structure
 * \param[in]       flags  sum of matrix assembler flag constants
 *                         (bitwise or)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_set_options(cs_matrix_assembler_t  *ma,
                                int                     flags);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the option flags for a given matrix assembler structure.
 *
 * Flags are defined as a sum (bitwise or) of constants, described in
 * \ref cs_matrix_assembler_set_options.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  option flags (sum of integer constants) used y this structure
 */
/*----------------------------------------------------------------------------*/

int
cs_matrix_assembler_get_options(const cs_matrix_assembler_t  *ma);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign specific MPI communicator to matrix assembler.
 *
 * This must be called before \ref cs_matrix_assembler_compute.
 *
 * \param[in, out]  ma    pointer to matrix assembler structure
 * \param[in]       comm  assigned MPI communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_set_comm(cs_matrix_assembler_t  *ma,
                             MPI_Comm                comm);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the halo structure associated with a
 *        matrix assembler.
 *
 * \param[in]  ma   pointer to matrix assembler structure
 *
 * \return  pointer to halo structure
 */
/*----------------------------------------------------------------------------*/

const cs_halo_t *
cs_matrix_assembler_get_halo(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if the matrix assembler is based on a separate diagonal.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  true if the associated structure has a separate diagonal,
 *          false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_matrix_assembler_get_separate_diag(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of rows associated with a matrix assembler.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  number of rows associated with matrix assembler
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_assembler_get_n_rows(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of columns associated with a matrix assembler.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  number of columns associated with matrix assembler
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_assembler_get_n_columns(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a row index associated with a matrix assembler.
 *
 * This index is a CSR type index relative to all columns, including or
 * excluding the diagonal depending on the \c separate_diag option used
 * when creating the matrix assembler. The matching column ids can be
 * obtained using \ref cs_matrix_assembler_get_col_ids.
 *
 * \warning the returned index is valid only as long as the matrix assembly
 * structure exists.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  pointer to matrix structure row index
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_matrix_assembler_get_row_index(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a column ids associated with a matrix assembler.
 *
 * These ids relative to all columns of successive rows, with columns
 * ordered by local id, so local columns appear first, distant columns
 * after. Depending on the \c separate_diag option used when creating the
 * matrix assembler, the columns may include the diagonal or not.
 * The matching index can be obtained using
 * \ref cs_matrix_assembler_get_row_index.
 *
 * \warning the returned index is valid only as long as the matrix assembly
 * structure exists.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  pointer to matrix structure row index
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_matrix_assembler_get_col_ids(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return info on the number of neighbor ranks a matrix assembler
 *        may communicate with.
 *
 * \param[in]   ma   pointer to matrix assembler structure
 * \param[out]  rc   rank counts; the 4 values are:
 *                   - 0 number of communicating ranks for vector update (halo)
 *                   - 1 number of communicating ranks for matrix assembly
 *                   - 2 maximum number of communicating ranks for halo
 *                       construction (assumed ranks algorithm)
 *                   - 3 maximum number of communicating ranks for matrix
 *                       assembly determination (assumed ranks algorithm)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_get_rank_counts(const cs_matrix_assembler_t  *ma,
                                    int                           rc[4]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log rank counts for a given matrix assembler.
 *
 * \param[in]  ma      pointer to matrix assembler structure
 * \param[in]  log_id  log type
 * \param[in]  name    name of this assembler
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_log_rank_counts(const cs_matrix_assembler_t  *ma,
                                    cs_log_t                      log_id,
                                    const char                   *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a matrix assembler values structure.
 *
 * The associated values will initially be set to zero.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * This is a low-level function, which should be called by a simpler
 * function (\ref cs_matrix_assembler_values_init) which provides
 * the necessary function pointers.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure. In a similar manner, the life cycle of the associated
 *           matrix assembler must also be at least as long as that
 *           of the assembler values structure.
 *
 * \param[in]       ma        associated matrix assembly structure
 * \param[in]       sep_diag  true if diagonal terms are stored separately
 * \param[in]       db_size   optional diagonal block sizes
 * \param[in]       eb_size   optional extra-diagonal block sizes
 * \param[in, out]  matrix    untyped pointer to matrix description structure
 * \param[in]       init      pointer to matrix coefficients
 *                            initialization function
 * \param[in]       add       pointer to matrix coefficients addition from
 *                            local ids function
 * \param[in]       add_g     pointer to matrix coefficients addition from
 *                            global ids function
 * \param[in]       begin     pointer to matrix coefficients assembly
 *                            start function (optional)
 * \param[in]       end       pointer to matrix coefficients assembly
 *                            end function (optional)
 *
 * \return  pointer to initialized matrix assembler structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_values_t *
cs_matrix_assembler_values_create(const cs_matrix_assembler_t          *ma,
                                  bool                                  sep_diag,
                                  const int                            *db_size,
                                  const int                            *eb_size,
                                  void                                 *matrix,
                                  cs_matrix_assembler_values_init_t    *init,
                                  cs_matrix_assembler_values_add_t     *add,
                                  cs_matrix_assembler_values_add_g_t   *add_g,
                                  cs_matrix_assembler_values_begin_t   *begin,
                                  cs_matrix_assembler_values_end_t     *end);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize matrix assembler values structure.
 *
 * When this function returns, the assembler values structure has been
 * destroyed, and the associated matrix is fully assembled (i.e. ready to
 * use).
 *
 * \param[in, out]  mav  pointer to matrix assembler values structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_finalize(cs_matrix_assembler_values_t  **mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using local
 *        row and column ids.
 *
 * If the matching matrix coefficients use a block structure, the
 * values passed to this function should use the same block size.
 *
 * Note also that if the matrix has different diagonal and extradiagonal
 * block sizes, separate calls to this function should be used for both
 * types of coefficients; in addition, in this case, diagonal coefficients
 * should only be provided by the owning rank (this also impacts how
 * the associated matrix assembler structure is defined).
 *
 * This function may be called by different threads, as long those threads
 * do not add contributions to the same rows (assuming caution is taken
 * in the case of external libraries so that their builder functions
 * have tht same property).
 *
 * \param[in, out]  mav     pointer to matrix assembler values structure
 * \param[in]       n       number of entries
 * \param[in]       col_id  local column ids associated with entries
 * \param[in]       row_id  local row ids associated with entries
 * \param[in]       val     values associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_add(cs_matrix_assembler_values_t  *mav,
                               cs_lnum_t                      n,
                               cs_lnum_t                      row_id[],
                               cs_lnum_t                      col_id[],
                               cs_real_t                      val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *
 * If the matching matrix coefficients use a block structure, the
 * values passed to this function should use the same block size.
 *
 * Note also that if the matrix has different diagonal and extradiagonal
 * block sizes, separate calls to this function should be used for both
 * types of coefficients; in addition, in this case, diagonal coefficients
 * should only be provided by the owning rank (this also impacts how
 * the associated matrix assembler structure is defined).
 *
 * This function may be called by different threads, as long those threads
 * do not add contributions to the same rows (assuming caution is taken
 * in the case of external libraries so that their builder functions
 * have tht same property).
 *
 * \param[in, out]  mav       pointer to matrix assembler values structure
 * \param[in]       n         number of entries
 * \param[in]       g_col_id  global column ids associated with entries
 * \param[in]       g_row_id  global row ids associated with entries
 * \param[in]       val       values associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_add_g(cs_matrix_assembler_values_t  *mav,
                                 cs_lnum_t                      n,
                                 cs_gnum_t                      g_row_id[],
                                 cs_gnum_t                      g_col_id[],
                                 cs_real_t                      val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Start assembly of matrix values structure.
 *
 * The call to this function is always optional, and indicates that assembly
 * may begin. Calling it before \ref cs_matrix_assembler_values_finalize
 * may be useful when the underlying libraries can overlap communication
 * (or other operations such as matrix trasformation) and computation.
 *
 * \remark When values have been added to rows belonging to another rank,
 *         communication will occur here. Splitting this function could
 *         add another opportunity for overlap of computation and
 *         communication, but is not done as of 2016, as it would add
 *         complexity, with probably limited returns, as the effective
 *         asynchonous progress of MPI libraries is usually limited.
 *
 * \param[in, out]  mav  pointer to matrix assembler data structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_done(cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_ASSEMBLER_H__ */
