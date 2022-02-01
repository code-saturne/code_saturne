#ifndef __CS_MATRIX_UTIL_H__
#define __CS_MATRIX_UTIL_H__

/*============================================================================
 * Utilitary functions for sparse matrixes.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_numbering.h"

#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute diagonal dominance metric.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   dd     --> diagonal dominance (normalized)
 *----------------------------------------------------------------------------*/

void
cs_matrix_diag_dominance(const cs_matrix_t  *matrix,
                         cs_real_t           dd[]);

/*----------------------------------------------------------------------------
 * Dump a linear system matrix and right-hand side to file.
 *
 * parameters:
 *   matrix <-- pointer to mat$rix structure
 *   rhs    <-- right hand side vector
 *   name   <-- identifier string used in file name
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump_linear_system(const cs_matrix_t  *matrix,
                             const cs_real_t     rhs[],
                             const char         *name);

/*----------------------------------------------------------------------------
 * Dump a vector to file.
 *
 * parameters:
 *   n_rows <-- number of local rows
 *   stride <-- associated stride
 *   x      <-- associated vector
 *   name   <-- identifier string used in file name
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump_vector(const cs_lnum_t     n_rows,
                      const cs_lnum_t     stride,
                      const cs_real_t     x[],
                      const char         *name);

/*----------------------------------------------------------------------------
 * Dump a matrix to file.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   name   <-- identifier string used in file name
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump(const cs_matrix_t  *matrix,
               const char         *name);

/*----------------------------------------------------------------------------
 * Log general info relative to matrix.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_matrix_log_info(const cs_matrix_t  *matrix,
                   int                 verbosity);

/*----------------------------------------------------------------------------
 * Test matrix dump operations.
 *
 * parameters:
 *   n_rows     <-- number of local rows
 *   n_cols_ext <-- number of colmuns including ghost columns (array size)
 *   n_edges    <-- local number of graph edges
 *   edges      <-- graph edges connectivity
 *   halo       <-- cell halo structure
 *   numbering  <-- vectorization or thread-related numbering info, or NULL
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump_test(cs_lnum_t              n_rows,
                    cs_lnum_t              n_cols_ext,
                    cs_lnum_t              n_edges,
                    const cs_lnum_2_t     *edges,
                    const cs_halo_t       *halo,
                    const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_UTIL_H__ */
