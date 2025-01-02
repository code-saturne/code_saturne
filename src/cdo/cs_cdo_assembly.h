#ifndef __CS_CDO_ASSEMBLY_H__
#define __CS_CDO_ASSEMBLY_H__

/*============================================================================
 * Set of functions and structures to handle the assembly of cellwise local CDO
 * systems into a cs_matrix_t structure through the cs_matrix_assembler_t and
 * its related structures
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_matrix.h"
#include "alge/cs_matrix_assembler.h"
#include "base/cs_param_types.h"
#include "base/cs_range_set.h"
#include "cdo/cs_sdm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Size of the buffer used to collect global ids for rows and columns when
 * assembling the values in the global matrix from the local cellwise matrices
 */

#define CS_CDO_ASSEMBLY_BUFSIZE  200

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

typedef struct _cs_cdo_assembly_t  cs_cdo_assembly_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Block or no block versions are handled
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_assembly_func_t)(const cs_sdm_t                    *m,
                         const cs_lnum_t                   *dof_ids,
                         const cs_range_set_t              *rset,
                         cs_cdo_assembly_t                 *eqa,
                         cs_matrix_assembler_values_t      *mav);

/* Metadata for the assembly of a row */

typedef struct {

  /* Row numberings */

  cs_gnum_t   g_id;             /* Global row numbering */
  cs_lnum_t   l_id;             /* Local range set numbering */
  int         i;                /* Cellwise system numbering */

  /* Column-related members */

  int                 n_cols;    /* Number of columns (cellwise system) */
  cs_gnum_t          *col_g_id;  /* Global numbering of columns */
  int                *col_idx;   /* Array to build and to give as parameter
                                    for *add_vals() functions */

  const cs_real_t    *val;       /* Row values */
  cs_real_t          *expval;    /* Expanded row values (when unrolling non
                                    scalar-valued blocks) */

} cs_cdo_assembly_row_t;

/* Cell-wise structure used for the assembly of local systems */

struct _cs_cdo_assembly_t {

  int      n_cw_dofs;   /* Number of DoFs in a cell */
  int      ddim;        /* Number of real values related to each diagonal
                           entry */
  int      edim;        /* Number of real values related to each
                           extra-diagonal entry */

  /* For matrices built by blocks (for instance in coupled systems with
     scalar-valued blocks), one may need to shift the row and/or the column
     local ids */

  cs_lnum_t   l_col_shift;
  cs_lnum_t   l_row_shift;

  cs_cdo_assembly_row_t    *row;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information for the setup.log file
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_setup_log(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a pointer to a cs_cdo_assembly_t structure related to a given
 *        thread
 *
 * \param[in] t_id    id in the array of pointer
 *
 * \return a pointer to a cs_cdo_assembly_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_assembly_t *
cs_cdo_assembly_get(int    t_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate cs_cdo_assembly_t structure (shared among schemes). Each
 *        thread has its own copy this structure to enable a multithreaded
 *        assembly process.
 *
 * \param[in] ddim          max number of dof values on the diagonal part
 * \param[in] edim          max number of dof values on the extra-diag. part
 * \param[in] n_cw_dofs     max number of DoFs in a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_init(int     ddim,
                     int     edim,
                     int     n_cw_dofs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free matrix-related structures used during the simulation.
 *        Display overall statistic about the assembly stage for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the current shift values to consider during the assembly stage
 *
 * \param[in, out] asb          pointer to a cs_cdo_assembly_t to update
 * \param[in]      l_row_shift  shift to apply to local row ids
 * \param[in]      l_col_shift  shift to apply to local col ids
 *
 * \return a function pointer cs_cdo_assembly_func_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_set_shift(cs_cdo_assembly_t    *asb,
                          cs_lnum_t             l_row_shift,
                          cs_lnum_t             l_col_shift);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Rely on the generic cs_matrix_assembler_values_add_g() function.
 *        Case of scalar-valued matrices.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_scal_generic(const cs_sdm_t                   *m,
                                    const cs_lnum_t                  *dof_ids,
                                    const cs_range_set_t             *rset,
                                    cs_cdo_assembly_t                *asb,
                                    cs_matrix_assembler_values_t     *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Rely on the generic cs_matrix_assembler_values_add_g() function.
 *        Case of vector-valued matrices with an expanded 3x3 block.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_e33_generic(const cs_sdm_t                  *m,
                                   const cs_lnum_t                 *dof_ids,
                                   const cs_range_set_t            *rset,
                                   cs_cdo_assembly_t               *asb,
                                   cs_matrix_assembler_values_t    *mav);

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Scalar-valued case. Parallel and with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_mpit(const cs_sdm_t                   *m,
                            const cs_lnum_t                  *dof_ids,
                            const cs_range_set_t             *rset,
                            cs_cdo_assembly_t                *asb,
                            cs_matrix_assembler_values_t     *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Scalar-valued case. Parallel without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_mpis(const cs_sdm_t                   *m,
                            const cs_lnum_t                  *dof_ids,
                            const cs_range_set_t             *rset,
                            cs_cdo_assembly_t                *asb,
                            cs_matrix_assembler_values_t     *mav);
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Scalar-valued case. Sequential and with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_seqt(const cs_sdm_t                  *m,
                            const cs_lnum_t                 *dof_ids,
                            const cs_range_set_t            *rset,
                            cs_cdo_assembly_t               *asb,
                            cs_matrix_assembler_values_t    *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix
 *        Scalar-valued case.
 *        Sequential and without openMP (single thread).
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_seqs(const cs_sdm_t                  *m,
                            const cs_lnum_t                 *dof_ids,
                            const cs_range_set_t            *rset,
                            cs_cdo_assembly_t               *asb,
                            cs_matrix_assembler_values_t    *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise (no-block) matrix into the global matrix
 *        corresponding to a system of coupled equations.
 *        Scalar-valued case. Sequential and without openMP.
 *        Block matrices assembled from cellwise scalar-valued matrices
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_seqs(const cs_sdm_t                  *m,
                                const cs_lnum_t                 *dof_ids,
                                const cs_range_set_t            *rset,
                                cs_cdo_assembly_t               *asb,
                                cs_matrix_assembler_values_t    *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise (no-block) matrix into the global matrix
 *        corresponding to a system of coupled equations.
 *        Scalar-valued case. Sequential and with openMP.
 *        Block matrices assembled from cellwise scalar-valued matrices
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_seqt(const cs_sdm_t                  *m,
                                const cs_lnum_t                 *dof_ids,
                                const cs_range_set_t            *rset,
                                cs_cdo_assembly_t               *asb,
                                cs_matrix_assembler_values_t    *mav);

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global coupled matrix.
 *        Scalar-valued case. Parallel without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_mpis(const cs_sdm_t                   *m,
                                const cs_lnum_t                  *dof_ids,
                                const cs_range_set_t             *rset,
                                cs_cdo_assembly_t                *asb,
                                cs_matrix_assembler_values_t     *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global coupled matrix.
 *        Scalar-valued case. Parallel with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_mpit(const cs_sdm_t                   *m,
                                const cs_lnum_t                  *dof_ids,
                                const cs_range_set_t             *rset,
                                cs_cdo_assembly_t                *asb,
                                cs_matrix_assembler_values_t     *mav);
#endif  /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_seqs(const cs_sdm_t               *m,
                                     const cs_lnum_t              *dof_ids,
                                     const cs_range_set_t         *rset,
                                     cs_cdo_assembly_t            *asb,
                                     cs_matrix_assembler_values_t *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Case of a block 3x3 entries. Expand each row.
 *        Sequential run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_seqt(const cs_sdm_t                *m,
                                     const cs_lnum_t               *dof_ids,
                                     const cs_range_set_t          *rset,
                                     cs_cdo_assembly_t             *asb,
                                     cs_matrix_assembler_values_t  *mav);

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Case of a block 3x3 entries. Expand each row.
 *        Parallel run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_mpis(const cs_sdm_t                *m,
                                     const cs_lnum_t               *dof_ids,
                                     const cs_range_set_t          *rset,
                                     cs_cdo_assembly_t             *asb,
                                     cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_mpit(const cs_sdm_t               *m,
                                     const cs_lnum_t              *dof_ids,
                                     const cs_range_set_t         *rset,
                                     cs_cdo_assembly_t            *asb,
                                     cs_matrix_assembler_values_t *mav);
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Case of a block 3x3 entries. Expand each row.
 *        Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_seqs(const cs_sdm_t                *m,
                                    const cs_lnum_t               *dof_ids,
                                    const cs_range_set_t          *rset,
                                    cs_cdo_assembly_t             *asb,
                                    cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Case of a block 3x3 entries. Expand each row.
 *        Sequential run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_seqt(const cs_sdm_t                *m,
                                    const cs_lnum_t               *dof_ids,
                                    const cs_range_set_t          *rset,
                                    cs_cdo_assembly_t             *asb,
                                    cs_matrix_assembler_values_t  *mav);

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Case of a block 3x3 entries. Expand each row.
 *        Parallel run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_mpis(const cs_sdm_t                *m,
                                    const cs_lnum_t               *dof_ids,
                                    const cs_range_set_t          *rset,
                                    cs_cdo_assembly_t             *asb,
                                    cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assemble a cellwise matrix into the global matrix.
 *        Case of a block 3x3 entries. Expand each row.
 *        Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_mpit(const cs_sdm_t                *m,
                                    const cs_lnum_t               *dof_ids,
                                    const cs_range_set_t          *rset,
                                    cs_cdo_assembly_t             *asb,
                                    cs_matrix_assembler_values_t  *mav);
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_seqs(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_seqt(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav);

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_mpis(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_mpit(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav);
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_ASSEMBLY_H__ */
