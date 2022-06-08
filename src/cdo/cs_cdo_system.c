/*============================================================================
 * Structure and functions used to manipulate elementary structures related to
 * the definition of a matrix, parallel synchronization/operations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_flag.h"
#include "cs_matrix_priv.h"

#if defined(HAVE_HYPRE)
#include "cs_matrix_hypre.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_CDO_SYSTEM_DBG        0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

static cs_mesh_t  *cs_shared_mesh = NULL;
static cs_cdo_connect_t  *cs_shared_connect = NULL;

static int  _n_cdo_block_structures = 0;
static cs_cdo_system_block_t  **_cdo_block_structures = NULL;

/*============================================================================
 * Static inline private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of scalar-valued matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_cdo_assembly_func_t *
_set_scalar_assembly_func(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_matrix_mpis;
    else                       /* With OpenMP */
      return cs_cdo_assembly_matrix_mpit;

  }
#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) { /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_matrix_seqs;
    else                       /* With OpenMP */
      return cs_cdo_assembly_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of scalar-valued matrices used inside a system of equations.
 *         This is system is seen as a slave w.r.t. the master one which has a
 *         view on the full system
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_cdo_assembly_func_t *
_set_scalar_slave_assembly_func(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {  /* Parallel */

    return cs_cdo_assembly_matrix_sys_mpis;
    /* TODO: Threaded version */

  }
#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) { /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_matrix_sys_seqs;
    else                       /* With OpenMP */
      return cs_cdo_assembly_matrix_sys_seqt;

  }

  return NULL; /* Case not handled */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of block 3x3 matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_cdo_assembly_func_t *
_set_eblock33_assembly_func(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_eblock33_matrix_mpis;
    else                      /* With OpenMP */
      return cs_cdo_assembly_eblock33_matrix_mpit;

  }
#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) {  /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_eblock33_matrix_seqs;
    else                      /* With OpenMP */
      return cs_cdo_assembly_eblock33_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of block 3x3 matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_cdo_assembly_func_t *
_set_block33_assembly_func(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_block33_matrix_mpis;
    else                      /* With OpenMP */
      return cs_cdo_assembly_block33_matrix_mpit;

  }
#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) {  /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_block33_matrix_seqs;
    else                      /* With OpenMP */
      return cs_cdo_assembly_block33_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of block NxN matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_cdo_assembly_func_t *
_set_block_assembly_func(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_eblock_matrix_mpis;
    else                      /* With OpenMP */
      return cs_cdo_assembly_eblock_matrix_mpit;

  }
#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) {  /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_cdo_assembly_eblock_matrix_seqs;
    else                      /* With OpenMP */
      return cs_cdo_assembly_eblock_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the assembly function associated to a block according to the
 *        the metadata of the block. Case of an assembly designed by coupled
 *        systems of equations
 *
 * \param[in] bi         block info structure to consider
 *
 * \return a pointer to a function or NULL if useless
 */
/*----------------------------------------------------------------------------*/

static cs_cdo_assembly_func_t *
_assign_slave_assembly_func(const cs_cdo_system_block_info_t   bi)
{
  CS_UNUSED(bi);

  /* Up to now, there is no choice */

  return _set_scalar_slave_assembly_func();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the assembly function associated to a block according to the
 *        metadata of the block
 *
 * \param[in] bi         block info structure to consider
 *
 * \return a pointer to a function or NULL if useless
 */
/*----------------------------------------------------------------------------*/

static cs_cdo_assembly_func_t *
_assign_assembly_func(const cs_cdo_system_block_info_t   bi)
{
  if (bi.stride == 1)  { /* Assemble a cell system with 1 DoF by element */

    if (bi.matrix_class == CS_CDO_SYSTEM_MATRIX_HYPRE)
      return cs_cdo_assembly_matrix_scal_generic;
    else
      return _set_scalar_assembly_func();

  }
  else if (bi.stride == 3) {

    if (bi.matrix_class == CS_CDO_SYSTEM_MATRIX_HYPRE)
      return cs_cdo_assembly_matrix_e33_generic;

    else {

      if (bi.unrolled)
        return _set_eblock33_assembly_func();
      else
        return _set_block33_assembly_func();

    }

  }
  else
    return _set_block_assembly_func();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set (and sometimes build) a cs_range_set_t structure and a
 *        cs_interface_set_t structure for the given block structure.
 *
 * \param[in]      forced    define a new structure even if already allocated
 * \param[in, out] b         block structure to update
 * \param[in]      scal_ifs  reference interface set in the scalar-valued case
 * \param[in]      scal_rset reference range set in the scalar-valued case
 */
/*----------------------------------------------------------------------------*/

static void
_assign_ifs_rset(bool                               forced,
                 cs_cdo_system_block_t             *b,
                 cs_interface_set_t                *scal_ifs,
                 cs_range_set_t                    *scal_rset)
{
  assert(b != NULL);
  assert(scal_ifs != NULL || cs_glob_n_ranks < 2);

  const cs_cdo_system_block_info_t  bi = b->info;

  switch (b->type) {

  case CS_CDO_SYSTEM_BLOCK_DEFAULT:
    /* --------------------------- */
    {
      cs_cdo_system_dblock_t  *db = b->block_pointer;

      if (bi.stride == 1) { /* scalar-valued DoF --> always shared */

        db->interface_set = scal_ifs;
        db->range_set = scal_rset;

      }
      else {

        if (forced && b->owner) { /* One forces the build of new structures */

          cs_interface_set_destroy(&(db->interface_set));
          cs_range_set_destroy(&(db->range_set));

          db->interface_set = NULL;
          db->range_set = NULL;

        }

        if (!bi.unrolled) {

            db->interface_set = scal_ifs;
            db->range_set = scal_rset;

        }

        if (db->interface_set == NULL) {

          if (bi.interlaced)
            db->interface_set = cs_interface_set_dup(scal_ifs, bi.stride);
          else
            db->interface_set = cs_interface_set_dup_blocks(scal_ifs,
                                                            bi.n_elements,
                                                            bi.stride);

        } /* interface set not defined */

        if (db->range_set == NULL) {

          db->range_set = cs_range_set_create(db->interface_set,
                                              NULL,
                                              bi.n_elements*bi.stride,
                                              false,  /* option for balance */
                                              1,      /* tr_ignore */
                                              0);     /* g_id_base */

        }

      } /* stride > 1 */
    }
    break; /* default block */

  case CS_CDO_SYSTEM_BLOCK_SPLIT:
    /* ------------------------- */
    {
      cs_cdo_system_sblock_t  *sb = b->block_pointer;

      /* Always treated as a scalar-valued block (if stride = 3, then the 9
         blocks are scalar-valued and thus, rely on shared interface and range
         sets */

      sb->interface_set = scal_ifs;
      sb->range_set = scal_rset;
    }
    break; /* split block */

  case CS_CDO_SYSTEM_BLOCK_UNASS:
    /* ------------------------- */
    {
      cs_cdo_system_ublock_t  *ub = b->block_pointer;

      if (bi.stride == 1) { /* scalar-valued DoF --> always shared */

        ub->interface_set = scal_ifs;
        ub->range_set = scal_rset;

      }
      else {

        if (forced && b->owner) { /* One forces the build of new structures */

          cs_interface_set_destroy(&(ub->interface_set));
          cs_range_set_destroy(&(ub->range_set));

          ub->interface_set = NULL;
          ub->range_set = NULL;

        }

        if (ub->interface_set == NULL) {

          if (bi.interlaced)
            ub->interface_set = cs_interface_set_dup(scal_ifs, bi.stride);
          else
            ub->interface_set = cs_interface_set_dup_blocks(scal_ifs,
                                                            bi.n_elements,
                                                            bi.stride);

        } /* interface set not defined */

        if (ub->range_set == NULL) {

          ub->range_set = cs_range_set_create(ub->interface_set,
                                              NULL,
                                              bi.n_elements*bi.stride,
                                              false,  /* option for balance */
                                              1,      /* tr_ignore */
                                              0);     /* g_id_base */

        }

      } /* stride > 1 */
    }
    break; /* unassembled block */

  default:
    break; /* do nothing */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and define a new cs_matrix_assembler_t structure.
 *         Case of structures with interlaced DoFs.
 *
 * \param[in]  stride   number of DoFs by "x" entity
 * \param[in]  x2x      pointer to a cs_adjacency_t structure
 * \param[in]  rs       pointer to a range set structure
 *
 * \return a pointer to a new allocated cs_matrix_assembler_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_assembler_t *
_build_interlaced_ma(int                      stride,
                     const cs_adjacency_t    *x2x,
                     const cs_range_set_t    *rs)
{
  cs_gnum_t  *grows = NULL, *gcols = NULL;

  /* The second parameter is set to "true" meaning that the diagonal is stored
     separately. This corresponds to a MSR matrix storage. This is always the
     case in CDO schemes */

  cs_matrix_assembler_t  *ma = cs_matrix_assembler_create(rs->l_range, true);

  /* First loop to count the max. size of the temporary buffers */

  cs_lnum_t  n_x = x2x->n_elts;
  cs_lnum_t  max_size = 0;
  for (cs_lnum_t i = 0; i < n_x; i++)
    max_size = CS_MAX(max_size, x2x->idx[i+1] - x2x->idx[i]);

  /* We increment the max. size to take into account the diagonal entry */

  int  buf_size = stride * stride * (max_size + 1);

  BFT_MALLOC(grows, buf_size, cs_gnum_t);
  BFT_MALLOC(gcols, buf_size, cs_gnum_t);

  if (stride == 1)  { /* Simplified version */

    for (cs_lnum_t row_id = 0; row_id < n_x; row_id++) {

      const cs_gnum_t  grow_id = rs->g_id[row_id];
      const cs_lnum_t  start = x2x->idx[row_id];
      const cs_lnum_t  end = x2x->idx[row_id+1];

      /* Diagonal term is excluded in the x2x connectivity. Add it. */

      grows[0] = grow_id, gcols[0] = grow_id;

      /* Extra diagonal couples */

      for (cs_lnum_t j = start, i = 1; j < end; j++, i++) {
        grows[i] = grow_id;
        gcols[i] = rs->g_id[x2x->ids[j]];
      }

      cs_matrix_assembler_add_g_ids(ma, end - start + 1, grows, gcols);

    } /* Loop on entities */

  }
  else { /* stride > 1 */

    for (cs_lnum_t row_id = 0; row_id < n_x; row_id++) {

      const cs_lnum_t  start = x2x->idx[row_id];
      const cs_lnum_t  end = x2x->idx[row_id+1];
      const int  n_entries = (end - start + 1) * stride * stride;
      const cs_gnum_t  *grow_ids = rs->g_id + row_id*stride;

      int shift = 0;

      /* Diagonal term is excluded in this connectivity. Add it "manually" */

      for (int dof_i = 0; dof_i < stride; dof_i++) {
        const cs_gnum_t  grow_id = grow_ids[dof_i];
        for (int dof_j = 0; dof_j < stride; dof_j++) {
          grows[shift] = grow_id;
          gcols[shift] = grow_ids[dof_j];
          shift++;
        }
      }

      /* Extra diagonal couples */

      for (cs_lnum_t j = start; j < end; j++) {

        const cs_lnum_t  col_id = x2x->ids[j];
        const cs_gnum_t  *gcol_ids = rs->g_id + col_id*stride;

        for (int dof_i = 0; dof_i < stride; dof_i++) {
          const cs_gnum_t  grow_id = grow_ids[dof_i];
          for (int dof_j = 0; dof_j < stride; dof_j++) {
            grows[shift] = grow_id;
            gcols[shift] = gcol_ids[dof_j];
            shift++;
          }
        }

      } /* Loop on number of DoFs by entity */

      assert(shift == n_entries);
      cs_matrix_assembler_add_g_ids(ma, n_entries, grows, gcols);

    } /* Loop on entities */

  }

  /* Now compute structure */

  cs_matrix_assembler_compute(ma);

  /* Free temporary buffers */

  BFT_FREE(grows);
  BFT_FREE(gcols);

  return ma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and define a new cs_matrix_assembler_t structure.
 *         Case of structures with no interlaced DoFs.
 *
 * \param[in]  stride   number of DoFs by "x" entity
 * \param[in]  x2x      pointer to a cs_adjacency_t structure
 * \param[in]  rs       pointer to a range set structure
 *
 * \return a pointer to a new allocated cs_matrix_assembler_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_assembler_t *
_build_no_interlaced_ma(int                      stride,
                        const cs_adjacency_t    *x2x,
                        const cs_range_set_t    *rs)
{
  cs_gnum_t  *grows = NULL, *gcols = NULL, *g_r_ids = NULL, *g_c_ids = NULL;

  /* The second parameter is set to "true" meaning that the diagonal is stored
     separately. This corresponds to a MSR matrix storage. This is always the
     case in CDO schemes */

  cs_matrix_assembler_t  *ma = cs_matrix_assembler_create(rs->l_range, true);

  /* First loop to count the max. size of the temporary buffers */

  cs_lnum_t  n_x = x2x->n_elts;
  cs_lnum_t  max_size = 0;
  for (cs_lnum_t i = 0; i < n_x; i++)
    max_size = CS_MAX(max_size, x2x->idx[i+1] - x2x->idx[i]);

  /* We increment the max. size to take into account the diagonal entry */

  int  buf_size = stride * stride * (max_size + 1);

  BFT_MALLOC(grows, buf_size, cs_gnum_t);
  BFT_MALLOC(gcols, buf_size, cs_gnum_t);

  if (stride == 1)  { /* Simplified version (equivalent to the interlaced
                         version) */

    for (cs_lnum_t row_id = 0; row_id < n_x; row_id++) {

      const cs_gnum_t  grow_id = rs->g_id[row_id];
      const cs_lnum_t  start = x2x->idx[row_id];
      const cs_lnum_t  end = x2x->idx[row_id+1];

      /* Diagonal term is excluded in the x2x connectivity. Add it. */

      grows[0] = grow_id, gcols[0] = grow_id;

      /* Extra diagonal couples */

      for (cs_lnum_t j = start, i = 1; j < end; j++, i++) {
        grows[i] = grow_id;
        gcols[i] = rs->g_id[x2x->ids[j]];
      }

      cs_matrix_assembler_add_g_ids(ma, end - start + 1, grows, gcols);

    } /* Loop on entities */

  }
  else { /* stride > 1 */

    BFT_MALLOC(g_r_ids, 2*stride, cs_gnum_t);
    g_c_ids = g_r_ids + stride;

    /*
     *   | A_00  | A_01  |  ...  | A_0n  |
     *   |-------|-------|-------|-------|
     *   | A_10  | A_11  |  ...  | A_1n  |
     *   |-------|-------|-------|-------|
     *   | A_n0  | A_n1  |  ...  | A_nn  |
     *
     *  Each block A_.. is n_x * n_x
     */

    for (cs_lnum_t row_id = 0; row_id < n_x; row_id++) {

      const cs_lnum_t  start = x2x->idx[row_id];
      const cs_lnum_t  end = x2x->idx[row_id+1];
      const int  n_entries = (end - start + 1) * stride * stride;

      for (int k = 0; k < stride; k++)
        g_r_ids[k] = rs->g_id[row_id + k*n_x];

      /* Diagonal term is excluded in this connectivity. Add it "manually" */

      int shift = 0;
      for (int ki = 0; ki < stride; ki++) {

        const cs_gnum_t  grow_id = g_r_ids[ki];

        for (int kj = 0; kj < stride; kj++) {

          grows[shift] = grow_id;
          gcols[shift] = g_r_ids[kj];
          shift++;

        }

      }

      /* Extra diagonal couples */

      for (cs_lnum_t j = start; j < end; j++) {

        const cs_lnum_t  col_id = x2x->ids[j];

        for (int k = 0; k < stride; k++)
          g_c_ids[k] = rs->g_id[col_id + k*n_x];

        for (int dof_i = 0; dof_i < stride; dof_i++) {

          const cs_gnum_t  grow_id = g_r_ids[dof_i];
          for (int dof_j = 0; dof_j < stride; dof_j++) {
            grows[shift] = grow_id;
            gcols[shift] = g_c_ids[dof_j];
            shift++;
          }

        }

      } /* Loop on number of DoFs by entity */

      assert(shift == n_entries);
      cs_matrix_assembler_add_g_ids(ma, n_entries, grows, gcols);

    } /* Loop on entities */

    BFT_FREE(g_r_ids);

  } /* stride > 1 */

  /* Now compute structure */

  cs_matrix_assembler_compute(ma);

  /* Free temporary buffers */

  BFT_FREE(grows);
  BFT_FREE(gcols);

  return ma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a matrix assembler and a matrix structure for the given
 *        block structure.
 *
 * \param[in]      forced    define a new structure even if already allocated
 * \param[in]      x2x       underpinning connectivity of the matrix stencil
 * \param[in, out] b         block structure to update
 */
/*----------------------------------------------------------------------------*/

static void
_assign_ma_ms(bool                           forced,
              const cs_adjacency_t          *x2x,
              cs_cdo_system_block_t         *b)
{
  assert(b != NULL);
  assert(x2x != NULL);

  const cs_cdo_system_block_info_t  bi = b->info;

  switch (b->type) {

  case CS_CDO_SYSTEM_BLOCK_DEFAULT:
    /* --------------------------- */
    {
      cs_cdo_system_dblock_t  *db = b->block_pointer;

      if (forced && b->owner) { /* One forces the build of new structures */

        cs_matrix_assembler_destroy(&(db->matrix_assembler));
        cs_matrix_structure_destroy(&(db->matrix_structure));
        db->matrix_assembler = NULL;
        db->matrix_structure = NULL;

      }

      if (db->matrix_assembler == NULL) {

        assert(db->matrix_structure == NULL);

        /* Define the matrix assembler */

        cs_matrix_assembler_t  *ma = NULL;

        if (bi.unrolled) {

          if (bi.interlaced)
            ma = _build_interlaced_ma(bi.stride, x2x, db->range_set);
          else
            ma = _build_no_interlaced_ma(bi.stride, x2x, db->range_set);

        }
        else
          ma = _build_interlaced_ma(1, x2x, db->range_set);

        db->matrix_assembler = ma;

        /* Define the matrix structure */

        db->matrix_structure =
          cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

      }
    }
    break; /* default type of block */

  case CS_CDO_SYSTEM_BLOCK_SPLIT:
    /* ------------------------- */
    {
      cs_cdo_system_sblock_t  *sb = b->block_pointer;

      if (forced && b->owner) { /* One forces the build of new structures */

        cs_matrix_assembler_destroy(&(sb->matrix_assembler));
        cs_matrix_structure_destroy(&(sb->matrix_structure));
        sb->matrix_assembler = NULL;
        sb->matrix_structure = NULL;

      }

      if (sb->matrix_assembler == NULL) {

        assert(sb->matrix_structure == NULL);

        /* Define the matrix assembler */

        cs_matrix_assembler_t  *ma = NULL;
        if (bi.interlaced)
          ma = _build_interlaced_ma(bi.stride, x2x, sb->range_set);
        else
          ma = _build_no_interlaced_ma(bi.stride, x2x, sb->range_set);

        sb->matrix_assembler = ma;

        /* Define the matrix structure */

        sb->matrix_structure =
          cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

      }
    }
    break; /* split type of block */

  default:
    break; /* Do nothing */

  } /* End of switch on the type of block */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find the block corresponding to the same set of metadata in the array
 *        of blocks
 *
 * \param[in]  b    pointer to a block structure to compare
 *
 * \return the id in the array or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static int
_find_in_block_array(cs_cdo_system_block_t   *b)
{
  int  id = -1;

  if (b == NULL)
    return id;

  for (int i = 0; i < _n_cdo_block_structures; i++) {

    cs_cdo_system_block_t  *b_store = _cdo_block_structures[i];

    if (b_store == NULL)        /* deleted block ? */
      continue;

    if (b->type != b_store->type)
      continue;
    if (b->info.matrix_class != b_store->info.matrix_class)
      continue;
    if (b->info.location != b_store->info.location)
      continue;
    if (b->info.n_elements != b_store->info.n_elements)
      continue;
    if (b->info.stride != b_store->info.stride)
      continue;
    if (b->info.unrolled != b_store->info.unrolled)
      continue;
    if (b->info.interlaced != b_store->info.interlaced)
      continue;

    return i; /* All metadata match */

  } /* Loop on blocks stored in the array */

  return id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_cdo_system_block_t structure
 *
 * \param[in, out] p_block    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

static void
_free_block(cs_cdo_system_block_t   **p_block)
{
  if (p_block == NULL)
    return;

  cs_cdo_system_block_t  *b = *p_block;
  if (b == NULL)
    return;

  switch (b->type) {

  case CS_CDO_SYSTEM_BLOCK_DEFAULT:
    /* --------------------------- */
    {
      cs_cdo_system_dblock_t  *db = b->block_pointer;
      assert(db != NULL);

      if (b->owner) {

        cs_matrix_assembler_destroy(&(db->matrix_assembler));
        cs_matrix_structure_destroy(&(db->matrix_structure));

        if (b->info.stride > 1 && b->info.unrolled) {
          cs_range_set_destroy(&(db->range_set));
          cs_interface_set_destroy(&(db->interface_set));
        }

        BFT_FREE(db);

      } /* Block is declared as owner of its structures */
    }
    break;

  case CS_CDO_SYSTEM_BLOCK_SPLIT:
    /* -------------------------- */
    {
      cs_cdo_system_sblock_t  *sb = b->block_pointer;
      assert(sb != NULL);

      if (b->owner) {

        if (sb->matrix_struct_ownership) {

          cs_matrix_assembler_destroy(&(sb->matrix_assembler));
          cs_matrix_structure_destroy(&(sb->matrix_structure));

        }

        BFT_FREE(sb);

      } /* Block is declared as owner of its structures */
    }
    break;

  case CS_CDO_SYSTEM_BLOCK_UNASS:
    /* ------------------------- */
    {
      cs_cdo_system_ublock_t  *ub = b->block_pointer;

      if (ub->_values != NULL) {

        BFT_FREE(ub->_values);
        ub->values = NULL;

      }

      if (!ub->shared_structures) {

        if (b->info.stride > 1) {
          cs_range_set_destroy(&(ub->range_set));
          cs_interface_set_destroy(&(ub->interface_set));
        }

      } /* Block is declared as owner of its structures */

      BFT_FREE(ub);

    }
    break;

  case CS_CDO_SYSTEM_BLOCK_EXT:
    /* ----------------------- */
    {
      cs_cdo_system_xblock_t  *xb = b->block_pointer;
      assert(xb != NULL);

      cs_matrix_assembler_destroy(&(xb->matrix_assembler));
      cs_matrix_structure_destroy(&(xb->matrix_structure));
      cs_range_set_destroy(&(xb->range_set));
      cs_interface_set_destroy(&(xb->interface_set));

      BFT_FREE(xb);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of block. Stop freeing the block structure.\n",
              __func__);

  } /* Switch on the type of block */

  if (b->owner) {

    /* Unset the pointer in the list of shared blocks */

    _cdo_block_structures[b->id] = NULL;

  }

  BFT_FREE(b);
  *p_block = NULL;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize matrix-related structures according to
 *         the type of discretization used for this simulation
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_init_sharing(cs_mesh_t           *mesh,
                                  cs_cdo_connect_t    *connect)
{
  assert(connect != NULL && mesh != NULL); /* Sanity check */

  cs_shared_connect = connect;
  cs_shared_mesh = mesh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a system_helper structure from its set of metadata.
 *        n_col_blocks and n_blocks may differ according to the settings. For
 *        instance, for a saddle-point system, n_col_blocks >= n_blocks
 *
 * \param[in] type             type of system to handle
 * \param[in] n_col_blocks     number of blocks in a row
 * \param[in] col_block_sizes  number of DoFs in each block of the row
 * \param[in] n_blocks         number of blocks associated to this system
 *
 * \return the pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_helper_t *
cs_cdo_system_helper_create(cs_cdo_system_type_t    type,
                            int                     n_col_blocks,
                            const cs_lnum_t        *col_block_sizes,
                            int                     n_blocks)
{
  cs_cdo_system_helper_t  *sh = NULL;

  if (n_col_blocks < 1 || n_blocks < 1)
    return sh;

  BFT_MALLOC(sh, 1, cs_cdo_system_helper_t);

  sh->type = type;
  sh->n_col_blocks = n_col_blocks;

  sh->col_block_sizes = NULL;
  BFT_MALLOC(sh->col_block_sizes, n_col_blocks, cs_lnum_t);

  sh->full_rhs_size = 0;
  for (int i = 0; i < n_col_blocks; i++) {

    sh->col_block_sizes[i] = col_block_sizes[i];
    sh->full_rhs_size +=  col_block_sizes[i];

  }

  sh->max_col_block_sizes = NULL;

  sh->rhs = NULL;
  sh->_rhs = NULL; /* private buffer => manage lifecycle */
  sh->rhs_array = NULL;

  assert(n_blocks > 0);
  sh->n_blocks = n_blocks;
  sh->blocks = NULL;

  BFT_MALLOC(sh->blocks, n_blocks, cs_cdo_system_block_t *);
  for (int i = 0; i < n_blocks; i++)
    sh->blocks[i] = NULL;

  return sh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_cdo_system_helper_t structure
 *
 * \param[in, out] p_helper    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_free(cs_cdo_system_helper_t   **p_helper)
{
  if (p_helper == NULL)
    return;

  cs_cdo_system_helper_t  *sh = *p_helper;
  if (sh == NULL)
    return;

  BFT_FREE(sh->col_block_sizes);
  BFT_FREE(sh->max_col_block_sizes);
  BFT_FREE(sh->rhs_array);      /* array of pointers */
  BFT_FREE(sh->_rhs);
  sh->rhs = NULL;               /* shared pointer */

  for (int i = 0; i < sh->n_blocks; i++)
    _free_block(&(sh->blocks[i]));

  BFT_FREE(sh->blocks);

  BFT_FREE(sh);
  *p_helper = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an unassembled block definition at position "block_id" in the
 *        helper structure Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      matclass    class of the matrix to handle
 * \param[in]      location    where DoFs are defined
 * \param[in]      n_elements  number of elements (support entities for DoFs)
 * \param[in]      stride      number of DoFs by element
 * \param[in]      interlaced  useful if stride > 1; way to store components
 * \param[in]      unrolled    useful if stride > 1; true=as scalar-valued
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_dblock(cs_cdo_system_helper_t       *sh,
                         int                           block_id,
                         cs_cdo_system_matrix_class_t  matclass,
                         cs_flag_t                     location,
                         cs_lnum_t                     n_elements,
                         int                           stride,
                         bool                          interlaced,
                         bool                          unrolled)
{
  if (sh == NULL)
    return NULL;

  if (block_id >= sh->n_blocks)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Block id (%d) is larger than the number of blocks (%d)\n",
              __func__, block_id, sh->n_blocks);

  cs_cdo_system_block_t  *b = NULL;
  BFT_MALLOC(b, 1, cs_cdo_system_block_t);

  b->type = CS_CDO_SYSTEM_BLOCK_DEFAULT;
  b->info.matrix_class = matclass;
  b->info.location = location;
  b->info.n_elements = n_elements;
  b->info.stride = stride;
  b->info.interlaced = interlaced;
  b->info.unrolled = unrolled;

  /* Try to share the most consuming part of the structure (range set,
   * interface set or matrix structures). Check if this type of block with
   * the given metadata already exists */

  int  id = _find_in_block_array(b);

  if (id > -1) { /* Already defined */

    b->block_pointer = _cdo_block_structures[id]->block_pointer;
    b->owner = false;
    b->id = id;

  }
  else { /* Allocate and initialize */

    cs_cdo_system_dblock_t  *db = NULL;
    BFT_MALLOC(db, 1, cs_cdo_system_dblock_t);

    db->matrix = NULL;
    db->mav = NULL;
    db->assembly_func = _assign_assembly_func(b->info);
    db->slave_assembly_func = _assign_slave_assembly_func(b->info);

    if (db->assembly_func == NULL && db->slave_assembly_func == NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: No assembly function set.\n", __func__);

    /* The following structures can be shared if the same block configuration
       is requested */

    db->range_set = NULL;
    db->interface_set = NULL;
    db->matrix_assembler = NULL;
    db->matrix_structure = NULL;

    b->block_pointer = db;
    b->owner = true;

    /* Update the array of reference blocks */

    b->id = _n_cdo_block_structures;
    _n_cdo_block_structures++;
    BFT_REALLOC(_cdo_block_structures, _n_cdo_block_structures,
                cs_cdo_system_block_t *);

    _cdo_block_structures[b->id] = b;

  }

  /* Update the system helper with this new block */

  sh->blocks[block_id] = b;

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a split block definition at position "block_id" in the helper
 *        structure. Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      matclass    class of the matrix to handle
 * \param[in]      location    where DoFs are defined
 * \param[in]      n_elements  number of elements (support entities for DoFs)
 * \param[in]      stride      number of DoFs by element
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_sblock(cs_cdo_system_helper_t       *sh,
                         int                           block_id,
                         cs_cdo_system_matrix_class_t  matclass,
                         cs_flag_t                     location,
                         cs_lnum_t                     n_elements,
                         int                           stride)
{
  if (sh == NULL)
    return NULL;

  if (block_id >= sh->n_blocks)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Block id (%d) is larger than the number of blocks (%d)\n",
              __func__, block_id, sh->n_blocks);

  cs_cdo_system_block_t  *b = NULL;
  BFT_MALLOC(b, 1, cs_cdo_system_block_t);

  b->type = CS_CDO_SYSTEM_BLOCK_SPLIT;
  b->info.matrix_class = matclass;
  b->info.location = location;
  b->info.n_elements = n_elements;
  b->info.stride = stride;
  b->info.interlaced = false;   /* only choice for this type of block */
  b->info.unrolled = false;     /* only choice for this type of block */

  /* Try to share the most consuming part of the structure.
   * Check if this type of block with these metadata already exists */

  int  id = _find_in_block_array(b);

  if (id > -1) { /* Already defined */

    b->block_pointer = _cdo_block_structures[id]->block_pointer;
    b->owner = false;
    b->id = id;

  }
  else { /* Allocate and initialize */

    cs_cdo_system_sblock_t  *sb = NULL;
    BFT_MALLOC(sb, 1, cs_cdo_system_sblock_t);

    sb->n_matrices = stride*stride;

    sb->matrices = NULL;
    BFT_MALLOC(sb->matrices, sb->n_matrices, cs_matrix_t *);
    for (int i = 0; i < sb->n_matrices; i++)
      sb->matrices[i] = NULL;

    sb->mav_array = NULL;
    BFT_MALLOC(sb->mav_array, sb->n_matrices, cs_matrix_assembler_values_t *);
    for (int i = 0; i < sb->n_matrices; i++)
      sb->mav_array[i] = NULL;

    /* The behavior of each sub block corresponds to the one described in
       b_tmp. Find if such block already exists to share the matrix structure
       and matrix assembler */

    cs_cdo_system_block_t  *b_tmp = NULL;
    BFT_MALLOC(b_tmp, 1, cs_cdo_system_block_t);

    b_tmp->type = CS_CDO_SYSTEM_BLOCK_DEFAULT;
    b->info.location = location;
    b->info.n_elements = n_elements;
    b->info.stride = 1;
    b->info.interlaced = false;
    b->info.unrolled = false;

    id = _find_in_block_array(b_tmp);
    if (id > -1)
      sb->matrix_struct_ownership = false;
    else
      sb->matrix_struct_ownership = true;

    sb->assembly_func = _assign_assembly_func(b_tmp->info);
    sb->slave_assembly_func = _assign_slave_assembly_func(b_tmp->info);

    if (sb->assembly_func == NULL && sb->slave_assembly_func == NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: No assembly function set.\n", __func__);

    BFT_FREE(b_tmp);

    /* Shared and private structures */

    sb->range_set = NULL;       /* always shared */
    sb->interface_set = NULL;   /* always shared */

    sb->matrix_assembler = NULL; /* shared status if matrix_struct_ownership */
    sb->matrix_structure = NULL; /* shared status if matrix_struct_ownership */

    b->block_pointer = sb;
    b->owner = true;

    /* Update the array of reference blocks */

    b->id = _n_cdo_block_structures;
    _n_cdo_block_structures++;
    BFT_REALLOC(_cdo_block_structures, _n_cdo_block_structures,
                cs_cdo_system_block_t *);

    _cdo_block_structures[b->id] = b;

  }

  /* Update the system helper with this new block */

  sh->blocks[block_id] = b;

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an unassembled block definition at position "block_id" in the
 *        helper structure Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      adjacency   shared adjacency structure
 * \param[in]      location    where DoFs are defined
 * \param[in]      n_elements  number of elements (support entities for DoFs)
 * \param[in]      stride      number of DoFs by element
 * \param[in]      interlaced  useful if stride > 1; way to store components
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_ublock(cs_cdo_system_helper_t   *sh,
                         int                       block_id,
                         const cs_adjacency_t     *adjacency,
                         cs_flag_t                 location,
                         cs_lnum_t                 n_elements,
                         int                       stride,
                         bool                      interlaced)
{
  if (sh == NULL)
    return NULL;

  if (block_id >= sh->n_blocks)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Block id (%d) is larger than the number of blocks (%d)\n",
              __func__, block_id, sh->n_blocks);

  cs_cdo_system_block_t  *b = NULL;
  BFT_MALLOC(b, 1, cs_cdo_system_block_t);

  b->type = CS_CDO_SYSTEM_BLOCK_UNASS;
  b->info.matrix_class = CS_CDO_SYSTEM_MATRIX_NONE;
  b->info.location = location;
  b->info.n_elements = n_elements;
  b->info.stride = stride;
  b->info.interlaced = interlaced;
  b->info.unrolled = false;   /* only choice for this type of block */

  /* Try to share the most consuming part of the structure.
   * Check if this type of block with these metadata already exists */

  int  id = _find_in_block_array(b);

  if (id > -1) { /* Already defined */

    b->block_pointer = _cdo_block_structures[id]->block_pointer;
    b->owner = false;
    b->id = id;

  }
  else { /* Allocate and initialize */

    cs_cdo_system_ublock_t  *ub = NULL;
    BFT_MALLOC(ub, 1, cs_cdo_system_ublock_t);

    ub->adjacency = adjacency;

    /* Shared and private structures */

    ub->values = ub->_values = NULL;
    ub->shared_structures = false;
    ub->range_set = NULL;
    ub->interface_set = NULL;

    b->block_pointer = ub;
    b->owner = true;

    /* Update the array of reference blocks */

    b->id = _n_cdo_block_structures;
    _n_cdo_block_structures++;
    BFT_REALLOC(_cdo_block_structures, _n_cdo_block_structures,
                cs_cdo_system_block_t *);

    _cdo_block_structures[b->id] = b;

  }

  /* Update the system helper with this new block */

  sh->blocks[block_id] = b;

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an external block definition at position "block_id" in the helper
 *        structure. Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      n_dofs      number of degrees of freedom
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_xblock(cs_cdo_system_helper_t   *sh,
                         int                       block_id,
                         cs_lnum_t                 n_dofs)
{
  if (sh == NULL)
    return NULL;

  if (block_id >= sh->n_blocks)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Block id (%d) is larger than the number of blocks (%d)\n",
              __func__, block_id, sh->n_blocks);

  cs_cdo_system_block_t  *b = NULL;
  BFT_MALLOC(b, 1, cs_cdo_system_block_t);

  b->type = CS_CDO_SYSTEM_BLOCK_EXT;
  b->info.matrix_class = CS_CDO_SYSTEM_MATRIX_CS;
  b->info.location = 0;         /* not defined */
  b->info.n_elements = n_dofs;
  b->info.stride = 1;
  b->info.interlaced = false;   /* not used */
  b->info.unrolled = false;     /* not used */

  cs_cdo_system_xblock_t  *xb = NULL;
  BFT_MALLOC(xb, 1, cs_cdo_system_xblock_t);

  xb->matrix = NULL;
  xb->mav = NULL;

  /* Private structures */

  xb->range_set = NULL;
  xb->interface_set = NULL;
  xb->matrix_assembler = NULL;
  xb->matrix_structure = NULL;

  b->block_pointer = xb;

  /* Always owner for this definition */

  b->owner = true;

  /* Update the array of reference blocks */

  b->id = _n_cdo_block_structures;
  _n_cdo_block_structures++;
  BFT_REALLOC(_cdo_block_structures, _n_cdo_block_structures,
              cs_cdo_system_block_t *);

  _cdo_block_structures[b->id] = b;

  /* Update the system helper with this new block */

  sh->blocks[block_id] = b;

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the range set structure associated to the given block_id
 *
 * \param[in]  sh         pointer to the system_helper structure to update
 * \param[in]  block_id   id of the block to consider
 *
 * \return a pointer to a range set structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_range_set_t *
cs_cdo_system_get_range_set(const cs_cdo_system_helper_t  *sh,
                            int                            block_id)
{
  if (sh == NULL)
    return NULL;
  if (block_id < 0 || block_id >= sh->n_blocks)
    return NULL;

  cs_range_set_t  *rs = NULL;

  cs_cdo_system_block_t  *b = sh->blocks[block_id];

  switch (b->type) {

  case CS_CDO_SYSTEM_BLOCK_DEFAULT:
    {
      cs_cdo_system_dblock_t  *db = b->block_pointer;
      rs = db->range_set;
    }
    break;

  case CS_CDO_SYSTEM_BLOCK_EXT:
    {
      cs_cdo_system_xblock_t  *xb = b->block_pointer;
      rs = xb->range_set;
    }
    break;

  case CS_CDO_SYSTEM_BLOCK_SPLIT:
    {
      cs_cdo_system_sblock_t  *sb = b->block_pointer;
      rs = sb->range_set;
    }
    break;

  case CS_CDO_SYSTEM_BLOCK_UNASS:
    {
      cs_cdo_system_ublock_t  *ub = b->block_pointer;
      rs = ub->range_set;
    }
    break;

  default:
    break;
  }

  return rs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the matrix associated to the given block_id. If the type of
 *        block is either CS_CDO_SYSTEM_BLOCK_DEFAULT or
 *        CS_CDO_SYSTEM_BLOCK_EXT. In other cases, a NULL pointer is
 *        returned. The unassembled block has no matrix and to get a matrix of
 *        a split block, one should use cs_cdo_system_get_sub_matrix(sh,
 *        block_id, sub_id)
 *
 * \param[in]  sh         pointer to the system_helper structure to update
 * \param[in]  block_id   id of the block to consider
 *
 * \return a pointer to a cs_matrix_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_cdo_system_get_matrix(const cs_cdo_system_helper_t  *sh,
                         int                            block_id)
{
  if (sh == NULL)
    return NULL;
  if (block_id < 0 || block_id >= sh->n_blocks)
    return NULL;

  cs_matrix_t  *m = NULL;

  cs_cdo_system_block_t  *b = sh->blocks[block_id];

  switch (b->type) {

  case CS_CDO_SYSTEM_BLOCK_DEFAULT:
    {
      cs_cdo_system_dblock_t  *db = b->block_pointer;
      m = db->matrix;
    }
    break;

  case CS_CDO_SYSTEM_BLOCK_EXT:
    {
      cs_cdo_system_xblock_t  *xb = b->block_pointer;
      m = xb->matrix;
    }
    break;

  default:
    break;
  }

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the (sub-)matrix associated to a split block with id equal
 *        to block_id. sub_id is the id in the list of matrices of size equal
 *        to stride*stride.
 *        If the type of the block is not CS_CDO_SYSTEM_BLOCK_SPLIT, then a
 *        NULL pointer is returned.
 *
 * \param[in, out]  sh         pointer to the system_helper structure to update
 * \param[in]       block_id   id of the block to consider
 *
 * \return a pointer to a cs_matrix_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_cdo_system_get_sub_matrix(cs_cdo_system_helper_t  *sh,
                             int                      block_id,
                             int                      sub_id)
{
  if (sh == NULL)
    return NULL;
  if (block_id < 0 || block_id >= sh->n_blocks)
    return NULL;

  cs_matrix_t  *m = NULL;

  cs_cdo_system_block_t  *b = sh->blocks[block_id];

  switch (b->type) {

  case CS_CDO_SYSTEM_BLOCK_SPLIT:
    {
      cs_cdo_system_sblock_t  *sb = b->block_pointer;
      if (sub_id > -1 && sub_id < sb->n_matrices)
        m = sb->matrices[sub_id];
    }
    break;

  default:
    break; /* Do nothing */
  }

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the associated structures for the given system_helper structure
 *        If a member is already allocated, one keeps the member as it is.
 *
 * \param[in, out]  sh         pointer to the system_helper structure to update
 * \param[in]       block_id   specific block to handle or -1 for all blocks
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_build_block(cs_cdo_system_helper_t  *sh,
                          int                      block_id)
{
  if (sh == NULL)
    return;

  if (block_id >= sh->n_blocks)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for the block_id (%d); n_blocks=%d\n",
              __func__, block_id, sh->n_blocks);

  cs_cdo_connect_t  *connect = cs_shared_connect;
  int _n_blocks = (block_id == -1) ? sh->n_blocks : 1;

  for (int i = 0; i < _n_blocks; i++) {

    cs_cdo_system_block_t  *b = NULL;
    if (i == 0 && block_id > -1)
      b = sh->blocks[block_id];
    else
      b = sh->blocks[i];
    assert(b != NULL);

    if (b->type == CS_CDO_SYSTEM_BLOCK_EXT)
      continue; /* All is defined by the calling code. Not generic */
    if (b->owner == false)
      continue;

    switch (b->info.location) {

    case CS_FLAG_LOCATION_PRIMAL_VTX:
    case CS_FLAG_LOCATION_DUAL_CELL:
      assert(connect->n_vertices == b->info.n_elements);
      assert(connect->n_vertices == connect->v2v->n_elts);
      _assign_ifs_rset(false, b, connect->vtx_ifs, connect->vtx_rset);
      if (b->type != CS_CDO_SYSTEM_BLOCK_UNASS)
        _assign_ma_ms(false, connect->v2v, b);
      break;

    case CS_FLAG_LOCATION_PRIMAL_EDGE:
    case CS_FLAG_LOCATION_DUAL_FACE:
      assert(connect->n_edges == b->info.n_elements);
      assert(connect->n_edges == connect->e2e->n_elts);
      _assign_ifs_rset(false, b, connect->edge_ifs, connect->edge_rset);
      if (b->type != CS_CDO_SYSTEM_BLOCK_UNASS)
        _assign_ma_ms(false, connect->e2e, b);
      break;

    case CS_FLAG_LOCATION_PRIMAL_FACE:
    case CS_FLAG_LOCATION_DUAL_EDGE:
      assert(connect->n_faces[CS_ALL_FACES] == b->info.n_elements);
      assert(connect->n_faces[CS_ALL_FACES] == connect->f2f->n_elts);
      _assign_ifs_rset(false, b, connect->face_ifs, connect->face_rset);
      if (b->type != CS_CDO_SYSTEM_BLOCK_UNASS)
        _assign_ma_ms(false, connect->f2f, b);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid location for a system_helper structure.\n",
                __func__);

    } /* End of switch */

  } /* Loop on blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize the matrix, rhs and the matrix assembler
 *        values. If p_rhs is NULL then one allocates the rhs inside this
 *        function. The ownership is transfered to this structure in that case.
 *
 * \param[in, out] sh       pointer to a system helper structure
 * \param[in, out] p_rhs    double pointer to the RHS array to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_init_system(cs_cdo_system_helper_t    *sh,
                                 cs_real_t                **p_rhs)
{
  if (sh == NULL)
    return;

  /* Right-hand side */

  cs_real_t *rhs = *p_rhs;
  if (rhs == NULL) {

    BFT_MALLOC(sh->_rhs, sh->full_rhs_size, cs_real_t);
    *p_rhs = sh->_rhs;
    sh->rhs = sh->_rhs;

    if (sh->n_col_blocks > 1) {

      if (sh->rhs_array == NULL)
        BFT_MALLOC(sh->rhs_array, sh->n_col_blocks, cs_real_t *);

      cs_lnum_t  shift = 0;
      for (int k = 0; k < sh->n_col_blocks; k++) {
        sh->rhs_array[k] = sh->rhs + shift;
        shift += sh->col_block_sizes[k];
      }

    }

  } /* rhs is managed by the helper */
  else
    sh->rhs = rhs;

#if defined(HAVE_OPENMP)
# pragma omp parallel for if  (sh->full_rhs_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < sh->full_rhs_size; i++) sh->rhs[i] = 0.0;
#else
  memset(sh->rhs, 0, sh->full_rhs_size*sizeof(cs_real_t));
#endif

  /* Initialize structures */

  for (int i = 0; i < sh->n_blocks; i++) {

    cs_cdo_system_block_t  *b = sh->blocks[i];

    switch (b->type) {

    case CS_CDO_SYSTEM_BLOCK_DEFAULT:
      {
        cs_cdo_system_dblock_t  *db = b->block_pointer;

        /* Matrix */

        if (db->matrix != NULL) {
          cs_matrix_release_coefficients(db->matrix);
          cs_matrix_destroy(&(db->matrix));
        }

        assert(db->matrix_structure != NULL);
        db->matrix = cs_matrix_create(db->matrix_structure);

#if defined(HAVE_HYPRE)
        if (b->info.matrix_class == CS_CDO_SYSTEM_MATRIX_HYPRE) {

          int device_id = cs_get_device_id();
          int use_device = (device_id < 0) ? 0 : 1;

          cs_matrix_set_type_hypre(db->matrix, use_device);

        }
#endif

        /* Matrix assembler values */

        if (db->mav != NULL)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Matrix assembler values has not been finalized.\n",
                    __func__);

        cs_lnum_t db_size = 1;
        cs_lnum_t eb_size = 1;

        if (!b->info.unrolled) {

          db_size = b->info.stride;
          eb_size = b->info.stride;

        }

        db->mav = cs_matrix_assembler_values_init(db->matrix, db_size, eb_size);
      }
      break;

    case CS_CDO_SYSTEM_BLOCK_SPLIT:
      {
        cs_cdo_system_sblock_t  *sb = b->block_pointer;
        assert(sb->matrices != NULL);
        assert(sb->mav_array != NULL);

        for (int k = 0; k < sb->n_matrices; k++) {

          /* Matrices */

          if (sb->matrices[k] != NULL) {
            cs_matrix_release_coefficients(sb->matrices[k]);
            cs_matrix_destroy(&(sb->matrices[k]));
          }

          assert(sb->matrix_structure != NULL);
          sb->matrices[k] = cs_matrix_create(sb->matrix_structure);

#if defined(HAVE_HYPRE)
        if (b->info.matrix_class == CS_CDO_SYSTEM_MATRIX_HYPRE) {

          int device_id = cs_get_device_id();
          int use_device = (device_id < 0) ? 0 : 1;

          cs_matrix_set_type_hypre(sb->matrices[k], use_device);

        }
#endif
          /* Matrix assembler values */

          if (sb->mav_array[k] != NULL)
            bft_error(__FILE__, __LINE__, 0,
                      "%s: Matrix assembler values has not been finalized.\n",
                      __func__);

          sb->mav_array[k] = cs_matrix_assembler_values_init(sb->matrices[k],
                                                             1, 1);

        } /* Loop on each matrix */
      }
      break;

    case CS_CDO_SYSTEM_BLOCK_EXT:
      {
        cs_cdo_system_xblock_t  *xb = b->block_pointer;
        assert(xb->matrix_structure != NULL);

        /* Matrix */

        if (xb->matrix != NULL) {
          cs_matrix_release_coefficients(xb->matrix);
          cs_matrix_destroy(&(xb->matrix));
        }

        assert(xb->matrix_structure != NULL);
        xb->matrix = cs_matrix_create(xb->matrix_structure);

        /* Matrix assembler values */

        if (xb->mav != NULL)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Matrix assembler values has not been finalized.\n",
                    __func__);

        xb->mav = cs_matrix_assembler_values_init(xb->matrix, 1, 1);
      }
      break;

    default:
      break; /* Do nothing for unassembled blocks */

    } /* End of switch on block type */

  } /* Loop on blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the assembly after the cellwise building and assembly
 *
 * \param[in, out]      sh       pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_finalize_assembly(cs_cdo_system_helper_t    *sh)
{
  if (sh == NULL)
    return;

  for (int i = 0; i < sh->n_blocks; i++) {

    cs_cdo_system_block_t  *b = sh->blocks[i];

    switch (b->type) {

    case CS_CDO_SYSTEM_BLOCK_DEFAULT:
      /* --------------------------- */
      {
        cs_cdo_system_dblock_t  *db = b->block_pointer;

        cs_matrix_assembler_values_done(db->mav);
        cs_matrix_assembler_values_finalize(&(db->mav));
      }
      break;

    case CS_CDO_SYSTEM_BLOCK_SPLIT:
      /* ------------------------- */
      {
        cs_cdo_system_sblock_t  *sb = b->block_pointer;

        for (int k = 0; k < sb->n_matrices; k++) {
          cs_matrix_assembler_values_done(sb->mav_array[k]);
          cs_matrix_assembler_values_finalize(&(sb->mav_array[k]));
        }
      }
      break;

    case CS_CDO_SYSTEM_BLOCK_EXT:
      /* ----------------------- */
      {
        cs_cdo_system_xblock_t  *xb = b->block_pointer;

        cs_matrix_assembler_values_done(xb->mav);
        cs_matrix_assembler_values_finalize(&(xb->mav));
      }
      break;

    default:
      /* CS_CDO_SYSTEM_BLOCK_UNASS */
      break; /* Do nothing */

    } /* Switch on the type of block */

  } /* Loop on blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free matrix and rhs after the solve step
 *
 * \param[in, out]      sh       pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_reset(cs_cdo_system_helper_t    *sh)
{
  if (sh == NULL)
    return;

  BFT_FREE(sh->_rhs);
  sh->rhs = NULL;

  /* Free matrix (or matrices) */

  for (int i = 0; i < sh->n_blocks; i++) {

    cs_cdo_system_block_t  *b = sh->blocks[i];

    switch (b->type) {

    case CS_CDO_SYSTEM_BLOCK_DEFAULT:
      /* --------------------------- */
      {
        cs_cdo_system_dblock_t  *db = b->block_pointer;

        cs_matrix_release_coefficients(db->matrix);
        cs_matrix_destroy(&(db->matrix));
      }
      break;

    case CS_CDO_SYSTEM_BLOCK_SPLIT:
      /* ------------------------- */
      {
        cs_cdo_system_sblock_t  *sb = b->block_pointer;

        for (int k = 0; k < sb->n_matrices; k++) {
          cs_matrix_release_coefficients(sb->matrices[k]);
          cs_matrix_destroy(&(sb->matrices[k]));
        }
      }
      break;

    case CS_CDO_SYSTEM_BLOCK_EXT:
      /* ----------------------- */
      {
        cs_cdo_system_xblock_t  *xb = b->block_pointer;

        cs_matrix_release_coefficients(xb->matrix);
        cs_matrix_destroy(&(xb->matrix));
      }
      break;

    default:
      /* CS_CDO_SYSTEM_BLOCK_UNASS */
      break; /* Do nothing */

    } /* Switch on the type of block */

  } /* Loop on blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize shared assembly structures from the existing helper
 *        structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_allocate_assembly(void)
{
  if (_n_cdo_block_structures < 1)
    return;

  const cs_cdo_connect_t  *connect = cs_shared_connect;

  int  n_max_cw_dofs = 0, max_ddim = 1, max_edim = 1;

  for (int i = 0; i < _n_cdo_block_structures; i++) {

    cs_cdo_system_block_t  *b = _cdo_block_structures[i];

    if (b == NULL) /* deleted block */
      continue;

    cs_cdo_system_block_info_t  bi = b->info;

    max_ddim = CS_MAX(max_ddim, bi.stride);
    max_edim = CS_MAX(max_edim, bi.stride);

    switch (bi.location) {

    case CS_FLAG_LOCATION_PRIMAL_VTX:
    case CS_FLAG_LOCATION_DUAL_CELL:
      n_max_cw_dofs = CS_MAX(n_max_cw_dofs, connect->n_max_vbyc);
      break;

    case CS_FLAG_LOCATION_PRIMAL_EDGE:
    case CS_FLAG_LOCATION_DUAL_FACE:
      n_max_cw_dofs = CS_MAX(n_max_cw_dofs, connect->n_max_ebyc);
      break;

    case CS_FLAG_LOCATION_PRIMAL_FACE:
    case CS_FLAG_LOCATION_DUAL_EDGE:
      n_max_cw_dofs = CS_MAX(n_max_cw_dofs, connect->n_max_fbyc);
      break;

    default:
      break; /* Do nothing */

    } /* End of switch */

  } /* Loop on system helper structures */

  /* Allocate shared buffers used during the assembly step */

  cs_cdo_assembly_init(max_ddim, max_edim, n_max_cw_dofs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all members associated to system helpers
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_destroy_all(void)
{
  /* Free remaining blocks if needed */

  for (int i = 0; i < _n_cdo_block_structures; i++)
    _free_block(&(_cdo_block_structures[i]));

  BFT_FREE(_cdo_block_structures);

  _n_cdo_block_structures = 0;
  _cdo_block_structures = NULL;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
