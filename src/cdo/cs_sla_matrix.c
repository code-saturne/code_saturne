/*============================================================================
 * Sparse linear algebra related to CDO techniques
 * Matrix operations (matvec, matmat, and other manipulations)
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_sort.h"
#include "cs_search.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sla.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define  SLA_MATRIX_DEBUG 0

struct _sla_mat_prop_t {

  double  spectral_radius;
  double  eps;
  double  iter_count;

};

/* Sparse Accumulator: see Gilbert etal. and Buluc (PhD) */
/* ===================================================== */

typedef struct { /* for CSR or MSR matrix operations */

  cs_lnum_t  *tag;  /* size n_cols */
  double     *val;  /* size n_cols */

  size_t      size; /* size of the allocated list */
  size_t      nnz;  /* number of non empty in the list */
  cs_lnum_t  *lst;

} _spa_t;

typedef struct { /* for DEC matrix operations */

  cs_lnum_t  *tag;  /* size n_cols */
  short int  *sgn;  /* size n_cols */

  size_t      size; /* size of the allocated list */
  size_t      nnz;  /* number of non empty in the list */
  cs_lnum_t  *lst;

} _spa_dec_t;

/*============================================================================
 * Private constant variables
 *============================================================================*/

static const char  _sla_err_stride[] =
  "  Incompatible stride value (>1).\n   Stop matrix computation.\n";

static const char _sla_matrix_type[CS_SLA_MAT_N_TYPES][CS_CDO_LEN_NAME] = {
  "None",
  "DEC",
  "CSR",
  "MSR"
};

/*============================================================================
 * Private function prototypes
 *============================================================================*/

static cs_sla_mat_prop_t *
_get_default_matrix_properties(void)
{
  cs_sla_mat_prop_t  *mp = NULL;

  BFT_MALLOC(mp, 1, cs_sla_mat_prop_t);

  mp->spectral_radius = -1;
  mp->eps = 1.e16;
  mp->iter_count = 0;

  return mp;
}

/*----------------------------------------------------------------------------*/
/*                         MATRIX MULTIPLICATION                              */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a Sparse Accumulator (SPA) for MSR and CSR type
 *
 * \param[in]  a_size     size of tag and val arrays
 * \param[in]  l_size     size of lst
 *
 * \return a pointer to a new allocated and initialize SPA struct.
 */
/*----------------------------------------------------------------------------*/

static _spa_t *
_spa_init(size_t   a_size,
          size_t   l_size)
{
  size_t  i;

  _spa_t  *spa = NULL;

  BFT_MALLOC(spa, 1, _spa_t);

  BFT_MALLOC(spa->tag, a_size, cs_lnum_t);
  BFT_MALLOC(spa->val, a_size, double);

  for (i = 0; i < a_size; i++) {
    spa->tag[i] = -1; /* not used by default */
    spa->val[i] = 0.0;
  }

  spa->size = l_size;
  spa->nnz = 0;
  BFT_MALLOC(spa->lst, l_size, cs_lnum_t);

  return spa;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a Sparse Accumulator (SPA) for MSR and CSR type
 *
 * \param[inout]   spa    pointer to a SPA struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static _spa_t *
_spa_free(_spa_t  *spa)
{
  if (spa == NULL)
    return spa;

  BFT_FREE(spa->tag);
  BFT_FREE(spa->val);
  BFT_FREE(spa->lst);

  BFT_FREE(spa);

  return spa;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a Sparse Accumulator (SPA) for DEC type
 *
 * \param[in]  a_size     size of tag and connect arrays
 * \param[in]  l_size     size of lst
 *
 * \return a pointer to a new allocated and initialize SPA struct.
 */
/*----------------------------------------------------------------------------*/

static _spa_dec_t *
_spa_dec_init(size_t   a_size,
              size_t   l_size)
{
  size_t  i;

  _spa_dec_t  *spa = NULL;

  BFT_MALLOC(spa, 1, _spa_dec_t);

  BFT_MALLOC(spa->tag, a_size, cs_lnum_t);
  BFT_MALLOC(spa->sgn, a_size, short int);

  for (i = 0; i < a_size; i++) {
    spa->tag[i] = -1; /* not used by default */
    spa->sgn[i] = 0;
  }

  spa->size = l_size;
  spa->nnz = 0;
  BFT_MALLOC(spa->lst, l_size, cs_lnum_t);

  return spa;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a Sparse Accumulator (SPA) for MSR and CSR type
 *
 * \param[inout]   spa    pointer to a SPA struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static _spa_dec_t *
_spa_dec_free(_spa_dec_t  *spa)
{
  if (spa == NULL)
    return spa;

  BFT_FREE(spa->tag);
  BFT_FREE(spa->sgn);
  BFT_FREE(spa->lst);

  BFT_FREE(spa);

  return spa;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of a new entry
 *
 * \param[inout]  spa     pointer to the SPA struct. to update
 * \param[in]     value   value to add
 * \param[in]     pos     position related to this value
 * \param[in]     row_id  current row
 *
 */
/*----------------------------------------------------------------------------*/

inline static void
_spa_add(_spa_t     *spa,
         double      value,
         cs_lnum_t   pos,
         cs_lnum_t   row_id)
{
  if (spa->tag[pos] != row_id) { /* Add a new item in lst */

    if (spa->nnz == spa->size) {  /* Increase size */
      spa->size = CS_MAX(spa->size + 1, 2*spa->size);
      BFT_REALLOC(spa->lst, spa->size, cs_lnum_t);
    }
    spa->lst[spa->nnz] = pos;
    spa->nnz += 1;

    spa->tag[pos] = row_id;
    spa->val[pos] = value;

  }
  else /* Entry already defined */
    spa->val[pos] += value;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of a new entry
 *
 * \param[inout]  spa     pointer to the SPA struct. to update
 * \param[in]     value   value to add
 * \param[in]     pos     position related to this value
 * \param[in]     row_id  current row
 *
 */
/*----------------------------------------------------------------------------*/

static inline void
_spa_dec_add(_spa_dec_t   *spa,
             short int     value,
             cs_lnum_t     pos,
             cs_lnum_t     row_id)
{
  if (spa->tag[pos] != row_id) { /* Add a new item in lst */

    if (spa->nnz == spa->size) {  /* Increase size */
      spa->size = CS_MAX(spa->size + 1, 2*spa->size);
      BFT_REALLOC(spa->lst, spa->size, cs_lnum_t);
    }
    spa->lst[spa->nnz] = pos;
    spa->nnz += 1;

    spa->tag[pos] = row_id;
    spa->sgn[pos] = value;

  }
  else /* Entry already defined */
    spa->sgn[pos] += value;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Gather data from the current SPA state to update matrix structure
 *          Do not add zero entries.
 *
 * \param[in]    spa      pointer to the SPA struct.
 * \param[in]    idx      current index position
 * \param[inout] col_id   pointer to matrix->col_id array
 * \param[inout] val      pointer to matrix->val array
 *
 * \return the number of entries added
 */
/*----------------------------------------------------------------------------*/

static size_t
_spa_gather(_spa_t     *spa,
            size_t      idx,
            cs_lnum_t  *col_id,
            double     *val)
{
  size_t  i, shift;

  size_t  n_entries = 0;

  for (i = 0; i < spa->nnz; i++) {

    cs_lnum_t  pos = spa->lst[i];
    double  value = spa->val[pos];

    /* Do not affect an entry with zero value */
    if (value > DBL_EPSILON || value < -DBL_EPSILON) {
      shift = idx + n_entries;
      col_id[shift] = pos;
      val[shift] = value;
      n_entries++;
    }

    /* Initialize for a future use */
    spa->val[pos] = 0.0;

  } /* End of loop on elements of SPA */

  /* Reset */
  spa->nnz = 0;

  return n_entries;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Gather data from the current SPA state to update matrix structure
 *          Do not add zero entries.
 *
 * \param[in]    spa     pointer to the SPA struct.
 * \param[in]    idx     current index position
 * \param[inout] col_id  pointer to matrix->col_id array
 * \param[inout] val     pointer to matrix->sgn array
 *
 * \return the number of entries added
 */
/*----------------------------------------------------------------------------*/

static size_t
_spa_dec_gather(_spa_dec_t  *spa,
                size_t       idx,
                cs_lnum_t   *col_id,
                short int   *connect)
{
  size_t  i, shift;

  size_t  n_entries = 0;

  for (i = 0; i < spa->nnz; i++) {

    cs_lnum_t  pos = spa->lst[i];
    short int  value = spa->sgn[pos];

    /* Do not affect an entry with zero value */
    if (value != 0) {
      shift = idx + n_entries;
      col_id[shift] = pos;
      connect[shift] = value;
      n_entries++;
    }

    /* Initialize for a future use */
    spa->sgn[pos] = 0;

  } /* End of loop on elements of SPA */

  /* Reset */
  spa->nnz = 0;

  return n_entries;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Estimate roughly the number of nnz of the matrix c = a*b and the
 *          stencil of c.
 *          Size estimation is given by E. Cohen (1997) in
 *          Structure Prediction and Computation of Sparse Matrix Product
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 * \param[inout] nnz      size_t element
 * \param[inout] stencil  size_t element
 *
 */
/*----------------------------------------------------------------------------*/

static void
_estimate_sizes(const cs_sla_matrix_t   *a,
                const cs_sla_matrix_t   *b,
                size_t                  *nnz,
                size_t                  *stencil)
{
  size_t  as = a->idx[a->n_rows], bs = b->idx[b->n_rows];
  size_t  ab = as*bs, bb = b->n_rows*b->n_cols;

  if (a->n_cols == 0)
    *nnz = 0;
  else
    *nnz = CS_MAX(5, ab/a->n_cols);

  if (bb == 0)
    *stencil = 0;
  else
    *stencil = CS_MAX(2, ab/bb);

  if (*nnz == 0 || *stencil == 0)
    printf(" << WARNING >> operation on an empty matrix. Check data !\n");
}

/*----------------------------------------------------------------------------
  Allocate a CSR matrix "c" coming from "c = a*b"
  ----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_mat(const cs_sla_matrix_t   *a,
          const cs_sla_matrix_t   *b,
          cs_sla_matrix_type_t     type,
          size_t                   guess_size)
{
  cs_sla_matrix_t  *c = NULL;
  _Bool  sym = false;

  if ((a->flag & CS_SLA_MATRIX_SYM) &&
      (b->flag & CS_SLA_MATRIX_SYM))
    sym = true;

  c = cs_sla_matrix_create(a->n_rows, b->n_cols, 1, type, sym);

  BFT_MALLOC(c->col_id, guess_size, cs_lnum_t);

  if (type == CS_SLA_MAT_CSR || type == CS_SLA_MAT_MSR)
    BFT_MALLOC(c->val, guess_size, double);
  else if (type == CS_SLA_MAT_DEC)
    BFT_MALLOC(c->sgn, guess_size, short int);

  return c;
}

/*----------------------------------------------------------------------------
  Resize a CSR matrix
  ----------------------------------------------------------------------------*/

static void
_resize_mat(cs_sla_matrix_t   *mat,
            size_t             cur_size,
            size_t            *max_size)
{
  size_t  ms = *max_size;

  if (cur_size + 1 > ms) { /* Realloc */

    ms = CS_MAX(cur_size + 1, floor(1.3*ms));
    BFT_REALLOC(mat->col_id, ms, cs_lnum_t);

    if (mat->type == CS_SLA_MAT_CSR || mat->type == CS_SLA_MAT_MSR)
      BFT_REALLOC(mat->val, ms, double);
    else if (mat->type == CS_SLA_MAT_DEC)
      BFT_REALLOC(mat->sgn, ms, short int);

  }
  *max_size = ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute c = a*b for a and b in DEC storage
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 *
 * \return a pointer to matrix struct. c
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_multiply_dec_matrices(const cs_sla_matrix_t  *a,
                       const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  short int  val;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_dec_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity check */
  assert(a->type == CS_SLA_MAT_DEC && b->type == CS_SLA_MAT_DEC);

  /* Initialize structures */
  _estimate_sizes(a, b, &size_max, &lst_size_init);
  c = _init_mat(a, b, CS_SLA_MAT_DEC, size_max);
  spa = _spa_dec_init(b->n_cols, lst_size_init);

  /* Compute c = a*b line by line */
  for (ii = 0; ii < a->n_rows; ii++) {

    for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {

      jj = a->col_id[j], val = a->sgn[j];
      for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
        kk = b->col_id[k], _spa_dec_add(spa, val*b->sgn[k], kk, ii);

    } /* End of loop on non-empty columnd of row ii of a */

    /* Fill a new row in c and prepare next step */
    if (spa->nnz + idx > size_max)
      _resize_mat(c, spa->nnz + idx, &size_max);
    n_entries = _spa_dec_gather(spa, idx, c->col_id, c->sgn);
    c->idx[ii+1] = idx + n_entries;
    idx = c->idx[ii+1];

  } /* End of loop on row ii of a */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->sgn, c->idx[c->n_rows], short int);
  spa = _spa_dec_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute c = a*b for a and b in CSR storage
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 *
 * \return a pointer to matrix struct. c
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_multiply_csr_matrices(const cs_sla_matrix_t  *a,
                       const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  double  val;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity check */
  assert(a->type == CS_SLA_MAT_CSR && b->type == CS_SLA_MAT_CSR);

  /* Initialize structures */
  _estimate_sizes(a, b, &size_max, &lst_size_init);
  c = _init_mat(a, b, CS_SLA_MAT_CSR, size_max);
  spa = _spa_init(b->n_cols, lst_size_init);

  /* Compute c = a*b line by line */
  for (ii = 0; ii < a->n_rows; ii++) {

    for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {

      jj = a->col_id[j], val = a->val[j];
      for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
        kk = b->col_id[k], _spa_add(spa, val*b->val[k], kk, ii);

    } /* End of loop on non-empty columns of row ii of a */

    /* Fill a new row in c and prepare next step */
    if (spa->nnz + idx > size_max)
      _resize_mat(c, spa->nnz + idx, &size_max);
    n_entries = _spa_gather(spa, idx, c->col_id, c->val);
    c->idx[ii+1] = idx + n_entries;
    idx = c->idx[ii+1];

  } /* End of loop on row ii of a */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute c = a*b for a in DEC storage and b in CSR storage
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 *
 * \return a pointer to matrix struct. c
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_multiply_deccsr_matrices(const cs_sla_matrix_t  *a,
                          const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  double  val;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity checks */
  assert(a->type == CS_SLA_MAT_DEC && b->type == CS_SLA_MAT_CSR);

  /* Initialize structures */
  _estimate_sizes(a, b, &size_max, &lst_size_init);
  c = _init_mat(a, b, CS_SLA_MAT_CSR, size_max);
  spa = _spa_init(b->n_cols, lst_size_init);

  /* Compute c = a*b line by line */
  for (ii = 0; ii < a->n_rows; ii++) {

    for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {

      jj = a->col_id[j], val = a->sgn[j];
      for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
        kk = b->col_id[k], _spa_add(spa, val*b->val[k], kk, ii);

    } /* End of loop on non-empty columns of row ii of a */

    /* Fill a new row in c and prepare next step */
    if (spa->nnz + idx > size_max)
      _resize_mat(c, spa->nnz + idx, &size_max);
    n_entries = _spa_gather(spa, idx, c->col_id, c->val);
    c->idx[ii+1] = idx + n_entries;
    idx = c->idx[ii+1];

  } /* End of loop on row ii of a */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute c = a*b for a in CSR storage and b in DEC storage
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 *
 * \return a pointer to matrix struct. c
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_multiply_csrdec_matrices(const cs_sla_matrix_t  *a,
                          const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  double  val;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity checks */
  assert(a->type == CS_SLA_MAT_CSR && b->type == CS_SLA_MAT_DEC);

  /* Initialize structures */
  _estimate_sizes(a, b, &size_max, &lst_size_init);
  c = _init_mat(a, b, CS_SLA_MAT_CSR, size_max);
  spa = _spa_init(b->n_cols, lst_size_init);

  /* Compute c = a*b line by line */
  for (ii = 0; ii < a->n_rows; ii++) {

    for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {

      jj = a->col_id[j], val = a->val[j];
      for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
        kk = b->col_id[k], _spa_add(spa, val*b->sgn[k], kk, ii);

    } /* End of loop on non-empty columns of row ii of a */

    /* Fill a new row in c and prepare next step */
    if (spa->nnz + idx > size_max)
      _resize_mat(c, spa->nnz + idx, &size_max);
    n_entries = _spa_gather(spa, idx, c->col_id, c->val);
    c->idx[ii+1] = idx + n_entries;
    idx = c->idx[ii+1];

  } /* End of loop on row ii of a */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute c = a*b for a in DEC storage and b in MSR storage
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 *
 * \return a pointer to matrix struct. c
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_multiply_decmsr_matrices(const cs_sla_matrix_t  *a,
                          const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  double  val;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity checks */
  assert(a->type == CS_SLA_MAT_DEC && b->type == CS_SLA_MAT_MSR);

  /* Initialize structures */
  _estimate_sizes(a, b, &size_max, &lst_size_init);
  c = _init_mat(a, b, CS_SLA_MAT_CSR, size_max);
  spa = _spa_init(b->n_cols, lst_size_init);

  for (ii = 0; ii < a->n_rows; ii++) {

    for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {

      jj = a->col_id[j], val = a->sgn[j];

      /* Diagonal term */
      _spa_add(spa, val*b->diag[jj], jj, ii);

      /* Extra diag term */
      for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
        kk = b->col_id[k], _spa_add(spa, val*b->val[k], kk, ii);

    } /* End of loop on non-empty columns of row ii of a */

    /* Fill a new row in c and prepare next step */
    if (spa->nnz + idx > size_max)
      _resize_mat(c, spa->nnz + idx, &size_max);
    n_entries = _spa_gather(spa, idx, c->col_id, c->val);
    c->idx[ii+1] = idx + n_entries;
    idx = c->idx[ii+1];

  } /* End of loop on row ii of a */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute c = a*b for a in MSR storage and b in DEC storage
 *
 * \param[in]    a        pointer to the matrix struct. a
 * \param[in]    b        pointer to the matrix struct. b
 *
 * \return a pointer to matrix struct. c
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_multiply_msrdec_matrices(const cs_sla_matrix_t  *a,
                          const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  double  val;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity checks */
  assert(a->type == CS_SLA_MAT_MSR && b->type == CS_SLA_MAT_DEC);

  /* Initialize structures */
  _estimate_sizes(a, b, &size_max, &lst_size_init);
  c = _init_mat(a, b, CS_SLA_MAT_CSR, size_max);
  spa = _spa_init(b->n_cols, lst_size_init);

  /* Compute c = a*b line by line */
  for (ii = 0; ii < a->n_rows; ii++) {

    /* Diagonal term */
    val = a->diag[ii];
    for (k = b->idx[ii]; k < b->idx[ii+1]; k++)
      kk = b->col_id[k], _spa_add(spa, val*b->sgn[k], kk, ii);

    /* Extra-diagonal terms */
    for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {

      jj = a->col_id[j], val = a->val[j];
      for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
        kk = b->col_id[k], _spa_add(spa, val*b->sgn[k], kk, ii);

    } /* End of loop on non-empty columns of row ii of a */

    /* Fill a new row in c and prepare next step */
    if (spa->nnz + idx > size_max)
      _resize_mat(c, spa->nnz + idx, &size_max);
    n_entries = _spa_gather(spa, idx, c->col_id, c->val);
    c->idx[ii+1] = idx + n_entries;
    idx = c->idx[ii+1];

  } /* End of loop on row ii of a */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the product C = At * Diag * A
 *          where A and At are DEC matrices
 *
 * \param[in]     At   pointer to a cs_sla_matrix_t struct. (DEC type)
 * \param[in]     D    array standing for a diagonal operator
 * \param[in]     A    pointer to a cs_sla_matrix_t struct. (DEC type)
 * \param[inout]  C    pointer to a cs_sla_matrix_t storing the result
 * \param[inout]  w    work buffer
 */
/*----------------------------------------------------------------------------*/

static void
_decdec_AtDA(const cs_sla_matrix_t  *At,
             const double            D[],
             const cs_sla_matrix_t  *A,
             cs_sla_matrix_t        *C,
             int                    *w)
{
  int  ii, j, jj, kk, k, shift;
  double  val;

  int  size = 0, size_max = At->n_rows;

  /* Sanity check */
  assert(A->type == CS_SLA_MAT_DEC);
  assert(At->type == CS_SLA_MAT_DEC);

  /* Default allocation */
  BFT_MALLOC(C->col_id, size_max, cs_lnum_t);
  BFT_MALLOC(C->val, size_max, double);

  /* Compute the double multiplication */
  for (ii = 0; ii < At->n_rows; ii++) {
    for (j = At->idx[ii]; j < At->idx[ii+1]; j++) {

      jj = At->col_id[j];
      val = D[jj-1] * At->sgn[j];

      for (k = A->idx[jj]; k < A->idx[jj+1]; k++) {

        kk = A->col_id[k];
        shift = w[kk];

        if (shift == -1) { /* Add a new entry */
          if (size + 1 > size_max) { /* Realloc buffers */
            size_max *= 1.5;
            assert(size_max > size);
            BFT_REALLOC(C->col_id, size_max, cs_lnum_t);
            BFT_REALLOC(C->val, size_max, double);
          }

          w[kk] = size;
          C->col_id[size] = kk;
          C->val[size] = val * A->sgn[k];
          size++;

        }
        else /* Entry is already defined */
          C->val[shift] += val * A->sgn[k];

      } /* End of loop on non-empty coloms of row jj of A */
    } /* End of loop on non-empty coloms of row ii of At */

    C->idx[ii+1] = size;
    shift = C->idx[ii];

    for (k = shift; k < C->idx[ii+1]; k++) {
      w[C->col_id[k]] = -1;
      if (fabs(C->val[k]) > cs_base_zthreshold) { /* Clean zero entry */
        if (k != shift) {
          C->col_id[shift] = C->col_id[k];
          C->val[shift] = C->val[k];
        }
        shift++;
      }
    } /* End of loop on new elements */

    C->idx[ii+1] = shift;
    size = shift;

  } /* End of loop on row ii of a */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the product C = At * Diag * A
 *          where A and At are CSR matrices
 *
 * \param[in]     At   pointer to a cs_sla_matrix_t struct. (CSR type)
 * \param[in]     D    array standing for a diagonal operator
 * \param[in]     A    pointer to a cs_sla_matrix_t struct. (CSR type)
 * \param[inout]  C    pointer to a cs_sla_matrix_t storing the result
 * \param[inout]  w    work buffer
 */
/*----------------------------------------------------------------------------*/

static void
_csrcsr_AtDA(const cs_sla_matrix_t  *At,
             const double            D[],
             const cs_sla_matrix_t  *A,
             cs_sla_matrix_t        *C,
             int                    *w)
{
  int  ii, j, jj, kk, k, shift;
  double  val;

  int  size = 0, size_max = A->n_rows;

  /* Sanity check */
  assert(A->type == CS_SLA_MAT_CSR);
  assert(At->type == CS_SLA_MAT_CSR);

  /* Default allocation */
  BFT_MALLOC(C->col_id, size_max, int);
  BFT_MALLOC(C->val, size_max, double);

  /* Compute the double multiplication */
  for (ii = 0; ii < At->n_rows; ii++) {
    for (j = At->idx[ii]; j < At->idx[ii+1]; j++) {

      jj = At->col_id[j];
      val = D[jj] * At->val[j];

      for (k = A->idx[jj]; k < A->idx[jj+1]; k++) {

        kk = A->col_id[k];
        shift = w[kk];

        if (shift == -1) { /* Add a new entry */
          if (size + 1 > size_max) { /* Realloc buffers */
            size_max *= 1.5;
            assert(size_max > size);
            BFT_REALLOC(C->col_id, size_max, int);
            BFT_REALLOC(C->val, size_max, double);
          }

          w[kk] = size;
          C->col_id[size] = kk;
          C->val[size] = val * A->val[k];
          size++;

        }
        else /* Entry is already defined */
          C->val[shift] += val * A->val[k];

      } /* End of loop on non-empty coloms of row jj of A */
    } /* End of loop on non-empty coloms of row ii of At */

    C->idx[ii+1] = size;
    shift = C->idx[ii];

    for (k = shift; k < C->idx[ii+1]; k++) {
      w[C->col_id[k]] = -1;
      if (fabs(C->val[k]) > cs_base_zthreshold) { /* Clean zero entry */
        if (k != shift) {
          C->col_id[shift] = C->col_id[k];
          C->val[shift] = C->val[k];
        }
        shift++;
      }
    } /* End of loop on new elements */

    C->idx[ii+1] = shift;
    size = shift;

  } /* End of loop on row ii of a */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Fill a matrix resulting from an extraction of row/column for
 *          matrix in DEC format
 *
 * \param[inout]  final        pointer to the final matrix struct.
 * \param[in]     init         init matrix to work with
 * \param[in]     row_z2i_ids  zipped -> init numbering for rows
 * \param[in]     col_i2z_ids  init   -> zipped numbering for columns
 *
 * \return a pointer to the new final matrix struct
 */
/*----------------------------------------------------------------------------*/

static void
_pack_dec(cs_sla_matrix_t         *final,
          const cs_sla_matrix_t   *init,
          const cs_lnum_t         *row_z2i_ids,
          const cs_lnum_t         *col_i2z_ids)
{
  int  j, shift, zrow_id, irow_id, icol_id;

  /* Fill entries */
  if (final->n_cols == init->n_cols) { /* No extraction of columns */

    for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
      irow_id = row_z2i_ids[zrow_id];
      shift = final->idx[zrow_id];
      for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++, shift++) {
        final->col_id[shift] = init->col_id[j];
        final->sgn[shift] = init->sgn[j];
      }
    }

  }
  else {

    for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
      irow_id = row_z2i_ids[zrow_id];
      shift = final->idx[zrow_id];
      for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++) {
        icol_id = init->col_id[j];
        if (col_i2z_ids[icol_id] > -1) { /* Not removed */
          final->col_id[shift] = col_i2z_ids[icol_id];
          final->sgn[shift] = init->sgn[j];
          shift++;
        }
      } /* Loop on columns in init */

    } /* Loop on rows in final */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Fill a matrix resulting from an extraction of row/column for
 *          matrix in CSR format
 *
 * \param[inout]  final        pointer to the block matrix struct.
 * \param[in]     init         init matrix to work with
 * \param[in]     row_z2i_ids  zipped -> init numbering for rows
 * \param[in]     col_i2z_ids  init -> zipped numbering for columns
 *
 * \return a pointer to the new final matrix struct
 */
/*----------------------------------------------------------------------------*/

static void
_pack_csr(cs_sla_matrix_t         *final,
          const cs_sla_matrix_t   *init,
          const cs_lnum_t         *row_z2i_ids,
          const cs_lnum_t         *col_i2z_ids)
{
  int  j, shift, zrow_id, irow_id, icol_id;

  /* Fill entries */
  if (final->n_cols == init->n_cols) { /* No extraction of columns */

    for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
      irow_id = row_z2i_ids[zrow_id];
      shift = final->idx[zrow_id];
      for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++, shift++) {
        final->col_id[shift] = init->col_id[j];
        final->val[shift] = init->val[j];
      }
    }

  }
  else {

    for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
      irow_id = row_z2i_ids[zrow_id];
      shift = final->idx[zrow_id];
      for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++) {
        icol_id = init->col_id[j];
        if (col_i2z_ids[icol_id] > -1) { /* Not removed */
          final->col_id[shift] = col_i2z_ids[icol_id];
          final->val[shift] = init->val[j];
          shift++;
        }
      } /* Loop on columns in init */

    } /* Loop on rows in final */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Fill a matrix resulting from an extraction of row/column for
 *          matrix in MSR format
 *
 * \param[inout]  final        pointer to the block matrix struct.
 * \param[in]     init         init matrix to work with
 * \param[in]     row_z2i_ids  zipped -> init numbering for rows
 * \param[in]     col_i2z_ids  init -> zipped numbering for columns
 * \param[in]     msr2csr      change matrix type
 *
 * \return a pointer to the new final matrix struct
 */
/*----------------------------------------------------------------------------*/

static void
_pack_msr(cs_sla_matrix_t         *final,
          const cs_sla_matrix_t   *init,
          const cs_lnum_t         *row_z2i_ids,
          const cs_lnum_t         *col_i2z_ids,
          _Bool                    msr2csr)
{
  int  j, shift, zrow_id, irow_id, icol_id;

  /* Fill entries */
  if (final->n_cols == init->n_cols) { /* No extraction of columns */

    if (msr2csr) {

      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        shift = final->idx[zrow_id];
        irow_id = row_z2i_ids[zrow_id];
        if (zrow_id < final->n_cols) { // Diag. entry
          final->col_id[shift] = zrow_id;
          final->val[shift] = init->diag[irow_id];
          shift++;
        }
        for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++, shift++) {
          final->col_id[shift] = init->col_id[j];
          final->val[shift] = init->val[j];
        }
      }

    }
    else { /* msr2csr = false */

      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        final->diag[zrow_id] = init->diag[irow_id];
        shift = final->idx[zrow_id];
        for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++, shift++) {
          final->col_id[shift] = init->col_id[j];
          final->val[shift] = init->val[j];
        }
      }

    } /* msr2csr */

  }
  else { /* Column and row extraction */

    if (msr2csr) {

      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        shift = final->idx[zrow_id];
        if (zrow_id < final->n_cols) {
          final->col_id[shift] = zrow_id;
          final->val[shift] = init->diag[irow_id];
          shift++;
        }
        for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++) {
          icol_id = init->col_id[j];
          if (col_i2z_ids[icol_id] > -1) { /* Not removed */
            final->col_id[shift] = col_i2z_ids[icol_id];
            final->val[shift] = init->val[j];
            shift++;
          }
        } /* Loop on columns in init */
      } /* Loop on rows in final */

    }
    else { /* msr2csr = false */

      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        final->diag[zrow_id] = init->diag[irow_id];
        shift = final->idx[zrow_id];
        for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++) {
          icol_id = init->col_id[j];
          if (col_i2z_ids[icol_id] > -1) { /* Not removed */
            final->col_id[shift] = col_i2z_ids[icol_id];
            final->val[shift] = init->val[j];
            shift++;
          }
        } /* Loop on columns in init */
      } /* Loop on rows in final */

    } /* test msr2csr */

  } /* test final->n_cols */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute a matrix norm
 *
 * Definition of these norms comes from the book "TEMPLATES for the solution
 * of linear systems" (R. Barret, M. Berry, T. Chan, J. Demmel, J. Donato,
 * J. Dongarra, V. Eijkhout, R. Polzo, C. Romine and H.K. Van der Vorst)
 *
 * CS_SLA_NORM_APPROX is not a true norm but rather a way of simply computing
 * an order of magnitude of what could be a norm of the matrix. (p.52 of the
 * same book)
 *----------------------------------------------------------------------------*/

double
cs_sla_get_matrix_norm(const cs_sla_matrix_t  *m)
{
  int  i, nnz;

  double  m_norm = 0.0;

  if (m == NULL)
    return m_norm;

  if (m->type == CS_SLA_MAT_NONE)
    return m_norm;

  assert(m->stride == 1);

  nnz = m->idx[m->n_rows];

  if (m->type == CS_SLA_MAT_CSR || m->type == CS_SLA_MAT_MSR)
    m_norm = cs_dp(nnz, m->val, m->val);
  else if (m->type == CS_SLA_MAT_DEC)
    for (i = 0; i < nnz; i++)
      m_norm += m->sgn[i]*m->sgn[i];

  assert(m_norm >= 0);
  m_norm = sqrt(m_norm);

  return m_norm;
}

/*----------------------------------------------------------------------------
 * Compute if necessary the spectral radius related to a matrix with a given
 * precision
 *----------------------------------------------------------------------------*/

int
cs_sla_get_mat_spectral_radius(const cs_sla_matrix_t  *matrix,
                               int                     iter_max,
                               double                  epsilon,
                               double                 *lambda)
{
  int  info;

  cs_sla_mat_prop_t  *mp = NULL;

  assert(matrix != NULL);
  assert(matrix->stride == 1);

  mp = matrix->properties;

  if (epsilon >= mp->eps || iter_max <= mp->iter_count)        {
    (*lambda) = mp->spectral_radius;
    info = -1;
    return info;
  }
  else        {

    int  ii, iter;
    double  test;

    int  n_rows = matrix->n_rows;
    double  *_aux_vectors = NULL, *x0 = NULL, *yk = NULL, *xk = NULL;

    BFT_MALLOC(_aux_vectors, n_rows*3, double);

    x0 = _aux_vectors;
    yk = _aux_vectors + n_rows;
    xk = _aux_vectors + 2*n_rows;

    for (ii = 1; ii < n_rows; ii++) x0[ii] = 0;

    x0[0] = 1;
    iter = 1;

    do        {
      cs_sla_matvec(matrix,x0, &yk, true);
      (*lambda) = sqrt(cs_dp(n_rows, yk, yk));

      for (ii = 0; ii<n_rows; ii++) {
        xk[ii] = yk[ii] / (*lambda);
        yk[ii] = xk[ii] - x0[ii];
        x0[ii] = xk[ii];
      }

      test = sqrt(cs_dp(n_rows, yk, yk));
      iter++;

    } while (test > epsilon && iter < iter_max);

    if (iter >= iter_max)
      info = 2;
    else
      info = 0;

    mp->spectral_radius = (*lambda);
    mp->iter_count = iter;
    mp->eps = test;

    BFT_FREE(_aux_vectors);

  } /* End of test */

  return info;
}

/*----------------------------------------------------------------------------
 * Create a new matrix from a block description (not done for all matrix
 * combinations, see below).
 *
 *     |A | B|
 * M = |--|--|
 *     |C | D|
 *----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_block2mat(cs_sla_matrix_t        *A,
                 const cs_sla_matrix_t  *B,
                 const cs_sla_matrix_t  *C,
                 const cs_sla_matrix_t  *D,
                 _Bool                   sym)
{
  int  i, j, n_rows, n_cols, nnz;

  int  shift = 0;
  _Bool  msr2csr = false;
  cs_sla_matrix_t  *M = NULL;

 /* Sanity check. D can be NULL */
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  /* Up to now it is the only case which is dealt */
  assert(A->type == CS_SLA_MAT_CSR || A->type == CS_SLA_MAT_MSR);
  assert(B->type == CS_SLA_MAT_DEC);
  assert(C->type == CS_SLA_MAT_DEC);
  assert(A->stride == 1);
  assert(B->stride == 1);
  assert(C->stride == 1);

  if (D != NULL) {
    assert(D->type == CS_SLA_MAT_CSR);
    assert(D->stride == 1);
  }

  if (A->type == CS_SLA_MAT_MSR) {
    msr2csr = true;
    cs_sla_matrix_msr2csr(A);
  }

  n_rows = A->n_rows + C->n_rows;
  n_cols = A->n_cols + B->n_cols;
  M = cs_sla_matrix_create(n_rows, n_cols, 1, CS_SLA_MAT_CSR, sym);

  /* Build index */
  for (i = 0; i < A->n_rows; i++)
    M->idx[i+1] =  M->idx[i] + A->idx[i+1] - A->idx[i]
      + B->idx[i+1] - B->idx[i];

  if (D != NULL) {
    for (j = A->n_rows, i = 0; i < C->n_rows; i++, j++)
      M->idx[j+1] =  M->idx[j] + C->idx[i+1] - C->idx[i]
        + D->idx[i+1] - D->idx[i];
  }
  else { /* D == NULL */
    for (j = A->n_rows, i = 0; i < C->n_rows; i++, j++)
      M->idx[j+1] = M->idx[j] + C->idx[i+1] - C->idx[i];
  }

  nnz = M->idx[M->n_rows];
  BFT_MALLOC(M->col_id, nnz, int);
  BFT_MALLOC(M->val, nnz, double);

  /* Build global matrix */
  for (i = 0; i < A->n_rows; i++) {

    shift = M->idx[i];
    for (j = A->idx[i]; j < A->idx[i+1]; j++) {
      M->col_id[shift] = A->col_id[j];
      M->val[shift] = A->val[j];
      shift++;
    }

    for (j = B->idx[i]; j < B->idx[i+1]; j++) {
      M->col_id[shift] = B->col_id[j] + A->n_cols;
      M->val[shift] = (double)B->sgn[j]; /* assume it's a DEC type */
      shift++;
    }

  }

  for (i = 0; i < C->n_rows; i++) {

    shift = M->idx[A->n_rows + i];
    for (j = C->idx[i]; j < C->idx[i+1]; j++) {
      M->col_id[shift] = C->col_id[j];
      M->val[shift] = (double)C->sgn[j];
      shift++;
    }

    if (D != NULL) {
      for (j = D->idx[i]; j < D->idx[i+1]; j++) {
        M->col_id[shift] = D->col_id[j] + C->n_cols;
        M->val[shift] = D->val[j];
        shift++;
      }
    }

  }

  assert(shift == nnz);

  if (msr2csr)
    cs_sla_matrix_csr2msr(A);

  return M;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two sparse matrices a and b. c = alpha*a + beta*b
 *
 * \param[in]  alpha   first coef.
 * \param[in]  a       first matrix to add
 * \param[in]  beta    second coef.
 * \param[in]  b       second matrix to add
 *
 * \return  a pointer to the resulting matrix
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_add(double                  alpha,
                  const cs_sla_matrix_t  *a,
                  double                  beta,
                  const cs_sla_matrix_t  *b)
{
  int  ii, j, jj;
  size_t  annz, bnnz;

  cs_sla_matrix_type_t  c_type = CS_SLA_MAT_CSR;
  size_t  idx = 0, n_entries = 0, size_max = 0, lst_max = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  const char errmsg[] =
    "  Incompatible matrix type.\n  This combination is not available yet.\n"
    "  Matrices addition is aborted.\n";

  /* Sanity check */
  assert(a != NULL && b != NULL);
  assert(a->n_rows == b->n_rows);
  assert(a->n_cols == b->n_cols);

  if (a->stride > 1 || b->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  /* Estimate sizes and allocate structures */
  annz = a->idx[a->n_rows];
  bnnz = b->idx[b->n_rows];
  size_max = CS_MAX((annz+bnnz)/2, (size_t)a->n_rows);

  if (a->type == CS_SLA_MAT_MSR && b->type == CS_SLA_MAT_MSR)
    c_type = CS_SLA_MAT_MSR;
  c = _init_mat(a, b, c_type, size_max);

  lst_max = CS_MAX(5, (annz + bnnz)/(2*a->n_rows));
  spa = _spa_init(b->n_cols, lst_max);

  if (a->type == b->type) {

    if (a->type == CS_SLA_MAT_DEC) {  /* DEC-DEC */

      /* Compute c = alpha*a + beta*b row by row */
      for (ii = 0; ii < a->n_rows; ii++) {

        /* Contribution of matrix a */
        for (j = a->idx[ii]; j < a->idx[ii+1]; j++)
          jj = a->col_id[j], _spa_add(spa, alpha*a->sgn[j], jj, ii);

        /* Contribution of matrix b */
        for (j = b->idx[ii]; j < b->idx[ii+1]; j++)
          jj = b->col_id[j], _spa_add(spa, beta*b->sgn[j], jj, ii);

        /* Fill a new row in c and prepare next step */
        if (spa->nnz + idx > size_max)
          _resize_mat(c, spa->nnz + idx, &size_max);
        n_entries = _spa_gather(spa, idx, c->col_id, c->val);
        c->idx[ii+1] = idx + n_entries;
        idx = c->idx[ii+1];

      } /* Loop on non-empty entry of row  */

    } /* a = b = DEC */
    else if (a->type == CS_SLA_MAT_CSR || a->type == CS_SLA_MAT_MSR) {

      if (a->type == CS_SLA_MAT_MSR)  /* Diagonal entries */
        for (ii = 0; ii < a->n_rows; ii++)
          c->diag[ii] = alpha*a->diag[ii] + beta*b->diag[ii];

      /* Compute c = alpha*a + beta*b row by row */
      for (ii = 0; ii < a->n_rows; ii++) {

        /* Contribution of matrix a */
        for (j = a->idx[ii]; j < a->idx[ii+1]; j++)
          jj = a->col_id[j], _spa_add(spa, alpha*a->val[j], jj, ii);

        /* Contribution of matrix b */
        for (j = b->idx[ii]; j < b->idx[ii+1]; j++)
          jj = b->col_id[j], _spa_add(spa, beta*b->val[j], jj, ii);

        /* Fill a new row in c and prepare next step */
        if (spa->nnz + idx > size_max)
          _resize_mat(c, spa->nnz + idx, &size_max);
        n_entries = _spa_gather(spa, idx, c->col_id, c->val);
        c->idx[ii+1] = idx + n_entries;
        idx = c->idx[ii+1];

      } /* Loop on non-empty entry of row  */

    }
    else
      bft_error(__FILE__, __LINE__, 0, errmsg);

  }
  else { /* a->type != b->type */

    if (a->type == CS_SLA_MAT_MSR && b->type == CS_SLA_MAT_CSR) {

      /* Compute c = alpha*a + beta*b row by row */
      for (ii = 0; ii < a->n_rows; ii++) {

        /* Contribution of matrix a */
        _spa_add(spa, alpha*a->diag[ii], ii, ii);
        for (j = a->idx[ii]; j < a->idx[ii+1]; j++)
          jj = a->col_id[j], _spa_add(spa, alpha*a->val[j], jj, ii);

        /* Contribution of matrix b */
        for (j = b->idx[ii]; j < b->idx[ii+1]; j++)
          jj = b->col_id[j], _spa_add(spa, beta*b->val[j], jj, ii);

        /* Fill a new row in c and prepare next step */
        if (spa->nnz + idx > size_max)
          _resize_mat(c, spa->nnz + idx, &size_max);
        n_entries = _spa_gather(spa, idx, c->col_id, c->val);
        c->idx[ii+1] = idx + n_entries;
        idx = c->idx[ii+1];

      } /* Loop on non-empty entry of row  */

    }
    else
      bft_error(__FILE__, __LINE__, 0, errmsg);

  }

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main subroutine to multiply two sparse matrices a and b. c= a*b
 *
 * \param[in]  a     first matrix to multiply
 * \param[in]  b     second matrix to multiply
 *
 * \return  a pointer to the resulting matrix
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_multiply(const cs_sla_matrix_t   *a,
                       const cs_sla_matrix_t   *b)
{
  cs_sla_matrix_t  *c = NULL;

  const char errmsg[] =
    "  Incompatible matrix type.\n  This combination is not available yet.\n"
    "  Matrices multiplication is aborted.\n";

  /* Sanity check */
  assert(a != NULL && b != NULL);
  assert(a->n_cols == b->n_rows);

  if (a->type == CS_SLA_MAT_NONE || b->type == CS_SLA_MAT_NONE)
    bft_error(__FILE__, __LINE__, 0,
              "  Incompatible matrix type.\n   Cannot mulitply matrices\n");
  if (a->stride > 1 || b->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  /* Switch to the right algorithm according to the kind of matrices */
  if (a->type == CS_SLA_MAT_DEC) {

    if (b->type == CS_SLA_MAT_DEC)        /* DEC-DEC */
      c = _multiply_dec_matrices(a, b);

    else if (b->type == CS_SLA_MAT_MSR)   /* DEC-MSR */
      c = _multiply_decmsr_matrices(a, b);

    else if(b->type == CS_SLA_MAT_CSR)    /* DEC-CSR */
      c = _multiply_deccsr_matrices(a, b);

    else
      bft_error(__FILE__, __LINE__, 0, errmsg);

  }
  else if (a->type == CS_SLA_MAT_CSR) {

    if (b->type == CS_SLA_MAT_CSR)       /* CSR-CSR */
      c = _multiply_csr_matrices(a, b);

    else if (b->type == CS_SLA_MAT_DEC)  /* CSR-DEC */
      c = _multiply_csrdec_matrices(a, b);

    else
      bft_error(__FILE__, __LINE__, 0, errmsg);

  }
  else if (a->type == CS_SLA_MAT_MSR) {

    if (b->type == CS_SLA_MAT_DEC)       /* MSR-DEC */
      c = _multiply_msrdec_matrices(a, b);

    else
      bft_error(__FILE__, __LINE__, 0, errmsg);

  } /* End of switch on a->type */
  else
    bft_error(__FILE__, __LINE__, 0, errmsg);

  /* Build diagonal index */
  cs_sla_matrix_diag_idx(c);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Specific matrix multiplication. Compute Bt * beta * B + alpha * A
 *          alpha and beta are scalar
 *
 * \param[in] alpha     real coefficient
 * \param[in] a         square sym. matrix
 * \param[in] beta      real coefficient
 * \param[in] b         matrix (CSR or DEC)
 * \param[in] bt        adjoint matrix of b
 *
 * \return the matrix corresponding to this operation in MSR storage
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_combine(double                  alpha,
                      const cs_sla_matrix_t  *a,
                      double                  beta,
                      const cs_sla_matrix_t  *bt,
                      const cs_sla_matrix_t  *b)
{
  int  ii, j, jj, k, kk;
  double  val, contrib;

  size_t  idx = 0, n_entries = 0, size_max = 0, lst_size_init = 0;
  _spa_t  *spa = NULL;
  cs_sla_matrix_t  *c = NULL;

  /* Sanity checks */
  assert(b != NULL && bt != NULL && a != NULL);
  assert(bt->n_cols == b->n_rows);
  assert(bt->n_rows == a->n_rows);
  assert(a->type == CS_SLA_MAT_MSR || a->type == CS_SLA_MAT_CSR);
  assert(bt->type == b->type);
  assert(bt->type != CS_SLA_MAT_MSR);

  if (b->stride > 1 || bt->stride > 1 || a->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  /* Initialize matrix c = alpha*a + beta*bt*b and structures */
  size_max = a->idx[a->n_rows];
  c = _init_mat(a, a, CS_SLA_MAT_CSR, size_max);
  lst_size_init = CS_MAX(1, a->idx[a->n_rows]/a->n_rows);
  spa = _spa_init(a->n_rows, lst_size_init);

  if (c->type == CS_SLA_MAT_MSR) { /* Define diagonal entries */

    if (bt->type == CS_SLA_MAT_DEC) {

      for (ii = 0; ii < a->n_rows; ii++) {
        contrib = 0.0;
        for (j = bt->idx[ii]; j < bt->idx[ii+1]; j++)
          val = bt->sgn[j], contrib += val*val;
        c->diag[ii] = alpha*a->diag[ii] + beta*contrib;
      }

    } /* DEC */
    else {
      assert(bt->type == CS_SLA_MAT_CSR);

      for (ii = 0; ii < a->n_rows; ii++) {
        contrib = 0.0;
        for (j = bt->idx[ii]; j < bt->idx[ii+1]; j++)
          val = bt->val[j], contrib += val*val;
        c->diag[ii] = alpha*a->diag[ii] + beta*contrib;
      }

    }

  } /* MSR type: define diagonal */

  /* Compute the remaining part line by line */
  if (bt->type == CS_SLA_MAT_DEC) {

    for (ii = 0; ii < a->n_rows; ii++) {

      /* Add contribution from alpha*a */
      for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {
        jj = a->col_id[j], val = alpha*a->val[j];
        _spa_add(spa, val, jj, ii);
      }
      if (a->type == CS_SLA_MAT_MSR)
        _spa_add(spa, alpha*a->diag[ii], ii, ii);

      for (j = bt->idx[ii]; j < bt->idx[ii+1]; j++) {

        jj = bt->col_id[j], val = beta*bt->sgn[j];

        for (k = b->idx[jj]; k < b->idx[jj+1]; k++) {
          kk = b->col_id[k];
          _spa_add(spa, val*b->sgn[k], kk, ii);
        }

      } /* End of loop on non-empty columns of row ii of bt */

      /* Fill a new row in c and prepare next step */
      if (spa->nnz + idx > size_max)
        _resize_mat(c, spa->nnz + idx, &size_max);
      n_entries = _spa_gather(spa, idx, c->col_id, c->val);
      c->idx[ii+1] = idx + n_entries;
      idx = c->idx[ii+1];

    } /* End of loop on rows */

  } /* DEC type */
  else  {

    for (ii = 0; ii < a->n_rows; ii++) {

      /* Add contribution from alpha*a */
      for (j = a->idx[ii]; j < a->idx[ii+1]; j++) {
        jj = a->col_id[j], val = alpha*a->val[j];
        _spa_add(spa, val, jj, ii);
      }
      if (a->type == CS_SLA_MAT_MSR)
        _spa_add(spa, alpha*a->diag[ii], ii, ii);

      for (j = bt->idx[ii]; j < bt->idx[ii+1]; j++) {

        jj = bt->col_id[j], val = beta * bt->val[j];

        for (k = b->idx[jj]; k < b->idx[jj+1]; k++)
          kk = b->col_id[k], _spa_add(spa, val*b->val[k], kk, ii);

      } /* End of loop on non-empty columns of row ii of bt */

      /* Fill a new row in c and prepare next step */
      if (spa->nnz + idx > size_max)
        _resize_mat(c, spa->nnz + idx, &size_max);
      n_entries = _spa_gather(spa, idx, c->col_id, c->val);
      c->idx[ii+1] = idx + n_entries;
      idx = c->idx[ii+1];

    } /* End of loop on rows */

  } /* CSR type */

  /* Memory management */
  BFT_REALLOC(c->col_id, c->idx[c->n_rows], cs_lnum_t);
  BFT_REALLOC(c->val, c->idx[c->n_rows], double);
  spa = _spa_free(spa);

  /* Build diagonal index */
  cs_sla_matrix_diag_idx(c);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the product C = At * Diag * A
 *
 * \param[in]     At   pointer to a cs_sla_matrix_t struct. (DEC type)
 * \param[in]     D    array standing for a diagonal operator
 * \param[in]     A    pointer to a cs_sla_matrix_t struct. (DEC type)
 * \param[inout]  w    work buffer
 *
 * \return a pointer to a new cs_sla_matrix_t
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_multiply_AtDA(const cs_sla_matrix_t  *At,
                     const double            D[],
                     const cs_sla_matrix_t  *A,
                     int                    *w)
{
  cs_lnum_t  i;

  int  use_work = 0;
  cs_sla_matrix_t  *C = NULL;

  /* Sanity check */
  assert(A != NULL && At != NULL && D != NULL);
  assert(At->n_cols == A->n_rows);

  if (A->stride > 1 || At->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  C = cs_sla_matrix_create(At->n_rows, A->n_cols, 1,
                           CS_SLA_MAT_CSR,
                           true); /* Symmetric by construction */

  /* Manage temporary buffer */
  if (w == NULL) {
    BFT_MALLOC(w, A->n_cols, int);
    use_work = 1;
  }
  for (i = 0; i < A->n_cols; i++)
    w[i] = -1;

  if (A->type == CS_SLA_MAT_CSR && At->type == CS_SLA_MAT_CSR)
    _csrcsr_AtDA(At, D, A, C, w);
  else if (A->type == CS_SLA_MAT_DEC && At->type == CS_SLA_MAT_DEC)
    _decdec_AtDA(At, D, A, C, w);

  /* Memory management */
  BFT_REALLOC(C->col_id, C->idx[C->n_rows], int);
  BFT_REALLOC(C->val, C->idx[C->n_rows], double);

  if (use_work == 1)
    BFT_FREE(w);

  /* Build diagonal index */
  cs_sla_matrix_diag_idx(C);

  return C;
}

/*----------------------------------------------------------------------------*/
/*                  Matrix by vector multiplication                           */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 * Compute matrix vector product.
 *   If inout is not allocated, allocation is done and initialization to 0
 *   If inout is allocated, add values to the given vector
 *----------------------------------------------------------------------------*/

void
cs_sla_matvec(const cs_sla_matrix_t  *m,
              const double            v[],
              double                 *inout[],
              _Bool                   reset)
{
  int  i, j;
  double  sum;

  double  *out = *inout;

  assert(m != NULL);

  if (m->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  if (out == NULL) {
    BFT_MALLOC(out, m->n_rows, double);
    reset = true;
  }

  if (reset == true)
    for (i = 0; i < m->n_rows; i++)
      out[i] = 0.0;

  switch (m->type) {
  case CS_SLA_MAT_DEC:
    for (i = 0; i < m->n_rows; i++) {
      sum = 0.0;
      for (j = m->idx[i]; j < m->idx[i+1]; j++)
        sum += m->sgn[j]*v[m->col_id[j]];
      out[i] += sum;
    }
    break;

  case CS_SLA_MAT_CSR:
    for (i = 0; i < m->n_rows; i++) {
      sum = 0.0;
      for (j = m->idx[i]; j < m->idx[i+1]; j++)
        sum += m->val[j] * v[m->col_id[j]];
      out[i] += sum;
    }
    break;

  case CS_SLA_MAT_MSR:
    assert(m->n_rows == m->n_cols); /* Should be a squared matrix */
    for (i = 0; i < m->n_rows; i++) {
      sum = m->diag[i] * v[i];
      for (j = m->idx[i]; j < m->idx[i+1]; j++)
        sum += m->val[j] * v[m->col_id[j]];
      out[i] += sum;
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "  Incompatible matrix type.\n"
              "  Cannot mulitply matrix by vector.\n");
    break;

  } // end of switch

  /* Return pointer */
  *inout = out;
}

/*----------------------------------------------------------------------------*
 * Compute out = alpha * Mx + beta * y
 *----------------------------------------------------------------------------*/

void
cs_sla_amxby(double                  alpha,
             const cs_sla_matrix_t  *m,
             const double            x[],
             double                  beta,
             const double            y[],
             double                 *inout[])
{
  int  i;

  double  *out = *inout;

  cs_sla_matvec(m, x, &out, true);

  for (i = 0; i < m->n_rows; i++) {
    out[i] *= alpha;
    out[i] += beta*y[i];
  }

  *inout = out;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Matrix block 2x2 multiply by a vector
 *
 *              | A  B | |X| = |F|= |AX + BY|
 *              | C  D | |Y|   |G|  |CX + DY|
 *
 * \param[in]    A     pointer to a cs_sla_matrix_t block (1,1)
 * \param[in]    B     pointer to a cs_sla_matrix_t block (1,2)
 * \param[in]    C     pointer to a cs_sla_matrix_t block (2,1)
 * \param[in]    D     pointer to a cs_sla_matrix_t block (2,2)
 * \param[in]    X     upper vector
 * \param[in]    Y     lower vector
 * \param[inout] F     upper result
 * \param[inout] G     lower result
 * \param[in]    reset reset before computation (true/false)
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matvec_block2(const cs_sla_matrix_t  *A,
                     const cs_sla_matrix_t  *B,
                     const cs_sla_matrix_t  *C,
                     const cs_sla_matrix_t  *D,
                     const double            X[],
                     const double            Y[],
                     double                 *F[],
                     double                 *G[],
                     _Bool                   reset)
{
  int  i;

  double  *_F = *F, *_G = *G;
  int  nx = 0, ny = 0;

  assert(A != NULL || B != NULL);
  assert(C != NULL || D != NULL);
  assert(X != NULL);
  assert(Y != NULL);

  if (A->stride > 1 || B->stride > 1 || C->stride > 1 || D->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  if (A != NULL)
    nx = A->n_rows;
  else
    nx = B->n_rows;
  if (C != NULL)
    ny = B->n_rows;
  else
    ny = D->n_rows;
  if (_F == NULL) {
    BFT_MALLOC(_F, nx, double);
    reset = true;
  }
  if (_G == NULL) {
    BFT_MALLOC(_G, ny, double);
    reset = true;
  }

  if (reset == true) {
    for (i = 0; i < nx; i++) _F[i] = 0.0;
    for (i = 0; i < ny; i++) _G[i] = 0.0;
  }

  /* Compute _F = A*X + B*Y */
  if (A != NULL)
    cs_sla_matvec(A, X, &_F, reset);
  if (B != NULL)
    cs_sla_matvec(B, Y, &_F, reset);

  /* Compute _G = A*X + B*Y */
  if (C != NULL)
    cs_sla_matvec(C, X, &_G, reset);
  if (D != NULL)
    cs_sla_matvec(D, Y, &_G, reset);

  /* Return pointers */
  *F = _F;
  *G = _G;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Change matrix representation from MSR to CSR
 *
 * \param[inout] a     matrix to transform
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_msr2csr(cs_sla_matrix_t  *a)
{
  int  i, j, s, e, new_nnz;

  cs_lnum_t  shift = 0;
  cs_lnum_t  *new_index = NULL, *new_col_id = NULL;
  double  *new_val = NULL;

  if (a->type == CS_SLA_MAT_CSR) /* Nothing to do */
    return;

  if (a->type != CS_SLA_MAT_MSR)
    bft_error(__FILE__, __LINE__, 0,
              "  Incompatible matrix type.\n"
              "  Cannot convert matrix from MSR -> CSR\n");

  if (a->stride > 1) /* Should not be useful */
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  /* Sanity checks */
  assert(a->diag != NULL);

  /* Allocate new buffers */
  BFT_MALLOC(new_index, a->n_rows + 1, cs_lnum_t);
  new_index[0] = 0;
  new_nnz = a->idx[a->n_rows] + a->n_rows;
  BFT_MALLOC(new_col_id, new_nnz, cs_lnum_t);
  BFT_MALLOC(new_val, new_nnz, double);

  /* Fill new buffers */
  for (i = 0; i < a->n_rows; i++) {

    s = a->idx[i], e = a->idx[i+1];

    /* Diagonal term */
    new_col_id[shift] = i;
    new_val[shift] = a->diag[i];
    shift++;

    /* Extra_diagonal terms */
    for (j = s; j < e; j++) {
      new_col_id[shift] = a->col_id[j];
      new_val[shift] = a->val[j];
      shift++;
    }
    new_index[i+1] = shift;

  } /* End of loop on rows */

  /* Change elements of the matrix structure */
  BFT_FREE(a->idx);
  BFT_FREE(a->col_id);
  BFT_FREE(a->val);
  BFT_FREE(a->diag);

  a->diag = NULL;
  a->idx = new_index;
  a->col_id = new_col_id;
  a->val = new_val;
  a->type = CS_SLA_MAT_CSR;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Change matrix representation from CSR to MSR
 *
 * \param[inout] a     matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_csr2msr(cs_sla_matrix_t  *a)
{
  int  row_id, j, s, e, new_nnz;

  cs_lnum_t  shift = 0;

  if (a->type == CS_SLA_MAT_MSR) /* Nothing to do */
    return;

  /* Sanity checks */
  assert(a->diag == NULL);
  assert(a->n_rows == a->n_cols); /* Only for squared matrices */

  if (a->type != CS_SLA_MAT_CSR)
    bft_error(__FILE__, __LINE__, 0,
              "  Incompatible matrix type.\n"
              "  Cannot convert matrix from CSR -> MSR\n");

  if (a->stride > 1) /* Should not be useful */
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  /* Allocate new buffers */
  BFT_MALLOC(a->diag, a->n_rows, double);
  for (row_id = 0; row_id < a->n_rows; row_id++)
    a->diag[row_id] = 0;

  /* Modify existing entries and fill diag */
  s = a->idx[0];
  for (row_id = 0; row_id < a->n_rows; row_id++) {

    /* Extra_diagonal terms */
    e = a->idx[row_id+1];
    for (j = s; j < e; j++) {

      if (a->col_id[j] == row_id)  /* Diagonal entry */
        a->diag[row_id] = a->val[j];
      else {                       /* Extra diag. entry */
        a->col_id[shift] = a->col_id[j];
        a->val[shift] = a->val[j];
        shift++;
      }

    }
    a->idx[row_id+1] = shift;
    s = e;

  } /* End of loop on rows */

  /* Change elements of the matrix structure */
  new_nnz = shift;
  BFT_REALLOC(a->col_id, new_nnz, cs_lnum_t);
  BFT_REALLOC(a->val, new_nnz, double);
  a->type = CS_SLA_MAT_MSR;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate its own pattern if shared
 *
 * \param[inout] a     matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_share2own(cs_sla_matrix_t  *a)
{
  int  *p = NULL;

  if (a == NULL)
    return;

  if (a->flag & CS_SLA_MATRIX_SHARED) { /* Pattern is shared */

    p = a->idx;
    BFT_MALLOC(a->idx, a->n_rows + 1, cs_lnum_t);
    memcpy(a->idx, p, (a->n_rows+1)*sizeof(cs_lnum_t));

    p = a->col_id;
    BFT_MALLOC(a->col_id, a->idx[a->n_rows], cs_lnum_t);
    memcpy(a->col_id, p, a->idx[a->n_rows]*sizeof(cs_lnum_t));

    if (a->didx != NULL) {
      p = a->didx;
      BFT_MALLOC(a->didx, a->n_rows, int);
      memcpy(a->didx, p, a->n_rows*sizeof(int));
    }

    a->flag = a->flag ^ CS_SLA_MATRIX_SHARED;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a new matrix resulting from the extraction of some listed
 *          rows and columns. The output is a new matrix packed (or zipped)
 *          w.r.t the initial matrix.
 *
 * \param[in]  n_final_rows   number of rows to extract
 * \param[in]  n_final_cols   number of columns to extract
 * \param[in]  init           init matrix to work with
 * \param[in]  row_z2i_ids    zipped -> initial ids for rows
 * \param[in]  col_i2z_ids    initial-> zipped ids for columns (-1 ==> remove)
 * \param[in]  keep_sym       true or false
 *
 * \return a pointer to the new (packed) matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_pack(cs_lnum_t                n_final_rows,
                   cs_lnum_t                n_final_cols,
                   const cs_sla_matrix_t   *init,
                   const cs_lnum_t         *row_z2i_ids,
                   const cs_lnum_t         *col_i2z_ids,
                   _Bool                    keep_sym)
{
  int  j, nnz, n_entries, zrow_id, irow_id;
  cs_sla_matrix_type_t  final_type;

  _Bool  msr2csr = false;
  cs_sla_matrix_t  *final = NULL;

  if (init == NULL)
    return final;

  /* Create a block matrix extracted from an initial (full) matrix */
  final_type = init->type;
  if (init->type == CS_SLA_MAT_MSR && n_final_rows != n_final_cols) {
    final_type = CS_SLA_MAT_CSR;
    msr2csr = true;
  }

  if (keep_sym && init->flag & CS_SLA_MATRIX_SYM)
    final = cs_sla_matrix_create(n_final_rows, n_final_cols, 1,
                                 final_type, true);
  else
    final = cs_sla_matrix_create(n_final_rows, n_final_cols, 1,
                                 final_type, false);

  if (init->type == CS_SLA_MAT_NONE)
    return final;

  /* Build index */
  if (!msr2csr) {

    if (n_final_cols == init->n_cols) { /* No size reduction for columns.
                                           Extract only rows */
      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        n_entries = init->idx[irow_id+1] - init->idx[irow_id];
        final->idx[zrow_id+1] = final->idx[zrow_id] + n_entries;
      }

    }
    else {

      /* Sanity check */
      assert(col_i2z_ids != NULL);

      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        n_entries = 0;
        for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++)
          if (col_i2z_ids[init->col_id[j]] > -1) /* Not removed */
            n_entries++;
        final->idx[zrow_id+1] = final->idx[zrow_id] + n_entries;

      } /* End of loop on final rows */

    } /* Columns reduction */

  }
  else { /* msr2csr = true */

    if (n_final_cols == init->n_cols) { /* No size reduction for columns.
                                           Extract only rows */
      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        n_entries = init->idx[irow_id+1] - init->idx[irow_id];
        if (zrow_id < n_final_cols) // Add diagonal term
          n_entries += 1;
        final->idx[zrow_id+1] = final->idx[zrow_id] + n_entries;
      }

    }
    else {

      /* Sanity check */
      assert(col_i2z_ids != NULL);

      for (zrow_id = 0; zrow_id < final->n_rows; zrow_id++) {
        irow_id = row_z2i_ids[zrow_id];
        n_entries = 0;
        if (zrow_id < n_final_cols)
          n_entries = 1;
        for (j = init->idx[irow_id]; j < init->idx[irow_id+1]; j++)
          if (col_i2z_ids[init->col_id[j]] > -1) /* Not removed */
            n_entries++;
        final->idx[zrow_id+1] = final->idx[zrow_id] + n_entries;

      } /* End of loop on final rows */

    } /* Columns reduction */

  } /* test msr2csr */

  /* Allocate arrays */
  nnz = final->idx[final->n_rows];
  BFT_MALLOC(final->col_id, nnz, int);

  if (init->type == CS_SLA_MAT_DEC) {

    BFT_MALLOC(final->sgn, nnz, short int);
    _pack_dec(final, init, row_z2i_ids, col_i2z_ids);

  }
  else { /* Not DEC */

    BFT_MALLOC(final->val, nnz, double);
    if (init->type == CS_SLA_MAT_CSR)
      _pack_csr(final, init, row_z2i_ids, col_i2z_ids);
    else if (init->type == CS_SLA_MAT_MSR)
      _pack_msr(final, init, row_z2i_ids, col_i2z_ids, msr2csr);

  } /* End test matrix type */

  return final;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build diagonal index
 *
 * \param[inout]  m   matrix to work with
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_diag_idx(cs_sla_matrix_t    *m)
{
  cs_lnum_t  row_id, pos;

  if (m == NULL)
    return;

  if (m->type == CS_SLA_MAT_CSR && m->n_rows == m->n_cols) {

    if (m->flag & CS_SLA_MATRIX_SHARED)
      bft_error(__FILE__, __LINE__, 0,
                _(" Cannot build diagonal index. Pattern is shared.\n"
                  " Stop execution.\n"));

    if (m->didx == NULL)
      BFT_MALLOC(m->didx, m->n_rows, int);

    for (row_id = 0; row_id < m->n_rows; row_id++) {
      m->didx[row_id] = -1; /* Default initialization: no diagonal entry */
      for (pos = m->idx[row_id]; pos < m->idx[row_id+1]; pos++) {
        if (m->col_id[pos] == row_id) {
          m->didx[row_id] = pos;
          break;
        }
      }
    } /* Loop on rows */

  } /* CSR + square */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the diagonal entries of a given matrix.
 *
 * \param[in]    m        matrix to work with
 * \param[inout] p_diag   pointer to diag array to define (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_get_diag(const cs_sla_matrix_t  *m,
                       double                 *p_diag[])
{
  int  i, j, k;

  double  *diag = *p_diag;

  /* Sanity checks */
  assert(m != NULL);
  assert(m->n_rows == m->n_cols); /* Only for squared matrices */
  assert(m->type == CS_SLA_MAT_MSR || m->type == CS_SLA_MAT_CSR);

  if (diag == NULL)
    BFT_MALLOC(diag, m->n_rows*m->stride, double);

  if (m->diag != NULL)
    memcpy(diag, m->diag, m->n_rows*m->stride*sizeof(double));

  else if (m->didx != NULL) { /* Diagonal is straightforward */

    if (m->stride > 1) {
      const int  stride = m->stride;

      for (i = 0; i < m->n_rows; i++) {
        if (m->didx[i] != -1) {
          for (k = 0; k < stride; k++)
            diag[stride*i+k] = m->val[m->didx[i]*stride+k];
        }
        else { /* No diagonal entry */
          for (k = 0; k < stride; k++)
            diag[stride*i+k] = 0;
        }
      } /* Loop on rows */

    }
    else { /* stride == 1 */

      for (i = 0; i < m->n_rows; i++) {
        if (m->didx[i] != -1)
          diag[i] = m->val[m->didx[i]];
        else /* No diagonal entry */
          diag[i] = 0;
      } /* Loop on rows */

    } /* stride == 1 */

  }
  else {

    /* Initialize diagonal entries */
    for (i = 0; i < m->n_rows*m->stride; i++)
      diag[i] = 0.0;

    /* Fill diagonal entries */
    if (m->stride > 1) {

      const int  stride = m->stride;

      for (i = 0; i < m->n_rows; i++) {
        for (j = m->idx[i]; j < m->idx[i+1]; j++) {
          if (m->col_id[j] == i) {
            for (k = 0; k < stride; k++)
              diag[i*stride+k] = m->val[j*stride+k];
            break;
          }
        }
      }

    }
    else {

      for (i = 0; i < m->n_rows; i++) {
        for (j = m->idx[i]; j < m->idx[i+1]; j++) {
          if (m->col_id[j] == i) {
            diag[i] = m->val[j];
            break;
          }
        }
      }

    } /* stride == 1 */

  } /* m->diag == NULL */

  /* return pointer */
  *p_diag = diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Sort each row by increasing colomn number
 *
 * \param[in]  mat     matrix to sort
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_sort(cs_sla_matrix_t  *m)
{
  int  ii;

  if (m == NULL)
    return;
  if (m->flag & CS_SLA_MATRIX_SORTED)
    return; /* Nothing to do */

  m->flag |= CS_SLA_MATRIX_SORTED;

  if (m->type == CS_SLA_MAT_CSR || m->type == CS_SLA_MAT_MSR)
    for (ii = 0; ii < m->n_rows; ii++)
      cs_sort_dcoupled_shell(m->idx[ii],
                             m->idx[ii + 1],
                             m->col_id,
                             m->val);

  else if (m->type == CS_SLA_MAT_DEC)
    for (ii = 0; ii < m->n_rows; ii++)
      cs_sort_sicoupled_shell(m->idx[ii],
                              m->idx[ii + 1],
                              m->col_id,
                              m->sgn);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read from a binary file a matrix in CSR format, its righ hand side
 *         and the solution. Matrix must have a stride equal to 1.
 *
 * \param[in]     name      name of the output file
 * \param[inout]  p_mat     system to solve
 * \param[inout]  p_rhs     right hand side
 * \param[inout]  p_sol       solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_bread(const char        *name,
             cs_sla_matrix_t  **p_mat,
             double            *p_rhs[],
             double            *p_sol[])
{
  int  n, n_rows, n_cols, nnz, flag;

  double  *rhs = NULL, *sol = NULL;
  cs_sla_matrix_t  *m = NULL;
  FILE *f = NULL;

  /* Sanity check */
  if (name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " No filename given. Can not read binary file!\n");

  f = fopen(name, "rb");

  /* Rhs : dim + values */
  fread(&n, sizeof(int), 1, f);
  BFT_MALLOC(rhs, n, double);
  BFT_MALLOC(sol, n, double);
  fread(rhs, n*sizeof(double), 1, f);
  fread(sol, n*sizeof(double), 1, f);

  /* Matrix flag 0 = CSR or 1 = MSR */
  fread(&flag, sizeof(int), 1, f);

  /* Matrix: dim + values */
  fread(&n_rows, sizeof(int), 1, f);
  fread(&n_cols, sizeof(int), 1, f);
  fread(&nnz, sizeof(int), 1, f);

  assert(n_rows == n_cols); /* This is not a real limitation */

  if (flag == 0)
    m = cs_sla_matrix_create(n_rows, n_cols, 1, CS_SLA_MAT_CSR, false);
  else {
    assert(flag == 1);
    m = cs_sla_matrix_create(n_rows, n_cols, 1, CS_SLA_MAT_MSR, false);

    fread(m->diag, n_rows*sizeof(double), 1, f);
  }

  /* index from 0 */
  fread(m->idx, (n_rows + 1)*sizeof(int), 1, f);

  if (nnz > 0) {

    BFT_MALLOC(m->col_id, nnz, int);
    BFT_MALLOC(m->val, nnz, double);

    /* cols from 1 to n_rows */
    fread(m->col_id, nnz*sizeof(int), 1, f);
    fread(m->val, nnz*sizeof(double), 1, f);

  }

  fclose(f);

  /* Return pointers */
  *p_mat = m;
  *p_rhs = rhs;
  *p_sol = sol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write in binary format a matrix in CSR format, its righ hand side
 *         and the solution
 *
 * \param[in]     name      name of the output file
 * \param[in]     A         system to solve
 * \param[in]     rhs       right hand side
 * \param[in]     x         solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_bwrite(const char              *name,
              const cs_sla_matrix_t   *m,
              const double            *rhs,
              const double            *sol)
{
  int  nnz, flag;

  FILE *f = NULL;

  /* Sanity check */
  assert(name != NULL);

  if (m == NULL) {
    fprintf(stdout, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);
    fprintf(stdout, " Stop file writing.\n");
    return;
  }
  else if (m->type == CS_SLA_MAT_NONE || m->type == CS_SLA_MAT_DEC) {
    fprintf(stdout, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);
    fprintf(stdout, "   type:   %8d\n", m->type);
    fprintf(stdout, " Stop file writing.\n");
    return;
  }
  else if (rhs == NULL) {
    fprintf(stdout, " Empty rhs array.\n");
    fprintf(stdout, " Stop file writing.\n");
    return;
  }
  else if (sol == NULL) {
    fprintf(stdout, " Empty sol array.\n");
    fprintf(stdout, " Stop file writing.\n");
    return;
  }

  if (m->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  f = fopen(name, "wb");

  /* array size + rhs + sol */
  fwrite(&(m->n_cols), sizeof(int), 1, f);
  fwrite(rhs, sizeof(double), m->n_cols, f);
  fwrite(sol, sizeof(double), m->n_cols, f);

  /* matrix type: 0 = CSR or 1 = MSR */
  flag = 0;
  if (m->type == CS_SLA_MAT_MSR)
    flag = 1;
  fwrite(&flag, sizeof(int), 1, f);

  /* matrix size and number of non-zeros */
  fwrite(&(m->n_rows), sizeof(int), 1, f);
  fwrite(&(m->n_cols), sizeof(int), 1, f);
  nnz = m->idx[m->n_rows];
  fwrite(&nnz, sizeof(int), 1, f);

  if (flag == 1) /* Dump diagonal */
    fwrite(m->diag, sizeof(double), m->n_rows, f);

  /* write index (start from 0) */
  fwrite(m->idx, sizeof(int), m->n_rows + 1, f);

  /* write cols (start from 1 to n_rows) and values */
  fwrite(m->col_id, sizeof(int), m->idx[m->n_rows], f);
  fwrite(m->val, sizeof(double), m->idx[m->n_rows], f);

  fclose(f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_sla_matrix_t structure
 *
 * \param[in]  n_rows     number of rows
 * \param[in]  n_cols     number of columns
 * \param[in]  stride     number of values related to each entry
 * \param[in]  type       kind of matrix
 * \param[in]  sym        true or false
 *
 * \return  a pointer to new allocated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_create(int                   n_rows,
                     int                   n_cols,
                     int                   stride,
                     cs_sla_matrix_type_t  type,
                     _Bool                 sym)
{
  int  i;

  cs_sla_matrix_t  *m = NULL;

  assert(n_rows > -1);
  assert(n_cols > -1);

  BFT_MALLOC(m, 1, cs_sla_matrix_t);

  m->properties = _get_default_matrix_properties();

  m->n_rows = n_rows;
  m->n_cols = n_cols;
  m->stride = stride;
  m->type = type;
  m->flag = 0;
  if (sym) m->flag |= CS_SLA_MATRIX_SYM;

  m->diag = NULL;
  m->idx = NULL;
  m->didx = NULL;
  m->col_id = NULL;
  m->sgn = NULL;
  m->val = NULL;

  if (m->type != CS_SLA_MAT_NONE) {

    BFT_MALLOC(m->idx, n_rows + 1, int);
    for (i = 0; i < n_rows + 1; i++)
      m->idx[i] = 0;

    if (m->type == CS_SLA_MAT_CSR && n_rows == n_cols) {
      BFT_MALLOC(m->didx, n_rows, int);
      for (i = 0; i < n_rows; i++)
        m->didx[i] = -1; /* Not set */
    }

    if (m->type == CS_SLA_MAT_MSR) {
      assert(n_rows == n_cols);
      BFT_MALLOC(m->diag, stride*n_rows, double);
      for (i = 0; i < stride*n_cols; i++)
        m->diag[i] = 0.0;
    }

  } /* Not an undefined type */

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_sla_matrix_t structure from an existing one
 *
 * \param[in]  ref      pointer to a reference matrix with the same pattern
 * \param[in]  type     type of the matrix to create
 * \param[in]  stride   true or false
 *
 * \return  a pointer to new allocated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_create_from_pattern(const cs_sla_matrix_t  *ref,
                                  cs_sla_matrix_type_t    type,
                                  int                     stride)
{
  int  i;

  cs_sla_matrix_t  *m = NULL;

  /* Sanity check */
  assert(ref != NULL);
  assert(ref->type != CS_SLA_MAT_NONE);

  BFT_MALLOC(m, 1, cs_sla_matrix_t);

  m->properties = _get_default_matrix_properties();

  m->type = type;
  m->n_rows = ref->n_rows;
  m->n_cols = ref->n_cols;
  m->stride = stride;
  m->diag = NULL;
  m->sgn = NULL;
  m->val = NULL;

  m->flag = ref->flag | CS_SLA_MATRIX_SHARED;
  m->idx = ref->idx;
  m->col_id = ref->col_id;
  m->didx = ref->didx;

  if (m->type == CS_SLA_MAT_MSR) {
    assert(m->n_rows == m->n_cols);
    BFT_MALLOC(m->diag, stride*m->n_rows, double);
    for (i = 0; i < stride*m->n_cols; i++)
      m->diag[i] = 0.0;
  }

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a new matrix structure from the copy of an existing one.
 *
 * \param[in] a       matrix to copy
 * \param[in] shared  true: only pointer are copied other allocated
 *
 * \return the matrix corresponding to this operation
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_copy(const cs_sla_matrix_t  *a,
                   _Bool                   shared)
{
  cs_sla_matrix_t  *b = NULL;

  if (a == NULL)
    return b;

  assert(a->stride == 1); /* Done only for this case up to now */

  if (shared)
    b = cs_sla_matrix_create_from_pattern(a, a->type, a->stride);

  else { /* Create and copy pattern */

    if (a->flag & CS_SLA_MATRIX_SYM)
      b = cs_sla_matrix_create(a->n_rows, a->n_cols, a->stride, a->type, true);
    else
      b = cs_sla_matrix_create(a->n_rows, a->n_cols, a->stride, a->type, false);

    if (a->type != CS_SLA_MAT_NONE) { /* Copy pattern */

      BFT_MALLOC(b->col_id, a->idx[a->n_rows], cs_lnum_t);
      memcpy(b->idx, a->idx, (a->n_rows+1)*sizeof(cs_lnum_t));
      memcpy(b->col_id, a->col_id, a->idx[a->n_rows]*sizeof(cs_lnum_t));

      if (a->didx != NULL) {
        BFT_MALLOC(b->didx, a->n_rows, int);
        memcpy(b->didx, a->didx, a->n_rows*sizeof(int));
      }

    } /* type != CS_SLA_MAT_NONE */

  } /* Shared or not */

  if (a->type != CS_SLA_MAT_NONE) {

    if (a->type == CS_SLA_MAT_DEC) {
      BFT_MALLOC(b->sgn, a->idx[a->n_rows], short int);
      memcpy(b->sgn, a->sgn, a->idx[a->n_rows]*sizeof(short int));
    }
    else { /* CSR or MSR */

      BFT_MALLOC(b->val, a->idx[a->n_rows], double);
      memcpy(b->val, a->val, a->idx[a->n_rows]*sizeof(double));

      if (a->type == CS_SLA_MAT_MSR) /* Copy also the diagonal */
        memcpy(b->diag, a->diag, a->n_rows*sizeof(double));

    } /* type */

  } /* type != CS_SLA_MAT_NONE */

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Transpose a cs_sla_matrix_t structure
 *
 * \param[in]  mat     matrix to transpose
 *
 * \return  a pointer to new allocated matrix which is transposed to mat
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_transpose(const cs_sla_matrix_t  *a)
{
  int  i, j, nnz, shift, row_id;

  int  *count = NULL;
  cs_sla_matrix_t  *at = NULL;

  if (a == NULL)
    return at;

  if (a->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  /* Simple cases: symmetric matrix */
  if (a->flag & CS_SLA_MATRIX_SYM) {
    at = cs_sla_matrix_copy(a, true); /* share pattern with a */
    return at;
  }

  at = cs_sla_matrix_create(a->n_cols, a->n_rows, a->stride, a->type, false);

  if (a->type == CS_SLA_MAT_NONE)
    return at;

  /* General case: 1) build index */
  BFT_MALLOC(at->col_id, a->idx[a->n_rows], cs_lnum_t);
  BFT_MALLOC(count, at->n_rows, int);

  for (i = 0; i < a->n_rows; i++)
    for (j = a->idx[i]; j < a->idx[i+1]; j++)
      at->idx[a->col_id[j]+1] += 1;

  for (i = 0; i < at->n_rows; i++) {
    count[i] = 0;
    at->idx[i+1] += at->idx[i];
  }

  nnz = a->idx[a->n_rows];
  assert(at->idx[at->n_rows] == nnz);

  switch (a->type) {

  case CS_SLA_MAT_DEC:

    BFT_MALLOC(at->sgn, nnz, short int);

    for (i = 0; i < a->n_rows; i++) {
      for (j = a->idx[i]; j < a->idx[i+1]; j++) {

        row_id = a->col_id[j];
        shift = at->idx[row_id] + count[row_id];
        at->col_id[shift] = i;
        at->sgn[shift] = a->sgn[j];
        count[row_id] += 1;

      }
    }
    break;

  case CS_SLA_MAT_MSR:

    assert(at->diag != NULL);
    for (i = 0; i < a->n_rows; i++)
      at->diag[i] = a->diag[i];

  case CS_SLA_MAT_CSR:

    BFT_MALLOC(at->val, nnz, double);

    for (i = 0; i < a->n_rows; i++) {
      for (j = a->idx[i]; j < a->idx[i+1]; j++) {

        row_id = a->col_id[j];
        shift = at->idx[row_id] + count[row_id];
        at->col_id[shift] = i;
        at->val[shift] = a->val[j];
        count[row_id] += 1;

      }
    }

    cs_sla_matrix_diag_idx(at);

    break;

  default:
    break;

  } /* End of switch */

  BFT_FREE(count);

  return at;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_sla_matrix_t structure
 *
 * \param[in]  m     matrix to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_free(cs_sla_matrix_t  *m)
{
  if (m == NULL)
    return NULL;

  if (m->type != CS_SLA_MAT_NONE) {

    switch (m->type) {

    case CS_SLA_MAT_DEC:
      BFT_FREE(m->sgn);
      break;
    case CS_SLA_MAT_CSR:
      BFT_FREE(m->val);
      if (m->diag != NULL)
        BFT_FREE(m->diag);
      break;
    case CS_SLA_MAT_MSR:
      BFT_FREE(m->val);
      BFT_FREE(m->diag);
      break;
    default:
      break;

    } /* End switch */

    /* Common pointers defining the pattern */
    if (!(m->flag & CS_SLA_MATRIX_SHARED)) {
      BFT_FREE(m->idx);
      BFT_FREE(m->col_id);
      if (m->didx != NULL) /* Should be only for square CSR matrix */
        BFT_FREE(m->didx);
    }

    if (m->properties != NULL)
      BFT_FREE(m->properties);

  } /* TYPE_NONE */

  BFT_FREE(m);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Remove entries in a cs_sla_matrix_t structure below a given
 *          threshold. |a(i,j)| < eps * max|a(i,j)|
 *                                       i
 *
 * \param[inout]  mat     matrix to clean
 * \param[in]     eps     value of the threshold
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_clean(cs_sla_matrix_t   *m,
                    double             eps)
{
  int  i, j, shift, start, init_nnz;
  double  max_val, limit;

  /* Sanity checks */
  assert(m != NULL);

  if (m->stride > 1)
    bft_error(__FILE__, __LINE__, 0, _sla_err_stride);

  init_nnz = m->idx[m->n_rows];

  /* Find and delete null entries */
  if (m->type == CS_SLA_MAT_MSR || m->type == CS_SLA_MAT_CSR) {

    shift = m->idx[0];
    start = m->idx[0];

    for (i = 0; i < m->n_rows; i++) {

      /* Get max absolute value and compute threshold */
      max_val = -DBL_MAX;
      for (j = start; j < m->idx[i+1]; j++)
        max_val = fmax(fabs(m->val[j]), max_val);
      limit = max_val*eps;

      /* Delete too small entries */
      for (j = start; j < m->idx[i+1]; j++) {

        if (fabs(m->val[j]) > limit) { /* Keep value */
          if (j != shift) { /* Need to copy */
            m->val[shift] = m->val[j];
            m->col_id[shift] = m->col_id[j];
          }
          shift++;
        }

      } /* End of loop on row entries */

      start = m->idx[i+1];
      m->idx[i+1] = shift;

    } /* End of loop on rows */

    if (init_nnz != shift) {
      BFT_REALLOC(m->val, shift, double);
      BFT_REALLOC(m->col_id, shift, int);
    }

    /* Output message */
    bft_printf(" -dbg- cs_sla_matrix_clean >>"
               " type: %s; n_rows: %6d; threshold: %6.3e; nnz: %d -> %d\n",
               _sla_matrix_type[m->type], m->n_rows, eps, init_nnz, shift);

  } /* MSR or CSR */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the number of non-zeros (nnz) elements in a matrix
 *
 * \param[in]  mat     matrix
 *
 * \return the number of nnz
 */
/*----------------------------------------------------------------------------*/

size_t
cs_sla_matrix_get_nnz(const cs_sla_matrix_t   *m)
{
  size_t  nnz = 0;

  if (m == NULL)
    return nnz;

  if (m->type == CS_SLA_MAT_NONE)
    return nnz;

  nnz = m->idx[m->n_rows];
  if (m->type == CS_SLA_MAT_MSR)
    nnz += m->n_rows;

  return nnz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute general information related to a cs_sla_matrix_t structure
 *
 * \param[in]  m           matrix to analyse
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_info_t
cs_sla_matrix_analyse(const cs_sla_matrix_t   *m)
{
  int  i;
  cs_sla_matrix_info_t  minfo;
  cs_data_info_t  dinfo;

  int  *row_size = NULL;

  /* Default initialization */
  minfo.stencil_min = 0;
  minfo.stencil_max = 0;
  minfo.stencil_mean = 0;
  minfo.nnz = 0;
  minfo.fillin = 0;

  if (m == NULL)
    return minfo;
  if (m->type == CS_SLA_MAT_NONE)
    return minfo;

  /* Compute information */
  minfo.nnz = cs_sla_matrix_get_nnz(m);
  minfo.fillin = (double)(minfo.nnz/m->n_rows)*(100./m->n_cols);

  /* Stencil */
  BFT_MALLOC(row_size, m->n_rows, int);
  for (i = 0; i < m->n_rows; i++)
    row_size[i] = m->idx[i+1] - m->idx[i];

  if (m->type == CS_SLA_MAT_MSR)
    for (i = 0; i < m->n_rows; i++)
      row_size[i] += 1;

  dinfo = cs_analysis_data(m->n_rows,  // n_elts
                           1,          // stride
                           CS_INT32,   // datatype
                           row_size,   // data
                           false);     // do_abs

  BFT_FREE(row_size);

  /* Stencil information */
  minfo.stencil_min = dinfo.min.number;
  minfo.stencil_max = dinfo.max.number;
  minfo.stencil_mean = dinfo.mean;

  return minfo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Synthesis of a cs_sla_matrix_t structure
 *
 * \param[in]  name        either name of the file if f is NULL or description
 * \param[in]  f           pointer to a FILE struct.
 * \param[in]  m           matrix to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_resume(const char              *name,
                     FILE                    *f,
                     const cs_sla_matrix_t   *m)
{
  cs_sla_matrix_info_t  minfo;
  cs_data_info_t  dinfo;

  char  *filename = NULL;
  FILE  *_f = f;
  _Bool close_file = false;

  if (_f == NULL) {
    if (name == NULL)
      _f = stdout;

    else {

      int len = strlen(name) + strlen("-resume.log") + 1;

      BFT_MALLOC(filename, len, char);
      sprintf(filename,"%s-resume.log", name);
      _f = fopen(filename,"w");
      close_file = true;
    }
  }

  fprintf(_f, "\n");
  if (m == NULL)
    fprintf(_f, " -sla-  SLA matrix structure: %p (%s)\n",
            (const void *)m, name);

  else if (m->type == CS_SLA_MAT_NONE) {
    fprintf(_f, " -sla-  SLA matrix structure: %p (%s)\n",
            (const void *)m, name);
    fprintf(_f, " -sla-  type:        %s\n", _sla_matrix_type[m->type]);
  }
  else {

    minfo = cs_sla_matrix_analyse(m);

    fprintf(_f, " -sla-  SLA matrix structure: %p (%s)\n",
            (const void *)m, name);
    fprintf(_f, " -sla-  type          %s\n", _sla_matrix_type[m->type]);
    fprintf(_f, " -sla-  n_rows        %d\n", m->n_rows);
    fprintf(_f, " -sla-  n_cols        %d\n", m->n_cols);
    fprintf(_f, " -sla-  stride        %d\n", m->stride);
    fprintf(_f, " -sla-  nnz           %lu\n", minfo.nnz);
    fprintf(_f, " -sla-  fill-in       %5.2e %%\n", minfo.fillin);
    fprintf(_f, " -sla-  stencil_min   %d\n", minfo.stencil_min);
    fprintf(_f, " -sla-  stencil_max   %d\n", minfo.stencil_max);
    fprintf(_f, " -sla-  stencil_mean  %5.2e\n", minfo.stencil_mean);
    if (m->flag & CS_SLA_MATRIX_SYM)
      fprintf(_f, " -sla-  sym           True\n");
    else
      fprintf(_f, " -sla-  sym           False\n");

    if (m->type == CS_SLA_MAT_MSR ) {
      dinfo = cs_analysis_data(m->n_rows,    // n_elts
                               1,            // stride
                               CS_DOUBLE,    // datatype
                               m->diag,      // data
                               false);       // work also with abs ?
      cs_data_info_dump("mat->diag", _f, m->n_rows, CS_DOUBLE, dinfo);
    }

    if (m->type == CS_SLA_MAT_CSR || m->type == CS_SLA_MAT_MSR) {
      dinfo = cs_analysis_data(minfo.nnz,    // n_elts
                               1,            // stride
                               CS_DOUBLE,    // datatype
                               m->val,       // data
                               false);       // work also with abs ?
      cs_data_info_dump("mat->val", _f, minfo.nnz, CS_DOUBLE, dinfo);
    }

  } /* m neither NULL nor CS_SLA_MAT_NONE */

  if (close_file) {
    BFT_FREE(filename);
    fclose(_f);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_sla_matrix_t structure
 *
 * \param[in]  name        either name of the file if f is NULL or description
 * \param[in]  f           pointer to a FILE struct.
 * \param[in]  m           matrix to dump
 * \param[in]  dump_level  level of information displayed
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_dump(const char              *name,
                   FILE                    *f,
                   const cs_sla_matrix_t   *m)
{
  int i, j, k;

  FILE  *_f = f;
  _Bool close_file = false;

  if (_f == NULL) {
    if (name == NULL)  _f = stdout;
    else
      _f = fopen(name,"w"), close_file = true;
  }

  if (m == NULL)
    fprintf(_f, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);

  else if (m->type == CS_SLA_MAT_NONE) {
    fprintf(_f, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);
    fprintf(_f, "   type:        %s\n", _sla_matrix_type[m->type]);
  }
  else {

    size_t  max_nnz = m->n_rows*m->n_cols;
    size_t  nnz = cs_sla_matrix_get_nnz(m);
    double  fill_in = nnz*100.0/max_nnz;

    fprintf(_f, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);
    fprintf(_f, "   n_rows: %8d\n", m->n_rows);
    fprintf(_f, "   n_cols: %8d\n", m->n_cols);
    fprintf(_f, "   nnz:    %lu\n", nnz);
    fprintf(_f, "   fill-in:     %.2f\n", fill_in);
    fprintf(_f, "   stride: %8d\n", m->stride);
    fprintf(_f, "   type:        %s\n", _sla_matrix_type[m->type]);
    if (m->flag & CS_SLA_MATRIX_SYM)
      fprintf(_f, "   sym:    True\n\n");
    else
      fprintf(_f, "   sym:    False\n\n");

    if (m->diag != NULL) { /* Dump diagonal */
      fprintf(_f, "  Diagonal entries:\n");
      for (i = 0; i < m->n_rows; i++) {
        fprintf(_f, " %5d >", i+1);
        for (k = 0; k < m->stride; k++)
          fprintf(_f, " % -8.4e", m->diag[i*m->stride+k]);
        fprintf(_f, "\n");
      }
    } /* diag != NULL */

    if (nnz > 0) {

      for (i = 0; i < m->n_rows; i++) {

        int  s = m->idx[i], e = m->idx[i+1];

        fprintf(_f, "%5d >", i+1);
        if (m->type == CS_SLA_MAT_DEC) { // DEC matrix
          for (j = s; j < e; j++) {
            for (k = 0; k < m->stride; k++)
              fprintf(_f, " %2d", m->sgn[j*m->stride+k]);
            fprintf(_f, " (%5d)", m->col_id[j]);
          }
        }
        else if (m->type == CS_SLA_MAT_CSR || m->type == CS_SLA_MAT_MSR) {
          for (j = s; j < e; j++) {
            for (k = 0; k < m->stride; k++)
              fprintf(_f, " % -8.4e", m->val[j*m->stride+k]);
            fprintf(_f, " (%5d)", m->col_id[j]);
          }
        }
        fprintf(_f, "\n");

      }
    } /* nnz > 0 */

  } /* m neither NULL nor CS_SLA_MAT_NONE */

  if (close_file)
    fclose(_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_sla_matrix_t structure and its related right-hand side
 *
 * \param[in]  name        either name of the file if f is NULL or description
 * \param[in]  f           pointer to a FILE struct.
 * \param[in]  m           matrix to dump
 * \param[in]  rhs         right-hand side to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_system_dump(const char              *name,
                   FILE                    *f,
                   const cs_sla_matrix_t   *m,
                   const double            *rhs)
{
  int i, j, k;

  FILE  *_f = f;
  _Bool close_file = false;

  if (_f == NULL) {
    if (name == NULL)  _f = stdout;
    else
      _f = fopen(name,"w"), close_file = true;
  }

  if (m == NULL)
    fprintf(_f, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);

  else if (m->type == CS_SLA_MAT_NONE) {
    fprintf(_f, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);
    fprintf(_f, "   type:        %s\n", _sla_matrix_type[m->type]);
  }
  else {

    size_t  max_nnz = m->n_rows*m->n_cols;
    size_t  nnz = cs_sla_matrix_get_nnz(m);
    double  fill_in = nnz*100.0/max_nnz;

    /* Sanity check */
    assert(rhs != NULL);

    fprintf(_f, "\n SLA matrix structure: %p (%s)\n", (const void *)m, name);
    fprintf(_f, "   n_rows: %8d\n", m->n_rows);
    fprintf(_f, "   n_cols: %8d\n", m->n_cols);
    fprintf(_f, "   nnz:    %lu\n", nnz);
    fprintf(_f, "   fill-in:     %.2f\n", fill_in);
    fprintf(_f, "   stride: %8d\n", m->stride);
    fprintf(_f, "   type:        %s\n", _sla_matrix_type[m->type]);
    if (m->flag & CS_SLA_MATRIX_SYM)
      fprintf(_f, "   sym:    True\n\n");
    else
      fprintf(_f, "   sym:    False\n\n");

    for (i = 0; i < m->n_rows; i++) {

      int  s = m->idx[i], e = m->idx[i+1];

      fprintf(_f, "\nrow: %3d >> rhs: % -8.4e", i+1, rhs[i]);

      if (m->type == CS_SLA_MAT_DEC) { // DEC matrix

        if (m->diag != NULL) {
          fprintf(_f, " <diag:");
          for (k = 0; k < m->stride; k++)
            fprintf(_f, " % -8.4e", m->diag[i*m->stride+k]);
          fprintf(_f, ">\t");
        }

        for (j = s; j < e; j++) {
          fprintf(_f, " <col: %3d;", m->col_id[j]);
          for (k = 0; k < m->stride; k++)
            fprintf(_f, " %2d", m->sgn[j*m->stride+k]);
          fprintf(_f,">");
        }

      }
      else if (m->type == CS_SLA_MAT_CSR || m->type == CS_SLA_MAT_MSR) {

        if (m->diag != NULL) {
          fprintf(_f, " <diag:");
          for (k = 0; k < m->stride; k++)
            fprintf(_f, " % -8.4e", m->diag[i*m->stride+k]);
          fprintf(_f, ">\t");
        }

        for (j = s; j < e; j++) {
          fprintf(_f, " <col: %3d;", m->col_id[j]);
          for (k = 0; k < m->stride; k++)
            fprintf(_f, " % -8.4e", m->val[j*m->stride+k]);
          fprintf(_f, ">");
        }
      }

    }

  } /* m neither NULL nor CS_SLA_MAT_NONE */

  if (close_file)
    fclose(_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assemble a MSR matrix from local contributions
 *          --> We assume that the local matrices are symmetric
 *          --> We assume that the assembled matrix has its columns sorted
 *
 * \param[in]     loc     pointer to a local matrix
 * \param[inout]  ass     pointer to a cs_sla_matrix_t struct. collecting data
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_assemble_msr(const cs_toolbox_locmat_t  *loc,
                    cs_sla_matrix_t            *ass)
{
  int  i, j, k_ij, k_ji;

  const int  n_ent = loc->n_ent;

  assert(ass->type == CS_SLA_MAT_MSR);
  assert(ass->flag & (CS_SLA_MATRIX_SYM | CS_SLA_MATRIX_SORTED));

  /* Assemble local contributions into ass */
  for (i = 0; i < n_ent; i++) {

    short int  pos_i = i*n_ent;
    cs_lnum_t  i_id = loc->ids[i];
    cs_lnum_t  start_i = ass->idx[i_id];
    size_t  n_i_ents = ass->idx[i_id+1] - start_i; // arg. binary search

    /* Add diagonal term : loc(i,i) */
    ass->diag[i_id] += loc->mat[pos_i+i];

    for (j = i + 1; j < n_ent; j++) {

      double  val_ij = loc->mat[pos_i+j];

      if (fabs(val_ij) > cs_base_zthreshold) { /* Not zero */

        int  j_id = loc->ids[j];
        cs_lnum_t  start_j = ass->idx[j_id];
        size_t  n_j_ents = ass->idx[j_id+1] - start_j;

        /* First add: loc(i,j) */
        k_ij = cs_search_binary(n_i_ents, j_id, &(ass->col_id[start_i]));
        assert(k_ij > -1);
        ass->val[start_i+k_ij] += val_ij;

        /* Second add: loc(j,i) */
        k_ji = cs_search_binary(n_j_ents, i_id, &(ass->col_id[start_j]));
        assert(k_ji > -1);
        ass->val[start_j+k_ji] += val_ij;

      } /* loc[ij] != 0.0 */

    } /* Loop on j (Add extra-diag terms) */

  } /* Loop on i (Add diagonal term) */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
