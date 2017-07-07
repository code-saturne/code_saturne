/*============================================================================
 * Basic operations: dot product, cross product, sum...
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo.h"
#include "cs_blas.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_sort.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm 2 of a vector of size len
 *         This algorithm tries to reduce round-off error thanks to
 *         intermediate sums.
 *
 *  \param[in] len     vector dimension
 *  \param[in] v       vector
 *
 * \return  the euclidean norm of a vector
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_euclidean_norm(int               len,
                  const cs_real_t   v[])
{
  cs_real_t  n2 = DBL_MAX;

  assert(v != NULL);

  if (len < 1 && cs_glob_n_ranks == 1)
    return 0.0;

  n2 = cs_dot(len, v, v);
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &n2);

  if (n2 > -DBL_MIN)
    n2 = sqrt(n2);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop norm computation. Norm value is < 0 !\n"));

  return n2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted sum of squared values of an array
 *
 *  \param[in] n       size of arrays x and weight
 *  \param[in] x       array of floating-point values
 *  \param[in] weight  floating-point values of weights
 *
 * \return the result of this operation
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_weighted_sum_squared(cs_lnum_t                   n,
                        const cs_real_t *restrict   x,
                        const cs_real_t *restrict   weight)
{
  if (n == 0 && cs_glob_n_ranks == 1)
    return 0.0;

  /* Sanity check */
  if (weight == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Weighted operation needs weigth array to be allocated.\n"
                " Stop execution.\n"));

  const cs_lnum_t  block_size = 60;
  const cs_lnum_t  n_blocks = n / block_size;
  const cs_lnum_t  n_sblocks = sqrt(n_blocks);
  const cs_lnum_t  blocks_in_sblocks = (n_sblocks > 0) ? n_blocks/n_sblocks : 0;

  cs_real_t result = 0.0;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

# pragma omp parallel for reduction(+:result) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    cs_real_t sresult = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cs_real_t _result = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++)
        _result += weight[i]*x[i]*x[i];
      sresult += _result;
    }

    result += sresult;

  }

  /* Remainder */
  cs_real_t _result = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++)
    _result += weight[i]*x[i]*x[i];
  result += _result;

  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_locmat_t structure
 *
 * \param[in]  n_max_ent    max number of entities
 *
 * \return  a new allocated cs_locmat_t structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_locmat_create(int   n_max_ent)
{
  int  i;

  cs_locmat_t  *lm = NULL;

  BFT_MALLOC(lm, 1, cs_locmat_t);

  lm->n_max_ent = n_max_ent;
  lm->n_ent = 0;
  lm->ids = NULL;
  lm->val = NULL;

  if (n_max_ent > 0) {

    int  msize = n_max_ent*n_max_ent;

    BFT_MALLOC(lm->ids, n_max_ent, cs_lnum_t);
    for (i = 0; i < n_max_ent; i++)
      lm->ids[i] = 0;

    BFT_MALLOC(lm->val, msize, double);
    for (i = 0; i < msize; i++)
      lm->val[i] = 0;

  }

  return lm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_locmat_t structure
 *
 * \param[in]  lm    pointer to a cs_locmat_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_locmat_free(cs_locmat_t  *lm)
{
  if (lm == NULL)
    return lm;

  if (lm->n_max_ent > 0) {
    BFT_FREE(lm->ids);
    BFT_FREE(lm->val);
  }

  BFT_FREE(lm);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Copy a cs_locmat_t structure into another cs_locmat_t structure
 *          which has been already allocated
 *
 * \param[in, out]  recv    pointer to a cs_locmat_t struct.
 * \param[in]       send    pointer to a cs_locmat_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_copy(cs_locmat_t        *recv,
               const cs_locmat_t  *send)
{
  /* Sanity check */
  assert(recv->n_max_ent >= send->n_max_ent);

  recv->n_ent = send->n_ent;

  /* Copy ids */
  for (int  i = 0; i < send->n_ent; i++)
    recv->ids[i] = send->ids[i];

  /* Copy values */
  memcpy(recv->val, send->val, sizeof(double)*send->n_ent*send->n_ent);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a local dense matrix-vector product
 *          matvec has been previously allocated
 *
 * \param[in]      loc    local matrix to use
 * \param[in]      vec    local vector to use
 * \param[in, out] matvec result of the local matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_matvec(const cs_locmat_t   *loc,
                 const cs_real_t     *vec,
                 cs_real_t           *matvec)
{
  /* Sanity checks */
  assert(loc != NULL && vec != NULL && matvec != NULL);

  const int  n = loc->n_ent;

  /* Init. matvec */
  cs_real_t  v = vec[0];
  for (int i = 0; i < n; i++)
    matvec[i] = v*loc->val[i*n];

  /* Increment matvec */
  for (int i = 0; i < n; i++) {
    int shift =  i*n;
    for (int j = 1; j < n; j++)
      matvec[i] += vec[j]*loc->val[shift + j];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two local dense matrices: loc += add
 *
 * \param[in, out] loc   local matrix storing the result
 * \param[in]      add   values to add to loc
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_add(cs_locmat_t        *loc,
              const cs_locmat_t  *add)
{
  /* Sanity checks */
  assert(loc != NULL && add != NULL);
  assert(loc->n_max_ent <= add->n_max_ent);

  for (int i = 0; i < loc->n_ent*loc->n_ent; i++)
    loc->val[i] += add->val[i];

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Give the result of the following operation: loc = loc + alpha*add
 *
 * \param[in, out] loc    local matrix storing the result
 * \param[in]      alpha  multiplicative coefficient
 * \param[in]      add    values to add to loc
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_mult_add(cs_locmat_t        *loc,
                   cs_real_t              alpha,
                   const cs_locmat_t  *add)
{
  /* Sanity checks */
  assert(loc != NULL && add != NULL);
  assert(loc->n_max_ent <= add->n_max_ent);

  for (int i = 0; i < loc->n_ent*loc->n_ent; i++)
    loc->val[i] += alpha*add->val[i];

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a new matrix by adding a local matrix with its transpose.
 *          Keep the transposed matrix for future use.
 *
 * \param[in, out] loc   local matrix to transpose and add
 * \param[in, out] tr    transposed of the local matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_add_transpose(cs_locmat_t  *loc,
                        cs_locmat_t  *tr)
{
  /* Sanity check */
  assert(loc != NULL && tr != NULL && tr->n_max_ent == loc->n_max_ent);

  if (loc->n_ent < 1)
    return;

  tr->n_ent = loc->n_ent;

  for (int i = 0; i < loc->n_ent; i++) {

    int  ii = i*loc->n_ent + i;

    tr->ids[i] = loc->ids[i];
    tr->val[ii] = loc->val[ii];
    loc->val[ii] *= 2;

    for (int j = i+1; j < loc->n_ent; j++) {

      int  ij = i*loc->n_ent + j;
      int  ji = j*loc->n_ent + i;

      tr->val[ji] = loc->val[ij];
      tr->val[ij] = loc->val[ji];
      loc->val[ij] += tr->val[ij];
      loc->val[ji] += tr->val[ji];

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the given matrix into its anti-symmetric part
 *
 * \param[in, out] loc   local matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_asymm(cs_locmat_t  *loc)
{
  /* Sanity check */
  assert(loc != NULL);

  if (loc->n_ent < 1)
    return;

  for (int i = 0; i < loc->n_ent; i++) {

    double  *mi = loc->val + i*loc->n_ent;

    mi[i] = 0;

    for (int j = i+1; j < loc->n_ent; j++) {

      int  ji = j*loc->n_ent + i;

      mi[j] = 0.5*(mi[j] - loc->val[ji]);
      loc->val[ji] = mi[j];

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local discrete Hodge operator
 *
 * \param[in]    parent_id  id of the related parent entity
 * \param[in]    lm         pointer to the cs_sla_locmat_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_dump(int                 parent_id,
               const cs_locmat_t  *lm)
{
  int  i, j;

  cs_log_printf(CS_LOG_DEFAULT, "\n  << parent id: %d >>\n", parent_id);

  /* List sub-entity ids */
  for (i = 0; i < lm->n_ent; i++)
    cs_log_printf(CS_LOG_DEFAULT, " %9d", lm->ids[i]);
  cs_log_printf(CS_LOG_DEFAULT, "\n");

  for (i = 0; i < lm->n_ent; i++) {
    cs_log_printf(CS_LOG_DEFAULT, " %5d", lm->ids[i]);
    for (j = 0; j < lm->n_ent; j++)
      cs_log_printf(CS_LOG_DEFAULT, " % .4e", lm->val[i*lm->n_ent+j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump an array in the listing (For DEBUG)
 */
/*----------------------------------------------------------------------------*/

void
cs_dump_array_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_real_t    array[],
                         int                n_cols)
{
  cs_log_printf(CS_LOG_DEFAULT, "\nDUMP>> %s\n", header);

  if (n_cols < 1) n_cols = 1;
  int  n_rows = size/n_cols;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    for (cs_lnum_t j = i*n_cols; j < (i+1)*n_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6.4e", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

  if (n_rows*n_cols < size) {
    for (cs_lnum_t j = n_rows*n_cols; j < size; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6.4e", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump an array in the listing (For DEBUG)
 */
/*----------------------------------------------------------------------------*/

void
cs_dump_integer_to_listing(const char        *header,
                           const cs_lnum_t    size,
                           const cs_lnum_t    array[],
                           int                n_cols)
{
  cs_log_printf(CS_LOG_DEFAULT, "\nDUMP>> %s\n", header);

  if (n_cols < 1) n_cols = 1;
  int  n_rows = size/n_cols;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    for (cs_lnum_t j = i*n_cols; j < (i+1)*n_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6d", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

  if (n_rows*n_cols < size) {
    for (cs_lnum_t j = n_rows*n_cols; j < size; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6d", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
