/*============================================================================
 * Basic operations: dot product, cross product, sum...
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

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo.h"
#include "cs_blas.h"
#include "cs_math.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Intialize by default a cs_data_info_t structure according to the
 *          datatype
 *
 * \param[in]  datatype
 *
 * \return a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_data_info_t
_init_dinfo(cs_datatype_t   datatype)
{
  cs_data_info_t  info;

  info.mean = 0.0;
  info.sigma = 0.0;
  info.euclidean_norm = 0.0;

  switch (datatype) {

  case CS_DOUBLE:
    info.min.value = DBL_MAX;
    info.max.value = -DBL_MAX;
    break;

  case CS_INT32:
    info.min.number = INT_MAX;
    info.max.number = -INT_MAX;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid datatype for analysing data.\n"));
  }

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute simple information about an array of data.
 *          >> Algorithm from Mark Hoemmen (U.C. Berkeley)
 *
 * \param[in]      n_elts    number of couples in data
 * \param[in]      data      buffer containing input data
 * \param[in, out] info      pointer to a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_info_double(cs_lnum_t         n_elts,
                     const cs_real_t   data[],
                     cs_data_info_t   *info)
{
  if (n_elts == 0)
    return;

  /* Compute statistics */
  for (cs_lnum_t i = 0; i < n_elts; i++) {

    const cs_real_t  val = data[i];
    cs_real_t  delta = val - info->mean;
    cs_real_t  chi = delta/(i+1);

    if (info->min.value > val)  info->min.value = val;
    if (info->max.value < val)  info->max.value = val;

    info->sigma += i*chi*delta;
    info->mean += chi;

  }

  info->sigma = sqrt(fabs(info->sigma)/n_elts);
  info->euclidean_norm = cs_euclidean_norm(n_elts, data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute simple information about an array of data.
 *          >> Algorithm from Mark Hoemmen (U.C. Berkeley)
 *
 * \param[in]      n_elts    number of couples in data
 * \param[in]      data      buffer containing input data
 * \param[in, out] info      pointer to a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_info_int32(cs_lnum_t         n_elts,
                    const cs_lnum_t   data[],
                    cs_data_info_t   *info)
{
  if (n_elts == 0)
    return;

  /* Compute statistics */
  for (cs_lnum_t i = 0; i < n_elts; i++) {

    const cs_lnum_t  val = data[i];
    double  delta = val - info->mean;
    double  chi = delta/(i+1);

    if (info->min.number > val)  info->min.number = val;
    if (info->max.number < val)  info->max.number = val;

    info->sigma += i*chi*delta;
    info->mean += chi;
    info->euclidean_norm += val*val;

  }

  info->sigma = sqrt(fabs(info->sigma)/n_elts);
  info->euclidean_norm = sqrt(info->euclidean_norm);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm 2 of a vector of size len
 *         This algorithm tries to reduce round-off error thanks to
 *          intermediate sums.
 *
 *  \param[in] len     vector dimension
 *  \param[in] v       vector
 *
 * \return  the euclidean norm of a vector
 */
/*----------------------------------------------------------------------------*/

double
cs_euclidean_norm(int            len,
                  const double   v[])
{
  double  n2 = DBL_MAX;

  if (len < 1 || v == NULL)
    return 0.0;

  n2 = cs_dot(len, v, v);
  if (n2 > -DBL_MIN)
    n2 = sqrt(n2);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop norm computation. Norm value is < 0 !\n"));

  return n2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted sum of square values of an array
 *
 *  \param[in] n       size of arrays x and weight
 *  \param[in] x       array of floating-point values
 *  \param[in] weight  floating-point values of weights
 *
 * \return the result of this operation
 */
/*----------------------------------------------------------------------------*/

double
cs_weighted_sum_square(cs_lnum_t                n,
                       const double *restrict   x,
                       const double *restrict   weight)
{
  if (n == 0)
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

  double result = 0.0;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

# pragma omp parallel for reduction(+:result) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sresult = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double _result = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++)
        _result += weight[i]*x[i]*x[i];
      sresult += _result;
    }

    result += sresult;

  }

  /* Remainder */
  double _result = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++)
    _result += weight[i]*x[i]*x[i];
  result += _result;

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate or reallocate a temporary buffer structure
 *
 * \param[in]       bufsize   reference size
 * \param[in, out]  p_tb      pointer to the temporary structure to allocate
 */
/*----------------------------------------------------------------------------*/

void
cs_tmpbuf_alloc(size_t           bufsize,
                cs_tmpbuf_t    **p_tb)
{
  cs_tmpbuf_t  *tb = *p_tb;

  if (bufsize == 0)
    return;

  if (tb == NULL) { /* Creation */
    BFT_MALLOC(tb, 1, cs_tmpbuf_t);
    tb->bufsize = bufsize;
    BFT_MALLOC(tb->buf, bufsize, char);
  }
  else { /* Reallocation if needed */
    if (tb->bufsize < bufsize) {
      tb->bufsize = bufsize;
      BFT_REALLOC(tb->buf, bufsize, char);
    }
  }

  /* Returns buffer structure */
  *p_tb = tb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a temporary buffer structure
 *
 * \param[in]  tb      pointer to the temporary structure to free
 *
 * \returns NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_tmpbuf_t *
cs_tmpbuf_free(cs_tmpbuf_t  *tb)
{
  if (tb == NULL)
    return tb;

  BFT_FREE(tb->buf);
  BFT_FREE(tb);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute some simple statistics from an array
 *
 * \param[in]  n_elts      number of couples in data
 * \param[in]  stride      size of a couple of data
 * \param[in]  datatype    datatype
 * \param[in]  indata      buffer containing input data
 * \param[in]  do_abs      analyse the absolute value of indata
 *
 * \return a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

cs_data_info_t
cs_analysis_data(cs_lnum_t       n_elts,
                 int             stride,
                 cs_datatype_t   datatype,
                 const void     *indata,
                 _Bool           do_abs)
{
  int  i, j;

  _Bool  deallocate = false;
  cs_data_info_t  info = _init_dinfo(datatype);

  if (indata == NULL)
    return info;

  switch (datatype) { // Only for double
  case CS_DOUBLE:
    {
      const cs_real_t  *data = (const cs_real_t *)indata;
      cs_real_t  *values = NULL;

      if (stride == 3) { // compute norm

        cs_real_3_t  v;

        BFT_MALLOC(values, n_elts, cs_real_t);
        deallocate = true;
        for (i = 0; i < n_elts; i++) {
          for (j = 0; j < 3; j++)
            v[j] = data[stride*i+j];
          values[i] = cs_math_3_norm(v);
        }
        _compute_info_double(n_elts, values, &info);

      }
      else {  /* stride = 1 */

        assert(stride==1);
        if (do_abs) {
          BFT_MALLOC(values, n_elts, cs_real_t);
          deallocate = true;
          for (i = 0; i < n_elts; i++)
            values[i] = fabs(data[i]);
          _compute_info_double(n_elts, values, &info);
        }
        else
          _compute_info_double(n_elts, data, &info);

      } // stride

      if (deallocate)  BFT_FREE(values);
    }
    break;

  case CS_INT32:
    {
      const cs_lnum_t  *data = (const cs_lnum_t *)indata;
      cs_lnum_t  *numbers = NULL;

      if (do_abs) {
        BFT_MALLOC(numbers, n_elts, cs_lnum_t);
        for (i = 0; i < n_elts; i++)
          numbers[i] = CS_ABS(data[i]);
        _compute_info_int32(n_elts, numbers, &info);
        BFT_FREE(numbers);
      }
      else
        _compute_info_int32(n_elts, data, &info);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid datatype for analysing data.\n"));
  }

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_data_info_t structure
 *
 * \param[in]  name        filename if not NULL
 * \param[in]  f           output file if not NULL
 * \param[in]  n_elts      number of couples in data
 * \param[in]  datatype    datatype
 * \param[in]  dinfo       cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_data_info_dump(const char              *name,
                  FILE                    *f,
                  cs_lnum_t                n_elts,
                  cs_datatype_t            datatype,
                  const cs_data_info_t     dinfo)
{
  FILE  *_f = f;
  _Bool close_file = false;

  if (f == NULL) {
    if (name == NULL) _f = stdout;
    else
      _f = fopen(name, "w"), close_file = true;
  }

  fprintf(_f, "\n");
  if (name != NULL)
    fprintf(_f, " -dat- name          %s\n", name);
  fprintf(_f, " -dat- n_elts        %d\n", n_elts);

  switch (datatype) {
  case CS_DOUBLE:
    fprintf(_f, " -dat- minimum    %- 9.6e\n", dinfo.min.value);
    fprintf(_f, " -dat- maximum    %- 9.6e\n", dinfo.max.value);
    break;
  case CS_INT32:
    fprintf(_f, " -dat- minimum    %10d\n", dinfo.min.number);
    fprintf(_f, " -dat- maximum    %10d\n", dinfo.max.number);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid datatype for analysing data.\n"));
  }

  fprintf(_f, " -dat- mean            %- 9.6e\n", dinfo.mean);
  fprintf(_f, " -dat- sigma           %- 9.6e\n", dinfo.sigma);
  fprintf(_f, " -dat- euclidean norm  %- 9.6e\n", dinfo.euclidean_norm);
  fprintf(_f, "\n");

  fflush(_f);
  if (close_file)  fclose(_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create an index structure of size n
 *
 * \param[in]  n     number of entries of the indexed list
 *
 * \return  a pointer to a cs_connect_index_t
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_create(int  n)
{
  int  i;

  cs_connect_index_t  *x = NULL;

  BFT_MALLOC(x, 1, cs_connect_index_t);

  x->n = n;
  x->owner = true;
  x->ids = NULL;

  BFT_MALLOC(x->idx, n+1, int);
  for (i = 0; i < x->n + 1; i++)  x->idx[i] = 0;

  return x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Map arrays into an index structure of size n (owner = false)
 *
 * \param[in]  n     number of entries of the indexed list
 * \param[in]  idx   array of size n+1
 * \param[in]  ids   array of size idx[n]
 *
 * \return  a pointer to a cs_connect_index_t
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_map(int    n,
             int   *idx,
             int   *ids)
{
  cs_connect_index_t  *x = NULL;

  BFT_MALLOC(x, 1, cs_connect_index_t);

  x->n = n;
  x->owner = false;
  x->idx = idx;
  x->ids = ids;

  return x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a cs_connect_index_t structure
 *
 * \param[in]  pidx     pointer of pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_free(cs_connect_index_t   **pidx)
{
  cs_connect_index_t  *x = *pidx;

  if (x == NULL)
    return;

  if (x->owner) {
    BFT_FREE(x->idx);
    BFT_FREE(x->ids);
  }

  BFT_FREE(x);
  *pidx = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   From 2 indexes : A -> B and B -> C create a new index A -> C
 *
 * \param[in]  nc      number of elements in C set
 * \param[in]  a2b     pointer to the index A -> B
 * \param[in]  b2c     pointer to the index B -> C
 *
 *\return  a pointer to the cs_connect_index_t structure A -> C
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_compose(int                        nc,
                 const cs_connect_index_t  *a2b,
                 const cs_connect_index_t  *b2c)
{
  int  i, pos_a, pos_b, a_id, b_id, c_id, shift;

  int  *ctag = NULL;
  cs_connect_index_t  *a2c = cs_index_create(a2b->n);

  BFT_MALLOC(ctag, nc, int);
  for (i = 0; i < nc; i++)
    ctag[i] = -1;

  /* Build index */
  for (a_id = 0; a_id < a2b->n; a_id++) {

    for (pos_a = a2b->idx[a_id]; pos_a < a2b->idx[a_id+1]; pos_a++) {

      b_id = a2b->ids[pos_a];
      for (pos_b = b2c->idx[b_id]; pos_b < b2c->idx[b_id+1]; pos_b++) {

        c_id = b2c->ids[pos_b];
        if (ctag[c_id] != a_id) { /* Not tagged yet */
          ctag[c_id] = a_id;
          a2c->idx[a_id+1] += 1;
        }

      } /* End of loop on C elements */
    } /* End of loop on B elements */
  } /* End of loop on A elements */

  for (i = 0; i < a2c->n; i++)
    a2c->idx[i+1] += a2c->idx[i];

  BFT_MALLOC(a2c->ids, a2c->idx[a2c->n], int);

  /* Reset ctag */
  for (i = 0; i < nc; i++)
    ctag[i] = -1;

  /* Fill ids */
  shift = 0;
  for (a_id = 0; a_id < a2b->n; a_id++) {

    for (pos_a = a2b->idx[a_id]; pos_a < a2b->idx[a_id+1]; pos_a++) {

      b_id = a2b->ids[pos_a];
      for (pos_b = b2c->idx[b_id]; pos_b < b2c->idx[b_id+1]; pos_b++) {

        c_id = b2c->ids[pos_b];
        if (ctag[c_id] != a_id) { /* Not tagged yet */
          ctag[c_id] = a_id;
          a2c->ids[shift++] = c_id;
        }

      } /* End of loop on C elements */
    } /* End of loop on B elements */
  } /* End of loop on A elements */

  assert(shift == a2c->idx[a2c->n]);

  BFT_FREE(ctag);

  return a2c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   From a cs_connect_index_t A -> B create a new index B -> A
 *
 * \param[in]  nb     size of the "b" set
 * \param[in]  a2b    pointer to the index A -> B
 *
 * \return  a new pointer to the cs_connect_index_t structure B -> A
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_transpose(int                        nb,
                   const cs_connect_index_t  *a2b)
{
  int  i, j, b_id, shift;
  int  *count = NULL;

  cs_connect_index_t  *b2a = cs_index_create(nb);

  if (nb == 0)
    return b2a;

  /* Build idx */
  for (i = 0; i < a2b->n; i++)
    for (j = a2b->idx[i]; j < a2b->idx[i+1]; j++)
      b2a->idx[a2b->ids[j]+1] += 1;

  for (i = 0; i < b2a->n; i++)
    b2a->idx[i+1] += b2a->idx[i];

  /* Allocate and initialize temporary buffer */
  BFT_MALLOC(count, nb, int);
  for (i = 0; i < nb; i++) count[i] = 0;

  /* Build ids */
  BFT_MALLOC(b2a->ids, b2a->idx[b2a->n], int);

  for (i = 0; i < a2b->n; i++) {
    for (j = a2b->idx[i]; j < a2b->idx[i+1]; j++) {
      b_id = a2b->ids[j];
      shift = count[b_id] + b2a->idx[b_id];
      b2a->ids[shift] = i;
      count[b_id] += 1;
    }
  }

  /* Free temporary buffer */
  BFT_FREE(count);

  return b2a;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Sort each sub-list related to an entry in a cs_connect_index_t
 *          structure
 *
 * \param[in]  x     pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_sort(cs_connect_index_t   *x)
{
  if (x == NULL)
    return;

  for (int i = 0; i < x->n; i++)
    cs_sort_shell(x->idx[i], x->idx[i+1], x->ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_connect_index_t structure to a file or into the
 *          standard output
 *
 * \param[in]  name  name of the dump file. Can be set to NULL
 * \param[in]  _f    pointer to a FILE structure. Can be set to NULL.
 * \param[in]  x     pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_dump(const char           *name,
              FILE                 *_f,
              cs_connect_index_t   *x)
{
  FILE  *f = _f;
  _Bool  close_file = false;

  if (f == NULL) {
    if (name == NULL)
      f = stdout;
    else {
      f = fopen(name,"w");
      close_file = true;
    }
  }

  fprintf(f, "\n Dump cs_connect_index_t struct: %p (%s)\n",
          (const void *)x, name);

  if (x == NULL) {
    if (close_file) fclose(f);
    return;
  }

  fprintf(f, "  owner:             %6d\n", x->owner);
  fprintf(f, "  n_elts:            %6d\n", x->n);
  fprintf(f, "  ids_size:          %6d\n", x->idx[x->n]);

  for (int i = 0; i < x->n; i++) {
    fprintf(f, "\n[%4d] ", i);
    for (int j = x->idx[i]; j < x->idx[i+1]; j++)
      fprintf(f, "%5d |", x->ids[j]);
  }

  if (close_file)
    fclose(f);
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
                   double              alpha,
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

  bft_printf("\n  << parent id: %d >>\n", parent_id);

  /* List sub-entity ids */
  for (i = 0; i < lm->n_ent; i++)
    bft_printf(" %9d", lm->ids[i]);
  bft_printf("\n");

  for (i = 0; i < lm->n_ent; i++) {
    bft_printf(" %5d", lm->ids[i]);
    for (j = 0; j < lm->n_ent; j++)
      bft_printf(" % .4e", lm->val[i*lm->n_ent+j]);
    bft_printf("\n");
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_locdec_t structure
 *
 * \param[in]  n_max_rows   max number of rows
 * \param[in]  n_max_cols   max number of columns
 *
 * \return  a new allocated cs_locdec_t structure
 */
/*----------------------------------------------------------------------------*/

cs_locdec_t *
cs_locdec_create(int   n_max_rows,
                 int   n_max_cols)
{
  int  i;

  cs_locdec_t  *m = NULL;

  BFT_MALLOC(m, 1, cs_locdec_t);

  m->n_max_rows = n_max_rows;
  m->n_max_cols = n_max_cols;
  m->n_rows = 0;
  m->n_cols = 0;
  m->row_ids = m->col_ids = NULL;
  m->sgn = NULL;

  int msize = n_max_rows * n_max_cols;

  if (msize > 0) {

    BFT_MALLOC(m->row_ids, n_max_rows, cs_lnum_t);
    for (i = 0; i < n_max_rows; i++)
      m->row_ids[i] = 0;
    BFT_MALLOC(m->col_ids, n_max_cols, cs_lnum_t);
    for (i = 0; i < n_max_cols; i++)
      m->col_ids[i] = 0;

    BFT_MALLOC(m->sgn, msize, short int);
    for (i = 0; i < msize; i++)
      m->sgn[i] = 0;

  }

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_locdec_t structure
 *
 * \param[in, out]  m    pointer to a cs_locdec_t structure to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_locdec_t *
cs_locdec_free(cs_locdec_t  *m)
{
  if (m == NULL)
    return m;

  if (m->n_max_rows > 0 && m->n_max_cols > 0) {
    BFT_FREE(m->col_ids);
    BFT_FREE(m->row_ids);
    BFT_FREE(m->sgn);
  }

  BFT_FREE(m);

  return NULL;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
