#ifndef __CS_CDO_TOOLBOX_H__
#define __CS_CDO_TOOLBOX_H__

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure for managing temporary buffers */
typedef struct {

  size_t    bufsize;
  void     *buf;

} cs_tmpbuf_t;

/* Index to scan a mesh connectivity */
typedef struct {

  bool     owner;
  int      n;
  int     *idx;   /* from 0, size = n+1 */
  int     *ids;   /* ids from 0 to n-1 (there is no multifold entry) */

} cs_connect_index_t;

/* Structure enabling the repeated usage of local dense square matrix
   associated  with a local set of entities */
typedef struct {

  int         n_max_ent;  // max number of entities by primal cells
  int         n_ent;      // current number of entities
  cs_lnum_t  *ids;        // list of entity ids (size = n_max_ent)
  double     *val;        // local matrix (size: n_max_ent*n_max_ent)

} cs_locmat_t;

/* Structure enabling the repeated usage of local dense DEC matrix associated
   with a local set of entities */
typedef struct {

  int         n_max_rows; // max number of entries in a row
  int         n_rows;     // current number of rows
  cs_lnum_t  *row_ids;    // list of entity ids in a row (size = n_max_rows)

  int         n_max_cols; // max number of entries in a column
  int         n_cols;     // current number of columns
  cs_lnum_t  *col_ids;    // list of entity ids in a col (size = n_max_rows)

  short int  *sgn;        // local matrix (size: n_max_rows*n_max_cols)

} cs_locdec_t;

/* ============= *
 * DATA ANALYSIS *
 * ============= */

typedef union {

  cs_lnum_t   number;
  double      value;

} cs_data_t;

typedef struct {

  cs_data_t    min;
  cs_data_t    max;
  double       mean;
  double       sigma;
  double       euclidean_norm;

} cs_data_info_t;

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

double
cs_euclidean_norm(int            len,
                  const double   v[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted sum of square values of an array
 *
 *  \param[in] n      size of arrays v and w
 *  \param[in] x      array of floating-point values
 *  \param[in] w      floating-point values of weights
 *
 * \return the result of this operation
 */
/*----------------------------------------------------------------------------*/

double
cs_weighted_sum_square(cs_lnum_t                n,
                       const double *restrict   x,
                       const double *restrict   weight);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate or reallocate a temporary buffer structure
 *
 * \param[in]       bufsize   reference size
 * \param[in, out]  p_tb      pointer to the temporary structure to allocate
 */
/*----------------------------------------------------------------------------*/

void
cs_tmpbuf_alloc(size_t          bufsize,
                cs_tmpbuf_t   **p_tb);

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
cs_tmpbuf_free(cs_tmpbuf_t    *tb);

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
                 bool            do_abs);

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
                  const cs_data_info_t     dinfo);

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
cs_index_create(int  n);

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
             int   *ids);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a cs_connect_index_t structure
 *
 * \param[in]  pidx     pointer of pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_free(cs_connect_index_t   **pidx);

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
                 const cs_connect_index_t  *b2c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   From a cs_connect_index_t struct. A -> B create a new index B -> A
 *
 * \param[in]  nb     size of the "b" set
 * \param[in]  a2b    pointer to the index A -> B
 *
 *\return  a new pointer to the cs_connect_index_t structure B -> A
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_transpose(int                         nb,
                   const cs_connect_index_t   *a2b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Sort each list related to an entry in a cs_connect_index_t structure
 *
 * \param[in]  x     pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_sort(cs_connect_index_t   *x);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_connect_index_t structure to a file or into the standard
 *          output
 *
 * \param[in]  name  name of the dump file. Can be set to NULL
 * \param[in]  _f    pointer to a FILE structure. Can be set to NULL.
 * \param[in]  x     pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_dump(const char          *name,
              FILE                *_f,
              cs_connect_index_t  *x);

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
cs_locmat_create(int   n_max_ent);

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
cs_locmat_free(cs_locmat_t  *lm);

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
               const cs_locmat_t  *send);

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
                 cs_real_t           *matvec);

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
              const cs_locmat_t  *add);

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
                   const cs_locmat_t  *add);

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
                        cs_locmat_t  *tr);

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
               const cs_locmat_t  *lm);

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
                 int   n_max_cols);

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
cs_locdec_free(cs_locdec_t  *m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TOOLBOX_H__ */
