/*============================================================================
 * General functions or variables for the CDO module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

#include "cs_log.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_dbg.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Local static variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(DEBUG) && !defined(NDEBUG)
/*----------------------------------------------------------------------------*/
/*!
 * \brief   Function used to select which element deserves a dump or specific
 *          treatment during a debugging stage
 *
 * \param[in]  eqp      pointer to a cs_equation_param_t structure
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  csys     pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

bool
cs_dbg_cw_test(const cs_equation_param_t   *eqp,
               const cs_cell_mesh_t        *cm,
               const cs_cell_sys_t         *csys)
{
#if 1 /* Example: Only search debug information for a givan equation */
  bool has_name = false;
  if (eqp != NULL) {
    if (strcmp(eqp->name, "Tracer1") == 0)
      has_name=true;
  }
#else
  bool has_name = true;
#endif

  if (has_name) {
#if 1 /* First example: Look for the cells which have the vertex 441 */
    short int _v = -1;
    for (int v = 0; v < cm->n_vc; v++)
      if (cm->v_ids[v] == 441)
        _v = v;
#else
    short int _v = 0;
#endif

#if 1 /* Second example: Look for the cells which have a previous DoF value
         greater than 0.06 and the vertex 441 */
    if (csys != NULL && _v > -1) {
      if (csys->val_n[_v] > 0.06)
        return true;
    }
#endif
  } /* The current equation has the requested name */

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a cs_sdm_t structure which is defined by block
 *          Print into the file f if given otherwise open a new file named
 *          fname if given otherwise print into the standard output
 *          The usage of threshold allows one to compare more easier matrices
 *          without taking into account numerical roundoff.
 *
 * \param[in]  fp         pointer to a file structure or NULL
 * \param[in]  fname      filename or NULL
 * \param[in]  thd        threshold (below this value --> set 0)
 * \param[in]  n_elts     size of the array
 * \param[in]  array      list of values to dump
 * \param[in]  n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_array_fprintf(FILE             *fp,
                     const char       *fname,
                     cs_real_t         thd,
                     cs_lnum_t         n_elts,
                     const cs_real_t   array[],
                     int               n_cols)
{
  FILE  *fout = stdout;
  if (fp != NULL)
    fout = fp;
  else if (fname != NULL) {
    fout = fopen(fname, "w");
  }

  fprintf(fout, "array %p\n", (const void *)array);

  if (array == NULL)
    return;

  if (n_cols < 1) n_cols = 1;
  int  n_rows = n_elts/n_cols;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    for (cs_lnum_t j = i*n_cols; j < (i+1)*n_cols; j++) {
      if (fabs(array[j]) < thd)
        fprintf(fout, "% -8.5e", 0.);
      else
        fprintf(fout, "% -8.5e", array[j]);
    }
    fprintf(fout, "\n");
  }

  if (n_rows*n_cols < n_elts) {
    for (cs_lnum_t j = n_rows*n_cols; j < n_elts; j++) {
      if (fabs(array[j]) < thd)
        fprintf(fout, "% -8.5e", 0.);
      else
        fprintf(fout, "% -8.5e", array[j]);
    }
    fprintf(fout, "\n");
  }

  if (fout != stdout && fout != fp)
    fclose(fout);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In debug mode, dump an array of double into the listing
 *
 * \param[in] header     header message to write
 * \param[in] size       number of elements in array
 * \param[in] array      pointer to the array of values
 * \param[in] n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_darray_to_listing(const char        *header,
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
 * \brief  In debug mode, dump an array of integer into the listing
 *
 * \param[in] header     header message to write
 * \param[in] size       number of elements in array
 * \param[in] array      pointer to the array of values
 * \param[in] n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_iarray_to_listing(const char        *header,
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
/*!
 * \brief  In debug mode, dump a linear system
 *
 * \param[in] eqname     name of the equation related to the current system
 * \param[in] size       number of elements in array
 * \param[in] x          solution array
 * \param[in] b          right-hand side
 * \param[in] row_index  index on row entries (column id and extra-diag values)
 * \param[in] col_id     list of column id
 * \param[in] xval       array of extra-diagonal values
 * \param[in] dval       array of diagonal values
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_dump_linear_system(const char        *eqname,
                          cs_lnum_t          size,
                          int                verbosity,
                          const cs_real_t    x[],
                          const cs_real_t    b[],
                          const cs_lnum_t    row_index[],
                          const cs_lnum_t    col_id[],
                          const cs_real_t    xval[],
                          const cs_real_t    dval[])
{
  cs_log_printf(CS_LOG_DEFAULT, "\nDUMP LINEAR SYSTEM FOR THE EQUATION >> %s\n",
                eqname);

  int  n_dump_cols = 8;
  int  n_dump_rows = size/n_dump_cols;

  cs_log_printf(CS_LOG_DEFAULT, " >> SOLUTION\n");
  for (cs_lnum_t i = 0; i < n_dump_rows; i++) {
    for (cs_lnum_t j = i*n_dump_cols; j < (i+1)*n_dump_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, "%4d % -6.4e |", j, x[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }
  if (n_dump_rows*n_dump_cols < size) {
    for (cs_lnum_t j = n_dump_rows*n_dump_cols; j < size; j++)
      cs_log_printf(CS_LOG_DEFAULT, "%4d % -6.4e |", j, x[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

  cs_log_printf(CS_LOG_DEFAULT, " >> RIGHT-HAND SIDE\n");
  for (cs_lnum_t i = 0; i < n_dump_rows; i++) {
    for (cs_lnum_t j = i*n_dump_cols; j < (i+1)*n_dump_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, "%4d % -6.4e |", j, b[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }
  if (n_dump_rows*n_dump_cols < size) {
    for (cs_lnum_t j = n_dump_rows*n_dump_cols; j < size; j++)
      cs_log_printf(CS_LOG_DEFAULT, "%4d % -6.4e |", j, b[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

  if (verbosity > 2) {

    cs_log_printf(CS_LOG_DEFAULT, " >> DIAGONAL ENTRIES\n");
    for (cs_lnum_t i = 0; i < n_dump_rows; i++) {
      for (cs_lnum_t j = i*n_dump_cols; j < (i+1)*n_dump_cols; j++)
        cs_log_printf(CS_LOG_DEFAULT, "%4d % -6.4e |", j, dval[j]);
      cs_log_printf(CS_LOG_DEFAULT, "\n");
    }
    if (n_dump_rows*n_dump_cols < size) {
      for (cs_lnum_t j = n_dump_rows*n_dump_cols; j < size; j++)
        cs_log_printf(CS_LOG_DEFAULT, "%4d % -6.4e |", j, dval[j]);
      cs_log_printf(CS_LOG_DEFAULT, "\n");
    }

    cs_log_printf(CS_LOG_DEFAULT, " >> EXTRA-DIAGONAL ENTRIES\n");
    for (cs_lnum_t i = 0; i < size; i++) {

      const cs_lnum_t  *idx = row_index + i;
      const cs_lnum_t  *_col = col_id + idx[0];
      const cs_real_t  *_val = xval + idx[0];
      const int  n_entries = idx[1] - idx[0];

      int  _n_cols = 6;
      int  _n_rows = n_entries/_n_cols;

      for (cs_lnum_t ii = 0; ii < _n_rows; ii++) {
        cs_log_printf(CS_LOG_DEFAULT, "ROW%4d >> ", i);
        for (cs_lnum_t jj = ii*_n_cols; jj < (ii+1)*_n_cols; jj++)
          cs_log_printf(CS_LOG_DEFAULT, "%4d: % -6.4e |", _col[jj], _val[jj]);
        cs_log_printf(CS_LOG_DEFAULT, "\n");
      }
      if (_n_rows*_n_cols < n_entries) {
        cs_log_printf(CS_LOG_DEFAULT, "ROW%4d >> ", i);
        for (cs_lnum_t jj = _n_rows*_n_cols; jj < n_entries; jj++)
          cs_log_printf(CS_LOG_DEFAULT, "%4d: % -6.4e |", _col[jj], _val[jj]);
        cs_log_printf(CS_LOG_DEFAULT, "\n");
      }

    }

  } /* verbosity */
}
#endif  /* Only in debug mode */

/*----------------------------------------------------------------------------*/

END_C_DECLS
