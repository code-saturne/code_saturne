/*============================================================================
 * Routines to handle common features for building algebraic system in CDO
 * schemes
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

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_local.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdofb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_common.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

static size_t  cs_equation_common_work_buffer_size = 0;
static cs_real_t  *cs_equation_common_work_buffer = NULL;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 * \param[in]  quant         pointer to additional mesh quantities struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_allocate_common_structures(const cs_cdo_connect_t     *connect,
                                       const cs_cdo_quantities_t  *quant,
                                       const cs_time_step_t       *time_step,
                                       cs_flag_t                   scheme_flag)
{
  assert(connect != NULL); // Sanity check

  /* Allocate cell-wise and face_wise view of a mesh */
  cs_cdo_local_initialize(connect);

  const cs_lnum_t  n_cells = connect->c_info->n_elts;
  const cs_lnum_t  n_faces = connect->f_info->n_elts;
  const cs_lnum_t  n_vertices = connect->v_info->n_elts;

  size_t  cwb_size = 2*n_cells;

  if (scheme_flag & CS_SCHEME_FLAG_CDOVB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    cwb_size = CS_MAX(cwb_size, (size_t)3*n_vertices);

    /* Initialize additional structures */
    cs_cdovb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdovb_scaleq_initialize();

  }

  if (scheme_flag & CS_SCHEME_FLAG_CDOVCB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    cwb_size = CS_MAX(cwb_size, (size_t)2*(n_vertices + n_cells));

    /* Initialize additional structures */
    cs_cdovcb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdovcb_scaleq_initialize();
  }

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    /* Initialize additional structures */
    cs_cdofb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdofb_scaleq_initialize();

    cwb_size = CS_MAX(cwb_size, (size_t)3*n_faces);

  }

  cs_equation_common_work_buffer_size = cwb_size;
  BFT_MALLOC(cs_equation_common_work_buffer, cwb_size, double);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_free_common_structures(cs_flag_t   scheme_flag)
{
  /* Free cell-wise and face-wise view of a mesh */
  cs_cdo_local_finalize();

  /* Free common structures specific to a numerical scheme */
  if ((scheme_flag & CS_SCHEME_FLAG_CDOVB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdovb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_CDOVCB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdovcb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_CDOFB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdofb_scaleq_finalize();

  BFT_FREE(cs_equation_common_work_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the 2*n_cells
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_tmpbuf(void)
{
  return cs_equation_common_work_buffer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_equation_get_tmpbuf_size(void)
{
  return cs_equation_common_work_buffer_size;
}
