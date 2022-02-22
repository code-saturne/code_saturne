/*============================================================================
 * Manage the list of solid cells and associated helper functions
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
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solid_selection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_solid_selection.c

  \brief Structure and functions handling the list of solid cells
         Useful for Navier-Stokes, thermal module or the solidification module

*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_SOLID_SELECTION_DBG       0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static variables
 *============================================================================*/

cs_solid_selection_t  *_cs_solid = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the information related to the list of solid cells
 *         If this structure does not exist, there is an initialization.
 *
 * \return a pointer to a cs_solid_selection_t structure
 */
/*----------------------------------------------------------------------------*/

cs_solid_selection_t *
cs_solid_selection_get(void)
{
  if (_cs_solid == NULL) {

    BFT_MALLOC(_cs_solid, 1, cs_solid_selection_t);

    _cs_solid->n_cells = 0;
    _cs_solid->n_g_cells = 0;
    _cs_solid->cell_ids = NULL;

    _cs_solid->cell_is_solid = NULL;
    _cs_solid->face_is_solid = NULL;

  }

  return _cs_solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the solid selection
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_selection_sync(const cs_cdo_connect_t   *connect)
{
  if (_cs_solid == NULL)
    cs_solid_selection_get();

  /* Parallel synchronization of the global number of solid cells */

  _cs_solid->n_g_cells = _cs_solid->n_cells;
  cs_parall_sum(1, CS_GNUM_TYPE, &_cs_solid->n_g_cells);

  if (_cs_solid->n_g_cells > 0) {

    /* Tag cells */

    if (_cs_solid->cell_is_solid == NULL)
      BFT_MALLOC(_cs_solid->cell_is_solid, connect->n_cells, bool);

    /* Set to false all cells */

    for (cs_lnum_t i = 0; i < connect->n_cells; i++)
      _cs_solid->cell_is_solid[i] = false;

    /* Tag cells associated to a solid cell */

    for (cs_lnum_t i = 0; i < _cs_solid->n_cells; i++)
      _cs_solid->cell_is_solid[_cs_solid->cell_ids[i]] = true;

    /* Tag faces */

    if (_cs_solid->face_is_solid == NULL)
      BFT_MALLOC(_cs_solid->face_is_solid, connect->n_faces[0], bool);

    /* Set to false all faces */

    for (cs_lnum_t i = 0; i < connect->n_faces[0]; i++)
      _cs_solid->face_is_solid[i] = false;

    /* Tag faces associated to a solid cell */

    const cs_adjacency_t  *c2f = connect->c2f;

    for (cs_lnum_t i = 0; i < _cs_solid->n_cells; i++) {

      const cs_lnum_t  c_id = _cs_solid->cell_ids[i];
      for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
        _cs_solid->face_is_solid[c2f->ids[j]] = true;

    } /* Loop on solid cells */

    /* Synchronize face tags */

    if (connect->face_ifs != NULL) {

      assert(sizeof(bool) == sizeof(char));
      cs_interface_set_max(connect->face_ifs,
                           connect->n_faces[0],
                           1,             /* stride */
                           false,         /* interlace (not useful here) */
                           CS_CHAR,       /* boolean */
                           _cs_solid->face_is_solid);

    }

  } /* n_g_cells > 0 */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the structure storing the information related to the list of
 *         solid cells.
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_selection_free(void)
{
  if (_cs_solid == NULL)
    return;

  BFT_FREE(_cs_solid->cell_is_solid);
  BFT_FREE(_cs_solid->face_is_solid);
  BFT_FREE(_cs_solid->cell_ids);
  BFT_FREE(_cs_solid);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
