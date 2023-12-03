/*============================================================================
 * Solid zones handling.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_flag_check.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_solid_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_solid_zone.c
        Volume zone handling.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Zone definitions */

static int  _n_solid_zones = -1;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check (and cache) number of solid zones
 */
/*----------------------------------------------------------------------------*/

static int
_solid_zone_n_zones(void)
{
  if (_n_solid_zones != cs_volume_zone_n_zones())
    _n_solid_zones = cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_SOLID);

  return _n_solid_zones;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief build solid flag for mesh cells.
 *
 * If no solid cells are present, NULL is returned.

 * If non-empty, the caller is responsible for freeing the flag
 *
 * \param[in]   m         pointer to mesh
 *
 * \return  solid cell flag array (0 is fluid, 1 if solid), or NULL
 */
/*----------------------------------------------------------------------------*/

int *
cs_solid_zone_flag(const cs_mesh_t   *m)
{
  int *c_is_solid = NULL;

  if (_solid_zone_n_zones() == 0)
    return c_is_solid;

  cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  BFT_MALLOC(c_is_solid, n_cells_ext, int);
  for (cs_lnum_t i = 0; i < n_cells_ext; i++)
    c_is_solid[i] = 0;
  cs_volume_zone_tag_cell_type(CS_VOLUME_ZONE_SOLID, 1, c_is_solid);

  return c_is_solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Zero an array on cells of a solid zone
 *
 * \param[in]   stride  array stride
 * \param[out]  a       array of cell values
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_zone_set_zero_on_cells(int         stride,
                                cs_real_t  *a)
{
  if (_solid_zone_n_zones() == 0)
    return;

  int n_zones = cs_volume_zone_n_zones();

  for (int i = 0; i < n_zones; i++) {
    const cs_zone_t *z = cs_volume_zone_by_id(i);

    if (z->type & CS_VOLUME_ZONE_SOLID) {

      const cs_lnum_t _n_cells = z->n_elts;
      const cs_lnum_t *_cell_ids = z->elt_ids;
      if (_cell_ids != NULL) {
        if (stride == 1) {
          for (cs_lnum_t j = 0; j < _n_cells; j++)
            a[_cell_ids[j]] = 0;
        }
        else {
          for (cs_lnum_t j = 0; j < _n_cells; j++) {
            cs_lnum_t k = _cell_ids[j];
            for (cs_lnum_t l = 0; l < stride; l++)
              a[k*stride + l] = 0;
          }
        }
      }
      else {
        if (stride == 1) {
          for (cs_lnum_t j = 0; j < _n_cells; j++)
            a[j] = 0;
        }
        else {
          for (cs_lnum_t j = 0; j < _n_cells; j++) {
            for (cs_lnum_t l = 0; l < stride; l++)
              a[l*stride + l] = 0;
          }
        }
      }

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant scalar value to cells of a solid zone.
 *
 * \param[in]   ref_val  reference value
 * \param[out]  a        array of cell values
 */
/*----------------------------------------------------------------------------*/

void
cs_solid_zone_set_scalar_on_cells(cs_real_t  ref_val,
                                  cs_real_t  a[])
{
  if (_solid_zone_n_zones() == 0)
    return;

  int n_zones = cs_volume_zone_n_zones();

  for (int i = 0; i < n_zones; i++) {
    const cs_zone_t *z = cs_volume_zone_by_id(i);

    if (z->type & CS_VOLUME_ZONE_SOLID) {

      const cs_lnum_t _n_cells = z->n_elts;
      const cs_lnum_t *_cell_ids = z->elt_ids;
      if (_cell_ids != NULL) {
        for (cs_lnum_t j = 0; j < _n_cells; j++)
          a[_cell_ids[j]] = ref_val;
      }
      else {
        for (cs_lnum_t j = 0; j < _n_cells; j++)
          a[j] = ref_val;
      }

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
