#ifndef __CS_MEDCOUPLING_REMAPPER_HXX__
#define __CS_MEDCOUPLING_REMAPPER_HXX__

/*============================================================================
 * Interpolation using MEDCoupling Remapper.
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

/*----------------------------------------------------------------------------
 * MED library headers
 *----------------------------------------------------------------------------*/

BEGIN_C_DECLS

#if defined(HAVE_MEDCOUPLING_LOADER)

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _cs_medcoupling_remapper_t cs_medcoupling_remapper_t;

/*============================================================================
 * Public C++ function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return remapper associated with a given id
 *
 * parameters:
 *   id <-- remapper id
 *
 * return:
 *   pointer to remapper
 *----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_id(int  r_id);

/*----------------------------------------------------------------------------
 * Return remapper associated with a given name
 *
 * parameters:
 *   name <-- remapper name
 *
 * return:
 *   pointer to remapper, or NULL
 *----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_name_try(const char  *name);

#if defined(HAVE_MEDCOUPLING_LOADER)

/*----------------------------------------------------------------------------
 * Create or update update the list of remappers in the case where
 * several remappers may be needed.
 *
 * parameters:
 *   name            <-- new remapper name
 *   elt_dim         <-- element dimension
 *   select_criteria <-- selection criteria
 *   medfile_path    <-- path of associated MED file
 *   mesh_name       <-- mesh name
 *   n_fields        <-- number of fields
 *   field_names     <-- associated field names
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *
 * return:
 *   id of the newly added remapper within the list
 *----------------------------------------------------------------------------*/

int
cs_medcoupling_remapper_initialize(const char   *name,
                                   int           elt_dim,
                                   const char   *select_criteria,
                                   const char   *medfile_path,
                                   const char   *mesh_name,
                                   int           n_fields,
                                   const char  **field_names,
                                   int           iteration,
                                   int           iteration_order);

#endif /* HAVE_MEDCOUPLING_LOADER */

/*----------------------------------------------------------------------------
 * Update field values (if several time steps are available in the MED file).
 *
 * parameters:
 *   r               <-- remapper object
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_set_iteration(cs_medcoupling_remapper_t  *r,
                                      int                         iteration,
                                      int                         iteration_order);

/*----------------------------------------------------------------------------
 * Create the interpolation matrix.
 *
 * This step is separated from the interpolation step since it only needs
 * to be done once per mesh, while interpolation can be done for several
 * fields.
 *
 * parameters:
 *   r               <-- remapper object
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_setup(cs_medcoupling_remapper_t  *r);

/*----------------------------------------------------------------------------
 * Copy interpolated values to a new array.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   field_id        <-- id of given field
 *   r               <-- pointer to remapper object
 *   default_val     <-- default value
 *
 * return:
 *   pointer to allocated values array
 *----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_remapper_copy_values(cs_medcoupling_remapper_t  *r,
                                    int                         field_id,
                                    double                      default_val);

/*----------------------------------------------------------------------------
 * Translate the mapped source mesh.
 *
 * Caution: cs_medcoupling_remapper_prepare() must to be called after this
 * function in order to update the interpolation matrix.
 *
 * parameters:
 *   r           <-- pointer to remapper object
 *   translation <-- translation vector
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_translate(cs_medcoupling_remapper_t  *r,
                                  cs_real_t                   translation[3]);

/*----------------------------------------------------------------------------
 * Rotate the mapped source mesh.
 *
 * Caution: cs_medcoupling_remapper_prepare() must to be called after this
 * function in order to update the interpolation matrix.
 *
 * parameters:
 *   r         <-- pointer to remapper object
 *   invariant <-- coordinates of invariant point
 *   axis      <-- rotation axis vector
 *   angle     <-- rotation angle
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_rotate(cs_medcoupling_remapper_t  *r,
                               cs_real_t                   invariant[3],
                               cs_real_t                   axis[3],
                               cs_real_t                   angle);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* MEDCOUPLING_LOADER */

#endif /* __CS_MEDCOUPLING_REMAPPER_HXX__ */
