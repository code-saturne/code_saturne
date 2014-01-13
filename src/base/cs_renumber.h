#ifndef __CS_RENUMBER_H__
#define __CS_RENUMBER_H__

/*============================================================================
 * Optional mesh renumbering
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Renumbering algorithm */

typedef enum {

  CS_RENUMBER_I_FACES_BLOCK,      /* No shared cell in block */
  CS_RENUMBER_I_FACES_MULTIPASS,  /* Use multipass face numbering */
  CS_RENUMBER_I_FACES_NONE        /* No interior face numbering */

} cs_renumber_i_faces_type_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set the target number of threads for mesh renumbering.
 *
 * By default, the target number of threads is set to cs_glob_n_threads,
 * but the value may be forced using this function. This is mainly useful
 * for testing purposes.
 *
 * parameters:
 *   n_threads <-- target number of threads for mesh numbering
 *----------------------------------------------------------------------------*/

void
cs_renumber_set_n_threads(int  n_threads);

/*----------------------------------------------------------------------------
 * Return the target number of threads for mesh renumbering.
 *
 * returns:
 *   the target number of threads for mesh numbering
 *----------------------------------------------------------------------------*/

int
cs_renumber_get_n_threads(void);

/*----------------------------------------------------------------------------
 * Set the minimum sunset sizes when renumbering for threads.
 *
 * parameters:
 *   min_i_subset_size <-- minimum number of interior faces per
 *                         thread per group
 *   min_b_subset_size <-- minimum number of boundary faces per
 *                         thread per group
 *----------------------------------------------------------------------------*/

void
cs_renumber_set_min_subset_size(cs_lnum_t  min_i_subset_size,
                                cs_lnum_t  min_b_subset_size);

/*----------------------------------------------------------------------------
 * Get the minimum sunset sizes when renumbering for threads.
 *
 *   min_i_subset_size --> minimum number of interior faces per
 *                         thread per group, or NULL
 *   min_b_subset_size --> minimum number of boundary faces per
 *                         thread per group, or NULL
 *----------------------------------------------------------------------------*/

void
cs_renumber_get_min_subset_size(cs_lnum_t  *min_i_subset_size,
                                cs_lnum_t  *min_b_subset_size);

/*----------------------------------------------------------------------------
 * Select the algorithm for interior faces renumbering.
 *
 * parameters:
 *   algorithm <-- algorithm type for interior faces renumbering
 *----------------------------------------------------------------------------*/

void
cs_renumber_set_i_face_algorithm(cs_renumber_i_faces_type_t  algorithm);

/*----------------------------------------------------------------------------
 * Return the algorithm for interior faces renumbering.
 *
 * returns:
 *   algorithm type for interior faces renumbering
 *----------------------------------------------------------------------------*/

cs_renumber_i_faces_type_t
cs_renumber_get_i_face_algorithm(void);

/*----------------------------------------------------------------------------
 * Renumber mesh elements for vectorization or OpenMP depending on code
 * options and target machine.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_mesh(cs_mesh_t             *mesh,
                 cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RENUMBER_H__ */
