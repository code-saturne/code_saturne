#ifndef __CS_MEDCOUPLING_UTILS_HXX__
#define __CS_MEDCOUPLING_UTILS_HXX__

/*============================================================================
 * Usage of MEDCoupling base components.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#include <MEDCoupling_version.h>
#include <MEDCouplingUMesh.hxx>

using namespace MEDCoupling;

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * MEDCoupling mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char               *sel_criteria;   /* Element selection criteria */

  int                 elt_dim;        /* Element dimension */

  cs_lnum_t           n_elts;         /* Number of coupled elements */
  cs_lnum_t          *elt_list;       /* List of associated elements
                                         (0 to n-1) */
  int                *new_to_old;     /* Connectivity used if only a section of
                                         the mesh is read */

  MEDCouplingUMesh   *med_mesh;       /* MED mesh structure */

  cs_real_t          *bbox;           /* Bounding box to optimize search */

} cs_medcoupling_mesh_t;

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a cs_medcoupling_mesh_t structure.
 *
 * parameters:
 *   name            <-- mesh name
 *   select_criteria <-- Selection criteria
 *   elt_dim         <-- Element type (3: cells, 2: faces)
 *----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_create(const char  *name,
                           const char  *select_criteria,
                           int          elt_dim);

/*----------------------------------------------------------------------------
 * Copy a base mesh to a medcoupling mesh structure.
 *
 * parameters:
 *   csmesh  <-- Code_Saturne FVM format mesh structure
 *   pmmesh  <-> partially ParaMEDMEM mesh coupling structure
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_mesh_copy_from_base(cs_mesh_t              *csmesh,
                                   cs_medcoupling_mesh_t  *pmmesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif //__CS_MEDCOUPLING_UTILS_HXX__
