#ifndef __CS_MEDCOUPLING_UTILS_HXX__
#define __CS_MEDCOUPLING_UTILS_HXX__

/*============================================================================
 * Usage of MEDCoupling base components.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#if defined(HAVE_MEDCOUPLING)

#include "MEDCouplingUMesh.hxx"

using namespace MEDCoupling;

#endif


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
  cs_lnum_t          *new_to_old;     /* Connectivity used if only a section of
                                         the mesh is read */

  cs_real_t          *bbox;           /* Bounding box to optimize search */

#if defined(HAVE_MEDCOUPLING)
  MEDCouplingUMesh   *med_mesh;       /* MED mesh structure */
#else
  void               *med_mesh;
#endif

} cs_medcoupling_mesh_t;

/*=============================================================================
 * Public functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   create a new cs_medcoupling_mesh_t instance
 *
 * \param[in] name                name of the mesh
 * \param[in] selection_criteria  selection criteria (entire mesh or part of it)
 * \param[in] elt_dim             dimension of elements.
 *                                2: faces
 *                                3: cells
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_create(const char  *name,
                           const char  *selection_criteria,
                           int          elt_dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief copy a cs_mesh_t into a cs_medcoupling_mesh_t
 *
 * \param[in] csmesh    pointer to the cs_mesh_t struct to copy data from
 * \param[in] pmmesh    pointer to the cs_medcoupling_mesh_t for copy
 * \param[in] use_bbox  flag indicating if a reduced bounding is used. Usefull
 *                      for interpolation to reduce the matrix sice.
 *                      0: Do not use a reduced bbox
 *                      1: Use a reduced bbox
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_mesh_copy_from_base(cs_mesh_t              *csmesh,
                                   cs_medcoupling_mesh_t  *pmmesh,
                                   int                     use_bbox);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_medcoupling_mesh_t
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_mesh_destroy(cs_medcoupling_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's spatial dimension
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return associated spatial dimension
 */
/*----------------------------------------------------------------------------*/

int
cs_medcoupling_mesh_get_dim(cs_medcoupling_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's number of elements
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return associated number of elements
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_medcoupling_mesh_get_n_elts(cs_medcoupling_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's (parent) elements list
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return ids of associated elements, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_medcoupling_mesh_get_elt_list(cs_medcoupling_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's (parent) elements list
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return ids of associated elements, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_medcoupling_mesh_get_connectivity(cs_medcoupling_mesh_t  *m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MEDCOUPLING_UTILS_HXX__ */
