#ifndef __CS_MEDCOUPLING_MESH_HXX__
#define __CS_MEDCOUPLING_MESH_HXX__

/*============================================================================
 * Usage of MEDCoupling base components.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_mesh.h"

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

  char               *sel_criteria;   /* Element selection criteria
                                         (if provided) */

  int                 elt_dim;        /* Element dimension */

  cs_lnum_t           n_elts;         /* Number of coupled elements */
  cs_lnum_t          *elt_list;       /* List of associated elements
                                         (0 to n-1) */

  cs_lnum_t           n_vtx;          /* Number of vertices */
  cs_lnum_t          *vtx_list;       /* List of associated vertices
                                         (0 to n-1), or NULL */

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
 * \brief   create a new cs_medcoupling_mesh_t instance based on cs_mesh_t
 *
 * \param[in] csmesh              pointer to cs_mesh_t instance
 * \param[in] name                name of the mesh
 * \param[in] selection_criteria  selection criteria string
 * \param[in] elt_dim             dimension of elements (2: faces, 3: cells)
 * \param[in] use_bbox            Use a reduced bounding box
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_from_base(cs_mesh_t   *csmesh,
                              const char  *name,
                              const char  *selection_criteria,
                              int          elt_dim,
                              int          use_bbox);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   create a new cs_medcoupling_mesh_t instance based on cs_mesh_t
 *
 * \param[in] csmesh      pointer to cs_mesh_t instance
 * \param[in] name        name of the mesh
 * \param[in] n_elts      local number of elements
 * \param[in] elt_ids     list of local elements
 * \param[in] elt_dim     dimension of elements (2: faces, 3: cells)
 * \param[in] use_bbox    use a reduced bounding box
 *
 * \return  pointer to the newly created cs_medcoupling_mesh_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_mesh_t *
cs_medcoupling_mesh_from_ids(cs_mesh_t       *csmesh,
                             const char      *name,
                             cs_lnum_t        n_elts,
                             const cs_lnum_t  elt_ids[],
                             int              elt_dim,
                             int              use_bbox);

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
 * \brief Destroy all cs_medcoupling_mesh_t instances.
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_mesh_destroy_all(void);

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
 * \brief Return a cs_medcoupling_mesh_t structure's number of vertices
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return associated number of vertices
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_medcoupling_mesh_get_n_vertices(cs_medcoupling_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a cs_medcoupling_mesh_t structure's (parent) vertices list
 *
 * \param[in] mesh  cs_medcoupling_mesh_t pointer
 *
 * \return ids of associated vertices, or NULL if all or no local vertices
 *         o parent mesh are present.
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_medcoupling_mesh_get_vertex_list(cs_medcoupling_mesh_t  *m);

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
/*!
 * \brief Returns a pointer to a MEDCouplingUMesh of a plane.
 *
 * \param[in] origin   Plane origin coordinates
 * \param[in] normal   Plane normal vector
 * \param[in] length1  Plane's edge length along first axis
 * \param[in] length2  Plane's edge length along second axis
 *
 * \return pointer to the MEDCouplingUMesh structure.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
MEDCouplingUMesh *
#else
void *
#endif
cs_medcoupling_create_plane_mesh(const cs_real_t origin[],
                                 const cs_real_t normal[],
                                 const cs_real_t length1,
                                 const cs_real_t length2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a pointer to a MEDCouplingUMesh of a disc.
 *
 * \param[in] origin     Disc origin coordinates
 * \param[in] normal     Disc normal vector
 * \param[in] radius     Disc radius
 * \param[in] n_sectors  Number of sectors for discretization. If negative,
 *                       default value of 36 is taken.
 *
 * \return pointer to the MEDCouplingUMesh structure.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
MEDCouplingUMesh *
#else
void *
#endif
cs_medcoupling_create_disc_mesh(const cs_real_t origin[],
                                const cs_real_t normal[],
                                const cs_real_t radius,
                                const int       n_sectors);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a pointer to a MEDCouplingUMesh of an annulus
 *
 * \param[in] origin     Annulus origin coordinates
 * \param[in] normal     Annulus normal vector
 * \param[in] radius1    Annulus inner radius
 * \param[in] radius2    Annulus outer radius
 * \param[in] n_sectors  Number of sectors for discretization. If negative,
 *                       default value of 36 is taken.
 *
 * \return pointer to the MEDCouplingUMesh structure.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)
MEDCouplingUMesh *
#else
void *
#endif
cs_medcoupling_create_annulus_mesh(const cs_real_t origin[],
                                   const cs_real_t normal[],
                                   const cs_real_t radius1,
                                   const cs_real_t radius2,
                                   const int       n_sectors);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MEDCOUPLING_MESH_HXX__ */
