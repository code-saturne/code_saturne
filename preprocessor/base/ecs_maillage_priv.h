#ifndef _ECS_MAILLAGE_PRIV_H_
#define _ECS_MAILLAGE_PRIV_H_

/*============================================================================
 *  Definition privee de la structure `_ecs_maillage_t' decrivant un maillage
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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


/*============================================================================
 *                                 Visibilite
 *============================================================================*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"
#include "ecs_famille.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*============================================================================
 *                       Definition de macros
 *============================================================================*/


/*============================================================================
 *                       Definition de la structure
 *============================================================================*/

typedef enum {

  ECS_MAILLAGE_CONNECT_NODALE,
  ECS_MAILLAGE_CONNECT_DESCENDANTE

} ecs_maillage_connect_typ_t;

struct _ecs_maillage_t {

  ecs_maillage_connect_typ_t typ_connect;   /* Connectivity type */

  size_t            n_vertices;             /* Number of vertices */
  ecs_coord_t      *vertex_coords;          /* Vertex coordinates */

  ecs_table_t      *table_def[2];           /* Face and cell connectivity */
  ecs_table_t      *table_att[2];           /* Face and cell attributes */

  int              *elt_fam[2];             /* Face and cell families
                                               (size: n_faces/n_cells) */

  /* Additionnal connectivity (information on connected faces or
     cells, defined by connected couples, for non-conforming meshes) */

  size_t            n_connect_couples[2];
  ecs_int_t        *connect_couples[2];

  /* Linked family lists, per entity (for family descriptor information) */

  ecs_famille_t    *famille[2];
};

/*============================================================================
 * Les valeurs des tableaux `famille' sont les numeros des familles
 *  numerotees a partir de `1' comme suit :
 *  - en commencant par les familles des cellules
 *  - puis          par les familles des faces
 *============================================================================*/

#endif /* _ECS_MAILLAGE_PRIV_H_ */
