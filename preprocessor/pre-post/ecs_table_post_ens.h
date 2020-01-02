#ifndef _ECS_CHAMP_POST_ENS_H_
#define _ECS_CHAMP_POST_ENS_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_champ_t' décrivant un champ
 *   et réalisant les sorties pour post-traitement EnSight
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
 *                                 Visibilité
 *============================================================================*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C ou BFT
 *----------------------------------------------------------------------------*/

#include <ecs_file.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*============================================================================
 *                       Définition de macro
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui écrit les connectivités des éléments
 *   selon leur type géometrique
 *
 *  Les éléments doivent avoir été triés suivant leur type géometrique
 *----------------------------------------------------------------------------*/

void
ecs_table_post_ens__ecr_part(const char            *nom_maillage,
                             size_t                 n_vertices,
                             const ecs_coord_t      vertex_coords[],
                             ecs_table_t           *table_def,
                             const ecs_tab_int_t   *tab_elt_typ_geo,
                             ecs_post_ens_t        *cas_ens);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_POST_ENS_H_ */

