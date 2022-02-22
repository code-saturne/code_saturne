#ifndef _ECS_TABLE_POST_H_
#define _ECS_TABLE_POST_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_table_t' décrivant une table
 *   et réalisant les sorties pour post-traitement
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

/*============================================================================
 *                                 Visibilité
 *============================================================================*/

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_post_ens.h"

#if defined(HAVE_MED)
#include "ecs_med.h"
#endif


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction ecrivant les elements d'une table donne pour le post traitement
 *
 *  Les elements doivent avoir ete tries suivant leur type geometrique
 *----------------------------------------------------------------------------*/

void
ecs_table_post__ecr_elt(const char            *nom_maillage,
                        int                    dim_entite_max,
                        size_t                 n_vertices,
                        ecs_coord_t            vertex_coords[],
                        ecs_table_t           *table_def,
                        const int              elt_fam[],
                        ecs_table_t           *table_def_inf,
                        const int              elt_fam_inf[],
                        const ecs_famille_t   *famille_elt,
                        const ecs_famille_t   *famille_inf,
                        ecs_post_type_t        type_post,
                        ecs_post_t            *cas_post);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_POST_H_ */
