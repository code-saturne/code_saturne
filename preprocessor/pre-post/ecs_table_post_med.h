#ifndef _ECS_TABLE_POST_MED_H_
#define _ECS_TABLE_POST_MED_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_table_t' décrivant une table
 *   et réalisant les sorties au format MED
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


/*============================================================================
 *                                 Visibilité
 *============================================================================*/

#include "cs_config.h"

#if defined(HAVE_MED)

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_med.h"


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
 *  Fonction écrivant les familles
 *----------------------------------------------------------------------------*/

void
ecs_table_post_med__ecr_famille(const char           *nom_maillage,
                                const ecs_famille_t  *famille_elt,
                                const ecs_famille_t  *famille_inf,
                                ecs_med_t            *cas_med);

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu des tables asociees aux sommets
 *----------------------------------------------------------------------------*/

void
ecs_table_post_med__ecr_som(const char         *nom_maillage,
                            size_t              n_vertices,
                            ecs_coord_t         vertex_coords[],
                            const ecs_med_t    *cas_med);

/*----------------------------------------------------------------------------
 *  Fonction qui écrit les connectivités des éléments
 *   selon leur type géometrique
 *
 *  Les éléments doivent avoir ete triés suivant leur type géometrique
 *----------------------------------------------------------------------------*/

void
ecs_table_post_med__ecr_elt(const char           *nom_maillage,
                            ecs_table_t          *table_def,
                            const int             elt_fam[],
                            const ecs_tab_int_t  *tab_elt_typ_geo,
                            const ecs_med_t      *cas_med);

#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_POST_MED_H_ */
