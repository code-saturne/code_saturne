#ifndef _ECS_CHAMP_POST_MED_H_
#define _ECS_CHAMP_POST_MED_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_champ_t' décrivant un champ
 *   et réalisant les sorties au format MED
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2009 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/


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

#include "ecs_tab_glob.h"


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

#include "ecs_champ.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction écrivant les familles
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_famille(const char           *nom_maillage,
                                const ecs_famille_t  *famille_elt,
                                const ecs_famille_t  *famille_inf,
                                ecs_med_t            *cas_med);

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu des champs asociees aux sommets
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_som(const char         *nom_maillage,
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
ecs_champ_post_med__ecr_elt(const char           *nom_maillage,
                            ecs_champ_t          *champ_def,
                            const int             elt_fam[],
                            const ecs_tab_int_t  *tab_elt_typ_geo,
                            const ecs_med_t      *cas_med);

/*----------------------------------------------------------------------------
 *  Fonction qui ajoute à une structure maillage_med les informations
 *   sur le nombre d'éléments de chaque type d'un maillage
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__cpt_elt_typ(const ecs_tab_int_t  *tab_elt_typ_geo,
                                const char           *nom_maillage,
                                ecs_med_t            *cas_med);

/*----------------------------------------------------------------------------
 *  Fonction qui écrit les valeurs par élément pour un tableau donné.
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_val(const ecs_tab_int_t  *tab_val,
                            const char           *nom_maillage,
                            const char           *nom_champ,
                            const ecs_med_t      *cas_med);

#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CHAMP_POST_MED_H_ */
