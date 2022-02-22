#ifndef _ECS_CHAMP_ATT_H_
#define _ECS_CHAMP_ATT_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_table_t' décrivant une table
 *   et propres aux tables auxiliaires de type "attribut"
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"
#include "ecs_descr.h"

/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui assemble une table donnée dans une table réceptrice donnée
 *
 *  L'assemblage consiste à :
 *  - regrouper sous la même position les valeurs des 2 tables
 *    (cela suppose donc que les 2 tables ont le même nombre de positions)
 *  - assembler les membres des descripteurs des 2 tables
 *    Les descripteurs des 2 tables peuvent être à `NULL'
 *    et si le descripteur de la table réceptrice est à `NULL',
 *    le descripteur de la table assemblée est celui de la table à assembler
 *
 *  La table à assembler est détruite apres assemblage
 * ----------------------------------------------------------------------------
 *
 *  Exemple :
 *  =======
 *
 *  Soit la table à assembler :
 *
 *                         .---.---..---..---.---.---.
 *     assemb->val         | 5 | 3 || 4 || 5 | 2 | 6 |
 *                         `---'---'`---'`---'---'---'
 *                           0   1    2    3   4   5
 *
 *
 *                         .---.---.---.---.---.
 *     assemb->pos         | 1 | 3 | 4 | 4 | 7 |
 *                         `---'---'---'---'---'
 *                           0   1   2   3   4
 *
 *
 *  dans la table réceptrice :
 *
 *                         .---..---..---.---..---.
 *     recept->val         | 4 || 5 || 6 | 6 || 1 |
 *                         `---'`---'`---'---'`---'
 *                           0    1    2   3    4
 *
 *
 *                         .---.---.---.---.---.
 *     recept->pos         | 1 | 2 | 3 | 5 | 6 |
 *                         `---'---'---'---'---'
 *                           0   1   2   3   4
 *
 *
 *  La table réceptrice devient :
 *
 *                         .---.---.---..---.---..---.---..---.---.---.---.
 *     recept->val         | 4 | 5 | 3 || 5 | 4 || 6 | 6 || 1 | 5 | 2 | 6 |
 *                         `---'---'---'`---'---'`---'---'`---'---'---'---'
 *                           0   1   2    3   4    5   6    7   8   9   10
 *
 *
 *                         .---.---.---.---.---.
 *     recept->pos         | 1 | 4 | 6 | 8 | 12|
 *                         `---'---'---'---'---'
 *                           0   1   2   3   4
 *
 *----------------------------------------------------------------------------*/

void
ecs_table_att__assemble(ecs_table_t  *table_recept,
                        ecs_table_t  *table_assemb);


/*----------------------------------------------------------------------------
 *  Fonction réalisant la mise à jour des valeurs d'un attribut des élément
 *   lorque ceux-ci sont renumérotés. Chaque élément d'origine a au
 *   maximum un élément correspondant, mais un nouvel élément peut résulter
 *   de la fusion de plusieurs éléments d'origine, et en conserver tous
 *   les attributs.
 *----------------------------------------------------------------------------*/

void
ecs_table_att__herite(ecs_table_t    *table_att_elt,
                      ecs_tab_int_t  *tab_old_new);


/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'une table
 *   en fusionnant les propriétés de ses éléments
 *   qui sont identiquement transformes par le vecteur de transformation donné
 *----------------------------------------------------------------------------*/

void
ecs_table_att__fusionne(ecs_table_t         *this_table_att,
                        size_t               nbr_elt_new,
                        const ecs_tab_int_t  vect_transf);

/*----------------------------------------------------------------------------
 *  Fonction qui construit les familles à partir
 *   de la liste chaînée de tous les tables de type "attribut"
 *   pour toutes des entités ; ces familles sont ajoutées à la liste
 *   chaînée fournie en argument ;
 *
 *  Elle remplace la table attribut par une table "famille" par entité
 *
 *  Elle détermine aussi :
 *   - le nombre de familles
 *----------------------------------------------------------------------------*/

int *
ecs_table_att__construit_fam(ecs_table_t     **table_att,
                             ecs_famille_t   **vect_fam_tete,
                             int               num_fam_deb,
                             int              *nbr_fam);

/*----------------------------------------------------------------------------
 *  Fonction qui crée le table "groupe" à partir
 *   de la table "famille" et de la liste chaînée des familles
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table_att__cree_att_fam(size_t           n_elts,
                            int             *elt_fam,
                            ecs_famille_t   *famille);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_ATT_H_ */
