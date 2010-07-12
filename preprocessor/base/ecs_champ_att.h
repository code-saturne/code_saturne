#ifndef _ECS_CHAMP_ATT_H_
#define _ECS_CHAMP_ATT_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_champ_t' décrivant un champ
 *   et propres aux champs auxiliaires de type "attribut"
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab_glob.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"
#include "ecs_descr.h"

/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_champ.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui assemble un champ donné dans un champ récepteur donné
 *
 *  L'assemblage consiste à :
 *  - regrouper sous la même position les valeurs des 2 champs
 *    (cela suppose donc que les 2 champs ont le même nombre de positions)
 *  - assembler les membres des descripteurs des 2 champs
 *    Les descripteurs des 2 champs peuvent être à `NULL'
 *    et si le descripteur du champ récepteur est à `NULL',
 *          le descripteur du champ assemblé est celui du champ à assembler
 *
 *  Le champ à assembler est détruit apres assemblage
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
ecs_champ_att__assemble(ecs_champ_t  *champ_recept,
                        ecs_champ_t  *champ_assemb);


/*----------------------------------------------------------------------------
 *  Fonction réalisant la mise à jour des valeurs d'un attribut des élément
 *   lorque ceux-ci sont renumérotés. Chaque élément d'origine a au
 *   maximum un élément correspondant, mais un nouvel élément peut résulter
 *   de la fusion de plusieurs éléments d'origine, et en conserver tous
 *   les attributs.
 *----------------------------------------------------------------------------*/

void
ecs_champ_att__herite(ecs_champ_t    *champ_att_elt,
                      ecs_tab_int_t  *tab_old_new);


/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un champ
 *   en fusionnant les propriétés de ses éléments
 *   qui sont identiquement transformes par le vecteur de transformation donné
 *----------------------------------------------------------------------------*/

void
ecs_champ_att__fusionne(ecs_champ_t         *this_champ_att,
                        size_t               nbr_elt_new,
                        const ecs_tab_int_t  vect_transf);


/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la liste des numéros de famille des éléments
 *
 *  Pour les éléments de famille 0 ou n'ayant pas de famille, on leur
 *   attribue le numéro de famille par défaut
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_champ_att__fam_elt(size_t          n_elts,
                       int            *elt_fam,
                       ecs_tab_int_t  *tab_nbr_elt_fam);


/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un tableau donnant pour chaque valeur du champ
 *   le nombre d'éléments ayant cette valeur
 *----------------------------------------------------------------------------*/

ecs_int_t *
ecs_champ_att__ret_nbr_elt_fam(size_t        n_elts,
                               int          *elt_fam,
                               size_t        nbr_val_fam);

/*----------------------------------------------------------------------------
 *  Fonction qui construit les familles à partir
 *   de la liste chaînée de tous les champs de type "attribut"
 *   pour toutes des entités ; ces familles sont ajoutées à la liste
 *   chaînée fournie en argument ;
 *
 *  Elle remplace le champ attribut par un champ "famille" par entité
 *
 *  Elle détermine aussi :
 *   - le nombre de familles
 *----------------------------------------------------------------------------*/

int *
ecs_champ_att__construit_fam(ecs_champ_t     **champ_att,
                             ecs_famille_t   **vect_fam_tete,
                             int               num_fam_deb,
                             int              *nbr_fam);

/*----------------------------------------------------------------------------
 *  Fonction qui crée le champ "groupe" à partir
 *   du champ "famille" et de la liste chaînée des familles
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ_att__cree_att_fam(size_t           n_elts,
                            int             *elt_fam,
                            ecs_famille_t   *famille);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CHAMP_ATT_H_ */
