#ifndef _ECS_MAILLAGE_H_
#define _ECS_MAILLAGE_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associees a la structure `ecs_maillage_t' decrivant un maillage
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
 *                                 Visibilite
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
 *  Fichiers `include' publics  du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                         Déclaration de la structure
 *============================================================================*/

typedef struct _ecs_maillage_t ecs_maillage_t;


typedef enum {

  ECS_ENTMAIL_NONE = -1,
  ECS_ENTMAIL_FAC,
  ECS_ENTMAIL_CEL,
  ECS_N_ENTMAIL

} ecs_entmail_t;

/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define a new empty mesh structure with nodal connectivity.
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_maillage__cree_nodal(void);

/*----------------------------------------------------------------------------
 * Free a mesh structure.
 *----------------------------------------------------------------------------*/

void
ecs_maillage__detruit(ecs_maillage_t  **this_maillage);

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_maillage_t' donnee
 *   dans le fichier preprocessor_dump.txt
 *----------------------------------------------------------------------------*/

void
ecs_maillage__imprime(const ecs_maillage_t  *maillage,
                      ecs_int_t              nbr_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui retourne le type d'éntité de plus grande dimension
 *   contenue dans une structure `ecs_maillage_t'
 *----------------------------------------------------------------------------*/

ecs_entmail_t
ecs_maillage__ret_entmail_max(const ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_maillage_t'
 *----------------------------------------------------------------------------*/

float
ecs_maillage__ret_taille(const ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Suppression des sommets ne participant pas à la connectivité
 *   et fusion des éléments surfaciques confondus éventuels
 *----------------------------------------------------------------------------*/

void
ecs_maillage__nettoie_nodal(ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation des éléments en
 *   connectivité nodale.
 *
 *  La liste de cellules avec erreur est optionnelle.
 *----------------------------------------------------------------------------*/

void
ecs_maillage__orient_nodal(ecs_maillage_t    *maillage,
                           ecs_tab_int_t     *liste_cel_err,
                           bool               correc_orient);

/*----------------------------------------------------------------------------
 *  Fonction qui assigne la tete de la liste chainee des familles donnee
 *   a la structure de maillage donnee
 *----------------------------------------------------------------------------*/

void
ecs_maillage__definit_famille(ecs_maillage_t   *maillage,
                              ecs_famille_t    *vect_famille[2]);

/*----------------------------------------------------------------------------
 *  Fonction realisant, a partir d'une connectivite de maillage donnee,
 *   la connectivite descendante du maillage
 *----------------------------------------------------------------------------*/

void
ecs_maillage__connect_descend(ecs_maillage_t * maillage);

/*----------------------------------------------------------------------------
 *  Fonction realisant le tri des elements suivant leur type geometrique
 *  La fonction affiche le nombre d'elements par type geometrique
 *----------------------------------------------------------------------------*/

void
ecs_maillage__trie_typ_geo(ecs_maillage_t  *maillage);


/*----------------------------------------------------------------------------
 *  Fonction qui définit un nouveau maillage
 *   par extraction d'une partie du maillage donné
 *
 *  Les éléments à extraire doivent être tous de même dimension :
 *  cellules ou faces ou arêtes ou sommets
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_maillage__extrait(ecs_maillage_t       *maillage,
                      ecs_entmail_t         entmail_sel,
                      const ecs_tab_int_t  *liste_filtre);

/*----------------------------------------------------------------------------
 *  Fonction qui concatène dans un maillage récepteur donné,
 *   un maillage à concaténer donné.
 *
 *  Le maillage à concaténer est détruit.
 *----------------------------------------------------------------------------*/

void
ecs_maillage__concatene_nodal(ecs_maillage_t  *maillage_recept,
                              ecs_maillage_t  *maillage_concat);

/*----------------------------------------------------------------------------
 *  Fonction qui construit la liste des cellules attachées à une liste
 *  de faces fournie en argument.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_maillage__liste_cel_fac(ecs_maillage_t       *maillage,
                            const ecs_tab_int_t   liste_fac);

/*----------------------------------------------------------------------------
 *  Fonction qui calcule les coordonnées min et max du domaine
 *----------------------------------------------------------------------------*/

void
ecs_maillage__calc_coo_ext(ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Fonction qui construit les familles
 *----------------------------------------------------------------------------*/

void
ecs_maillage__cree_famille(ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Fonction qui detruit les familles
 *----------------------------------------------------------------------------*/

void
ecs_maillage__detruit_famille(ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Fonction qui construit les attributs "groupe" a partir des familles
 *----------------------------------------------------------------------------*/

void
ecs_maillage__cree_attributs(ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Fonction qui supprime les attributs "groupe"
 *----------------------------------------------------------------------------*/

void
ecs_maillage__supprime_attributs(ecs_maillage_t  *maillage);

/*----------------------------------------------------------------------------
 *  Vérification d'un maillage et calcul de critères de qualité
 *----------------------------------------------------------------------------*/

bool
ecs_maillage__verif(ecs_maillage_t  *maillage,
                    ecs_post_t      *cas_post);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_MAILLAGE_H_ */
