#ifndef _ECS_DESCR_CHAINE_H_
#define _ECS_DESCR_CHAINE_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à une liste chaînée de structures `ecs_descr_t' décrivant
 *   un descripteur de table
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

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C ou BFT
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_tab.h"

/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"

/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction libérant la portion d'une liste chaînée de descripteurs
 *   à partir d'un noeud dont le pointeur est donné en argument.
 *  Le noeud est à NULL au retour de la fonction
 *----------------------------------------------------------------------------*/

void
ecs_descr_chaine__detruit(ecs_descr_t  **descr_noeud);

/*----------------------------------------------------------------------------
 *  Fonction imprimant à partir d'un noeud `ecs_descr_t' donné
 *   une liste chaînée de tables
 *   sur le flux décrit par la structure `ecs_file_t'
 *----------------------------------------------------------------------------*/

void ecs_descr_chaine__imprime(const ecs_descr_t  *descr_noeud,
                               int                 imp_col,
                               FILE               *fic_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets
 *   d'une chaîné de structures `ecs_descr_t'
 *----------------------------------------------------------------------------*/

float
ecs_descr_chaine__ret_taille(const ecs_descr_t  *descr_noeud);

/*----------------------------------------------------------------------------
 *  Fonction qui ajoute à la fin d'une liste chaînée de descripteurs
 *   réceptrice dont la tête est donnée, une liste chaînée de descripteurs
 *   à concaténer dont la tête est donnée
 *
 *  Les numéros des descripteurs de la liste à concaténer sont incrementés
 *   à partir du nombre de descripteur de la liste réceptrice
 *
 *  Remarque: cette fonction se contente d'ajouter des descripteurs sans
 *            vérifier si le descripteur ajoute a le même contenu qu'un autre
 *            descripteur déjà présent dans la liste.
 *            Pour une vérification, utiliser `ecs_descr_chaine__concatene()'
 *----------------------------------------------------------------------------*/

void
ecs_descr_chaine__ajoute(ecs_descr_t  **descr_tete,
                         ecs_descr_t   *descr_concat_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre de descripteurs
 *   de la liste chaînée des descripteurs dont la tête est donnée
 *----------------------------------------------------------------------------*/

int
ecs_descr_chaine__ret_nbr(const ecs_descr_t  *descr_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui copie une liste chaînée de descripteurs
 *   dont la tête est donnée
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr_chaine__copie(ecs_descr_t  *descr_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui concatène,
 *   à la fin d'une liste chaînée de descripteurs dont la tête est donnée,
 *   une autre liste chaînée de descripteurs dont la tête est donnée,
 *   en supprimant les descripteurs déjà présents dans la 1ère liste
 *   et en décalant la renumérotation des descripteurs de la 2nde liste
 *
 *  La fonction renvoie la renumérotation des descripteurs de la 2nde liste
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_descr_chaine__concatene(ecs_descr_t  **descr_recept_tete,
                            ecs_descr_t  **descr_concat_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui affiche les contenus des descripteurs
 *   de la liste chaînée des descripteurs dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_descr_chaine__affiche(ecs_descr_t  *descr_tete,
                          int           decal);

/*----------------------------------------------------------------------------
 *  Fonction qui recherche dans une liste chaînée de descripteurs
 *   dont la tête est donnée,
 *   un numéro de descripteur donné
 *
 *  La fonction renvoie :
 *  -    le pointeur du descripteur si le numéro de descripteur a été trouve
 *  - ou NULL                       sinon
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr_chaine__cherche_num(ecs_descr_t  *descr_tete,
                              int           num);

/*----------------------------------------------------------------------------
 *  Fonction qui recherche dans une liste chaînée de descripteurs
 *   dont la tête est donnée,
 *   un descripteur ayant les mêmes type, identificateur et nom
 *   que le descripteur donné
 *
 *  La fonction renvoie :
 *  -    le numéro du descripteur si le descripteur   a     été trouve
 *  - ou ECS_DESCR_NUM_NUL        si le descripteur n'a pas été trouve
 *----------------------------------------------------------------------------*/

int
ecs_descr_chaine__trouve_num(ecs_descr_t        *descr_tete,
                             const ecs_descr_t  *descr_rech);

/*----------------------------------------------------------------------------
 *  Fonction qui crée une nouvelle chaîne de descripteurs
 *   à partir d'une chaîne de descripteurs dont la tête est donnée
 *  Un descripteur est copié dans la nouvelle chaîne si son numéro
 *   ne se transforme pas par le vecteur de transformation donné
 *   en `ECS_DESCR_NUM_NUL'
 *  Les membres du descripteur sont copies dans le nouveau sans modification
 *   sauf le numéro qui devient celui transformé par le vecteur
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr_chaine__renumerote(ecs_descr_t          *descr_tete ,
                             const ecs_tab_int_t   vect_transf);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre et la liste des pointeurs sur les noms
 *   des descripteurs de type groupe d'une liste chaînée dont la tête est
 *   donnée en argument
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_descr_chaine__ret_nom(ecs_descr_t   *descr_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la liste des références des descripteurs
 *   de la liste chaînée des descripteurs dont la tête est donnée
 *----------------------------------------------------------------------------*/

ecs_descr_t **
ecs_descr_chaine__ret_ref(ecs_descr_t  *descr_tete,
                          int          *nbr_descr);

/*----------------------------------------------------------------------------
 *  Fonction qui crée une nouvelle chaîne de descripteurs
 *   à partir d'une chaîne de descripteurs dont la tête est donnée
 *  Un descripteur est copié dans la nouvelle chaîne si son numéro
 *   ne se transforme pas par le vecteur de transformation donné
 *   en `ECS_DESCR_NUM_NUL'
 *  Les membres du descripteur sont copies dans le nouveau sans modification
 *   sauf le numéro qui devient celui transformé par le vecteur
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_descr_chaine__trie(ecs_descr_t  *descr_tete);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_DESCR_CHAINE_H_ */
