#ifndef _ECS_TAB_H_
#define _ECS_TAB_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `tab_t' décrivant un tableau
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "ecs_def.h"

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {
  size_t        nbr;
  ecs_int_t    *val;
} ecs_tab_int_t;

typedef struct {
  size_t       nbr;
  char       **val;
} ecs_tab_char_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui cree un tableau de dimension donnée
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__cree(size_t  dim_tab);

/*----------------------------------------------------------------------------
 *  Fonction qui crée un tableau de dimension donnée
 *   et initialise les valeurs avec une constante
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__cree_init(size_t      dim_tab,
                       ecs_int_t   val);

/*----------------------------------------------------------------------------
 *  Fonction qui transforme le tableau donne en son inverse
 *----------------------------------------------------------------------------*/

void
ecs_tab_int__inverse(ecs_tab_int_t  *this_tab);

/*----------------------------------------------------------------------------
 *  Fonction de tri lexicographique d'un vecteur d'entiers.
 *
 *  La liste n'est pas modifiée directement,
 *   mais on construit un vecteur de renumérotation,
 *   afin de pouvoir appliquer cette renumérotation à d'autres tableaux
 *
 *  Le tri utilisé est de type "heapsort", de complexité O(nlog(n)).
 *  Les éléments sont rangés en ordre croissant.
 *----------------------------------------------------------------------------*/

void
ecs_tab_int__trie(const ecs_tab_int_t  this_vect,
                  ecs_tab_int_t        vect_renum);

/*----------------------------------------------------------------------------
 *  Fonction de tri lexicographique d'un vecteur de chaînes de caractères
 *
 *  La liste n'est pas modifiée directement,
 *   mais on construit un vecteur de renumérotation,
 *   afin de pouvoir appliquer cette renumérotation à d'autres tableaux
 *
 *  Le tri utilisé est de type "heapsort", de complexité O(nlog(n)).
 *  Les éléments sont rangés en ordre croissant.
 *----------------------------------------------------------------------------*/

void
ecs_tab_char__trie(const ecs_tab_char_t  this_vect,
                   ecs_tab_int_t         vect_renum);

/*----------------------------------------------------------------------------
 *  Fonction qui trie un vecteur d'entiers donné
 *   en renvoyant le vecteur trié
 *
 *  La fonction détermine aussi le vecteur de renumérotation des indices
 *  (pour des indices commençant à `0')
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__trie_et_renvoie(const ecs_tab_int_t  this_vect,
                             ecs_tab_int_t        vect_renum);

/*----------------------------------------------------------------------------
 *  Fonction qui trie un vecteur de chaînes de caractères donné
 *   en renvoyant le vecteur trié. Les chaînes ne sont pas dupliquées,
 *   seuls les pointeurs sont copiés.
 *
 *  La fonction détermine aussi le vecteur de renumérotation des indices
 *  (pour des indices commençant à `0')
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_tab_char__trie_et_renvoie(const ecs_tab_char_t  this_vect,
                              ecs_tab_int_t         vect_renum);

/*----------------------------------------------------------------------------
 *  Fonction qui compacte un vecteur de chaînes de caractères donné
 *   en renvoyant le vecteur compacté; les chaînes ne sont pas dupliquées,
 *   seuls les pointeurs sont copiés.
 *
 *  Le vecteur d'origine doit être trié.
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_tab_char__compacte(const ecs_tab_char_t  this_vect);

/*----------------------------------------------------------------------------
 *  Fonction de recherche d'une collection d'entiers
 *   dans une autre collection d'entiers strictement ordonnée.
 *   (Méthode de recherche dichotomique)
 *
 *  La fonction retourne un vecteur  d'indices correspondant
 *   à la position des entiers dans le vecteur ou est faite la recherche
 *  Si un entier n'est pas contenu dans le vecteur,
 *   on lui adresse un "indice" `-1'
 *----------------------------------------------------------------------------*/

void
ecs_tab_int__recherche(ecs_tab_int_t  this_vect_rec,
                       ecs_tab_int_t  vect_ord,
                       ecs_tab_int_t  vect_ind);

/*----------------------------------------------------------------------------
 *  Fonction de construction du tableau de remplacement référence -> indice
 *   Si bool_copie est à true, on alloue et on renvoie une copie de
 *   tab_att_reference, qui n'est pas modifié; sinon, tab_att_reference est
 *   transformé.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__ref_en_indice(ecs_tab_int_t        tab_att_reference,
                           const ecs_tab_int_t  tab_val_idx,
                           bool                 bool_copie);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* _ECS_TAB_H_ */
