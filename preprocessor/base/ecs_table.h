#ifndef _ECS_CHAMP_H_
#define _ECS_CHAMP_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à la structure `ecs_table_t' décrivant une table
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

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                         Déclaration de la structure
 *============================================================================*/

typedef struct _ecs_table_t  ecs_table_t;


/*============================================================================
 *                         Définitions d'énumérations
 *============================================================================*/

typedef enum {

  ECS_TABLE_NUL = -1,
  ECS_TABLE_DEF,
  ECS_TABLE_ATT,
  ECS_TABLE_FAM,
  ECS_TABLE_CNN,
  ECS_TABLE_FIN

} ECS_TABLE_E;


#define ECS_TABLE_DEB ECS_TABLE_DEF


/*============================================================================
 *                         Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_table_t'
 *
 *  La structure devient propriétaire des tableaux tab_pos et tab_val
 *   fournis en argument.
 *
 *   nbr      : Nombre d'éléments à remplir
 *   pas      : Pas des positions  si REGLE
 *   pos      : Positions des éléments si non REGLE
 *   val      : Valeurs des éléments
 *   descr    : Pointeur sur le descripteur
 *   statut_e : Statut dans une transformation
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__cree(size_t         nbr,
                size_t         pas,
                ecs_size_t    *pos,
                ecs_int_t     *val,
                ecs_descr_t   *descr);

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_table_t'
 *
 *   nbr      : Nombre d'éléments à remplir
 *   nbr_val  : Nombre de valeurs à remplir
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__alloue(size_t  nbr,
                  size_t  nbr_val);

/*----------------------------------------------------------------------------
 *  Fonction libérant une structure `ecs_table_t' donnée en argument.
 *----------------------------------------------------------------------------*/

void
ecs_table__detruit(ecs_table_t  **this_table);

/*----------------------------------------------------------------------------
 *  Fonction qui convertit, si possible,
 *   le tableau des positions d'une table en REGLE
 *----------------------------------------------------------------------------*/

void
ecs_table__pos_en_regle(ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction qui construit, si nécessaire, un tableau des positions à
 *   partir d'une REGLE.
 *----------------------------------------------------------------------------*/

void
ecs_table__regle_en_pos(ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction qui libère, si possible, le tableau des positions d'un table.
 *  Ce tableau ne doit pas avoir été modifié.
 *----------------------------------------------------------------------------*/

void
ecs_table__libere_pos(ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_table_t' donnée
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_table__imprime(const ecs_table_t  *this_table,
                   size_t              imp_col,
                   size_t              nbr_imp,
                   FILE               *fic_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_table_t'
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_taille(const ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie une table entièrement réallouée
 *   dont le contenu est copié à partir de la table donnée
 *
 *  Le membre donnant le lien sur une table suivante `l_table_sui'
 *   n'est pas copié et est mis à `NULL'
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__copie(ecs_table_t  *table_init);

/*----------------------------------------------------------------------------
 *  Fonction qui créé une structure `ecs_table_t'
 *   à partir d'un tableau `tab_elt' contenant les valeurs du table.
 *
 *  Si un élément n'a pas de valeur associée, la valeur correspondante
 *   dans `tab_elt' est `0'
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__transforme_tableau(size_t                nbr_elt,
                              const ecs_int_t      *tab_elt,
                              ecs_descr_t          *descr);

/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre d'éléments associés à une table donnée
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_elt_nbr(const ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre de valeurs associées à une table donnée
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_val_nbr(const ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction retournant le nombre de descripteurs d'une table donnée
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_descr_nbr(const ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction retournant le type des valeurs d'une table donnée
 *----------------------------------------------------------------------------*/

ecs_type_t
ecs_table__ret_val_typ(const ecs_table_t  *this_table);

/*----------------------------------------------------------------------------
 *  Fonction libérant un pointeur sur le tableau des positions d'une
 *   structure `ecs_table_t' donnée.
 *
 *  Si les positions correspondent à une REGLE, le tableau est libéré.
 *   Sinon, il est conservé par la structure ecs_table_t.
 *----------------------------------------------------------------------------*/

void
ecs_table__libere_pos_tab(const ecs_table_t  *this_table,
                          ecs_size_t         *pos_tab);

/*----------------------------------------------------------------------------
 *  Fonction qui concatène deux tables, et supprime la table à concaténer
 *----------------------------------------------------------------------------*/

void
ecs_table__concatene(ecs_table_t  **this_table,
                     ecs_table_t  **concat_table,
                     size_t         nbr_elt_init,
                     size_t         nbr_elt_ent_concat);

/*----------------------------------------------------------------------------
 *  Fonction qui concatène deux tables de type connectivité,
 *   et supprime la table à concaténer
 *----------------------------------------------------------------------------*/

void
ecs_table__concatene_connect(ecs_table_t  **this_table,
                             ecs_table_t  **concat_table,
                             size_t         nbr_elt_init,
                             size_t         nbr_elt_concat);

/*----------------------------------------------------------------------------
 *  Fonction qui prolonge une table réceptrice donnée
 *
 *  Il s'agit en fait de concaténer le table avec une table vide. Seule la
 *  table des positions est modifiée. Les autres membres de la structure du
 *  table récepteur ne sont pas modifiés.
 *----------------------------------------------------------------------------*/

void
ecs_table__prolonge(ecs_table_t  *this_table,
                    size_t        nbr_elt_prec,
                    size_t        nbr_elt_suiv);

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'une table
 *   en appliquant directement le vecteur de transformation donné
 *   sur ses positions
 *
 *  Le nombre de valeurs transformées doit être égal
 *  au nombre de valeurs avant transformation
 *----------------------------------------------------------------------------*/

void
ecs_table__transforme_pos(ecs_table_t          *this_table,
                          size_t                nbr_elt_ref,
                          const ecs_tab_int_t   vect_transf);

/*----------------------------------------------------------------------------
 *  Fonction qui incrémente les valeurs d'une table donnée
 *   d'une constante donnée
 *----------------------------------------------------------------------------*/

void
ecs_table__incremente_val(ecs_table_t      *this_table,
                          const ecs_int_t   increment);

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un vecteur indexé
 *   en appliquant directement le vecteur de transformation donné
 *   sur les valeurs associées à ses éléments
 *----------------------------------------------------------------------------*/

void
ecs_table__renumerote(ecs_table_t          *this_table,
                      const ecs_tab_int_t   vect_transf,
                      const ecs_tab_int_t   signe_elt);

/*----------------------------------------------------------------------------
 *  Fonction qui détermine une nouvelle table à partir d'une table de référence
 *   en extrayant de ce dernier les éléments sélectionnés
 *   par le tableau de booléens
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__extrait(ecs_table_t  *table_ref,
                   bool          elt_select[]);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_H_ */
