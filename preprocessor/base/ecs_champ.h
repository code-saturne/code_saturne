#ifndef _ECS_CHAMP_H_
#define _ECS_CHAMP_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à la structure `ecs_champ_t' décrivant un champ
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

typedef struct _ecs_champ_t  ecs_champ_t;


/*============================================================================
 *                         Définitions d'énumérations
 *============================================================================*/

typedef enum {

  ECS_CHAMP_NUL = -1,
  ECS_CHAMP_DEF,
  ECS_CHAMP_ATT,
  ECS_CHAMP_FAM,
  ECS_CHAMP_CNN,
  ECS_CHAMP_FIN

} ECS_CHAMP_E;


#define ECS_CHAMP_DEB ECS_CHAMP_DEF


/*============================================================================
 *                         Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_champ_t'
 *
 *  La structure devient propriétaire des tableaux tab_pos et tab_val
 *   fournis en argument.
 *
 *   nbr      : Nombre d'éléments à remplir
 *   pas      : Pas des positions  si REGLE
 *   pos      : Positions du champ si non REGLE
 *   val      : Valeurs du champ
 *   descr    : Pointeur sur le descripteur
 *   statut_e : Statut dans une transformation
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__cree(size_t         nbr,
                size_t         pas,
                ecs_size_t    *pos,
                ecs_int_t     *val,
                ecs_descr_t   *descr);

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_champ_t'
 *
 *   nbr      : Nombre d'éléments à remplir
 *   nbr_val  : Nombre de valeurs à remplir
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__alloue(size_t  nbr,
                  size_t  nbr_val);

/*----------------------------------------------------------------------------
 *  Fonction libérant une structure `ecs_champ_t' donnée en argument.
 *----------------------------------------------------------------------------*/

void
ecs_champ__detruit(ecs_champ_t  **this_champ);

/*----------------------------------------------------------------------------
 *  Fonction qui convertit, si possible,
 *   le tableau des positions d'un champ en REGLE
 *----------------------------------------------------------------------------*/

void
ecs_champ__pos_en_regle(ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction qui construit, si nécessaire, un tableau des positions à
 *   partir d'une REGLE.
 *----------------------------------------------------------------------------*/

void
ecs_champ__regle_en_pos(ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction qui libère, si possible, le tableau des positions d'un champ.
 *  Ce tableau ne doit pas avoir été modifié.
 *----------------------------------------------------------------------------*/

void
ecs_champ__libere_pos(ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_champ_t' donnée
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_champ__imprime(const ecs_champ_t  *this_champ,
                   size_t              imp_col,
                   size_t              nbr_imp,
                   FILE               *fic_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_champ_t'
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_taille(const ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un champ entièrement réalloué
 *   dont le contenu est copié à partir du champ donné
 *
 *  Le membre donnant le lien sur un champ suivant `l_champ_sui'
 *   n'est pas copié et est mis à `NULL'
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__copie(ecs_champ_t  *champ_init);

/*----------------------------------------------------------------------------
 *  Fonction qui créé une structure `ecs_champ_t'
 *   à partir d'un tableau `tab_elt' contenant les valeurs du champ.
 *
 *  Si un élément n'a pas de valeur associée, la valeur correspondante
 *   dans `tab_elt' est `0'
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__transforme_tableau(size_t                nbr_elt,
                              const ecs_int_t      *tab_elt,
                              ecs_descr_t          *descr);

/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre d'éléments associés à un champ donné
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_elt_nbr(const ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre de valeurs associées à un champ donné
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_val_nbr(const ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction retournant le nombre de descripteurs d'un champ donné
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_descr_nbr(const ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction retournant le type des valeurs d'un champ donné
 *----------------------------------------------------------------------------*/

ecs_type_t
ecs_champ__ret_val_typ(const ecs_champ_t  *this_champ);

/*----------------------------------------------------------------------------
 *  Fonction libérant un pointeur sur le tableau des positions d'une
 *   structure `ecs_champ_t' donnée.
 *
 *  Si les positions correspondent à une REGLE, le tableau est libéré.
 *   Sinon, il est conservé par la structure ecs_champ_t.
 *----------------------------------------------------------------------------*/

void
ecs_champ__libere_pos_tab(const ecs_champ_t  *this_champ,
                          ecs_size_t         *pos_tab);

/*----------------------------------------------------------------------------
 *  Fonction qui concatène deux champs, et supprime le champ à concaténer
 *----------------------------------------------------------------------------*/

void
ecs_champ__concatene(ecs_champ_t  **this_champ,
                     ecs_champ_t  **concat_champ,
                     size_t         nbr_elt_init,
                     size_t         nbr_elt_ent_concat);

/*----------------------------------------------------------------------------
 *  Fonction qui concatène deux champs de type connectivité,
 *   et supprime le champ à concaténer
 *----------------------------------------------------------------------------*/

void
ecs_champ__concatene_connect(ecs_champ_t  **this_champ,
                             ecs_champ_t  **concat_champ,
                             size_t         nbr_elt_init,
                             size_t         nbr_elt_concat);

/*----------------------------------------------------------------------------
 *  Fonction qui prolonge un champ récepteur donné
 *
 *  Il s'agit en fait de concaténer le champ avec un champ vide. Seule la
 *  table des positions est modifiée. Les autres membres de la structure du
 *  champ récepteur ne sont pas modifiés.
 *----------------------------------------------------------------------------*/

void
ecs_champ__prolonge(ecs_champ_t  *this_champ,
                    size_t        nbr_elt_prec,
                    size_t        nbr_elt_suiv);

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un champ
 *   en appliquant directement le vecteur de transformation donné
 *   sur ses positions
 *
 *  Le nombre de valeurs transformées doit être égal
 *  au nombre de valeurs avant transformation
 *----------------------------------------------------------------------------*/

void
ecs_champ__transforme_pos(ecs_champ_t          *this_champ,
                          size_t                nbr_elt_ref,
                          const ecs_tab_int_t   vect_transf);

/*----------------------------------------------------------------------------
 *  Fonction qui incrémente les valeurs d'un champ donné
 *   d'une constante donnée
 *----------------------------------------------------------------------------*/

void
ecs_champ__incremente_val(ecs_champ_t      *this_champ,
                          const ecs_int_t   increment);

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un vecteur indexé
 *   en appliquant directement le vecteur de transformation donné
 *   sur les valeurs associées à ses éléments
 *----------------------------------------------------------------------------*/

void
ecs_champ__renumerote(ecs_champ_t          *this_champ,
                      const ecs_tab_int_t   vect_transf,
                      const ecs_tab_int_t   signe_elt);

/*----------------------------------------------------------------------------
 *  Fonction qui détermine un nouveau champ à partir d'un champ de référence
 *   en extrayant de ce dernier les éléments sélectionnés
 *   par le tableau de booléens
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__extrait(ecs_champ_t            *champ_ref,
                   const ecs_tab_bool_t    bool_elt_select);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CHAMP_H_ */
