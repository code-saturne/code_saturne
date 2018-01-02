/*============================================================================
 *  Définitions des fonctions
 *   associées à la structure `tab_t' décrivant un tableau
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

#include "ecs_def.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_mem.h"

/*----------------------------------------------------------------------------
 *  Headers for the current file
 *----------------------------------------------------------------------------*/

#include "ecs_tab.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local function defintions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *    Fonction de descente d'un arbre binaire pour le tri lexicographique
 *  d'un vecteur d'entiers (voir ecs_int_table_ord).
 *----------------------------------------------------------------------------*/

static void
ecs_loc_tab_int__desc_arbre(const ecs_int_t  *val_tab,
                            ecs_int_t         ltree,
                            ecs_int_t         ntree,
                            ecs_int_t        *elem_ord)
{
  ecs_int_t ktree;
  ecs_int_t i_save, p1, p2;

  i_save = *(elem_ord+ltree);

  while (ltree <= (ntree/2)) {

    ktree = (2*ltree)+1;

    if (ktree < ntree - 1) {

      p1 = *(elem_ord+ktree+1);
      p2 = *(elem_ord+ktree);

      if (*(val_tab+p1) > *(val_tab+p2)) ktree++;
    }

    if (ktree >= ntree) break;

    p1 = i_save;
    p2 = *(elem_ord+ktree);

    if (*(val_tab+p1) >= *(val_tab+p2)) break;

    *(elem_ord+ltree) = *(elem_ord+ktree);
    ltree = ktree;

  }

  *(elem_ord+ltree) = i_save;
}

/*----------------------------------------------------------------------------
 *    Fonction de descente d'un arbre binaire pour le tri lexicographique
 *  d'un vecteur de chaînes de caractères (voir ecs_char_table_ord).
 *----------------------------------------------------------------------------*/

static void
ecs_loc_tab_char__desc_arbre(char              **val_tab,
                             ecs_int_t           ltree,
                             const ecs_int_t     ntree,
                             ecs_int_t          *elem_ord)
{
  ecs_int_t ktree;
  ecs_int_t i_save, p1, p2;

  i_save = *(elem_ord+ltree);

  while (ltree <= (ntree/2)) {

    ktree = (2*ltree)+1;

    if (ktree < ntree - 1) {

      p1 = *(elem_ord+ktree+1);
      p2 = *(elem_ord+ktree);

      if (strcmp(*(val_tab+p1), *(val_tab+p2)) > 0) ktree++;
    }

    if (ktree >= ntree) break;

    p1 = i_save;
    p2 = *(elem_ord+ktree);

    if (strcmp(*(val_tab+p1), *(val_tab+p2)) >= 0) break;

    *(elem_ord+ltree) = *(elem_ord+ktree);
    ltree = ktree;

  }

  *(elem_ord+ltree) = i_save;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui crée un tableau de dimension donnée
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__cree(size_t  dim_tab)
{
  ecs_tab_int_t tab;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tab.nbr = dim_tab;
  ECS_MALLOC(tab.val, tab.nbr, ecs_int_t);

  return tab;
}

/*----------------------------------------------------------------------------
 *  Fonction qui crée un tableau de dimension donnée
 *   et initialise les valeurs avec une constante
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__cree_init(size_t     dim_tab,
                       ecs_int_t  val_init)
{
  size_t      itab;
  ecs_tab_int_t tab;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tab = ecs_tab_int__cree(dim_tab);

  for (itab = 0; itab < tab.nbr; itab++)
    tab.val[itab] = val_init;

  return tab;
}

/*----------------------------------------------------------------------------
 *  Fonction qui transforme le tableau donné en son inverse
 *----------------------------------------------------------------------------
 *
 * Exemple :
 * -------
 *                          .---.---.---.---.---.---.
 *     Tableau en entrée :  | 2 | 5 | 4 | 3 | 0 | 1 |
 *                          `---'---'---'---'---'---'
 *                            0   1   2   3   4   5
 *
 *                          .---.---.---.---.---.---.
 *     Tableau en sortie :  | 4 | 5 | 0 | 3 | 2 | 1 |
 *                          `---'---'---'---'---'---'
 *                            0   1   2   3   4   5
 *
 *----------------------------------------------------------------------------*/

void
ecs_tab_int__inverse(ecs_tab_int_t  *this_tab)
{
  size_t      itab   ;
  ecs_int_t * val_tab;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_tab->nbr != 0);
  assert(this_tab->val != NULL);

  ECS_MALLOC(val_tab, this_tab->nbr, ecs_int_t);

  for (itab = 0; itab < this_tab->nbr; itab++)
    val_tab[itab] = this_tab->val[itab];

  for (itab = 0; itab < this_tab->nbr; itab++)
    this_tab->val[val_tab[itab]] = itab;

  ECS_FREE(val_tab);
}

/*----------------------------------------------------------------------------
 *  Fonction de tri lexicographique d'un vecteur d'entiers.
 *
 *  La liste n'est pas modifiée directement,
 *   mais on construit un vecteur de renumérotation,
 *   afin de pouvoir appliquer cette renumérotation à d'autres tableaux
 *
 *  Le tri utilise est de type "heapsort", de complexité O(nlog(n)).
 *  Les éléments sont rangés en ordre croissant.
 *----------------------------------------------------------------------------*/

void
ecs_tab_int__trie(const ecs_tab_int_t  this_vect,
                  ecs_tab_int_t        vect_renum)
{
  ecs_int_t    i, i_save;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Création de l'arbre binaire vec_renum->val_tab[vect_renum->nbr] */

  for (i = (this_vect.nbr / 2) - 1; i >= 0; i--) {

    ecs_loc_tab_int__desc_arbre(this_vect.val,
                                i,
                                this_vect.nbr,
                                vect_renum.val);
  }

  /* Tri de l'arbre binaire */

  for (i = this_vect.nbr - 1; i > 0; i--) {

    i_save                = *(vect_renum.val    );
    *(vect_renum.val    ) = *(vect_renum.val + i);
    *(vect_renum.val + i) = i_save               ;


    ecs_loc_tab_int__desc_arbre(this_vect.val,
                                0,
                                i,
                                vect_renum.val);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction de tri lexicographique d'un vecteur de chaînes de caractères
 *
 *  La liste n'est pas modifiée directement,
 *   mais on construit un vecteur de renumérotation,
 *   afin de pouvoir appliquer cette renumérotation à d'autres tableaux
 *
 *  Le tri utilise est de type "heapsort", de complexité O(nlog(n)).
 *  Les éléments sont rangés en ordre croissant.
 *----------------------------------------------------------------------------*/

void
ecs_tab_char__trie(const ecs_tab_char_t  this_vect,
                   ecs_tab_int_t         vect_renum)
{
  ecs_int_t    i, i_save;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Création de l'arbre binaire vec_renum->val_tab[vect_renum->nbr] */

  for (i = (this_vect.nbr / 2) - 1; i >= 0; i--) {

    ecs_loc_tab_char__desc_arbre(this_vect.val,
                                 i,
                                 this_vect.nbr,
                                 vect_renum.val);
  }

  /* Tri de l'arbre binaire */

  for (i = this_vect.nbr - 1; i > 0; i--) {

    i_save                = *(vect_renum.val    );
    *(vect_renum.val    ) = *(vect_renum.val + i);
    *(vect_renum.val + i) = i_save               ;

    ecs_loc_tab_char__desc_arbre(this_vect.val,
                                 0,
                                 i,
                                 vect_renum.val);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui trie un vecteur d'entiers donné
 *   en renvoyant le vecteur trié
 *
 *  La fonction détermine aussi le vecteur de renumérotation des indices
 *  (pour des indices commençant à `0')
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__trie_et_renvoie(const ecs_tab_int_t  this_vect,
                             ecs_tab_int_t        vect_renum)
{
  size_t      ival;
  ecs_tab_int_t vect_trie;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ival = 0; ival < this_vect.nbr; ival++)
    vect_renum.val[ival] = ival;

  ecs_tab_int__trie(this_vect ,
                    vect_renum );

  ECS_MALLOC(vect_trie.val, this_vect.nbr, ecs_int_t);
  vect_trie.nbr = this_vect.nbr;

  for (ival = 0; ival < this_vect.nbr; ival++)
    vect_trie.val[ival] = this_vect.val[vect_renum.val[ival]];

  return vect_trie;
}

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
                              ecs_tab_int_t         vect_renum)
{
  size_t      ival;
  ecs_tab_char_t vect_trie;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ival = 0; ival < this_vect.nbr; ival++)
    vect_renum.val[ival] = ival;

  ecs_tab_char__trie(this_vect ,
                     vect_renum );

  ECS_MALLOC(vect_trie.val, this_vect.nbr, char *);
  vect_trie.nbr = this_vect.nbr;

  for (ival = 0; ival < this_vect.nbr; ival++)
    vect_trie.val[ival] = this_vect.val[vect_renum.val[ival]];

  return vect_trie;
}

/*----------------------------------------------------------------------------
 *  Fonction qui compacte un vecteur de chaînes de caractères donné
 *   en renvoyant le vecteur compacté; les chaînes ne sont pas dupliquées,
 *   seuls les pointeurs sont copiés.
 *
 *  Le vecteur d'origine doit être trié.
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_tab_char__compacte(const ecs_tab_char_t  this_vect)
{
  size_t         ival;
  size_t         ival_last;

  ecs_tab_char_t vect_compact ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ECS_MALLOC(vect_compact.val, this_vect.nbr, char *);
  vect_compact.nbr = this_vect.nbr;

  if (this_vect.nbr == 0)
    return vect_compact;

  ival_last = 0;
  vect_compact.val[0] = this_vect.val[0];

  for (ival = 1; ival < this_vect.nbr; ival++) {

    if (strcmp(this_vect.val[ival], vect_compact.val[ival_last]) != 0) {

      vect_compact.val[++ival_last] = this_vect.val[ival];

    }
  }

  vect_compact.nbr = ival_last + 1;
  ECS_REALLOC(vect_compact.val,
              vect_compact.nbr, char *);

  return vect_compact;
}

/*----------------------------------------------------------------------------
 *  Fonction de recherche d'une collection d'entiers
 *   dans une autre collection d'entiers strictement ordonnée.
 *   (Méthode de recherche dichotomique)
 *
 *  La fonction détermine un vecteur d'indices correspondant
 *   à la position des entiers dans le vecteur ou est faite la recherche
 *  Si un entier n'est pas contenu dans le vecteur,
 *   on lui adresse un "indice" `-1'
 *----------------------------------------------------------------------------*/

void
ecs_tab_int__recherche(ecs_tab_int_t  this_vect_rec,
                       ecs_tab_int_t  vect_ord,
                       ecs_tab_int_t  vect_ind)
{
  ecs_int_t   val_rec;   /* Valeur de l'entier  à rechercher                 */
  size_t      irec;      /* Indice de boucle sur les entiers à rechercher    */

  ecs_int_t   mid;       /* Valeur au centre dans l'intervalle de recherche  */
  ecs_int_t   left;      /* Valeur à gauche dans l'intervalle de recherche   */
  ecs_int_t   right;     /* Valeur à droite dans la recherche dichotomique   */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (irec = 0; irec < this_vect_rec.nbr; irec++) {

    val_rec = *(this_vect_rec.val + irec);

    if (val_rec <= vect_ord.val[0])
      *(vect_ind.val + irec) = 0;

    else if (val_rec > vect_ord.val[vect_ord.nbr - 1])
      *(vect_ind.val + irec) = vect_ord.nbr;

    else {

      left  = 0;
      right = vect_ord.nbr - 1;

      while ((right - left) > 1) {

        mid = (right + left) / 2;  /* Division entière */
        if (val_rec <= vect_ord.val[mid])
          right = mid;
        else
          left  = mid;

      }

      *(vect_ind.val + irec) = right;

    }

    if (   *(vect_ind.val + irec) >= (ecs_int_t)(vect_ord.nbr)
        || val_rec != vect_ord.val[*(vect_ind.val + irec)] ) {

      /* Échec de la recherche :                                 */
      /*  l'élément `val_rec' n'a pas été trouvé dans le vecteur */
      /* => on affecte une valeur `-1'                           */

      *(vect_ind.val + irec) = -1;

    }

  } /* Fin de la boucle sur les éléments de la liste à rechercher */
}

/*----------------------------------------------------------------------------
 *  Fonction de construction du tableau de remplacement référence -> indice
 *   Si bool_copie est à true, on alloue et on renvoie une copie de
 *   tab_att_reference, qui n'est pas modifie; sinon, tab_att_reference est
 *   transforme.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_tab_int__ref_en_indice(ecs_tab_int_t        tab_att_reference,
                           const ecs_tab_int_t  tab_val_idx,
                           bool                 bool_copie)
{
  size_t         ival;

  ecs_tab_int_t  tab_ind_ref_ord;
  ecs_tab_int_t  tab_renum;
  ecs_tab_int_t  tab_trie;

  ecs_tab_int_t  tab_ind_reference;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocations et initialisations */

  tab_ind_reference.nbr = tab_att_reference.nbr;

  if (bool_copie == true)
    ECS_MALLOC(tab_ind_reference.val, tab_ind_reference.nbr, ecs_int_t);

  else
    tab_ind_reference.val = tab_att_reference.val;

  tab_ind_ref_ord.nbr = tab_ind_reference.nbr;
  ECS_MALLOC(tab_ind_ref_ord.val, tab_ind_ref_ord.nbr, ecs_int_t);

  tab_renum.nbr = tab_val_idx.nbr;
  ECS_MALLOC(tab_renum.val, tab_renum.nbr,  ecs_int_t);

  /* renumérotation de vec_idx_sous_ent->val_tab  */

  tab_trie = ecs_tab_int__trie_et_renvoie(tab_val_idx,
                                          tab_renum);

  ecs_tab_int__recherche(tab_att_reference,
                         tab_trie,
                         tab_ind_ref_ord);

  ECS_FREE(tab_trie.val);

  for (ival = 0; ival < tab_ind_reference.nbr; ival++)
    tab_ind_reference.val[ival] = tab_renum.val[tab_ind_ref_ord.val[ival]];

  /* Libération de mémoire */

  ECS_FREE(tab_renum.val);
  ECS_FREE(tab_ind_ref_ord.val);

  /* Retour */

  return tab_ind_reference;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
