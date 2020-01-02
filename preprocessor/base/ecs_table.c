/*============================================================================
 *  Définitions des fonctions de base
 *   associées à la structure `ecs_table_t' décrivant une table
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <assert.h>
#include <stdio.h>
#include <string.h> /* strcpy() */


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr_chaine.h"
#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction d'impression d'une table avec position réglée en ASCII
 *----------------------------------------------------------------------------*/

static void
_imprime_pos_pas(FILE             *fic_imp,
                 size_t            nbr_ent,
                 size_t            pos_pas,
                 const ecs_int_t   val_int[],
                 size_t            nbr_imp)
{
  /* Variables locales */

  size_t  ient;
  size_t  iloc;

  size_t  ind_ent_1 = 0;
  size_t  ind_ent_2;

  char pos_str[32]; /* Largement surdimensionné pour contenir une
                       chaîne de type [%10d], un entier "long" pouvant
                       nécessiter un format plus large */

  /* Instructions */

  assert(val_int != NULL);

  ind_ent_2 = ECS_MIN(nbr_ent, nbr_imp);


  /* Impression des valeurs */
  /*========================*/

  while (1) {

    for (ient = ind_ent_1; ient < ind_ent_2; ient++) {

      if (pos_pas == 1) {

        fprintf(fic_imp, "%50s %12lu %12ld" "\n",
                " ", (unsigned long)(ient+1),
                (long)val_int[ient]);

      }
      else if (pos_pas > 1) {

        sprintf(pos_str, "[%d]", (int)(pos_pas*ient + 1));

        fprintf(fic_imp, "%37s %12lu %12s %12ld" "\n",
                " ", (unsigned long)(ient+1), pos_str,
                (long)val_int[pos_pas*ient]);

        for (iloc = 1; iloc < pos_pas; iloc++)
          fprintf(fic_imp, "%63s %12ld" "\n",
                  " ", (long)val_int[pos_pas*ient + iloc]);

      }
    }

    if (ind_ent_2 == nbr_ent)
      break;

    ind_ent_1 = ECS_MAX(nbr_ent - nbr_imp, nbr_imp);

    if (ind_ent_1 > ind_ent_2)

      fprintf(fic_imp,
              "%77s", "............\n");

    ind_ent_2 = nbr_ent;

  }
}


/*----------------------------------------------------------------------------
 *  Fonction d'impression d'une table avec position non réglée en ASCII
 *  (Table entier uniquement)
 *----------------------------------------------------------------------------*/

static void
_imprime_pos_tab(FILE             *fic_imp,
                 size_t            nbr_ent,
                 ecs_size_t        pos_tab[],
                 const ecs_int_t   val_tab[],
                 size_t            nbr_imp)
{
  /* Variables locales */

  size_t  ient;
  size_t  iloc;

  size_t  ind_ent_1 = 0;
  size_t  ind_ent_2;

  size_t  nbr_loc;


  /* Instructions */

  assert(pos_tab != NULL);

  ind_ent_2 = ECS_MIN(nbr_ent, nbr_imp);



  /* Impression des valeurs */
  /*========================*/

  while (1) {

    for (ient = ind_ent_1; ient < ind_ent_2; ient++) {

      nbr_loc = pos_tab[ient + 1] - pos_tab[ient];

      if (nbr_loc > 0)
        fprintf(fic_imp, "%37s %12lu %12lu %12ld\n",
                " ", (unsigned long)(ient+1),
                (unsigned long)pos_tab[ient],
                (long)val_tab[pos_tab[ient] - 1]);
      else
        fprintf(fic_imp, "%37s %12lu %12lu\n",
                " ", (unsigned long)(ient+1),
                (unsigned long)pos_tab[ient]);

      for (iloc = 1; iloc < nbr_loc; iloc++)
        fprintf(fic_imp, "%63s %12ld\n",
                " ", (long)val_tab[pos_tab[ient] + iloc - 1]);

    }

    if (ind_ent_2 == nbr_ent)
      break;

    ind_ent_1 = ECS_MAX(nbr_ent - nbr_imp, nbr_imp);

    if (ind_ent_1 > ind_ent_2)

      fprintf(fic_imp,
              "%77s", "............\n");

    ind_ent_2 = nbr_ent;

  }

  fprintf(fic_imp, "%50s %12lu\n",
          " ", (unsigned long)pos_tab[nbr_ent]);
}

/*----------------------------------------------------------------------------
 *  Fonction qui concatène dans une table réceptrice donnée,
 *   une table à concaténer donnée
 *
 *  La concaténation de 2 tables consiste à concaténer :
 *  - les tables des positions des 2 tables;
 *  - les tables des valeurs   des 2 tables;
 *  - les listes chaînées des descripteurs des 2 tables
 *    (nécessaire uniquement pour des tables de type "attribut")
 *
 *  Les autres membres de la table réceptrice ne sont pas modifiés
 *----------------------------------------------------------------------------*/

static void
_table__concatene(ecs_table_t  *table_recept,
                  ecs_table_t  *table_concat)
{
  size_t    ipos;
  size_t    ival;
  size_t    nbr_elt_recept;
  size_t    nbr_val_recept;
  size_t    nbr_val_concat;
  size_t    pos_recept_fin;

  ecs_descr_t  *descr_concat_copie;
  ecs_tab_int_t tab_renum_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_recept != NULL);
  assert(table_concat != NULL);

  /* Dimensions avant concaténation */

  nbr_elt_recept = table_recept->nbr;
  nbr_val_recept = ecs_table__ret_val_nbr(table_recept);
  nbr_val_concat = ecs_table__ret_val_nbr(table_concat);

  /* Pour les attributs,                                 */
  /*  il faut renuméroter les descripteurs               */
  /*  et propager ces nouvelles valeurs sur les éléments */

  tab_renum_descr.nbr = 0;
  tab_renum_descr.val = NULL;

  if (table_concat->descr != NULL) {

    /* Unification des descripteurs d'attributs */

    descr_concat_copie = ecs_descr_chaine__copie(table_concat->descr);

    tab_renum_descr = ecs_descr_chaine__concatene(&table_recept->descr,
                                                  &descr_concat_copie);
  }

  /* Membres restant à modifier pour la table receptrice : */
  /* - `nbr_elt'                                           */
  /* - `pos*'                                              */
  /* - `val*'                                              */

  table_recept->nbr += table_concat->nbr;

  /* Traitement de la table des positions; si l'on n'a pas
     une REGLE identique de part et d'autre, on doit la reconstruire */

  /* Si l'on a un pas identique de part et d'autre, on n'a rien à faire */

  if (   table_recept->pos != NULL
      || table_concat->pos != NULL
      || table_recept->pas != table_concat->pas) {

    /* 1ère étape : construire ou agrandir le tableau des positions,
       et le remplir des valeurs de la table initiale */

    if (table_recept->pos == NULL) {

      ECS_MALLOC(table_recept->pos, table_recept->nbr + 1, ecs_size_t);
      table_recept->pos[0] = 1;
      for (ipos = 0; ipos <= nbr_elt_recept; ipos++)
        table_recept->pos[ipos] = (table_recept->pas * ipos) + 1;

    }
    else
      ECS_REALLOC(table_recept->pos, table_recept->nbr + 1, ecs_size_t);

    /* 2ème étape : ajouter les positions à concaténer */

    pos_recept_fin = table_recept->pos[nbr_elt_recept];

    if (table_concat->pos == NULL) {

      for (ipos = 1; ipos <= table_concat->nbr; ipos++)
        table_recept->pos[nbr_elt_recept + ipos]
          = pos_recept_fin + (ipos * table_concat->pas);

    }
    else { /* if (table_concat->pos != NULL) */

      for (ipos = 1; ipos <= table_concat->nbr; ipos++)
        table_recept->pos[nbr_elt_recept + ipos]
          = pos_recept_fin - 1 + table_concat->pos[ipos];

    }

  }

  if (table_recept->pos != NULL)
    table_recept->pas = 0;

  /* Traitement de la table des valeurs */
  /*------------------------------------*/

  /* On concatène les tables de valeurs,
     en renumérotant éventuellement des attributs */

  if (table_recept->nbr > 0) {

    /* 1ère étape : construire ou agrandir le tableau des valeurs,
       et le remplir des valeurs de la table initiale */

    if (nbr_elt_recept == 0) {
      assert(table_recept->val == NULL);
      ECS_MALLOC(table_recept->val,
                 nbr_val_recept + nbr_val_concat,
                 ecs_int_t);
    }
    else {
      assert(table_recept->val != NULL);
      ECS_REALLOC(table_recept->val,
                  nbr_val_recept + nbr_val_concat,
                  ecs_int_t);
    }

    /* 2ème étape : ajouter les valeurs à concaténer, en les renumérotant,
       si nécessaire */

    if (tab_renum_descr.nbr == 0 && nbr_val_concat > 0)

      memcpy(((ecs_int_t *)table_recept->val) + nbr_val_recept,
             table_concat->val,
             nbr_val_concat * sizeof(ecs_int_t));

    else {

      ecs_int_t *tab_val_recept = table_recept->val;
      ecs_int_t *tab_val_concat = table_concat->val;

      for (ival = 0; ival < nbr_val_concat; ival++)
        tab_val_recept[nbr_val_recept + ival]
          = tab_renum_descr.val[ECS_ABS(tab_val_concat[ival]) - 1] + 1;

    }

  } /* Fin de la concaténation des tables de valeurs */


  /* Suppression du tableau de renumérotation des descripteurs */

  if (tab_renum_descr.nbr > 0)
    ECS_FREE(tab_renum_descr.val);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_table_t'
 *
 *  La structure devient propriétaire des tableaux tab_pos et tab_val
 *   fournis en argument.
 *
 *   nbr      : Nombre d'éléments à remplir
 *   pas      : Pas des positions  si REGLE
 *   pos      : Positions de la table si non REGLE
 *   val      : Valeurs de la table
 *   descr    : Pointeur sur le descripteur
 *   statut_e : Statut dans une transformation
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__cree(size_t                nbr,
                size_t                pas,
                ecs_size_t           *pos,
                ecs_int_t            *val,
                ecs_descr_t          *descr)
{
  ecs_table_t  *this_table;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocation de la structure `ecs_table_t' */
  /*------------------------------------------*/

  ECS_MALLOC(this_table, 1, ecs_table_t);

  this_table->nbr = nbr;

  /* Définition des positions des valeurs (itérateur) */
  /*--------------------------------------------------*/

  this_table->pas = pas;
  this_table->pos = pos;

  ecs_table__pos_en_regle(this_table);

  /* Définition de la table des valeurs (conteneur) */
  /*------------------------------------------------*/

  this_table->val = val;

  /* Affectation du descripteur de table */
  /*-------------------------------------*/

  this_table->descr = descr;

  return this_table;
}

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_table_t'
 *
 *   nbr      : Nombre d'éléments à remplir
 *   nbr_val  : Nombre de valeurs à remplir
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__alloue(size_t  nbr,
                  size_t  nbr_val)
{
  ecs_table_t  *this_table;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocation de la structure `ecs_table_t' */
  /*------------------------------------------*/

  ECS_MALLOC(this_table, 1, ecs_table_t);

  this_table->nbr = nbr;

  /* Définition des positions des valeurs (itérateur) */
  /*--------------------------------------------------*/

  this_table->pas = 0;

  ECS_MALLOC(this_table->pos, nbr + 1, ecs_size_t);

  /* Définition de la table des valeurs (conteneur) */
  /*------------------------------------------------*/

  ECS_MALLOC(this_table->val, nbr_val, ecs_int_t);

  /* Affectation du descripteur de table */
  /*-------------------------------------*/

  this_table->descr = NULL;

  return this_table;
}

/*----------------------------------------------------------------------------
 *  Fonction libérant une structure `ecs_table_t' donnée en argument.
 *----------------------------------------------------------------------------*/

void
ecs_table__detruit(ecs_table_t  **this_table)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table == NULL)
    return;

  if (*this_table == NULL)
    return;

  /* Libération du contenu de la structure `ecs_table_t' */
  /*=================================================*/

  /* Libération de la structure des positions */
  /*------------------------------------------*/

  if ((*this_table)->pos != NULL)
    ECS_FREE((*this_table)->pos);

  /* Libération de la structure des valeurs */
  /*----------------------------------------*/

  if ((*this_table)->val != NULL)
    ECS_FREE((*this_table)->val);

  /* Libération du descripteur de table */
  /*------------------------------------*/

  /* Appel à la fonction de libération d'un descripteur de table */

  if ((*this_table)->descr != NULL)
    ecs_descr_chaine__detruit(&((*this_table)->descr));

  /* Libération de la structure `ecs_table_t' */
  /*==========================================*/

  ECS_FREE(*this_table);
}

/*----------------------------------------------------------------------------
 *  Fonction qui convertit, si possible,
 *   le tableau des positions d'une table en REGLE
 *----------------------------------------------------------------------------*/

void
ecs_table__pos_en_regle(ecs_table_t  *this_table)
{
  size_t      ipos;
  size_t      pos_pas;
  bool        bool_regle;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);

  bool_regle = true;

  if (this_table->pos != NULL && this_table->nbr > 0) {

    this_table->pas = 0;

    pos_pas = this_table->pos[1] - this_table->pos[0];

    for (ipos = 1; ipos < this_table->nbr; ipos++) {

      if (this_table->pos[ipos + 1] - this_table->pos[ipos] != pos_pas) {
        bool_regle = false;
        break;
      }

    }

    if (bool_regle == true) {

      this_table->pas = pos_pas;

      ECS_FREE(this_table->pos);

    }
  }

}

/*----------------------------------------------------------------------------
 *  Fonction qui construit, si nécessaire, un tableau des positions à
 *   partir d'une REGLE.
 *----------------------------------------------------------------------------*/

void
ecs_table__regle_en_pos(ecs_table_t  *this_table)
{
  size_t ipos;
  ecs_size_t *tab_pos;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);

  if (this_table->pos == NULL) {

    ECS_MALLOC(tab_pos, this_table->nbr + 1, ecs_size_t);

    for (ipos = 0; ipos <= this_table->nbr; ipos++)
      tab_pos[ipos] = (ipos * this_table->pas) + 1;

    this_table->pos = tab_pos;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui libère, si possible, le tableau des positions d'une table.
 *  Ce tableau ne doit pas avoir été modifié.
 *----------------------------------------------------------------------------*/

void
ecs_table__libere_pos(ecs_table_t  *this_table)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);

  if (this_table->pas != 0 && this_table->pos != NULL)
    ECS_FREE(this_table->pos);
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_table_t' donnée
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_table__imprime(const ecs_table_t  *this_table,
                   size_t              imp_col,
                   size_t              nbr_imp,
                   FILE               *fic_imp)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table   != NULL);
  assert(fic_imp      != NULL);

  imp_col++;

  /* Impression des tables d'information de la table */
  /*-------------------------------------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, "nbr_elt", ECS_TYPE_size_t,
                       &this_table->nbr);

  ecs_fic__imprime_val(fic_imp, imp_col, "pos_pas", ECS_TYPE_size_t,
                       &this_table->pas);

  ecs_fic__imprime_ptr(fic_imp, imp_col, "pos_tab", this_table->pos);
  ecs_fic__imprime_ptr(fic_imp, imp_col, "val_tab", this_table->val);


  /* Impression des positions et des valeurs */
  /*-----------------------------------------*/

  if (this_table->pos == NULL && this_table->val != NULL)
    _imprime_pos_pas(fic_imp,
                     this_table->nbr,
                     this_table->pas,
                     this_table->val,
                     nbr_imp);

  else if (this_table->pos != NULL)
    _imprime_pos_tab(fic_imp,
                     this_table->nbr,
                     this_table->pos,
                     this_table->val,
                     nbr_imp);

  /* Impression de la liste chaînée des descripteurs */
  /*-------------------------------------------------*/

  /* Impression du pointeur sur le descripteur de tête */

  ecs_fic__imprime_val(fic_imp, imp_col, "descr_tete", ECS_TYPE_void,
                       this_table->descr);

  if (this_table->descr != NULL) {

    /* Appel à la fonction d'impression d'une chaîne de descripteurs de table */

    ecs_descr_chaine__imprime(this_table->descr,
                              imp_col + 1,
                              fic_imp);
  }

}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_table_t'
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_taille(const ecs_table_t  *this_table)
{
  size_t        nbr_val;
  size_t        taille;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table == NULL)
    return 0;

  taille = sizeof(*this_table);

  if (this_table->pos != NULL)
    taille += (sizeof(ecs_int_t) * (this_table->nbr + 1));

  if (this_table->val != NULL) {
    nbr_val = ecs_table__ret_val_nbr(this_table);
    taille += (sizeof(ecs_int_t) * nbr_val);
  }

  if (this_table->descr != NULL)
    taille += ecs_descr_chaine__ret_taille(this_table->descr);

  return taille;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie une table entièrement réallouée
 *   dont le contenu est copié à partir de la table donnée
 *
 *  Le membre donnant le lien sur une table suivant `l_table_sui'
 *   n'est pas copié et est mis à `NULL'
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__copie(ecs_table_t  *table_init)
{
  size_t        nbr_val;
  ecs_table_t * this_table;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_init != NULL);

  ECS_MALLOC(this_table, 1, ecs_table_t);

  this_table->nbr = table_init->nbr;

  this_table->pas = table_init->pas;

  if (table_init->pos != NULL) {
    ECS_MALLOC(this_table->pos, table_init->nbr + 1, ecs_size_t);
    memcpy(this_table->pos,
           table_init->pos,
           (table_init->nbr + 1) * sizeof(ecs_size_t));
  }
  else
    this_table->pos = NULL;

  if (table_init->val != NULL) {
    nbr_val = ecs_table__ret_val_nbr(table_init);
    ECS_MALLOC(this_table->val, nbr_val, ecs_int_t);
    memcpy(this_table->val,
           table_init->val,
           sizeof(ecs_int_t) * nbr_val);
  }
  else
    this_table->val = NULL;

  this_table->descr  = ecs_descr_chaine__copie(table_init->descr);

  return this_table;
}

/*----------------------------------------------------------------------------
 *  Fonction qui créé une structure `ecs_table_t'
 *   à partir d'un tableau `tab_elt' contenant les valeurs de la table.
 *
 *  Si un élément n'a pas de valeur associée, la valeur correspondante
 *   dans `tab_elt' est `0'
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__transforme_tableau(size_t                nbr_elt,
                              const ecs_int_t      *tab_elt,
                              ecs_descr_t          *descr)
{
  ecs_size_t *pos_tab;
  ecs_int_t *val_tab;
  size_t cpt_val;
  size_t ielt;

  ecs_table_t * this_table;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(nbr_elt != 0);
  assert(tab_elt != NULL);

  /* Allocation de la structure `ecs_table_t' */
  /*--------------------------------------*/

  ECS_MALLOC(this_table, 1, ecs_table_t);

  this_table->nbr = nbr_elt;

  /* Construction des tableaux de positions et de valeurs */
  /*------------------------------------------------------*/

  ECS_MALLOC(pos_tab, nbr_elt + 1, ecs_size_t);
  ECS_MALLOC(val_tab, nbr_elt    , ecs_int_t);

  pos_tab[0] = 1;
  cpt_val    = 0;

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    if (tab_elt[ielt] != 0) {
      pos_tab[ielt + 1]  = pos_tab[ielt] + 1;
      val_tab[cpt_val++] = tab_elt[ielt];
    }
    else
      pos_tab[ielt + 1]  = pos_tab[ielt];

  }

  ECS_REALLOC(val_tab, cpt_val, ecs_int_t);

  /* Création de la table des positions des valeurs (itérateur) */
  /*------------------------------------------------------------*/

  this_table->pas = 0;
  this_table->pos = pos_tab;

  ecs_table__pos_en_regle(this_table);

  /* Création de la table des valeurs (conteneur) */
  /*----------------------------------------------*/

  this_table->val = val_tab;

  /* Affectation du descripteur de table */
  /*-------------------------------------*/

  this_table->descr = descr;

  return this_table;
}

/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre d'éléments associés à une table donnée
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_elt_nbr(const ecs_table_t  *this_table)
{
  size_t retval = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table != NULL)
    retval = this_table->nbr;

  return retval;
}


/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre de valeurs associées à une table donnée
 *----------------------------------------------------------------------------*/

size_t
 ecs_table__ret_val_nbr(const ecs_table_t  *this_table)
{
  size_t  nbr_val;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table   != NULL);

  if (this_table->pos != NULL)
    nbr_val = this_table->pos[this_table->nbr] - 1;
  else
    nbr_val = (this_table->pas * this_table->nbr);

  return nbr_val;
}

/*----------------------------------------------------------------------------
 *  Fonction retournant le nombre de descripteurs d'une table donnée
 *----------------------------------------------------------------------------*/

size_t
ecs_table__ret_descr_nbr(const ecs_table_t  *this_table)
{
  assert(this_table != NULL);

  return ecs_descr_chaine__ret_nbr(this_table->descr);
}

/*----------------------------------------------------------------------------
 *  Fonction libérant un pointeur sur le tableau des positions d'une
 *   structure `ecs_table_t' donnée.
 *
 *  Si les positions correspondent à une REGLE, le tableau est libéré.
 *   Sinon, il est conservé par la structure ecs_table_t.
 *----------------------------------------------------------------------------*/

void
ecs_table__libere_pos_tab(const ecs_table_t  *this_table,
                          ecs_size_t         *pos_tab)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);
  assert(this_table->pos == NULL || this_table->pos == pos_tab);

  if (this_table->pos == NULL)
    ECS_FREE(pos_tab);
}

/*----------------------------------------------------------------------------
 *  Fonction qui concatène deux tables, et supprime la table à concaténer
 *----------------------------------------------------------------------------*/

void
ecs_table__concatene(ecs_table_t  **this_table,
                     ecs_table_t  **concat_table,
                     size_t         nbr_elt_init,
                     size_t         nbr_elt_concat)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);
  assert(concat_table != NULL);

  /* ------------------------------------------------------------------ */
  /* PARTIE REALISEE PENDANT LA FUSION                                  */
  /* ------------------------------------------------------------------ */
  /*    Si la table n'existe pas dans l'entité de maillage à concaténer */
  /* et si c'est une table de type "attribut"                           */
  /* (pour lequel tous les éléments doivent être renseignes)            */
  /* on initialise les valeurs de l'attribut pour la table à concaténer */
  /* ------------------------------------------------------------------ */

  /* La table existe dans l'entité de maillage à concaténer */

  if (*concat_table != NULL) {

    /* On le concatène */

    if (*this_table != NULL) {
      _table__concatene(*this_table, *concat_table);
      ecs_table__detruit(concat_table);
    }
    else {
      ecs_table__prolonge(*concat_table,
                          nbr_elt_init,
                          0);
      *this_table = *concat_table;
      *concat_table = NULL;
    }
  }

  /* La table n'existe pas dans l'entité de maillage à concaténer */

  else {
    ecs_table__prolonge(*this_table,
                        0,
                        nbr_elt_concat);

  } /* Fin : boucle sur les tables de l'entité de maillage receptrice */
}

/*----------------------------------------------------------------------------
 *  Fonction qui prolonge une table réceptrice donné
 *
 *  Il s'agit en fait de concaténer la table avec une table vide. Seule la
 *  table des positions est modifiée. Les autres membres de la structure de
 *  la table réceptrice ne sont pas modifiés.
 *----------------------------------------------------------------------------*/

void
ecs_table__prolonge(ecs_table_t  *this_table,
                    size_t        nbr_elt_prec,
                    size_t        nbr_elt_suiv)
{
  size_t      ipos;
  size_t      nbr_elt_ini;
  ecs_int_t   pos_fin;
  bool        bool_pos_ini;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table == NULL)
    return;

  ecs_table__regle_en_pos(this_table);

  nbr_elt_ini = this_table->nbr;

#if defined(SX) && defined(_SX) /* NEC SX compiler may bug on : operator */
  if (this_table->pos != NULL)
    bool_pos_ini = true;
  else
    bool_pos_ini = false;
#else
  bool_pos_ini = (this_table->pos != NULL) ? true : false;
#endif

  /* Mise à jour du nombre d'éléments */

  this_table->nbr += nbr_elt_prec + nbr_elt_suiv;

  /* Si la table n'est pas déjà vide, la table des positions ne
     correspondra pas à une REGLE, et devra être construite */

  if (this_table->pos != NULL || this_table->pas != 0) {

    ECS_REALLOC(this_table->pos, this_table->nbr + 1, ecs_size_t);

    if (bool_pos_ini == true) {
      memmove(this_table->pos + (nbr_elt_prec * sizeof(ecs_int_t)),
              this_table->pos,
              (this_table->nbr + 1) * sizeof(ecs_int_t));
    }
    else {
      for (ipos = 0; ipos <= nbr_elt_ini; ipos++)
        this_table->pos[nbr_elt_prec + ipos]
          = (ipos * this_table->pas) + 1;
    }

    this_table->pos[0] = 1;

    for (ipos = 0; ipos < nbr_elt_prec; ipos++)
      this_table->pos[ipos + 1] = 1;

    pos_fin = this_table->pos[nbr_elt_prec + nbr_elt_ini];
    for (ipos = 0; ipos < nbr_elt_suiv; ipos++)
      this_table->pos[nbr_elt_prec + nbr_elt_ini + ipos + 1] = pos_fin;

  }

  ecs_table__pos_en_regle(this_table);
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un table
 *   en appliquant directement le vecteur de transformation donné
 *   sur ses positions
 *
 *  Le nombre de valeurs transformées doit être égal
 *  au nombre de valeurs avant transformation
 *----------------------------------------------------------------------------*/

void
ecs_table__transforme_pos(ecs_table_t          *this_table,
                          size_t                nbr_elt_ref,
                          const ecs_tab_int_t   vect_transf)
{
  ecs_table_t  *table_ref;
  size_t       nbr_val_ref;
  ecs_int_t    ival;
  size_t       ielt_ref;
  size_t       ielt_transf;
  size_t       ival_ref;
  ecs_int_t    elt_nbr_val;
  ecs_int_t    pre_pos;
  ecs_int_t    cpt_val;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table == NULL)
    return;

  ecs_table__regle_en_pos(this_table);

  assert(vect_transf.nbr == nbr_elt_ref);
  assert(this_table->nbr == nbr_elt_ref);

  nbr_val_ref = ecs_table__ret_val_nbr(this_table);

  table_ref = ecs_table__alloue(nbr_elt_ref, nbr_val_ref);

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref + 1; ielt_ref++)
    table_ref->pos[ielt_ref] = this_table->pos[ielt_ref];

  for (ival_ref = 0; ival_ref < nbr_val_ref; ival_ref++)
    table_ref->val[ival_ref] = this_table->val[ival_ref];

  this_table->pos[0] = 1;

  cpt_val = 0;

  for (ielt_transf = 0; ielt_transf < nbr_elt_ref; ielt_transf++) {

    pre_pos = table_ref->pos[vect_transf.val[ielt_transf]];

    elt_nbr_val
      = table_ref->pos[vect_transf.val[ielt_transf] + 1] - pre_pos;

    this_table->pos[ielt_transf + 1]
      = this_table->pos[ielt_transf] + elt_nbr_val;

    for (ival = 0; ival < elt_nbr_val; ival++) {

      this_table->val[cpt_val++]
        = table_ref->val[pre_pos - 1 + ival];

    }
  }

  ecs_table__detruit(&table_ref);

  ecs_table__pos_en_regle(this_table);
}

/*----------------------------------------------------------------------------
 *  Fonction qui incrémente les valeurs d'une table donnée
 *   d'une constante donnée
 *----------------------------------------------------------------------------*/

void
ecs_table__incremente_val(ecs_table_t      *this_table,
                          const ecs_int_t   increment)
{
  size_t     nbr_val;
  size_t     ival;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table == NULL)
    return;

  nbr_val = ecs_table__ret_val_nbr(this_table);

  for (ival = 0; ival < nbr_val; ival++) {

    if (this_table->val[ival] > 0)
      this_table->val[ival] += increment;
    else if (this_table->val[ival] < 0)
      this_table->val[ival] -= increment;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un vecteur indexé
 *   en appliquant directement le vecteur de transformation donné
 *   sur les valeurs associées à ses éléments
 *----------------------------------------------------------------------------*/

void
ecs_table__renumerote(ecs_table_t          *this_table,
                      const ecs_tab_int_t   vect_transf,
                      const ecs_tab_int_t   signe_elt)
{
  ecs_table_t *table_ref;
  size_t       nbr_val_ref;
  ecs_int_t    val_ref;
  ecs_int_t    sgn_ref;
  size_t       ival_ref;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table == NULL)
    return;

  nbr_val_ref = ecs_table__ret_val_nbr(this_table);

  table_ref = ecs_table__alloue(0, nbr_val_ref);

  for (ival_ref = 0; ival_ref < nbr_val_ref; ival_ref++)
    table_ref->val[ival_ref] = this_table->val[ival_ref];

  for (ival_ref = 0; ival_ref < nbr_val_ref; ival_ref++) {

    val_ref = ECS_ABS(table_ref->val[ival_ref]);
    sgn_ref = table_ref->val[ival_ref] / val_ref;

    this_table->val[ival_ref]  = vect_transf.val[val_ref - 1] + 1;

    this_table->val[ival_ref] *= signe_elt.val[val_ref - 1] * sgn_ref;

  }

  ecs_table__detruit(&table_ref);
}

/*----------------------------------------------------------------------------
 *  Fonction qui détermine une nouvelle table à partir d'une table de
 *   référence en extrayant de ce dernier les éléments sélectionnés
 *   par le tableau de booléens
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table__extrait(ecs_table_t  *table_ref,
                   bool          elt_select[])
{
  size_t        cpt_elt_new;
  size_t        cpt_val_new;
  size_t        nbr_elt_ref;
  size_t        nbr_val_ref;
  size_t        pos_ref_inf;
  size_t        pos_ref_sup;

  size_t        ielt_ref;
  size_t        ipos_ref;

  ecs_table_t  *table_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_table__regle_en_pos(table_ref);

  nbr_elt_ref = table_ref->nbr;
  nbr_val_ref = ecs_table__ret_val_nbr(table_ref);

  table_new = ecs_table__alloue(nbr_elt_ref,
                                nbr_val_ref);

  /* Extraction de la table */

  nbr_elt_ref = table_ref->nbr;
  nbr_val_ref = ecs_table__ret_val_nbr(table_ref);

  table_new->pos[0] = 1;
  cpt_elt_new           = 0;
  cpt_val_new           = 0;

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref; ielt_ref++) {

    if (elt_select[ielt_ref] == true) {

      /* L'élément est à extraire */

      pos_ref_inf = table_ref->pos[ielt_ref    ] - 1;
      pos_ref_sup = table_ref->pos[ielt_ref + 1] - 1;

      for (ipos_ref = pos_ref_inf; ipos_ref < pos_ref_sup; ipos_ref++)
        table_new->val[cpt_val_new++] = table_ref->val[ipos_ref];

      table_new->pos[cpt_elt_new + 1] = cpt_val_new + 1;

      cpt_elt_new++;

    }
  }

  ECS_REALLOC(table_new->pos, cpt_elt_new + 1, ecs_size_t);
  ECS_REALLOC(table_new->val, cpt_val_new,     ecs_int_t);

  table_new->nbr = cpt_elt_new;

  ecs_table__pos_en_regle(table_ref);

  return table_new;
}

/*----------------------------------------------------------------------------*/
