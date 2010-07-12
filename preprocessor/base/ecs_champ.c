/*============================================================================
 *  Définitions des fonctions de base
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

#include "ecs_champ.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_champ_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction d'impression d'un champ avec position réglée en ASCII
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
 *  Fonction d'impression d'un champ avec position non réglée en ASCII
 *  (Champ entier uniquement)
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
 *  Fonction qui concatène dans un champ récepteur donné,
 *   un champ à concaténer donné
 *
 *  La concaténation de 2 champs consiste à concaténer :
 *  - les tables des positions des 2 champs;
 *  - les tables des valeurs   des 2 champs;
 *  - les listes chaînées des descripteurs des 2 champs
 *    (nécessaire uniquement pour des champs de type "attribut")
 *
 *  Les autres membres de la structure du champ récepteur ne sont pas modifiés
 *----------------------------------------------------------------------------*/

static void
_champ__concatene(ecs_champ_t  *champ_recept,
                  ecs_champ_t  *champ_concat)
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

  assert(champ_recept != NULL);
  assert(champ_concat != NULL);

  /* Dimensions avant concaténation */

  nbr_elt_recept = champ_recept->nbr;
  nbr_val_recept = ecs_champ__ret_val_nbr(champ_recept);
  nbr_val_concat = ecs_champ__ret_val_nbr(champ_concat);

  /* Pour les attributs,                                 */
  /*  il faut renuméroter les descripteurs               */
  /*  et propager ces nouvelles valeurs sur les éléments */

  tab_renum_descr.nbr = 0;
  tab_renum_descr.val = NULL;

  if (champ_concat->descr != NULL) {

    /* Unification des descripteurs d'attributs */

    descr_concat_copie = ecs_descr_chaine__copie(champ_concat->descr);

    tab_renum_descr = ecs_descr_chaine__concatene(&champ_recept->descr,
                                                  &descr_concat_copie);
  }

  /* Membres restant à modifier pour le champ recepteur : */
  /* - `nbr_elt'                                          */
  /* - `pos*'                                             */
  /* - `val*'                                             */

  champ_recept->nbr += champ_concat->nbr;

  /* Traitement de la table des positions; si l'on n'a pas
     une REGLE identique de part et d'autre, on doit la reconstruire */

  /* Si l'on a un pas identique de part et d'autre, on n'a rien à faire */

  if (   champ_recept->pos != NULL
      || champ_concat->pos != NULL
      || champ_recept->pas != champ_concat->pas) {

    /* 1ère étape : construire ou agrandir le tableau des positions,
       et le remplir des valeurs du champ initial */

    if (champ_recept->pos == NULL) {

      ECS_MALLOC(champ_recept->pos, champ_recept->nbr + 1, ecs_size_t);
      champ_recept->pos[0] = 1;
      for (ipos = 0; ipos <= nbr_elt_recept; ipos++)
        champ_recept->pos[ipos] = (champ_recept->pas * ipos) + 1;

    }
    else
      ECS_REALLOC(champ_recept->pos, champ_recept->nbr + 1, ecs_size_t);

    /* 2ème étape : ajouter les positions à concaténer */

    pos_recept_fin = champ_recept->pos[nbr_elt_recept];

    if (champ_concat->pos == NULL) {

      for (ipos = 1; ipos <= champ_concat->nbr; ipos++)
        champ_recept->pos[nbr_elt_recept + ipos]
          = pos_recept_fin + (ipos * champ_concat->pas);

    }
    else { /* if (champ_concat->pos != NULL) */

      for (ipos = 1; ipos <= champ_concat->nbr; ipos++)
        champ_recept->pos[nbr_elt_recept + ipos]
          = pos_recept_fin - 1 + champ_concat->pos[ipos];

    }

  }

  if (champ_recept->pos != NULL)
    champ_recept->pas = 0;

  /* Traitement de la table des valeurs */
  /*------------------------------------*/

  /* On concatène les tables de valeurs,
     en renumérotant éventuellement des attributs */

  if (champ_recept->nbr > 0) {

    /* 1ère étape : construire ou agrandir le tableau des valeurs,
       et le remplir des valeurs du champ initial */

    if (nbr_elt_recept == 0) {
      assert(champ_recept->val == NULL);
      ECS_MALLOC(champ_recept->val,
                 nbr_val_recept + nbr_val_concat,
                 ecs_int_t);
    }
    else {
      assert(champ_recept->val != NULL);
      ECS_REALLOC(champ_recept->val,
                  nbr_val_recept + nbr_val_concat,
                  ecs_int_t);
    }

    /* 2ème étape : ajouter les valeurs à concaténer, en les renumérotant,
       si nécessaire */

    if (tab_renum_descr.nbr == 0 && nbr_val_concat > 0)

      memcpy(((ecs_int_t *)champ_recept->val) + nbr_val_recept,
             champ_concat->val,
             nbr_val_concat * sizeof(ecs_int_t));

    else {

      ecs_int_t *tab_val_recept = champ_recept->val;
      ecs_int_t *tab_val_concat = champ_concat->val;

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
ecs_champ__cree(size_t                nbr,
                size_t                pas,
                ecs_size_t           *pos,
                ecs_int_t            *val,
                ecs_descr_t          *descr)
{
  ecs_champ_t  *this_champ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocation de la structure `ecs_champ_t' */
  /*------------------------------------------*/

  ECS_MALLOC(this_champ, 1, ecs_champ_t);

  this_champ->nbr = nbr;

  /* Définition des positions des valeurs (itérateur) */
  /*--------------------------------------------------*/

  this_champ->pas = pas;
  this_champ->pos = pos;

  ecs_champ__pos_en_regle(this_champ);

  /* Définition de la table des valeurs (conteneur) */
  /*------------------------------------------------*/

  this_champ->val = val;

  /* Affectation du descripteur de champ */
  /*-------------------------------------*/

  this_champ->descr = descr;

  return this_champ;
}

/*----------------------------------------------------------------------------
 *  Fonction qui crée une structure `ecs_champ_t'
 *
 *   nbr      : Nombre d'éléments à remplir
 *   nbr_val  : Nombre de valeurs à remplir
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__alloue(size_t  nbr,
                  size_t  nbr_val)
{
  ecs_champ_t  *this_champ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocation de la structure `ecs_champ_t' */
  /*------------------------------------------*/

  ECS_MALLOC(this_champ, 1, ecs_champ_t);

  this_champ->nbr = nbr;

  /* Définition des positions des valeurs (itérateur) */
  /*--------------------------------------------------*/

  this_champ->pas = 0;

  ECS_MALLOC(this_champ->pos, nbr + 1, ecs_size_t);

  /* Définition de la table des valeurs (conteneur) */
  /*------------------------------------------------*/

  ECS_MALLOC(this_champ->val, nbr_val, ecs_int_t);

  /* Affectation du descripteur de champ */
  /*-------------------------------------*/

  this_champ->descr = NULL;

  return this_champ;
}

/*----------------------------------------------------------------------------
 *  Fonction libérant une structure `ecs_champ_t' donnée en argument.
 *----------------------------------------------------------------------------*/

void
ecs_champ__detruit(ecs_champ_t  **this_champ)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_champ == NULL)
    return;

  if (*this_champ == NULL)
    return;

  /* Libération du contenu de la structure `ecs_champ_t' */
  /*=================================================*/

  /* Libération de la structure des positions */
  /*------------------------------------------*/

  if ((*this_champ)->pos != NULL)
    ECS_FREE((*this_champ)->pos);

  /* Libération de la structure des valeurs */
  /*----------------------------------------*/

  if ((*this_champ)->val != NULL)
    ECS_FREE((*this_champ)->val);

  /* Libération du descripteur de champ */
  /*------------------------------------*/

  /* Appel à la fonction de libération d'un descripteur de champ */

  if ((*this_champ)->descr != NULL)
    ecs_descr_chaine__detruit(&((*this_champ)->descr));

  /* Libération de la structure `ecs_champ_t' */
  /*==========================================*/

  ECS_FREE(*this_champ);
}

/*----------------------------------------------------------------------------
 *  Fonction qui convertit, si possible,
 *   le tableau des positions d'un champ en REGLE
 *----------------------------------------------------------------------------*/

void
ecs_champ__pos_en_regle(ecs_champ_t  *this_champ)
{
  size_t      ipos;
  size_t      pos_pas;
  bool        bool_regle;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ != NULL);

  bool_regle = true;

  if (this_champ->pos != NULL && this_champ->nbr > 0) {

    this_champ->pas = 0;

    pos_pas = this_champ->pos[1] - this_champ->pos[0];

    for (ipos = 1; ipos < this_champ->nbr; ipos++) {

      if (this_champ->pos[ipos + 1] - this_champ->pos[ipos] != pos_pas) {
        bool_regle = false;
        break;
      }

    }

    if (bool_regle == true) {

      this_champ->pas = pos_pas;

      ECS_FREE(this_champ->pos);

    }
  }

}

/*----------------------------------------------------------------------------
 *  Fonction qui construit, si nécessaire, un tableau des positions à
 *   partir d'une REGLE.
 *----------------------------------------------------------------------------*/

void
ecs_champ__regle_en_pos(ecs_champ_t  *this_champ)
{
  size_t ipos;
  ecs_size_t *tab_pos;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ != NULL);

  if (this_champ->pos == NULL) {

    ECS_MALLOC(tab_pos, this_champ->nbr + 1, ecs_size_t);

    for (ipos = 0; ipos <= this_champ->nbr; ipos++)
      tab_pos[ipos] = (ipos * this_champ->pas) + 1;

    this_champ->pos = tab_pos;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui libère, si possible, le tableau des positions d'un champ.
 *  Ce tableau ne doit pas avoir été modifié.
 *----------------------------------------------------------------------------*/

void
ecs_champ__libere_pos(ecs_champ_t  *this_champ)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ != NULL);

  if (this_champ->pas != 0 && this_champ->pos != NULL)
    ECS_FREE(this_champ->pos);
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_champ_t' donnée
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_champ__imprime(const ecs_champ_t  *this_champ,
                   size_t              imp_col,
                   size_t              nbr_imp,
                   FILE               *fic_imp)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ   != NULL);
  assert(fic_imp      != NULL);

  imp_col++;

  /* Impression des champs d'information du champ */
  /*----------------------------------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, "nbr_elt", ECS_TYPE_size_t,
                       &this_champ->nbr);

  ecs_fic__imprime_val(fic_imp, imp_col, "pos_pas", ECS_TYPE_size_t,
                       &this_champ->pas);

  ecs_fic__imprime_ptr(fic_imp, imp_col, "pos_tab", this_champ->pos);
  ecs_fic__imprime_ptr(fic_imp, imp_col, "val_tab", this_champ->val);


  /* Impression des positions et des valeurs */
  /*-----------------------------------------*/

  if (this_champ->pos == NULL && this_champ->val != NULL)
    _imprime_pos_pas(fic_imp,
                     this_champ->nbr,
                     this_champ->pas,
                     this_champ->val,
                     nbr_imp);

  else if (this_champ->pos != NULL)
    _imprime_pos_tab(fic_imp,
                     this_champ->nbr,
                     this_champ->pos,
                     this_champ->val,
                     nbr_imp);

  /* Impression de la liste chaînée des descripteurs */
  /*-------------------------------------------------*/

  /* Impression du pointeur sur le descripteur de tête */

  ecs_fic__imprime_val(fic_imp, imp_col, "descr_tete", ECS_TYPE_void,
                       this_champ->descr);

  if (this_champ->descr != NULL) {

    /* Appel à la fonction d'impression d'une chaîne de descripteurs de champ */

    ecs_descr_chaine__imprime(this_champ->descr,
                              imp_col + 1,
                              fic_imp);
  }

}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_champ_t'
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_taille(const ecs_champ_t  *this_champ)
{
  size_t        nbr_val;
  size_t        taille;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_champ == NULL)
    return 0;

  taille = sizeof(*this_champ);

  if (this_champ->pos != NULL)
    taille += (sizeof(ecs_int_t) * (this_champ->nbr + 1));

  if (this_champ->val != NULL) {
    nbr_val = ecs_champ__ret_val_nbr(this_champ);
    taille += (sizeof(ecs_int_t) * nbr_val);
  }

  if (this_champ->descr != NULL)
    taille += ecs_descr_chaine__ret_taille(this_champ->descr);

  return taille;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un champ entièrement réalloué
 *   dont le contenu est copié à partir du champ donné
 *
 *  Le membre donnant le lien sur un champ suivant `l_champ_sui'
 *   n'est pas copié et est mis à `NULL'
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__copie(ecs_champ_t  *champ_init)
{
  size_t        nbr_val;
  ecs_champ_t * this_champ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(champ_init != NULL);

  ECS_MALLOC(this_champ, 1, ecs_champ_t);

  this_champ->nbr = champ_init->nbr;

  this_champ->pas = champ_init->pas;

  if (champ_init->pos != NULL) {
    ECS_MALLOC(this_champ->pos, champ_init->nbr + 1, ecs_size_t);
    memcpy(this_champ->pos,
           champ_init->pos,
           (champ_init->nbr + 1) * sizeof(ecs_size_t));
  }
  else
    this_champ->pos = NULL;

  if (champ_init->val != NULL) {
    nbr_val = ecs_champ__ret_val_nbr(champ_init);
    ECS_MALLOC(this_champ->val, nbr_val, ecs_int_t);
    memcpy(this_champ->val,
           champ_init->val,
           sizeof(ecs_int_t) * nbr_val);
  }
  else
    this_champ->val = NULL;

  this_champ->descr  = ecs_descr_chaine__copie(champ_init->descr);

  return this_champ;
}

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
                              ecs_descr_t          *descr)
{
  ecs_size_t *pos_tab;
  ecs_int_t *val_tab;
  size_t cpt_val;
  size_t ielt;

  ecs_champ_t * this_champ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(nbr_elt != 0);
  assert(tab_elt != NULL);

  /* Allocation de la structure `ecs_champ_t' */
  /*--------------------------------------*/

  ECS_MALLOC(this_champ, 1, ecs_champ_t);

  this_champ->nbr = nbr_elt;

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

  this_champ->pas = 0;
  this_champ->pos = pos_tab;

  ecs_champ__pos_en_regle(this_champ);

  /* Création de la table des valeurs (conteneur) */
  /*----------------------------------------------*/

  this_champ->val = val_tab;

  /* Affectation du descripteur de champ */
  /*-------------------------------------*/

  this_champ->descr = descr;

  return this_champ;
}

/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre d'éléments associés à un champ donné
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_elt_nbr(const ecs_champ_t  *this_champ)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ   != NULL);

  return this_champ->nbr;
}


/*----------------------------------------------------------------------------
 *  Fonction renvoyant le nombre de valeurs associées à un champ donné
 *----------------------------------------------------------------------------*/

size_t
 ecs_champ__ret_val_nbr(const ecs_champ_t  *this_champ)
{
  size_t  nbr_val;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ   != NULL);

  if (this_champ->pos != NULL)
    nbr_val = this_champ->pos[this_champ->nbr] - 1;
  else
    nbr_val = (this_champ->pas * this_champ->nbr);

  return nbr_val;
}

/*----------------------------------------------------------------------------
 *  Fonction retournant le nombre de descripteurs d'un champ donné
 *----------------------------------------------------------------------------*/

size_t
ecs_champ__ret_descr_nbr(const ecs_champ_t  *this_champ)
{
  assert(this_champ != NULL);

  return ecs_descr_chaine__ret_nbr(this_champ->descr);
}

/*----------------------------------------------------------------------------
 *  Fonction libérant un pointeur sur le tableau des positions d'une
 *   structure `ecs_champ_t' donnée.
 *
 *  Si les positions correspondent à une REGLE, le tableau est libéré.
 *   Sinon, il est conservé par la structure ecs_champ_t.
 *----------------------------------------------------------------------------*/

void
ecs_champ__libere_pos_tab(const ecs_champ_t  *this_champ,
                          ecs_size_t         *pos_tab)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ != NULL);
  assert(this_champ->pos == NULL || this_champ->pos == pos_tab);

  if (this_champ->pos == NULL)
    ECS_FREE(pos_tab);
}

/*----------------------------------------------------------------------------
 *  Fonction qui concatène deux champs, et supprime le champ à concaténer
 *----------------------------------------------------------------------------*/

void
ecs_champ__concatene(ecs_champ_t  **this_champ,
                     ecs_champ_t  **concat_champ,
                     size_t         nbr_elt_init,
                     size_t         nbr_elt_concat)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_champ != NULL);
  assert(concat_champ != NULL);

  /* ------------------------------------------------------------------ */
  /* PARTIE REALISEE PENDANT LA FUSION                                  */
  /* ------------------------------------------------------------------ */
  /*    Si le champ n'existe pas dans l'entité de maillage à concaténer */
  /* et si c'est un champ de type "attribut"                            */
  /* (pour lequel tous les éléments doivent être renseignes)            */
  /* on initialise les valeurs de l'attribut pour le champ à concaténer */
  /* ------------------------------------------------------------------ */

  /* Le champ existe dans l'entité de maillage à concaténer */

  if (*concat_champ != NULL) {

    /* On le concatène */

    if (*this_champ != NULL) {
      _champ__concatene(*this_champ, *concat_champ);
      ecs_champ__detruit(concat_champ);
    }
    else {
      ecs_champ__prolonge(*concat_champ,
                          nbr_elt_init,
                          0);
      *this_champ = *concat_champ;
      *concat_champ = NULL;
    }
  }

  /* Le champ n'existe pas dans l'entité de maillage à concaténer */

  else {
    ecs_champ__prolonge(*this_champ,
                        0,
                        nbr_elt_concat);

  } /* Fin : boucle sur les champs de l'entité de maillage receptrice */
}

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
                    size_t        nbr_elt_suiv)
{
  size_t      ipos;
  size_t      nbr_elt_ini;
  ecs_int_t   pos_fin;
  bool        bool_pos_ini;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_champ == NULL)
    return;

  ecs_champ__regle_en_pos(this_champ);

  nbr_elt_ini = this_champ->nbr;

#if defined(SX) && defined(_SX) /* NEC SX compiler may bug on : operator */
  if (this_champ->pos != NULL)
    bool_pos_ini = true;
  else
    bool_pos_ini = false;
#else
  bool_pos_ini = (this_champ->pos != NULL) ? true : false;
#endif

  /* Mise à jour du nombre d'éléments */

  this_champ->nbr += nbr_elt_prec + nbr_elt_suiv;

  /* Si le champ n'est pas déjà vide, la table des positions ne
     correspondra pas à une REGLE, et devra être construite */

  if (this_champ->pos != NULL || this_champ->pas != 0) {

    ECS_REALLOC(this_champ->pos, this_champ->nbr + 1, ecs_size_t);

    if (bool_pos_ini == true) {
      memmove(this_champ->pos + (nbr_elt_prec * sizeof(ecs_int_t)),
              this_champ->pos,
              (this_champ->nbr + 1) * sizeof(ecs_int_t));
    }
    else {
      for (ipos = 0; ipos <= nbr_elt_ini; ipos++)
        this_champ->pos[nbr_elt_prec + ipos]
          = (ipos * this_champ->pas) + 1;
    }

    this_champ->pos[0] = 1;

    for (ipos = 0; ipos < nbr_elt_prec; ipos++)
      this_champ->pos[ipos + 1] = 1;

    pos_fin = this_champ->pos[nbr_elt_prec + nbr_elt_ini];
    for (ipos = 0; ipos < nbr_elt_suiv; ipos++)
      this_champ->pos[nbr_elt_prec + nbr_elt_ini + ipos + 1] = pos_fin;

  }

  ecs_champ__pos_en_regle(this_champ);
}

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
                          const ecs_tab_int_t   vect_transf)
{
  ecs_champ_t  *champ_ref;
  size_t       nbr_val_ref;
  ecs_int_t    ival;
  size_t       ielt_ref;
  size_t       ielt_transf;
  size_t       ival_ref;
  ecs_int_t    elt_nbr_val;
  ecs_int_t    pre_pos;
  ecs_int_t    cpt_val;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_champ == NULL)
    return;

  ecs_champ__regle_en_pos(this_champ);

  assert(vect_transf.nbr == nbr_elt_ref);
  assert(this_champ->nbr == nbr_elt_ref);

  nbr_val_ref = ecs_champ__ret_val_nbr(this_champ);

  champ_ref = ecs_champ__alloue(nbr_elt_ref, nbr_val_ref);

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref + 1; ielt_ref++)
    champ_ref->pos[ielt_ref] = this_champ->pos[ielt_ref];

  for (ival_ref = 0; ival_ref < nbr_val_ref; ival_ref++)
    champ_ref->val[ival_ref] = this_champ->val[ival_ref];

  this_champ->pos[0] = 1;

  cpt_val = 0;

  for (ielt_transf = 0; ielt_transf < nbr_elt_ref; ielt_transf++) {

    pre_pos = champ_ref->pos[vect_transf.val[ielt_transf]];

    elt_nbr_val
      = champ_ref->pos[vect_transf.val[ielt_transf] + 1] - pre_pos;

    this_champ->pos[ielt_transf + 1]
      = this_champ->pos[ielt_transf] + elt_nbr_val;

    for (ival = 0; ival < elt_nbr_val; ival++) {

      this_champ->val[cpt_val++]
        = champ_ref->val[pre_pos - 1 + ival];

    }
  }

  ecs_champ__detruit(&champ_ref);

  ecs_champ__pos_en_regle(this_champ);
}

/*----------------------------------------------------------------------------
 *  Fonction qui incrémente les valeurs d'un champ donné
 *   d'une constante donnée
 *----------------------------------------------------------------------------*/

void
ecs_champ__incremente_val(ecs_champ_t      *this_champ,
                          const ecs_int_t   increment)
{
  size_t     nbr_val;
  size_t     ival;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_champ == NULL)
    return;

  nbr_val = ecs_champ__ret_val_nbr(this_champ);

  for (ival = 0; ival < nbr_val; ival++) {

    if (this_champ->val[ival] > 0)
      this_champ->val[ival] += increment;
    else if (this_champ->val[ival] < 0)
      this_champ->val[ival] -= increment;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un vecteur indexé
 *   en appliquant directement le vecteur de transformation donné
 *   sur les valeurs associées à ses éléments
 *----------------------------------------------------------------------------*/

void
ecs_champ__renumerote(ecs_champ_t          *this_champ,
                      const ecs_tab_int_t   vect_transf,
                      const ecs_tab_int_t   signe_elt)
{
  ecs_champ_t *champ_ref;
  size_t       nbr_val_ref;
  ecs_int_t    val_ref;
  ecs_int_t    sgn_ref;
  size_t       ival_ref;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_champ == NULL)
    return;

  nbr_val_ref = ecs_champ__ret_val_nbr(this_champ);

  champ_ref = ecs_champ__alloue(0, nbr_val_ref);

  for (ival_ref = 0; ival_ref < nbr_val_ref; ival_ref++)
    champ_ref->val[ival_ref] = this_champ->val[ival_ref];

  for (ival_ref = 0; ival_ref < nbr_val_ref; ival_ref++) {

    val_ref = ECS_ABS(champ_ref->val[ival_ref]);
    sgn_ref = champ_ref->val[ival_ref] / val_ref;

    this_champ->val[ival_ref]  = vect_transf.val[val_ref - 1] + 1;

    this_champ->val[ival_ref] *= signe_elt.val[val_ref - 1] * sgn_ref;

  }

  ecs_champ__detruit(&champ_ref);
}

/*----------------------------------------------------------------------------
 *  Fonction qui détermine un nouveau champ à partir d'un champ de référence
 *   en extrayant de ce dernier les éléments sélectionnés
 *   par le tableau de booléens
 *----------------------------------------------------------------------------*/

ecs_champ_t *
ecs_champ__extrait(ecs_champ_t            *champ_ref,
                   const ecs_tab_bool_t    bool_elt_select)
{
  size_t        cpt_elt_new;
  size_t        cpt_val_new;
  size_t        cpt_val_old_new;
  size_t        nbr_elt_ref;
  size_t        nbr_val_ref;
  size_t        pos_ref_inf;
  size_t        pos_ref_sup;

  size_t        ielt_ref;
  size_t        ipos_ref;

  ecs_champ_t  *champ_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_champ__regle_en_pos(champ_ref);

  nbr_elt_ref = champ_ref->nbr;
  nbr_val_ref = ecs_champ__ret_val_nbr(champ_ref);

  champ_new = ecs_champ__alloue(nbr_elt_ref,
                                nbr_val_ref);

  /* Extraction du champ */

  nbr_elt_ref = champ_ref->nbr;
  nbr_val_ref = ecs_champ__ret_val_nbr(champ_ref);

  champ_new->pos[0] = 1;
  cpt_elt_new           = 0;
  cpt_val_new           = 0;
  cpt_val_old_new       = 0;

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref; ielt_ref++) {

    if (bool_elt_select.val[ielt_ref] == true) {

      /* L'élément est à extraire */

      pos_ref_inf = champ_ref->pos[ielt_ref    ] - 1;
      pos_ref_sup = champ_ref->pos[ielt_ref + 1] - 1;

      for (ipos_ref = pos_ref_inf; ipos_ref < pos_ref_sup; ipos_ref++)
        champ_new->val[cpt_val_new++] = champ_ref->val[ipos_ref];

      champ_new->pos[cpt_elt_new + 1] = cpt_val_new + 1;

      cpt_elt_new++;

    }
  }

  ECS_REALLOC(champ_new->pos, cpt_elt_new + 1, ecs_size_t);
  ECS_REALLOC(champ_new->val, cpt_val_new,     ecs_int_t);

  champ_new->nbr = cpt_elt_new;

  ecs_champ__pos_en_regle(champ_ref);

  return champ_new;
}

/*----------------------------------------------------------------------------*/
