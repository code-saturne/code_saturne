/*============================================================================
 *  Définitions des fonctions de base
 *   associées à une liste chaînée de structures `ecs_descr_t' décrivant
 *   un descripteur de table
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
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr_chaine.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui supprime un noeud donné
 *   dans une liste chaînée de descripteurs dont la tête est donnée
 *
 *  La tête de la liste qui n'est pas modifiée
 *   sauf si le noeud à supprimer est aussi la tête de la liste !
 *----------------------------------------------------------------------------*/

static int
_ecs_descr_chaine_compare(const void  *descr1,
                          const void  *descr2)
{
  const ecs_descr_t *d1 = descr1;
  const ecs_descr_t *d2 = descr2;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (d1->nom != NULL && d2->nom != NULL)
    return strcmp(d1->nom, d2->nom);

  else if (d1->nom != NULL && d2->nom == NULL)
    return 1;

  else if (d1->nom == NULL && d2->nom != NULL)
    return -1;

  return 0;
}

/*----------------------------------------------------------------------------
 *  Fonction qui supprime un noeud donné
 *   dans une liste chaînée de descripteurs dont la tête est donnée
 *
 *  La tête de la liste qui n'est pas modifiée
 *   sauf si le noeud à supprimer est aussi la tête de la liste !
 *----------------------------------------------------------------------------*/

static void
ecs_loc_descr_chaine__supprime(ecs_descr_t  **descr_noeud,
                               ecs_descr_t  **descr_tete)
{
  ecs_descr_t *loc_descr_sui;
  ecs_descr_t *ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(*descr_tete != NULL);
  assert(*descr_noeud != NULL);

  if ((*descr_noeud) == *descr_tete) {

    /* Le descripteur noeud à supprimer est la tête de liste des descripteurs */
    /* Le descripteur suivant est la nouvelle tête                            */

    *descr_tete = (*descr_tete)->l_descr_sui;

  }
  else {

    /* On recherche le noeud qui précède le noeud à supprimer */

    loc_descr_sui = (*descr_noeud)->l_descr_sui;

    for (ptr_descr = *descr_tete;
         ptr_descr != NULL && ptr_descr->l_descr_sui != (*descr_noeud);
         ptr_descr = ptr_descr->l_descr_sui);

    /* Le noeud à supprimer doit être contenu dans la liste */
    /* `*descr_tete' est la tête                            */
    assert(ptr_descr != NULL);

    ptr_descr->l_descr_sui = loc_descr_sui;

  } /* Fin else : le descripteur à supprimer n'est pas le descripteur de tête */

  /* Libération du descripteur correspondant au noeud à supprimer */

  *descr_noeud = ecs_descr__detruit(*descr_noeud);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction libérant la portion d'une liste chaînée de descripteurs
 *   à partir d'un noeud dont le pointeur est donné en argument.
 *  Le noeud est à NULL au retour de la fonction
 *----------------------------------------------------------------------------*/

void
ecs_descr_chaine__detruit(ecs_descr_t  **descr_noeud)
{
  if (*descr_noeud != NULL) {

    ecs_descr_chaine__detruit(&(*descr_noeud)->l_descr_sui);

    *descr_noeud = ecs_descr__detruit(*descr_noeud);

  }
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant à partir d'un noeud `ecs_descr_t' donné
 *   une liste chaînée de tables
 *   sur le flux décrit par la structure `ecs_file_t'
 *----------------------------------------------------------------------------*/

void ecs_descr_chaine__imprime(const ecs_descr_t  *descr_noeud,
                               int                 imp_col,
                               FILE               *fic_imp)
{
#define ECS_FCT_IMP_DESCR_NOEUD       "descr"

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (descr_noeud != NULL) {

    ecs_descr_chaine__imprime(descr_noeud->l_descr_sui,
                              imp_col,
                              fic_imp);

    ecs_fic__imprime_ptr(fic_imp, imp_col, ECS_FCT_IMP_DESCR_NOEUD,
                         (const void *)descr_noeud);

    ecs_descr__imprime(descr_noeud,
                       imp_col,
                       fic_imp);

  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets
 *   d'une chaîné de structures `ecs_descr_t'
 *----------------------------------------------------------------------------*/

float
ecs_descr_chaine__ret_taille(const ecs_descr_t  *descr_noeud)
{

  float         taille;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  taille = 0.;

  if (descr_noeud != NULL) {

    taille +=
      ecs_descr_chaine__ret_taille(descr_noeud->l_descr_sui);

    taille += ecs_descr__ret_taille(descr_noeud);

  }

  return taille;
}

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
                         ecs_descr_t   *descr_concat_tete)
{
  ecs_int_t   inum;

  ecs_descr_t  *loc_descr_prec;
  ecs_descr_t  *ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_tete   != NULL);

  if (*descr_tete != NULL) {

    /* On va à la fin de la chaîne réceptrice */

    for (ptr_descr = *descr_tete;
         ptr_descr != NULL;
         ptr_descr = ptr_descr->l_descr_sui  )
      loc_descr_prec = ptr_descr;

    /* On ajoute le lien avec le début de la chaîne à concaténer */

    loc_descr_prec->l_descr_sui = descr_concat_tete;

    /* Les numéros des descripteurs de la liste à concaténer sont incrémentés */
    /*  à partir du nombre de descripteur de la liste réceptrice              */

    for (ptr_descr = descr_concat_tete, inum = 1;
         ptr_descr != NULL;
         ptr_descr = ptr_descr->l_descr_sui, inum++    )
      ptr_descr->num = loc_descr_prec->num + inum;

  }
  else {

    *descr_tete = descr_concat_tete;

  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre de descripteurs
 *   de la liste chaînée des descripteurs dont la tête est donnée
 *----------------------------------------------------------------------------*/

int
ecs_descr_chaine__ret_nbr(const ecs_descr_t  *descr_tete)
{
  const ecs_descr_t  *ptr_descr;
  ecs_int_t           nbr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_tete != NULL);

  nbr_descr = 0;

  for (ptr_descr = descr_tete;
       ptr_descr != NULL;
       ptr_descr = ptr_descr->l_descr_sui) {

    nbr_descr++;

  }

  return nbr_descr;
}

/*----------------------------------------------------------------------------
 *  Fonction qui copie une liste chaînée de descripteurs
 *   dont la tête est donnée
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr_chaine__copie(ecs_descr_t  *descr_tete)
{
  ecs_descr_t  *descr_copie;
  ecs_descr_t  *descr_tete_copie;
  ecs_descr_t  *ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  descr_tete_copie = NULL;

  for (ptr_descr  = descr_tete;
       ptr_descr != NULL;
       ptr_descr  = ptr_descr->l_descr_sui  ) {

    descr_copie = ecs_descr__copie(ptr_descr);

    ecs_descr_chaine__ajoute(&descr_tete_copie,
                             descr_copie);

  }

  return descr_tete_copie;
}

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
                            ecs_descr_t  **descr_concat_tete)
{
  ecs_int_t     cpt_descr;

  ecs_descr_t    *ptr_descr_concat;
  ecs_descr_t    *ptr_descr_recept;

  ecs_tab_int_t   tab_renum_descr_concat;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(*descr_concat_tete != NULL);

  /* On repère le plus grand numéro de descripteur de la chaîne à concaténer */

  for (ptr_descr_concat = *descr_concat_tete, cpt_descr = 0;
       ptr_descr_concat != NULL;
       ptr_descr_concat = ptr_descr_concat->l_descr_sui, cpt_descr++);

  ECS_MALLOC(tab_renum_descr_concat.val, cpt_descr, ecs_int_t);
  tab_renum_descr_concat.nbr = cpt_descr;

  /* On repère le plus grand numéro de descripteur de la chaîne réceptrice */

  for (ptr_descr_recept = *descr_recept_tete, cpt_descr = 0;
       ptr_descr_recept != NULL;
       ptr_descr_recept = ptr_descr_recept->l_descr_sui, cpt_descr++);

  /* On parcourt les descripteurs de la chaîne réceptrice   :               */
  /* on parcourt les descripteurs de la chaîne à concaténer :               */
  /*  - si un descripteur n'a pas un descripteur identique                  */
  /*     dans la chaîne à concaténer                                        */
  /*     -> le descripteur de la liste à concaténer prend un nouveau numéro */
  /*  - sinon                                                               */
  /*     -> il est supprime de la chaîne à concaténer                       */

  for (ptr_descr_recept = *descr_recept_tete;
       ptr_descr_recept != NULL;
       ptr_descr_recept = ptr_descr_recept->l_descr_sui  ) {

    ptr_descr_concat = *descr_concat_tete;
    while (ptr_descr_concat != NULL                          &&
           ecs_descr__compare(ptr_descr_recept,
                              ptr_descr_concat) == false    )
      ptr_descr_concat = ptr_descr_concat->l_descr_sui;

    if (ptr_descr_concat != NULL) {

      /* Il y a un descripteur identique dans la chaîne à concaténer */

      tab_renum_descr_concat.val[ptr_descr_concat->num - 1]
        = ptr_descr_recept->num - 1;

      ecs_loc_descr_chaine__supprime(&ptr_descr_concat,
                                 descr_concat_tete);
    }
  }

  if (*descr_concat_tete != NULL) {

    /* Les numéros des descripteurs restants de la chaîne à concaténer */
    /*  sont modifies                                                  */

    for (ptr_descr_concat = *descr_concat_tete;
         ptr_descr_concat != NULL;
         ptr_descr_concat = ptr_descr_concat->l_descr_sui) {

      tab_renum_descr_concat.val[ptr_descr_concat->num - 1] = cpt_descr;
      cpt_descr++;
    }

    /* On ajoute à la fin de la chaîne réceptrice,          */
    /*  les descripteurs restants de la chaîne à concaténer */

    ecs_descr_chaine__ajoute(descr_recept_tete,
                             *descr_concat_tete);
  }

  return tab_renum_descr_concat;
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche les contenus des descripteurs
 *   de la liste chaînée des descripteurs dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_descr_chaine__affiche(ecs_descr_t  *descr_tete,
                          int           decal)
{
  ecs_descr_t  *ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ptr_descr  = descr_tete;
       ptr_descr != NULL;
       ptr_descr  = ptr_descr->l_descr_sui) {

    ecs_descr__affiche(ptr_descr,
                       decal);
  }
}

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
                              int           num)
{
  ecs_descr_t  *ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_tete != NULL);

  ptr_descr = descr_tete;

  while (ptr_descr != NULL && ptr_descr->num != num)
    ptr_descr = ptr_descr->l_descr_sui;

  return ptr_descr;
}

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
                             const ecs_descr_t  *descr_rech)
{
  ecs_descr_t  *ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_tete != NULL && descr_rech != NULL);

  ptr_descr = descr_tete;

  while (   ptr_descr != NULL
         && ecs_descr__compare(ptr_descr, descr_rech ) == false)
    ptr_descr = ptr_descr->l_descr_sui;

  if (ptr_descr != NULL) {

    return ptr_descr->num;

  }
  else {

    return ECS_DESCR_NUM_NUL;

  }
}

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
                             const ecs_tab_int_t   vect_transf)
{
  size_t        inum;

  ecs_descr_t  *descr_new;
  ecs_descr_t  *descr_tete_new;
  ecs_descr_t  *ptr_descr;        /* Pointeur de boucle sur les descripteurs */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  descr_tete_new = NULL;

  for (inum = 0; inum < vect_transf.nbr; inum++) {

    ptr_descr = descr_tete;
    while (ptr_descr != NULL && vect_transf.val[inum] != (ptr_descr->num - 1))
      ptr_descr = ptr_descr->l_descr_sui;

    assert(ptr_descr != NULL);

    descr_new = ecs_descr__cree(ECS_DESCR_IDE_NUL,
                                ptr_descr->nom);

    ecs_descr_chaine__ajoute(&descr_tete_new,
                             descr_new);

  }

  return descr_tete_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre et la liste des pointeurs sur les noms
 *   des descripteurs de type groupe d'une liste chaînée dont la tête est
 *   donnée en argument
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_descr_chaine__ret_nom(ecs_descr_t   *descr_tete)
{
  ecs_int_t       cpt_descr;
  ecs_descr_t    *ptr_descr;

  ecs_tab_char_t  tab_nom_descr_chaine;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_tete != NULL);

  /* Comptage des descripteurs de type groupe */

  cpt_descr = 0;

  for (ptr_descr  = descr_tete;
       ptr_descr != NULL;
       ptr_descr  = ptr_descr->l_descr_sui) {
    assert(ptr_descr->nom != NULL);
    cpt_descr += 1;
  }

  ECS_MALLOC(tab_nom_descr_chaine.val, cpt_descr, char *);
  tab_nom_descr_chaine.nbr = cpt_descr;

  /* Construction de la liste des descripteurs de type groupe */

  cpt_descr = 0;

  for (ptr_descr  = descr_tete;
       ptr_descr != NULL;
       ptr_descr  = ptr_descr->l_descr_sui) {
    tab_nom_descr_chaine.val[cpt_descr] = ptr_descr->nom;
    cpt_descr += 1;
  }

  return tab_nom_descr_chaine;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la liste des références des descripteurs
 *   de la liste chaînée des descripteurs dont la tête est donnée
 *----------------------------------------------------------------------------*/

ecs_descr_t **
ecs_descr_chaine__ret_ref(ecs_descr_t  *descr_tete,
                          int          *nbr_descr)
{
  ecs_int_t     idescr;

  ecs_descr_t * * liste_ref_descr;
  ecs_descr_t   * ptr_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Détermination du nombre de descripteurs de la liste chaînée */
  /*-------------------------------------------------------------*/

  *nbr_descr = ecs_descr_chaine__ret_nbr(descr_tete);

  ECS_MALLOC(liste_ref_descr, *nbr_descr, ecs_descr_t *);

  idescr = 0;

  for (ptr_descr  = descr_tete;
       ptr_descr != NULL;
       ptr_descr  = ptr_descr->l_descr_sui) {

    liste_ref_descr[idescr++] = ptr_descr;

  }

  return liste_ref_descr;
}

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
ecs_descr_chaine__trie(ecs_descr_t  *descr_tete)
{
  size_t        nbr_descr;
  size_t        i, j;

  ecs_descr_t  *ptr_descr;
  ecs_descr_t  *tab_descr;

  ecs_tab_int_t  tab_renum;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nbr_descr = ecs_descr_chaine__ret_nbr(descr_tete);

  ECS_MALLOC(tab_descr, nbr_descr, ecs_descr_t);

  for (ptr_descr = descr_tete, i = 0;
       ptr_descr != NULL;
       ptr_descr = ptr_descr->l_descr_sui, i++)
    tab_descr[i] = *ptr_descr;

  /* Sort descriptors */

  qsort(tab_descr,
        nbr_descr,
        sizeof(ecs_descr_t),
        _ecs_descr_chaine_compare);

  tab_renum.nbr = nbr_descr;
  ECS_MALLOC(tab_renum.val, tab_renum.nbr, ecs_int_t);

  for (ptr_descr = descr_tete, i = 0;
       ptr_descr != NULL;
       ptr_descr = ptr_descr->l_descr_sui, i++) {
    for (j = 0;
         j < nbr_descr && ptr_descr->num != tab_descr[j].num;
         j++);
    assert(j < nbr_descr);
    assert(ptr_descr->num <= (int)nbr_descr);
    tab_renum.val[tab_descr[j].num - 1] = j + 1;
    ptr_descr->num = j + 1;
    ptr_descr->nom = tab_descr[j].nom;
  }

  ECS_FREE(tab_descr);

  return tab_renum;
}

/*----------------------------------------------------------------------------*/

