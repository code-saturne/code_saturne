/*============================================================================
 *  Définitions des fonctions de base
 *   associées à une liste chaînée de structures `ecs_famille_t' décrivant
 *   une famille
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
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h> /* strlen() */


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

#include "ecs_descr.h"
#include "ecs_descr_chaine.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille_chaine.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille_priv.h"



/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui crée une liste chaînée de familles à partir :
 *   - des définitions de chaque famille en fonction des numéros de descripteur
 *   - de la liste chaînée des descripteurs
 *  La fonction renvoie la tête de la liste chaînée
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille_chaine__cree(ecs_int_t        **def_fam_descr,
                         const ecs_int_t   *nbr_descr_fam,
                         int                num_fam_deb,
                         int                nbr_fam,
                         ecs_descr_t       *descr_tete)
{
  ecs_int_t  ifam;
  ecs_int_t  idescr;

  ecs_descr_t  *descr_fam;
  ecs_descr_t  *descr_tete_fam;
  ecs_descr_t  *ptr_descr;

  ecs_famille_t  *fam ;
  ecs_famille_t  *fam_tete;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  fam_tete = NULL;

  for (ifam = 0; ifam < nbr_fam; ifam++) {

    descr_tete_fam = NULL;

    for (idescr = 0; idescr < nbr_descr_fam[ifam]; idescr++) {

      ptr_descr = ecs_descr_chaine__cherche_num(descr_tete,
                                                def_fam_descr[ifam][idescr]);

      assert(ptr_descr != NULL);

      descr_fam = ecs_descr__copie(ptr_descr);

      ecs_descr_chaine__ajoute(&descr_tete_fam,
                               descr_fam);
    }

    fam = ecs_famille__cree(ifam + num_fam_deb, descr_tete_fam);

    ecs_famille_chaine__ajoute(&fam_tete, fam);

  }

  return fam_tete;
}

/*----------------------------------------------------------------------------
 *  Fonction libérant la portion d'une liste chaînée de familles
 *   à partir d'un noeud dont le pointeur est donné en argument.
 *  Le noeud est à NULL au retour de la fonction
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__detruit(ecs_famille_t  **this_fam_noeud)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (*this_fam_noeud != NULL) {

    ecs_famille_chaine__detruit(&(*this_fam_noeud)->l_famille_sui);

    *this_fam_noeud = ecs_famille__detruit(*this_fam_noeud);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant à partir d'un noeud `ecs_famille_t' donné
 *   une liste chaînée de tables
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__imprime(const ecs_famille_t  *this_fam_noeud,
                            ecs_int_t             imp_col,
                            FILE                 *fic_imp)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_fam_noeud != NULL) {

    ecs_famille_chaine__imprime(this_fam_noeud->l_famille_sui,
                                imp_col,
                                fic_imp);

    ecs_fic__imprime_ptr(fic_imp, imp_col, "famille",
                         (const void *)this_fam_noeud);

    ecs_famille__imprime(this_fam_noeud,
                         imp_col,
                         fic_imp);

  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets
 *   d'une chaîne de structures `ecs_famille_t'
 *----------------------------------------------------------------------------*/

float
ecs_famille_chaine__ret_taille(const ecs_famille_t  *this_fam_noeud)
{
  float  taille = 0.;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_fam_noeud != NULL) {

    taille += ecs_famille_chaine__ret_taille(this_fam_noeud->l_famille_sui);

    taille += ecs_famille__ret_taille(this_fam_noeud);
  }

  return taille;
}

/*----------------------------------------------------------------------------
 *  Fonction qui ajoute à la fin d'une liste chaînée de familles
 *   réceptrice dont la tête est donnée,
 *   une liste chaînée de familles à concaténer dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__ajoute(ecs_famille_t  **this_fam_tete,
                           ecs_famille_t   *fam_concat_tete)
{
  ecs_famille_t  *loc_fam_prec;
  ecs_famille_t  *ptr_fam ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_fam_tete != NULL);

  if (*this_fam_tete != NULL) {

    /* On va à la fin de la chaîne réceptrice */

    for (ptr_fam = *this_fam_tete;
         ptr_fam != NULL;
         ptr_fam = ptr_fam->l_famille_sui  )
      loc_fam_prec = ptr_fam;

    /* On ajoute le lien avec le début de la chaîne à concaténer */

    loc_fam_prec->l_famille_sui = fam_concat_tete;

  }
  else {

    *this_fam_tete = fam_concat_tete;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche la définition de la famille de numéro donné
 *   à partir de la liste chaînée des familles dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__affiche(const ecs_int_t   num_fam,
                            ecs_famille_t    *fam_tete)
{
  int   ifam;
  ecs_famille_t  *ptr_fam;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(fam_tete != NULL);

  ifam = 1;
  ptr_fam = fam_tete;
  while (ptr_fam != NULL && ifam++ != num_fam)
    ptr_fam = ptr_fam->l_famille_sui;

  assert(ptr_fam != NULL);

  ecs_famille__affiche(ptr_fam);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre de familles
 *   de la liste chaînée des familles dont la tête est donnée
 *----------------------------------------------------------------------------*/

int
ecs_famille_chaine__ret_nbr(const ecs_famille_t  *fam_tete)
{
  int  nbr_fam;

  const ecs_famille_t  *ptr_fam;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Détermination du nombre de familles de la liste chaînée */
  /*---------------------------------------------------------*/

  nbr_fam = 0;

  for (ptr_fam  = fam_tete;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    nbr_fam++;
  }

  return nbr_fam;
}

/*----------------------------------------------------------------------------
 *  Fonction qui copie une liste chaînée de familles
 *   dont la tête est donnée
 *----------------------------------------------------------------------------*/

ecs_famille_t  *
ecs_famille_chaine__copie(ecs_famille_t  *famille_tete)
{
  ecs_famille_t  *famille_copie;
  ecs_famille_t  *famille_tete_copie;
  ecs_famille_t  *ptr_famille;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  famille_tete_copie = NULL;

  for (ptr_famille  = famille_tete;
       ptr_famille != NULL;
       ptr_famille  = ptr_famille->l_famille_sui   ) {

    famille_copie = ecs_famille__copie(ptr_famille);

    ecs_famille_chaine__ajoute(&famille_tete_copie,
                               famille_copie);

  }

  return famille_tete_copie;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie pour chaque numéro de famille
 *   le nombre et une liste de pointeurs sur les noms des identificateurs
 *   de type groupe des descripteurs de la famille
 *----------------------------------------------------------------------------*/

ecs_tab_char_t *
ecs_famille_chaine__ret_nom(ecs_famille_t   *fam_tete)
{
  int  ifam;
  int  nbr_fam;

  ecs_famille_t   *ptr_fam;

  ecs_tab_char_t  *tab_propr_fam;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(fam_tete    != NULL);

  /* Détermination du nombre de familles de la liste chaînée */
  /*---------------------------------------------------------*/

  nbr_fam = ecs_famille_chaine__ret_nbr(fam_tete);

  /* Détermination des groupes pour chaque famille de la liste chaînée */
  /*-------------------------------------------------------------------*/

  ECS_MALLOC(tab_propr_fam, nbr_fam, ecs_tab_char_t);

  for (ptr_fam  = fam_tete, ifam = 0;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui, ifam++) {

    tab_propr_fam[ifam] = ecs_famille__ret_nom(ptr_fam);

  }

  return tab_propr_fam;
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit une liste chaînée de descripteurs
 *   pour chaque numéro de famille contenu dans le tableau donné
 *   et à partir de la liste chaînée des familles
 *
 *  Cette fonction détermine aussi le tableau donnant pour chaque famille
 *   la liste des numéros de descripteurs  associes à la famille
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__cree_descr(ecs_famille_t   *famille,
                               ecs_tab_int_t    tab_fam,
                               ecs_descr_t    **descr_tete_att,
                               ecs_tab_int_t   *tab_att_fam,
                               int             *nbr_max_att_fam)
{
  size_t        ifam;
  size_t        iatt;
  ecs_int_t     ind_fam;

  ecs_descr_t    *descr_tete_att_fam;

  ecs_famille_t  *ptr_fam;

  ecs_tab_int_t   tab_renum_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(famille != NULL);

  *nbr_max_att_fam = 0;
  *descr_tete_att  = NULL;

  for (ifam = 0; ifam < tab_fam.nbr; ifam++) {

    ind_fam = tab_fam.val[ifam] - 1;

    ptr_fam = famille;
    while (ptr_fam != NULL && ECS_ABS(ptr_fam->num) != tab_fam.val[ifam])
      ptr_fam = ptr_fam->l_famille_sui;

    assert(ptr_fam != NULL);

    descr_tete_att_fam = ecs_descr_chaine__copie(ptr_fam->descr);

    if (descr_tete_att_fam != NULL) {

      tab_renum_descr = ecs_descr_chaine__concatene(descr_tete_att,
                                                    &descr_tete_att_fam);

      tab_att_fam[ind_fam].nbr = tab_renum_descr.nbr;
      ECS_MALLOC(tab_att_fam[ind_fam].val,
                 tab_att_fam[ind_fam].nbr, ecs_int_t);

      for (iatt = 0; iatt < tab_renum_descr.nbr; iatt++)
        tab_att_fam[ind_fam].val[iatt] = tab_renum_descr.val[iatt] + 1;

      ECS_FREE(tab_renum_descr.val);

      *nbr_max_att_fam = ECS_MAX(*nbr_max_att_fam,
                                 (int)(tab_renum_descr.nbr));

    }
    else {

      tab_att_fam[ind_fam].val = NULL;
      tab_att_fam[ind_fam].nbr = 0;

    }
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un tableau marquant chaque numéro de famille.
 *
 *  La libération du tableau est à la charge du code appelant
 *----------------------------------------------------------------------------*/

bool *
ecs_famille_chaine__indic_fam_att(const ecs_famille_t  *fam_tete)
{
  int       ifam;
  int       num_fam_max;

  const ecs_famille_t  *ptr_fam;

  bool  *tab_select = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Détermination du numéro de famille maximal de la liste chaînée */
  /*----------------------------------------------------------------*/

  num_fam_max = 0;

  for (ptr_fam  = fam_tete;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    if (ptr_fam->num > num_fam_max)
      num_fam_max = ptr_fam->num;
  }

  ECS_MALLOC(tab_select, num_fam_max + 1, bool);

  for (ifam = 0; ifam < num_fam_max + 1; ifam++)
    tab_select[ifam] = false;

  /* Marquage des familles correspondant aux critères de sélection */
  /*---------------------------------------------------------------*/

  for (ptr_fam  = fam_tete;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    tab_select[ptr_fam->num] = true;
  }

  return tab_select;
}

/*----------------------------------------------------------------------------*/


