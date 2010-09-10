/*============================================================================
 *  Définitions des fonctions de base
 *   associées à une liste chaînée de structures `ecs_famille_t' décrivant
 *   une famille
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

/*----------------------------------------------------------------------------
 *  Fonction qui détermine la liste des attributs référencés par une
 *   liste chaînée de familles dont la tête est donnée, ainsi que
 *   la liste des familles correspondant à chaque attribut
 *
 *  Le tableau optionnel nbr_elt_fam permet de limiter éventuellement
 *   la liste des attributs à ceux qu sont effectivement référencés
 *----------------------------------------------------------------------------*/

static void
ecs_loc_famille_chaine__fam_att(const ecs_famille_t   *fam_tete,
                                const ecs_int_t       *nbr_elt_fam,
                                ecs_tab_int_t         *tab_descr,
                                ecs_tab_int_t        **tab_fam_descr,
                                ecs_descr_t         ***liste_ref_descr)
{
  ecs_descr_typ_t     descr_typ;

  ecs_int_t           cpt_descr;
  ecs_int_t           cpt_fam;
  ecs_int_t           icoul;
  ecs_int_t           idescr;
  ecs_int_t           ifam;
  ecs_int_t           igrp;
  ecs_int_t           nbr_fam;
  ecs_int_t           nbr_max_descr;
  ecs_int_t           nbr_ref_descr;
  int                *liste_nbr_descr_fam;
  ecs_int_t          *liste_num_descr_coul;
  ecs_int_t          *liste_num_descr_grp;

  ecs_descr_t      ***liste_ref_descr_fam;

  ecs_tab_int_t       tab_couleur;
  ecs_tab_int_t       tab_groupe ;
  ecs_tab_int_t       tab_renum_coul ;
  ecs_tab_int_t       tab_renum_grp;

  const ecs_famille_t  *ptr_fam;

  ecs_tab_int_t    _tab_descr;
  ecs_tab_int_t   *_tab_fam_descr;
  ecs_descr_t    **_liste_ref_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Détermination du nombre de familles de la liste chaînée */
  /*---------------------------------------------------------*/

  for (ptr_fam  = fam_tete, cpt_fam = 0;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    if (nbr_elt_fam != NULL && nbr_elt_fam[ptr_fam->num] == 0)
      continue;

    cpt_fam++;
  }

  nbr_fam = cpt_fam;

  ECS_MALLOC(liste_ref_descr_fam, nbr_fam, ecs_descr_t * *);
  ECS_MALLOC(liste_nbr_descr_fam, nbr_fam, int);

  nbr_max_descr = 0;

  for (ptr_fam  = fam_tete, cpt_fam = 0;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    if (nbr_elt_fam != NULL && nbr_elt_fam[ptr_fam->num] == 0)
      continue;

    liste_ref_descr_fam[cpt_fam]
      = ecs_descr_chaine__ret_ref(ptr_fam->descr,
                                  &liste_nbr_descr_fam[cpt_fam]);

    nbr_max_descr += liste_nbr_descr_fam[cpt_fam];

    cpt_fam++;

  }

  /* Construction du tableau donnant pour chaque descripteur */
  /*  les numéros de famille auxquelles il appartient        */
  /*---------------------------------------------------------*/

  ECS_MALLOC(_tab_fam_descr, nbr_max_descr, ecs_tab_int_t);
  _liste_ref_descr = NULL;
  nbr_ref_descr = 0;

  for (ptr_fam  = fam_tete, cpt_fam = 0;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    if (nbr_elt_fam != NULL && nbr_elt_fam[ptr_fam->num] == 0)
      continue;

    for (idescr = 0; idescr < liste_nbr_descr_fam[cpt_fam]; idescr++) {

      cpt_descr = 0;
      while (   cpt_descr < nbr_ref_descr
             && (   ecs_descr__compare(_liste_ref_descr[cpt_descr],
                                       liste_ref_descr_fam[cpt_fam][idescr])
                 == false))
        cpt_descr++;

      if (cpt_descr == nbr_ref_descr) {

        /* Le descripteur n'a pas déjà été rencontre */
        /* On le rajoute à la liste */

        ECS_REALLOC(_liste_ref_descr, nbr_ref_descr + 1, ecs_descr_t *);

        _liste_ref_descr[cpt_descr] = liste_ref_descr_fam[cpt_fam][idescr];

        nbr_ref_descr++;

        /* On stocke le numéro de la famille à laquelle il appartient */

        _tab_fam_descr[cpt_descr].nbr = 1;
        ECS_MALLOC(_tab_fam_descr[cpt_descr].val,
                   _tab_fam_descr[cpt_descr].nbr, ecs_int_t);

        _tab_fam_descr[cpt_descr].val[0] = ptr_fam->num;

      }
      else {

        /* Le descripteur a déjà été rencontre */
        /* On rajoute le numéro de famille */

        ECS_REALLOC(_tab_fam_descr[cpt_descr].val,
                    _tab_fam_descr[cpt_descr].nbr + 1, ecs_int_t);

        _tab_fam_descr[cpt_descr].val[_tab_fam_descr[cpt_descr].nbr]
          = ptr_fam->num;

        _tab_fam_descr[cpt_descr].nbr++;

      }

    } /* Fin : boucle sur les descripteurs de la famille */

    cpt_fam++;

  } /* Fin : boucle sur les familles */

  for (ifam = 0; ifam < nbr_fam; ifam++)
    ECS_FREE(liste_ref_descr_fam[ifam]);
  ECS_FREE(liste_ref_descr_fam);
  ECS_FREE(liste_nbr_descr_fam);

  /* Séparation des descripteurs de type "couleur" */
  /*         et des descripteurs de type "groupe"  */
  /*-----------------------------------------------*/

  tab_couleur.nbr = 0;
  tab_groupe.nbr  = 0;

  ECS_MALLOC(tab_couleur.val,      nbr_ref_descr, ecs_int_t);
  ECS_MALLOC(liste_num_descr_coul, nbr_ref_descr, ecs_int_t);
  ECS_MALLOC(tab_groupe.val,       nbr_ref_descr, ecs_int_t);
  ECS_MALLOC(liste_num_descr_grp,  nbr_ref_descr, ecs_int_t);

  for (idescr = 0; idescr < nbr_ref_descr; idescr++) {

    descr_typ = ecs_descr__ret_typ(_liste_ref_descr[idescr]);

    switch (descr_typ) {

    case ECS_DESCR_COULEUR:

      tab_couleur.val[tab_couleur.nbr]
        = ecs_descr__ret_ide(_liste_ref_descr[idescr]);
      liste_num_descr_coul[tab_couleur.nbr++] = idescr;

      break;

    case ECS_DESCR_GROUPE:

      tab_groupe.val[tab_groupe.nbr]
        = ecs_descr__ret_ide(_liste_ref_descr[idescr]);
      liste_num_descr_grp[tab_groupe.nbr++] = idescr;

      break;

    default:

      assert(descr_typ == ECS_DESCR_COULEUR ||
             descr_typ == ECS_DESCR_GROUPE    );

    }
  }

  if (tab_couleur.nbr != 0)
    ECS_REALLOC(tab_couleur.val, tab_couleur.nbr, ecs_int_t);
  else
    ECS_FREE(tab_couleur.val);

  if (tab_groupe.nbr != 0)
    ECS_REALLOC(tab_groupe.val, tab_groupe.nbr, ecs_int_t);
  else
    ECS_FREE(tab_groupe.val);

  /* Ordination des couleurs et groupes par leur numéro */
  /*----------------------------------------------------*/

  assert(nbr_ref_descr = tab_couleur.nbr + tab_groupe.nbr);

  _tab_descr.nbr = nbr_ref_descr;
  ECS_MALLOC(_tab_descr.val, _tab_descr.nbr, ecs_int_t);

  for (idescr = 0; idescr < (int)(_tab_descr.nbr); idescr++)
    _tab_descr.val[idescr] = -1;

  idescr = 0;

  if (tab_couleur.nbr != 0) {

    tab_renum_coul.nbr = tab_couleur.nbr;
    ECS_MALLOC(tab_renum_coul.val, tab_renum_coul.nbr, ecs_int_t);
    for (icoul = 0; icoul < (int)(tab_renum_coul.nbr); icoul++)
      tab_renum_coul.val[icoul] = icoul;

    ecs_tab_int__trie(tab_couleur,
                      tab_renum_coul);

    ECS_FREE(tab_couleur.val);

    for (icoul = 0; icoul < (int)(tab_renum_coul.nbr); icoul++)
      _tab_descr.val[idescr++]
        = liste_num_descr_coul[tab_renum_coul.val[icoul]];

    ECS_FREE(tab_renum_coul.val);

  }

  if (tab_groupe.nbr != 0) {

    tab_renum_grp.nbr = tab_groupe.nbr;
    ECS_MALLOC(tab_renum_grp.val, tab_renum_grp.nbr, ecs_int_t);
    for (igrp = 0; igrp < (int)(tab_renum_grp.nbr); igrp++)
      tab_renum_grp.val[igrp] = igrp;

    ecs_tab_int__trie(tab_groupe,
                      tab_renum_grp);

    ECS_FREE(tab_groupe.val);

    for (igrp = 0; igrp < (int)(tab_renum_grp.nbr); igrp++)
      _tab_descr.val[idescr++]
        = liste_num_descr_grp[tab_renum_grp.val[igrp]];

    ECS_FREE(tab_renum_grp.val);

  }

  ECS_FREE(liste_num_descr_coul);
  ECS_FREE(liste_num_descr_grp);

  assert(idescr = _tab_descr.nbr);

  /* Valeurs de retour */
  /*-------------------*/

  *tab_descr       = _tab_descr;
  *tab_fam_descr   = _tab_fam_descr;
  *liste_ref_descr = _liste_ref_descr;
}

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
 *   une liste chaînée de champs
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
 *  Fonction qui affiche par attribut (couleur ou groupe),
 *   les numéros des familles auxquelles l'attribut appartient
 *   à partir d'une liste chaînée de familles dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__aff_fam_att(ecs_famille_t  *fam_tete)
{
  size_t  idescr;
  size_t  ifam;

  ecs_tab_int_t     tab_descr;
  ecs_tab_int_t    *tab_fam_descr;
  ecs_descr_t     **liste_ref_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_loc_famille_chaine__fam_att(fam_tete,
                                  NULL,
                                  &tab_descr,
                                  &tab_fam_descr,
                                  &liste_ref_descr);

  /* Affichage des attributs triés en fonction des */
  /* familles auxquelles ils appartiennent         */
  /*-----------------------------------------------*/

  for (idescr = 0; idescr < tab_descr.nbr; idescr++) {

    ecs_descr__affiche(liste_ref_descr[idescr], 0);

    for (ifam = 0; ifam < tab_fam_descr[idescr].nbr; ifam++) {

      printf("  %*s%s %d\n",
             (int)(strlen(_("Color")) + 1), "", _("Family"),
             (int)(tab_fam_descr[idescr].val[ifam]));

    }
  }

  /* Libération mémoire */

  ECS_FREE(liste_ref_descr);

  for (idescr = 0; idescr < tab_descr.nbr; idescr++)
    ECS_FREE(tab_fam_descr[idescr].val);
  ECS_FREE(tab_fam_descr);

  ECS_FREE(tab_descr.val);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie pour chaque numéro de famille
 *   le nombre et la liste des identificateurs de type couleur
 *   des descripteurs de la famille
 *----------------------------------------------------------------------------*/

ecs_tab_int_t *
ecs_famille_chaine__ret_ide(ecs_famille_t   *fam_tete)
{
  int  ifam;
  int  nbr_fam;

  ecs_famille_t  *ptr_fam;

  ecs_tab_int_t  *tab_propr_fam;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(fam_tete    != NULL);

  /* Détermination du nombre de familles de la liste chaînée */
  /*---------------------------------------------------------*/

  nbr_fam = ecs_famille_chaine__ret_nbr(fam_tete);

  /* Détermination des couleurs pour chaque famille de la liste chaînée */
  /*--------------------------------------------------------------------*/

  ECS_MALLOC(tab_propr_fam, nbr_fam, ecs_tab_int_t);

  for (ptr_fam  = fam_tete, ifam = 0;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui, ifam++) {

    tab_propr_fam[ifam] = ecs_famille__ret_ide(ptr_fam);

  }

  return tab_propr_fam;
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
 *  Fonction qui construit une liste chaînée de descripteurs de type
 *   "couleur" et "groupe"
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

ecs_tab_bool_t
ecs_famille_chaine__indic_fam_att(const ecs_famille_t  *fam_tete)
{
  size_t    ifam;
  int       num_fam_max;

  const ecs_famille_t  *ptr_fam;

  ecs_tab_bool_t   tab_select;

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

  tab_select.nbr = num_fam_max + 1;
  ECS_MALLOC(tab_select.val, tab_select.nbr, bool      );

  for (ifam = 0; ifam < tab_select.nbr; ifam++)
    tab_select.val[ifam] = false;

  /* Marquage des familles correspondant aux critères de sélection */
  /*---------------------------------------------------------------*/

  for (ptr_fam  = fam_tete;
       ptr_fam != NULL;
       ptr_fam  = ptr_fam->l_famille_sui) {

    tab_select.val[ptr_fam->num] = true;
  }

  return tab_select;
}

/*----------------------------------------------------------------------------*/


