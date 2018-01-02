#ifndef _ECS_FAMILLE_CHAINE_H_
#define _ECS_FAMILLE_CHAINE_H_

/*============================================================================
 *  Prototypes des fonctions de base
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

#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
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
                         ecs_descr_t       *descr_tete);

/*----------------------------------------------------------------------------
 *  Fonction libérant la portion d'une liste chaînée de familles
 *   à partir d'un noeud dont le pointeur est donné en argument.
 *  Le noeud est à NULL au retour de la fonction
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__detruit(ecs_famille_t  **this_fam_noeud);

/*----------------------------------------------------------------------------
 *  Fonction imprimant à partir d'un noeud `ecs_famille_t' donné
 *   une liste chaînée de tables
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__imprime(const ecs_famille_t  *this_fam_noeud,
                            ecs_int_t             imp_col,
                            FILE                 *fic_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets
 *   d'une chaîne de structures `ecs_famille_t'
 *----------------------------------------------------------------------------*/

float
ecs_famille_chaine__ret_taille(const ecs_famille_t  *this_fam_noeud);

/*----------------------------------------------------------------------------
 *  Fonction qui ajoute à la fin d'une liste chaînée de familles
 *   réceptrice dont la tête est donnée,
 *   une liste chaînée de familles à concaténer dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__ajoute(ecs_famille_t  **this_fam_tete,
                           ecs_famille_t   *fam_concat_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui affiche la définition de la famille de numéro donné
 *   à partir de la liste chaînée des familles dont la tête est donnée
 *----------------------------------------------------------------------------*/

void
ecs_famille_chaine__affiche(const ecs_int_t   num_fam,
                            ecs_famille_t    *fam_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre de familles
 *   de la liste chaînée des familles dont la tête est donnée
 *----------------------------------------------------------------------------*/

int
ecs_famille_chaine__ret_nbr(const ecs_famille_t  *fam_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui copie une liste chaînée de familles
 *   dont la tête est donnée
 *----------------------------------------------------------------------------*/

ecs_famille_t  *
ecs_famille_chaine__copie(ecs_famille_t  *famille_tete);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie pour chaque numéro de famille
 *   le nombre et une liste de pointeurs sur les noms des groupes
 *   des descripteurs de la famille
 *----------------------------------------------------------------------------*/

ecs_tab_char_t *
ecs_famille_chaine__ret_nom(ecs_famille_t   *fam_tete);

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
                               int             *nbr_max_att_fam);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un tableau marquant chaque numéro de famille.
 *
 *  La libération du tableau est à la charge du code appelant
 *----------------------------------------------------------------------------*/

bool *
ecs_famille_chaine__indic_fam_att(const ecs_famille_t  *fam_tete);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_FAMILLE_CHAINE_H_ */
