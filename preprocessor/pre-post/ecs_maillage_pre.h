#ifndef _ECS_MAILLAGE_PRE_H_
#define _ECS_MAILLAGE_PRE_H_

/*============================================================================
 *  Definition des fonctions de base
 *   de remplissage de la structure de maillage a partir des donnees lues
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
 *                                 Visibilite
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Assign vertex coordinates to a mesh structure.
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__cree_som(ecs_maillage_t  *maillage,
                           size_t           nbr_som,
                           ecs_coord_t     *som_val_coord);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un vecteur de structures "entité de maillage"
 *   pour les éléments (aretes, faces et cellules),
 *   construites à partir de tableaux remplis lors de la lecture des données
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__cree_elt(ecs_maillage_t   *maillage,
                           const size_t      nbr_elt_ent[],
                           ecs_size_t       *elt_pos_som_ent[],
                           ecs_int_t        *elt_val_som_ent[],
                           int              *elt_val_famille_ent[],
                           ecs_int_t        *elt_val_couleur_ent[],
                           const ecs_int_t   nbr_coul_ent[],
                           ecs_int_t        *val_coul_ent[],
                           ecs_size_t       *nbr_elt_coul_ent[]);

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation de références à des labels en indices
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__label_en_indice(size_t      nbr_label,
                                  size_t      nbr_ref,
                                  ecs_int_t  *val_label,
                                  ecs_int_t  *val_ref);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie l'identificateur de l'entité
 *   auquel appartient le type géométrique donné
 *--------------------------------------------------------------------------- */

ecs_entmail_t
ecs_maillage_pre__ret_typ_geo(int  typ_geo);

/*----------------------------------------------------------------------------
 * Fonction qui affiche le nombre d'éléments par entité
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__aff_nbr_par_ent(size_t        nbr_som,
                                  const size_t  nbr_elt_ent[],
                                  int           lng_imp);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_MAILLAGE_PRE_H_ */
