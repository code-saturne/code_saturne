#ifndef _ECS_POST_ENS_PRIV_H_
#define _ECS_POST_ENS_PRIV_H_

/*============================================================================
 *  Définition de la structure `ecs_post_ens_t' pour post-traitement EnSight
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------*
 *  Fichiers `include' librairie standard C ou BFT
 *----------------------------------------------------------------------------*/

#include <ecs_file.h>


/*----------------------------------------------------------------------------*
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------*
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens.h"


/*============================================================================
 *  Définition de macros
 *============================================================================*/


/*============================================================================
 *  Définitions d'énumerations
 *============================================================================*/


/*============================================================================
 *  Définitions de types
 *============================================================================*/

/* Structure conservant des informations sur les parts EnSight */

typedef struct {

  char        *nom_part;       /* Nom du "part" EnSight */
  int          num_part;       /* Numéro du "part" */
  ecs_int_t    nbr_som;        /* Nombre de sommets associés */
  int          nbr_typ_ele;    /* Nombre de type d'éléments considérés */
  ecs_int_t   *nbr_ele_typ;    /* Nombre d'éléments par type */
  char       **nom_typ_ele;    /* Noms des types d'éléments */
  ecs_int_t   *lst_parents;    /* Liste des éléments parents */

} ecs_post_ens_part_t ;


/* Structure définissant un cas EnSight */

struct _ecs_post_ens_t {

  char                  *nom_cas;        /* Nom du cas */
  char                  *prefixe_rep;    /* Préfixe du répertoire EnSight  */
  char                  *prefixe_fic;    /* Préfixe des fichiers EnSight */
  char                  *nom_fic_case;   /* Nom du fichier "case" */

  int                    nbr_part;       /* Nombre de ``parts'' géométrie */
  ecs_post_ens_part_t  **tab_part;       /* Descripteurs des ``parts'' */

  ecs_file_t            *fic_geo;        /* Pointeur sur fichier géométrie */

  bool                   modifie;        /* Modification depuis dernière
                                            écriture du fichier case */

};

/*----------------------------------------------------------------------------*/

#endif /* _ECS_POST_ENS_PRIV_H_ */
