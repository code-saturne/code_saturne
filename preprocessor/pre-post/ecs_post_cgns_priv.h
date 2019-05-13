#ifndef _ECS_POST_CGNS_PRIV_H_
#define _ECS_POST_CGNS_PRIV_H_

/*============================================================================
 *  Définition de la structure `_ecs_cgns_t' pour les sorties au format CGNS
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_config.h"

#if defined(HAVE_CGNS)

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"

/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_cgns.h"

/*============================================================================
 *                       Définition de macros
 *============================================================================*/

#define ECS_CGNS_TAILLE_NOM    32      /* Longueur nom max (la documentation
                                          ne précise rien, mais les exemples
                                          CGNS et le fichier ADF.h prévoient
                                          des chaînes de 32 caractères) */

/* Compatibilité avec diverses versions de CGNS */

#if !defined(CG_OK)
#define CG_OK  ALL_OK
#endif

/*============================================================================
 *                          Définitions de types
 *============================================================================*/


/* Structure d'information sur les maillages pour le post traitement */

typedef struct {

  char                *nom_maillage;   /* Nom du maillage CGNS */

  char                *nom_fic;        /* Nom du fichier CGNS */
  int                  num_fic;        /* Identificateur de fichier CGNS */

  /* Informations sur état du cas en post traitement */

  int                  dim_espace;     /* Dimension de l'espace associé */
  int                  dim_entite;     /* Dimension entités maillage */

  bool                 fic_ouvert;     /* Fichier ouvert ou non ? */

} ecs_post_cgns_base_t;


/* Structure définissant un cas de sortie CGNS */

struct _ecs_post_cgns_t {

  /* Informations principales (post traitement) */

  char                   *nom_cas;        /* Nom du cas */

  ecs_int_t               nbr_bases;      /* Nombre de ``bases'' */
  ecs_post_cgns_base_t  **tab_bases;      /* Descripteurs des ``bases'' */

};

#endif /* HAVE_CGNS */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_POST_CGNS_PRIV_H_ */

