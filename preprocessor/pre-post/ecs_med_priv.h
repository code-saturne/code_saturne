#ifndef _ECS_MED_PRIV_H_
#define _ECS_MED_PRIV_H_

/*============================================================================
 *  Définition de la structure `_ecs_med_t' pour les entrées ou sorties
 *   au format MED
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

#include "cs_config.h"

#if defined(HAVE_MED)


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "MED"
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MED)

#ifdef __cplusplus
extern "C" {
#endif

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include <med.h>

#ifdef __cplusplus
}
#endif

#endif /* HAVE_MED */


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_med.h"


/*============================================================================
 *                       Définition de macros
 *============================================================================*/

/* Definition des éléments */
/*=========================*/

/* Tableau donnant la liste des éléments `paraboliques' ou `cubiques'
   qui sont transformés en leur equivalent `lineaire' */

#define ECS_MED_ORDER_LINEAR                               1
#define ECS_MED_ORDER_PARABOLIC                            2

#define ECS_MED_NBR_TYP_ELT                               14
#define ECS_MED_NBR_MAX_SOM                                8

/*============================================================================
 *                         Définitions d'énumerations
 *============================================================================*/


/*============================================================================
 *                          Définitions de types
 *============================================================================*/

/* Structure d'information sur les maillages pour le post traitement */

typedef struct {

  char                   *nom_maillage;      /* Nom du maillage MED */
  char                    nom_maillage_med[MED_TAILLE_NOM + 1];  /* Nom MED */

  ecs_int_t               dim_entite;        /* Dimension entité max. */
  ecs_int_t               nbr_typ_ele;       /* Nombre de types d'éléments */
  ecs_int_t              *nbr_ele_typ;       /* Nombre d'éléments par type */
  med_geometrie_element  *med_typ;           /* Types MED des élements */

} ecs_med_maillage_t;

/* Structure définissant un cas MED */

struct _ecs_med_t {

  char                  *nom_cas;         /* Nom du cas */
  char                  *nom_fic;         /* Nom du fichier MED */
  med_idt                fid;             /* Identificateur de fichier MED */

  ecs_int_t              nbr_var;         /* Nombre de variables */
  char                 **nom_var;         /* Noms des variables */

  ecs_int_t              nbr_maillages;   /* Nombre de maillages */
  ecs_med_maillage_t   **tab_maillages;   /* Descripteurs des maillages */

  bool                   no_poly;         /* Ne pas écrire les polygones
                                             ou polyèdres */

};


typedef struct {

  med_geometrie_element  med_typ;     /* Type MED de l'element */
  ecs_elt_typ_t          ecs_typ;     /* Type ECS de l'element */
  ecs_int_t              order;       /* Ordre    de l'element */
                                      /* Liste des numeros de sommet ECS */
  ecs_int_t              num_som[ECS_MED_NBR_MAX_SOM];

} ecs_fic_med_init_elt_t;

/*============================================================================
 *                Définitions de variables globales statiques
 *============================================================================*/

extern const ecs_fic_med_init_elt_t
ecs_fic_med_init_elt_liste_c[ECS_MED_NBR_TYP_ELT];


#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_MED_PRIV_H_ */

