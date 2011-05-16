#ifndef _ECS_MED_PRIV_H_
#define _ECS_MED_PRIV_H_

/*============================================================================
 *  Définition de la structure `_ecs_med_t' pour les entrées ou sorties
 *   au format MED
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2011 EDF S.A., France

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

#undef HAVE_MPI /* For MED 2.9 */

#include <med.h>

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

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

/* Med version */

#if !defined(MED_NUM_MAJEUR)
#define MED_NUM_MAJEUR 2
#define MED_NUM_MINEUR 3
#endif

#if !defined(MED_MAJOR_NUM)
#define MED_MAJOR_NUM MED_NUM_MAJEUR
#define MED_MINOR_NUM MED_NUM_MINEUR
#endif

#if MED_MAJOR_NUM == 2 && MED_MINOR_NUM < 9
#define ECS_MED_VERSION 2
#else
#define ECS_MED_VERSION 3
#endif

/* Map MED 2 to MED3 names */

#if ECS_MED_VERSION == 2

#define MED_NAME_SIZE MED_TAILLE_NOM
#define MED_SNAME_SIZE MED_TAILLE_PNOM
#define MED_LNAME_SIZE MED_TAILLE_LNOM
#define MED_COMMENT_SIZE MED_TAILLE_DESC

#define MED_CELL MED_MAILLE
#define MED_POLYGON MED_POLYGONE
#define MED_POLYHEDRON MED_POLYEDRE
#define MED_DESCENDING_EDGE MED_ARETE
#define MED_DESCENDING_FACE MED_FACE

#define med_axis_type med_repere
#define med_bool med_booleen
#define med_mesh_type med_maillage
#define med_geometry_type med_geometrie_element
#define med_entity_type med_entite_maillage
#define med_field_type med_type_champ
#endif

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

  char               *nom_maillage;      /* Nom du maillage MED */
  char                nom_maillage_med[MED_NAME_SIZE + 1];  /* Nom MED */

} ecs_med_maillage_t;

/* Structure définissant un cas MED */

struct _ecs_med_t {

  char                  *nom_cas;         /* Nom du cas */
  char                  *nom_fic;         /* Nom du fichier MED */

  med_idt                fid;             /* Identificateur de fichier MED */
  med_int                version[3];      /* MED version used to write file */

  ecs_int_t              nbr_maillages;   /* Nombre de maillages */
  ecs_med_maillage_t   **tab_maillages;   /* Descripteurs des maillages */

};


typedef struct {

  med_geometry_type   med_typ;     /* Type MED de l'element */
  ecs_elt_typ_t       ecs_typ;     /* Type ECS de l'element */
  ecs_int_t           order;       /* Ordre    de l'element */
                                   /* Liste des numeros de sommet ECS */
  ecs_int_t           num_som[ECS_MED_NBR_MAX_SOM];

} ecs_fic_med_init_elt_t;

/*============================================================================
 *                Définitions de variables globales statiques
 *============================================================================*/

extern const ecs_fic_med_init_elt_t
ecs_fic_med_init_elt_liste_c[ECS_MED_NBR_TYP_ELT];


#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_MED_PRIV_H_ */

