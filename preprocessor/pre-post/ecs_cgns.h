#ifndef _ECS_CGNS_PUBL_H_
#define _ECS_CGNS_PUBL_H_

/*============================================================================
 *  Définitions pour le format CGNS
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

#if defined(HAVE_CGNS)

/*----------------------------------------------------------------------------
 * Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 * Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "CGNS"
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#include <cgnslib.h>

#ifdef __cplusplus
}
#endif

/*----------------------------------------------------------------------------
 * Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 * Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 * Définition de macros
 *============================================================================*/

/* Tableau donnant la liste des éléments              */
/* `paraboliques' ou `cubiques'                       */
/* qui sont transformés en leur equivalent `linéaire' */
/* -------------------------------------------------- */

#define ECS_CGNS_NBR_TYP_ELT                               20
#define ECS_CGNS_NBR_MAX_SOM                                8

#if defined(CGNS_SCOPE_ENUMS)
#define CS_CG_ENUM( e ) CG_ ## e
#else
#define CS_CG_ENUM( e ) e
#endif

/*============================================================================
 * Définitions de types
 *============================================================================*/

/* Structure définissant un cas CGNS */

typedef struct _ecs_cgns_t ecs_cgns_t;


/* Structure définissant un type d'élément CGNS */

typedef struct {

#if defined(CGNS_SCOPE_ENUMS)
  CG_ElementType_t  cgns_type;   /* Type CGNS de l'élément */
#else
  ElementType_t     cgns_type;   /* Type CGNS de l'élément */
#endif
  ecs_elt_typ_t     ecs_type;    /* Type ECS  de l'élément */
  ecs_int_t         nbr_som ;    /* Nombre de sommets associés */
  ecs_int_t         num_som[ECS_CGNS_NBR_MAX_SOM] ; /* Sommets ECS */

} ecs_cgns_elt_t;


/*============================================================================
 * Définitions de variables globales statiques
 *============================================================================*/

/* Le tableau suivant est donné dans le même ordre que les définitions
   de cgnslib.h (mais ne contient pas les deux premières entrées) */

extern const ecs_cgns_elt_t
ecs_cgns_elt_liste_c[ECS_CGNS_NBR_TYP_ELT];

#endif /* HAVE_CGNS */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CGNS_PUBL_H_ */

