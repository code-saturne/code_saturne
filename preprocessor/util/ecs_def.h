#ifndef _ECS_DEF_H_
#define _ECS_DEF_H_

/*============================================================================
 * Définitions, variables globales, et fonctions de base
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


/*----------------------------------------------------------------------------
 *  Fichier  `include' associé à la configuration du cas
 *  (généré par le Makefile)
 *----------------------------------------------------------------------------*/

#include "cs_config.h"


/*============================================================================
 * Définition de l'architecture
 *============================================================================*/

#if defined(__sgi__) || defined(__sgi) || defined(sgi)
#define ECS_ARCH_IRIX_64

#elif defined(__hpux__) || defined(__hpux) || defined(hpux)
#define ECS_ARCH_HP_UX

#elif defined(__linux__) || defined(__linux) || defined(linux)
#define ECS_ARCH_Linux

#elif defined(__sun__) || defined(__sun) || defined(sun)
#define ECS_ARCH_SunOS

#elif defined(__uxpv__) || defined(__uxpv) || defined(uxpv)
#define ECS_ARCH_UNIX_System_V

#elif defined(__osf__)
#define ECS_ARCH_OSF1

#endif


/*============================================================================
 * Définitions de type C99 qui ne sont pas toujurs fournies par des
 * compilateurs ou environnements plus anciens.
 *============================================================================*/

/*
 * En général, stdint.h est inclus par inttypes.h, mais seulement inttypes.h
 * existe sur certains systèmes, tels que Tru64Unix.
 */

#if HAVE_STDINT_H
# include <stdint.h>
#elif HAVE_INTTYPES_H
# include <inttypes.h>
#endif

/* _Bool */

#if HAVE_STDBOOL_H
# include <stdbool.h>
#else
# if !HAVE__BOOL
#  ifdef __cplusplus
typedef bool _Bool;
#  else
typedef unsigned char _Bool;
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
# define __bool_true_false_are_defined 1
#endif

/* int32_t */

#if !defined(HAVE_INT32_T)
# if (SIZEOF_INT == 4)
typedef int int32_t;
# elif (SIZEOF_SHORT == 4)
typedef short int32_t;
# else
#  error
# endif
#endif

/* int64_t */

#if !defined(HAVE_INT64_T)
# if (SIZEOF_INT == 8)
typedef int int64_t;
# elif (SIZEOF_LONG == 8)
typedef long int64_t;
# elif (HAVE_LONG_LONG == 8)
typedef long long int64_t;
# else
#  error
# endif
#endif

/* uint32_t */

#if !defined(HAVE_UINT32_T)
# if (SIZEOF_INT == 4)
typedef unsigned uint32_t;
# elif (SIZEOF_SHORT == 4)
typedef unsigned short uint32_t;
# else
#  error
# endif
#endif

/* uint64_t */

#if !defined(HAVE_UINT64_T)
# if (SIZEOF_INT == 8)
typedef unsigned uint64_t;
# elif (SIZEOF_LONG == 8)
typedef unsigned long uint64_t;
# elif (HAVE_LONG_LONG)
typedef unsigned long long uint64_t;
# else
#  error
# endif
#endif


/*============================================================================
 * Définitions de types
 *============================================================================*/

/*  Types de dimension fixée. */

#if (defined ECS_ARCH_Linux)
#include <stdint.h>

typedef int32_t         ecs_int_32_t;   /* Entier sur 4 octets */

#else

typedef int             ecs_int_32_t;   /* Entier sur 4 octets */

#endif


/* Types usuels */

#if defined(USE_LONG_INT)

typedef long            ecs_int_t;      /* Entier */
typedef unsigned long   ecs_size_t;     /* Taille pour les index */

#else

typedef int             ecs_int_t;      /* Entier */
typedef unsigned        ecs_size_t;     /* Taille pour les index */

#endif

typedef double          ecs_coord_t;    /* Réel (virgule flottante) */
typedef char            ecs_byte_t;     /* Octet (unité de mémoire non typée) */

/* Énumération de type ("type de type") pour transmettre le type d'une donnée */

typedef enum {
  ECS_TYPE_char,
  ECS_TYPE_bool,
  ECS_TYPE_ecs_int_t,
  ECS_TYPE_ecs_int_32_t,
  ECS_TYPE_ecs_coord_t,
  ECS_TYPE_ecs_size_t,
  ECS_TYPE_size_t,
  ECS_TYPE_void
} ecs_type_t;

/* Énumération liée au type d'élément */

typedef enum {

  ECS_ELT_TYP_NUL,         /*  Pas de type */

  ECS_ELT_TYP_FAC_TRIA,    /*  Triangle */
  ECS_ELT_TYP_FAC_QUAD,    /*  Quadrangle */

  ECS_ELT_TYP_CEL_TETRA,   /*  Tétraèdre */
  ECS_ELT_TYP_CEL_PYRAM,   /*  Pyramide */
  ECS_ELT_TYP_CEL_PRISM,   /*  Prisme */
  ECS_ELT_TYP_CEL_HEXA,    /*  Hexaedre */

  ECS_ELT_TYP_FAC_POLY,    /*  Polygone quelconque */
  ECS_ELT_TYP_CEL_POLY,    /*  Polyedre quelconque */

  ECS_ELT_TYP_FIN

} ecs_elt_typ_t ;


/* Qualificateurs de types restrict (existe en standard en C99) */

#if defined(__GNUC__)
#define restrict __restrict
#else
#define restrict
#endif


/*=============================================================================
 * Définitions de macros
 *============================================================================*/

/* Macros "classiques" */

#define ECS_ABS(a)     ((a) <  0  ? -(a) : (a))  /* Valeur absolue de a */
#define ECS_MIN(a,b)   ((a) > (b) ?  (b) : (a))  /* Minimum de a et b */
#define ECS_MAX(a,b)   ((a) < (b) ?  (b) : (a))  /* Maximum de a et b */


/* Définition du caractère de séparation de répertoires */

#define ECS_PATH_SEP             '/'

#define ECS_REAL_PRECISION    1.e-13

#define ECS_STR_SIZE       80

#define ECS_PAS_NUL         0
#define ECS_PAS_UNITE       1

#define ECS_LNG_AFF_STR           43
#define ECS_LNG_AFF_ENT            8
#define ECS_LNG_AFF_REE_MANTIS    11
#define ECS_LNG_AFF_REE_PRECIS     2

#define ECS_FMT_AFF_REE_PARAM     "%.15E"

/*
 * Macros pour internationalisation éventuelle via gettext() ou une fonction
 * semblable (pour encadrer les chaînes de caractères imprimables)
 */

#if defined(ENABLE_NLS)

#include <libintl.h>
#define _(String) gettext(String)
#define gettext_noop(String) String
#define N_(String) gettext_noop(String)

#else

#define _(String) String
#define N_(String) String
#define textdomain(Domain)
#define bindtextdomain(Package, Directory)

#endif

/*=============================================================================
 * Définitions de variables globales
 *============================================================================*/

extern char      ecs_glob_build_date[]; /* Date de compilation */

extern int       ecs_glob_have_cgns;    /* Support du format CGNS */
extern int       ecs_glob_cgns_ver_maj;
extern int       ecs_glob_cgns_ver_min;
extern int       ecs_glob_cgns_ver_rel;

extern int       ecs_glob_have_med;     /* Support de la librairie d'échange de
                                           maillages et champs MED */

/* Détermination du type d'élément en connectivité nodale en fonction du
   nombre de ses sommets, pour des éléments de dimension 2 (faces)
   ou 3 (cellules) ;
   Au dessus de 8 sommets, on a toujours le type ECS_ELT_TYP_FAC_POLY pour
   les faces et ECS_ELT_TYP_FAC_POLY pour les cellules */

extern  const ecs_elt_typ_t  ecs_glob_typ_elt[2][9];

/*=============================================================================
 * Prototypes de fonctions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction d'initialisation de la gestion des erreurs
 *----------------------------------------------------------------------------*/

void
ecs_init_gestion_erreur(void);

/*----------------------------------------------------------------------------
 * Fonction d'arret
 *----------------------------------------------------------------------------*/

void
ecs_exit(int  statut);

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un avertissement
 *----------------------------------------------------------------------------*/

void
ecs_warn(void);

/*----------------------------------------------------------------------------
 * Fonction d'arrêt sur erreur
 *----------------------------------------------------------------------------*/

void
ecs_error(const char  *file_name,
          const int    line_num,
          const int    sys_error_code,
          const char  *format,
          ...);

/*----------------------------------------------------------------------------
 *  Fonction qui imprime une chaîne de caractères avec une largeur
 *  de colonne donnée.
 *----------------------------------------------------------------------------*/

void
ecs_print_padded_str(const char  *str,
                     int          width);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_DEF_H_ */
