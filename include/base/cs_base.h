/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
*
*     contact: saturne-support@edf.fr
*
*     The Code_Saturne Kernel is free software; you can redistribute it
*     and/or modify it under the terms of the GNU General Public License
*     as published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*
*     The Code_Saturne Kernel is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
*     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public License
*     along with the Code_Saturne Kernel; if not, write to the
*     Free Software Foundation, Inc.,
*     51 Franklin St, Fifth Floor,
*     Boston, MA  02110-1301  USA
*
*============================================================================*/

#ifndef __CS_BASE_H__
#define __CS_BASE_H__

/*============================================================================
 * Définitions, variables globales, et fonctions de base
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#if defined(_CS_HAVE_MPI)

#include <mpi.h>

#if defined(_CS_HAVE_MPE)
#include <mpe.h>
#endif

#endif


/*=============================================================================
 * Définitions de macros
 *============================================================================*/

/* Nom du système */

#if defined(__sgi__) || defined(__sgi) || defined(sgi)
#define _CS_ARCH_IRIX_64

#elif defined(__hpux__) || defined(__hpux) || defined(hpux)
#define _CS_ARCH_HP_UX

#elif defined(__linux__) || defined(__linux) || defined(linux)
#define _CS_ARCH_Linux

#elif defined(__sun__) || defined(__sun) || defined(sun)
#define _CS_ARCH_SunOS

#elif defined(__uxpv__) || defined(__uxpv) || defined(uxpv)
#define _CS_ARCH_UNIX_System_V

#endif

/*
 * Macro utile pour gérér les differences de noms de symboles (underscore ou
 * non, miniscules ou majuscules entre C et FORTRAN) pour l'édition de liens
 */

#if !defined (__hpux)
#define CS_PROCF(x, y) x##_
#else
#define CS_PROCF(x, y) x
#endif

/*
 * Macro utile pour gérér arguments 'longueur de chaîne de caractères Fortran,
 * inutilisés lors des appel mais placés par de nombreux compilateurs)
 * Le compilateur Fujitsu VPP 5000 ne supporte pas ces listes  de longueur
 * variables dans les appels entre C et FORTRAN (mais ca marche pour les
 * appels C-C et FORTRAN-FORTRAN)
 */

#if defined (__uxpv__)  /* Cas Fujitsu VPP 5000 */
#define CS_ARGF_SUPP_CHAINE
#else
#define CS_ARGF_SUPP_CHAINE , ...
#endif

/* Sur certaines machines tells que IBM Blue Gene/L, certaines opérations
 * peuvent être mieux optimisées sur des données respectant un certain
 * alignement en mémoire; (si 0, aucun alignement exploité) */

#if defined(__blrts__)
#define CS_MEM_ALIGN 16
#else
#define CS_MEM_ALIGN 0
#endif

#define CS_DIM_3              3                 /* Dimension de l'espace */

/* Macros "classiques" */

#define CS_ABS(a)     ((a) <  0  ? -(a) : (a))  /* Valeur absolue de a */
#define CS_MIN(a,b)   ((a) > (b) ?  (b) : (a))  /* Minimum de a et b */
#define CS_MAX(a,b)   ((a) < (b) ?  (b) : (a))  /* Maximum de a et b */

/*
 * Macros pour internationalisation éventuelle via gettext() ou une fonction
 * semblable (pour encadrer les chaînes de caractères imprimables)
 */

#undef _
#define _(String) String

#undef N_
#define N_(String) String

/* Définition de la version du langage C utilisé (C89 ou C99) */

#if defined(__STDC_VERSION__)
#  define _CS_STDC_VERSION __STDC_VERSION__
#else
#  define _CS_STDC_VERSION 1989
#endif

/*
 * Redéfinition des commandes "inline" et "restrict" incompatible avec
 * certains compilateurs C89 (standard en C99)
 */

#if (_CS_STDC_VERSION < 199901L)

#  if defined(__GNUC__)
#    define inline __inline
#    define restrict __restrict
#  else
#    define inline
#    define restrict
#  endif

#else

/* Même en C99, le compilateur IRIX64 (ancien) ne semble pas accepter
 * inline (à vérifier) */

#  if defined(_CS_ARCH_IRIX_64) && !defined(__GNUC__)
#    define inline
#    define restrict
#  endif

#endif


/*============================================================================
 * Définitions de types
 *============================================================================*/

typedef int              cs_int_t;      /* Entier */
typedef double           cs_real_t;     /* Réel (virgule flottante) */
typedef char             cs_byte_t;     /* Octet (unité de mémoire non typée) */

typedef cs_real_t        cs_point_t[3];

typedef enum {                          /* Booléen */
  CS_FALSE ,
  CS_TRUE
} cs_bool_t;

#if !defined(false)
#define false CS_FALSE
#endif

#if !defined(true)
#define true CS_TRUE
#endif

/* Définitions pour opérations collectives (min, max, somme) sous MPI */

#if defined(_CS_HAVE_MPI)

#define CS_MPI_INT       MPI_INT         /* Si cs_real_t est un double ;
                                            sinon redéfinir en MPI_xxx */
#define CS_MPI_REAL      MPI_DOUBLE      /* Si cs_real_t est un double ;
                                            sinon redéfinir en MPI_REAL */
#define CS_MPI_REAL_INT  MPI_DOUBLE_INT  /* Si cs_real_t est un double ;
                                            sinon redéfinir en MPI_REAL_INT */

typedef struct
{
  cs_real_t val;
  cs_int_t  rang;
} cs_mpi_real_int_t;

#endif /* defined(_CS_HAVE_MPI) */

/* Énumération de type ("type de type") pour transmettre le type d'une donnée */

typedef enum {
  CS_TYPE_char,
  CS_TYPE_cs_int_t,
  CS_TYPE_cs_real_t,
  CS_TYPE_cs_bool_t,
  CS_TYPE_cs_point_t,
  CS_TYPE_void
} cs_type_t;


/*=============================================================================
 * Définitions de variables globales
 *============================================================================*/

extern cs_int_t  cs_glob_base_rang;     /* Rang du processus dans le groupe   */
extern cs_int_t  cs_glob_base_nbr;      /* Nombre de processus dans le groupe */

#if defined(_CS_HAVE_MPI)
extern MPI_Comm      cs_glob_base_mpi_comm;            /* Intra-communicateur */
#endif


/* Variables globales associées à l'instrumentation */

#if defined(_CS_HAVE_MPI) && defined(_CS_HAVE_MPE)
extern int  cs_glob_mpe_broadcast_a;
extern int  cs_glob_mpe_broadcast_b;
extern int  cs_glob_mpe_synchro_a;
extern int  cs_glob_mpe_synchro_b;
extern int  cs_glob_mpe_send_a;
extern int  cs_glob_mpe_send_b;
extern int  cs_glob_mpe_rcv_a;
extern int  cs_glob_mpe_rcv_b;
extern int  cs_glob_mpe_reduce_a;
extern int  cs_glob_mpe_reduce_b;
extern int  cs_glob_mpe_compute_a;
extern int  cs_glob_mpe_compute_b;
#endif


/*============================================================================
 *  Prototypes de fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction d'arret depuis du code Fortran
 *
 * Interface Fortran :
 *
 * SUBROUTINE CSEXIT (STATUT)
 * *****************
 *
 * INTEGER          STATUT      : --> : 0 pour succès, 1 ou + pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (csexit, CSEXIT)
(
  const cs_int_t  *const statut
);


/*----------------------------------------------------------------------------
 * Temps CPU écoulé depuis le début de l'exécution
 *
 * Interface Fortran :
 *
 * SUBROUTINE DMTMPS (TCPU)
 * *****************
 *
 * DOUBLE PRECISION TCPU        : --> : temps CPU (utilisateur + système)
 *----------------------------------------------------------------------------*/

void CS_PROCF (dmtmps, DMTMPS)
(
  cs_real_t  *const tcpu
);


/*=============================================================================
 * Prototypes de fonctions
 *============================================================================*/


#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Initialisation MPI ; les variables globales `cs_glob_base_nbr' indiquant
 *  le nombre de processus Code_Saturne et `cs_glob_base_rang' indiquant le
 *  rang du processus courant parmi les processus Code_Saturne sont
 * (re)positionnées par cette fonction.
 *----------------------------------------------------------------------------*/

void cs_base_mpi_init
(
 int         *argc,      /* --> Nombre d'arguments ligne de commandes        */
 char      ***argv,      /* --> Tableau des arguments ligne de commandes     */
 cs_int_t     rang_deb   /* --> Rang du premier processus du groupe
                          *     dans MPI_COMM_WORLD                          */
);


/*----------------------------------------------------------------------------
 *  Finalisation MPI
 *----------------------------------------------------------------------------*/

void cs_base_mpi_fin
(
 void
);


#endif /* defined(_CS_HAVE_MPI) */


/*----------------------------------------------------------------------------
 * Fonction d'arret
 *----------------------------------------------------------------------------*/

void cs_exit
(
  const cs_int_t  statut
);


/*----------------------------------------------------------------------------
 * Fonction initialisant la gestion des erreurs et des signaux
 *----------------------------------------------------------------------------*/

void cs_base_erreur_init
(
 void
);


/*----------------------------------------------------------------------------
 * Fonction initialisant la gestion de contrôle de la mémoire allouée
 *----------------------------------------------------------------------------*/

void cs_base_mem_init
(
 void
);


/*----------------------------------------------------------------------------
 * Fonction terminant la gestion de contrôle de la mémoire allouée
 * et affichant le bilan de la mémoire consommée.
 *----------------------------------------------------------------------------*/

void cs_base_mem_fin
(
 void
);


/*----------------------------------------------------------------------------
 * Fonction affichant le bilan du temps de calcul et temps écoulé.
 *----------------------------------------------------------------------------*/

void cs_base_bilan_temps
(
 void
);


/*----------------------------------------------------------------------------
 * Fonction affichant le bilan du temps de calcul et temps écoulé.
 *----------------------------------------------------------------------------*/

void cs_base_info_systeme
(
 void
);


/*----------------------------------------------------------------------------
 * Modification du comportement des fonctions bft_printf() par défaut
 *----------------------------------------------------------------------------*/

void cs_base_bft_printf_set
(
 void
);


/*----------------------------------------------------------------------------
 * Fonction d'impression d'un message "avertissement"
 *----------------------------------------------------------------------------*/

void cs_base_warn
(
 const char  *file_name,
 const int    line_num
);


/*----------------------------------------------------------------------------
 * Conversion d'une chaîne de l'API Fortran vers l'API C,
 * (avec suppression des blancs en début ou fin de chaîne).
 *----------------------------------------------------------------------------*/

char  * cs_base_chaine_f_vers_c_cree
(
 const char      *const chaine,             /* --> Chaîne Fortran             */
 const cs_int_t         longueur            /* --> Longueur de la chaîne      */
);


/*----------------------------------------------------------------------------
 *  Libération d'une chaîne convertie de l'API Fortran vers l'API C
 *----------------------------------------------------------------------------*/

char  * cs_base_chaine_f_vers_c_detruit
(
 char  * chaine                             /* --> Chaîne C                   */
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_BASE_H__ */
