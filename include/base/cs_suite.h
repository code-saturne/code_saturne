/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_SUITE_H__
#define __CS_SUITE_H__

/*============================================================================
 *  Gestion des fichiers suite
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "cs_base.h"


/*============================================================================
 *  Definitions d'enumerations
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Type de fichier suite
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SUITE_TYPE_ASCII,           /* Fichier suite ASCII                       */
  CS_SUITE_TYPE_BINAIRE          /* Fichier suite binaire                     */

} cs_suite_type_t;


/*----------------------------------------------------------------------------
 *  Fichier en lecture ou ecriture
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SUITE_MODE_LECTURE,         /* Communication en reception                */
  CS_SUITE_MODE_ECRITURE         /* Communication en emission                 */

} cs_suite_mode_t;


/*----------------------------------------------------------------------------
 *  Type de support de maillage associe a une rubrique
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SUITE_SUPPORT_SCAL,         /* Scalaire (sans support)                   */
  CS_SUITE_SUPPORT_CEL,          /* Valeurs associees aux cellules            */
  CS_SUITE_SUPPORT_FAC_INT,      /* Valeurs associees aux faces internes      */
  CS_SUITE_SUPPORT_FAC_BRD,      /* Valeurs associees aux faces de bord       */
  CS_SUITE_SUPPORT_SOM           /* Valeurs associees aux sommets             */

} cs_suite_support_t;


/*============================================================================
 *  Definition de macros
 *============================================================================*/

/* Codes d'erreur */

#define CS_SUITE_SUCCES          0 /* Reussite */
#define CS_SUITE_ERR_NUM_FIC    -1 /* Pas de suite du numero indique */
#define CS_SUITE_ERR_TYPE_FIC   -2 /* Type de fichier incorrect */
#define CS_SUITE_ERR_SUPPORT    -3 /* Support indefini/de dimension differente */
#define CS_SUITE_ERR_TYPE_VAL   -4 /* Type de valeur inconnu ou imprevu */
#define CS_SUITE_ERR_NBR_VAL    -5 /* Nombre de valeurs ne correspond pas */
#define CS_SUITE_ERR_MODE       -6 /* Mode d'ouverure du fichier incompatible */
#define CS_SUITE_ERR_EXISTE     -7 /* Enregistrement non disponible */


/*============================================================================
 *  Declaration de structures
 *============================================================================*/

/*
  Pointeur associe a un fichier suite. La structure elle-meme est declaree
  dans le fichier "cs_suite.c", car elle n'est pas necessaire ailleurs.
*/

typedef struct _cs_suite_t cs_suite_t;


/*=============================================================================
 * Definitions de variables globales
 *============================================================================*/


/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ouverture d'un fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE OPNSUI (NOMSUI, LNGNOM, IREAWR, IFORMA, NUMSUI)
 * *****************
 *
 * CHARACTER*       NOMSUI      : --> : Nom du fichier suite
 * INTEGER          LNGNOM      : --> : Longueur du nom du fichier suite
 * INTEGER          IREAWR      : --> : 1 pour lecture, 2 pour ecriture
 * INTEGER          IFORMA      : --> : 0 pour binaire, 1 pour formate
 * INTEGER          NUMSUI      : <-- : Numero du fichier suite ouvert
 * INTEGER          IERROR      : <-- : 0 pour succes, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (opnsui, OPNSUI)
(
 const char       *const nomsui,  /* --> Nom du fichier                       */
 const cs_int_t   *const lngnom,  /* --> Longueur du nom                      */
 const cs_int_t   *const ireawr,  /* --> 1 pour lecture, 2 pour ecriture      */
 const cs_int_t   *const iforma,  /* --> 0 pour binaire, 1 pour formate       */
       cs_int_t   *const numsui,  /* <-- Numero du ficher suite ouvert        */
       cs_int_t   *const ierror   /* <-- 0 pour succes, < 0 pour erreur       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' eventuels F77, */
                                  /*     inutilises lors de l'appel mais      */
                                  /*     places par de nombreux compilateurs) */
);


/*----------------------------------------------------------------------------
 * Fermeture d'un fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE CLSSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : numero du fichier suite a fermer
 * INTEGER          IERROR      : <-- : 0 pour succes, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (clssui, CLSSUI)
(
 const cs_int_t   *const numsui,  /* --> Numero du ficher suite a fermer      */
       cs_int_t   *const ierror   /* <-- Numero du ficher suite ouvert        */
);


/*----------------------------------------------------------------------------
 *  Verification du support associe a un fichier suite ;
 *  On renvoie pour chaque type d'entite 1 si le nombre d'entites associees
 *  au fichier suite correspond au nombre d'entites en cours (et donc que
 *  l'on considere que le support est bien le meme), 0 sinon.
 *
 * Interface Fortran :
 *
 * SUBROUTINE TSTSUI (NUMSUI, INDCEL, INDFAC, INDFBR, INDSOM)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numero du fichier suite
 * INTEGER          INDCEL      : <-- : Indicateur corresp. cellules
 * INTEGER          INDFAC      : <-- : Indicateur corresp. faces internes
 * INTEGER          INDFBR      : <-- : Indicateur corresp. faces de bord
 * INTEGER          INDSOM      : <-- : Indicateur corresp. sommets
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstsui, TSTSUI)
(
 const cs_int_t  *const numsui,   /* --> Numero du fichier suite              */
       cs_int_t  *const indcel,   /* <-- Indicateur corresp. cellules         */
       cs_int_t  *const indfac,   /* <-- Indicateur corresp. faces internes   */
       cs_int_t  *const indfbr,   /* <-- Indicateur corresp. faces de bord    */
       cs_int_t  *const indsom    /* <-- Indicateur corresp. sommets          */
);


/*----------------------------------------------------------------------------
 *  Affichage de l'index associe a un fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE INFSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numero du fichier suite
 *----------------------------------------------------------------------------*/

void CS_PROCF (infsui, INFSUI)
(
 const cs_int_t  *const numsui    /* --> Numero du fichier suite              */
);


/*----------------------------------------------------------------------------
 * Lecture d'une rubrique sur fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE LECSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numero du fichier suite
 * CHARACTER*       NOMRUB      : --> : Nom de la rubrique
 * INTEGER          LNGNOM      : --> : Longueur du nom de la rubrique
 * INTEGER          ITYSUP      : --> : Type de support :
 *                              :     :  0 : scalaire (pas de support)
 *                              :     :  1 : cellules
 *                              :     :  2 : faces internes
 *                              :     :  3 : faces de bord
 *                              :     :  4 : sommets (si disponibles)
 * INTEGER          NBVENT      : --> : Nb. valeurs par entite de support
 * INTEGER          IRTYPE      : --> : 1 pour entiers, 2 pour double precision
 * (?)              TABVAR      : <-> : Tableau des valeurs a lire
 * INTEGER          IERROR      : <-- : 0 pour succes, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (lecsui, LECSUI)
(
 const cs_int_t   *const numsui,  /* --> Numero du fichier suite              */
 const char       *const nomrub,  /* --> Nom de la rubrique                   */
 const cs_int_t   *const lngnom,  /* --> Longueur du nom de la rubrique       */
 const cs_int_t   *const itysup,  /* --> Type de support (voir ci-dessus)     */
 const cs_int_t   *const nbvent,  /* --> Nb. valeurs par entite du support    */
 const cs_int_t   *const irtype,  /* --> 1 pour entiers, 2 pour double prec.  */
       void       *const tabvar,  /* <-- Tableur des valeurs a lire           */
       cs_int_t   *const ierror   /* <-- 0 pour succes, < 0 pour erreur       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' eventuels F77, */
                                  /*     inutilises lors de l'appel mais      */
                                  /*     places par de nombreux compilateurs) */
);


/*----------------------------------------------------------------------------
 * ecriture d'une rubrique sur fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE ECRSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numero du fichier suite
 * CHARACTER*       NOMRUB      : --> : Nom de la rubrique
 * INTEGER          LNGNOM      : --> : Longueur du nom de la rubrique
 * INTEGER          ITYSUP      : --> : Type de support :
 *                              :     :  0 : scalaire (pas de support)
 *                              :     :  1 : cellules
 *                              :     :  2 : faces internes
 *                              :     :  3 : faces de bord
 *                              :     :  4 : sommets (si disponibles)
 * INTEGER          NBVENT      : --> : Nb. valeurs par entite de support
 * INTEGER          IRTYPE      : --> : 1 pour entiers, 2 pour double precision
 * (?)              TABVAR      : --> : Tableau des valeurs fournies
 * INTEGER          IERROR      : <-- : 0 pour succes, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrsui, ECRSUI)
(
 const cs_int_t   *const numsui,  /* --> Numero du fichier suite              */
 const char       *const nomrub,  /* --> Nom de la rubrique                   */
 const cs_int_t   *const lngnom,  /* --> Longueur du nom de la rubrique       */
 const cs_int_t   *const itysup,  /* --> Type de support (voir ci-dessus)     */
 const cs_int_t   *const nbvent,  /* --> Nb. valeurs par entite du support    */
 const cs_int_t   *const irtype,  /* --> 1 pour entiers, 2 pour double prec.  */
 const void       *const tabvar,  /* --> Tableur des valeurs fournies         */
       cs_int_t   *const ierror   /* <-- 0 pour succes, < 0 pour erreur       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' eventuels F77, */
                                  /*     inutilises lors de l'appel mais      */
                                  /*     places par de nombreux compilateurs) */
);


/*============================================================================
 *  Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui initialise un fichier suite
 *----------------------------------------------------------------------------*/

cs_suite_t * cs_suite_cree
(
 const char             *const nom,         /* --> nom de base du fichier     */
 const cs_suite_mode_t         mode,        /* --> Lecture ou ecriture        */
 const cs_suite_type_t         type         /* --> ASCII ou binaire           */
);


/*----------------------------------------------------------------------------
 *  Fonction qui detruit la structure associee a un fichier suite (et ferme
 *  le fichier associe) ; elle renvoie un pointeur NULL.
 *----------------------------------------------------------------------------*/

cs_suite_t * cs_suite_detruit
(
 cs_suite_t * suite                         /* --> Fichier suite              */
);


/*----------------------------------------------------------------------------
 *  Fonction qui verifie le support associe a un fichier suite ;
 *  On renvoie pour chaque type d'entite CS_TRUE si le nombre d'entites
 *  associees au fichier suite correspond au nombre d'entites en cours (et
 *  donc que l'on considere que le support est bien le meme), CS_FALSE sinon.
 *----------------------------------------------------------------------------*/

void cs_suite_verif_support
(
 const cs_suite_t  *const suite,            /* --> Fichier suite              */
       cs_bool_t   *const corresp_cel,      /* <-- Corresp. cellules          */
       cs_bool_t   *const corresp_fac,      /* <-- Corresp. faces internes    */
       cs_bool_t   *const corresp_fbr,      /* <-- Corresp. faces de bord     */
       cs_bool_t   *const corresp_som       /* <-- Corresp. sommets           */
);


/*----------------------------------------------------------------------------
 *  Fonction qui affiche l'index genere lors de l'analyse du fichier
 *----------------------------------------------------------------------------*/

void cs_suite_affiche_index
(
 const cs_suite_t  *const  suite          /* --> Structure suite              */
);


/*----------------------------------------------------------------------------
 *  Fonction qui lit un enregistrement sur fichier suite ; On renvoie 0
 *  (CS_SUITE_SUCCES) en cas de succes, une valeur negative (de type
 *  CS_SUITE_ERR_xxx) en cas d'echec.
 *----------------------------------------------------------------------------*/

cs_int_t cs_suite_lit_rub
(
 const cs_suite_t          *const suite,       /* --> Ptr. structure suite    */
 const char                *const nom_rub,     /* --> Nom de la rubrique      */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent, /* --> Nb. val/point support   */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
       void                *const val          /* <-- Valeurs a lire          */
);


/*----------------------------------------------------------------------------
 *  Fonction qui ecrit un enregistrement sur fichier suite
 *----------------------------------------------------------------------------*/

void cs_suite_ecr_rub
(
       cs_suite_t          *const suite,       /* --> Ptr. structure suite    */
 const char                *const nom_rub,     /* --> Nom de la rubrique      */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent, /* --> Nb. val/point support   */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
 const void                *const val          /* --> Valeurs a ecrire        */
);


/*----------------------------------------------------------------------------
 *  Fonction qui initialise l'API Fortran
 *----------------------------------------------------------------------------*/

void cs_suite_f77_api_init
(
 void
);


/*----------------------------------------------------------------------------
 *  Fonction qui termine l'API Fortran
 *----------------------------------------------------------------------------*/

void cs_suite_f77_api_finalize
(
 void
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_SUITE_H__ */
