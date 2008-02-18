/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
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

#ifndef __CS_COUPLAGE_H__
#define __CS_COUPLAGE_H__

/*============================================================================
 * Définitions, variables globales, et fonctions associées aux couplages
 * du code avec lui-même ou avec des modules reconnus.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "cs_base.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Définitions d'énumerations
 *============================================================================*/


/*============================================================================
 * Définition de macros
 *============================================================================*/


/*============================================================================
 * Déclaration de structures
 *============================================================================*/

/*
  Pointeur associé à un couplage. La structure elle-même est déclarée
  dans le fichier "cs_couplage.c", car elle n'est pas nécessaire ailleurs.
*/

typedef struct _cs_couplage_t cs_couplage_t;


/*============================================================================
 * Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Récupération du nombre de cas de couplage
 *
 * Interface Fortran :
 *
 * SUBROUTINE NBCCPL
 * *****************
 *
 * INTEGER          NBRCPL         : <-- : nombre de couplages
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbccpl, NBCCPL)
(
 cs_int_t  *const nbrcpl              /* <-- nombre de couplages              */
);


/*----------------------------------------------------------------------------
 * Affectation des listes de cellules et de faces de bord associées à
 * un couplage, ainsi que d'un ensemble de points.
 *
 * Les cellules et faces de bord locales "support" servent de base de
 * localisation des valeurs aux cellules et faces "couplée" distantes.
 * Selon le rôle émetteur et/ou destinataire du processus courant dans le
 * couplage, certains de ces ensembles peuvent être vides ou non.
 *
 * Les valeurs aux cellules seront toujours localisées et interpolées
 * par rapport au support "cellules" distant. Les valeurs aux faces
 * de bord seront localisées et interpolées par rapport au support
 * "faces de bord" s'il existe, et par rapport au support "cellules"
 * sinon. Vu du processeur local, on affecte (généralement par
 * interpolation) des valeurs à 0 à 2 ensembles de points distants,
 * dont l'un prendra les valeurs basées sur les cellules, et l'autre
 * soit sur les cellules, soit sur les faces de bord (selon si l'on
 * a défini les faces de bord comme support ou non).
 *
 * Si les tableaux LCESUP et LFBSUP ne sont pas triés en entrée, ils
 * le seront en sortie
 *
 * Interface Fortran :
 *
 * SUBROUTINE DEFCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NCESUP         : --> : nombre de cellules support
 * INTEGER          NFBSUP         : --> : nombre de faces de bord support
 * INTEGER          NCECPL         : --> : nombre de cellules couplées
 * INTEGER          NFBCPL         : --> : nombre de faces de bord couplées
 * INTEGER          LCESUP(NCESUP) : --> : liste des cellules support
 * INTEGER          LFBSUP(NFBSUP) : --> : liste des faces de bord support
 * INTEGER          LCECPL(NCECPL) : --> : liste des cellules couplées
 * INTEGER          LFBCPL(NFBCPL) : --> : liste des faces de bord couplées
 *----------------------------------------------------------------------------*/

void CS_PROCF (defcpl, DEFCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
 const cs_int_t  *const ncesup,       /* --> nombre de cellules support       */
 const cs_int_t  *const nfbsup,       /* --> nombre de faces de bord support  */
 const cs_int_t  *const ncecpl,       /* --> nombre de cellules couplées      */
 const cs_int_t  *const nfbcpl,       /* --> nombre de faces de bord couplées */
       cs_int_t         lcesup[],     /* <-> liste des cellules support       */
       cs_int_t         lfbsup[],     /* <-> liste des faces de bord support  */
 const cs_int_t         lcecpl[],     /* --> liste des cellules couplées      */
 const cs_int_t         lfbcpl[]      /* --> liste des faces de bord couplées */
);


/*----------------------------------------------------------------------------
 * Récupération des nombres de cellules et faces de bord support, couplées,
 * et non localisées associées à un couplage
 *
 * Interface Fortran :
 *
 * SUBROUTINE NBECPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NCESUP         : <-- : nombre de cellules support
 * INTEGER          NFBSUP         : <-- : nombre de faces de bord support
 * INTEGER          NCECPL         : <-- : nombre de cellules couplées
 * INTEGER          NFBCPL         : <-- : nombre de faces de bord couplées
 * INTEGER          NCENCP         : <-- : nombre de cellules non couplées
 *                                 :     : car non localisées
 * INTEGER          NFBNCP         : <-- : nombre de faces de bord non
 *                                 :     : couplées car non localisées
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbecpl, NBECPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
       cs_int_t  *const ncesup,       /* <-- nombre de cellules support       */
       cs_int_t  *const nfbsup,       /* <-- nombre de faces de bord support  */
       cs_int_t  *const ncecpl,       /* <-- nombre de cellules couplées      */
       cs_int_t  *const nfbcpl,       /* <-- nombre de faces de bord couplées */
       cs_int_t  *const ncencp,       /* <-- nombre de cellules non couplées
                                       *     car non localisées               */
       cs_int_t  *const nfbncp        /* <-- nombre de faces de bord non
                                       *     couplées car non localisées      */
);


/*----------------------------------------------------------------------------
 * Récupération des listes de cellules et de faces de bord couplées
 * (i.e. réceptrices) associées à un couplage.
 *
 * Le nombre de cellules et faces de bord, obtenus via NBECPL(), sont
 * fournis à des fins de vérification de cohérence des arguments.
 *
 * Interface Fortran :
 *
 * SUBROUTINE LELCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NCECPL         : --> : nombre de cellules couplées
 * INTEGER          NFBCPL         : --> : nombre de faces de bord couplées
 * INTEGER          LCECPL(*)      : <-- : liste des cellules couplées
 * INTEGER          LFBCPL(*)      : <-- : liste des faces de bord couplées
 *----------------------------------------------------------------------------*/

void CS_PROCF (lelcpl, LELCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du cas de couplage        */
 const cs_int_t  *const ncecpl,       /* --> nombre de cellules couplées      */
 const cs_int_t  *const nfbcpl,       /* --> nombre de faces de bord couplées */
       cs_int_t  *const lcecpl,       /* <-- liste des cellules couplées      */
       cs_int_t  *const lfbcpl        /* <-- liste des faces de bord couplées */
);


/*----------------------------------------------------------------------------
 * Récupération des listes de cellules et de faces de bord non couplées
 * (i.e. réceptrices mais non localisées) associées à un couplage
 *
 * Le nombre de cellules et faces de bord, obtenus via NBECPL(), sont
 * fournis à des fins de vérification de cohérence des arguments.
 *
 * Interface Fortran :
 *
 * SUBROUTINE LENCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NCENCP         : <-- : nombre de cellules non couplées
 *                                 :     : car non localisées
 * INTEGER          NFBNCP         : <-- : nombre de faces de bord non
 *                                 :     : couplées car non localisées
 * INTEGER          LCENCP(*)      : <-- : liste des cellules non couplées
 * INTEGER          LFBNCP(*)      : <-- : liste des faces de bord non couplées
 *----------------------------------------------------------------------------*/

void CS_PROCF (lencpl, LENCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du cas de couplage        */
 const cs_int_t  *const ncencp,       /* --> nombre de cellules non couplées
                                       *     car non localisées               */
 const cs_int_t  *const nfbncp,       /* --> nombre de faces de bord non
                                       *     couplées car non localisées      */
       cs_int_t  *const lcencp,       /* <-- liste des cellules non couplées  */
       cs_int_t  *const lfbncp        /* <-- liste des faces de bord non
                                       *     couplées                         */
);


/*----------------------------------------------------------------------------
 * Récupération du nombre de points distants associés à un couplage
 * et localisés par rapport au domaine local
 *
 * Interface Fortran :
 *
 * SUBROUTINE NPDCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NCEDIS         : <-- : nombre de cellules distantes
 * INTEGER          NFBDIS         : <-- : nombre de faces de bord distantes
 *----------------------------------------------------------------------------*/

void CS_PROCF (npdcpl, NPDCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
       cs_int_t  *const ncedis,       /* <-- nombre de cellules distantes     */
       cs_int_t  *const nfbdis        /* <-- nombre de faces de bord dist.    */
);


/*----------------------------------------------------------------------------
 * Récupération des coordonnées des points distants affectés à un
 * couplage et une liste de points, ainsi que les numéros et le
 * type d'élément (cellules ou faces) "contenant" ces points.
 *
 * Le nombre de points distants NBRPTS doit être égal à l'un des arguments
 * NCEDIS ou NFBDIS retournés par NPDCPL(), et est fourni ici à des fins
 * de vérification de cohérence avec les arguments NUMCPL et ITYSUP.
 *
 * Interface Fortran :
 *
 * SUBROUTINE COOCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NBRPTS         : --> : nombre de points distants
 * INTEGER          ITYDIS         : --> : 1 : accès aux points affectés
 *                                 :     :     aux cellules distantes
 *                                 :     : 2 : accès aux points affectés
 *                                 :     :     aux faces de bord distantes
 * INTEGER          ITYLOC         : <-- : 1 : localisation par rapport
 *                                 :     :     aux cellules locales
 *                                 :     : 2 : localisation par rapport
 *                                 :     :     aux faces de bord locales
 * INTEGER          LOCPTS(*)      : <-- : numéro du "contenant" associé à
 *                                 :     :   chaque point
 * DOUBLE PRECISION COOPTS(3,*)    : <-- : coordonnées des points distants
 *----------------------------------------------------------------------------*/

void CS_PROCF (coocpl, COOCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
 const cs_int_t  *const nbrpts,       /* --> nombre de points distants        */
 const cs_int_t  *const itydis,       /* --> 1 : accès aux points affectés
                                       *         aux cellules distantes
                                       *     2 : accès aux points affectés
                                       *         aux faces de bord distantes  */
       cs_int_t  *const ityloc,       /* <-- 1 : localisation par rapport
                                       *         aux cellules locales
                                       *     2 : localisation par rapport
                                       *         aux faces de bord locales    */
       cs_int_t  *const locpts,       /* <-- liste des mailles associées      */
       cs_real_t *const coopts        /* <-- coord. des points à localiser    */
);


/*----------------------------------------------------------------------------
 * Echange d'une variable associée à un ensemble de points et à un couplage.
 *
 * Interface Fortran :
 *
 * SUBROUTINE VARCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NBRDIS         : --> : Nombre de valeurs à envoyer
 * INTEGER          NBRLOC         : --> : Nombre de valeurs à recevoir
 * INTEGER          ITYVAR         : --> : 1 : variables aux cellules
 *                                 :     : 2 : variables aux faces de bord
 * DOUBLE PRECISION VARDIS(*) )    : --> : variable distante (à envoyer)
 * DOUBLE PRECISION VARLOC(*) )    : <-- : variable locale (à recevoir)
 *----------------------------------------------------------------------------*/

void CS_PROCF (varcpl, VARCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
 const cs_int_t  *const nbrdis,       /* --> nombre de valeurs à envoyer      */
 const cs_int_t  *const nbrloc,       /* --> nombre de valeurs à recevoir     */
 const cs_int_t  *const ityvar,       /* --> 1 : variables aux cellules
                                       *     2 : variables aux faces de bord  */
       cs_real_t *const vardis,       /* --> variable distante (à envoyer)    */
       cs_real_t *const varloc        /* <-- variable locale (à recevoir)     */
);


/*----------------------------------------------------------------------------
 * Echange de tableaux d'entiers associés à un couplage.
 *
 * On suppose que les tableaux à échanger sont de même taille et contiennent
 * les mêmes valeurs sur chaque groupe de processus (locaux et distants).
 *
 * Interface Fortran :
 *
 * SUBROUTINE TBICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NBRDIS         : --> : Nombre de valeurs à envoyer
 * INTEGER          NBRLOC         : --> : Nombre de valeurs à recevoir
 * INTEGER          TABDIS(*) )    : --> : valeurs distantes (à envoyer)
 * INTEGER          TABLOC(*) )    : --> : valeurs locales (à recevoir)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbicpl, TBICPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
 const cs_int_t  *const nbrdis,       /* --> nombre de valeurs à envoyer      */
 const cs_int_t  *const nbrloc,       /* --> nombre de valeurs à recevoir     */
       cs_int_t  *const vardis,       /* --> variable distante (à envoyer)    */
       cs_int_t  *const varloc        /* <-- variable locale (à recevoir)     */
);


/*----------------------------------------------------------------------------
 * Echange de tableaux de réels associés à un couplage.
 *
 * On suppose que les tableaux à échanger sont de même taille et contiennent
 * les mêmes valeurs sur chaque groupe de processus (locaux et distants).
 *
 * Interface Fortran :
 *
 * SUBROUTINE TBRCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numéro du couplage
 * INTEGER          NBRDIS         : --> : Nombre de valeurs à envoyer
 * INTEGER          NBRLOC         : --> : Nombre de valeurs à recevoir
 * DOUBLE PRECISION TABDIS(*) )    : --> : valeurs distantes (à envoyer)
 * DOUBLE PRECISION TABLOC(*) )    : --> : valeurs locales (à recevoir)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbrcpl, TBRCPL)
(
 const cs_int_t  *const numcpl,       /* --> numéro du couplage               */
 const cs_int_t  *const nbrdis,       /* --> nombre de valeurs à envoyer      */
 const cs_int_t  *const nbrloc,       /* --> nombre de valeurs à recevoir     */
       cs_real_t *const vardis,       /* --> variable distante (à envoyer)    */
       cs_real_t *const varloc        /* <-- variable locale (à recevoir)     */
);


/*============================================================================
 *  Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ajout d'un couplage.
 *
 * On autorise les couplages soit avec des groupes de processus totalement
 * distincts du groupe principal (correspondant à cs_glob_base_mpi_comm),
 * soit avec ce même groupe.
 *----------------------------------------------------------------------------*/

void  cs_couplage_ajoute
(
 const cs_int_t   rang_deb            /* --> rang du premier processus couplé */
);


/*----------------------------------------------------------------------------
 * Suppression des couplages.
 *----------------------------------------------------------------------------*/

void cs_couplage_detruit_tout
(
 void
);


#endif /* __CS_COUPLAGE_H__ */
