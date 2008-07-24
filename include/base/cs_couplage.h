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

#ifndef __CS_COUPLAGE_H__
#define __CS_COUPLAGE_H__

/*============================================================================
 * Definitions, variables globales, et fonctions associees aux couplages
 * du code avec lui-meme ou avec des modules reconnus.
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
 * Definitions d'enumerations
 *============================================================================*/


/*============================================================================
 * Definition de macros
 *============================================================================*/


/*============================================================================
 * Declaration de structures
 *============================================================================*/

/*
  Pointeur associe a un couplage. La structure elle-meme est declaree
  dans le fichier "cs_couplage.c", car elle n'est pas necessaire ailleurs.
*/

typedef struct _cs_couplage_t cs_couplage_t;


/*============================================================================
 * Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Recuperation du nombre de cas de couplage
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
 * Affectation des listes de cellules et de faces de bord associees a
 * un couplage, ainsi que d'un ensemble de points.
 *
 * Les cellules et faces de bord locales "support" servent de base de
 * localisation des valeurs aux cellules et faces "couplee" distantes.
 * Selon le r^ole emetteur et/ou destinataire du processus courant dans le
 * couplage, certains de ces ensembles peuvent etre vides ou non.
 *
 * Les valeurs aux cellules seront toujours localisees et interpolees
 * par rapport au support "cellules" distant. Les valeurs aux faces
 * de bord seront localisees et interpolees par rapport au support
 * "faces de bord" s'il existe, et par rapport au support "cellules"
 * sinon. Vu du processeur local, on affecte (generalement par
 * interpolation) des valeurs a 0 a 2 ensembles de points distants,
 * dont l'un prendra les valeurs basees sur les cellules, et l'autre
 * soit sur les cellules, soit sur les faces de bord (selon si l'on
 * a defini les faces de bord comme support ou non).
 *
 * Si les tableaux LCESUP et LFBSUP ne sont pas tries en entree, ils
 * le seront en sortie
 *
 * Interface Fortran :
 *
 * SUBROUTINE DEFCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NCESUP         : --> : nombre de cellules support
 * INTEGER          NFBSUP         : --> : nombre de faces de bord support
 * INTEGER          NCECPL         : --> : nombre de cellules couplees
 * INTEGER          NFBCPL         : --> : nombre de faces de bord couplees
 * INTEGER          LCESUP(NCESUP) : --> : liste des cellules support
 * INTEGER          LFBSUP(NFBSUP) : --> : liste des faces de bord support
 * INTEGER          LCECPL(NCECPL) : --> : liste des cellules couplees
 * INTEGER          LFBCPL(NFBCPL) : --> : liste des faces de bord couplees
 *----------------------------------------------------------------------------*/

void CS_PROCF (defcpl, DEFCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
 const cs_int_t  *const ncesup,       /* --> nombre de cellules support       */
 const cs_int_t  *const nfbsup,       /* --> nombre de faces de bord support  */
 const cs_int_t  *const ncecpl,       /* --> nombre de cellules couplees      */
 const cs_int_t  *const nfbcpl,       /* --> nombre de faces de bord couplees */
       cs_int_t         lcesup[],     /* <-> liste des cellules support       */
       cs_int_t         lfbsup[],     /* <-> liste des faces de bord support  */
 const cs_int_t         lcecpl[],     /* --> liste des cellules couplees      */
 const cs_int_t         lfbcpl[]      /* --> liste des faces de bord couplees */
);


/*----------------------------------------------------------------------------
 * Recuperation des nombres de cellules et faces de bord support, couplees,
 * et non localisees associees a un couplage
 *
 * Interface Fortran :
 *
 * SUBROUTINE NBECPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NCESUP         : <-- : nombre de cellules support
 * INTEGER          NFBSUP         : <-- : nombre de faces de bord support
 * INTEGER          NCECPL         : <-- : nombre de cellules couplees
 * INTEGER          NFBCPL         : <-- : nombre de faces de bord couplees
 * INTEGER          NCENCP         : <-- : nombre de cellules non couplees
 *                                 :     : car non localisees
 * INTEGER          NFBNCP         : <-- : nombre de faces de bord non
 *                                 :     : couplees car non localisees
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbecpl, NBECPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
       cs_int_t  *const ncesup,       /* <-- nombre de cellules support       */
       cs_int_t  *const nfbsup,       /* <-- nombre de faces de bord support  */
       cs_int_t  *const ncecpl,       /* <-- nombre de cellules couplees      */
       cs_int_t  *const nfbcpl,       /* <-- nombre de faces de bord couplees */
       cs_int_t  *const ncencp,       /* <-- nombre de cellules non couplees
                                       *     car non localisees               */
       cs_int_t  *const nfbncp        /* <-- nombre de faces de bord non
                                       *     couplees car non localisees      */
);


/*----------------------------------------------------------------------------
 * Recuperation des listes de cellules et de faces de bord couplees
 * (i.e. receptrices) associees a un couplage.
 *
 * Le nombre de cellules et faces de bord, obtenus via NBECPL(), sont
 * fournis a des fins de verification de coherence des arguments.
 *
 * Interface Fortran :
 *
 * SUBROUTINE LELCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NCECPL         : --> : nombre de cellules couplees
 * INTEGER          NFBCPL         : --> : nombre de faces de bord couplees
 * INTEGER          LCECPL(*)      : <-- : liste des cellules couplees
 * INTEGER          LFBCPL(*)      : <-- : liste des faces de bord couplees
 *----------------------------------------------------------------------------*/

void CS_PROCF (lelcpl, LELCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du cas de couplage        */
 const cs_int_t  *const ncecpl,       /* --> nombre de cellules couplees      */
 const cs_int_t  *const nfbcpl,       /* --> nombre de faces de bord couplees */
       cs_int_t  *const lcecpl,       /* <-- liste des cellules couplees      */
       cs_int_t  *const lfbcpl        /* <-- liste des faces de bord couplees */
);


/*----------------------------------------------------------------------------
 * Recuperation des listes de cellules et de faces de bord non couplees
 * (i.e. receptrices mais non localisees) associees a un couplage
 *
 * Le nombre de cellules et faces de bord, obtenus via NBECPL(), sont
 * fournis a des fins de verification de coherence des arguments.
 *
 * Interface Fortran :
 *
 * SUBROUTINE LENCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NCENCP         : <-- : nombre de cellules non couplees
 *                                 :     : car non localisees
 * INTEGER          NFBNCP         : <-- : nombre de faces de bord non
 *                                 :     : couplees car non localisees
 * INTEGER          LCENCP(*)      : <-- : liste des cellules non couplees
 * INTEGER          LFBNCP(*)      : <-- : liste des faces de bord non couplees
 *----------------------------------------------------------------------------*/

void CS_PROCF (lencpl, LENCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du cas de couplage        */
 const cs_int_t  *const ncencp,       /* --> nombre de cellules non couplees
                                       *     car non localisees               */
 const cs_int_t  *const nfbncp,       /* --> nombre de faces de bord non
                                       *     couplees car non localisees      */
       cs_int_t  *const lcencp,       /* <-- liste des cellules non couplees  */
       cs_int_t  *const lfbncp        /* <-- liste des faces de bord non
                                       *     couplees                         */
);


/*----------------------------------------------------------------------------
 * Recuperation du nombre de points distants associes a un couplage
 * et localises par rapport au domaine local
 *
 * Interface Fortran :
 *
 * SUBROUTINE NPDCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NCEDIS         : <-- : nombre de cellules distantes
 * INTEGER          NFBDIS         : <-- : nombre de faces de bord distantes
 *----------------------------------------------------------------------------*/

void CS_PROCF (npdcpl, NPDCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
       cs_int_t  *const ncedis,       /* <-- nombre de cellules distantes     */
       cs_int_t  *const nfbdis        /* <-- nombre de faces de bord dist.    */
);


/*----------------------------------------------------------------------------
 * Recuperation des coordonnees des points distants affectes a un
 * couplage et une liste de points, ainsi que les numeros et le
 * type d'element (cellules ou faces) "contenant" ces points.
 *
 * Le nombre de points distants NBRPTS doit etre egal a l'un des arguments
 * NCEDIS ou NFBDIS retournes par NPDCPL(), et est fourni ici a des fins
 * de verification de coherence avec les arguments NUMCPL et ITYSUP.
 *
 * Interface Fortran :
 *
 * SUBROUTINE COOCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NBRPTS         : --> : nombre de points distants
 * INTEGER          ITYDIS         : --> : 1 : acces aux points affectes
 *                                 :     :     aux cellules distantes
 *                                 :     : 2 : acces aux points affectes
 *                                 :     :     aux faces de bord distantes
 * INTEGER          ITYLOC         : <-- : 1 : localisation par rapport
 *                                 :     :     aux cellules locales
 *                                 :     : 2 : localisation par rapport
 *                                 :     :     aux faces de bord locales
 * INTEGER          LOCPTS(*)      : <-- : numero du "contenant" associe a
 *                                 :     :   chaque point
 * DOUBLE PRECISION COOPTS(3,*)    : <-- : coordonnees des points distants
 *----------------------------------------------------------------------------*/

void CS_PROCF (coocpl, COOCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
 const cs_int_t  *const nbrpts,       /* --> nombre de points distants        */
 const cs_int_t  *const itydis,       /* --> 1 : acces aux points affectes
                                       *         aux cellules distantes
                                       *     2 : acces aux points affectes
                                       *         aux faces de bord distantes  */
       cs_int_t  *const ityloc,       /* <-- 1 : localisation par rapport
                                       *         aux cellules locales
                                       *     2 : localisation par rapport
                                       *         aux faces de bord locales    */
       cs_int_t  *const locpts,       /* <-- liste des mailles associees      */
       cs_real_t *const coopts        /* <-- coord. des points a localiser    */
);


/*----------------------------------------------------------------------------
 * Echange d'une variable associee a un ensemble de points et a un couplage.
 *
 * Interface Fortran :
 *
 * SUBROUTINE VARCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NBRDIS         : --> : Nombre de valeurs a envoyer
 * INTEGER          NBRLOC         : --> : Nombre de valeurs a recevoir
 * INTEGER          ITYVAR         : --> : 1 : variables aux cellules
 *                                 :     : 2 : variables aux faces de bord
 * DOUBLE PRECISION VARDIS(*) )    : --> : variable distante (a envoyer)
 * DOUBLE PRECISION VARLOC(*) )    : <-- : variable locale (a recevoir)
 *----------------------------------------------------------------------------*/

void CS_PROCF (varcpl, VARCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
 const cs_int_t  *const nbrdis,       /* --> nombre de valeurs a envoyer      */
 const cs_int_t  *const nbrloc,       /* --> nombre de valeurs a recevoir     */
 const cs_int_t  *const ityvar,       /* --> 1 : variables aux cellules
                                       *     2 : variables aux faces de bord  */
       cs_real_t *const vardis,       /* --> variable distante (a envoyer)    */
       cs_real_t *const varloc        /* <-- variable locale (a recevoir)     */
);


/*----------------------------------------------------------------------------
 * Echange de tableaux d'entiers associes a un couplage.
 *
 * On suppose que les tableaux a echanger sont de meme taille et contiennent
 * les memes valeurs sur chaque groupe de processus (locaux et distants).
 *
 * Interface Fortran :
 *
 * SUBROUTINE TBICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NBRDIS         : --> : Nombre de valeurs a envoyer
 * INTEGER          NBRLOC         : --> : Nombre de valeurs a recevoir
 * INTEGER          TABDIS(*) )    : --> : valeurs distantes (a envoyer)
 * INTEGER          TABLOC(*) )    : --> : valeurs locales (a recevoir)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbicpl, TBICPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
 const cs_int_t  *const nbrdis,       /* --> nombre de valeurs a envoyer      */
 const cs_int_t  *const nbrloc,       /* --> nombre de valeurs a recevoir     */
       cs_int_t  *const vardis,       /* --> variable distante (a envoyer)    */
       cs_int_t  *const varloc        /* <-- variable locale (a recevoir)     */
);


/*----------------------------------------------------------------------------
 * Echange de tableaux de reels associes a un couplage.
 *
 * On suppose que les tableaux a echanger sont de meme taille et contiennent
 * les memes valeurs sur chaque groupe de processus (locaux et distants).
 *
 * Interface Fortran :
 *
 * SUBROUTINE TBRCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : numero du couplage
 * INTEGER          NBRDIS         : --> : Nombre de valeurs a envoyer
 * INTEGER          NBRLOC         : --> : Nombre de valeurs a recevoir
 * DOUBLE PRECISION TABDIS(*) )    : --> : valeurs distantes (a envoyer)
 * DOUBLE PRECISION TABLOC(*) )    : --> : valeurs locales (a recevoir)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbrcpl, TBRCPL)
(
 const cs_int_t  *const numcpl,       /* --> numero du couplage               */
 const cs_int_t  *const nbrdis,       /* --> nombre de valeurs a envoyer      */
 const cs_int_t  *const nbrloc,       /* --> nombre de valeurs a recevoir     */
       cs_real_t *const vardis,       /* --> variable distante (a envoyer)    */
       cs_real_t *const varloc        /* <-- variable locale (a recevoir)     */
);


/*============================================================================
 *  Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ajout d'un couplage.
 *
 * On autorise les couplages soit avec des groupes de processus totalement
 * distincts du groupe principal (correspondant a cs_glob_base_mpi_comm),
 * soit avec ce meme groupe.
 *----------------------------------------------------------------------------*/

void  cs_couplage_ajoute
(
 const cs_int_t   rang_deb            /* --> rang du premier processus couple */
);


/*----------------------------------------------------------------------------
 * Suppression des couplages.
 *----------------------------------------------------------------------------*/

void cs_couplage_detruit_tout
(
 void
);


#endif /* __CS_COUPLAGE_H__ */
