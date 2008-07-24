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

#ifndef __CS_POST_H__
#define __CS_POST_H__

/*============================================================================
 * Definitions, variables globales, et fonctions associees au post traitement
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include` librairies BFT et FVM
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_writer.h>

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

/* enumeration pour transmettre le type d'une donnee */

typedef enum {
  CS_POST_TYPE_cs_int_t,
  CS_POST_TYPE_cs_real_t,
  CS_POST_TYPE_int,
  CS_POST_TYPE_float,
  CS_POST_TYPE_double
} cs_post_type_t;


/*============================================================================
 * Definition de macros
 *============================================================================*/


/*============================================================================
 * Declaration de structures et types
 *============================================================================*/

/* Pointeur associe a un "writer" : cet objet correspond au choix d'un
 * nom de cas, de repertoire, et de format, ainsi qu'un indicateur precisant
 * si les maillages associes doivent dependre ou non du temps, et la
 * frequence de sortie par defaut pour les variables associees. */

typedef struct _cs_post_writer_t cs_post_writer_t;

/* Pointeur associe a un maillage de post traitement ; cet objet
 * gere le lien entre un tel maillage et les "writers" associes. */

typedef struct _cs_post_maillage_t cs_post_maillage_t;

/* Pointeur de fonction associe a un post-traitement particulier ;
 * on enregistre de telles fonctions via la fonction
 * cs_post_ajoute_var_temporelle(), et toutes les fonctions enregistrees
 * de la sorte sont appellees automatiquement par PSTVAR. */

typedef void
(cs_post_var_temporelle_t) (cs_int_t     id_instance,
                            cs_int_t     nt_cur_abs,
                            cs_real_t    t_cur_abs);


/*============================================================================
 * Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation d'un "writer" a partir des donnees du Fortran ; cet objet
 * correspond au choix d'un nom de cas, de repertoire, et de format, ainsi
 * qu'un indicateur precisant si les maillages associes doivent dependre ou
 * non du temps, et la frequence de sortie par defaut pour les
 * variables associees.
 *
 * Interface Fortran : utiliser PSTCWR (voir cs_post_util.F)
 *
 * SUBROUTINE PSTCWR (NUMGEP, NOMCAS, NOMREP, NOMFMT, OPTFMT,
 * *****************
 *                    LNMCAS, LNMFMT, LNMREP, LOPFMT,
 *                    INDMOD, NTCHR)
 *
 * INTEGER          NUMGEP      : --> : Numero du filtre a creer (< 0 pour
 *                              :     : filtre standard ou developpeur,
 *                              :     : > 0 pour filtre utilisateur)
 * CHARACTER        NOMCAS      : --> : Nom du cas associe
 * CHARACTER        NOMREP      : --> : Nom du repertoire associe
 * INTEGER          NOMFMT      : --> : Nom de format associe
 * INTEGER          OPTFMT      : --> : Options associees au format
 * INTEGER          LNMCAS      : --> : Longueur du nom du cas
 * INTEGER          LNMREP      : --> : Longueur du nom du repertoire
 * INTEGER          LNMFMT      : --> : Longueur du nom du format
 * INTEGER          LOPFMT      : --> : Longueur des options du format
 * INTEGER          INDMOD      : --> : 0 si fige, 1 si deformable,
 *                              :     : 2 si la topologie change
 * INTEGER          NTCHR       : --> : Frequence de sortie par defaut
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcw1, PSTCW1)
(
 const cs_int_t   *const numwri,  /* --> numero du writer a creer
                                   *     < 0 pour writer reserve,
                                   *     > 0 pour writer utilisateur)         */
 const char       *const nomcas,  /* --> nom du cas associe                   */
 const char       *const nomrep,  /* --> nom de repertoire associe            */
 const char       *const nomfmt,  /* --> nom de format associe                */
 const char       *const optfmt,  /* --> options associees au format          */
 const cs_int_t   *const lnmcas,  /* --> longueur du nom du cas               */
 const cs_int_t   *const lnmrep,  /* --> longueur du nom du repertoire        */
 const cs_int_t   *const lnmfmt,  /* --> longueur du nom du format            */
 const cs_int_t   *const lopfmt,  /* --> longueur des options du format       */
 const cs_int_t   *const indmod,  /* --> 0 si fige, 1 si deformable,
                                   *     2 si topologie change                */
 const cs_int_t   *const ntchr    /* --> frequence de sortie par defaut       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' eventuels,
                                          Fortran, inutilises lors de
                                          l'appel mais places par de
                                          nombreux compilateurs)              */
);


/*----------------------------------------------------------------------------
 * Creation d'un maillage de post traitement ; les listes de cellules ou
 * faces a extraire sont triees en sortie, qu'elles le soient deja en entree
 * ou non.
 *
 * La liste des cellules associees n'est necessaire que si le nombre
 * de cellules a extraire est strictement superieur a 0 et inferieur au
 * nombre de cellules du maillage.
 *
 * Les listes de faces ne sont prises en compte que si le nombre de cellules
 * a extraire est nul ; si le nombre de faces de bord a extraire est egal au
 * nombre de faces de bord du maillage global, et le nombre de faces internes
 * a extraire est nul, alors on extrait par defaut le maillage de bord, et la
 * liste des faces de bord associees n'est donc pas necessaire.
 *
 * Interface Fortran : utiliser PSTCMA (voir cs_post_util.F)
 *
 * SUBROUTINE PSTCM1 (NUMMAI, NOMMAI, LNMMAI,
 * *****************
 *                    NBRCEL, NBRFAC, NBRFBR, LSTCEL, LSTFAC, LSTFBR)
 *
 * INTEGER          NUMMAI      : --> : Numero du maillage externe a creer
 *                              :     : (< 0 pour maillage standard ou
 *                              :     : developpeur, > 0 pour maillage
 *                              :     : utilisateur)
 * CHARACTER        NOMMAI      : --> : Nom du maillage externe associe
 * INTEGER          LNMMAI      : --> : Longueur du nom de maillage
 * INTEGER          NBRCEL      : --> : Nombre de cellules associees
 * INTEGER          NBRFAC      : --> : Nombre de faces internes associees
 * INTEGER          NBRFBR      : --> : Nombre de faces de bord associees
 * INTEGER          LSTCEL      : <-> : Liste des cellules associees
 * INTEGER          LSTFAC      : <-> : Liste des faces internes associees
 * INTEGER          LSTFBR      : <-> : Liste des faces de bord associees
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcm1, PSTCM1)
(
 const cs_int_t   *const nummai,    /* --> numero du maillage a creer (< 0 pour
                                     *     maillage standard ou developpeur,
                                     *     > 0 pour maillage utilisateur)     */
 const char       *const nommai,    /* --> nom du maillage externe            */
 const cs_int_t   *const lnmmai,    /* --> longueur du nom du maillage        */
 const cs_int_t   *const nbrcel,    /* --> nombre de cellules                 */
 const cs_int_t   *const nbrfac,    /* --> nombre de faces internes           */
 const cs_int_t   *const nbrfbr,    /* --> nombre de faces de bord            */
       cs_int_t          lstcel[],  /* <-> liste des cellules                 */
       cs_int_t          lstfac[],  /* <-> liste des faces internes           */
       cs_int_t          lstfbr[]   /* <-> liste des faces de bord            */
 CS_ARGF_SUPP_CHAINE                /*     (arguments 'longueur' eventuels,
                                           Fortran, inutilises lors de
                                           l'appel mais places par de
                                           nombreux compilateurs)             */
);


/*----------------------------------------------------------------------------
 * Creation d'un alias sur un maillage de post traitement.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTALM (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : --> : Numero de l'alias a creer
 * INTEGER          NUMREF      : --> : Numero du maillage externe associe
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstalm, PSTALM)
(
 const cs_int_t   *nummai,      /* --> numero de l'alias a creer              */
 const cs_int_t   *numref       /* --> numero du maillage associe             */
);


/*----------------------------------------------------------------------------
 * Association d'un "writer" a un maillage pour le post traitement.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTASS (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : --> : Numero du maillage externe associe
 * INTEGER          NUMWRI      : --> : Numero du "writer"
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstass, PSTASS)
(
 const cs_int_t   *nummai,      /* --> numero du maillage externe associe     */
 const cs_int_t   *numwri       /* --> numero du "writer"                     */
);


/*----------------------------------------------------------------------------
 * Mise a jour de l'indicateur "actif" ou "inactif" des "writers" en
 * fonction du pas de temps et de leur frequence de sortie par defaut.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTNTC (NTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : --> : Numero du pas de temps
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t   *ntcabs         /* --> numero de pas de temps associe       */
);


/*----------------------------------------------------------------------------
 * Forcer de l'indicateur "actif" ou "inactif" d'un "writers" specifique
 * ou de l'ensemble des "writers" pour le pas de temps en cours.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTNTC (NTCABS, TTCABS)
 * *****************
 *
 * INTEGER          NUMWRI      : --> : Numero du writer, ou 0 pour forcer
 *                              :     : simultanement tous les writers
 * INTEGER          INDACT      : --> : 0 pour desactiver, 1 pour activer
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstact, PSTACT)
(
 const cs_int_t   *numwri,     /* --> numero du writer, ou 0 pour forcer
                                *     simultanement tous les writers          */
 const cs_int_t   *indact      /* --> 0 pour desactiver, 1 pour activer       */
);


/*----------------------------------------------------------------------------
 * Ecriture des maillages de post traitement en fonction des writers
 * associes.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTEMA (NTCABS, TTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : --> : Numero du pas de temps
 * DOUBLE PRECISION TTCABS      : --> : Temps physique associe
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstema, PSTEMA)
(
 const cs_int_t   *ntcabs,        /* --> numero de pas de temps associe       */
 const cs_real_t  *ttcabs         /* --> valeur du pas de temps associe       */
);


/*----------------------------------------------------------------------------
 * Boucle sur les maillages de post traitement pour ecriture  des variables
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *const idbia0,      /* --> numero 1ere case libre dans IA   */
 const cs_int_t   *const idbra0,      /* --> numero 1ere case libre dans RA   */
 const cs_int_t   *const ndim,        /* --> dimension de l'espace            */
 const cs_int_t   *const ntcabs,      /* --> numero de pas de temps courant   */
 const cs_int_t   *const ncelet,      /* --> nombre de cellules etendu        */
 const cs_int_t   *const ncel,        /* --> nombre de cellules               */
 const cs_int_t   *const nfac,        /* --> nombre de faces internes         */
 const cs_int_t   *const nfabor,      /* --> nombre de faces de bord          */
 const cs_int_t   *const nfml,        /* --> nombre de familles               */
 const cs_int_t   *const nprfml,      /* --> nombre de proprietes des familles*/
 const cs_int_t   *const nnod,        /* --> nombre de noeuds                 */
 const cs_int_t   *const lndfac,      /* --> longueur de nodfac               */
 const cs_int_t   *const lndfbr,      /* --> longueur de nodfbr               */
 const cs_int_t   *const ncelbr,      /* --> nombre de cellules de bord       */
 const cs_int_t   *const nvar,        /* --> nombre de variables              */
 const cs_int_t   *const nscal,       /* --> nombre de scalaires              */
 const cs_int_t   *const nphas,       /* --> nombre de phases                 */
 const cs_int_t   *const nvlsta,      /* --> nombre de variables stat. (lagr) */
 const cs_int_t   *const nvisbr,      /* --> nombre de variables stat. (lagr) */
 const cs_int_t   *const nideve,      /* --> longueur du tableau idevel[]     */
 const cs_int_t   *const nrdeve,      /* --> longueur du tableau rdevel[]     */
 const cs_int_t   *const nituse,      /* --> longueur du tableau ituser[]     */
 const cs_int_t   *const nrtuse,      /* --> longueur du tableau rtuser[]     */
 const cs_int_t          ifacel[],    /* --> liste des faces internes         */
 const cs_int_t          ifabor[],    /* --> liste des faces de bord          */
 const cs_int_t          ifmfbr[],    /* --> liste des familles des faces bord*/
 const cs_int_t          ifmcel[],    /* --> liste des familles des cellules  */
 const cs_int_t          iprfml[],    /* --> liste des proprietes des familles*/
 const cs_int_t          ipnfac[],    /* --> rg ds nodfac 1er sommet faces int*/
 const cs_int_t          nodfac[],    /* --> numero des sommets des faces int.*/
 const cs_int_t          ipnfbr[],    /* --> rg ds nodfbr 1er sommet faces brd*/
 const cs_int_t          nodfbr[],    /* --> numero des sommets des faces bord*/
 const cs_int_t          idevel[],    /* --> tab. complementaire developpeur  */
 const cs_int_t          ituser[],    /* --> tab. complementaire utilisateur  */
 const cs_int_t          ia[],        /* --> macro-tableau entier             */
 const cs_real_t  *const ttcabs,      /* --> temps courant absolu             */
 const cs_real_t         xyzcen[],    /* --> c.d.g. des cellules              */
 const cs_real_t         surfac[],    /* --> surfaces des faces internes      */
 const cs_real_t         surfbo[],    /* --> surfaces des faces de bord       */
 const cs_real_t         cdgfac[],    /* --> c.d.g. des faces internes        */
 const cs_real_t         cdgfbo[],    /* --> c.d.g. des faces de bord         */
 const cs_real_t         xyznod[],    /* --> coordonnees des sommets          */
 const cs_real_t         volume[],    /* --> volumes des cellules             */
 const cs_real_t         dt[],        /* --> pas de temps                     */
 const cs_real_t         rtpa[],      /* --> variables aux cellules (prec.)   */
 const cs_real_t         rtp[],       /* --> variables aux cellules           */
 const cs_real_t         propce[],    /* --> proprietes physiques cellules    */
 const cs_real_t         propfa[],    /* --> proprietes physiques aux faces   */
 const cs_real_t         propfb[],    /* --> proprietes physiques faces bord  */
 const cs_real_t         coefa[],     /* --> cond. limites aux faces de bord  */
 const cs_real_t         coefb[],     /* --> cond. limites aux faces de bord  */
 const cs_real_t         statce[],    /* --> moyennes statistiques (Lagrangien*/
 const cs_real_t         stativ[],    /* --> variances statistiques (Lagrangie*/
 const cs_real_t         statfb[],    /* --> moyennes statistiques (Lagrangien*/
 const cs_real_t         rdevel[],    /* --> tab. complementaire developpeur  */
 const cs_real_t         rtuser[],    /* --> tab. complementaire utilisateur  */
 const cs_real_t         ra[]         /* --> macro-tableau reel               */
);


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement defini sur les cellules ou faces
 * d'un maillage en fonction des "writers" associes.
 *
 * Interface Fortran : utiliser PSTEVA (voir cs_post_util.F)
 *
 * SUBROUTINE PSTEVA (NUMMAI, NOMVAR, IDIMT,  IENTLA, IVARPR,
 * *****************
 *                    NTCABS, TTCABS, VARCEL, VARFAC, VARFBR)
 *
 * INTEGER          NUMMAI      : --> : Numero du maillage associe
 * CHARACTER        NOMVAR      : --> : Nom de la variable
 * INTEGER          IDIMT       : --> : 1 pour scalaire, 3 pour vecteur
 * INTEGER          IENTLA      : --> : Si vecteur, 1 si valeurs entrelacees
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 sinon (x1, x2, ...xn, y1, y2, ...)
 * INTEGER          IVARPR      : --> : 1 si variable definie sur maillage
 *                              :     : "parent", 2 si variable restreinte
 *                              :     : au maillage post
 * INTEGER          NTCABS      : --> : Numero du pas de temps
 * DOUBLE PRECISION TTCABS      : --> : Temps physique associe
 * DOUBLE PRECISION VARCEL(*)   : --> : Valeurs associees aux cellules
 * DOUBLE PRECISION VARFAC(*)   : --> : Valeurs associees aux faces internes
 * DOUBLE PRECISION VARFBO(*)   : --> : Valeurs associees aux faces de bord
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstev1, PSTEV1)
(
 const cs_int_t   *const nummai,      /* --> numero du maillage associe       */
 const char       *const nomvar,      /* --> nom de la variable               */
 const cs_int_t   *const lnmvar,      /* --> longueur du nom de la variable   */
 const cs_int_t   *const idimt,       /* --> 1 pour scalaire, 3 pour vecteur  */
 const cs_int_t   *const ientla,      /* --> si vecteur, 1 si valeurs
                                       *     entrelacees, 0 sinon             */
 const cs_int_t   *const ivarpr,      /* --> 1 si variable definie sur
                                       *     maillage "parent", 2 si variable
                                       *     restreinte au maillage post      */
 const cs_int_t   *const ntcabs,      /* --> numero de pas de temps associe   */
 const cs_real_t  *const ttcabs,      /* --> valeur du pas de temps associe   */
 const cs_real_t         varcel[],    /* --> valeurs aux cellules             */
 const cs_real_t         varfac[],    /* --> valeurs aux faces internes       */
 const cs_real_t         varfbr[]     /* --> valeurs aux faces de bord        */
 CS_ARGF_SUPP_CHAINE                  /*     (arguments 'longueur' eventuels,
                                             Fortran, inutilises lors de
                                             l'appel mais places par de
                                             nombreux compilateurs)           */
);


/*----------------------------------------------------------------------------
 * Prise en compte de la renumerotation des faces et faces de bord
 * dans les liens de "parente" des maillages post.
 *
 * Cette fonction ne doit etre appellee qu'une fois, apres la renumerotation
 * evuentuelle des faces, pour adapter les maillages post existants.
 * Des nouveaux maillages post seront automatiquement bases sur la
 * "bonne" numerotation, par construction.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTRNM(IVECTI, IVECTB, INUMFI, INUMFB)
 * *****************
 *
 * INTEGER IVECTI               : --> : Indicateur de renum. faces internes
 * INTEGER IVECTB               : --> : Indicateur de renum. faces de bord
 * INTEGER INUMFI(NFAC)         : --> : Table de renum. des faces internes
 * INTEGER INUMFB(NFABOR)       : --> : Table de renum. des faces de bord
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstrnm, PSTRNM)
(
 cs_int_t  *ivecti,           /* --> vectorisation des faces internes         */
 cs_int_t  *ivectb,           /* --> vectorisation des faces de bord          */
 cs_int_t  *inumfi,           /* --> numerotation initiale des faces internes */
 cs_int_t  *inumfb            /* --> numerotation initiale des faces de bord  */
);


/*============================================================================
 *  Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation d'un "writer" ; cet objet correspond au choix d'un nom de cas,
 * de repertoire, et de format, ainsi qu'un indicateur precisant si les
 * maillages associes doivent dependre ou non du temps, et la frequence de
 * sortie par defaut pour les variables associees.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_writer
(
       cs_int_t          id_writer,  /* --> numero du writer a creer
                                      *     (< 0 pour writer reserve,
                                      *      > 0 pour writer utilisateur)     */
 const char       *const nom_cas,    /* --> nom du cas associe                */
 const char       *const nom_rep,    /* --> nom de repertoire associe         */
 const char       *const nom_fmt,    /* --> nom de format associe             */
 const char       *const opt_fmt,    /* --> options associees au format       */
       cs_int_t          ind_mod,    /* --> 0 si fige, 1 si deformable,
                                      *     2 si topologie change, +10 pour
                                      *     ajouter un champ deplacement      */
       cs_int_t          frequence   /* --> frequence de sortie par defaut    */
);


/*----------------------------------------------------------------------------
 * Creation d'un maillage de post traitement ; les listes de cellules ou
 * faces a extraire sont triees en sortie, qu'elles le soient deja en entree
 * ou non.
 *
 * La liste des cellules associees n'est necessaire que si le nombre
 * de cellules a extraire est strictement superieur a 0 et inferieur au
 * nombre de cellules du maillage.
 *
 * Les listes de faces ne sont prises en compte que si le nombre de cellules
 * a extraire est nul ; si le nombre de faces de bord a extraire est egal au
 * nombre de faces de bord du maillage global, et le nombre de faces internes
 * a extraire est nul, alors on extrait par defaut le maillage de bord, et la
 * liste des faces de bord associees n'est donc pas necessaire.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_maillage
(
 const cs_int_t          id_maillage,  /* --> numero du maillage a creer
                                        *     (< 0 pour maillage reserve,
                                        *      > 0 pour maillage utilisateur) */
 const char       *const nom_maillage, /* --> nom du maillage externe         */
 const cs_int_t          nbr_cel,      /* --> nombre de cellules              */
 const cs_int_t          nbr_fac,      /* --> nombre de faces internes        */
 const cs_int_t          nbr_fbr,      /* --> nombre de faces de bord         */
       cs_int_t          liste_cel[],  /* <-> liste des cellules              */
       cs_int_t          liste_fac[],  /* <-> liste des faces internes        */
       cs_int_t          liste_fbr[]   /* <-> liste des faces de bord         */
);


/*----------------------------------------------------------------------------
 * Creation d'un maillage de post traitement par association d'un maillage
 * externe existant.
 *
 * Si le maillage externe n'est plus destine a etre utilise par ailleurs,
 * on peut choisir d'en transferer la propriete au maillage de post traitement,
 * qui gerera alors son cycle de vie selon ses seuls besoins.
 *
 * Si le maillage externe doit continuer a etre partage, on devra veiller
 * a maintenir la coherence entre ce maillage et le posttraitement au cours
 * du temps.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_maillage_existant
(
 cs_int_t            id_maillage,      /* --> numero du maillage a creer
                                        *     (< 0 pour maillage reserve,
                                        *      > 0 pour maillage utilisateur) */
 fvm_nodal_t  *const maillage_ext,     /* --> maillage externe */
 cs_bool_t           transferer        /* --> indique si l'on transfere la
                                        *     propriete du maillage externe
                                              au maillage de post traitement  */
);


/*----------------------------------------------------------------------------
 * Creation d'un alias sur un maillage de post traitement.
 *
 * Un alias permet d'associer un numero supplementaire a un maillage de
 * post traitement deja defini, et donc de lui associer d'autres
 * "writers" qu'au maillage initial ; ceci permet par exemple d'ecrire
 * un jeu de variables principales tous les n1 pas de temps dans un
 * jeu de donnees de post traitement, et de sortir quelques variables
 * specifiques tous les n2 pas de temps dans un autre jeu de donnees
 * de post traitement, sans necessiter de duplication du maillage support.
 *
 * Un alias est donc traite en tout point comme le maillage principal
 * associe ; en particulier, si la definition de l'un est modifie, celle
 * de l'autre l'est aussi.
 *
 * Il est impossible d'associer un alias a un autre alias (cela n'aurait
 * pas d'utilite), mais on peut associer plusieurs alias a un maillage.
 *----------------------------------------------------------------------------*/

void cs_post_alias_maillage
(
 const cs_int_t          id_alias,     /* --> numero de l'alias a creer
                                        *     (< 0 pour alias reserve,
                                        *      > 0 pour alias utilisateur)    */
 const cs_int_t          id_maillage   /* --> numero du maillage  associe     */
);


/*----------------------------------------------------------------------------
 * Verifie l'existence d'un "writer" associe a un numero donne.
 *----------------------------------------------------------------------------*/

cs_bool_t cs_post_existe_writer
(
 const cs_int_t   numwri        /* --> numero du writer associe               */
);


/*----------------------------------------------------------------------------
 * Get the writer associated to a writer_id.
 *
 * writer_id       -->  id of the writer in cs_glob_post_writers
 *
 * Returns:
 *  a pointer to a fvm_writer_t structure
 *----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(cs_int_t   writer_id);

/*----------------------------------------------------------------------------
 * Verifie l'existence d'un maillage de post traitement associe a un
 * numero donne.
 *----------------------------------------------------------------------------*/

cs_bool_t cs_post_existe_maillage
(
 const cs_int_t   nummai        /* --> numero du maillage externe associe     */
);


/*----------------------------------------------------------------------------
 * Modification d'un maillage de post traitement existant.
 *
 * Il s'agit ici de modifier les listes de cellules ou faces du maillage,
 * par exemple pour faire evoluer une coupe en fonction des zones
 * "interessantes (il n'est pas necessaire de recourir a cette fonction
 * si le maillage se deforme simplement).
 *----------------------------------------------------------------------------*/

void cs_post_modifie_maillage
(
 const cs_int_t          id_maillage,  /* --> numero du writer a creer
                                        *     (< 0 pour maillage reserve,
                                        *      > 0 pour maillage utilisateur) */
 const cs_int_t          nbr_cel,      /* --> nombre de cellules              */
 const cs_int_t          nbr_fac,      /* --> nombre de faces internes        */
 const cs_int_t          nbr_fbr,      /* --> nombre de faces de bord         */
       cs_int_t          liste_cel[],  /* <-> liste des cellules              */
       cs_int_t          liste_fac[],  /* <-> liste des faces internes        */
       cs_int_t          liste_fbr[]   /* <-> liste des faces de bord         */
);


/*----------------------------------------------------------------------------
 * Recuperation du prochain numero de maillage standard ou developpeur
 * disponible (base sur le plus petit numero negatif present -1).
 *----------------------------------------------------------------------------*/

cs_int_t cs_post_ret_num_maillage_libre
(
 void
);


/*----------------------------------------------------------------------------
 * Association d'un "writer" a un maillage pour le post traitement.
 *----------------------------------------------------------------------------*/

void cs_post_associe
(
 const cs_int_t   id_maillage,  /* --> numero du maillage externe associe     */
 const cs_int_t   id_writer     /* --> numero du writer                       */
);


/*----------------------------------------------------------------------------
 * Mise a jour de l'indicateur "actif" ou "inactif" des "writers" en
 * fonction du pas de temps et de leur frequence de sortie par defaut.
 *----------------------------------------------------------------------------*/

void cs_post_activer_selon_defaut
(
 const cs_int_t   nt_cur_abs    /* --> numero de pas de temps courant         */
);


/*----------------------------------------------------------------------------
 * Forcer de l'indicateur "actif" ou "inactif" d'un "writers" specifique
 * ou de l'ensemble des "writers" pour le pas de temps en cours.
 *----------------------------------------------------------------------------*/

void cs_post_activer_writer
(
 const cs_int_t   id_writer,    /* --> numero du writer,ou 0 pour forcer
                                 *     simultanement tous les writers         */
 const cs_int_t   activer       /* --> 0 pour desactiver, 1 pour activer      */
);


/*----------------------------------------------------------------------------
 * Ecriture des maillages de post traitement en fonction des "writers"
 * associes.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_maillages
(
 const cs_int_t   nt_cur_abs,         /* --> numero de pas de temps courant   */
 const cs_real_t  t_cur_abs           /* --> valeur du temps physique associe */
);


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement defini sur les cellules ou faces
 * d'un maillage en fonction des "writers" associes.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_var
(
       cs_int_t          id_maillage,  /* --> numero du maillage post associe */
 const char             *nom_var,      /* --> nom de la variable              */
       cs_int_t          dim_var,      /* --> 1 pour scalaire, 3 pour vecteur */
       cs_bool_t         entrelace,    /* --> si vecteur, vrai si valeurs
                                        *     entrelacees, faux sinon         */
       cs_bool_t         var_parent,   /* --> vrai si valeurs definies sur
                                        *     maillage "parent", faux si
                                        *     restreintes au maillage post    */
       cs_post_type_t    var_type,     /* --> type de donnees associe         */
       cs_int_t          nt_cur_abs,   /* --> numero de pas de temps courant  */
       cs_real_t         t_cur_abs,    /* --> valeur du temps physique        */
 const void             *var_cel,      /* --> valeurs aux cellules            */
 const void             *var_fac,      /* --> valeurs aux faces internes      */
 const void             *var_fbr       /* --> valeurs aux faces de bord       */
);


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement defini sur les sommets
 * d'un maillage en fonction des "writers" associes.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_var_som
(
       cs_int_t          id_maillage,  /* --> numero du maillage post associe */
 const char             *nom_var,      /* --> nom de la variable              */
       cs_int_t          dim_var,      /* --> 1 pour scalaire, 3 pour vecteur */
       cs_bool_t         entrelace,    /* --> si vecteur, vrai si valeurs
                                        *     entrelacees, faux sinon         */
       cs_bool_t         var_parent,   /* --> vrai si valeurs definies sur
                                        *     maillage "parent", faux si
                                        *     restreintes au maillage post    */
       cs_post_type_t    var_type,     /* --> type de donnees associe         */
       cs_int_t          nt_cur_abs,   /* --> numero de pas de temps courant  */
       cs_real_t         t_cur_abs,    /* --> valeur du temps physique        */
 const void             *var_som       /* --> valeurs aux sommets             */
);


/*----------------------------------------------------------------------------
 * Prise en compte de la renumerotation des faces et faces de bord
 * dans les liens de "parente" des maillages post.
 *
 * Cette fonction ne doit etre appellee qu'une fois, apres la renumerotation
 * evuentuelle des faces, pour adapter les maillages post existants.
 * Des nouveaux maillages post seront automatiquement bases sur la
 * "bonne" numerotation, par construction.
 *----------------------------------------------------------------------------*/

void cs_post_renum_faces
(
 cs_int_t  *init_i_face_num,  /* --> numerotation initiale des faces internes */
 cs_int_t  *init_b_face_num   /* --> numerotation initiale des faces de bord  */
);

/*----------------------------------------------------------------------------
 * Destruction des structures associees aux post traitements
 *----------------------------------------------------------------------------*/

void cs_post_detruit
(
 void
);


/*----------------------------------------------------------------------------
 * Initialisation du post-traitement principal
 *----------------------------------------------------------------------------*/

void cs_post_init_pcp
(
 void
);


/*----------------------------------------------------------------------------
 * Ajout d'un traitement de variable temporelle a l'appel de PSTVAR.
 *
 * L'identificateur d'instance associe a la fonction permet d'ajouter
 * une meme fonction plusieurs fois, avec un identificateur different
 * permettant a la fonction de selectionner un sous-traitement.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_var_temporelle
(
 cs_post_var_temporelle_t  *fonction,    /* Fonction associee                 */
 cs_int_t                   id_instance  /* Indentificateur d'instance
                                            associe a la fonction             */
 );


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_POST_H__ */
