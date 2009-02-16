/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Post-processing management
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Définitions d'énumerations
 *============================================================================*/

/* Énumération pour transmettre le type d'une donnée */

typedef enum {
  CS_POST_TYPE_cs_int_t,
  CS_POST_TYPE_cs_real_t,
  CS_POST_TYPE_int,
  CS_POST_TYPE_float,
  CS_POST_TYPE_double
} cs_post_type_t;


/*============================================================================
 * Définition de macros
 *============================================================================*/


/*============================================================================
 * Déclaration de structures et types
 *============================================================================*/

/* Pointeur associé à un "writer" : cet objet correspond au choix d'un
 * nom de cas, de répertoire, et de format, ainsi qu'un indicateur précisant
 * si les maillages associés doivent dépendre ou non du temps, et la
 * fréquence de sortie par défaut pour les variables associées. */

typedef struct _cs_post_writer_t cs_post_writer_t;

/* Pointeur associé à un maillage de post traitement ; cet objet
 * gère le lien entre un tel maillage et les "writers" associés. */

typedef struct _cs_post_maillage_t cs_post_maillage_t;

/* Pointeur de fonction associé à un post-traitement particulier ;
 * on enregistre de telles fonctions via la fonction
 * cs_post_ajoute_var_temporelle(), et toutes les fonctions enregistrées
 * de la sorte sont appellées automatiquement par PSTVAR. */

typedef void
(cs_post_var_temporelle_t) (int          id_instance,
                            int          nt_cur_abs,
                            cs_real_t    t_cur_abs);


/*============================================================================
 * Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Création d'un "writer" à partir des données du Fortran ; cet objet
 * correspond au choix d'un nom de cas, de répertoire, et de format, ainsi
 * qu'un indicateur précisant si les maillages associés doivent dépendre ou
 * non du temps, et la fréquence de sortie par défaut pour les
 * variables associées.
 *
 * Interface Fortran : utiliser PSTCWR (voir cs_post_util.F)
 *
 * SUBROUTINE PSTCWR (NUMGEP, NOMCAS, NOMREP, NOMFMT, OPTFMT,
 * *****************
 *                    LNMCAS, LNMFMT, LNMREP, LOPFMT,
 *                    INDMOD, NTCHR)
 *
 * INTEGER          NUMGEP      : --> : Numéro du filtre à créer (< 0 pour
 *                              :     : filtre standard ou développeur,
 *                              :     : > 0 pour filtre utilisateur)
 * CHARACTER        NOMCAS      : --> : Nom du cas associé
 * CHARACTER        NOMREP      : --> : Nom du répertoire associé
 * INTEGER          NOMFMT      : --> : Nom de format associé
 * INTEGER          OPTFMT      : --> : Options associées au format
 * INTEGER          LNMCAS      : --> : Longueur du nom du cas
 * INTEGER          LNMREP      : --> : Longueur du nom du répertoire
 * INTEGER          LNMFMT      : --> : Longueur du nom du format
 * INTEGER          LOPFMT      : --> : Longueur des options du format
 * INTEGER          INDMOD      : --> : 0 si figé, 1 si déformable,
 *                              :     : 2 si la topologie change
 * INTEGER          NTCHR       : --> : Fréquence de sortie par défaut
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcw1, PSTCW1)
(
 const cs_int_t   *const numwri,  /* --> numéro du writer à créer
                                   *     < 0 pour writer réservé,
                                   *     > 0 pour writer utilisateur)         */
 const char       *const nomcas,  /* --> nom du cas associé                   */
 const char       *const nomrep,  /* --> nom de répertoire associé            */
 const char       *const nomfmt,  /* --> nom de format associé                */
 const char       *const optfmt,  /* --> options associées au format          */
 const cs_int_t   *const lnmcas,  /* --> longueur du nom du cas               */
 const cs_int_t   *const lnmrep,  /* --> longueur du nom du répertoire        */
 const cs_int_t   *const lnmfmt,  /* --> longueur du nom du format            */
 const cs_int_t   *const lopfmt,  /* --> longueur des options du format       */
 const cs_int_t   *const indmod,  /* --> 0 si figé, 1 si déformable,
                                   *     2 si topologie change                */
 const cs_int_t   *const ntchr    /* --> fréquence de sortie par défaut       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' éventuels,
                                          Fortran, inutilisés lors de
                                          l'appel mais placés par de
                                          nombreux compilateurs)              */
);


/*----------------------------------------------------------------------------
 * Création d'un maillage de post traitement ; les listes de cellules ou
 * faces à extraire sont triées en sortie, qu'elles le soient déjà en entrée
 * ou non.
 *
 * La liste des cellules associées n'est nécessaire que si le nombre
 * de cellules à extraire est strictement supérieur à 0 et inférieur au
 * nombre de cellules du maillage.
 *
 * Les listes de faces ne sont prises en compte que si le nombre de cellules
 * à extraire est nul ; si le nombre de faces de bord à extraire est égal au
 * nombre de faces de bord du maillage global, et le nombre de faces internes
 * à extraire est nul, alors on extrait par défaut le maillage de bord, et la
 * liste des faces de bord associées n'est donc pas nécessaire.
 *
 * Interface Fortran : utiliser PSTCMA (voir cs_post_util.F)
 *
 * SUBROUTINE PSTCM1 (NUMMAI, NOMMAI, LNMMAI,
 * *****************
 *                    NBRCEL, NBRFAC, NBRFBR, LSTCEL, LSTFAC, LSTFBR)
 *
 * INTEGER          NUMMAI      : --> : Numéro du maillage externe à créer
 *                              :     : (< 0 pour maillage standard ou
 *                              :     : développeur, > 0 pour maillage
 *                              :     : utilisateur)
 * CHARACTER        NOMMAI      : --> : Nom du maillage externe associé
 * INTEGER          LNMMAI      : --> : Longueur du nom de maillage
 * INTEGER          NBRCEL      : --> : Nombre de cellules associées
 * INTEGER          NBRFAC      : --> : Nombre de faces internes associées
 * INTEGER          NBRFBR      : --> : Nombre de faces de bord associées
 * INTEGER          LSTCEL      : <-> : Liste des cellules associées
 * INTEGER          LSTFAC      : <-> : Liste des faces internes associées
 * INTEGER          LSTFBR      : <-> : Liste des faces de bord associées
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcm1, PSTCM1)
(
 const cs_int_t   *const nummai,    /* --> numéro du maillage à créer (< 0 pour
                                     *     maillage standard ou développeur,
                                     *     > 0 pour maillage utilisateur)     */
 const char       *const nommai,    /* --> nom du maillage externe            */
 const cs_int_t   *const lnmmai,    /* --> longueur du nom du maillage        */
 const cs_int_t   *const nbrcel,    /* --> nombre de cellules                 */
 const cs_int_t   *const nbrfac,    /* --> nombre de faces internes           */
 const cs_int_t   *const nbrfbr,    /* --> nombre de faces de bord            */
       cs_int_t          lstcel[],  /* <-> liste des cellules                 */
       cs_int_t          lstfac[],  /* <-> liste des faces internes           */
       cs_int_t          lstfbr[]   /* <-> liste des faces de bord            */
 CS_ARGF_SUPP_CHAINE                /*     (arguments 'longueur' éventuels,
                                           Fortran, inutilisés lors de
                                           l'appel mais placés par de
                                           nombreux compilateurs)             */
);


/*----------------------------------------------------------------------------
 * Création d'un alias sur un maillage de post traitement.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTALM (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : --> : Numéro de l'alias à créer
 * INTEGER          NUMREF      : --> : Numéro du maillage externe associé
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstalm, PSTALM)
(
 const cs_int_t   *nummai,      /* --> numéro de l'alias à créer              */
 const cs_int_t   *numref       /* --> numéro du maillage associe             */
);


/*----------------------------------------------------------------------------
 * Association d'un "writer" à un maillage pour le post traitement.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTASS (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : --> : Numéro du maillage externe associé
 * INTEGER          NUMWRI      : --> : Numéro du "writer"
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstass, PSTASS)
(
 const cs_int_t   *nummai,      /* --> numéro du maillage externe associé     */
 const cs_int_t   *numwri       /* --> numéro du "writer"                     */
);


/*----------------------------------------------------------------------------
 * Mise à jour de l'indicateur "actif" ou "inactif" des "writers" en
 * fonction du pas de temps et de leur fréquence de sortie par défaut.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTNTC (NTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : --> : Numéro du pas de temps
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t   *ntcabs         /* --> numéro de pas de temps associé       */
);


/*----------------------------------------------------------------------------
 * Forcer de l'indicateur "actif" ou "inactif" d'un "writers" spécifique
 * ou de l'ensemble des "writers" pour le pas de temps en cours.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTNTC (NTCABS, TTCABS)
 * *****************
 *
 * INTEGER          NUMWRI      : --> : Numéro du writer, ou 0 pour forcer
 *                              :     : simultanément tous les writers
 * INTEGER          INDACT      : --> : 0 pour désactiver, 1 pour activer
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstact, PSTACT)
(
 const cs_int_t   *numwri,     /* --> numéro du writer, ou 0 pour forcer
                                *     simultanément tous les writers          */
 const cs_int_t   *indact      /* --> 0 pour désactiver, 1 pour activer       */
);


/*----------------------------------------------------------------------------
 * Ecriture des maillages de post traitement en fonction des writers
 * associés.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTEMA (NTCABS, TTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : --> : Numéro du pas de temps
 * DOUBLE PRECISION TTCABS      : --> : Temps physique associé
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstema, PSTEMA)
(
 const cs_int_t   *ntcabs,        /* --> numéro de pas de temps associé       */
 const cs_real_t  *ttcabs         /* --> valeur du pas de temps associé       */
);


/*----------------------------------------------------------------------------
 * Boucle sur les maillages de post traitement pour écriture  des variables
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *const idbia0,      /* --> numéro 1ère case libre dans IA   */
 const cs_int_t   *const idbra0,      /* --> numéro 1ère case libre dans RA   */
 const cs_int_t   *const ndim,        /* --> dimension de l'espace            */
 const cs_int_t   *const ntcabs,      /* --> numéro de pas de temps courant   */
 const cs_int_t   *const ncelet,      /* --> nombre de cellules étendu        */
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
 const cs_int_t          nodfbr[],    /* --> numéro des sommets des faces bord*/
 const cs_int_t          idevel[],    /* --> tab. complémentaire développeur  */
 const cs_int_t          ituser[],    /* --> tab. complémentaire utilisateur  */
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
 const cs_real_t         rtpa[],      /* --> variables aux cellules (préc.)   */
 const cs_real_t         rtp[],       /* --> variables aux cellules           */
 const cs_real_t         propce[],    /* --> propriétés physiques cellules    */
 const cs_real_t         propfa[],    /* --> propriétés physiques aux faces   */
 const cs_real_t         propfb[],    /* --> propriétés physiques faces bord  */
 const cs_real_t         coefa[],     /* --> cond. limites aux faces de bord  */
 const cs_real_t         coefb[],     /* --> cond. limites aux faces de bord  */
 const cs_real_t         statce[],    /* --> moyennes statistiques (Lagrangien*/
 const cs_real_t         stativ[],    /* --> variances statistiques (Lagrangie*/
 const cs_real_t         statfb[],    /* --> moyennes statistiques (Lagrangien*/
 const cs_real_t         rdevel[],    /* --> tab. complémentaire développeur  */
 const cs_real_t         rtuser[],    /* --> tab. complémentaire utilisateur  */
 const cs_real_t         ra[]         /* --> macro-tableau réel               */
);


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement défini sur les cellules ou faces
 * d'un maillage en fonction des "writers" associés.
 *
 * Interface Fortran : utiliser PSTEVA (voir cs_post_util.F)
 *
 * SUBROUTINE PSTEVA (NUMMAI, NOMVAR, IDIMT,  IENTLA, IVARPR,
 * *****************
 *                    NTCABS, TTCABS, VARCEL, VARFAC, VARFBR)
 *
 * INTEGER          NUMMAI      : --> : Numéro du maillage associé
 * CHARACTER        NOMVAR      : --> : Nom de la variable
 * INTEGER          IDIMT       : --> : 1 pour scalaire, 3 pour vecteur
 * INTEGER          IENTLA      : --> : Si vecteur, 1 si valeurs entrelacées
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 sinon (x1, x2, ...xn, y1, y2, ...)
 * INTEGER          IVARPR      : --> : 1 si variable définie sur maillage
 *                              :     : "parent", 2 si variable restreinte
 *                              :     : au maillage post
 * INTEGER          NTCABS      : --> : Numéro du pas de temps
 * DOUBLE PRECISION TTCABS      : --> : Temps physique associé
 * DOUBLE PRECISION VARCEL(*)   : --> : Valeurs associées aux cellules
 * DOUBLE PRECISION VARFAC(*)   : --> : Valeurs associées aux faces internes
 * DOUBLE PRECISION VARFBO(*)   : --> : Valeurs associées aux faces de bord
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstev1, PSTEV1)
(
 const cs_int_t   *const nummai,      /* --> numéro du maillage associé       */
 const char       *const nomvar,      /* --> nom de la variable               */
 const cs_int_t   *const lnmvar,      /* --> longueur du nom de la variable   */
 const cs_int_t   *const idimt,       /* --> 1 pour scalaire, 3 pour vecteur  */
 const cs_int_t   *const ientla,      /* --> si vecteur, 1 si valeurs
                                       *     entrelacées, 0 sinon             */
 const cs_int_t   *const ivarpr,      /* --> 1 si variable définie sur
                                       *     maillage "parent", 2 si variable
                                       *     restreinte au maillage post      */
 const cs_int_t   *const ntcabs,      /* --> numéro de pas de temps associé   */
 const cs_real_t  *const ttcabs,      /* --> valeur du pas de temps associé   */
 const cs_real_t         varcel[],    /* --> valeurs aux cellules             */
 const cs_real_t         varfac[],    /* --> valeurs aux faces internes       */
 const cs_real_t         varfbr[]     /* --> valeurs aux faces de bord        */
 CS_ARGF_SUPP_CHAINE                  /*     (arguments 'longueur' éventuels,
                                             Fortran, inutilisés lors de
                                             l'appel mais placés par de
                                             nombreux compilateurs)           */
);


/*============================================================================
 *  Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Création d'un "writer" ; cet objet correspond au choix d'un nom de cas,
 * de répertoire, et de format, ainsi qu'un indicateur précisant si les
 * maillages associés doivent dépendre ou non du temps, et la fréquence de
 * sortie par défaut pour les variables associées.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_writer
(
       cs_int_t          id_writer,  /* --> numéro du writer à créer
                                      *     (< 0 pour writer réservé,
                                      *      > 0 pour writer utilisateur)     */
 const char       *const nom_cas,    /* --> nom du cas associé                */
 const char       *const nom_rep,    /* --> nom de répertoire associé         */
 const char       *const nom_fmt,    /* --> nom de format associé             */
 const char       *const opt_fmt,    /* --> options associées au format       */
       cs_int_t          ind_mod,    /* --> 0 si figé, 1 si déformable,
                                      *     2 si topologie change, +10 pour
                                      *     ajouter un champ déplacement      */
       cs_int_t          frequence   /* --> fréquence de sortie par défaut    */
);


/*----------------------------------------------------------------------------
 * Création d'un maillage de post traitement ; les listes de cellules ou
 * faces à extraire sont triées en sortie, qu'elles le soient déjà en entrée
 * ou non.
 *
 * La liste des cellules associées n'est nécessaire que si le nombre
 * de cellules à extraire est strictement supérieur à 0 et inférieur au
 * nombre de cellules du maillage.
 *
 * Les listes de faces ne sont prises en compte que si le nombre de cellules
 * à extraire est nul ; si le nombre de faces de bord à extraire est égal au
 * nombre de faces de bord du maillage global, et le nombre de faces internes
 * à extraire est nul, alors on extrait par défaut le maillage de bord, et la
 * liste des faces de bord associées n'est donc pas nécessaire.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_maillage
(
 const cs_int_t          id_maillage,  /* --> numéro du maillage à créer
                                        *     (< 0 pour maillage réservé,
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
 * Création d'un maillage de post traitement par association d'un maillage
 * externe existant.
 *
 * Si le maillage externe n'est plus destiné à être utilisé par ailleurs,
 * on peut choisir d'en transférer la propriété au maillage de post traitement,
 * qui gèrera alors son cycle de vie selon ses seuls besoins.
 *
 * Si le maillage externe doit continuer à être partagé, on devra veiller
 * à maintenir la cohérence entre ce maillage et le posttraitement au cours
 * du temps.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_maillage_existant
(
 cs_int_t            id_maillage,      /* --> numéro du maillage à créer
                                        *     (< 0 pour maillage réservé,
                                        *      > 0 pour maillage utilisateur) */
 fvm_nodal_t  *const maillage_ext,     /* --> maillage externe */
 cs_bool_t           transferer        /* --> indique si l'on transfère la
                                        *     propriété du maillage externe
                                              au maillage de post traitement  */
);


/*----------------------------------------------------------------------------
 * Création d'un alias sur un maillage de post traitement.
 *
 * Un alias permet d'associer un numéro supplémentaire à un maillage de
 * post traitement déjà défini, et donc de lui associer d'autres
 * "writers" qu'au maillage initial ; ceci permet par exemple d'écrire
 * un jeu de variables principales tous les n1 pas de temps dans un
 * jeu de données de post traitement, et de sortir quelques variables
 * spécifiques tous les n2 pas de temps dans un autre jeu de données
 * de post traitement, sans nécessiter de duplication du maillage support.
 *
 * Un alias est donc traité en tout point comme le maillage principal
 * associé ; en particulier, si la définition de l'un est modifié, celle
 * de l'autre l'est aussi.
 *
 * Il est impossible d'associer un alias à un autre alias (cela n'aurait
 * pas d'utilité), mais on peut associer plusieurs alias à un maillage.
 *----------------------------------------------------------------------------*/

void cs_post_alias_maillage
(
 const cs_int_t          id_alias,     /* --> numéro de l'alias à créer
                                        *     (< 0 pour alias réservé,
                                        *      > 0 pour alias utilisateur)    */
 const cs_int_t          id_maillage   /* --> numéro du maillage  associé     */
);


/*----------------------------------------------------------------------------
 * Vérifie l'existence d'un "writer" associé à un numéro donné.
 *----------------------------------------------------------------------------*/

cs_bool_t cs_post_existe_writer
(
 const cs_int_t   numwri        /* --> numéro du writer associé               */
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
 * Vérifie l'existence d'un maillage de post traitement associé à un
 * numéro donné.
 *----------------------------------------------------------------------------*/

cs_bool_t cs_post_existe_maillage
(
 const cs_int_t   nummai        /* --> numéro du maillage externe associé     */
);


/*----------------------------------------------------------------------------
 * Modification d'un maillage de post traitement existant.
 *
 * Il s'agit ici de modifier les listes de cellules ou faces du maillage,
 * par exemple pour faire évoluer une coupe en fonction des zones
 * "intéressantes (il n'est pas nécessaire de recourir à cette fonction
 * si le maillage se déforme simplement).
 *----------------------------------------------------------------------------*/

void cs_post_modifie_maillage
(
 const cs_int_t          id_maillage,  /* --> numéro du writer à créer
                                        *     (< 0 pour maillage réservé,
                                        *      > 0 pour maillage utilisateur) */
 const cs_int_t          nbr_cel,      /* --> nombre de cellules              */
 const cs_int_t          nbr_fac,      /* --> nombre de faces internes        */
 const cs_int_t          nbr_fbr,      /* --> nombre de faces de bord         */
       cs_int_t          liste_cel[],  /* <-> liste des cellules              */
       cs_int_t          liste_fac[],  /* <-> liste des faces internes        */
       cs_int_t          liste_fbr[]   /* <-> liste des faces de bord         */
);


/*----------------------------------------------------------------------------
 * Récupération du prochain numéro de maillage standard ou développeur
 * disponible (basé sur le plus petit numéro négatif présent -1).
 *----------------------------------------------------------------------------*/

cs_int_t cs_post_ret_num_maillage_libre
(
 void
);


/*----------------------------------------------------------------------------
 * Association d'un "writer" à un maillage pour le post traitement.
 *----------------------------------------------------------------------------*/

void cs_post_associe
(
 const cs_int_t   id_maillage,  /* --> numéro du maillage externe associé     */
 const cs_int_t   id_writer     /* --> numéro du writer                       */
);


/*----------------------------------------------------------------------------
 * Mise à jour de l'indicateur "actif" ou "inactif" des "writers" en
 * fonction du pas de temps et de leur fréquence de sortie par défaut.
 *----------------------------------------------------------------------------*/

void cs_post_activer_selon_defaut
(
 const cs_int_t   nt_cur_abs    /* --> numéro de pas de temps courant         */
);


/*----------------------------------------------------------------------------
 * Forcer de l'indicateur "actif" ou "inactif" d'un "writers" spécifique
 * ou de l'ensemble des "writers" pour le pas de temps en cours.
 *----------------------------------------------------------------------------*/

void cs_post_activer_writer
(
 const cs_int_t   id_writer,    /* --> numéro du writer,ou 0 pour forcer
                                 *     simultanément tous les writers         */
 const cs_int_t   activer       /* --> 0 pour désactiver, 1 pour activer      */
);


/*----------------------------------------------------------------------------
 * Ecriture des maillages de post traitement en fonction des "writers"
 * associés.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_maillages
(
 const cs_int_t   nt_cur_abs,         /* --> numéro de pas de temps courant   */
 const cs_real_t  t_cur_abs           /* --> valeur du temps physique associé */
);


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement défini sur les cellules ou faces
 * d'un maillage en fonction des "writers" associés.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_var
(
       cs_int_t          id_maillage,  /* --> numéro du maillage post associé */
 const char             *nom_var,      /* --> nom de la variable              */
       cs_int_t          dim_var,      /* --> 1 pour scalaire, 3 pour vecteur */
       cs_bool_t         entrelace,    /* --> si vecteur, vrai si valeurs
                                        *     entrelacées, faux sinon         */
       cs_bool_t         var_parent,   /* --> vrai si valeurs définies sur
                                        *     maillage "parent", faux si
                                        *     restreintes au maillage post    */
       cs_post_type_t    var_type,     /* --> type de données associé         */
       cs_int_t          nt_cur_abs,   /* --> numéro de pas de temps courant  */
       cs_real_t         t_cur_abs,    /* --> valeur du temps physique        */
 const void             *var_cel,      /* --> valeurs aux cellules            */
 const void             *var_fac,      /* --> valeurs aux faces internes      */
 const void             *var_fbr       /* --> valeurs aux faces de bord       */
);


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement défini sur les sommets
 * d'un maillage en fonction des "writers" associés.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_var_som
(
       cs_int_t          id_maillage,  /* --> numéro du maillage post associé */
 const char             *nom_var,      /* --> nom de la variable              */
       cs_int_t          dim_var,      /* --> 1 pour scalaire, 3 pour vecteur */
       cs_bool_t         entrelace,    /* --> si vecteur, vrai si valeurs
                                        *     entrelacées, faux sinon         */
       cs_bool_t         var_parent,   /* --> vrai si valeurs définies sur
                                        *     maillage "parent", faux si
                                        *     restreintes au maillage post    */
       cs_post_type_t    var_type,     /* --> type de données associé         */
       cs_int_t          nt_cur_abs,   /* --> numéro de pas de temps courant  */
       cs_real_t         t_cur_abs,    /* --> valeur du temps physique        */
 const void             *var_som       /* --> valeurs aux sommets             */
);


/*----------------------------------------------------------------------------
 * Prise en compte de la renumérotation des cellules
 * dans les liens de "parenté" des maillages post.
 *
 * Cette fonction ne doit être appellée qu'une fois, après la renumérotation
 * évuentuelle des cellules, pour adapter les maillages post existants.
 * Des nouveaux maillages post seront automatiquement basés sur la
 * "bonne" numérotation, par construction.
 *----------------------------------------------------------------------------*/

void cs_post_renum_cells
(
 const cs_int_t  *init_cell_num   /* --> numérotation initiale des cellules   */
);


/*----------------------------------------------------------------------------
 * Prise en compte de la renumérotation des faces et faces de bord
 * dans les liens de "parenté" des maillages post.
 *
 * Cette fonction ne doit être appellée qu'une fois, après la renumérotation
 * évuentuelle des faces, pour adapter les maillages post existants.
 * Des nouveaux maillages post seront automatiquement basés sur la
 * "bonne" numérotation, par construction.
 *----------------------------------------------------------------------------*/

void cs_post_renum_faces
(
 const cs_int_t  *init_i_face_num,  /* --> numérotation init. faces internes  */
 const cs_int_t  *init_b_face_num   /* --> numérotation init. faces de bord   */
);

/*----------------------------------------------------------------------------
 * Destruction des structures associées aux post traitements
 *----------------------------------------------------------------------------*/

void cs_post_detruit
(
 void
);


/*----------------------------------------------------------------------------
 * Initialisation du "writer" du post-traitement principal
 *----------------------------------------------------------------------------*/

void cs_post_init_pcp_writer
(
 void
);


/*----------------------------------------------------------------------------
 * Initialisation des maillages du post-traitement principal
 *----------------------------------------------------------------------------*/

void cs_post_init_pcp_maillages
(
 void
);


/*----------------------------------------------------------------------------
 * Ajout d'un traitement de variable temporelle à l'appel de PSTVAR.
 *
 * L'identificateur d'instance associé à la fonction permet d'ajouter
 * une même fonction plusieurs fois, avec un identificateur différent
 * permettant à la fonction de sélectionner un sous-traitement.
 *----------------------------------------------------------------------------*/

void cs_post_ajoute_var_temporelle
(
 cs_post_var_temporelle_t  *fonction,    /* Fonction associée                 */
 cs_int_t                   id_instance  /* Indentificateur d'instance
                                            associé à la fonction             */
 );


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_H__ */
