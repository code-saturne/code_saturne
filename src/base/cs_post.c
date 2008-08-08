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

/*============================================================================
 * Management of the post-processing
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 *  Définitions d'énumerations
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Type de support de maillage
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_POST_SUPPORT_CEL,           /* Valeurs associées aux cellules            */
  CS_POST_SUPPORT_FAC_INT,       /* Valeurs associées aux faces internes      */
  CS_POST_SUPPORT_FAC_BRD,       /* Valeurs associées aux faces de bord       */
  CS_POST_SUPPORT_SOM            /* Valeurs associées aux sommets             */

} cs_post_support_t;



/*============================================================================
 *  Structures locales
 *============================================================================*/

/* Structure définissant un "writer" : cet objet correspond au choix d'un
 * nom de cas, de répertoire, et de format, ainsi qu'un indicateur précisant
 * si les maillages associés doivent dépendre ou non du temps, et la
 * fréquence de sortie par défaut pour les variables associées. */

struct _cs_post_writer_t {

  cs_int_t                 id;           /* Identificateur (< 0 pour writer
                                            standard ou développeur, > 0 pour
                                            writer utilisateur */
  cs_int_t                 freq_sortie;  /* Fréquence de sortie par défaut
                                            associée */
  cs_bool_t                ecr_depl;     /* Ecriture d'un champ déplacement
                                            si vrai */

  int                      actif;        /* 0 si pas de sortie au pas de
                                            temps courant, 1 en cas de sortie */

  fvm_writer_t            *writer;       /* Gestionnaire d'écriture associé */

};


/* Structure définissant un maillage de post traitement ; cet objet
 * gère le lien entre un tel maillage et les "writers" associés. */

struct _cs_post_maillage_t {

  int                     id;            /* Identificateur (< 0 pour maillage
                                            standard ou développeur, > 0 pour
                                            maillage utilisateur */

  int                     ind_ent[3];    /* Présence de cellules (ind_ent[0],
                                            de faces internes (ind_ent[1]),
                                            ou de faces de bord (ind_ent[2])
                                            sur un processeur au moins */

  int                     alias;         /* Si > -1, indice dans le tableau
                                            des maillages de post traitement
                                            du premier maillage partageant
                                            le même maillage externe */

  int                     nbr_writers;   /* Nombre de gestionnaires de sortie */
  int                    *ind_writer;    /* Tableau des indices des
                                          * gestionnaires de sortie associés */
  int                     nt_ecr;        /* Numéro du pas de temps de la
                                            dernière écriture (-1 avant la
                                            première écriture) */

  cs_int_t                nbr_fac;       /* Nombre de faces internes associées */
  cs_int_t                nbr_fbr;       /* Nombre de faces de bord associées */

  const fvm_nodal_t      *maillage_ext;  /* Maillage externe associé */
  fvm_nodal_t            *_maillage_ext; /* Maillage externe associé, si
                                            propriétaire */

  fvm_writer_time_dep_t   ind_mod_min;   /* Indicateur de possibilité de
                                            modification au cours du temps */
  fvm_writer_time_dep_t   ind_mod_max;   /* Indicateur de possibilité de
                                            modification au cours du temps */
};


/*============================================================================
 *  Variables globales statiques
 *============================================================================*/

/* Tableau de sauvegarde des coordonnées initiales des sommets */

static cs_bool_t           cs_glob_post_deformable = CS_FALSE;
static cs_real_t          *cs_glob_post_coo_som_ini = NULL;

/* Indicateur de sortie du domaine en parallèle */

static cs_bool_t           cs_glob_post_domaine = CS_TRUE;


/* Tableau des maillages externes associés aux post-traitements */
/* (maillages, -1 et -2 réservés, donc numérotation commence  -2)*/
static int                 cs_glob_post_num_maillage_min = -2;
static int                 cs_glob_post_nbr_maillages = 0;
static int                 cs_glob_post_nbr_maillages_max = 0;
static cs_post_maillage_t *cs_glob_post_maillages = NULL;


/* Tableau des writers d'écriture associés aux post-traitements */

static int                 cs_glob_post_nbr_writers = 0;
static int                 cs_glob_post_nbr_writers_max = 0;
static cs_post_writer_t   *cs_glob_post_writers = NULL;


/* Tableau des fonctions et instances enregistrés */

static int                 cs_glob_post_nbr_var_tp = 0;
static int                 cs_glob_post_nbr_var_tp_max = 0;

static cs_post_var_temporelle_t  **cs_glob_post_f_var_tp = NULL;
static int                        *cs_glob_post_i_var_tp = NULL;


/*============================================================================
 *  Définitions de macros
 *============================================================================*/


/*============================================================================
 * Prototypes de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Conversion d'un type de données cs_post_type_t en type fvm_datatype_t.
 *----------------------------------------------------------------------------*/

static fvm_datatype_t _cs_post_cnv_datatype
(
 cs_post_type_t type_cs
);


/*----------------------------------------------------------------------------
 * Recherche de l'indice d'un writer associé à un numéro donné.
 *----------------------------------------------------------------------------*/

static int _cs_post_ind_writer
(
 const cs_int_t   id_writer     /* --> numéro du writer                       */
);


/*----------------------------------------------------------------------------
 * Recherche de l'indice d'un maillage de post traitement associé à un
 * numéro donné.
 *----------------------------------------------------------------------------*/

static int _cs_post_ind_maillage
(
 const cs_int_t   nummai        /* --> numéro du maillage externe associé     */
);


/*----------------------------------------------------------------------------
 * Ajout d'un maillage de post traitement et initialisation de base,
 * et renvoi d'un pointeur sur la structure associée
 *----------------------------------------------------------------------------*/

static cs_post_maillage_t * _cs_post_ajoute_maillage
(
 const cs_int_t  id_maillage           /* --> numéro du maillage  demandé     */
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
 * liste des faces de bord associée n'est donc pas nécessaire.
 *----------------------------------------------------------------------------*/

static void _cs_post_definit_maillage
(
 cs_post_maillage_t  *const maillage_post, /* <-> maillage post à compléter   */
 const char          *const nom_maillage,  /* --> nom du maillage externe     */
 const cs_int_t             nbr_cel,       /* --> nombre de cellules          */
 const cs_int_t             nbr_fac,       /* --> nombre de faces internes    */
 const cs_int_t             nbr_fbr,       /* --> nombre de faces de bord     */
       cs_int_t             liste_cel[],   /* <-> liste des cellules          */
       cs_int_t             liste_fac[],   /* <-> liste des faces internes    */
       cs_int_t             liste_fbr[]    /* <-> liste des faces de bord     */
);


/*----------------------------------------------------------------------------
 * Mise à jour en cas d'alias des critères de modification au cours du temps
 * des maillages en fonction des propriétés des writers qui lui
 * sont associés :
 *
 * La topologie d'un maillage ne pourra pas être modifiée si le critère
 * de modification minimal résultant est trop faible (i.e. si l'un des
 * writers associés ne permet pas la redéfinition de la topologie du maillage).
 *
 * Les coordonnées des sommets et la connectivité ne pourront être libérés
 * de la mémoire si la critère de modification maximal résultant est trop
 * élevé (i.e. si l'un des writers associés permet l'évolution du maillage
 * en temps, et nécessite donc sa réécriture).
 *----------------------------------------------------------------------------*/

static void _cs_post_ind_mod_alias
(
 const int  indmai              /* --> indice du maillage post en cours       */
);


/*----------------------------------------------------------------------------
 * Découpage des polygones ou polyèdres en éléments simples si nécessaire.
 *----------------------------------------------------------------------------*/

static void _cs_post_divise_poly
(
 cs_post_maillage_t      *const maillage_post,  /* --> maillage ext. associé  */
 const cs_post_writer_t  *const writer          /* --> writer associé         */
);


/*----------------------------------------------------------------------------
 * Assemblage des valeurs d'une variable définie sur une combinaison de
 * faces de bord et de faces internes (sans indirection) en un tableau
 * défini sur un ensemble unique de faces.
 *
 * La variable résultante n'est pas entrelacée.
 *----------------------------------------------------------------------------*/

static void _cs_post_assmb_var_faces
(
 const fvm_nodal_t      *const maillage_ext,     /* --> maillage externe   */
 const cs_int_t                nbr_fac,          /* --> nb faces internes  */
 const cs_int_t                nbr_fbr,          /* --> nb faces bord      */
 const int                     dim_var,          /* --> dim. variable      */
 const fvm_interlace_t         interlace,        /* --> indic. entrelacage */
 const cs_real_t               var_fac[],        /* --> valeurs faces int. */
 const cs_real_t               var_fbr[],        /* --> valeurs faces bord */
       cs_real_t               var_tmp[]         /* <-- valeurs assemblées */
);


/*----------------------------------------------------------------------------
 * Ecriture d'un maillage de post traitement en fonction des "writers".
 *----------------------------------------------------------------------------*/

static void _cs_post_ecrit_maillage
(
       cs_post_maillage_t  *const maillage_post,
 const cs_int_t                   nt_cur_abs,   /* --> numéro de pas de temps */
 const cs_real_t                  t_cur_abs     /* --> temps physique courant */
);


/*----------------------------------------------------------------------------
 * Transformation d'un tableau d'indicateurs (marqueurs) en liste ;
 * renvoie la taille effective de la liste.
 *----------------------------------------------------------------------------*/

static cs_int_t _cs_post_ind_vers_liste
(
 cs_int_t  nbr,                       /* <-> taille indicateur                */
 cs_int_t  liste[]                    /* <-> indicateur, puis liste           */
);


/*----------------------------------------------------------------------------
 * Boucle sur les maillages de post traitement pour écriture des variables
 *----------------------------------------------------------------------------*/

static void _cs_post_ecrit_deplacements
(
 const cs_int_t   nt_cur_abs,         /* --> numéro de pas de temps courant   */
 const cs_real_t  t_cur_abs           /* --> valeur du temps physique associé */
);


/*----------------------------------------------------------------------------
 * Écriture du domaine sur maillage de post traitement
 *----------------------------------------------------------------------------*/

static void _cs_post_ecrit_domaine
(
       fvm_writer_t   *writer,        /* --> writer FVM                       */
 const fvm_nodal_t    *maillage_ext,  /* --> maillage externe                 */
       cs_int_t        nt_cur_abs,    /* --> numéro pas de temps courant      */
       cs_real_t       t_cur_abs      /* --> valeur du temps physique associé */
);


/*============================================================================
 * Prototypes de fonctions Fortran appellées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction développeur pour la sortie de variables sur un maillage post
 *----------------------------------------------------------------------------*/

void CS_PROCF (dvvpst, DVVPST)
(
 const cs_int_t  *const idbia0,    /* --> numéro 1ère case libre dans IA      */
 const cs_int_t  *const idbra0,    /* --> numéro 1ère case libre dans RA      */
 const cs_int_t  *const nummai,    /* --> numéro du maillage post             */
 const cs_int_t  *const ndim,      /* --> dimension de l'espace               */
 const cs_int_t  *const ncelet,    /* --> nombre de cellules étendu           */
 const cs_int_t  *const ncel,      /* --> nombre de cellules                  */
 const cs_int_t  *const nfac,      /* --> nombre de faces internes            */
 const cs_int_t  *const nfabor,    /* --> nombre de faces de bord             */
 const cs_int_t  *const nfml,      /* --> nombre de familles                  */
 const cs_int_t  *const nprfml,    /* --> nombre de proprietes des familles   */
 const cs_int_t  *const nnod,      /* --> nombre de sommets                   */
 const cs_int_t  *const lndfac,    /* --> longueur de nodfac                  */
 const cs_int_t  *const lndfbr,    /* --> longueur de nodfbr                  */
 const cs_int_t  *const ncelbr,    /* --> nombre de cellules de bord          */
 const cs_int_t  *const nvar,      /* --> nombre de variables                 */
 const cs_int_t  *const nscal,     /* --> nombre de scalaires                 */
 const cs_int_t  *const nphas,     /* --> nombre de phases                    */
 const cs_int_t  *const nvlsta,    /* --> nombre de variables stat. (lagr)    */
 const cs_int_t  *const nvisbr,    /* --> nombre de variables stat. (lagr)    */
 const cs_int_t  *const ncelps,    /* --> nombre de cellules post             */
 const cs_int_t  *const nfacps,    /* --> nombre de faces internes post       */
 const cs_int_t  *const nfbrps,    /* --> nombre de faces de bord post        */
 const cs_int_t  *const nideve,    /* --> longueur du tableau idevel[]        */
 const cs_int_t  *const nrdeve,    /* --> longueur du tableau rdevel[]        */
 const cs_int_t  *const nituse,    /* --> longueur du tableau ituser[]        */
 const cs_int_t  *const nrtuse,    /* --> longueur du tableau rtuser[]        */
 const cs_int_t         itypps[],  /* --> indicateur (0 ou 1) de présence     */
                                   /*     de cellules, faces internes et bord */
 const cs_int_t         ifacel[],  /* --> connect. faces internes / cellules  */
 const cs_int_t         ifabor[],  /* --> connect. faces de bord / cellules   */
 const cs_int_t         ifmfbr[],  /* --> liste des familles des faces bord   */
 const cs_int_t         ifmcel[],  /* --> liste des familles des cellules     */
 const cs_int_t         iprfml[],  /* --> liste des propriétés des familles   */
 const cs_int_t         ipnfac[],  /* --> rang ds nodfac 1er sommet faces int */
 const cs_int_t         nodfac[],  /* --> numéro des sommets des faces int    */
 const cs_int_t         ipnfbr[],  /* --> rang ds nodfbr 1er sommt faces brd  */
 const cs_int_t         nodfbr[],  /* --> numéro des sommets des faces bord   */
 const cs_int_t         lstcel[],  /* --> liste des cellules post             */
 const cs_int_t         lstfac[],  /* --> liste des faces internes post       */
 const cs_int_t         lstfbr[],  /* --> liste des faces de bord post        */
 const cs_int_t         idevel[],  /* --> tab. complémentaire développeur     */
 const cs_int_t         ituser[],  /* --> tab. complémentaire utilisateur     */
 const cs_int_t         ia[],      /* --> macro-tableau entier                */
 const cs_real_t        xyzcen[],  /* --> c.d.g. des cellules                 */
 const cs_real_t        surfac[],  /* --> surfaces des faces internes         */
 const cs_real_t        surfbo[],  /* --> surfaces des faces de bord          */
 const cs_real_t        cdgfac[],  /* --> c.d.g. des faces internes           */
 const cs_real_t        cdgfbo[],  /* --> c.d.g. des faces de bord            */
 const cs_real_t        xyznod[],  /* --> coordonnées des sommets             */
 const cs_real_t        volume[],  /* --> volumes des cellules                */
 const cs_real_t        dt[],      /* --> pas de temps                        */
 const cs_real_t        rtpa[],    /* --> variables aux cellules (préc.)      */
 const cs_real_t        rtp[],     /* --> variables aux cellules              */
 const cs_real_t        propce[],  /* --> propriétés physiques cellules       */
 const cs_real_t        propfa[],  /* --> propriétés physiques aux faces      */
 const cs_real_t        propfb[],  /* --> propriétés physiques faces bord     */
 const cs_real_t        coefa[],   /* --> cond. limites aux faces de bord     */
 const cs_real_t        coefb[],   /* --> cond. limites aux faces de bord     */
 const cs_real_t        statce[],  /* --> moy. statistiques (Lagrangien)      */
 const cs_real_t        stativ[],  /* --> var. statistiques (Lagrangien)      */
 const cs_real_t        statfb[],  /* --> moy. statistiques (Lagrangien)      */
 const cs_real_t        valcel[],  /* --- vals. aux cellules post             */
 const cs_real_t        valfac[],  /* --- vals. aux faces internes post       */
 const cs_real_t        valfbr[],  /* --- vals. aux faces de bord post        */
 const cs_real_t        rdevel[],  /* --> tab. complémentaire développeur     */
 const cs_real_t        rtuser[],  /* --> tab. complémentaire utilisateur     */
 const cs_real_t        ra[]       /* --> macro-tableau réel                  */
);


/*----------------------------------------------------------------------------
 * Fonction utilisateur pour la modification d'un maillage post
 *----------------------------------------------------------------------------*/

void CS_PROCF (usmpst, USMPST)
(
 const cs_int_t  *const idbia0,    /* --> numéro 1ère case libre dans IA      */
 const cs_int_t  *const idbra0,    /* --> numéro 1ère case libre dans RA      */
 const cs_int_t  *const nummai,    /* --> numéro du maillage post             */
 const cs_int_t  *const ndim,      /* --> dimension de l'espace               */
 const cs_int_t  *const ncelet,    /* --> nombre de cellules étendu           */
 const cs_int_t  *const ncel,      /* --> nombre de cellules                  */
 const cs_int_t  *const nfac,      /* --> nombre de faces internes            */
 const cs_int_t  *const nfabor,    /* --> nombre de faces de bord             */
 const cs_int_t  *const nfml,      /* --> nombre de familles                  */
 const cs_int_t  *const nprfml,    /* --> nombre de proprietes des familles   */
 const cs_int_t  *const nnod,      /* --> nombre de sommets                   */
 const cs_int_t  *const lndfac,    /* --> longueur de nodfac                  */
 const cs_int_t  *const lndfbr,    /* --> longueur de nodfbr                  */
 const cs_int_t  *const ncelbr,    /* --> nombre de cellules de bord          */
 const cs_int_t  *const nvar,      /* --> nombre de variables                 */
 const cs_int_t  *const nscal,     /* --> nombre de scalaires                 */
 const cs_int_t  *const nphas,     /* --> nombre de phases                    */
 const cs_int_t  *const nvlsta,    /* --> nombre de variables stat. (lagr)    */
 const cs_int_t  *const ncelps,    /* --> nombre de cellules post             */
 const cs_int_t  *const nfacps,    /* --> nombre de faces internes post       */
 const cs_int_t  *const nfbrps,    /* --> nombre de faces de bord post        */
 const cs_int_t  *const nideve,    /* --> longueur du tableau idevel[]        */
 const cs_int_t  *const nrdeve,    /* --> longueur du tableau rdevel[]        */
 const cs_int_t  *const nituse,    /* --> longueur du tableau ituser[]        */
 const cs_int_t  *const nrtuse,    /* --> longueur du tableau rtuser[]        */
       cs_int_t  *const imodif,    /* <-- 0 si maillage non modifié, 1 sinon  */
 const cs_int_t         itypps[],  /* --> indicateur (0 ou 1) de présence     */
                                   /*     de cellules, faces internes et bord */
 const cs_int_t         ifacel[],  /* --> connect. faces internes / cellules  */
 const cs_int_t         ifabor[],  /* --> connect. faces de bord / cellules   */
 const cs_int_t         ifmfbr[],  /* --> liste des familles des faces bord   */
 const cs_int_t         ifmcel[],  /* --> liste des familles des cellules     */
 const cs_int_t         iprfml[],  /* --> liste des propriétés des familles   */
 const cs_int_t         ipnfac[],  /* --> rang ds nodfac 1er sommet faces int */
 const cs_int_t         nodfac[],  /* --> numéro des sommets des faces int    */
 const cs_int_t         ipnfbr[],  /* --> rang ds nodfbr 1er sommt faces brd  */
 const cs_int_t         nodfbr[],  /* --> numéro des sommets des faces bord   */
 const cs_int_t         lstcel[],  /* --> liste des cellules post             */
 const cs_int_t         lstfac[],  /* --> liste des faces internes post       */
 const cs_int_t         lstfbr[],  /* --> liste des faces de bord post        */
 const cs_int_t         idevel[],  /* --> tab. complémentaire développeur     */
 const cs_int_t         ituser[],  /* --> tab. complémentaire utilisateur     */
 const cs_int_t         ia[],      /* --> macro-tableau entier                */
 const cs_real_t        xyzcen[],  /* --> c.d.g. des cellules                 */
 const cs_real_t        surfac[],  /* --> surfaces des faces internes         */
 const cs_real_t        surfbo[],  /* --> surfaces des faces de bord          */
 const cs_real_t        cdgfac[],  /* --> c.d.g. des faces internes           */
 const cs_real_t        cdgfbo[],  /* --> c.d.g. des faces de bord            */
 const cs_real_t        xyznod[],  /* --> coordonnées des sommets             */
 const cs_real_t        volume[],  /* --> volumes des cellules                */
 const cs_real_t        dt[],      /* --> pas de temps                        */
 const cs_real_t        rtpa[],    /* --> variables aux cellules (préc.)      */
 const cs_real_t        rtp[],     /* --> variables aux cellules              */
 const cs_real_t        propce[],  /* --> propriétés physiques cellules       */
 const cs_real_t        propfa[],  /* --> propriétés physiques aux faces      */
 const cs_real_t        propfb[],  /* --> propriétés physiques faces bord     */
 const cs_real_t        coefa[],   /* --> cond. limites aux faces de bord     */
 const cs_real_t        coefb[],   /* --> cond. limites aux faces de bord     */
 const cs_real_t        statce[],  /* --> moy. statistiques (Lagrangien)      */
 const cs_real_t        valcel[],  /* --- vals. aux cellules post             */
 const cs_real_t        valfac[],  /* --- vals. aux faces internes post       */
 const cs_real_t        valfbr[],  /* --- vals. aux faces de bord post        */
 const cs_real_t        rdevel[],  /* --> tab. complémentaire développeur     */
 const cs_real_t        rtuser[],  /* --> tab. complémentaire utilisateur     */
 const cs_real_t        ra[]       /* --> macro-tableau réel                  */
);


/*----------------------------------------------------------------------------
 * Fonction utilisateur pour la sortie de variables sur un maillage post
 *----------------------------------------------------------------------------*/

void CS_PROCF (usvpst, USVPST)
(
 const cs_int_t  *const idbia0,    /* --> numéro 1ère case libre dans IA      */
 const cs_int_t  *const idbra0,    /* --> numéro 1ère case libre dans RA      */
 const cs_int_t  *const nummai,    /* --> numéro du maillage post             */
 const cs_int_t  *const ndim,      /* --> dimension de l'espace               */
 const cs_int_t  *const ncelet,    /* --> nombre de cellules étendu           */
 const cs_int_t  *const ncel,      /* --> nombre de cellules                  */
 const cs_int_t  *const nfac,      /* --> nombre de faces internes            */
 const cs_int_t  *const nfabor,    /* --> nombre de faces de bord             */
 const cs_int_t  *const nfml,      /* --> nombre de familles                  */
 const cs_int_t  *const nprfml,    /* --> nombre de proprietes des familles   */
 const cs_int_t  *const nnod,      /* --> nombre de sommets                   */
 const cs_int_t  *const lndfac,    /* --> longueur de nodfac                  */
 const cs_int_t  *const lndfbr,    /* --> longueur de nodfbr                  */
 const cs_int_t  *const ncelbr,    /* --> nombre de cellules de bord          */
 const cs_int_t  *const nvar,      /* --> nombre de variables                 */
 const cs_int_t  *const nscal,     /* --> nombre de scalaires                 */
 const cs_int_t  *const nphas,     /* --> nombre de phases                    */
 const cs_int_t  *const nvlsta,    /* --> nombre de variables stat. (lagr)    */
 const cs_int_t  *const ncelps,    /* --> nombre de cellules post             */
 const cs_int_t  *const nfacps,    /* --> nombre de faces internes post       */
 const cs_int_t  *const nfbrps,    /* --> nombre de faces de bord post        */
 const cs_int_t  *const nideve,    /* --> longueur du tableau idevel[]        */
 const cs_int_t  *const nrdeve,    /* --> longueur du tableau rdevel[]        */
 const cs_int_t  *const nituse,    /* --> longueur du tableau ituser[]        */
 const cs_int_t  *const nrtuse,    /* --> longueur du tableau rtuser[]        */
 const cs_int_t         itypps[],  /* --> indicateur (0 ou 1) de présence     */
                                   /*     de cellules, faces internes et bord */
 const cs_int_t         ifacel[],  /* --> connect. faces internes / cellules  */
 const cs_int_t         ifabor[],  /* --> connect. faces de bord / cellules   */
 const cs_int_t         ifmfbr[],  /* --> liste des familles des faces bord   */
 const cs_int_t         ifmcel[],  /* --> liste des familles des cellules     */
 const cs_int_t         iprfml[],  /* --> liste des propriétés des familles   */
 const cs_int_t         ipnfac[],  /* --> rang ds nodfac 1er sommet faces int */
 const cs_int_t         nodfac[],  /* --> numéro des sommets des faces int    */
 const cs_int_t         ipnfbr[],  /* --> rang ds nodfbr 1er sommt faces brd  */
 const cs_int_t         nodfbr[],  /* --> numéro des sommets des faces bord   */
 const cs_int_t         lstcel[],  /* --> liste des cellules post             */
 const cs_int_t         lstfac[],  /* --> liste des faces internes post       */
 const cs_int_t         lstfbr[],  /* --> liste des faces de bord post        */
 const cs_int_t         idevel[],  /* --> tab. complémentaire développeur     */
 const cs_int_t         ituser[],  /* --> tab. complémentaire utilisateur     */
 const cs_int_t         ia[],      /* --> macro-tableau entier                */
 const cs_real_t        xyzcen[],  /* --> c.d.g. des cellules                 */
 const cs_real_t        surfac[],  /* --> surfaces des faces internes         */
 const cs_real_t        surfbo[],  /* --> surfaces des faces de bord          */
 const cs_real_t        cdgfac[],  /* --> c.d.g. des faces internes           */
 const cs_real_t        cdgfbo[],  /* --> c.d.g. des faces de bord            */
 const cs_real_t        xyznod[],  /* --> coordonnées des sommets             */
 const cs_real_t        volume[],  /* --> volumes des cellules                */
 const cs_real_t        dt[],      /* --> pas de temps                        */
 const cs_real_t        rtpa[],    /* --> variables aux cellules (préc.)      */
 const cs_real_t        rtp[],     /* --> variables aux cellules              */
 const cs_real_t        propce[],  /* --> propriétés physiques cellules       */
 const cs_real_t        propfa[],  /* --> propriétés physiques aux faces      */
 const cs_real_t        propfb[],  /* --> propriétés physiques faces bord     */
 const cs_real_t        coefa[],   /* --> cond. limites aux faces de bord     */
 const cs_real_t        coefb[],   /* --> cond. limites aux faces de bord     */
 const cs_real_t        statce[],  /* --> moy. statistiques (Lagrangien)      */
 const cs_real_t        valcel[],  /* --- vals. aux cellules post             */
 const cs_real_t        valfac[],  /* --- vals. aux faces internes post       */
 const cs_real_t        valfbr[],  /* --- vals. aux faces de bord post        */
 const cs_real_t        rdevel[],  /* --> tab. complémentaire développeur     */
 const cs_real_t        rtuser[],  /* --> tab. complémentaire utilisateur     */
 const cs_real_t        ra[]       /* --> macro-tableau réel                  */
);


/*----------------------------------------------------------------------------
 * Fonction récupèrant les paramètres contenus dans les commons FORTRAN qui
 * sont utiles à l'initialisation du post-traitement.
 *----------------------------------------------------------------------------*/

void CS_PROCF (inipst, INIPST)
(
 const cs_int_t  *const ichrvl,    /* --> indic. de post du volume fluide     */
 const cs_int_t  *const ichrbo,    /* --> indic. de post des faces de bord    */
 const cs_int_t  *const ichrsy,    /* --> indic. de post des faces de bord    */
                                   /*     couplées avec Syrthes               */
 const cs_int_t  *const ipstmd,    /* --> indic. de maillage déformable       */
                                   /*     0 : pas de déformation              */
                                   /*     1 : déformation des maillages post  */
                                   /*     2 : écriture d'un champ déplacement */
 const cs_int_t  *const ntchr,     /* --> fréquence des sorties post          */
 const char      *const fmtchr,    /* --> nom du format de post-traitement    */
 const char      *const optchr     /* --> options du format de post           */
);


/*============================================================================
 * Fonctions publiques pour l'API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Création d'un "writer" à partir des données du Fortran.
 *
 * Interface Fortran : utiliser PSTCWR (voir cs_post_util.F)
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
)
{
  /* variables locales */

  char  *nom_cas;
  char  *nom_rep;
  char  *nom_format;
  char  *opt_format;

  /* Conversion des chaînes de caractères Fortran en chaînes C */

  nom_cas    = cs_base_chaine_f_vers_c_cree(nomcas, *lnmcas);
  nom_rep    = cs_base_chaine_f_vers_c_cree(nomrep, *lnmrep);
  nom_format = cs_base_chaine_f_vers_c_cree(nomfmt, *lnmfmt);
  opt_format = cs_base_chaine_f_vers_c_cree(optfmt, *lopfmt);

  /* Traitement principal */

  cs_post_ajoute_writer(*numwri,
                        nom_cas,
                        nom_rep,
                        nom_format,
                        opt_format,
                        *indmod,
                        *ntchr);

  /* Libération des chaînes C temporaires */

  nom_cas = cs_base_chaine_f_vers_c_detruit(nom_cas);
  nom_rep = cs_base_chaine_f_vers_c_detruit(nom_rep);
  nom_format = cs_base_chaine_f_vers_c_detruit(nom_format);
  opt_format = cs_base_chaine_f_vers_c_detruit(opt_format);

}


/*----------------------------------------------------------------------------
 * Création d'un maillage de post traitement à partir des données du
 * Fortran.
 *
 * Interface Fortran : utiliser PSTCMA (voir cs_post_util.F)
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
)
{
  /* variables locales */

  char  *nom_maillage = NULL;


  /* Conversion des chaînes de caractères Fortran en chaînes C */

  nom_maillage = cs_base_chaine_f_vers_c_cree(nommai, *lnmmai);

  /* Traitement principal */

  cs_post_ajoute_maillage(*nummai,
                          nom_maillage,
                          *nbrcel,
                          *nbrfac,
                          *nbrfbr,
                          lstcel,
                          lstfac,
                          lstfbr);

  /* Libération des chaînes C temporaires */

  nom_maillage = cs_base_chaine_f_vers_c_detruit(nom_maillage);

}


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
)
{
  cs_post_alias_maillage(*nummai, *numref);
}


/*----------------------------------------------------------------------------
 * Association d'un "writer" à un maillage pour le post traitement.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTASS (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : --> : Numéro du maillage externe associé
 * INTEGER          NUMWRI      : --> : Numéro du writer
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstass, PSTASS)
(
 const cs_int_t   *nummai,      /* --> numéro du maillage externe associé     */
 const cs_int_t   *numwri       /* --> numéro du writer                       */
)
{
  cs_post_associe(*nummai, *numwri);
}


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
)
{
  cs_post_activer_selon_defaut(*ntcabs);
}


/*----------------------------------------------------------------------------
 * Forcer de l'indicateur "actif" ou "inactif" d'un "writers" spécifique
 * ou de l'ensemble des "writers" pour le pas de temps en cours.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTACT (NTCABS, TTCABS)
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
)
{
  cs_post_activer_writer(*numwri,
                         *indact);
}


/*----------------------------------------------------------------------------
 * Ecriture des maillages de post traitement en fonction des gestionnaires
 * de sortie associés.
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
)
{
  cs_post_ecrit_maillages(*ntcabs, *ttcabs);
}


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
)
{
  /* variables locales */

  int i, j, k;
  int dim_ent;
  cs_int_t   itypps[3];
  cs_int_t   ind_cel, ind_fac, dec_num_fbr;
  cs_int_t   nbr_ent, nbr_ent_max;
  cs_int_t   nummai, imodif;

  cs_bool_t  actif;

  cs_post_maillage_t  *maillage_post;
  cs_post_writer_t  *writer;

  cs_int_t    nbr_cel, nbr_fac, nbr_fbr;
  cs_int_t    *liste_cel, *liste_fac, *liste_fbr;

  cs_int_t    *num_ent_parent = NULL;
  cs_real_t   *var_trav = NULL;
  cs_real_t   *var_cel = NULL;
  cs_real_t   *var_fac = NULL;
  cs_real_t   *var_fbr = NULL;


  /* Boucle sur les writers pour vérifier si l'on a quelque chose à faire */
  /*----------------------------------------------------------------------*/

  for (j = 0 ; j < cs_glob_post_nbr_writers ; j++) {
    writer = cs_glob_post_writers + j;
    if (writer->actif == 1)
      break;
  }
  if (j == cs_glob_post_nbr_writers)
    return;


  /* Modification éventuelle des définitions des maillages post */
  /*------------------------------------------------------------*/

  nbr_ent_max = 0;

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    actif = CS_FALSE;

    for (j = 0 ; j < maillage_post->nbr_writers ; j++) {
      writer = cs_glob_post_writers + maillage_post->ind_writer[j];
      if (writer->actif == 1)
        actif = CS_TRUE;
    }

    /* Maillage utilisateur modifiable, non alias, actif à ce pas de temps */

    if (   actif == CS_TRUE
        && maillage_post->alias < 0
        && maillage_post->id > 0
        && maillage_post->ind_mod_min == FVM_WRITER_TRANSIENT_CONNECT) {

      const fvm_nodal_t * maillage_ext = maillage_post->maillage_ext;

      dim_ent = fvm_nodal_get_max_entity_dim(maillage_ext);
      nbr_ent = fvm_nodal_get_n_entities(maillage_ext, dim_ent);

      if (nbr_ent > nbr_ent_max) {
        nbr_ent_max = nbr_ent;
        BFT_REALLOC(num_ent_parent, nbr_ent_max, cs_int_t);
      }

      nummai = maillage_post->id;

      /* Récupération des listes d'entités correspondantes */

      fvm_nodal_get_parent_num(maillage_ext, dim_ent, num_ent_parent);

      for (k = 0 ; k < 3 ; k++)
        itypps[k] = maillage_post->ind_ent[k];


      /* On surdimensionne ici les listes, l'utilisateur pouvant
         en remplir une partie arbitraire */

      BFT_MALLOC(liste_cel, cs_glob_mesh->n_cells, cs_int_t);
      BFT_MALLOC(liste_fac, cs_glob_mesh->n_i_faces, cs_int_t);
      BFT_MALLOC(liste_fbr, cs_glob_mesh->n_b_faces, cs_int_t);

      nbr_cel = 0;
      nbr_fac = 0;
      nbr_fbr = 0;

      /* Mise à zéro des listes */

      if (dim_ent == 3)
        for (ind_cel = 0 ; ind_cel < cs_glob_mesh->n_cells ; ind_cel++)
          liste_cel[ind_cel] = 0;
      else if (dim_ent == 2) {
        for (ind_fac = 0 ; ind_fac < cs_glob_mesh->n_b_faces ; ind_fac++)
          liste_fbr[ind_fac] = 0;
        for (ind_fac = 0 ; ind_fac < cs_glob_mesh->n_i_faces ; ind_fac++)
          liste_fac[ind_fac] = 0;
      }

      /* Si les éléments du maillage FVM sont découpés, un même numéro
         parent peut apparaître plusieurs fois ; On utilise donc une
         logique par indicateur. */

      if (dim_ent == 3) {
        for (ind_cel = 0 ; ind_cel < nbr_ent ; ind_cel++)
          liste_cel[num_ent_parent[ind_cel] - 1] = 1;
      }

      /* Pour les faces, les numéros de faces internes "parentes"
         connus par FVM décalés du nombre total de faces de bord
         (c.f. construction dans cs_maillage_extrait...()) */

      else if (dim_ent == 2) {
        dec_num_fbr = cs_glob_mesh->n_b_faces;
        for (ind_fac = 0 ; ind_fac < nbr_ent ; ind_fac++) {
          if (num_ent_parent[ind_fac] > dec_num_fbr)
            liste_fac[num_ent_parent[ind_fac] - dec_num_fbr - 1] = 1;
          else
            liste_fbr[num_ent_parent[ind_fac] - 1] = 1;
        }
      }

      /* Transformation des indicateurs en listes */

      if (dim_ent == 3) {
        nbr_cel = _cs_post_ind_vers_liste(cs_glob_mesh->n_cells,
                                          liste_cel);
      }
      else if (dim_ent == 2) {
        nbr_fac = _cs_post_ind_vers_liste(cs_glob_mesh->n_i_faces,
                                          liste_fac);
        nbr_fbr = _cs_post_ind_vers_liste(cs_glob_mesh->n_b_faces,
                                          liste_fbr);
      }

      /* Modification de la définition du maillage par l'utilisateur */

      imodif = 0;

      CS_PROCF(usmpst, USMPST) (idbia0, idbra0, &nummai,
                                ndim, ncelet, ncel, nfac, nfabor, nfml, nprfml,
                                nnod, lndfac, lndfbr, ncelbr,
                                nvar, nscal, nphas, nvlsta,
                                &nbr_cel, &nbr_fac, &nbr_fbr,
                                nideve, nrdeve, nituse, nrtuse, &imodif,
                                itypps, ifacel, ifabor, ifmfbr, ifmcel, iprfml,
                                ipnfac, nodfac, ipnfbr, nodfbr,
                                liste_cel, liste_fac, liste_fbr,
                                idevel, ituser, ia,
                                xyzcen, surfac, surfbo, cdgfac, cdgfbo, xyznod,
                                volume, dt, rtpa, rtp, propce, propfa, propfb,
                                coefa, coefb, statce,
                                var_cel, var_fac, var_fbr,
                                rdevel, rtuser, ra);

      if (imodif > 0)
        cs_post_modifie_maillage(maillage_post->id,
                                 nbr_cel,
                                 nbr_fac,
                                 nbr_fbr,
                                 liste_cel,
                                 liste_fac,
                                 liste_fbr);

      BFT_FREE(liste_cel);
      BFT_FREE(liste_fac);
      BFT_FREE(liste_fbr);

    }

  }

  /* On s'assure maintenant de la synchronisation des alias */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    if (maillage_post->alias > -1) {

      const cs_post_maillage_t  *maillage_ref;

      maillage_ref = cs_glob_post_maillages + maillage_post->alias;

      for (j = 0 ; j < 3 ; j++)
        maillage_post->ind_ent[j] = maillage_ref->ind_ent[j];

      maillage_post->nbr_fac = maillage_ref->nbr_fac;
      maillage_post->nbr_fbr = maillage_ref->nbr_fbr;

    }

  }

  /* Sortie des maillages ou champs de déplacement des sommets si nécessaire */
  /*-------------------------------------------------------------------------*/

  cs_post_ecrit_maillages(*ntcabs, *ttcabs);

  if (cs_glob_post_deformable == CS_TRUE)
    _cs_post_ecrit_deplacements(*ntcabs, *ttcabs);


  /* Sorties des variables réalisés par des traitements enregistrés */
  /*----------------------------------------------------------------*/

  for (i = 0; i < cs_glob_post_nbr_var_tp; i++) {
    cs_glob_post_f_var_tp[i](cs_glob_post_i_var_tp[i],
                             *ntcabs,
                             *ttcabs);
  }

  /* Sortie des variables associées aux maillages de post traitement */
  /*-----------------------------------------------------------------*/

  /* nbr_ent_max déja initialisé avant et au cours de la
     modification éventuelle des définitions de maillages post,
     et num_ent_parent alloué si nbr_ent_max > 0 */

  BFT_MALLOC(var_trav, nbr_ent_max * 3, cs_real_t);

  /* Boucle principale sur les maillages de post traitment */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    actif = CS_FALSE;

    for (j = 0 ; j < maillage_post->nbr_writers ; j++) {
      writer = cs_glob_post_writers + maillage_post->ind_writer[j];
      if (writer->actif == 1)
        actif = CS_TRUE;
    }

    /* Si le maillage est actif à ce pas de temps */
    /*--------------------------------------------*/

    if (actif == CS_TRUE) {

      const fvm_nodal_t * maillage_ext = maillage_post->maillage_ext;

      dim_ent = fvm_nodal_get_max_entity_dim(maillage_ext);
      nbr_ent = fvm_nodal_get_n_entities(maillage_ext, dim_ent);

      if (nbr_ent > nbr_ent_max) {
        nbr_ent_max = nbr_ent;
        BFT_REALLOC(var_trav, nbr_ent_max * 3, cs_real_t);
        BFT_REALLOC(num_ent_parent, nbr_ent_max, cs_int_t);
      }

      nummai = maillage_post->id;


      /* Récupération des listes d'entités correspondantes */

      fvm_nodal_get_parent_num(maillage_ext, dim_ent, num_ent_parent);

      for (k = 0 ; k < 3 ; k++)
        itypps[k] = maillage_post->ind_ent[k];


      /* On peut sortir des variables pour ce pas de temps */
      /*---------------------------------------------------*/

      nbr_cel = 0;
      nbr_fac = 0;
      nbr_fbr = 0;
      liste_cel = NULL;
      liste_fac = NULL;
      liste_fbr = NULL;

      /* Ici, les listes sont dimensionnées au plus juste, et on pointe
         sur le tableau rempli par fvm_nodal_get_parent_num() si possible. */

      if (dim_ent == 3) {
        nbr_cel = nbr_ent;
        liste_cel = num_ent_parent;
      }

      /* Les numéros de faces internes "parentes" connus par FVM
         sont décalés du nombre total de faces de bord */

      else if (dim_ent == 2 && nbr_ent > 0) {

        dec_num_fbr = cs_glob_mesh->n_b_faces;

        for (ind_fac = 0 ; ind_fac < nbr_ent ; ind_fac++) {
          if (num_ent_parent[ind_fac] > dec_num_fbr)
            nbr_fac++;
          else
            nbr_fbr++;
        }

        /* faces de bord seulement : numéros faces parentes FVM adaptés */
        if (nbr_fac == 0) {
          liste_fbr = num_ent_parent;
        }

        /* faces internes seulement : numéros faces parentes FVM décalés */
        else if (nbr_fbr == 0) {
          for (ind_fac = 0 ; ind_fac < nbr_ent ; ind_fac++)
            num_ent_parent[ind_fac] -= dec_num_fbr;
          liste_fac = num_ent_parent;
        }

        /* faces internes et de bord : numéros à séparer */

        else {

          BFT_MALLOC(liste_fac, nbr_fac, cs_int_t);
          BFT_MALLOC(liste_fbr, nbr_fbr, cs_int_t);

          nbr_fac = 0, nbr_fbr = 0;

          for (ind_fac = 0 ; ind_fac < nbr_ent ; ind_fac++) {
            if (num_ent_parent[ind_fac] > dec_num_fbr)
              liste_fac[nbr_fac++] = num_ent_parent[ind_fac] - dec_num_fbr;
            else
              liste_fbr[nbr_fbr++] = num_ent_parent[ind_fac];
          }

        }

        /* Dans tous les cas, mise à jour du nombre de faces internes
           et faces de bord (utile en cas de découpage du maillage FVM)
           pour les fonctions appellées par celle-ci */

        maillage_post->nbr_fac = nbr_fac;
        maillage_post->nbr_fbr = nbr_fbr;

      }

      /* Pointeurs sur tableaux d'assemblage des variables,
         mis à NULL si inutiles (afin de provoquer si possible
         une erreur immédiate en cas de mauvaise utilisation) */

      var_cel = var_trav;
      var_fac = var_cel + (nbr_cel * 3);
      var_fbr = var_fac + (nbr_fac * 3);

      if (nbr_cel == 0)
        var_cel = NULL;
      if (nbr_fac == 0)
        var_fac = NULL;
      if (nbr_fbr == 0)
        var_fbr = NULL;

      /* Post traitement automatique des variables */

      if (nummai < 0)
        CS_PROCF(dvvpst, DVVPST) (idbia0, idbra0, &nummai,
                                  ndim, ncelet, ncel, nfac, nfabor, nfml, nprfml,
                                  nnod, lndfac, lndfbr, ncelbr,
                                  nvar, nscal, nphas, nvlsta, nvisbr,
                                  &nbr_cel, &nbr_fac, &nbr_fbr,
                                  nideve, nrdeve, nituse, nrtuse,
                                  itypps, ifacel, ifabor, ifmfbr, ifmcel, iprfml,
                                  ipnfac, nodfac, ipnfbr, nodfbr,
                                  liste_cel, liste_fac, liste_fbr,
                                  idevel, ituser, ia,
                                  xyzcen, surfac, surfbo, cdgfac, cdgfbo, xyznod,
                                  volume, dt, rtpa, rtp, propce, propfa, propfb,
                                  coefa, coefb, statce, stativ , statfb ,
                                  var_cel, var_fac, var_fbr,
                                  rdevel, rtuser, ra);

      /* Appel Fortran utilisateur pour post traitement des variables */

      CS_PROCF(usvpst, USVPST) (idbia0, idbra0, &nummai,
                                ndim, ncelet, ncel, nfac, nfabor, nfml, nprfml,
                                nnod, lndfac, lndfbr, ncelbr,
                                nvar, nscal, nphas, nvlsta,
                                &nbr_cel, &nbr_fac, &nbr_fbr,
                                nideve, nrdeve, nituse, nrtuse,
                                itypps, ifacel, ifabor, ifmfbr, ifmcel, iprfml,
                                ipnfac, nodfac, ipnfbr, nodfbr,
                                liste_cel, liste_fac, liste_fbr,
                                idevel, ituser, ia,
                                xyzcen, surfac, surfbo, cdgfac, cdgfbo, xyznod,
                                volume, dt, rtpa, rtp, propce, propfa, propfb,
                                coefa, coefb, statce,
                                var_cel, var_fac, var_fbr,
                                rdevel, rtuser, ra);


      /* En cas de mélange de faces internes et de bord, tableaux
         supplémentaires alloués, à libérer */

      if (liste_fac != NULL && liste_fbr != NULL) {
        BFT_FREE(liste_fac);
        BFT_FREE(liste_fbr);
      }

    }

  }


  /* Libération mémoire */

  BFT_FREE(num_ent_parent);
  BFT_FREE(var_trav);

}


/*----------------------------------------------------------------------------
 * Sortie d'un champ de post traitement défini sur les cellules ou faces
 * d'un maillage en fonction des "writers" associés.
 *
 * Interface Fortran : utiliser PSTEVA (voir cs_post_util.F)
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
                                       *     maillage "parent", 0 si variable
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
)
{
  cs_bool_t  var_parent;
  cs_bool_t  entrelace;

  char  *nom_var = NULL;

  if (*ivarpr == 1)
    var_parent = CS_TRUE;
  else if (*ivarpr == 0)
    var_parent = CS_FALSE;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("L'argument IVARPR du sous-programme PSTEVA doit être\n"
                "égal à 0 ou 1, et non %d.\n"), (int)(*ivarpr));

  if (*ientla == 0)
    entrelace = CS_FALSE;
  else if (*ientla == 1)
    entrelace = CS_TRUE;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("L'argument IENTLA du sous-programme PSTEVA doit être\n"
                "égal à 0 ou 1, et non %d.\n"), (int)(*ientla));


  /* Conversion des chaînes de caractères Fortran en chaînes C */

  nom_var = cs_base_chaine_f_vers_c_cree(nomvar, *lnmvar);


  /* Traitement principal */

  cs_post_ecrit_var(*nummai,
                    nom_var,
                    *idimt,
                    entrelace,
                    var_parent,
                    CS_POST_TYPE_cs_real_t,
                    *ntcabs,
                    *ttcabs,
                    varcel,
                    varfac,
                    varfbr);

  /* Libération des tableaux C temporaires */

  nom_var = cs_base_chaine_f_vers_c_detruit(nom_var);

}


/*----------------------------------------------------------------------------
 * Prise en compte de la renumérotation des faces et faces de bord
 * dans les liens de "parenté" des maillages post.
 *
 * Cette fonction ne doit être appellée qu'une fois, après la renumérotation
 * évuentuelle des faces, pour adapter les maillages post existants.
 * Des nouveaux maillages post seront automatiquement basés sur la
 * "bonne" numérotation, par construction.
 *
 * Interface Fortran :
 *
 * SUBROUTINE PSTRNM
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstrnm, PSTRNM)
(
 void
)
{
  cs_post_renum_faces();
}

/*============================================================================
 * Fonctions publiques
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
)
{
  /* variables locales */

  int    i;

  cs_post_writer_t  *writer = NULL;
  fvm_writer_time_dep_t  dep_temps = FVM_WRITER_FIXED_MESH;


  /* Vérification que le numéro demandé est disponible */

  if (id_writer == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Le numéro de gestionnaire de post traitement demandé\n"
                "doit être < 0 (réservé) ou > 0 (utilisateur).\n"));

  for (i = 0 ; i < cs_glob_post_nbr_writers ; i++) {
    if ((cs_glob_post_writers + i)->id == id_writer)
      bft_error(__FILE__, __LINE__, 0,
                _("Le numéro de gestionnaire de post traitement demandé\n"
                  "(%d) a déjà été affecté.\n"), (int)id_writer);
  }


  /* Redimensionnement du tableau global des writers */

  if (cs_glob_post_nbr_writers == cs_glob_post_nbr_writers_max) {

    if (cs_glob_post_nbr_writers_max == 0)
      cs_glob_post_nbr_writers_max = 4;
    else
      cs_glob_post_nbr_writers_max *= 2;

    BFT_REALLOC(cs_glob_post_writers,
                cs_glob_post_nbr_writers_max,
                cs_post_writer_t);

  }

  cs_glob_post_nbr_writers += 1;


  /* Affectation du writer nouvellement créé à la structure */

  writer = cs_glob_post_writers + cs_glob_post_nbr_writers - 1;

  writer->id = id_writer;
  writer->freq_sortie = frequence;
  writer->ecr_depl = CS_FALSE;
  writer->actif = 0;

  if (ind_mod >= 10) {
    writer->ecr_depl = CS_TRUE;
    ind_mod -= 10;
  }

  if (ind_mod == 1)
    dep_temps = FVM_WRITER_TRANSIENT_COORDS;
  else if (ind_mod >= 2)
    dep_temps = FVM_WRITER_TRANSIENT_CONNECT;

  writer->writer = fvm_writer_init(nom_cas,
                                   nom_rep,
                                   nom_fmt,
                                   opt_fmt,
                                   dep_temps);
}


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
)
{
  /* variables locales */

  cs_post_maillage_t  *maillage_post = NULL;


  /* Ajout et initialisation de la structure de base */

  maillage_post = _cs_post_ajoute_maillage(id_maillage);


  /* Création du maillage et affectation à la structure */

  _cs_post_definit_maillage(maillage_post,
                            nom_maillage,
                            nbr_cel,
                            nbr_fac,
                            nbr_fbr,
                            liste_cel,
                            liste_fac,
                            liste_fbr);
}


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
)
{
  /* variables locales */

  int       i;
  int       indic_glob[3];
  cs_int_t  dec_num_fbr, ind_fac;

  int    indic_loc[3] = {1, 1, 1};  /* Indicateurs 0 à 2 "inversés" par
                                       rapport aux autres pour pouvoir
                                       utiliser un même appel à
                                       MPI_Allreduce(..., MPI_MIN, ...) */

  int         dim_ent = 0;
  cs_bool_t   maj_ind_ent = CS_FALSE;
  fvm_lnum_t  nbr_ent = 0;

  fvm_lnum_t          *num_ent_parent = NULL;
  cs_post_maillage_t  *maillage_post = NULL;


  /* Ajout et initialisation de la structure de base */

  maillage_post = _cs_post_ajoute_maillage(id_maillage);


  /* Affectation du maillage à la structure */

  maillage_post->maillage_ext = maillage_ext;

  if (transferer == CS_TRUE)
    maillage_post->_maillage_ext = maillage_ext;


  /* Calcul du nombre de cellules et/ou de faces */

  dim_ent = fvm_nodal_get_max_entity_dim(maillage_ext);
  nbr_ent = fvm_nodal_get_n_entities(maillage_ext, dim_ent);

  if (dim_ent == 3 && nbr_ent > 0)
    indic_loc[0] = 0;

  else if (dim_ent == 2 && nbr_ent > 0) {

    BFT_MALLOC(num_ent_parent, nbr_ent, cs_int_t);

    fvm_nodal_get_parent_num(maillage_ext, dim_ent, num_ent_parent);

    dec_num_fbr = cs_glob_mesh->n_b_faces;
    for (ind_fac = 0 ; ind_fac < nbr_ent ; ind_fac++) {
      if (num_ent_parent[ind_fac] > dec_num_fbr)
        maillage_post->nbr_fac += 1;
      else
        maillage_post->nbr_fbr += 1;
    }

    BFT_FREE(num_ent_parent);

    if (maillage_post->nbr_fac > 0)
      indic_loc[1] = 0;
    else if (maillage_post->nbr_fbr > 0)
      indic_loc[2] = 0;

  }

  for (i = 0 ; i < 3 ; i++)
    indic_glob[i] = indic_loc[i];

#if defined(_CS_HAVE_MPI)
  if (cs_glob_base_nbr > 1)
    MPI_Allreduce (indic_loc, indic_glob, 3, MPI_INT, MPI_MIN,
                   cs_glob_base_mpi_comm);
#endif


  /* Indicateurs globaux de présence de types de mailles
     on ne met à jour que si le maillage n'est pas totalement vide
     (pour des maillages dépendant du temps, vides à certains instants,
     on veut pouvoir connaître le dernier type de maille contenu
     dans USMPST) */

  for (i = 0 ; i < 3 ; i++) {
    if (indic_glob[i] == 0)
      maj_ind_ent = CS_TRUE;
  }

  if (maj_ind_ent == CS_TRUE) {
    for (i = 0 ; i < 3 ; i++) {
      if (indic_glob[i] == 0)           /* Logique indic_glob 0 à 2 inversée */
        maillage_post->ind_ent[i] = 1;  /* (c.f. remarque ci-dessus) */
      else
        maillage_post->ind_ent[i] = 0;
    }
  }

  /* Indicateurs de modification min et max inversés initialement,
     seront recalculés lors des associations maillages - post-traitements */

  maillage_post->ind_mod_min = FVM_WRITER_TRANSIENT_CONNECT;
  maillage_post->ind_mod_max = FVM_WRITER_FIXED_MESH;

}


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
)
{
  /* variables locales */

  int    indref, j;

  cs_post_maillage_t  *maillage_post = NULL;
  cs_post_maillage_t  *maillage_ref = NULL;


  /* Vérifications initiales */

  indref = _cs_post_ind_maillage(id_maillage);
  maillage_ref = cs_glob_post_maillages + indref;

  if (maillage_ref->alias > -1)
    bft_error(__FILE__, __LINE__, 0,
              _("Le maillage %d ne peut être un alias du maillage %d,\n"
                "qui est lui-même déjà un alias du maillage %d.\n"),
              (int)id_alias, (int)id_maillage,
              (int)((cs_glob_post_maillages + maillage_ref->alias)->id));

  /* Ajout et initialisation de la structure de base */

  maillage_post = _cs_post_ajoute_maillage(id_alias);

  /* On réactualise maillage_ref, car on a pu déplacer l'adresse
     de cs_glob_pos_maillages en le réallouant */

  maillage_ref = cs_glob_post_maillages + indref;

  /* Liens avec le maillage de référence */

  maillage_post->alias = indref;

  maillage_post->maillage_ext = maillage_ref->maillage_ext;

  maillage_post->ind_mod_min = maillage_ref->ind_mod_min;
  maillage_post->ind_mod_max = maillage_ref->ind_mod_max;

  for (j = 0 ; j < 3 ; j++)
    maillage_post->ind_ent[j] = maillage_ref->ind_ent[j];

  maillage_post->nbr_fac = maillage_ref->nbr_fac;
  maillage_post->nbr_fbr = maillage_ref->nbr_fbr;

}


/*----------------------------------------------------------------------------
 * Vérifie l'existence d'un "writer" associé à un numéro donné.
 *----------------------------------------------------------------------------*/

cs_bool_t cs_post_existe_writer
(
 const cs_int_t   numwri        /* --> numéro du writer associé               */
)
{
  /* variables locales */

  int indwri;
  cs_post_writer_t  *writer = NULL;


  /* Recherche du writer demandé */

  for (indwri = 0 ; indwri < cs_glob_post_nbr_writers ; indwri++) {
    writer = cs_glob_post_writers + indwri;
    if (writer->id == numwri)
      return CS_TRUE;
  }

  return CS_FALSE;
}


/*----------------------------------------------------------------------------
 * Vérifie l'existence d'un maillage de post traitement associé à un
 * numéro donné.
 *----------------------------------------------------------------------------*/

cs_bool_t cs_post_existe_maillage
(
 const cs_int_t   nummai        /* --> numéro du maillage externe associé     */
)
{
  /* variables locales */

  int indmai;
  cs_post_maillage_t  *maillage_post = NULL;


  /* Recherche du maillage demandé */

  for (indmai = 0 ; indmai < cs_glob_post_nbr_maillages ; indmai++) {
    maillage_post = cs_glob_post_maillages + indmai;
    if (maillage_post->id == nummai)
      return CS_TRUE;
  }

  return CS_FALSE;
}


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
)
{
  /* variables locales */

  int i, indmai;
  char  *nom_maillage = NULL;
  cs_post_maillage_t  *maillage_post = NULL;
  cs_post_writer_t    *writer = NULL;


  /* Récupération de la structure de base
     (sortie si l'on n'est pas propriétaire du maillage) */

  indmai = _cs_post_ind_maillage(id_maillage);
  maillage_post = cs_glob_post_maillages + indmai;

  if (maillage_post->_maillage_ext == NULL)
    return;


  /* Remplacement de la structure de base */

  BFT_MALLOC(nom_maillage,
             strlen(fvm_nodal_get_name(maillage_post->maillage_ext)) + 1,
             char);
  strcpy(nom_maillage, fvm_nodal_get_name(maillage_post->maillage_ext));

  fvm_nodal_destroy(maillage_post->_maillage_ext);
  maillage_post->maillage_ext = NULL;

  _cs_post_definit_maillage(maillage_post,
                            nom_maillage,
                            nbr_cel,
                            nbr_fac,
                            nbr_fbr,
                            liste_cel,
                            liste_fac,
                            liste_fbr);

  BFT_FREE(nom_maillage);


  /* Mise à jour des alias éventuels */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {
    if ((cs_glob_post_maillages + i)->alias == indmai)
      (cs_glob_post_maillages + i)->maillage_ext
        = maillage_post->maillage_ext;
  }

  /* Découpage des polygones ou polyèdres en éléments simples */

  for (i = 0 ; i < maillage_post->nbr_writers ; i++) {

    writer = cs_glob_post_writers + maillage_post->ind_writer[i];
    _cs_post_divise_poly(maillage_post, writer);

  }

}


/*----------------------------------------------------------------------------
 * Récupération du prochain numéro de maillage standard ou développeur
 * disponible (basé sur le plus petit numéro négatif présent -1).
 *----------------------------------------------------------------------------*/

cs_int_t cs_post_ret_num_maillage_libre
(
 void
)
{
  return (cs_glob_post_num_maillage_min - 1);
}


/*----------------------------------------------------------------------------
 * Association d'un "writer" à un maillage pour le post traitement.
 *----------------------------------------------------------------------------*/

void cs_post_associe
(
 const cs_int_t   id_maillage,  /* --> numéro du maillage externe associé     */
 const cs_int_t   id_writer     /* --> numéro du writer                       */
)
{
  int  i;
  int  indmai, indgep;
  fvm_writer_time_dep_t ind_mod;

  cs_post_maillage_t  *maillage_post = NULL;
  cs_post_writer_t  *writer = NULL;


  /* Recherche du maillage et writer demandés */

  indmai = _cs_post_ind_maillage(id_maillage);
  indgep = _cs_post_ind_writer(id_writer);

  maillage_post = cs_glob_post_maillages + indmai;

  /* On vérifie que le writer n'est pas déjà associé */

  for (i = 0 ; i < maillage_post->nbr_writers ; i++) {
    if (maillage_post->ind_writer[i] == indgep)
      break;
  }

  /* Si le writer n'est pas déjà associé, on l'associe */

  if (i >= maillage_post->nbr_writers) {

    maillage_post->nbr_writers += 1;
    BFT_REALLOC(maillage_post->ind_writer,
                maillage_post->nbr_writers,
                cs_int_t);

    maillage_post->ind_writer[maillage_post->nbr_writers - 1] = indgep;
    maillage_post->nt_ecr = - 1;

    /* Mise à jour de la structure */

    writer = cs_glob_post_writers + indgep;
    ind_mod = fvm_writer_get_time_dep(writer->writer);

    if (ind_mod < maillage_post->ind_mod_min)
      maillage_post->ind_mod_min = ind_mod;
    if (ind_mod > maillage_post->ind_mod_max)
      maillage_post->ind_mod_max = ind_mod;

    _cs_post_ind_mod_alias(indmai);

    /* Si l'on doit calculer le champ déplacement des sommets, on devra
       sauvegarder les coordonnées initiales des sommets */

    if (   cs_glob_post_deformable == CS_FALSE
        && cs_glob_post_coo_som_ini == NULL
        && writer->ecr_depl == CS_TRUE) {

      cs_mesh_t *maillage = cs_glob_mesh;

      if (maillage->n_vertices > 0) {
        BFT_MALLOC(cs_glob_post_coo_som_ini,
                   maillage->n_vertices * 3,
                   cs_real_t);
        memcpy(cs_glob_post_coo_som_ini,
               maillage->vtx_coord,
               maillage->n_vertices * 3 * sizeof(cs_real_t));
      }

      cs_glob_post_deformable = CS_TRUE;

    }

    /* Découpage des polygones ou polyèdres en éléments simples */

    _cs_post_divise_poly(maillage_post, writer);

  }

}



/*----------------------------------------------------------------------------
 * Mise à jour de l'indicateur "actif" ou "inactif" des "writers" en
 * fonction du pas de temps et de leur fréquence de sortie par défaut.
 *----------------------------------------------------------------------------*/

void cs_post_activer_selon_defaut
(
 const cs_int_t   nt_cur_abs      /* --> numéro de pas de temps courant       */
)
{
  int  i;
  cs_post_writer_t  *writer;

  for (i = 0 ; i < cs_glob_post_nbr_writers ; i++) {

    writer = cs_glob_post_writers + i;

    if (writer->freq_sortie > 0) {
      if (nt_cur_abs % (writer->freq_sortie) == 0)
        writer->actif = 1;
      else
        writer->actif = 0;
    }
    else
      writer->actif = 0;

  }
}


/*----------------------------------------------------------------------------
 * Forcer de l'indicateur "actif" ou "inactif" d'un "writers" spécifique
 * ou de l'ensemble des "writers" pour le pas de temps en cours.
 *----------------------------------------------------------------------------*/

void cs_post_activer_writer
(
 const cs_int_t   id_writer,    /* --> numéro du writer,ou 0 pour forcer
                                 *     simultanément tous les writers         */
 const cs_int_t   activer       /* --> 0 pour désactiver, 1 pour activer      */
)
{
  int i;
  cs_post_writer_t  *writer;

  if (id_writer != 0) {
    i = _cs_post_ind_writer(id_writer);
    writer = cs_glob_post_writers + i;
    writer->actif = (activer > 0) ? 1 : 0;
  }
  else {
    for (i = 0 ; i < cs_glob_post_nbr_writers ; i++) {
      writer = cs_glob_post_writers + i;
      writer->actif = (activer > 0) ? 1 : 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Get the writer associated to a writer_id.
 *
 * writer_id       -->  id of the writer in cs_glob_post_writers
 *
 * Returns:
 *  a pointer to a fvm_writer_t structure
 *----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(cs_int_t   writer_id)
{
  int  id;
  const cs_post_writer_t  *writer = NULL;

  id = _cs_post_ind_writer(writer_id);
  writer = cs_glob_post_writers + id;

  return writer->writer;
}

/*----------------------------------------------------------------------------
 * Ecriture des maillages de post traitement en fonction des "writers"
 * associés.
 *----------------------------------------------------------------------------*/

void cs_post_ecrit_maillages
(
 const cs_int_t   nt_cur_abs,         /* --> numéro de pas de temps courant   */
 const cs_real_t  t_cur_abs           /* --> valeur du temps physique associé */
)
{
  int  i;
  cs_post_maillage_t  *maillage_post;

  /* Boucles sur les maillages et "writers" */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    _cs_post_ecrit_maillage(maillage_post,
                            nt_cur_abs,
                            t_cur_abs);

  }

}


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
)
{
  cs_int_t  i;
  int       indmai;


  fvm_interlace_t      interlace;
  fvm_datatype_t       datatype;

  size_t       dec_ptr = 0;
  int          nbr_listes_parents = 0;
  fvm_lnum_t   dec_num_parent[2]  = {0, 0};
  cs_real_t   *var_tmp = NULL;
  cs_post_maillage_t  *maillage_post = NULL;
  cs_post_writer_t    *writer = NULL;

  const void  *var_ptr[2*9] = {NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL};

  /* Initialisations */

  indmai = _cs_post_ind_maillage(id_maillage);
  maillage_post = cs_glob_post_maillages + indmai;

  if (entrelace == CS_TRUE)
    interlace = FVM_INTERLACE;
  else
    interlace = FVM_NO_INTERLACE;

  datatype =  _cs_post_cnv_datatype(var_type);


  /* Affectation du tableau approprié à FVM pour la sortie */

  /* Cas des cellules */
  /*------------------*/

  if (maillage_post->ind_ent[CS_POST_SUPPORT_CEL] == 1) {

    if (var_parent == CS_TRUE) {
      nbr_listes_parents = 1;
      dec_num_parent[0] = 0;
    }
    else
      nbr_listes_parents = 0;

    var_ptr[0] = var_cel;
    if (entrelace == CS_FALSE) {
      if (var_parent == CS_TRUE)
        dec_ptr = cs_glob_mesh->n_cells_with_ghosts;
      else
        dec_ptr = fvm_nodal_get_n_entities(maillage_post->maillage_ext, 3);
      dec_ptr *= fvm_datatype_size[datatype];
      for (i = 1 ; i < dim_var ; i++)
        var_ptr[i] = ((const char *)var_cel) + i*dec_ptr;
    }
  }

  /* Cas des faces */
  /*---------------*/

  else if (   maillage_post->ind_ent[CS_POST_SUPPORT_FAC_INT] == 1
           || maillage_post->ind_ent[CS_POST_SUPPORT_FAC_BRD] == 1) {

    /* En cas d'indirection, il suffit de positionner les pointeurs */

    if (var_parent == CS_TRUE) {

      nbr_listes_parents = 2;
      dec_num_parent[0] = 0;
      dec_num_parent[1] = cs_glob_mesh->n_b_faces;

      if (maillage_post->ind_ent[CS_POST_SUPPORT_FAC_BRD] == 1) {
        if (entrelace == CS_FALSE) {
          dec_ptr = cs_glob_mesh->n_b_faces * fvm_datatype_size[datatype];
          for (i = 0 ; i < dim_var ; i++)
            var_ptr[i] = ((const char *)var_fbr) + i*dec_ptr;
        }
        else
          var_ptr[0] = var_fbr;
      }

      if (maillage_post->ind_ent[CS_POST_SUPPORT_FAC_INT] == 1) {
        if (entrelace == CS_FALSE) {
          dec_ptr = cs_glob_mesh->n_i_faces * fvm_datatype_size[datatype];
          for (i = 0 ; i < dim_var ; i++)
            var_ptr[dim_var + i] = ((const char *)var_fac) + i*dec_ptr;
        }
        else
          var_ptr[1] = var_fac;
      }

    }

    /* Sans indirection, on doit repasser d'une variable définie sur deux
       listes de faces à une variable définie sur une liste */

    else {

      nbr_listes_parents = 0;

      if (maillage_post->ind_ent[CS_POST_SUPPORT_FAC_BRD] == 1) {

        /* Cas où la variable est définie à la fois sur des faces de
           bord et des faces internes : on doit repasser à une liste
           unique, étant donné que l'on n'utilise pas l'indirection */

        if (maillage_post->ind_ent[CS_POST_SUPPORT_FAC_INT] == 1) {

          BFT_MALLOC(var_tmp,
                     (   maillage_post->nbr_fac
                      +  maillage_post->nbr_fbr) * dim_var,
                     cs_real_t);

          _cs_post_assmb_var_faces(maillage_post->maillage_ext,
                                   maillage_post->nbr_fac,
                                   maillage_post->nbr_fbr,
                                   dim_var,
                                   interlace,
                                   var_fac,
                                   var_fbr,
                                   var_tmp);

          interlace = FVM_NO_INTERLACE;

          dec_ptr = fvm_datatype_size[datatype] * (  maillage_post->nbr_fac
                                                   + maillage_post->nbr_fbr);

          for (i = 0 ; i < dim_var ; i++)
            var_ptr[i] = ((char *)var_tmp) + i*dec_ptr;

        }

        /* Cas où l'on a que des faces de bord */

        else {

          if (entrelace == CS_FALSE) {
            dec_ptr = fvm_datatype_size[datatype] * maillage_post->nbr_fbr;
            for (i = 0 ; i < dim_var ; i++)
              var_ptr[i] = ((const char *)var_fbr) + i*dec_ptr;
          }
          else
            var_ptr[0] = var_fbr;
        }

      }

      /* Cas où l'on a que des faces internes */

      else if (maillage_post->ind_ent[CS_POST_SUPPORT_FAC_INT] == 1) {

        if (entrelace == CS_FALSE) {
          dec_ptr = fvm_datatype_size[datatype] * maillage_post->nbr_fac;
          for (i = 0 ; i < dim_var ; i++)
            var_ptr[i] = ((const char *)var_fac) + i*dec_ptr;
        }
        else
          var_ptr[0] = var_fac;
      }

    }

  }


  /* Sortie effective : boucle sur les writers */
  /*-------------------------------------------*/

  for (i = 0 ; i < maillage_post->nbr_writers ; i++) {

    writer = cs_glob_post_writers + maillage_post->ind_writer[i];

    if (writer->actif == 1)
      fvm_writer_export_field(writer->writer,
                              maillage_post->maillage_ext,
                              nom_var,
                              FVM_WRITER_PER_ELEMENT,
                              dim_var,
                              interlace,
                              nbr_listes_parents,
                              dec_num_parent,
                              datatype,
                              (int)nt_cur_abs,
                              (double)t_cur_abs,
                              (const void * *)var_ptr);

  }

  /* Libération mémoire (si faces internes et de bord simultanées) */

  if (var_tmp != NULL)
    BFT_FREE(var_tmp);

}


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
)
{
  cs_int_t  i;
  int       indmai;


  cs_post_maillage_t  *maillage_post;
  cs_post_writer_t    *writer;
  fvm_interlace_t      interlace;
  fvm_datatype_t       datatype;

  size_t       dec_ptr = 0;
  int          nbr_listes_parents = 0;
  fvm_lnum_t   dec_num_parent[1]  = {0};

  const void  *var_ptr[9] = {NULL, NULL, NULL,
                             NULL, NULL, NULL,
                             NULL, NULL, NULL};


  /* Initialisations */

  indmai = _cs_post_ind_maillage(id_maillage);
  maillage_post = cs_glob_post_maillages + indmai;

  if (entrelace == CS_TRUE)
    interlace = FVM_INTERLACE;
  else
    interlace = FVM_NO_INTERLACE;

  assert(   sizeof(cs_real_t) == sizeof(double)
         || sizeof(cs_real_t) == sizeof(float));

  datatype =  _cs_post_cnv_datatype(var_type);


  /* Affectation du tableau approprié à FVM pour la sortie */

  if (var_parent == CS_TRUE)
    nbr_listes_parents = 1;
  else
    nbr_listes_parents = 0;

  var_ptr[0] = var_som;
  if (entrelace == CS_FALSE) {
    if (var_parent == CS_TRUE)
      dec_ptr = cs_glob_mesh->n_vertices;
    else
      dec_ptr =   fvm_nodal_get_n_entities(maillage_post->maillage_ext, 0)
                * fvm_datatype_size[datatype];
    for (i = 1 ; i < dim_var ; i++)
      var_ptr[i] = ((const char *)var_som) + i*dec_ptr;
  }


  /* Sortie effective : boucle sur les writers */
  /*-------------------------------------------*/

  for (i = 0 ; i < maillage_post->nbr_writers ; i++) {

    writer = cs_glob_post_writers + maillage_post->ind_writer[i];

    if (writer->actif == 1)
      fvm_writer_export_field(writer->writer,
                              maillage_post->maillage_ext,
                              nom_var,
                              FVM_WRITER_PER_NODE,
                              dim_var,
                              interlace,
                              nbr_listes_parents,
                              dec_num_parent,
                              datatype,
                              (int)nt_cur_abs,
                              (double)t_cur_abs,
                              (const void * *)var_ptr);

  }

}


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
 void
)
{
  int       i;
  cs_int_t  ifac;
  cs_int_t  nbr_ent;

  cs_int_t  *renum_ent_parent = NULL;

  cs_bool_t  a_traiter = CS_FALSE;

  cs_post_maillage_t   *maillage_post;
  const cs_mesh_t  *maillage = cs_glob_mesh;


  /* Boucles sur les maillages */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    if (   maillage_post->ind_ent[CS_POST_SUPPORT_FAC_INT] > 0
        || maillage_post->ind_ent[CS_POST_SUPPORT_FAC_BRD] > 0) {
      a_traiter = CS_TRUE;
    }

  }

  if (a_traiter == CS_TRUE) {

    /* Préparation de la renumérotation */

    nbr_ent = maillage->n_i_faces + maillage->n_b_faces;

    BFT_MALLOC(renum_ent_parent, nbr_ent, cs_int_t);

    if (maillage->init_b_face_num == NULL) {
      for (ifac = 0 ; ifac < maillage->n_b_faces ; ifac++)
        renum_ent_parent[ifac] = ifac + 1;
    }
    else {
      for (ifac = 0 ; ifac < maillage->n_b_faces ; ifac++)
        renum_ent_parent[maillage->init_b_face_num[ifac] - 1] = ifac + 1;
    }

    if (maillage->init_i_face_num == NULL) {
      for (ifac = 0, i = maillage->n_b_faces ;
           ifac < maillage->n_i_faces ;
           ifac++, i++)
        renum_ent_parent[maillage->n_b_faces + ifac]
          = maillage->n_b_faces + ifac + 1;
    }
    else {
      for (ifac = 0, i = maillage->n_b_faces ;
           ifac < maillage->n_i_faces ;
           ifac++, i++)
        renum_ent_parent[  maillage->n_b_faces
                         + maillage->init_i_face_num[ifac] - 1]
          = maillage->n_b_faces + ifac + 1;
    }

    /* Modification effective */

    for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

      maillage_post = cs_glob_post_maillages + i;

      if (   maillage_post->_maillage_ext != NULL
          && (   maillage_post->ind_ent[CS_POST_SUPPORT_FAC_INT] > 0
              || maillage_post->ind_ent[CS_POST_SUPPORT_FAC_BRD] > 0)) {

        fvm_nodal_change_parent_num(maillage_post->_maillage_ext,
                                    renum_ent_parent,
                                    2);

      }

    }

    BFT_FREE(renum_ent_parent);

  }

}


/*----------------------------------------------------------------------------
 * Destruction des structures associées aux post traitements
 *----------------------------------------------------------------------------*/

void cs_post_detruit
(
 void
)
{
  int i;
  cs_post_maillage_t  *maillage_post = NULL;

  /* Chronométrages */

  for (i = 0 ; i < cs_glob_post_nbr_writers ; i++) {
    double m_wtime = 0.0, m_cpu_time = 0.0, c_wtime = 0.0, c_cpu_time = 0.;
    fvm_writer_get_times((cs_glob_post_writers + i)->writer,
                         &m_wtime, &m_cpu_time, &c_wtime, &c_cpu_time);
    bft_printf(_("\n"
                 "Bilan des écritures de \"%s\" (%s) :\n\n"
                 "  Temps CPU pour les maillages :    %12.3f\n"
                 "  Temps CPU pour les champs :       %12.3f\n\n"
                 "  Temps écoulé pour les maillages : %12.3f\n"
                 "  Temps écoulé pour les champs :    %12.3f\n"),
               fvm_writer_get_name((cs_glob_post_writers + i)->writer),
               fvm_writer_get_format((cs_glob_post_writers + i)->writer),
               m_cpu_time, c_cpu_time, m_wtime, c_wtime);
  }

  /* Coordonnées initiales si maillage déformable) */

  if (cs_glob_post_coo_som_ini != NULL)
    BFT_FREE(cs_glob_post_coo_som_ini);

  /* Maillages externes */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {
    maillage_post = cs_glob_post_maillages + i;
    if (maillage_post->_maillage_ext != NULL)
      fvm_nodal_destroy(maillage_post->_maillage_ext);
    BFT_FREE(maillage_post->ind_writer);
  }

  BFT_FREE(cs_glob_post_maillages);

  cs_glob_post_num_maillage_min = -2;
  cs_glob_post_nbr_maillages = 0;
  cs_glob_post_nbr_maillages_max = 0;

  /* Writers */

  for (i = 0 ; i < cs_glob_post_nbr_writers ; i++)
    fvm_writer_finalize((cs_glob_post_writers + i)->writer);

  BFT_FREE(cs_glob_post_writers);

  cs_glob_post_nbr_writers = 0;
  cs_glob_post_nbr_writers_max = 0;

  /* Traitements enregistrés si nécessaire */

  if (cs_glob_post_nbr_var_tp_max > 0) {
    BFT_FREE(cs_glob_post_f_var_tp);
    BFT_FREE(cs_glob_post_i_var_tp);
  }

}


/*----------------------------------------------------------------------------
 * Initialisation du "writer" du post-traitement principal
 *----------------------------------------------------------------------------*/

void cs_post_init_pcp_writer
(
 void
)
{
  /* Valeurs par défaut */

  cs_int_t  indic_vol = -1, indic_brd = -1, indic_syr = -1;
  cs_int_t  indic_mod = -1;
  char  fmtchr[32 + 1] = "";
  char  optchr[96 + 1] = "";
  cs_int_t  ntchr = -1;

  const char  nomcas[] = "chr";
  const char  nomrep_ens[] = "chr.ensight";
  const char  nomrep_def[] = ".";
  const char *nomrep = NULL;

  const cs_int_t  id_writer = -1; /* Numéro du writer par défaut */

  /* Récupération des paramètres contenus dans les COMMONS Fortran */

  CS_PROCF(inipst, INIPST)(&indic_vol,
                           &indic_brd,
                           &indic_syr,
                           &indic_mod,
                           &ntchr,
                           fmtchr,
                           optchr);

  fmtchr[32] = '\0';
  optchr[96] = '\0';

  if (indic_vol == 0 && indic_brd == 0 && indic_syr == 0)
    return;

  /* Création du writer par défaut */

  if (fmtchr[0] == 'e' || fmtchr[0] == 'E')
    nomrep = nomrep_ens;
  else
    nomrep = nomrep_def;

  cs_post_ajoute_writer(id_writer,
                        nomcas,
                        nomrep,
                        fmtchr,
                        optchr,
                        indic_mod,
                        ntchr);

}


/*----------------------------------------------------------------------------
 * Initialisation des maillages du post-traitement principal
 *----------------------------------------------------------------------------*/

void cs_post_init_pcp_maillages
(
 void
)
{
  /* Valeurs par défaut */

  cs_int_t  indic_vol = -1, indic_brd = -1, indic_syr = -1;
  cs_int_t  indic_mod = -1;
  char  fmtchr[32 + 1] = "";
  char  optchr[96 + 1] = "";
  cs_int_t  ntchr = -1;

  cs_int_t  id_maillage = -1;

  const cs_int_t  id_writer = -1; /* Numéro du writer par défaut */

  /* Récupération des paramètres contenus dans les COMMONS Fortran */

  CS_PROCF(inipst, INIPST)(&indic_vol,
                           &indic_brd,
                           &indic_syr,
                           &indic_mod,
                           &ntchr,
                           fmtchr,
                           optchr);

  fmtchr[32] = '\0';
  optchr[96] = '\0';

  /* Définition des maillages de post-traitement */

  if (cs_glob_mesh->n_i_faces > 0 || cs_glob_mesh->n_b_faces > 0) {

    /*
      Si on dispose de la connectivité faces -> sommets, on
      peut reconstruire la connectivité nodale pour le post
      traitement (mécanisme usuel).
    */

    if (indic_vol > 0) { /* Maillage volumique */

      id_maillage = -1; /* Numéro de maillage réservé */

      cs_post_ajoute_maillage(id_maillage,
                              _("Volume fluide"),
                              cs_glob_mesh->n_cells,
                              0,
                              0,
                              NULL,
                              NULL,
                              NULL);

      cs_post_associe(id_maillage, id_writer);

    }

    if (indic_brd > 0) { /* Maillage de peau */

      id_maillage = -2;  /* Numéro de maillage réservé */

      cs_post_ajoute_maillage(id_maillage,
                              _("Bord"),
                              0,
                              0,
                              cs_glob_mesh->n_b_faces,
                              NULL,
                              NULL,
                              NULL);

      cs_post_associe(id_maillage, id_writer);

    }

  } /* Fin if cs_glob_mesh->n_i_faces > 0 || cs_glob_mesh->n_b_faces > 0 */

  /*
    Si on ne dispose pas de la connectivité faces -> sommets, on
    ne peut pas reconstruire la connectivité nodale, donc on doit
    l'obtenir par un autre moyen.

    Ceci ne doit se produire que lorsqu'on a lu directement certains
    maillages au format solcom, dans quel cas la connectivité
    nodale a déja été lue et affectée à un maillage post
    (voir LETGEO et cs_maillage_solcom_lit).
  */

  else if (indic_vol > 0) {

    id_maillage = -1;

    if (cs_post_existe_maillage(id_maillage))
      cs_post_associe(id_maillage, id_writer);

  }

}


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
)
{
  /* Redimensionnement du tableau des traitements enregistrés si nécessaire */

  if (cs_glob_post_nbr_var_tp <= cs_glob_post_nbr_var_tp_max) {
    if (cs_glob_post_nbr_var_tp_max == 0)
      cs_glob_post_nbr_var_tp_max = 8;
    else
      cs_glob_post_nbr_var_tp_max *= 2;
    BFT_REALLOC(cs_glob_post_f_var_tp,
                cs_glob_post_nbr_var_tp_max,
                cs_post_var_temporelle_t *);
    BFT_REALLOC(cs_glob_post_i_var_tp, cs_glob_post_nbr_var_tp_max, cs_int_t);
  }

  /* Ajout du traitement */

  cs_glob_post_f_var_tp[cs_glob_post_nbr_var_tp] = fonction;
  cs_glob_post_i_var_tp[cs_glob_post_nbr_var_tp] = id_instance;

  cs_glob_post_nbr_var_tp += 1;

}


/*============================================================================
 * Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Conversion d'un type de données cs_post_type_t en type fvm_datatype_t.
 *----------------------------------------------------------------------------*/

static fvm_datatype_t _cs_post_cnv_datatype
(
 cs_post_type_t type_cs
)
{
  fvm_datatype_t type_fvm = FVM_DATATYPE_NULL;

  switch(type_cs) {

  case CS_POST_TYPE_cs_int_t:
    if (sizeof(cs_int_t) == 4)
      type_fvm = FVM_INT32;
    else if (sizeof(cs_int_t) == 8)
      type_fvm = FVM_INT64;
    break;

  case CS_POST_TYPE_cs_real_t:
    if (sizeof(cs_real_t) == sizeof(double))
      type_fvm = FVM_DOUBLE;
    else if (sizeof(cs_real_t) == sizeof(float))
      type_fvm = FVM_FLOAT;
    break;

  case CS_POST_TYPE_int:
    if (sizeof(int) == 4)
      type_fvm = FVM_INT32;
    else if (sizeof(int) == 8)
      type_fvm = FVM_INT64;
    break;

  case CS_POST_TYPE_float:
    type_fvm = FVM_FLOAT;
    break;

  case CS_POST_TYPE_double:
    type_fvm = FVM_DOUBLE;
    break;

  default:
    assert(0);
  }

  return type_fvm;
}


/*----------------------------------------------------------------------------
 * Recherche de l'indice d'un writer associé à un numéro donné.
 *----------------------------------------------------------------------------*/

static int _cs_post_ind_writer
(
 const cs_int_t   id_writer     /* --> numéro du writer                       */
)
{
  cs_int_t  indgep;

  cs_post_writer_t  *writer = NULL;


  /* Recherche du writer demandé */

  for (indgep = 0 ; indgep < cs_glob_post_nbr_writers ; indgep++) {
    writer = cs_glob_post_writers + indgep;
    if (writer->id == id_writer)
      break;
  }
  if (indgep >= cs_glob_post_nbr_writers)
    bft_error(__FILE__, __LINE__, 0,
              _("Le gestionnaire de post traitement numéro %d demandé\n"
                "n'est pas défini.\n"), (int)(id_writer));

  return indgep;

}


/*----------------------------------------------------------------------------
 * Recherche de l'indice d'un maillage de post traitement associé à un
 * numéro donné.
 *----------------------------------------------------------------------------*/

static int _cs_post_ind_maillage
(
 const cs_int_t   nummai        /* --> numéro du maillage externe associé     */
)
{
  int indmai;
  cs_post_maillage_t  *maillage_post = NULL;


  /* Recherche du maillage demandé */

  for (indmai = 0 ; indmai < cs_glob_post_nbr_maillages ; indmai++) {
    maillage_post = cs_glob_post_maillages + indmai;
    if (maillage_post->id == nummai)
      break;
  }
  if (indmai >= cs_glob_post_nbr_maillages)
    bft_error(__FILE__, __LINE__, 0,
              _("Le maillage de post traitement numéro %d demandé\n"
                "n'est pas défini.\n"), (int)nummai);

  return indmai;

}


/*----------------------------------------------------------------------------
 * Ajout d'un maillage de post traitement et initialisation de base,
 * et renvoi d'un pointeur sur la structure associée
 *----------------------------------------------------------------------------*/

static cs_post_maillage_t * _cs_post_ajoute_maillage
(
 const cs_int_t          id_maillage   /* --> numéro du maillage  demandé     */
)
{
  /* variables locales */

  int    i, j;

  cs_post_maillage_t  *maillage_post = NULL;


  /* Vérification que le numéro demandé est disponible */

  if (id_maillage == 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Le numéro de maillage de post traitement demandé\n"
                  "doit être < 0 (réservé) ou > 0 (utilisateur).\n"));

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {
    if ((cs_glob_post_maillages + i)->id == id_maillage)
      bft_error(__FILE__, __LINE__, 0,
                _("Le numéro de maillage de post traitement demandé\n"
                  "(%d) a déjà été affecté.\n"), (int)id_maillage);
  }


  /* Redimensionnement du tableau global des maillages externes */

  if (cs_glob_post_nbr_maillages == cs_glob_post_nbr_maillages_max) {

    if (cs_glob_post_nbr_maillages_max == 0)
      cs_glob_post_nbr_maillages_max = 8;
    else
      cs_glob_post_nbr_maillages_max *= 2;

    BFT_REALLOC(cs_glob_post_maillages,
                cs_glob_post_nbr_maillages_max,
                cs_post_maillage_t);

  }

  cs_glob_post_nbr_maillages += 1;

  if (id_maillage < cs_glob_post_num_maillage_min)
    cs_glob_post_num_maillage_min = id_maillage;

  /* Affectation du maillage nouvellement crée à la structure */

  maillage_post = cs_glob_post_maillages + cs_glob_post_nbr_maillages - 1;

  maillage_post->id = id_maillage;
  maillage_post->alias = -1;

  maillage_post->nbr_writers = 0;
  maillage_post->ind_writer = NULL;

  maillage_post->nt_ecr = -1;

  for (j = 0 ; j < 3 ; j++)
    maillage_post->ind_ent[j] = 0;

  maillage_post->nbr_fac = 0;
  maillage_post->nbr_fbr = 0;

  maillage_post->maillage_ext = NULL;
  maillage_post->_maillage_ext = NULL;

  /* Indicateurs de modification min et max inversés initialement,
     seront recalculés lors des associations maillages - post-traitements */

  maillage_post->ind_mod_min = FVM_WRITER_TRANSIENT_CONNECT;
  maillage_post->ind_mod_max = FVM_WRITER_FIXED_MESH;

  return maillage_post;
}


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
 * liste des faces de bord associée n'est donc pas nécessaire.
 *----------------------------------------------------------------------------*/

static void _cs_post_definit_maillage
(
 cs_post_maillage_t  *const maillage_post, /* <-> maillage post à compléter   */
 const char          *const nom_maillage,  /* --> nom du maillage externe     */
 const cs_int_t             nbr_cel,       /* --> nombre de cellules          */
 const cs_int_t             nbr_fac,       /* --> nombre de faces internes    */
 const cs_int_t             nbr_fbr,       /* --> nombre de faces de bord     */
       cs_int_t             liste_cel[],   /* <-> liste des cellules          */
       cs_int_t             liste_fac[],   /* <-> liste des faces internes    */
       cs_int_t             liste_fbr[]    /* <-> liste des faces de bord     */
)
{
  /* variables locales */

  int    i;
  int    indic_glob[5];

  int    indic_loc[5] = {1, 1, 1, 0, 0};  /* Indicateurs 0 à 2 "inversés" par
                                             rapport aux autres pour pouvoir
                                             utiliser un même appel à
                                             MPI_Allreduce(..., MPI_MIN, ...) */

  fvm_nodal_t  *maillage_ext = NULL;
  cs_bool_t     maj_ind_ent = CS_FALSE;


  /* Indicateurs:
     0 : 0 si cellules, 1 si pas de cellules,
     1 : 0 si faces internes, 1 si pas de faces internes,
     2 : 0 si faces de bord, 1 si pas de faces de bord,
     3 : 1 si toutes les cellules sont sélectionnées,
     4 : 1 si toutes les faces de bord et aucune face interne sélectionnées */

  if (nbr_cel > 0)
    indic_loc[0] = 0;
  else {
    if (nbr_fac > 0)
      indic_loc[1] = 0;
    if (nbr_fbr > 0)
      indic_loc[2] = 0;
  }

  if (nbr_cel >= cs_glob_mesh->n_cells)
    indic_loc[3] = 1;
  else
    indic_loc[3] = 0;

  if (   nbr_fbr >= cs_glob_mesh->n_b_faces
      && nbr_fac == 0)
    indic_loc[4] = 1;
  else
    indic_loc[4] = 0;

  for (i = 0 ; i < 5 ; i++)
    indic_glob[i] = indic_loc[i];

#if defined(_CS_HAVE_MPI)
  if (cs_glob_base_nbr > 1)
    MPI_Allreduce (indic_loc, indic_glob, 5, MPI_INT, MPI_MIN,
                   cs_glob_base_mpi_comm);
#endif

  /* Création de la structure associée */

  if (indic_glob[0] == 0) {

    if (indic_glob[3] == 1)
      maillage_ext = cs_maillage_extrait_cel_nodal(cs_glob_mesh,
                                                   nom_maillage,
                                                   cs_glob_mesh->n_cells,
                                                   NULL);
    else
      maillage_ext = cs_maillage_extrait_cel_nodal(cs_glob_mesh,
                                                   nom_maillage,
                                                   nbr_cel,
                                                   liste_cel);

  }
  else {

    if (indic_glob[4] == 1)
      maillage_ext = cs_maillage_extrait_fac_nodal(cs_glob_mesh,
                                                   nom_maillage,
                                                   0,
                                                   cs_glob_mesh->n_b_faces,
                                                   NULL,
                                                   NULL);
    else
      maillage_ext = cs_maillage_extrait_fac_nodal(cs_glob_mesh,
                                                   nom_maillage,
                                                   nbr_fac,
                                                   nbr_fbr,
                                                   liste_fac,
                                                   liste_fbr);

  }

  /* Indicateurs globaux de présence de types de mailles ;
     on ne met à jour que si le maillage n'est pas totalement vide
     (pour des maillages dépendant du temps, vides à certains instants,
     on veut pouvoir connaître le dernier type de maille contenu
     dans USMPST) */

  for (i = 0 ; i < 3 ; i++) {
    if (indic_glob[i] == 0)
      maj_ind_ent = CS_TRUE;
  }

  if (maj_ind_ent == CS_TRUE) {
    for (i = 0 ; i < 3 ; i++) {
      if (indic_glob[i] == 0)           /* Logique indic_glob 0 à 2 inversée */
        maillage_post->ind_ent[i] = 1;  /* (c.f. remarque ci-dessus) */
      else
        maillage_post->ind_ent[i] = 0;
    }
  }

  /* Dimensions locales */

  maillage_post->nbr_fac = nbr_fac;
  maillage_post->nbr_fbr = nbr_fbr;

  /* Lien sur le maillage nouvellement créé */

  maillage_post->maillage_ext = maillage_ext;
  maillage_post->_maillage_ext = maillage_ext;

}


/*----------------------------------------------------------------------------
 * Mise à jour en cas d'alias des critères de modification au cours du temps
 * des maillages en fonction des propriétés des writers qui lui
 * sont associés :
 *
 * La topologie d'un maillage ne pourra pas être modifiée si le critère
 * de modification minimal résultant est trop faible (i.e. si l'un des
 * writers associés ne permet pas la redéfinition de la topologie du maillage).
 *
 * Les coordonnées des sommets et la connectivité ne pourront être libérés
 * de la mémoire si la critère de modification maximal résultant est trop
 * élevé (i.e. si l'un des writers associés permet l'évolution du maillage
 * en temps, et nécessite donc sa réécriture).
 *----------------------------------------------------------------------------*/

static void _cs_post_ind_mod_alias
(
 const int  indmai              /* --> indice du maillage post en cours       */
)
{
  int  i;

  cs_post_maillage_t  *maillage_post = NULL;
  cs_post_maillage_t  *maillage_ref = NULL;


  /* Mise à jour référence */

  maillage_post = cs_glob_post_maillages + indmai;

  if (maillage_post->alias > -1) {

    maillage_ref = cs_glob_post_maillages + maillage_post->alias;

    if (maillage_post->ind_mod_min < maillage_ref->ind_mod_min)
      maillage_ref->ind_mod_min = maillage_post->ind_mod_min;

    if (maillage_post->ind_mod_max < maillage_ref->ind_mod_max)
      maillage_ref->ind_mod_max = maillage_post->ind_mod_max;

  }


  /* Mise  à jour alias) */

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    if (maillage_post->alias > -1) {

      maillage_ref = cs_glob_post_maillages + maillage_post->alias;

      if (maillage_post->ind_mod_min > maillage_ref->ind_mod_min)
        maillage_post->ind_mod_min = maillage_ref->ind_mod_min;

      if (maillage_post->ind_mod_max > maillage_ref->ind_mod_max)
        maillage_post->ind_mod_max = maillage_ref->ind_mod_max;
    }

  }

}


/*----------------------------------------------------------------------------
 * Découpage des polygones ou polyèdres en éléments simples si nécessaire.
 *----------------------------------------------------------------------------*/

static void _cs_post_divise_poly
(
 cs_post_maillage_t      *const maillage_post,  /* --> maillage ext. associé  */
 const cs_post_writer_t  *const writer          /* --> writer associé         */
)
{
  /* Découpage des polygones ou polyèdres en éléments simples */

  if (fvm_writer_needs_tesselation(writer->writer,
                                   maillage_post->maillage_ext,
                                   FVM_CELL_POLY) > 0)
    fvm_nodal_tesselate(maillage_post->_maillage_ext,
                        FVM_CELL_POLY,
                        NULL);

  if (fvm_writer_needs_tesselation(writer->writer,
                                   maillage_post->maillage_ext,
                                   FVM_FACE_POLY) > 0)
    fvm_nodal_tesselate(maillage_post->_maillage_ext,
                        FVM_FACE_POLY,
                        NULL);
}


/*----------------------------------------------------------------------------
 * Assemblage des valeurs d'une variable définie sur une combinaison de
 * faces de bord et de faces internes (sans indirection) en un tableau
 * défini sur un ensemble unique de faces.
 *
 * La variable résultante n'est pas entrelacée.
 *----------------------------------------------------------------------------*/

static void _cs_post_assmb_var_faces
(
 const fvm_nodal_t      *const maillage_ext,     /* --> maillage externe   */
 const cs_int_t                nbr_fac,          /* --> nb faces internes  */
 const cs_int_t                nbr_fbr,          /* --> nb faces bord      */
 const int                     dim_var,          /* --> dim. variable      */
 const fvm_interlace_t         interlace,        /* --> indic. entrelacage */
 const cs_real_t               var_fac[],        /* --> valeurs faces int. */
 const cs_real_t               var_fbr[],        /* --> valeurs faces bord */
       cs_real_t               var_tmp[]         /* <-- valeurs assemblées */
)
{
  cs_int_t  i, j, pas_1, pas_2;

  cs_int_t  nbr_ent = nbr_fac + nbr_fbr;

  assert(maillage_ext != NULL);

  /* La variable est définie sur les faces internes et faces de bord
     du maillage de post traitement, et a été construite à partir
     des listes de faces internes et de bord correspondantes */

  /* Contribution des faces de bord */

  if (interlace == FVM_INTERLACE) {
    pas_1 = dim_var;
    pas_2 = 1;
  }
  else {
    pas_1 = 1;
    pas_2 = nbr_fbr;
  }

  for (i = 0 ; i < nbr_fbr ; i++) {
    for (j = 0 ; j < dim_var ; j++)
      var_tmp[i + j*nbr_ent] = var_fbr[i*pas_1 + j*pas_2];
  }

  /* Contribution des faces internes */

  if (interlace == FVM_INTERLACE) {
    pas_1 = dim_var;
    pas_2 = 1;
  }
  else {
    pas_1 = 1;
    pas_2 = nbr_fac;
  }

  for (i = 0 ; i < nbr_fac ; i++) {
    for (j = 0 ; j < dim_var ; j++)
      var_tmp[i + nbr_fbr + j*nbr_ent] = var_fac[i*pas_1 + j*pas_2];
  }

}


/*----------------------------------------------------------------------------
 * Ecriture d'un maillage de post traitement en fonction des "writers".
 *----------------------------------------------------------------------------*/

static void _cs_post_ecrit_maillage
(
       cs_post_maillage_t  *const maillage_post,
 const cs_int_t                   nt_cur_abs,   /* --> numéro de pas de temps */
 const cs_real_t                  t_cur_abs     /* --> temps physique courant */
)
{
  int  j;
  cs_bool_t  ecrire_maillage;
  fvm_writer_time_dep_t  dep_temps;

  cs_post_writer_t *writer = NULL;


  /* Boucles sur les "writers" */

  for (j = 0 ; j < maillage_post->nbr_writers ; j++) {

    writer = cs_glob_post_writers + maillage_post->ind_writer[j];

    dep_temps = fvm_writer_get_time_dep(writer->writer);

    ecrire_maillage = CS_FALSE;

    if (dep_temps == FVM_WRITER_FIXED_MESH) {
      if (maillage_post->nt_ecr < 0)
        ecrire_maillage = CS_TRUE;
    }
    else {
      if (maillage_post->nt_ecr < nt_cur_abs && writer->actif == 1)
        ecrire_maillage = CS_TRUE;
    }

    if (ecrire_maillage == CS_TRUE) {
      fvm_writer_set_mesh_time(writer->writer, nt_cur_abs, t_cur_abs);
      fvm_writer_export_nodal(writer->writer, maillage_post->maillage_ext);
    }

    if (ecrire_maillage == CS_TRUE && maillage_post->id == -1)
      _cs_post_ecrit_domaine(writer->writer,
                             maillage_post->maillage_ext,
                             nt_cur_abs,
                             t_cur_abs);

  }

  if (ecrire_maillage == CS_TRUE)
    maillage_post->nt_ecr = nt_cur_abs;

  if (   maillage_post->ind_mod_max == FVM_WRITER_FIXED_MESH
      && maillage_post->_maillage_ext != NULL)
    fvm_nodal_reduce(maillage_post->_maillage_ext, 0);
}


/*----------------------------------------------------------------------------
 * Transformation d'un tableau d'indicateurs (marqueurs) en liste ;
 * renvoie la taille effective de la liste.
 *----------------------------------------------------------------------------*/

static cs_int_t _cs_post_ind_vers_liste
(
 cs_int_t  nbr,         /* <-> taille indicateur                */
 cs_int_t  liste[]      /* <-> indicateur, puis liste           */
)
{
  cs_int_t   cpt, ind;

  for (cpt = 0, ind = 0 ; ind < nbr ; ind++) {
    if (liste[ind] != 0) {
      liste[ind] = 0;
      liste[cpt++] = ind + 1;
    }
  }

  return cpt;
}


/*----------------------------------------------------------------------------
 * Boucle sur les maillages de post traitement pour écriture  des variables
 *----------------------------------------------------------------------------*/

static void _cs_post_ecrit_deplacements
(
 const cs_int_t   nt_cur_abs,         /* --> numéro de pas de temps courant   */
 const cs_real_t  t_cur_abs           /* --> valeur du temps physique associé */
)
{
  int i, j;
  cs_int_t  k, nbr_val;
  fvm_datatype_t datatype;

  fvm_lnum_t   dec_num_parent[1]  = {0};
  cs_post_maillage_t  *maillage_post = NULL;
  cs_post_writer_t  *writer = NULL;
  cs_real_t   *deplacements = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_real_t   *var_ptr[1] = {NULL};


  /* Boucle sur les writers pour vérifier si l'on a quelque chose à faire */
  /*----------------------------------------------------------------------*/

  if (cs_glob_post_deformable == CS_FALSE)
    return;

  for (j = 0 ; j < cs_glob_post_nbr_writers ; j++) {
    writer = cs_glob_post_writers + j;
    if (writer->actif == 1 && writer->ecr_depl == CS_TRUE)
      break;
  }
  if (j == cs_glob_post_nbr_writers)
    return;


  /* Calcul du champ de déformation principal */
  /*------------------------------------------*/

  nbr_val = mesh->n_vertices * 3;

  BFT_MALLOC(deplacements, nbr_val, cs_real_t);

  assert(mesh->n_vertices == 0 || cs_glob_post_coo_som_ini != NULL);

  for (k = 0; k < nbr_val; k++)
    deplacements[k] = mesh->vtx_coord[k] - cs_glob_post_coo_som_ini[k];


  /* Préparation du post traitement */
  /*--------------------------------*/

  if (sizeof(cs_real_t) == sizeof(double))
    datatype = FVM_DOUBLE;
  else if (sizeof(cs_real_t) == sizeof(float))
    datatype = FVM_FLOAT;

  var_ptr[0] = deplacements;


  /* Boucle sur les maillages pour l'écriture des déplacements */
  /*-----------------------------------------------------------*/

  for (i = 0 ; i < cs_glob_post_nbr_maillages ; i++) {

    maillage_post = cs_glob_post_maillages + i;

    for (j = 0 ; j < maillage_post->nbr_writers ; j++) {

      writer = cs_glob_post_writers + maillage_post->ind_writer[j];

      if (writer->actif == 1 && writer->ecr_depl == CS_TRUE) {

        fvm_writer_export_field(writer->writer,
                                maillage_post->maillage_ext,
                                _("displacement"),
                                FVM_WRITER_PER_NODE,
                                3,
                                FVM_INTERLACE,
                                1,
                                dec_num_parent,
                                datatype,
                                (int)nt_cur_abs,
                                (double)t_cur_abs,
                                (const void * *)var_ptr);

      }

    }

  }

  /* Libération mémoire */

  BFT_FREE(deplacements);
}


/*----------------------------------------------------------------------------
 * Écriture du domaine sur maillage de post traitement
 *---------------------------------------------------------------------------*/

static void _cs_post_ecrit_domaine
(
       fvm_writer_t     *const writer,         /* --> writer FVM              */
 const fvm_nodal_t      *const maillage_ext,   /* --> maillage externe        */
       cs_int_t                nt_cur_abs,     /* --> numéro pas de temps     */
       cs_real_t               t_cur_abs       /* --> temps physique associé  */
)
{
  int  dim_ent;
  fvm_lnum_t  i, nbr_ent;
  fvm_datatype_t  datatype;

  fvm_lnum_t   dec_num_parent[1]  = {0};
  cs_int_t *domaine = NULL;

  int _nt_cur_abs = -1;
  double _t_cur_abs = 0.;

  const cs_int_t   *var_ptr[1] = {NULL};

  if (cs_glob_base_nbr < 2 || cs_glob_post_domaine == CS_FALSE)
    return;

  dim_ent = fvm_nodal_get_max_entity_dim(maillage_ext);
  nbr_ent = fvm_nodal_get_n_entities(maillage_ext, dim_ent);


  /* Préparation du numéro de domaine */

  BFT_MALLOC(domaine, nbr_ent, cs_int_t);

  for (i = 0; i < nbr_ent; i++)
    domaine[i] = cs_glob_mesh->domain_num;

  /* Préparation du post traitement */

  if (sizeof(cs_int_t) == 4)
    datatype = FVM_INT32;
  else if (sizeof(cs_real_t) == 8)
    datatype = FVM_INT64;

  var_ptr[0] = domaine;

  if (fvm_writer_get_time_dep(writer) != FVM_WRITER_FIXED_MESH) {
    _nt_cur_abs = nt_cur_abs;
    _t_cur_abs = t_cur_abs;
  }

  fvm_writer_export_field(writer,
                          maillage_ext,
                          _("parallel domain"),
                          FVM_WRITER_PER_ELEMENT,
                          1,
                          FVM_INTERLACE,
                          1,
                          dec_num_parent,
                          datatype,
                          _nt_cur_abs,
                          _t_cur_abs,
                          (const void * *)var_ptr);

  /* Libération mémoire */

  BFT_FREE(domaine);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
