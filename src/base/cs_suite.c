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
 *  Gestion des fichiers suite
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie BFT
 *----------------------------------------------------------------------------*/

#include <bft_file.h>
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 *  Fichiers  `include' associés au fichier courant
 *----------------------------------------------------------------------------*/

#include "cs_suite.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 *  Structures locales
 *============================================================================*/

typedef struct _cs_suite_rec_t {

  char            *nom;          /* Nom de l'enregistrement */
  cs_int_t         ind_support;  /* Indice du support (code cs_suite_support_t
                                    ou utilisateur dans version ultérieure) */
  cs_int_t         nbr_val_ent;  /* Nombre de valeurs/entité support */
  cs_type_t        typ_val;      /* Type de variable */

  cs_int_t         ind_fic;      /* Indice du fichier correspondant */
  bft_file_off_t   pos_fic;      /* Position du début dans le fichier */

} cs_suite_rec_t;

struct _cs_suite_t {

  char            *nom;          /* Nom du fichier suite */

  cs_int_t         nbr_cel;      /* Nombre de cellules du support en lecture */
  cs_int_t         nbr_fac;      /* Nb. faces internes du support en lecture */
  cs_int_t         nbr_fbr;      /* Nb. faces de bord du support en lecture */
  cs_int_t         nbr_som;      /* Nombre de sommets du support en lecture */

  cs_int_t         nbr_rec;      /* Nombre d'enregistrements courant */
  cs_int_t         nbr_rec_max;  /* Nombre d'enregistrements maximal avant
                                    redimensionnement de tab_rec */
  cs_suite_rec_t  *tab_rec;      /* Index des enregistrements */

  cs_int_t         nbr_fic;      /* Nombre de fichiers associés (un seul en
                                    écriture, plusieurs possibles en lecture) */
  bft_file_t     **fic;          /* Pointeurs sur fichiers associés ; tableau
                                    de taille 1 en écriture (un fichier ouvert
                                    réutilise la position 0), et de taille
                                    nbr_fic en lecture (plusieurs fichiers
                                    peuvent être ouverts si données trop
                                    volumineuses pour un seul */

  cs_suite_type_t  type;         /* Type de codage des données */
  cs_suite_mode_t  mode;         /* Lecture ou écriture */

};


/*============================================================================
 *  Constantes et Macros
 *============================================================================*/

/* Taille maximale d'un fichier : 0.95 * entier 32 bits max (2147483647) */

#define CS_SUITE_TAILLE_FIC_MAX   2040109464

#define CS_SUITE_LNG_NOM_TYPE_ELT      2    /* Longueur du nom de type     */
#define CS_SUITE_LNG_BUF_ASCII       120    /* Longueur tampon ligne ASCII */

/*
 * Label MPI (pour distinguer messages associés aux suites des autres,
 * sans passer par un communicateur distinct)
 */

#define CS_SUITE_MPI_TAG     'C'+'S'+'_'+'S'+'U'+'I'+'T'+'E'

/* Formats pour mode texte (voire chaînes globales plus bas) */

#define CS_SUITE_FMT_ASCII_DIM_COL_INT                11 /* %10d    + 1 col */
#define CS_SUITE_FMT_ASCII_DIM_COL_REAL               24 /* %23.15E + 1 col */
#define CS_SUITE_FMT_ASCII_NBR_COL_INT                 7
#define CS_SUITE_FMT_ASCII_NBR_COL_REAL                3

/* API Fortran */
/* ----------- */

/*
 * (longueur max 'usuelle' de nom ; un nom plus long est possible
 * mais provoquera une allocation de mémoire dynamique).
 */

#define CS_SUITE_LNG_NOM                              64


/*============================================================================
 *  Variables globales statiques
 *============================================================================*/

/* Taille minimale d'un buffer sur le rang 0 (pour limiter le
   nombre de blocs lorsqu'on a un grand nombre de processeurs) */

static int cs_suite_taille_buf_def = 1024*1024*8;

/* Chaîne indiquant que la suite continue sur un autre fichier */

static char const cs_suite_nom_fin_fic[]    = "reprise : fin";
static char const cs_suite_nom_decoup_fic[] = "reprise : fic suivant";
static char const cs_suite_nom_partie_fic[] = "reprise : partie num ";

/* Type de valeur (dans l'ordre de définition de cs_type_t dans cs_base.h,
   seuls les types "principaux étant" traités ici) */

static char const cs_suite_nom_typ_elt[3][4]  = {"c  ",
                                                 "i32",
                                                 "r64"};

/* Nom de support (dans l'ordre de définition de cs_suite_support_t) */

static char const *cs_suite_nom_support[5] = {"---",
                                              "cel",
                                              "fac",
                                              "fbr",
                                              "som"};

/* Formats d'écriture de valeurs */

static char  cs_suite_fmt_ascii_tab_int[]  = " %10d";
static char  cs_suite_fmt_ascii_tab_real[] = " %23.15E";

/* Tableau pour API Fortran */

static cs_int_t     cs_glob_suite_ptr_nbr = 0;
static cs_suite_t **cs_glob_suite_ptr_tab = NULL;

#if defined(_CS_HAVE_MPI)

/* Mode de rapatriement sûr (légèrement plus lent) */

/* IBM Blue Gene ou Cray XT */
#if   defined(__blrts__) || defined(__bgp__) \
   || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
static cs_bool_t cs_glob_suite_sync_gather = CS_TRUE;
#else
static cs_bool_t cs_glob_suite_sync_gather = CS_FALSE;
#endif

#endif /* defined(_CS_HAVE_MPI) */


/*============================================================================
 *  Prototypes de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Calcul du nombre de valeurs d'un enregistrement
 *----------------------------------------------------------------------------*/

static cs_int_t cs_loc_suite_calc_nbr_ent
(
 const cs_suite_t          *const suite,       /* --> Suite associée          */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent  /* --> Nb. val/point support   */
);


/*----------------------------------------------------------------------------
 *  Calcul du déplacement dans un fichier correspondant aux valeurs
 *  d'une rubrique.
 *----------------------------------------------------------------------------*/

static bft_file_off_t cs_loc_suite_calc_avance
(
 const cs_suite_t          *const suite,       /* --> Suite associée          */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent, /* --> Nb. val/point support   */
 const cs_type_t                  typ_val      /* --> Type de valeurs         */
);


/*----------------------------------------------------------------------------
 *  Fonction qui initialise une structure fichier àssociée à une suite ;
 *  Le nom de fichier associé correspond au nom de la suite pour le
 *  premier fichier, prolongé éventuellement par '_pxx' pour les fichiers
 *  successifs en cas de découpage du fichier en plusieurs morceaux pour
 *  des gros volumes de données.
 *  On indique en retour si les données se poursuivent sur un autre fichier.
 *----------------------------------------------------------------------------*/

static cs_bool_t cs_loc_suite_fic_ajoute
(
 cs_suite_t       *const suite              /* --> Structure suite            */
);


/*----------------------------------------------------------------------------
 *  Fonction qui lit une liste de valeurs sur un fichier ; on suppose que
 *  l'on est déjà à la bonne position de départ dans le fichier.
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_lit_val
(
 const cs_suite_type_t         type,           /* --> Structure suite         */
 const bft_file_t       *const fic,            /* --> Fichier associé         */
 const cs_int_t                nbr_val,        /* --> Nb. valeurs à lire      */
 const cs_type_t               typ_val,        /* --> Type de valeurs         */
       void             *const val,            /* <-- Valeurs à lire          */
       char                    asc_buffer[]    /* <-> Buffer texte            */
);


/*----------------------------------------------------------------------------
 *  Fonction qui écrit une liste de valeurs sur un fichier
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_ecr_val
(
 const cs_suite_t          *const suite,       /* --> Ptr. fichier suite      */
 const cs_int_t                   nbr_val,     /* --> Nb. valeurs à écrire    */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
 const void                *const val,         /* --> Valeurs à écrire        */
       cs_int_t            *const ind_col      /* <-> Colonne départ/fin      */
);


/*----------------------------------------------------------------------------
 *  Fonction qui termine l'écriture d'une liste de valeurs sur un fichier
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_ecr_val_fin
(
 const cs_suite_t          *const suite,       /* --> Ptr. fichier suite      */
       cs_int_t            *const ind_col      /* <-> Colonne départ/fin      */
);


/*----------------------------------------------------------------------------
 *  Fonction qui analyse le contenu d'un fichier suite en mode texte
 *  On démare en haut de deuxième ligne, ayant lu l'en-tête, et on
 *  indique en retour si la suite continue sur un autre fichier
 *----------------------------------------------------------------------------*/

static cs_bool_t cs_loc_suite_analyse_txt
(
 cs_suite_t  *const  suite                     /* --> Structure suite         */
);


/*----------------------------------------------------------------------------
 *  Fonction qui analyse le contenu d'un fichier suite en mode binaire
 *  On démarre à la position suivantl'en-tête, et on indique en retour
 *  si la suite continue sur un autre fichier
 *----------------------------------------------------------------------------*/

static cs_bool_t cs_loc_suite_analyse_bin
(
 cs_suite_t  *const  suite                     /* --> Structure suite         */
);


/*----------------------------------------------------------------------------
 *  Fonction qui prépare l'index généré lors de l'analyse du fichier à
 *  l'utilisation par les fonctions de lecture d'enregistrements
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_prepare_index
(
 cs_suite_t  *const  suite                /* --> Structure suite              */
);

/*----------------------------------------------------------------------------
 *  Conversion d'arguments de lecture/écriture de l'API Fortran vers l'API C
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_rub_f77_vers_C
(
 const cs_int_t             *const numsui,  /* --> Numéro du fichier suite    */
 const cs_int_t             *const itysup,  /* --> Code type de support       */
 const cs_int_t             *const irtype,  /* --> Entiers ou réels ?         */
       cs_suite_t          **const suite,   /* <-- Pointeur structure suite   */
       cs_suite_support_t   *const support, /* <-- Type de support            */
       cs_type_t            *const typ_val, /* <-- Entiers ou réels           */
       cs_int_t             *const ierror   /* <-- 0 = succès, < 0 = erreur   */
);


#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Fonction qui créé des listes locales->globales associées à chaque bloc.
 *
 *  Pour le processus "maître E/S", le tableau "nbr_ent_bloc" indiquant le
 *  nombre d'entités appartenant à chaque bloc, et ceci pour chaque domaine :
 *  nbr_ent(bloc, domaine) = nbr_ent_bloc[nbr_bloc * ind_dom + ind_bloc] ;
 *  Pour les autres processus, le tableau "nbr_ent_loc" indique le nombre
 *  d'entités locales appartenant à chaque bloc.
 *  Les tableaux "lst_ent_loc" et "lst_ent_glob" donnent les indices locaux
 *  et globaux des entités locales appartenant à chaque bloc (nbr_ent_loc[0]
 *  premières valeurs pour le bloc 0, nbr_ent_loc[1] suivantes pour le
 *  second, etc.).
 *
 *  Ces tableaux sont alloués ici et devront être libérés par la suite.
 *
 *  On renvoie la plus grande valeur (globale) de nbr_ent_loc.
 *----------------------------------------------------------------------------*/

static cs_int_t cs_loc_suite_cree_listes_ent
(
 const cs_int_t            nbr_bloc,      /* --> Nombre de blocs              */
 const cs_int_t            nbr_ent_glob,  /* --> Nombre global d'entités      */
 const cs_int_t            nbr_ent_loc,   /* --> Nombre local d'entités       */
 const fvm_gnum_t   *const num_glob_ent,  /* --> Numéros globaux des entités  */
       cs_int_t     *const pas_bloc,      /* <-- Taille blocs (sauf dernier)  */
       cs_int_t   * *const nbr_ent_bloc,  /* <-- Nombre d'entités par bloc    */
       cs_int_t   * *const lst_ent_loc,   /* <-- Liste des entités par bloc   */
       cs_int_t   * *const lst_ent_glob   /* <-- Liste des entités par bloc   */
);


/*----------------------------------------------------------------------------
 *  Fonction qui lit les valeurs associées à une variable définie sur
 *  une entité de maillage.
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_lit_val_ent
(
 const cs_suite_t   *const suite,         /* --> Structure suite              */
 const bft_file_t   *const fic,           /* --> Fichier associé aux valeurs  */
 const cs_int_t            nbr_bloc,      /* --> Nombre de blocs              */
 const cs_int_t            nbr_ent_glob,  /* --> Nombre global d'entités      */
 const cs_int_t            nbr_ent_loc,   /* --> Nombre local d'entités       */
 const fvm_gnum_t   *const num_glob_ent,  /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
       cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
);


/*----------------------------------------------------------------------------
 *  Fonction qui écrit les valeurs associées à une variable définie sur
 *  une entité de maillage.
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_ecr_val_ent
(
 const cs_suite_t  *const  suite,         /* --> Structure suite              */
 const cs_int_t            nbr_bloc,      /* --> Nombre de blocs              */
 const cs_int_t            nbr_ent_glob,  /* --> Nombre global d'entités      */
 const cs_int_t            nbr_ent_loc,   /* --> Nombre local d'entités       */
 const fvm_gnum_t   *const num_glob_ent,  /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
 const cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
);


/*----------------------------------------------------------------------------
 *  Fonction qui distribue une liste de variables
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_distr_val
(
 const cs_int_t                   nbr_val,     /* --> Nb. valeurs             */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
       void                *const val          /* <-> Valeurs à lire          */
) ;

#endif /* #if defined(_CS_HAVE_MPI) */


/*----------------------------------------------------------------------------
 *  Permutation des valeurs d'un tableau renuméroté en lecture
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_permute_lec
(
 const cs_int_t            nbr_ent,       /* --> Nombre d'entités             */
 const cs_int_t     *const num_ent_ini,   /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
       cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
);


/*----------------------------------------------------------------------------
 *  Permutation des valeurs d'un tableau renuméroté en écriture
 *----------------------------------------------------------------------------*/

static cs_byte_t * cs_loc_suite_permute_ecr
(
 const cs_int_t            nbr_ent,       /* --> Nombre d'entités             */
 const cs_int_t     *const num_ent_ini,   /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
 const cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
);


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
 * INTEGER          IREAWR      : --> : 1 pour lecture, 2 pour écriture
 * INTEGER          IFORMA      : --> : 0 pour binaire, 1 pour formaté
 * INTEGER          NUMSUI      : <-- : Numéro du fichier suite ouvert
 * INTEGER          IERROR      : <-- : 0 pour succès, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (opnsui, OPNSUI)
(
 const char       *const nomsui,  /* --> Nom du fichier                       */
 const cs_int_t   *const lngnom,  /* --> Longueur du nom                      */
 const cs_int_t   *const ireawr,  /* --> 1 pour lecture, 2 pour écriture      */
 const cs_int_t   *const iforma,  /* --> 0 pour binaire, 1 pour formaté      */
       cs_int_t   *const numsui,  /* <-- Numéro du ficher suite ouvert        */
       cs_int_t   *const ierror   /* <-- 0 pour succès, < 0 pour erreur       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' éventuels F77, */
                                  /*     inutilisés lors de l'appel mais      */
                                  /*     placés par de nombreux compilateurs) */
)
{
  char    *nombuf;

  cs_int_t ind;

  cs_suite_mode_t suite_mode;
  cs_suite_type_t suite_type;


  /* Initialisation */

  *ierror = CS_SUITE_SUCCES;

  /* Traitement du nom pour l'API C */

  nombuf = cs_base_chaine_f_vers_c_cree(nomsui, *lngnom);

  /* Options de création du fichier */

  {
    switch(*ireawr) {
    case 1:
      suite_mode = CS_SUITE_MODE_LECTURE;
      break;
    case 2:
      suite_mode = CS_SUITE_MODE_ECRITURE;
      break;
    default:
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The access mode of the restart file <%s>\n"
                   "must be equal to 1 (read) or 2 (write) and not <%d>."),
                 nombuf, (int)(*ireawr));

      *ierror = CS_SUITE_ERR_MODE;
    }

  }

  if (*ierror == CS_SUITE_SUCCES) {

    if (suite_mode == CS_SUITE_MODE_ECRITURE) {

      switch (*iforma) {
      case 0:
        suite_type = CS_SUITE_TYPE_BINAIRE;
        break;
      case 1:
        suite_type = CS_SUITE_TYPE_ASCII;
        break;
      default:
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(_("The type of the restart file <%s>\n"
                     "must be equal to 0 (binary) or 1 (formatted) and not <%d>\n"
                     "(default is binary)."),
                   nombuf, (int)(*iforma));
        *ierror = CS_SUITE_ERR_TYPE_FIC;
      }

    }
    else
      suite_type = CS_SUITE_TYPE_BINAIRE; /* Inutile, mais évite warning
                                             sous Purify */

  }

  /* Recherche d'un emplacement disponible */

  if (*ierror == CS_SUITE_SUCCES) {

    for (ind = 0 ;
         ind < cs_glob_suite_ptr_nbr && cs_glob_suite_ptr_tab[ind] != NULL ;
         ind++);

    /* Si aucun emplacement disponible, on autorise plus de fichiers suite */

    if (ind == cs_glob_suite_ptr_nbr) {

      BFT_REALLOC(cs_glob_suite_ptr_tab, cs_glob_suite_ptr_nbr * 2,
                  cs_suite_t *);
      for (ind = cs_glob_suite_ptr_nbr ;
           ind < cs_glob_suite_ptr_nbr * 2 ;
           ind++)
        cs_glob_suite_ptr_tab[ind] = NULL;
      cs_glob_suite_ptr_nbr *= 2;

    }

  }

  /* Création du fichier suite */

  if (*ierror == CS_SUITE_SUCCES)
    cs_glob_suite_ptr_tab[ind] = cs_suite_cree(nombuf,
                                               suite_mode,
                                               suite_type);

  /* Libération de mémoire si nécessaire */

  nombuf = cs_base_chaine_f_vers_c_detruit(nombuf);

  /*
   * On renvoie le numéro de position du pointeur dans le tableau
   * (indice + 1 pour avoir un numéro de 1 à n, plus classique en F77)
  */

  if (*ierror == CS_SUITE_SUCCES)
    *numsui = ind + 1;
  else
    *numsui = -1;

}


/*----------------------------------------------------------------------------
 * Fermeture d'un fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE CLSSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : numéro du fichier suite à fermer
 * INTEGER          IERROR      : <-- : 0 pour succès, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (clssui, CLSSUI)
(
 const cs_int_t   *const numsui,  /* --> Numéro du ficher suite à fermer      */
       cs_int_t   *const ierror   /* <-- Numéro du ficher suite ouvert        */
)
{
  cs_int_t indsui = *numsui - 1;

  *ierror = CS_SUITE_SUCCES;

  /* Vérification que le fichier est valide */

  if (   indsui < 0
      || indsui > cs_glob_suite_ptr_nbr
      || cs_glob_suite_ptr_tab[indsui] == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The restart file number <%d> cannot be closed\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *ierror = CS_SUITE_ERR_NUM_FIC;
    return;
  }

  /* Fermeture du fichier */

  cs_suite_detruit(cs_glob_suite_ptr_tab[indsui]);
  cs_glob_suite_ptr_tab[indsui] = NULL;

}


/*----------------------------------------------------------------------------
 *  Vérification du support associé à un fichier suite ;
 *  On renvoie pour chaque type d'entité 1 si le nombre d'entités associées
 *  au fichier suite correspond au nombre d'entités en cours (et donc que
 *  l'on considère que le support est bien le même), 0 sinon.
 *
 * Interface Fortran :
 *
 * SUBROUTINE TSTSUI (NUMSUI, INDCEL, INDFAC, INDFBR, INDSOM)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numéro du fichier suite
 * INTEGER          INDCEL      : <-- : Indicateur corresp. cellules
 * INTEGER          INDFAC      : <-- : Indicateur corresp. faces internes
 * INTEGER          INDFBR      : <-- : Indicateur corresp. faces de bord
 * INTEGER          INDSOM      : <-- : Indicateur corresp. sommets
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstsui, TSTSUI)
(
 const cs_int_t  *const numsui,   /* --> Numéro du fichier suite              */
       cs_int_t  *const indcel,   /* <-- Indicateur corresp. cellules         */
       cs_int_t  *const indfac,   /* <-- Indicateur corresp. faces internes   */
       cs_int_t  *const indfbr,   /* <-- Indicateur corresp. faces de bord    */
       cs_int_t  *const indsom    /* <-- Indicateur corresp. sommets          */
)
{
  cs_bool_t  corresp_cel, corresp_fac, corresp_fbr, corresp_som;

  cs_int_t   indsui   = *numsui - 1;

  /* Pointeur de structure suite associé */

  if (   indsui < 0
      || indsui > cs_glob_suite_ptr_nbr
      || cs_glob_suite_ptr_tab[indsui] == NULL) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Information on the restart file number <%d> unavailable\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *indcel = 0;
    *indfac = 0;
    *indfbr = 0;
    *indsom = 0;
    return;
  }

  else {

    cs_suite_verif_support(cs_glob_suite_ptr_tab[indsui],
                           &corresp_cel, &corresp_fac,
                           &corresp_fbr, &corresp_som);

    *indcel = (corresp_cel == CS_TRUE ? 1 : 0);
    *indfac = (corresp_fac == CS_TRUE ? 1 : 0);
    *indfbr = (corresp_fbr == CS_TRUE ? 1 : 0);
    *indsom = (corresp_som == CS_TRUE ? 1 : 0);

  }

}


/*----------------------------------------------------------------------------
 *  Affichage de l'index associé à un fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE INFSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numéro du fichier suite
 *----------------------------------------------------------------------------*/

void CS_PROCF (infsui, INFSUI)
(
 const cs_int_t  *const numsui    /* --> Numéro du fichier suite              */
)
{
  cs_int_t   indsui   = *numsui - 1;

  /* Pointeur de structure suite associé */

  if (   indsui < 0
      || indsui > cs_glob_suite_ptr_nbr
      || cs_glob_suite_ptr_tab[indsui] == NULL) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Information on the restart file number <%d> unavailable\n"
                 "(file already closed or invalid number)."), (int)(*numsui));
  }
  else {

    cs_suite_affiche_index(cs_glob_suite_ptr_tab[indsui]);

  }
}


/*----------------------------------------------------------------------------
 * Lecture d'une rubrique sur fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE LECSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numéro du fichier suite
 * CHARACTER*       NOMRUB      : --> : Nom de la rubrique
 * INTEGER          LNGNOM      : --> : Longueur du nom de la rubrique
 * INTEGER          ITYSUP      : --> : Type de support :
 *                              :     :  0 : scalaire (pas de support)
 *                              :     :  1 : cellules
 *                              :     :  2 : faces internes
 *                              :     :  3 : faces de bord
 *                              :     :  4 : sommets (si disponibles)
 * INTEGER          NBVENT      : --> : Nb. valeurs par entité de support
 * INTEGER          IRTYPE      : --> : 1 pour entiers, 2 pour double précision
 * (?)              TABVAR      : <-> : Tableau des valeurs à lire
 * INTEGER          IERROR      : <-- : 0 pour succès, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (lecsui, LECSUI)
(
 const cs_int_t   *const numsui,  /* --> Numéro du fichier suite              */
 const char       *const nomrub,  /* --> Nom de la rubrique                   */
 const cs_int_t   *const lngnom,  /* --> Longueur du nom de la rubrique       */
 const cs_int_t   *const itysup,  /* --> Type de support (voir ci-dessus)     */
 const cs_int_t   *const nbvent,  /* --> Nb. valeurs par entité du support    */
 const cs_int_t   *const irtype,  /* --> 1 pour entiers, 2 pour double préc.  */
       void       *const tabvar,  /* <-- Tableur des valeurs à lire           */
       cs_int_t   *const ierror   /* <-- 0 pour succès, < 0 pour erreur       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' éventuels F77, */
                                  /*     inutilisés lors de l'appel mais      */
                                  /*     placés par de nombreux compilateurs) */
)
{
  char    *nombuf;

  cs_type_t          typ_val;

  cs_suite_t         *suite;
  cs_suite_support_t  support;


  *ierror = CS_SUITE_SUCCES;

  /* Traitement du nom pour l'API C */

  nombuf = cs_base_chaine_f_vers_c_cree(nomrub,
                                        *lngnom);

  /* Traitement des autres arguments pour l'API C */

  cs_loc_suite_rub_f77_vers_C(numsui,
                              itysup,
                              irtype,
                              &suite,
                              &support,
                              &typ_val,
                              ierror);

  if (*ierror < CS_SUITE_SUCCES)
    return;

  /* Lecture de la rubrique */

  *ierror = cs_suite_lit_rub(suite,
                             nombuf,
                             support,
                             *nbvent,
                             typ_val,
                             tabvar);

  /* Libération de mémoire si nécessaire */

  nombuf = cs_base_chaine_f_vers_c_detruit(nombuf);

}


/*----------------------------------------------------------------------------
 * Écriture d'une rubrique sur fichier suite
 *
 * Interface Fortran :
 *
 * SUBROUTINE ECRSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Numéro du fichier suite
 * CHARACTER*       NOMRUB      : --> : Nom de la rubrique
 * INTEGER          LNGNOM      : --> : Longueur du nom de la rubrique
 * INTEGER          ITYSUP      : --> : Type de support :
 *                              :     :  0 : scalaire (pas de support)
 *                              :     :  1 : cellules
 *                              :     :  2 : faces internes
 *                              :     :  3 : faces de bord
 *                              :     :  4 : sommets (si disponibles)
 * INTEGER          NBVENT      : --> : Nb. valeurs par entité de support
 * INTEGER          IRTYPE      : --> : 1 pour entiers, 2 pour double précision
 * (?)              TABVAR      : --> : Tableau des valeurs fournies
 * INTEGER          IERROR      : <-- : 0 pour succès, < 0 pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrsui, ECRSUI)
(
 const cs_int_t   *const numsui,  /* --> Numéro du fichier suite              */
 const char       *const nomrub,  /* --> Nom de la rubrique                   */
 const cs_int_t   *const lngnom,  /* --> Longueur du nom de la rubrique       */
 const cs_int_t   *const itysup,  /* --> Type de support (voir ci-dessus)     */
 const cs_int_t   *const nbvent,  /* --> Nb. valeurs par entité du support    */
 const cs_int_t   *const irtype,  /* --> 1 pour entiers, 2 pour double préc.  */
 const void       *const tabvar,  /* --> Tableur des valeurs fournies         */
       cs_int_t   *const ierror   /* <-- 0 pour succès, < 0 pour erreur       */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' éventuels F77, */
                                  /*     inutilisés lors de l'appel mais      */
                                  /*     placés par de nombreux compilateurs) */
)
{
  char    *nombuf;

  cs_type_t          typ_val;

  cs_suite_t         *suite;
  cs_suite_support_t  support;


  *ierror = CS_SUITE_SUCCES;

  /* Traitement du nom pour l'API C */

  nombuf = cs_base_chaine_f_vers_c_cree(nomrub,
                                        *lngnom);

  /* Traitement des autres arguments pour l'API C */

  cs_loc_suite_rub_f77_vers_C(numsui,
                              itysup,
                              irtype,
                              &suite,
                              &support,
                              &typ_val,
                              ierror);

  if (*ierror < CS_SUITE_SUCCES)
    return;

  /* Écriture de la rubrique */

  cs_suite_ecr_rub(suite,
                   nombuf,
                   support,
                   *nbvent,
                   typ_val,
                   tabvar);

  /* Libération de mémoire si nécessaire */

  nombuf = cs_base_chaine_f_vers_c_detruit(nombuf);

}


/*============================================================================
 *  Définitions de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui initialise un fichier suite
 *----------------------------------------------------------------------------*/

cs_suite_t * cs_suite_cree
(
 const char             *const nom,         /* --> nom de base du fichier     */
 const cs_suite_mode_t         mode,        /* --> Lecture ou écriture        */
 const cs_suite_type_t         type         /* --> ASCII ou binaire           */
)
{
  cs_bool_t     lit_fichier;
  cs_suite_t  * suite;

  cs_mesh_t  *mesh = cs_glob_mesh;

  /* Instructions */

  BFT_MALLOC(suite, 1, cs_suite_t);

  /* Construction du nom du fichier suite */

  BFT_MALLOC(suite->nom, strlen(nom) + 1, char);

  strcpy(suite->nom, nom);


  /* Initialisation des autres champs */

  suite->mode = mode;
  suite->type
    = (mode == CS_SUITE_MODE_LECTURE) ? CS_SUITE_TYPE_BINAIRE : type ;

  suite->nbr_fic = 0;
  suite->fic     = NULL;

  /* Initialisation des données du support */

  if (suite->mode == CS_SUITE_MODE_LECTURE) {
    suite->nbr_cel = 0;
    suite->nbr_fac = 0;
    suite->nbr_fbr = 0;
    suite->nbr_som = 0;
  }
  else {
    suite->nbr_cel = mesh->n_g_cells;
    suite->nbr_fac = mesh->n_g_i_faces;
    suite->nbr_fbr = mesh->n_g_b_faces;
    suite->nbr_som = mesh->n_g_vertices;
  }

  /* Initialisation de l'index (pour lecture) */

  suite->nbr_rec = 0;

  suite->tab_rec = NULL;
  suite->nbr_rec_max = 0;

  if (suite->mode == CS_SUITE_MODE_LECTURE && cs_glob_base_rang <= 0) {
    suite->nbr_rec_max = 1;
    BFT_MALLOC(suite->tab_rec, suite->nbr_rec_max, cs_suite_rec_t);
  }

  /*
    Ouverture du ou des fichier(s) associé(s) et analyse du contenu
    en lecture ou écriture de l'en-tête en écriture.
  */

  if (cs_glob_base_rang <= 0) {

    lit_fichier = CS_TRUE;

    while (lit_fichier == CS_TRUE)
      lit_fichier = cs_loc_suite_fic_ajoute(suite);

  }

  /* En lecture, nettoyage et distribution de l'index */

  if (suite->mode == CS_SUITE_MODE_LECTURE)
    cs_loc_suite_prepare_index(suite);

  return suite;

}


/*----------------------------------------------------------------------------
 *  Fonction qui détruit la structure associée à un fichier suite (et ferme
 *  le fichier associé) ; elle renvoie un pointeur NULL.
 *----------------------------------------------------------------------------*/

cs_suite_t * cs_suite_detruit
(
 cs_suite_t * suite                         /* --> Fichier suite              */
)
{

  assert(suite != NULL);

  if (suite->fic != NULL) {

    cs_int_t ind, nbr_fic;

    if (suite->mode == CS_SUITE_MODE_ECRITURE) {

      nbr_fic = 1;

      /* Marquage numéro de suite sur fichier suivant */

      switch(suite->type) {
      case CS_SUITE_TYPE_ASCII:
        bft_file_printf(suite->fic[0], "[%s]\n",
                        cs_suite_nom_fin_fic);
        break;
      case CS_SUITE_TYPE_BINAIRE:
        {
          cs_int_t  buf[4] = {0, 0, 0, 0};
          buf[0] = strlen(cs_suite_nom_fin_fic) + 1;
          bft_file_write(buf, sizeof(cs_int_t), 4, suite->fic[0]);
          bft_file_write(cs_suite_nom_fin_fic, 1, buf[0], suite->fic[0]);
        }
        break;
      }
    }
    else /* if (suite->mode == CS_SUITE_MODE_LECTURE) */

      nbr_fic = suite->nbr_fic;

    /* Fermeture des fichiers et libération structures associées */

    for (ind = 0 ; ind < nbr_fic ; ind++)
      bft_file_free(suite->fic[ind]);

    BFT_FREE(suite->fic);

  }

  /* Libérations de l'index */

  if (suite->nbr_rec > 0) {
    cs_int_t ind;
    for (ind = 0 ; ind < suite->nbr_rec ; ind++)
      BFT_FREE((suite->tab_rec[ind]).nom);
  }
  if (suite->tab_rec != NULL)
    BFT_FREE(suite->tab_rec);

  /* Dernières libérations mémoire */

  BFT_FREE(suite->nom);
  BFT_FREE(suite);

  return NULL;

}


/*----------------------------------------------------------------------------
 *  Fonction qui vérifie le support associé à un fichier suite ;
 *  On renvoie pour chaque type d'entité CS_TRUE si le nombre d'entités
 *  associées au fichier suite correspond au nombre d'entités en cours (et
 *  donc que l'on considère que le support est bien le même), CS_FALSE sinon.
 *----------------------------------------------------------------------------*/

void cs_suite_verif_support
(
 const cs_suite_t  *const suite,            /* --> Fichier suite              */
       cs_bool_t   *const corresp_cel,      /* <-- Corresp. cellules          */
       cs_bool_t   *const corresp_fac,      /* <-- Corresp. faces internes    */
       cs_bool_t   *const corresp_fbr,      /* <-- Corresp. faces de bord     */
       cs_bool_t   *const corresp_som       /* <-- Corresp. sommets           */
)
{
  cs_mesh_t  *mesh = cs_glob_mesh;

  assert(suite != NULL);

  *corresp_cel
    = (mesh->n_g_cells == (fvm_gnum_t)suite->nbr_cel ? CS_TRUE : CS_FALSE);
  *corresp_fac
    = (mesh->n_g_i_faces == (fvm_gnum_t)suite->nbr_fac ? CS_TRUE : CS_FALSE);
  *corresp_fbr
    = (mesh->n_g_b_faces == (fvm_gnum_t)suite->nbr_fbr ? CS_TRUE : CS_FALSE);
  *corresp_som
    = (mesh->n_g_vertices == (fvm_gnum_t)suite->nbr_som ? CS_TRUE : CS_FALSE);

  /* Messages d'avertissement pour listing "maître E/S" */

  if (cs_glob_base_rang <= 0) {

    /*
     * Dimensions du support
     */
    if ((fvm_gnum_t)suite->nbr_cel != mesh->n_g_cells) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The number of cells associated with the restart file\n"
                   "<%s> is %d and does not correspond to the current mesh.\n"),
                 suite->nom, (int)suite->nbr_cel);
    }
    if ((fvm_gnum_t)suite->nbr_fac != mesh->n_g_i_faces) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The number of interior faces associated with the restart file\n"
                   "<%s> is %d and does not correspond to the current mesh.\n"),
                 suite->nom, (int)suite->nbr_fac);
    }
    if ((fvm_gnum_t)suite->nbr_fbr != mesh->n_g_b_faces) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The number of boundary faces associated with the restart file\n"
                   "<%s> is %d and does not correspond to the current mesh\n"),
                 suite->nom, (int)suite->nbr_fbr);
    }
    if ((fvm_gnum_t)suite->nbr_som != mesh->n_g_vertices) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The number of vertices associated with the restart file\n"
                   "<%s> is %d and does not correspond to the current mesh\n"),
                 suite->nom, (int)suite->nbr_som);
    }

  }

}


/*----------------------------------------------------------------------------
 *  Fonction qui affiche l'index généré lors de l'analyse du fichier
 *----------------------------------------------------------------------------*/

void cs_suite_affiche_index
(
 const cs_suite_t  *const  suite          /* --> Structure suite              */
)
{
  cs_int_t  ind_rec;

  assert(suite != NULL);

  bft_printf(_("  Index associated to the restart: %s\n"
               "  (location, n_val/ent, val_type, [file_ind, file_pos], name):\n"),
             suite->nom);

  for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++)
    bft_printf("    %s  %d  %s  [%2d %10d]  %s\n",
               cs_suite_nom_support[(suite->tab_rec[ind_rec]).ind_support],
               (suite->tab_rec[ind_rec]).nbr_val_ent,
               cs_suite_nom_typ_elt[(suite->tab_rec[ind_rec]).typ_val],
               (suite->tab_rec[ind_rec]).ind_fic,
               (int)(suite->tab_rec[ind_rec]).pos_fic,
               (suite->tab_rec[ind_rec]).nom);

  bft_printf("\n");

}


/*----------------------------------------------------------------------------
 *  Fonction qui lit un enregistrement sur fichier suite ; On renvoie 0
 *  (CS_SUITE_SUCCES) en cas de succès, une valeur négative (de type
 *  CS_SUITE_ERR_xxx) en cas d'échec.
 *----------------------------------------------------------------------------*/

cs_int_t cs_suite_lit_rub
(
 const cs_suite_t          *const suite,       /* --> Ptr. structure suite    */
 const char                *const nom_rub,     /* --> Nom de la rubrique      */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent, /* --> Nb. val/point support   */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
       void                *const val          /* <-- Valeurs à lire          */
)
{

  char         buffer_ascii[CS_SUITE_LNG_BUF_ASCII + 1];

  cs_int_t     nbr_val_tot, nbr_ent_glob, nbr_ent_loc;
  fvm_gnum_t  *num_glob_ent;

  cs_int_t     ind_rec;
  bft_file_t  *fic;

  cs_mesh_t  *mesh = cs_glob_mesh;

  assert(suite != NULL);

  buffer_ascii[0] = '\0';

  /* Vérification du support associé */

  switch (support) {

  case CS_SUITE_SUPPORT_SCAL:
    break;
  case CS_SUITE_SUPPORT_CEL:
    if (mesh->n_g_cells != (fvm_gnum_t)suite->nbr_cel)
      return CS_SUITE_ERR_SUPPORT;
    break;
  case CS_SUITE_SUPPORT_FAC_INT:
    if (mesh->n_g_i_faces != (fvm_gnum_t)suite->nbr_fac)
      return CS_SUITE_ERR_SUPPORT;
    break;
  case CS_SUITE_SUPPORT_FAC_BRD:
    if (mesh->n_g_b_faces != (fvm_gnum_t)suite->nbr_fbr)
      return CS_SUITE_ERR_SUPPORT;
    break;
  case CS_SUITE_SUPPORT_SOM:
    if (mesh->n_g_vertices != (fvm_gnum_t)suite->nbr_som)
      return CS_SUITE_ERR_SUPPORT;
    break;
  default:
    assert(support <= CS_SUITE_SUPPORT_SOM);

  }

  /* On recherche l'enregistrement correspondant dans l'index */

  for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++) {
    if (strcmp((suite->tab_rec[ind_rec]).nom, nom_rub) == 0)
      break;
  }

  /* Si on n'a pas trouvé l'enregistrement */

  if (ind_rec >= suite->nbr_rec)
    return CS_SUITE_ERR_EXISTE;

  /*
    Si le support ne correspond pas ; on recherche un enregistrement
    de même nom mais de support approprié.
  */

  if (   (cs_suite_support_t)((suite->tab_rec[ind_rec]).ind_support)
      != support) {

    for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++) {
      if (   (strncmp((suite->tab_rec[ind_rec]).nom,
                      nom_rub,
                      strlen((suite->tab_rec[ind_rec]).nom)) == 0)
          && ((cs_suite_support_t)((suite->tab_rec[ind_rec]).ind_support)
              == support))
      break;
    }

    if (ind_rec >= suite->nbr_rec)
      return CS_SUITE_ERR_SUPPORT;

  }

  /* Si le nombre de valeurs par entité ne correspond pas */

  if ((suite->tab_rec[ind_rec]).nbr_val_ent != nbr_val_ent)
     return CS_SUITE_ERR_NBR_VAL;

  /* Si le type de valeurs ne correspond pas */

  if ((cs_type_t)((suite->tab_rec[ind_rec]).typ_val) != typ_val)
    return CS_SUITE_ERR_TYPE_VAL;

   /* Sinon, on se positionne pour effectuer la lecture */

  if (cs_glob_base_rang <= 0) {

    fic = suite->fic[(suite->tab_rec[ind_rec]).ind_fic];

    bft_file_seek(fic, (suite->tab_rec[ind_rec]).pos_fic, BFT_FILE_SEEK_SET);

  }

  /* Contenu de la rubrique */
  /*------------------------*/

  nbr_val_tot = cs_loc_suite_calc_nbr_ent(suite,
                                          support,
                                          nbr_val_ent);

  /* En mode non-parallèle */

  if (cs_glob_base_rang < 0)
    cs_loc_suite_lit_val(suite->type, fic, nbr_val_tot, typ_val, val,
                         buffer_ascii);

#if defined(_CS_HAVE_MPI)

  /* Variable indépendante du support et commune aux processus */

  else if (support == CS_SUITE_SUPPORT_SCAL) {

    if (cs_glob_base_rang == 0)
      cs_loc_suite_lit_val(suite->type, fic, nbr_val_tot, typ_val, val,
                           buffer_ascii);

    cs_loc_suite_distr_val(nbr_val_tot, typ_val, val);

  }

  /* En mode parallèle sur une entité de maillage */

  else {

    cs_int_t  nbr_bloc;

    switch(support) {

    case CS_SUITE_SUPPORT_CEL:
      nbr_ent_glob = mesh->n_g_cells;
      nbr_ent_loc  = mesh->n_cells;
      num_glob_ent = mesh->global_cell_num;
      break;
    case CS_SUITE_SUPPORT_FAC_INT:
      nbr_ent_glob = mesh->n_g_i_faces;
      nbr_ent_loc  = mesh->n_i_faces;
      num_glob_ent = mesh->global_i_face_num;
      break;
    case CS_SUITE_SUPPORT_FAC_BRD:
      nbr_ent_glob = mesh->n_g_b_faces;
      nbr_ent_loc  = mesh->n_b_faces;
      num_glob_ent = mesh->global_b_face_num;
      break;
    case CS_SUITE_SUPPORT_SOM:
      nbr_ent_glob = mesh->n_g_vertices;
      nbr_ent_loc  = mesh->n_vertices;
      num_glob_ent = mesh->global_vtx_num;
      break;
    default:
      assert(   support > CS_SUITE_SUPPORT_SCAL
             && support <= CS_SUITE_SUPPORT_SOM);

    }

    nbr_bloc = (  ((sizeof(cs_real_t) * nbr_ent_glob * nbr_val_ent) - 1)
                / cs_suite_taille_buf_def) + 1;
    if (nbr_bloc > cs_glob_base_nbr)
      nbr_bloc = cs_glob_base_nbr;
    if (nbr_bloc == 0 )
      nbr_bloc = 1;

    cs_loc_suite_lit_val_ent(suite,
                             fic,
                             nbr_bloc,
                             nbr_ent_glob,
                             nbr_ent_loc,
                             num_glob_ent,
                             nbr_val_ent,
                             typ_val,
                             (cs_byte_t *)val);

  }

#endif /* #if defined(_CS_HAVE_MPI) */

  /* Traitement des renumérotations éventuelles */

  switch (support) {

  case CS_SUITE_SUPPORT_FAC_INT:
    cs_loc_suite_permute_lec(mesh->n_i_faces,
                             mesh->init_i_face_num,
                             nbr_val_ent,
                             typ_val,
                             val);
    break;
  case CS_SUITE_SUPPORT_FAC_BRD:
    cs_loc_suite_permute_lec(mesh->n_b_faces,
                             mesh->init_b_face_num,
                             nbr_val_ent,
                             typ_val,
                             val);
    break;
  default:
    break;

  }

  /* Retour */

  return CS_SUITE_SUCCES;

}


/*----------------------------------------------------------------------------
 *  Fonction qui écrit un enregistrement sur fichier suite
 *----------------------------------------------------------------------------*/

void cs_suite_ecr_rub
(
       cs_suite_t          *const suite,       /* --> Ptr. structure suite    */
 const char                *const nom_rub,     /* --> Nom de la rubrique      */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent, /* --> Nb. val/point support   */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
 const void                *const val          /* --> Valeurs à écrire        */
)
{

  cs_int_t         nbr_val_tot, nbr_ent_glob, nbr_ent_loc;
  fvm_gnum_t      *num_glob_ent;
  bft_file_off_t   taille_rub, pos_fic_cur;
  cs_byte_t       *val_tmp;

  cs_int_t         ind_col_0 = 0;

  cs_mesh_t   *mesh = cs_glob_mesh;


  assert(suite != NULL);
  assert(   suite->type == CS_SUITE_TYPE_ASCII
         || suite->type == CS_SUITE_TYPE_BINAIRE);


  /* En fonction de la taille max, on bascule ou non sur un nouveau fichier. */
  /* ----------------------------------------------------------------------- */

  if (cs_glob_base_rang <= 0) {

    /* Taille données */
    taille_rub = cs_loc_suite_calc_avance(suite,
                                          support,
                                          nbr_val_ent,
                                          typ_val);

    /*
     * Taille en-tête (surestimée : on doit pouvoir contenir l'entête,
     * et il doit rester assez de place après l'enregistrement pour
     * écrire un indicateur de poursuite sur fichier suivant, sans
     * dépasser la taille indiquée)
     */

    taille_rub += strlen(nom_rub);
    taille_rub += 200;

    /* Traitement */

    pos_fic_cur = bft_file_tell(suite->fic[0]);

    if (CS_SUITE_TAILLE_FIC_MAX - taille_rub <= pos_fic_cur) {

      switch(suite->type) {
      case CS_SUITE_TYPE_ASCII:
        bft_file_printf(suite->fic[0], "[%s]\n", cs_suite_nom_decoup_fic);
        break;
      case CS_SUITE_TYPE_BINAIRE:
        {
          cs_int_t  buf[4] = {0, 0, 0, 0};
          buf[0] = strlen(cs_suite_nom_decoup_fic) + 1;
          bft_file_write(buf, sizeof(cs_int_t), 4, suite->fic[0]);
          bft_file_write(cs_suite_nom_decoup_fic, 1, buf[0], suite->fic[0]);
        }
        break;
      }

      /* Bascule effective sur fichier suivant */

      cs_loc_suite_fic_ajoute(suite);

      /* Marquage numéro de suite sur fichier suivant */

      switch(suite->type) {
      case CS_SUITE_TYPE_ASCII:
        bft_file_printf(suite->fic[0], "[%s%2d]\n",
                        cs_suite_nom_partie_fic, (int)suite->nbr_fic);
        break;
      case CS_SUITE_TYPE_BINAIRE:
        {
          cs_int_t  buf[4] = {0, 0, 0, 0};
          buf[0] = strlen(cs_suite_nom_partie_fic) + 1;
          buf[1] = suite->nbr_fic;
          bft_file_write(buf, sizeof(cs_int_t), 4, suite->fic[0]);
          bft_file_write(cs_suite_nom_partie_fic, 1, buf[0], suite->fic[0]);
        }
        break;
      }

    }

  }

  /* Entête de la rubrique */
  /* --------------------- */

  if (cs_glob_base_rang <= 0) {

    switch(suite->type) {

    case CS_SUITE_TYPE_ASCII:

      bft_file_printf(suite->fic[0], "[%s\n", nom_rub);
      bft_file_printf(suite->fic[0],
                      " support = %3s, nbr_val_ent = %3d, typ_val = %3s]\n",
                      cs_suite_nom_support[(int)support],
                      (int)nbr_val_ent,
                      cs_suite_nom_typ_elt[(int)typ_val]);
      break;

    case CS_SUITE_TYPE_BINAIRE:

      {
        cs_int_t  buf[4];

        buf[0] = strlen(nom_rub) + 1;
        buf[1] = (int)support;
        buf[2] = (int)nbr_val_ent;
        buf[3] = (int)typ_val;

        bft_file_write(buf, sizeof(cs_int_t), 4, suite->fic[0]);
        bft_file_write(nom_rub, 1, buf[0], suite->fic[0]);
      }
      break;

    }

  }

  /* Traitement des renumérotations éventuelles */

  switch (support) {

  case CS_SUITE_SUPPORT_FAC_INT:
    val_tmp = cs_loc_suite_permute_ecr(mesh->n_i_faces,
                                       mesh->init_i_face_num,
                                       nbr_val_ent,
                                       typ_val,
                                       val);
    break;
  case CS_SUITE_SUPPORT_FAC_BRD:
    val_tmp = cs_loc_suite_permute_ecr(mesh->n_b_faces,
                                       mesh->init_b_face_num,
                                       nbr_val_ent,
                                       typ_val,
                                       val);
    break;
  default:
    val_tmp = NULL;

  }

  /* Contenu de la rubrique */
  /*------------------------*/

  nbr_val_tot = cs_loc_suite_calc_nbr_ent(suite,
                                          support,
                                          nbr_val_ent);

  /* En mode non-parallèle */

  if (cs_glob_base_rang < 0) {

    cs_loc_suite_ecr_val(suite, nbr_val_tot, typ_val,
                         (val_tmp != NULL) ? val_tmp : val,
                         &ind_col_0);
    cs_loc_suite_ecr_val_fin(suite, &ind_col_0);

  }

#if defined(_CS_HAVE_MPI)

  /* Variable indépendante du support et commune aux processus */

  else if (support == CS_SUITE_SUPPORT_SCAL) {

    if (cs_glob_base_rang == 0) {
      cs_loc_suite_ecr_val(suite, nbr_val_tot, typ_val,
                           (val_tmp != NULL) ? val_tmp : val,
                           &ind_col_0);
      cs_loc_suite_ecr_val_fin(suite, &ind_col_0);
    }

  }

  /* En mode parallèle sur une entité de maillage */

  else {

    cs_int_t  nbr_bloc;

    switch(support) {

    case CS_SUITE_SUPPORT_CEL:
      nbr_ent_glob = mesh->n_g_cells;
      nbr_ent_loc  = mesh->n_cells;
      num_glob_ent = mesh->global_cell_num;
      break;
    case CS_SUITE_SUPPORT_FAC_INT:
      nbr_ent_glob = mesh->n_g_i_faces;
      nbr_ent_loc  = mesh->n_i_faces;
      num_glob_ent = mesh->global_i_face_num;
      break;
    case CS_SUITE_SUPPORT_FAC_BRD:
      nbr_ent_glob = mesh->n_g_b_faces;
      nbr_ent_loc  = mesh->n_b_faces;
      num_glob_ent = mesh->global_b_face_num;
      break;
    case CS_SUITE_SUPPORT_SOM:
      nbr_ent_glob = mesh->n_g_vertices;
      nbr_ent_loc  = mesh->n_vertices;
      num_glob_ent = mesh->global_vtx_num;
      break;
    default:
      assert(   support > CS_SUITE_SUPPORT_SCAL
             && support <= CS_SUITE_SUPPORT_SOM);

    }

    nbr_bloc = (  ((sizeof(cs_real_t) * nbr_ent_glob * nbr_val_ent) - 1)
                / cs_suite_taille_buf_def) + 1;
    if (nbr_bloc > cs_glob_base_nbr)
      nbr_bloc = cs_glob_base_nbr;
    if (nbr_bloc == 0 )
      nbr_bloc = 1;

    cs_loc_suite_ecr_val_ent
      (suite,
       nbr_bloc,
       nbr_ent_glob,
       nbr_ent_loc,
       num_glob_ent,
       nbr_val_ent,
       typ_val,
       (val_tmp != NULL) ? val_tmp : (const cs_byte_t *)val);

  }

#endif /* #if defined(_CS_HAVE_MPI) */

  if (val_tmp != NULL)
    BFT_FREE (val_tmp);

}


/*----------------------------------------------------------------------------
 *  Fonction qui initialise l'API Fortran
 *----------------------------------------------------------------------------*/

void cs_suite_f77_api_init
(
 void
)
{
  cs_int_t ind;

  /* Allocation du tableau des pointeurs */

  cs_glob_suite_ptr_nbr = 10;
  BFT_MALLOC(cs_glob_suite_ptr_tab, cs_glob_suite_ptr_nbr, cs_suite_t *);

  /* Mise à zéro du tableau des pointeurs */

  for (ind = 0 ; ind < cs_glob_suite_ptr_nbr ; ind++)
    cs_glob_suite_ptr_tab[ind] = NULL;

}


/*----------------------------------------------------------------------------
 *  Fonction qui termine l'API Fortran
 *----------------------------------------------------------------------------*/

void cs_suite_f77_api_finalize
(
 void
)
{
  cs_int_t ind;

  /* Fermeture des fichiers qui ne le seraient pas */

  for (ind = 0 ; ind < cs_glob_suite_ptr_nbr ; ind++) {
    if (cs_glob_suite_ptr_tab[ind] != NULL)
      cs_suite_detruit (cs_glob_suite_ptr_tab[ind]);
  }

  /* Libération du tableau des pointeurs */

  cs_glob_suite_ptr_nbr = 0;
  BFT_FREE(cs_glob_suite_ptr_tab);

}


/*============================================================================
 *  Définitions de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Calcul du nombre de valeurs d'un enregistrement
 *----------------------------------------------------------------------------*/

static cs_int_t cs_loc_suite_calc_nbr_ent
(
 const cs_suite_t          *const suite,       /* --> Suite associée          */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent  /* --> Nb. val/point support   */
)
{
  switch(support) {

  case CS_SUITE_SUPPORT_SCAL:
    return (nbr_val_ent);
  case CS_SUITE_SUPPORT_CEL:
    return (suite->nbr_cel * nbr_val_ent);
  case CS_SUITE_SUPPORT_FAC_INT:
    return (suite->nbr_fac * nbr_val_ent);
  case CS_SUITE_SUPPORT_FAC_BRD:
    return (suite->nbr_fbr * nbr_val_ent);
  case CS_SUITE_SUPPORT_SOM:
    return (suite->nbr_som * nbr_val_ent);
  default:
    assert(support <= CS_SUITE_SUPPORT_SOM);

  }

  return 0;

}


/*----------------------------------------------------------------------------
 *  Calcul du déplacement dans un fichier correspondant aux valeurs
 *  d'une rubrique.
 *----------------------------------------------------------------------------*/

static bft_file_off_t cs_loc_suite_calc_avance
(
 const cs_suite_t          *const suite,       /* --> Suite associée          */
 const cs_suite_support_t         support,     /* --> Support de la variable  */
 const cs_int_t                   nbr_val_ent, /* --> Nb. val/point support   */
 const cs_type_t                  typ_val      /* --> Type de valeurs         */
)
{
  cs_int_t        nbr_val;
  bft_file_off_t  taille;

  assert(   suite->type == CS_SUITE_TYPE_ASCII
         || suite->type == CS_SUITE_TYPE_BINAIRE);
  assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);

  nbr_val = cs_loc_suite_calc_nbr_ent(suite,
                                      support,
                                      nbr_val_ent);

  if (suite->type == CS_SUITE_TYPE_ASCII) {

    cs_int_t     nbr_val_par_ligne, largeur_col;
    cs_int_t     nbr_ligne_complet, nbr_val_reste;

    switch (typ_val) {
    case CS_TYPE_cs_int_t:
      nbr_val_par_ligne = CS_SUITE_FMT_ASCII_NBR_COL_INT;
      largeur_col       = CS_SUITE_FMT_ASCII_DIM_COL_INT;
      break;
    case CS_TYPE_cs_real_t:
      nbr_val_par_ligne = CS_SUITE_FMT_ASCII_NBR_COL_REAL;
      largeur_col       = CS_SUITE_FMT_ASCII_DIM_COL_REAL;
      break;
    default:
      return 0;
    }

    nbr_ligne_complet = nbr_val / nbr_val_par_ligne;
    nbr_val_reste     = nbr_val % nbr_val_par_ligne;

    /* Nombre valeurs * largeur associée + fins des lignes entières */
    taille = (nbr_val * largeur_col) + nbr_ligne_complet;

    if (nbr_val_reste > 0)
    taille += 1;                           /* dernière fin de ligne */

  }
  else if (suite->type == CS_SUITE_TYPE_BINAIRE) {

    switch (typ_val) {
    case CS_TYPE_cs_int_t:
      taille = nbr_val * sizeof(cs_int_t);
      break;
    case CS_TYPE_cs_real_t:
      taille = nbr_val * sizeof(cs_real_t);
      break;
    default:
      return 0;
    }

  }

  return (taille);

}


/*----------------------------------------------------------------------------
 *  Fonction qui initialise une structure fichier associée à une suite ;
 *  Le nom de fichier associé correspond au nom de la suite pour le
 *  premier fichier, prolongé éventuellement par '_pxx' pour les fichiers
 *  successifs en cas de découpage du fichier en plusieurs morceaux pour
 *  des gros volumes de données.
 *  On indique en retour si les données se poursuivent sur un autre fichier.
 *----------------------------------------------------------------------------*/

static cs_bool_t cs_loc_suite_fic_ajoute
(
 cs_suite_t       *const suite              /* --> Structure suite            */
)
{
/* On ne change le numero que si les fichiers suite des
 * versions precedentes ne sont plus compatibles. Ce numero est en
 * fait deconnecte du numero de version */
  char entete_txt[] = "Code_Saturne 1.1 (reprise)\n";
  char entete_bin[] = "Code_Saturne_1.1_bin_reprise\n";
  char buffer[CS_SUITE_LNG_BUF_ASCII + 1];

  char *nom_fic = NULL;

  cs_int_t ind_fic = 0;

  bft_file_mode_t fic_mode;
  bft_file_type_t fic_type;

  cs_bool_t autre_partie = CS_FALSE;


  /* Instructions */

  if (cs_glob_base_rang > 0)
    return CS_FALSE;

  /* Incrémentation du compteur */

  suite->nbr_fic += 1;

  /* Construction du nom du fichier suite */

  if (suite->nbr_fic == 1)
    nom_fic = suite->nom;
  else {
    BFT_MALLOC(nom_fic, strlen(suite->nom) + strlen("_pxx") + 1, char);
    sprintf(nom_fic, "%s_p%02d", suite->nom, suite->nbr_fic);
  }


  /* Création du descripteur de fichier d'interface */
  /*------------------------------------------------*/

  /*
    En lecture, on détecte le type binaire/texte automatiquement
    (en commencant par une tentative en mode binaire) ; en cas de
    découpage en plusieurs fichiers pour gros volumes de données,
    tous les fichiers doivent être de même type.
    En écriture, on choisit le type.
  */

  if (suite->mode == CS_SUITE_MODE_LECTURE) {
    fic_mode = BFT_FILE_MODE_READ;
    fic_type = BFT_FILE_TYPE_BINARY;
  }

  else {
    fic_mode = BFT_FILE_MODE_WRITE;
    if (suite->type == BFT_FILE_TYPE_TEXT)
      fic_type = BFT_FILE_TYPE_TEXT;
    else
      fic_type = BFT_FILE_TYPE_BINARY;
  }


  /* Création du fichier */

  if (suite->mode == CS_SUITE_MODE_LECTURE) {
    ind_fic = suite->nbr_fic - 1;
    BFT_REALLOC(suite->fic, suite->nbr_fic, bft_file_t *);
  }

  else if (suite->mode == CS_SUITE_MODE_ECRITURE) {
    ind_fic = 0;
    if (suite->fic == NULL)
      BFT_MALLOC(suite->fic, 1, bft_file_t *);
    else if (suite->fic[0] != NULL)
      bft_file_free(suite->fic[0]);
  }

  suite->fic[ind_fic] = bft_file_open(nom_fic,
                                      fic_mode,
                                      fic_type);

  if (fic_type == BFT_FILE_TYPE_BINARY)
    bft_file_set_big_endian(suite->fic[ind_fic]);


  /* Gestion des en-têtes selon le mode d'ouverture */
  /*------------------------------------------------*/

  /* Lecture de l'en-tête */

  if (suite->mode == CS_SUITE_MODE_LECTURE) {

    cs_int_t nbr_rec, nbr_lus;
    cs_int_t num_ligne = 1;

    char     *str;

    cs_suite_type_t suite_type_detect = CS_SUITE_TYPE_BINAIRE;


    /* Détection automatique de type de fichier */

    nbr_rec = strlen(entete_bin);
    nbr_lus = bft_file_read_try(buffer, 1, nbr_rec, suite->fic[ind_fic]);

    if (strncmp(entete_bin, buffer, nbr_rec) != 0) {

      assert(   CS_SUITE_LNG_BUF_ASCII > strlen(entete_txt)
             && CS_SUITE_LNG_BUF_ASCII > strlen(entete_bin));

      bft_file_free(suite->fic[ind_fic]);
      fic_type = BFT_FILE_TYPE_TEXT;
      suite->fic[ind_fic] = bft_file_open(nom_fic,
                                          fic_mode,
                                          fic_type);

      nbr_rec = strlen(entete_txt);
      str = bft_file_gets_try(buffer,
                              CS_SUITE_LNG_BUF_ASCII,
                              suite->fic[ind_fic],
                              &num_ligne);

      if (str != NULL) {
        nbr_rec = strlen(entete_txt);
        if (strncmp(entete_txt, buffer, nbr_rec) == 0)
          suite_type_detect = CS_SUITE_TYPE_ASCII;
      }

      if (suite_type_detect != CS_SUITE_TYPE_ASCII)
        bft_error(__FILE__, __LINE__, 0,
                  _("The file <%s> is not a valid restart file.\n"),
                  nom_fic);

    }

    if (suite->type != suite_type_detect) {
      if (suite->nbr_fic > 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("The restart file <%s> and its extension\n"
                    "<%s> are not of the same type (text/binary);\n"
                    "they do not correspond to the same data set.\n"),
                  suite->nom, nom_fic);
      else
        suite->type = suite_type_detect;
    }

    if (suite->type == CS_SUITE_TYPE_ASCII)
      autre_partie = cs_loc_suite_analyse_txt(suite);
    else if (suite->type == CS_SUITE_TYPE_BINAIRE)
      autre_partie = cs_loc_suite_analyse_bin(suite);

  }

  /* Écriture de l'en-tête */

  else if (suite->mode == CS_SUITE_MODE_ECRITURE) {

    switch(suite->type) {

    case CS_SUITE_TYPE_ASCII:

      bft_file_printf(suite->fic[0], entete_txt);
      bft_file_printf(suite->fic[0], "[nbr_cel     = %d]\n",
                      (int)cs_glob_mesh->n_g_cells);
      bft_file_printf(suite->fic[0], "[nbr_fac_int = %d]\n",
                      (int)cs_glob_mesh->n_g_i_faces);
      bft_file_printf(suite->fic[0], "[nbr_fac_brd = %d]\n",
                      (int)cs_glob_mesh->n_g_b_faces);
      bft_file_printf(suite->fic[0], "[nbr_som     = %d]\n",
                      (int)cs_glob_mesh->n_g_vertices);
      break;

    case CS_SUITE_TYPE_BINAIRE:

      {
        cs_int_t  buf[4];

        buf[0] = cs_glob_mesh->n_g_cells;
        buf[1] = cs_glob_mesh->n_g_i_faces;
        buf[2] = cs_glob_mesh->n_g_b_faces;
        buf[3] = cs_glob_mesh->n_g_vertices;

        bft_file_write(entete_bin, 1, strlen(entete_bin),suite->fic[0]);
        bft_file_write(buf, sizeof(cs_int_t), 4, suite->fic[0]);
      }
      break;

    default:
      assert (   suite->type == CS_SUITE_TYPE_ASCII
              || suite->type == CS_SUITE_TYPE_BINAIRE);

    }

  }

  /* Libération mémoire */

  if (nom_fic != suite->nom)
    BFT_FREE(nom_fic);

  /* On indique en retour si le fichier continue */

  return autre_partie;

}


/*----------------------------------------------------------------------------
 *  Fonction qui lit une liste de valeurs sur un fichier ; on suppose que
 *  l'on est déjà à la bonne position de départ dans le fichier.
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_lit_val
(
 const cs_suite_type_t         type,           /* --> Structure suite         */
 const bft_file_t       *const fic,            /* --> Fichier associé         */
 const cs_int_t                nbr_val,        /* --> Nb. valeurs à lire      */
 const cs_type_t               typ_val,        /* --> Type de valeurs         */
       void             *const val,            /* <-- Valeurs à lire          */
       char                    asc_buffer[]    /* <-> Buffer texte            */
)
{

  /* Traitement des fichiers binaires (simple) */

  if (type == CS_SUITE_TYPE_BINAIRE) {

    cs_int_t  typ_taille = 1;

    switch (typ_val) {
    case CS_TYPE_cs_int_t:
      typ_taille = sizeof(cs_int_t);
      break;
    case CS_TYPE_cs_real_t:
      typ_taille = sizeof(cs_real_t);
      break;
    default:
      assert (   typ_val != CS_TYPE_cs_int_t
              || typ_val != CS_TYPE_cs_real_t);
    }

    bft_file_read(val, typ_taille, nbr_val, fic);

  }

  /* Traitement des fichiers formatés (plus complexe) */

  else if (type == CS_SUITE_TYPE_ASCII) {

    char      fmt_lec[4];
    cs_int_t  ind_val;

    char      *str = asc_buffer;
    cs_int_t  num_ligne = 1;

    fmt_lec[0] = '%';

    switch (typ_val) {

    case CS_TYPE_cs_int_t:
      {
        const cs_int_t *const tab_val_int = (const cs_int_t *)val;

        /* Attention au format de lecture (longueur des types) ! */

        if (sizeof(cs_int_t) == sizeof(int))
          fmt_lec[1] = 'd', fmt_lec[2] = '\0';
        else if(sizeof(cs_int_t) == sizeof(long))
          fmt_lec[1] = 'l', fmt_lec[2] = 'd', fmt_lec[3] = '\0';
        else
          fmt_lec[1] = 'L', fmt_lec[2] = 'd', fmt_lec[3] = '\0';

        /* Lecture effective */

        for (ind_val = 0 ; ind_val < nbr_val ; ind_val++) {

          char col_deb = *str;

          if (col_deb == '\n' || col_deb == '\0')
            str = bft_file_gets(asc_buffer,
                                CS_SUITE_LNG_BUF_ASCII,
                                fic,
                                &num_ligne);

          else {
            if (col_deb != ' ' && col_deb != '\t')
              bft_error(__FILE__, __LINE__, 0,
                        _("The restart file <%s> and its extension\n"
                          "<%s> are not of the same type (text/binary);\n"
                          "they do not correspond to the same data set.\n"),
                        CS_SUITE_FMT_ASCII_DIM_COL_INT);
          }

          if (str == NULL)
            break;

          if (sscanf(str, fmt_lec, &(tab_val_int[ind_val])) != 1)
            bft_error(__FILE__, __LINE__, 0,
                      _("Error while reading an array of integers in a restart file:\n"
                        "%ul on %ul read values"),
                      (unsigned long)ind_val, (unsigned long)nbr_val);

          str += CS_SUITE_FMT_ASCII_DIM_COL_INT;

        }
        if (ind_val < nbr_val)
          bft_error(__FILE__, __LINE__, 0,
                    _("Error while reading an array of integers in a restart file:\n"
                      "%ul on %ul read values"),
                    (unsigned long)ind_val, (unsigned long)nbr_val);

      }

      break;

    case CS_TYPE_cs_real_t:
      {
        const cs_real_t  *const tab_val_real = (const cs_real_t *)val;

        /* Attention au format de lecture (longueur des types) ! */

        if (sizeof(cs_real_t) == sizeof(float))
          fmt_lec[1] = 'g', fmt_lec[2] = '\0';
        else if (sizeof(cs_real_t) == sizeof(double))
          fmt_lec[1] = 'l', fmt_lec[2] = 'g', fmt_lec[3] = '\0';
        else
          fmt_lec[1] = 'L', fmt_lec[2] = 'g', fmt_lec[3] = '\0';

        /* Lecture effective */

        for (ind_val = 0 ; ind_val < nbr_val ; ind_val++) {

          char col_deb = *str;

          if (col_deb == '\n' || col_deb == '\0')
            str = bft_file_gets(asc_buffer,
                                CS_SUITE_LNG_BUF_ASCII,
                                fic,
                                &num_ligne);

          else {
            if (col_deb != ' ' && col_deb != '\t')
              bft_error(__FILE__, __LINE__, 0,
                        _("Error while reading an array of reals in a restart file:\n"
                          "field of different size than %d"),
                        CS_SUITE_FMT_ASCII_DIM_COL_REAL);
          }

          if (str == NULL)
            break;

          if (sscanf(str, fmt_lec, &(tab_val_real[ind_val])) != 1)
            bft_error(__FILE__, __LINE__, 0,
                      _("Error while reading an array of reals in a restart file:\n"
                        "%ul on %ul read values"),
                      (unsigned long)ind_val, (unsigned long)nbr_val);

          str += CS_SUITE_FMT_ASCII_DIM_COL_REAL;

        }
        if (ind_val < nbr_val)
          bft_error(__FILE__, __LINE__, 0,
                    _("Error while reading an array of reals in a restart file:\n"
                      "%ul on %ul read values"),
                    (unsigned long)ind_val, (unsigned long)nbr_val);

      }
      break;

    default:

      assert(   typ_val != CS_TYPE_cs_int_t
             || typ_val != CS_TYPE_cs_real_t);

    } /* Fin : switch (typ_val) */

  }

}


/*----------------------------------------------------------------------------
 *  Fonction qui écrit une liste de valeurs sur un fichier
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_ecr_val
(
 const cs_suite_t          *const suite,       /* --> Ptr. fichier suite      */
 const cs_int_t                   nbr_val,     /* --> Nb. valeurs à écrire    */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
 const void                *const val,         /* --> Valeurs à écrire        */
       cs_int_t            *const ind_col      /* <-> Colonne départ/fin      */
)
{

  /* Traitement des fichiers binaires (simple) */

  if (suite->type == CS_SUITE_TYPE_BINAIRE) {

    cs_int_t  typ_taille = 1;

    switch (typ_val) {
    case CS_TYPE_cs_int_t:
      typ_taille = sizeof(cs_int_t);
      break;
    case CS_TYPE_cs_real_t:
      typ_taille = sizeof(cs_real_t);
      break;
    default:
      assert(   typ_val != CS_TYPE_cs_int_t
             || typ_val != CS_TYPE_cs_real_t);
    }

    bft_file_write(val, typ_taille, nbr_val, suite->fic[0]);

  }


  /* Traitement des fichiers formatés (plus complexe) */

  else if (suite->type == CS_SUITE_TYPE_ASCII) {

    cs_int_t ind_val;

    switch (typ_val) {

    case CS_TYPE_cs_int_t:
      {
        const cs_int_t *const tab_val_int = (const cs_int_t *)val;

        for (ind_val = 0 ; ind_val < nbr_val ; ind_val++) {

          if (*ind_col == CS_SUITE_FMT_ASCII_NBR_COL_INT) {
            *ind_col = 0;
            bft_file_printf(suite->fic[0], "\n");
          }
          *ind_col += 1;

          bft_file_printf(suite->fic[0], cs_suite_fmt_ascii_tab_int,
                          tab_val_int[ind_val]);

        }

      }
      break;

    case CS_TYPE_cs_real_t:
      {
        const cs_real_t  *const tab_val_real = (const cs_real_t *)val;

        for (ind_val = 0 ; ind_val < nbr_val ; ind_val++) {

          if (*ind_col == CS_SUITE_FMT_ASCII_NBR_COL_REAL) {
            *ind_col = 0;
            bft_file_printf(suite->fic[0], "\n");
          }
          *ind_col += 1;

          bft_file_printf(suite->fic[0], cs_suite_fmt_ascii_tab_real,
                          tab_val_real[ind_val]);

        }

      }
      break;

    default:

      assert(   typ_val != CS_TYPE_cs_int_t
             || typ_val != CS_TYPE_cs_real_t);

    } /* Fin : switch (typ_val) */

  }

}


/*----------------------------------------------------------------------------
 *  Fonction qui termine l'écriture d'une liste de valeurs sur un fichier
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_ecr_val_fin
(
 const cs_suite_t          *const suite,       /* --> Ptr. fichier suite      */
       cs_int_t            *const ind_col      /* <-> Colonne départ/fin      */
)
{
  /* Rien à faire pour les fichiers binaires */

  /* Traitement des fichiers formatés */

  if (suite->type == CS_SUITE_TYPE_ASCII) {

    bft_file_printf(suite->fic[0], "\n");

    *ind_col = 0;

  }

}


#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Fonction qui créé des listes locales->globales associées à chaque bloc.
 *
 *  Pour le processus "maître E/S", le tableau "nbr_ent_bloc" indiquant le
 *  nombre d'entités appartenant à chaque bloc, et ceci pour chaque domaine :
 *  nbr_ent(bloc, domaine) = nbr_ent_bloc[nbr_bloc * ind_dom + ind_bloc] ;
 *  Pour les autres processus, le tableau "nbr_ent_loc" indique le nombre
 *  d'entités locales appartenant à chaque bloc.
 *  Les tableaux "lst_ent_loc" et "lst_ent_glob" donnent les indices locaux
 *  et globaux des entités locales appartenant à chaque bloc (nbr_ent_loc[0]
 *  premières valeurs pour le bloc 0, nbr_ent_loc[1] suivantes pour le
 *  second, etc.)
 *
 *  Ces tableaux sont alloués ici et devront être libérés par la suite.
 *
 *  On renvoie la plus grande valeur (globale) de nbr_ent_loc.
 *----------------------------------------------------------------------------*/

static cs_int_t cs_loc_suite_cree_listes_ent
(
 const cs_int_t            nbr_bloc,      /* --> Nombre de blocs              */
 const cs_int_t            nbr_ent_glob,  /* --> Nombre global d'entités      */
 const cs_int_t            nbr_ent_loc,   /* --> Nombre local d'entités       */
 const fvm_gnum_t   *const num_glob_ent,  /* --> Numéros globaux des entités  */
       cs_int_t     *const pas_bloc,      /* <-- Taille blocs (sauf dernier)  */
       cs_int_t   * *const nbr_ent_bloc,  /* <-- Nombre d'entités par bloc    */
       cs_int_t   * *const lst_ent_loc,   /* <-- Liste des entités par bloc   */
       cs_int_t   * *const lst_ent_glob   /* <-- Liste des entités par bloc   */
)
{

  cs_int_t  ind, ind_bloc, ind_glob_ent;
  cs_int_t  taille_max_loc, taille_max_glob;

  cs_int_t  *pos_ent_bloc;
  cs_int_t  *nbr_ent_bloc_loc;

  /* pas_bloc = ceil(nbr_ent_glob/nbr_bloc) */

  *pas_bloc = nbr_ent_glob / nbr_bloc;
  if (nbr_ent_glob % nbr_bloc > 0)
    *pas_bloc += 1;

  /* Allocation */

  BFT_MALLOC(*lst_ent_loc,  nbr_ent_loc, cs_int_t);
  BFT_MALLOC(*lst_ent_glob, nbr_ent_loc, cs_int_t);

  if (cs_glob_base_rang == 0) {
    BFT_MALLOC(*nbr_ent_bloc, nbr_bloc * cs_glob_base_nbr, cs_int_t);
    BFT_MALLOC(nbr_ent_bloc_loc, nbr_bloc, cs_int_t);
  }
  else {
    BFT_MALLOC(*nbr_ent_bloc, nbr_bloc, cs_int_t);
    nbr_ent_bloc_loc = *nbr_ent_bloc;
  }

  BFT_MALLOC(pos_ent_bloc, nbr_bloc, cs_int_t);


  /* Comptage */

  for (ind_bloc = 0 ; ind_bloc < nbr_bloc ; ind_bloc++)
    nbr_ent_bloc_loc[ind_bloc] = 0;

  for (ind = 0 ; ind < nbr_ent_loc ; ind++) {
    ind_bloc = (num_glob_ent[ind] - 1) / *pas_bloc;
    nbr_ent_bloc_loc[ind_bloc] += 1;
  }

  /* Affectation */

  pos_ent_bloc[0] = 0;
  for (ind_bloc = 1 ; ind_bloc < nbr_bloc ; ind_bloc++)
    pos_ent_bloc[ind_bloc]
      = pos_ent_bloc[ind_bloc - 1] + nbr_ent_bloc_loc[ind_bloc - 1];

  for (ind = 0 ; ind < nbr_ent_loc ; ind++) {
    ind_glob_ent = num_glob_ent[ind] - 1;
    ind_bloc = ind_glob_ent / *pas_bloc;
    (*lst_ent_loc) [pos_ent_bloc[ind_bloc]] = ind;
    (*lst_ent_glob)[pos_ent_bloc[ind_bloc]] = ind_glob_ent % (*pas_bloc);
    pos_ent_bloc[ind_bloc] += 1;
  }

  BFT_FREE(pos_ent_bloc);

  /*
   * Dimensionnement des messages ; le processus "maître E/S" stocke
   * le nombre d'entités par bloc et par domaine externe :
   * nbr_ent_bloc[ind_dom][ind_bloc]
   */

  MPI_Gather(nbr_ent_bloc_loc, nbr_bloc, CS_MPI_INT,
             *nbr_ent_bloc, nbr_bloc, CS_MPI_INT, 0,
             cs_glob_base_mpi_comm);

  if (cs_glob_base_rang == 0)
    BFT_FREE(nbr_ent_bloc_loc);
  else
    nbr_ent_bloc_loc = NULL;

  /* Tailles des divers messages */

  taille_max_loc = (*nbr_ent_bloc)[0];
  for (ind_bloc = 1 ; ind_bloc < nbr_bloc ; ind_bloc++)
    taille_max_loc = CS_MAX(taille_max_loc, (*nbr_ent_bloc)[ind_bloc]);

  MPI_Allreduce(&taille_max_loc, &taille_max_glob, 1, CS_MPI_INT, MPI_MAX,
                cs_glob_base_mpi_comm);

  return taille_max_glob;

}


/*----------------------------------------------------------------------------
 *  Fonction qui lit les valeurs associées à une variable définie sur
 *  une entité de maillage.
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_lit_val_ent
(
 const cs_suite_t  *const  suite,         /* --> Structure suite              */
 const bft_file_t  *const  fic,           /* --> Fichier associé aux valeurs  */
 const cs_int_t            nbr_bloc,      /* --> Nombre de blocs              */
 const cs_int_t            nbr_ent_glob,  /* --> Nombre global d'entités      */
 const cs_int_t            nbr_ent_loc,   /* --> Nombre local d'entités       */
 const fvm_gnum_t   *const num_glob_ent,  /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
       cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
)
{

  char      buffer_ascii[CS_SUITE_LNG_BUF_ASCII + 1];

  cs_int_t  pas_bloc, ind, ind_dom, ind_bloc, ind_loc;
  cs_int_t  ind_deb_glob_bloc, ind_fin_glob_bloc, nbr_ent_glob_bloc;
  cs_int_t  ind_deb_loc_bloc, nbr_ent_bloc_cur;
  cs_int_t  nbr_bloc_max;

  cs_int_t  *nbr_ent_bloc;
  cs_int_t  *lst_ent_loc;
  cs_int_t  *lst_ent_glob;
  cs_int_t  *ind_ent_bloc;

  cs_int_t   *buffer_ent_bloc = NULL;

  cs_byte_t  *buffer_fic = NULL, *buffer_msg = NULL;

  size_t      ind_byte, nbr_byte_ent;

  MPI_Datatype  mpi_typ;
  MPI_Status    status;

  /* Initialisations */

  buffer_ascii[0] = '\0';

  switch (typ_val) {
  case CS_TYPE_cs_int_t:
    mpi_typ      = CS_MPI_INT;
    nbr_byte_ent = nbr_val_ent * sizeof(cs_int_t);
    break;
  case CS_TYPE_cs_real_t:
    mpi_typ      = CS_MPI_REAL;
    nbr_byte_ent = nbr_val_ent * sizeof(cs_real_t);
    break;
    break;
  default:
    assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);
  }

  /* Création des listes d'entités associées à la redistribution */
  /*-------------------------------------------------------------*/

  nbr_bloc_max = cs_loc_suite_cree_listes_ent(nbr_bloc,
                                              nbr_ent_glob,
                                              nbr_ent_loc,
                                              num_glob_ent,
                                              &pas_bloc,
                                              &nbr_ent_bloc,
                                              &lst_ent_loc,
                                              &lst_ent_glob);

  /* Allocation des tableaux de travail */

  if (cs_glob_base_rang == 0)
    BFT_MALLOC(buffer_ent_bloc, nbr_bloc_max, cs_int_t);

  BFT_MALLOC(buffer_fic, pas_bloc     * nbr_byte_ent, cs_byte_t);
  BFT_MALLOC(buffer_msg, nbr_bloc_max * nbr_byte_ent, cs_byte_t);

  /* Boucle sur les blocs */
  /*----------------------*/

  ind_deb_glob_bloc = 0;
  ind_deb_loc_bloc  = 0;

  for (ind_bloc = 0 ; ind_bloc < nbr_bloc ; ind_bloc++) {

    ind_fin_glob_bloc = CS_MIN(ind_deb_glob_bloc + pas_bloc, nbr_ent_glob);

    /* Traitement processus "maître E/S" */

    if (cs_glob_base_rang == 0) {

      nbr_ent_glob_bloc = ind_fin_glob_bloc - ind_deb_glob_bloc;

      /* Lecture sur fichier */

      cs_loc_suite_lit_val(suite->type, fic, nbr_ent_glob_bloc * nbr_val_ent,
                           typ_val, (void *)buffer_fic, buffer_ascii);

      for (ind_dom = 0 ; ind_dom < cs_glob_base_nbr ; ind_dom++) {

        /* Index des valeurs */

        nbr_ent_bloc_cur = nbr_ent_bloc[ind_dom * nbr_bloc + ind_bloc];

        if (ind_dom > 0) {
          if (nbr_ent_bloc_cur > 0) {
            MPI_Recv((void *)buffer_ent_bloc, nbr_ent_bloc_cur, CS_MPI_INT,
                     ind_dom, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm,
                     &status);
            ind_ent_bloc = buffer_ent_bloc;
          }
        }
        else /* if (ind_dom == 0) */
          ind_ent_bloc = lst_ent_glob + ind_deb_loc_bloc;

        /* Envoi des valeurs indexées */

        if (ind_dom > 0) {
          if (nbr_ent_bloc_cur > 0) {
            ind = 0;
            for (ind_loc = 0 ; ind_loc < nbr_ent_bloc_cur ; ind_loc++) {
              for (ind_byte = 0 ; ind_byte < nbr_byte_ent ; ind_byte++)
                buffer_msg[ind++]
                  = buffer_fic[ind_ent_bloc[ind_loc] * nbr_byte_ent + ind_byte];
            }
            MPI_Send((void *)buffer_msg,
                     nbr_ent_bloc_cur * nbr_val_ent, mpi_typ,
                     ind_dom, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm);
          }
        }
        else { /* if (ind_dom == 0) */
          for (ind_loc = 0 ; ind_loc < nbr_ent_bloc_cur ; ind_loc++) {
            for (ind_byte = 0 ; ind_byte < nbr_byte_ent ; ind_byte++)
              tab_val[((lst_ent_loc + ind_deb_loc_bloc)[ind_loc])
                        * nbr_byte_ent + ind_byte]
                = buffer_fic[ind_ent_bloc[ind_loc] * nbr_byte_ent + ind_byte];
          }
        }

      }

    }

    /* Traitement autres processus */

    else if (nbr_ent_bloc[ind_bloc] > 0) {

      /* Index des valeurs */

      ind_ent_bloc = lst_ent_glob + ind_deb_loc_bloc;

      MPI_Send((void *)ind_ent_bloc,
               nbr_ent_bloc[ind_bloc], CS_MPI_INT,
               0, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm);

      /* Remplissage et envoi des valeurs indexées */

      ind_ent_bloc = lst_ent_loc + ind_deb_loc_bloc;

      MPI_Recv((void *)buffer_msg, nbr_ent_bloc[ind_bloc] * nbr_val_ent,
               mpi_typ, 0, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm, &status);
      ind = 0;
      for (ind_loc = 0 ; ind_loc < nbr_ent_bloc[ind_bloc] ; ind_loc++) {
        for (ind_byte = 0 ; ind_byte < nbr_byte_ent ; ind_byte++)
          tab_val[ind_ent_bloc[ind_loc] * nbr_byte_ent + ind_byte]
            = buffer_msg[ind++];
      }

    }

    ind_deb_glob_bloc = ind_fin_glob_bloc;
    ind_deb_loc_bloc += nbr_ent_bloc[ind_bloc];

  }

  /* Libération des tableaux de travail */

  BFT_FREE(buffer_fic);
  BFT_FREE(buffer_msg);

  if (cs_glob_base_rang == 0)
    BFT_FREE(buffer_ent_bloc);

  BFT_FREE(nbr_ent_bloc);
  BFT_FREE(lst_ent_loc);
  BFT_FREE(lst_ent_glob);

}


/*----------------------------------------------------------------------------
 *  Fonction qui écrit les valeurs associées à une variable définie sur
 *  une entité de maillage.
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_ecr_val_ent
(
 const cs_suite_t  *const  suite,         /* --> Structure suite              */
 const cs_int_t            nbr_bloc,      /* --> Nombre de blocs              */
 const cs_int_t            nbr_ent_glob,  /* --> Nombre global d'entités      */
 const cs_int_t            nbr_ent_loc,   /* --> Nombre local d'entités       */
 const fvm_gnum_t   *const num_glob_ent,  /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
 const cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
)
{

  cs_int_t  pas_bloc, ind, ind_dom, ind_bloc, ind_loc;
  cs_int_t  ind_deb_glob_bloc, ind_fin_glob_bloc, nbr_ent_glob_bloc;
  cs_int_t  ind_deb_loc_bloc, nbr_ent_bloc_cur;
  cs_int_t  nbr_bloc_max;

  cs_int_t  *nbr_ent_bloc;
  cs_int_t  *lst_ent_loc;
  cs_int_t  *lst_ent_glob;
  cs_int_t  *ind_ent_bloc;

  cs_int_t    ind_col = 0;
  cs_int_t   *buffer_ent_bloc = NULL;

  cs_byte_t  *buffer_fic = NULL, *buffer_msg = NULL;

  size_t      ind_byte, nbr_byte_ent;

  MPI_Datatype  mpi_typ;
  MPI_Status    status;

  /* Initialisations */

  switch (typ_val) {
  case CS_TYPE_cs_int_t:
    mpi_typ      = CS_MPI_INT;
    nbr_byte_ent = nbr_val_ent * sizeof(cs_int_t);
    break;
  case CS_TYPE_cs_real_t:
    mpi_typ      = CS_MPI_REAL;
    nbr_byte_ent = nbr_val_ent * sizeof(cs_real_t);
    break;
  default:
    assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);
  }

  /* Création des listes d'entités associées à la redistribution */
  /*-------------------------------------------------------------*/

  nbr_bloc_max = cs_loc_suite_cree_listes_ent(nbr_bloc,
                                              nbr_ent_glob,
                                              nbr_ent_loc,
                                              num_glob_ent,
                                              &pas_bloc,
                                              &nbr_ent_bloc,
                                              &lst_ent_loc,
                                              &lst_ent_glob);

  /* Allocation des tableaux de travail */

  if (cs_glob_base_rang == 0)
    BFT_MALLOC(buffer_ent_bloc, nbr_bloc_max, cs_int_t);

  BFT_MALLOC(buffer_fic, pas_bloc     * nbr_byte_ent, cs_byte_t);
  BFT_MALLOC(buffer_msg, nbr_bloc_max * nbr_byte_ent, cs_byte_t);

  /* Boucle sur les blocs */
  /*----------------------*/

  ind_deb_glob_bloc = 0;
  ind_deb_loc_bloc  = 0;

  for (ind_bloc = 0 ; ind_bloc < nbr_bloc ; ind_bloc++) {

    ind_fin_glob_bloc = CS_MIN(ind_deb_glob_bloc + pas_bloc, nbr_ent_glob);

    /* Traitement processus "maître E/S" */

    if (cs_glob_base_rang == 0) {

      nbr_ent_glob_bloc = ind_fin_glob_bloc - ind_deb_glob_bloc;

      for (ind_dom = 0 ; ind_dom < cs_glob_base_nbr ; ind_dom++) {

        nbr_ent_bloc_cur = nbr_ent_bloc[ind_dom * nbr_bloc + ind_bloc];

        /* Synchronisation forcée pour éviter que tous les autres
           rangs ne postent des envois en même temps (pouvant
           mener à un problème d'allocation de buffers) */

        if (   cs_glob_suite_sync_gather == CS_TRUE
            && ind_dom > 0 && nbr_ent_bloc_cur > 0) {
          int _ind_dom = ind_dom;
          MPI_Send(&_ind_dom, 1, MPI_INT, ind_dom,
                   CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm);
        }

        /* Index des valeurs */

        if (ind_dom > 0) {
          if (nbr_ent_bloc_cur > 0) {
            MPI_Recv((void *)buffer_ent_bloc, nbr_ent_bloc_cur, CS_MPI_INT,
                     ind_dom, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm,
                     &status);
            ind_ent_bloc = buffer_ent_bloc;
          }
        }
        else /* if (ind_dom == 0) */
          ind_ent_bloc = lst_ent_glob + ind_deb_loc_bloc;

        /* Récupération des valeurs indexées et préparation écriture */

        if (ind_dom > 0) {
          if (nbr_ent_bloc_cur > 0) {
            MPI_Recv((void *)buffer_msg,
                     nbr_ent_bloc_cur * nbr_val_ent, mpi_typ,
                     ind_dom, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm,
                     &status);
            ind = 0;
            for (ind_loc = 0 ; ind_loc < nbr_ent_bloc_cur ; ind_loc++) {
              for (ind_byte = 0 ; ind_byte < nbr_byte_ent ; ind_byte++)
                buffer_fic[ind_ent_bloc[ind_loc] * nbr_byte_ent + ind_byte]
                  = buffer_msg[ind++];
            }
          }
        }
        else { /* if (ind_dom == 0) */
          for (ind_loc = 0 ; ind_loc < nbr_ent_bloc_cur ; ind_loc++) {
            for (ind_byte = 0 ; ind_byte < nbr_byte_ent ; ind_byte++)
              buffer_fic[ind_ent_bloc[ind_loc] * nbr_byte_ent + ind_byte]
                = tab_val[((lst_ent_loc + ind_deb_loc_bloc)[ind_loc])
                          * nbr_byte_ent + ind_byte];
          }
        }

      }

      /* Écriture sur fichier */

      cs_loc_suite_ecr_val(suite, nbr_ent_glob_bloc * nbr_val_ent,
                           typ_val, (void *)buffer_fic, &ind_col);

    }

    /* Traitement autres processus */

    else if (nbr_ent_bloc[ind_bloc] > 0) {

      /* Remplissage des valeurs indexées */

      ind_ent_bloc = lst_ent_loc + ind_deb_loc_bloc;

      ind = 0;
      for (ind_loc = 0 ; ind_loc < nbr_ent_bloc[ind_bloc] ; ind_loc++) {
        for (ind_byte = 0 ; ind_byte < nbr_byte_ent ; ind_byte++)
          buffer_msg[ind++]
            = tab_val[ind_ent_bloc[ind_loc] * nbr_byte_ent + ind_byte];
      }

      /* Synchronisation forcée pour éviter que tous les autres
         rangs ne postent des envois en même temps (pouvant
         mener à un problème d'allocation de buffers) */

      if (cs_glob_suite_sync_gather == CS_TRUE) {
        int _ind_sync;
        MPI_Recv(&_ind_sync, 1, MPI_INT, 0,
                 CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm, &status);
      }

      /* Index des valeurs */

      ind_ent_bloc = lst_ent_glob + ind_deb_loc_bloc;

      MPI_Send((void *)ind_ent_bloc,
               nbr_ent_bloc[ind_bloc], CS_MPI_INT,
               0, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm);

      /* Envoi des valeurs indexées */

      MPI_Send((void *)buffer_msg, nbr_ent_bloc[ind_bloc] * nbr_val_ent,
               mpi_typ, 0, CS_SUITE_MPI_TAG, cs_glob_base_mpi_comm);

    }

    ind_deb_glob_bloc = ind_fin_glob_bloc;
    ind_deb_loc_bloc += nbr_ent_bloc[ind_bloc];

  }

  /* Retour à la ligne pour fichier ASCII */

  if (cs_glob_base_rang == 0)
    cs_loc_suite_ecr_val_fin(suite, &ind_col);

  /* Libération des tableaux de travail */

  BFT_FREE(buffer_fic);
  BFT_FREE(buffer_msg);

  if (cs_glob_base_rang == 0)
    BFT_FREE(buffer_ent_bloc);

  BFT_FREE(nbr_ent_bloc);
  BFT_FREE(lst_ent_loc);
  BFT_FREE(lst_ent_glob);

}


/*----------------------------------------------------------------------------
 *  Fonction qui distribue une liste de variables
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_distr_val
(
 const cs_int_t                   nbr_val,     /* --> Nb. valeurs             */
 const cs_type_t                  typ_val,     /* --> Type de valeurs         */
       void                *const val          /* <-> Valeurs à lire          */
)
{

  cs_int_t      ind_dom;

  MPI_Datatype  mpi_typ;
  MPI_Status    status;

  /* Initialisations */

  switch (typ_val) {
  case CS_TYPE_cs_int_t:
    mpi_typ      = CS_MPI_INT;
    break;
  case CS_TYPE_cs_real_t:
    mpi_typ      = CS_MPI_REAL;
    break;
  default:
    assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);
  }

  if (cs_glob_base_rang == 0) {

    for (ind_dom = 1 ; ind_dom < cs_glob_base_nbr ; ind_dom++)
      MPI_Send((void *)val, nbr_val, mpi_typ, ind_dom, CS_SUITE_MPI_TAG,
               cs_glob_base_mpi_comm);

  }

  else if (cs_glob_base_rang > 0) {

    MPI_Recv((void *)val, nbr_val, mpi_typ, 0, CS_SUITE_MPI_TAG,
             cs_glob_base_mpi_comm, &status);

  }

}

#endif /* #if defined(_CS_HAVE_MPI) */


/*----------------------------------------------------------------------------
 *  Fonction qui analyse le contenu d'un fichier suite en mode texte
 *  On démare en haut de deuxième ligne, ayant lu l'en-tête, et on
 *  indique en retour si la suite continue sur un autre fichier
 *----------------------------------------------------------------------------*/

static cs_bool_t cs_loc_suite_analyse_txt
(
 cs_suite_t  *const  suite                /* --> Structure suite              */
)
{
  char buffer [CS_SUITE_LNG_BUF_ASCII + 1];
  char sub    [CS_SUITE_LNG_BUF_ASCII + 1];
  char str_typ[CS_SUITE_LNG_BUF_ASCII + 1];
  char *str;

  cs_int_t  ind, ind_rec, ind_rub, ind_fin, lng, etape;
  cs_int_t  nbr_cel, nbr_fac, nbr_fbr, nbr_som;

  int nbr_ent, nbr_val_ent;

  cs_int_t  num_ligne = 1;
  cs_int_t  ind_fic   = suite->nbr_fic - 1;

  cs_bool_t  err_fmt = CS_FALSE;
  cs_bool_t  fin_fic = CS_FALSE;

  /* Initialisations */

  for (ind = CS_SUITE_SUPPORT_CEL ; ind <= CS_SUITE_SUPPORT_SOM ; ind++) {

    str = bft_file_gets_try(buffer,
                            CS_SUITE_LNG_BUF_ASCII,
                            suite->fic[ind_fic],
                            &num_ligne);

    if (str != NULL) {

      /* Décodage en-tête globale */

      if (sscanf(buffer, "[%s = %d]", sub, &nbr_ent) == 2) {

        if (strncmp("nbr_cel", sub, strlen("nbr_cel")) == 0)
          nbr_cel = nbr_ent;
        else if (strncmp ("nbr_fac_int", sub, strlen("nbr_fac_int")) == 0)
          nbr_fac = nbr_ent;
        else if (strncmp("nbr_fac_brd", sub, strlen("nbr_fac_brd")) == 0)
          nbr_fbr = nbr_ent;
        else if (strncmp("nbr_som", sub, strlen("nbr_som")) == 0)
          nbr_som = nbr_ent;
        else {
          err_fmt = CS_TRUE;
          break;
        }

      }

      else {
        err_fmt = CS_TRUE;
        break;
      }

    }
    else {               /* if (str != NULL) */

      err_fmt = CS_TRUE;
      /* break; */

    }

  }

  if (err_fmt == CS_TRUE) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The restart file <%s> number <%d> does not comply\n"),
               suite->nom, suite->nbr_fic);

    return CS_FALSE;
  }

  /* Dimensions du support */

  if (ind_fic == 0) {
    suite->nbr_cel = nbr_cel;
    suite->nbr_fac = nbr_fac;
    suite->nbr_fbr = nbr_fbr;
    suite->nbr_som = nbr_som;
  }
  else if (ind_fic > 0 && (   nbr_cel != suite->nbr_cel
                           || nbr_fac != suite->nbr_fac
                           || nbr_fbr != suite->nbr_fbr
                           || nbr_som != suite->nbr_som)) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The dimensions of the mesh associated to the restart file <%s>\n"
                 "number <%d> does not correspond to the file number <1>\n"),
               suite->nom, suite->nbr_fic);

    return CS_FALSE;
  }

  /* Tableaux et champs contenus dans le fichier */

  while (fin_fic == CS_FALSE) {

    str = bft_file_gets_try(buffer,
                            CS_SUITE_LNG_BUF_ASCII,
                            suite->fic[ind_fic],
                            &num_ligne);

    if (str == NULL) {
      fin_fic = CS_TRUE;
      break;
    }

    /* Traitement d'une rubrique */

    else if (buffer[0] == '[') {

      lng = strlen(buffer);
      while (   buffer[lng - 1] == '\n' || buffer[lng - 1] == '\r'
                || buffer[lng - 1] == '\t' || buffer[lng - 1] == ' ') {
        buffer[lng - 1] = '\0';
        lng --;
      }
      strcpy (sub, buffer + 1);

      /* Indicateur de fin de fichier */

      if (strncmp(sub, cs_suite_nom_fin_fic,
                  strlen(cs_suite_nom_fin_fic)) == 0)
        return CS_FALSE;

      /* Indicateur de découpage de fichier */

      else if (strncmp(sub, cs_suite_nom_decoup_fic,
                       strlen(cs_suite_nom_decoup_fic)) == 0)
        return CS_TRUE;

      /* Numéro de partie de fichier (en cas de découpage) */

      else if (strncmp(sub, cs_suite_nom_partie_fic,
                       strlen(cs_suite_nom_partie_fic)) == 0) {
        int num_part;
        ind = strlen(cs_suite_nom_partie_fic);
        for (ind_fin = lng - 1 ;
             ind_fin >= ind && sub[ind_fin] != ']' ;
             ind_fin--);
        sub[ind_fin] = '\0';
        if (   sscanf(sub + ind, "%d", &num_part) !=1
            || num_part != suite->nbr_fic) {
          cs_base_warn(__FILE__, __LINE__);
          bft_printf(_("The restart file <%s_p%02d> does not correspond to\n"
                       "part <%d> of the original restart file\n"),
                     suite->nom, (int)suite->nbr_fic, (int)suite->nbr_fic);

          return CS_FALSE;
        }
        continue;
      }

      /* Traitement d'un enregistrement standard */

      if (suite->nbr_rec == suite->nbr_rec_max) {
        suite->nbr_rec_max *= 2;
        BFT_REALLOC(suite->tab_rec, suite->nbr_rec_max, cs_suite_rec_t);
      }

      ind_rec = suite->nbr_rec;

      BFT_MALLOC((suite->tab_rec[ind_rec]).nom, strlen(sub) + 1, char);
      strcpy((suite->tab_rec[ind_rec]).nom, sub);

      str = bft_file_gets_try(buffer,
                              CS_SUITE_LNG_BUF_ASCII,
                              suite->fic[ind_fic],
                              &num_ligne);

      if (str == NULL) {
        fin_fic = CS_TRUE;
        break;
      }

      /*
        Décodage de la deuxième ligne d'entête (avec vérification
        de la syntaxe) ; on part de la fin, qui doit être caractérisée
        par le caractère ']', et on progresse vers le début
      */

      ind_fin = strlen(buffer);
      for (ind = ind_fin - 1 ; ind > 0 && buffer[ind] != ']' ; ind--);
      buffer[ind] = '\0';
      ind_fin = ind - 1;

      for (etape = 0 ; etape < 3 ;etape++) {

        if (ind_fin > 0 && err_fmt == CS_FALSE) {
          for ( ; ind > 0 && buffer[ind] != '=' ; ind--);
          for (ind_rub = ind ;
               ind_rub > 0 && buffer[ind_rub] != ',' ;
               ind_rub--);
          ind++;
          for ( ; ind < ind_fin && buffer[ind] == ' ' ; ind++);

          switch(etape) {

          case 0:
            if (strncmp(", typ_val = ", buffer+ind_rub,
                        strlen(", typ_val = ")) == 0)
              strcpy(str_typ, buffer+ind);
            else
              err_fmt = CS_TRUE;
            break;

          case 1:
            if (strncmp(", nbr_val_ent = ", buffer+ind_rub,
                        strlen(", nbr_val_ent = ")) == 0) {
              if (sscanf(buffer+ind, "%d", &nbr_val_ent) != 1)
                err_fmt = CS_TRUE;
            }
            else
              err_fmt = CS_TRUE;
            break;

          case 2:
            if (strncmp(" support = ", buffer+ind_rub,
                        strlen(" support = ")) == 0)
              strcpy(sub, buffer+ind);
            else
              err_fmt = CS_TRUE;
            break;

          }

          ind = ind_rub;
          buffer[ind_rub] = '\0';
          ind_fin = ind_rub - 1;
        }
        else
          err_fmt = CS_TRUE;

      }

      /* Remplissage de la structure */

      for (ind = CS_SUITE_SUPPORT_SCAL ;
           ind <= CS_SUITE_SUPPORT_SOM ; ind++) {
        if (strcmp(sub, cs_suite_nom_support[ind]) == 0) {
          (suite->tab_rec[ind_rec]).ind_support = (cs_suite_support_t)ind;
            break;
        }
      }
      if (ind > CS_SUITE_SUPPORT_SOM)
        err_fmt = CS_TRUE;
      (suite->tab_rec[ind_rec]).nbr_val_ent = (cs_int_t) nbr_val_ent;
      for (ind = 0 ; ind < 3 ; ind++) {
        if (strcmp(str_typ, cs_suite_nom_typ_elt[ind]) == 0) {
          (suite->tab_rec[ind_rec]).typ_val = (cs_type_t)ind;
            break;
        }
      }
      if (ind == 3)
        err_fmt = CS_TRUE;

      /* Position dans le fichier */

      (suite->tab_rec[ind_rec]).ind_fic = ind_fic;
      (suite->tab_rec[ind_rec]).pos_fic = bft_file_tell(suite->fic[ind_fic]);

      /* Libération du nom si erreur de format */

      if (   (fin_fic == CS_TRUE || err_fmt == CS_TRUE)
          && (suite->tab_rec[ind_rec]).nom != NULL)
        BFT_FREE((suite->tab_rec[ind_rec]).nom);

      else {

        /* Incrémentation compteur */

        suite->nbr_rec += 1;

        /* Passage rapide à la suite */

        bft_file_seek
          (suite->fic[ind_fic],
           (  (suite->tab_rec[ind_rec]).pos_fic
            + cs_loc_suite_calc_avance
                (suite,
                 (cs_suite_support_t)((suite->tab_rec[ind_rec]).ind_support),
                 (suite->tab_rec[ind_rec]).nbr_val_ent,
                 (suite->tab_rec[ind_rec]).typ_val)),
           BFT_FILE_SEEK_SET);

      }

    }

  }

  return CS_FALSE;

}


/*----------------------------------------------------------------------------
 *  Fonction qui analyse le contenu d'un fichier suite en mode binaire
 *  On démarre à la position suivant l'en-tête, et on indique en retour
 *  si la suite continue sur un autre fichier
 *----------------------------------------------------------------------------*/

static cs_bool_t cs_loc_suite_analyse_bin
(
 cs_suite_t  *const  suite                /* --> Structure suite              */
)
{
  char       buf_nom[CS_SUITE_LNG_NOM + 1];
  char      *nom = buf_nom;
  cs_int_t   lng_nom_max = CS_SUITE_LNG_NOM;

  cs_int_t  buf[4], nbr_lus, ind_rec;

  cs_int_t  ind_fic   = suite->nbr_fic - 1;

  cs_bool_t  fic_suiv = CS_FALSE;
  cs_bool_t  fin_fic  = CS_FALSE;

  /* Initialisations */

  nbr_lus = bft_file_read_try(buf,
                              sizeof(cs_int_t),
                              4,
                              suite->fic[ind_fic]);

  if (nbr_lus == 4) {

    /* Dimensions du support */

    if (ind_fic == 0) {
      suite->nbr_cel = buf[0];
      suite->nbr_fac = buf[1];
      suite->nbr_fbr = buf[2];
      suite->nbr_som = buf[3];
    }

    else if (ind_fic > 0 && (   buf[0] != suite->nbr_cel
                             || buf[1] != suite->nbr_fac
                             || buf[2] != suite->nbr_fbr
                             || buf[3] != suite->nbr_som)) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The dimensions of the mesh associated to the restart file <%s>\n"
                   "number <%d> does not correspond to the file number <1>\n"),
                 suite->nom, suite->nbr_fic);
      return CS_FALSE;
    }

  }
  else {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The restart file <%s> number <%d> does not comply\n"),
               suite->nom, suite->nbr_fic);
    return CS_FALSE;
  }


  /* Tableaux et champs contenus dans le fichier */

  while (fin_fic == CS_FALSE) {

    nbr_lus = bft_file_read_try(buf,
                                sizeof(cs_int_t),
                                4,
                                suite->fic[ind_fic]);

    if (nbr_lus < 4) {
      fin_fic = CS_TRUE;
      break;
    }

    if (buf[0] > lng_nom_max) {
      lng_nom_max = buf[0];
      if (nom == buf_nom)
        BFT_MALLOC(nom, lng_nom_max, char);
      else
        BFT_REALLOC(nom, lng_nom_max, char);
    }

    bft_file_read(nom,
                  1,
                  buf[0],
                  suite->fic[ind_fic]);

    /* Traitement d'une rubrique */

    /* Indicateur de fin de fichier */

    if (strncmp(nom, cs_suite_nom_fin_fic,
                strlen(cs_suite_nom_fin_fic)) == 0) {
      fin_fic = CS_TRUE;
      break;
    }

    /* Indicateur de découpage de fichier */

    else if (strncmp(nom, cs_suite_nom_decoup_fic,
                     strlen(cs_suite_nom_decoup_fic)) == 0) {
      fin_fic  = CS_TRUE;
      fic_suiv = CS_TRUE;
      break;
    }

    /* Numéro de partie de fichier (en cas de découpage) */

    else if (strncmp(nom, cs_suite_nom_partie_fic,
                     strlen(cs_suite_nom_partie_fic)) == 0) {
      if (buf[1] != suite->nbr_fic) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(_("The restart file <%s_p%02d> does not correspond to\n"
                     "part <%d> of the original restart file\n"),
                   suite->nom, (int)suite->nbr_fic, (int)suite->nbr_fic);
        fin_fic = CS_TRUE;
        break;
      }
      continue;
    }

    /* Traitement d'un enregistrement standard */

    if (suite->nbr_rec == suite->nbr_rec_max) {
      suite->nbr_rec_max *= 2;
      BFT_REALLOC(suite->tab_rec, suite->nbr_rec_max, cs_suite_rec_t);
    }

    ind_rec = suite->nbr_rec;

    BFT_MALLOC((suite->tab_rec[ind_rec]).nom, strlen(nom) + 1, char);
    strcpy((suite->tab_rec[ind_rec]).nom, nom);

    (suite->tab_rec[ind_rec]).ind_support = buf[1];
    (suite->tab_rec[ind_rec]).nbr_val_ent = buf[2];
    (suite->tab_rec[ind_rec]).typ_val     = (cs_type_t)buf[3];

    /* Position dans le fichier */

    (suite->tab_rec[ind_rec]).ind_fic = ind_fic;
    (suite->tab_rec[ind_rec]).pos_fic = bft_file_tell(suite->fic[ind_fic]);

    /* Incrémentation compteur */

    suite->nbr_rec += 1;

    /* Passage rapide à la suite */

    bft_file_seek
      (suite->fic[ind_fic],
       (  (suite->tab_rec[ind_rec]).pos_fic
        + cs_loc_suite_calc_avance
            (suite,
             (cs_suite_support_t)((suite->tab_rec[ind_rec]).ind_support),
             (suite->tab_rec[ind_rec]).nbr_val_ent,
             (suite->tab_rec[ind_rec]).typ_val)),
       BFT_FILE_SEEK_SET);

  }

  if (nom != buf_nom)
    BFT_FREE(nom);

  return fic_suiv;

}


/*----------------------------------------------------------------------------
 *  Fonction qui prépare l'index généré lors de l'analyse du fichier à
 *  l'utilisation par les fonctions de lecture d'enregistrements
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_prepare_index
(
 cs_suite_t  *const  suite                /* --> Structure suite              */
)
{
  cs_int_t  ind_rec;

  /* Traitement pour le parallélisme */

#if defined(_CS_HAVE_MPI)

  if (cs_glob_base_rang >= 0) {

    cs_int_t   buf[5];
    cs_int_t   lng_nom;
    cs_int_t  *pos_nom = NULL, *buf_idx = NULL;
    char      *buf_nom = NULL;

    /* Distribution dimensions et index */

    buf[0] = suite->nbr_cel;
    buf[1] = suite->nbr_fac;
    buf[2] = suite->nbr_fbr;
    buf[3] = suite->nbr_som;
    buf[4] = suite->nbr_rec;

    MPI_Bcast((void *)buf, 5, CS_MPI_INT, 0,
              cs_glob_base_mpi_comm);

    if (cs_glob_base_rang > 0) {

      suite->nbr_cel = buf[0];
      suite->nbr_fac = buf[1];
      suite->nbr_fbr = buf[2];
      suite->nbr_som = buf[3];
      suite->nbr_rec = buf[4];

      suite->nbr_rec_max = suite->nbr_rec;
      BFT_MALLOC(suite->tab_rec, suite->nbr_rec, cs_suite_rec_t);

    }

    /* Noms des enregistrements */

    BFT_MALLOC(pos_nom, suite->nbr_rec + 1, cs_int_t);

    if (cs_glob_base_rang == 0) {
      pos_nom[0] = 0;
      for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++)
        pos_nom[ind_rec + 1] =   pos_nom[ind_rec]
                               + strlen((suite->tab_rec[ind_rec]).nom);
    }

    MPI_Bcast((void *)pos_nom, suite->nbr_rec + 1, CS_MPI_INT, 0,
              cs_glob_base_mpi_comm);

    BFT_MALLOC(buf_nom, pos_nom[suite->nbr_rec], char);

    if (cs_glob_base_rang == 0) {
      for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++) {
        lng_nom = pos_nom[ind_rec + 1] - pos_nom[ind_rec];
        strncpy(buf_nom + pos_nom[ind_rec],
                (suite->tab_rec[ind_rec]).nom, lng_nom);
      }
    }

    MPI_Bcast((void *)buf_nom, pos_nom[suite->nbr_rec], MPI_CHAR, 0,
              cs_glob_base_mpi_comm);

    if (cs_glob_base_rang > 0) {
      for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++) {
        lng_nom = pos_nom[ind_rec + 1] - pos_nom[ind_rec];
        BFT_MALLOC ((suite->tab_rec[ind_rec]).nom, lng_nom + 1, char);
        strncpy((suite->tab_rec[ind_rec]).nom,
                buf_nom + pos_nom[ind_rec], lng_nom);
        (suite->tab_rec[ind_rec]).nom[lng_nom] = '\0';
      }
    }

    BFT_FREE(buf_nom);
    BFT_FREE(pos_nom);

    /* Autres parties de l'index */

    BFT_MALLOC(buf_idx, suite->nbr_rec * 3, cs_int_t);

    if (cs_glob_base_rang == 0) {
      for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++) {
        buf_idx[ind_rec * 3]
          = (suite->tab_rec[ind_rec]).ind_support;
        buf_idx[ind_rec * 3 + 1]
          = (suite->tab_rec[ind_rec]).nbr_val_ent;
        buf_idx[ind_rec * 3 + 2]
          = (cs_int_t)(suite->tab_rec[ind_rec]).typ_val;
      }
    }

    MPI_Bcast((void *)buf_idx, suite->nbr_rec * 3, CS_MPI_INT, 0,
              cs_glob_base_mpi_comm);

    if (cs_glob_base_rang > 0) {
      for (ind_rec = 0 ; ind_rec < suite->nbr_rec ; ind_rec++) {
        (suite->tab_rec[ind_rec]).ind_support = buf_idx[ind_rec * 3];
        (suite->tab_rec[ind_rec]).nbr_val_ent = buf_idx[ind_rec * 3 + 1];
        (suite->tab_rec[ind_rec]).typ_val
          = (cs_type_t)(buf_idx[ind_rec * 3 + 2]);
        (suite->tab_rec[ind_rec]).ind_fic = -1;
        (suite->tab_rec[ind_rec]).pos_fic = -1;
      }
    }

    BFT_FREE(buf_idx);

  }

#endif /* defined(_CS_HAVE_MPI) */

}


/*----------------------------------------------------------------------------
 *  Conversion d'arguments de lecture/écriture de l'API Fortran vers l'API C
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_rub_f77_vers_C
(
 const cs_int_t             *const numsui,  /* --> Numéro du fichier suite    */
 const cs_int_t             *const itysup,  /* --> Code type de support       */
 const cs_int_t             *const irtype,  /* --> Entiers ou réels ?         */
       cs_suite_t          **const suite,   /* <-- Pointeur structure suite   */
       cs_suite_support_t   *const support, /* <-- Type de support            */
       cs_type_t            *const typ_val, /* <-- Entiers ou réels           */
       cs_int_t             *const ierror   /* <-- 0 = succès, < 0 = erreur   */
)
{
  cs_int_t indsui = *numsui - 1;


  *ierror = CS_SUITE_SUCCES;

  /* Pointeur de structure suite associé */

  if (   indsui < 0
      || indsui > cs_glob_suite_ptr_nbr
      || cs_glob_suite_ptr_tab[indsui] == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The restart file number <%d> cannot be closed\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *ierror = CS_SUITE_ERR_NUM_FIC;
    return;
  }

  else
    *suite = cs_glob_suite_ptr_tab[indsui];


  /* Support associé à la rubrique */

  switch (*itysup) {

  case 0:
    *support = CS_SUITE_SUPPORT_SCAL;
    break;

  case 1:
    *support = CS_SUITE_SUPPORT_CEL;
    break;

  case 2:
    *support = CS_SUITE_SUPPORT_FAC_INT;
    break;

  case 3:
    *support = CS_SUITE_SUPPORT_FAC_BRD;
    break;

  case 4:
    *support = CS_SUITE_SUPPORT_SOM;
    break;

  default:
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The type of location <%d> indicated for a restart file\n"
                 "section is invalid"), (int)(*itysup));
    *ierror = CS_SUITE_ERR_SUPPORT;
    return;

  }


  /* Type associé à la rubrique */

  switch (*irtype) {

  case 1:
    *typ_val = CS_TYPE_cs_int_t;
    break;

  case 2:
    *typ_val = CS_TYPE_cs_real_t;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("The type of value <%d> indicated for a restart file\n"
                "section is invalid"), (int)(*irtype));
    *ierror = CS_SUITE_ERR_TYPE_VAL;
    return;

  }

}


/*----------------------------------------------------------------------------
 *  Permutation des valeurs d'un tableau renuméroté en lecture
 *----------------------------------------------------------------------------*/

static void cs_loc_suite_permute_lec
(
 const cs_int_t            nbr_ent,       /* --> Nombre d'entités             */
 const cs_int_t     *const num_ent_ini,   /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
       cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
)
{
  cs_int_t  ind_ent, ind_loc;

  cs_int_t  ind = 0;

  /* Instructions */

  if (num_ent_ini == NULL)
    return;

  switch (typ_val) {

  case CS_TYPE_cs_int_t:
    {
      cs_int_t  *val_ord;
      cs_int_t  *val_cur = (cs_int_t *)tab_val;

      BFT_MALLOC(val_ord, nbr_ent * nbr_val_ent, cs_int_t);

      for (ind_ent = 0 ; ind_ent < nbr_ent ; ind_ent++) {
        for (ind_loc = 0 ; ind_loc < nbr_val_ent ; ind_loc++)
          val_ord[ind++]
            = val_cur[(num_ent_ini[ind_ent] - 1) * nbr_val_ent + ind_loc];
      }

      for (ind = 0 ; ind < nbr_ent * nbr_val_ent ; ind++)
        val_cur[ind] = val_ord[ind];

      BFT_FREE(val_ord);
    }
    break;

  case CS_TYPE_cs_real_t:
    {
      cs_real_t  *val_ord;
      cs_real_t  *val_cur = (cs_real_t *)tab_val;

      BFT_MALLOC (val_ord, nbr_ent * nbr_val_ent, cs_real_t);

      for (ind_ent = 0 ; ind_ent < nbr_ent ; ind_ent++) {
        for (ind_loc = 0 ; ind_loc < nbr_val_ent ; ind_loc++)
          val_ord[ind++]
            = val_cur[(num_ent_ini[ind_ent] - 1) * nbr_val_ent + ind_loc];
      }

      for (ind = 0 ; ind < nbr_ent * nbr_val_ent ; ind++)
        val_cur[ind] = val_ord[ind];

      BFT_FREE(val_ord);
    }
    break;

  default:
    assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);

  }

}


/*----------------------------------------------------------------------------
 *  Permutation des valeurs d'un tableau renuméroté en écriture
 *----------------------------------------------------------------------------*/

static cs_byte_t * cs_loc_suite_permute_ecr
(
 const cs_int_t            nbr_ent,       /* --> Nombre d'entités             */
 const cs_int_t     *const num_ent_ini,   /* --> Numéros globaux des entités  */
 const cs_int_t            nbr_val_ent,   /* --> Nombre de valeurs/entité     */
 const cs_type_t           typ_val,       /* --> Type de valeur               */
 const cs_byte_t    *const tab_val        /* --> Tableau des valeurs          */
)
{
  cs_int_t  ind_ent, ind_loc;

  cs_int_t  ind = 0;

  /* Instructions */

  if (num_ent_ini == NULL)
    return NULL;

  switch (typ_val) {

  case CS_TYPE_cs_int_t:
    {
      cs_int_t  *val_ord;
      const cs_int_t  *val_cur = (const cs_int_t *)tab_val;

      BFT_MALLOC(val_ord, nbr_ent * nbr_val_ent, cs_int_t);

      for (ind_ent = 0 ; ind_ent < nbr_ent ; ind_ent++) {
        for (ind_loc = 0 ; ind_loc < nbr_val_ent ; ind_loc++)
          val_ord[(num_ent_ini[ind_ent] - 1) * nbr_val_ent + ind_loc]
            = val_cur[ind++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  case CS_TYPE_cs_real_t:
    {
      cs_real_t  *val_ord;
      const cs_real_t  *val_cur = (const cs_real_t *)tab_val;

      BFT_MALLOC(val_ord, nbr_ent * nbr_val_ent, cs_real_t);

      for (ind_ent = 0 ; ind_ent < nbr_ent ; ind_ent++) {
        for (ind_loc = 0 ; ind_loc < nbr_val_ent ; ind_loc++)
          val_ord[(num_ent_ini[ind_ent] - 1) * nbr_val_ent + ind_loc]
            = val_cur[ind++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  default:
    assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);
    return NULL;

  }

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
