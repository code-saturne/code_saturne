/*============================================================================
 *  Modification d'un fichier NOPO (Simail)
 *============================================================================*/

/*
  This file is part of the code_saturne Preprocessor, element of the
  code_saturne CFD tool.

  Copyright (C) 1999-2022 EDF S.A., France

  contact: saturne-support@edf.fr

  The code_saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The code_saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the code_saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/*============================================================================
 *                  Définitions de paramètres et macros
 *============================================================================*/

/*  Types de dimension fixée. */

#if defined(__linux__) || defined(__linux) || defined(linux)
#include <stdint.h>

typedef int32_t         int_32_t;   /* Entier sur 4 octets */
typedef float           real_32_t;  /* Réel sur 4 octets */

#else

typedef int             int_32_t;   /* Entier sur 4 octets */
typedef float           real_32_t;  /* Réel sur 4 octets */

#endif


/* Indices premier et dernier groupe de 4 octets à ne pas permuter */

#define NOPO_CHAR_DEB     (1 + 14 + 1) + (1 + 1)
#define NOPO_CHAR_FIN     (1 + 14 + 1) + (1 + 30)

#define NOPO_TAILLE_E     56        /* Longueur entête (octets)   */
#define NOPO_TAILLE_0     (32+1)*4  /* Longueur NOP0 (octets)     */
#define NOPO_TAILLE_2     (27+1)*4  /* Longueur NOP2 (mots 4 o.)  */

#define NOPO_NBR_MAX_TAB_REF      26        /* Valeur max de NMAE + 1
                                                   (nb. sommets + nb. arêtes
                                                   +nb. faces : 8 + 12 + 6)   */

#define NOPO_NUL                   0        /* Inexistant                 */
#define NOPO_NODE                  1        /* Noeud                      */
#define NOPO_SEGMENT               2        /* Segment                    */
#define NOPO_TRIANGLE              3        /* Triangle                   */
#define NOPO_QUADRANGLE            4        /* Quadrangle                 */
#define NOPO_TETRAHEDRON           5        /* Tétraèdre                  */
#define NOPO_PENTAHEDRON           6        /* Prisme                     */
#define NOPO_HEXAHEDRON            7        /* Hexaèdre                   */
#define NOPO_SUPER_ELEMENT         8        /* Super-élément (inutilisé)  */

typedef struct {

  int_32_t  dim_e;      /* Dimension de l'espace */
  int_32_t  ndsr;       /* numéro de référence maximal */
  int_32_t  ndsd;       /* numéro de sous-domaine maximal */
  int_32_t  ncopnp;     /* 1 si sommets = noeuds, 0 sinon */
  int_32_t  np;         /* nombre de points */
  int_32_t  ne;         /* nombre d'éléments */
  int_32_t  ne_typ[9];  /* nombre d'éléments par type */
  int_32_t  nepo;       /* nombre d'éléments point */
  int_32_t  nesg;       /* nombre de segments */
  int_32_t  ntri;       /* nombre de triangles */
  int_32_t  nqua;       /* nombre de quadrangles */
  int_32_t  ntet;       /* nombre de tétraèdres */
  int_32_t  npen;       /* nombre de pentaèdres */
  int_32_t  nhex;       /* nombre d'hexaèdres */
  int_32_t  noe;        /* nombre de noeuds */

  int_32_t  lnop5;      /* Taille du tableau nop5 */
  real_32_t  *coord;    /* Coordonnées des sommets */
  int_32_t   *nop5;     /* Connectivités et références des éléments */

} nopo_maillage_t;

/*============================================================================
 *                  Définition de structures locales
 *============================================================================*/


/* Définition des éléments */
/*=========================*/

typedef struct {
  int_32_t elt_typ;
  int_32_t som[4];
} nopo_sous_elt_t ;

typedef struct {
  int_32_t         nopo_typ;    /* Type NOPO de l'élément  */
  int_32_t         num_som[8];  /* Numéros de sommets      */
  int_32_t         nbr_som;     /* Nombre de sommets */
  int_32_t         nbr_are_vol; /* Nombre arêtes si volume */
  int_32_t         nbr_sselt;   /* Nombre de sous-éléments */
  nopo_sous_elt_t  sous_elt[6] ;
} nopo_elt_t ;


static const nopo_elt_t  nopo_elt_liste_c[8] = {

  {                        /* 1 */
    NOPO_NUL,
    { 0 },
    0,
    0,
    0,
    {
      {0,{0}}
    }
  },
  {                        /* 1 */
    NOPO_NODE,
    { 1 },
    1,
    0,
    0,
    {
      {0,{0}}
    }
  },
  {                        /* 2 */
    NOPO_SEGMENT,
    { 1, 2 },
    2,
    0,
    2,
    {                                          /*    1       2            */
      {NOPO_NODE, { 1 }},                      /*    x-------x            */
      {NOPO_NODE, { 2 }}
    }
  },
  {                        /* 3 */
    NOPO_TRIANGLE,
    { 1, 2, 3 },
    3,
    0,
    3,
    {                                          /*        x 3              */
      {NOPO_SEGMENT, { 1 , 2 }} ,              /*       / \               */
      {NOPO_SEGMENT, { 2 , 3 }} ,              /*      /   \              */
      {NOPO_SEGMENT, { 3 , 1 }}                /*     /     \             */
    }                                          /*  1 x-------x 2          */
  },
  {                        /* 4 */
    NOPO_QUADRANGLE,
    { 1, 2, 3, 4 },
    4,
    0,
    4,
    {                                          /*  4 x-------x 3          */
      {NOPO_SEGMENT, { 1 , 2 }} ,              /*    |       |            */
      {NOPO_SEGMENT, { 2 , 3 }} ,              /*    |       |            */
      {NOPO_SEGMENT, { 3 , 4 }} ,              /*    |       |            */
      {NOPO_SEGMENT, { 4 , 1 }}                /*  1 x-------x 2          */
    }                                          /*                         */
  },
  {                        /* 5 */
    NOPO_TETRAHEDRON,
    { 1, 2, 3, 4 },
    4,                                         /*                         */
    6,                                         /*        x 4              */
    4,                                         /*       /|\               */
    {                                          /*      / | \              */
      {NOPO_TRIANGLE, { 1 , 3 , 2 }},          /*     /  |  \             */
      {NOPO_TRIANGLE, { 1 , 4 , 3 }},          /*  1 x- -|- -x 3          */
      {NOPO_TRIANGLE, { 1 , 2 , 4 }},          /*     \  |  /             */
      {NOPO_TRIANGLE, { 2 , 3 , 4 }}           /*      \ | /              */
    }                                          /*       \|/               */
  },                                           /*        x 2              */
  {                       /*  6 */
    NOPO_PENTAHEDRON,
    { 1 , 2 , 3 , 4 , 5 , 6 },
    6,                                         /*                         */
    9,                                         /*  4 x-------x 6          */
    5,                                         /*    |\     /|            */
    {                                          /*    | \   / |            */
      {NOPO_TRIANGLE,    { 1 , 3 , 2 }    },   /*  1 x- \-/ -x 3          */
      {NOPO_QUADRANGLE,  { 1 , 4 , 6 , 3 }},   /*     \ 5x  /             */
      {NOPO_QUADRANGLE,  { 1 , 2 , 5 , 4 }},   /*      \ | /              */
      {NOPO_TRIANGLE,    { 4 , 5 , 6 }    },   /*       \|/               */
      {NOPO_QUADRANGLE,  { 2 , 3 , 6 , 5 }}    /*        x 2              */
    }
  },
  {                       /*  7 */
    NOPO_HEXAHEDRON,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
    8,
    12,
    6,                                         /*     8 x-------x 7       */
    {                                          /*      /|      /|         */
      {NOPO_QUADRANGLE, { 1 , 4 , 3 , 2 }},    /*     / |     / |         */
      {NOPO_QUADRANGLE, { 1 , 5 , 8 , 4 }},    /*  5 x-------x6 |         */
      {NOPO_QUADRANGLE, { 1 , 2 , 6 , 5 }},    /*    | 4x----|--x 3       */
      {NOPO_QUADRANGLE, { 5 , 6 , 7 , 8 }},    /*    | /     | /          */
      {NOPO_QUADRANGLE, { 2 , 3 , 7 , 6 }},    /*    |/      |/           */
      {NOPO_QUADRANGLE, { 3 , 4 , 8 , 7 }}     /*  1 x-------x 2          */
    }
  }
} ;


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Modification des coordonnées d'un sommet ou suppression du sommet
 *  (les éléments s'appuyant sur les sommets supprimés seront supprimés.
 *  On renvoie 0 en cas normal, et -1 si le sommet est supprimé.
 *----------------------------------------------------------------------------*/

static int nopo_modif_som
(
 real_32_t coord[3]
)
{

  /* Exemple : translation [1,1,-1] ; attention aux symétries,
     qui nécessitent une modification des connectivités
     des cellules correspondantes */

  coord[0] = coord[0] + 1.0;
  coord[1] = coord[2] + 1.0;
  coord[2] = coord[0] - 1.0;

  /* Suppression du sommet en dessous d'un demi-sphère de
     centre [0,0,-10] et de rayon 5 */

  if (coord[2] < 10.0) {
    if (  coord[0]*coord[0]
        + coord[1]*coord[1]
        + (coord[2]-10.0)*(coord[2]-10.0) > 5.0)
      return -1;
  }

  return 0;

}


/*----------------------------------------------------------------------------
 * Modification de la référence d'une face (ou arête en 2D) en fonction
 * des coordonnées du centre de la face et de sa référence précédente.
 * La nouvelle référence est renvoyée.
 *----------------------------------------------------------------------------*/

static int_32_t nopo_modif_ref
(
 const real_32_t coord[3],
 const int_32_t  ref_old
)
{

  /* Exemple : pas de modification */

  return ref_old;

}


/*----------------------------------------------------------------------------
 *  Conversion big-endian/little-endian d'un enregistrement d'un fichier
 *   NOPO (Format INRIA utilisé par Simail)
 *----------------------------------------------------------------------------*/

static void nopo_permut_4
(
 char *rec,
 const int nbr
)
{

  int i;
  char c0, c1, c2, c3;
  char *p = rec;

  for (i = 0 ; i < nbr ; i++) {

    c0 = p[0];
    c1 = p[1];
    c2 = p[2];
    c3 = p[3];

    p[0] = c3;
    p[1] = c2;
    p[2] = c1;
    p[3] = c0;

    p += 4;

  }

}


/*----------------------------------------------------------------------------
 * Message d'erreur pour fichier NOPO
 *----------------------------------------------------------------------------*/

static void nopo_err_fic
(
 FILE *fic
)
{
  int  ind_err ;

  /* Vérification erreur éventuelle */

  if (fic == NULL)
    ind_err = errno;
  else
    ind_err = ferror(fic);

  if (ind_err != 0) {
    fprintf(stderr,
            "Erreur sur le fichier :\n%s\n",
            strerror(errno));
    exit(EXIT_FAILURE);
  }
  else if (feof(fic) != 0) {
    fprintf(stderr,
            "Fin de fichier prématurée\n");
    exit(EXIT_FAILURE);
  }

}
                  

/*----------------------------------------------------------------------------
 *  Lecture d'un enregistrement d'un fichier NOPO
 *----------------------------------------------------------------------------*/

static void nopo_lit_rec
(
 FILE *fic,
 const int swap_end_lec,
 const size_t taille,
 char **const rec
)
{
  int_32_t nbr_cmp;

  /* Marqueur début Fortran */

  if (fread((char *)(&nbr_cmp), 4, 1, fic) != 1)
     nopo_err_fic(fic);

  if (swap_end_lec == 1)
    nopo_permut_4((char *)(&nbr_cmp), 1);

  if (nbr_cmp == taille) {

    /* allocation mémoire */

    *rec = malloc(taille);
    if (*rec == NULL) {
      fprintf(stderr, "Impossible d'allouer %d octets\n", taille);
      exit(EXIT_FAILURE);
    }

    /* Contenu de l'enregistrement */

    if (fread((char *)(*rec), 1, taille, fic) != taille)
      nopo_err_fic(fic);

    if (swap_end_lec == 1)
      nopo_permut_4((char *)(*rec), taille/4);

    /* Marqueur fin Fortran */

    if (fread((char *)(&nbr_cmp), 4, 1, fic) != 1)
      nopo_err_fic(fic);

    if (swap_end_lec == 1)
      nopo_permut_4((char *)(&nbr_cmp), 1);

  }

  if (nbr_cmp != taille) {
    fprintf(stderr, "Erreur à la lecture d'un enregistrement : "
            "On attendait une taille\nde %d octets, il semble y "
            "en avoir %d\n", (int)taille, (int)nbr_cmp);
    exit(EXIT_FAILURE);
  }

}
                  

/*----------------------------------------------------------------------------
 *  Écriture d'un enregistrement d'un fichier NOPO
 *----------------------------------------------------------------------------*/

static void nopo_ecr_rec
(
 FILE *fic,
 const size_t taille,
 char *const rec
)
{
  int_32_t nbr_cmp = taille ;

  /* Marqueur début Fortran */

  if (fwrite((char *)(&nbr_cmp), 4, 1, fic) != 1)
    nopo_err_fic(fic);

  /* Contenu de l'enregistrement */

  if (fwrite((char *)rec, 1, taille, fic) != taille)
    nopo_err_fic(fic);

  /* Marqueur fin Fortran */

  if (fwrite((char *)(&nbr_cmp), 4, 1, fic) != 1)
    nopo_err_fic(fic);
}
                  

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier NOPO (Format INRIA utilisé par Simail)
 *   et affectation des donnees dans la structure de maillage
 *----------------------------------------------------------------------------*/

static nopo_maillage_t nopo_lit_maillage
(                                       /* <-- Renvoie un pointeur sur */
                                        /*     une structure de maillage */
 const char *const nom_fic_maillage     /* --> Nom du fichier a lire */
)
{
  FILE        *fic_maillage;         /* Descripteur du fichier */
  int_32_t    dim_e;                 /* Dimension spatiale */

  int          swap_end_lec = 0;
  int_32_t     ind;
  int_32_t     ind_test;

  int_32_t     lnop1; /* nombre de mots dans NOP1 */
  int_32_t     lnop3; /* nombre de mots dans NOP3 */
  int_32_t     lnop4; /* nombre de mots dans NOP4 */

  int_32_t  *nopo_e;

  int_32_t  *nop0 = NULL;
  int_32_t  *nop1 = NULL;
  int_32_t  *nop3 = NULL;
  int_32_t  *nop2 = NULL;
  real_32_t *nop4 = NULL;

  int_32_t  ntacoo; /* type de système de coordonnées :
                       1 : cartésien ; 2 : cyl ; 3 : sphérique */

  nopo_maillage_t m; /* Structure maillage associciée */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Affichage du titre */

  printf("\n\n"
         "Lecture du fichier de maillage au format NOPO (Simail)\n"
         "------------------------------\n") ;

  printf("  Fichier de maillage : %s\n\n\n", nom_fic_maillage) ;


  /* Ouverture du fichier NOPO en lecture */

  fic_maillage = fopen(nom_fic_maillage, "rb");

  if (fic_maillage == NULL)
     nopo_err_fic(fic_maillage);

  /* Test si le fichier est au format natif ou non */
  
  if (fread((char *)(&ind_test), 4, 1, fic_maillage) != 1)
     nopo_err_fic(fic_maillage);

  if (ind_test != NOPO_TAILLE_E) {

    nopo_permut_4((char *)(&ind_test), 1);

    if (ind_test == NOPO_TAILLE_E)
      swap_end_lec = 1;
    else {
      fprintf(stderr,
              "Erreur de format du fichier \"%s\" :\n"
              "Ce fichier ne semble pas être au format NOPO\n"
              "(la longueur de l'entête est différente de %d).",
              nom_fic_maillage, NOPO_TAILLE_E);
      exit(EXIT_FAILURE);
    }

  }

  /* Repositionnement en début de fichier */

  if (fseek(fic_maillage, 0, SEEK_SET) < 0)
     nopo_err_fic(fic_maillage);

  /*  On lit l'entête du fichier donnant les dimensions des tableaux */
  /*-----------------------------------------------------------------*/

  /*
    Attention, la première valeur d'un tableau NOP* donne la
    dimension du tableau, et le tableau est donc de taille n+1 ;
    les valeurs utiles étant aux positions 1 à n : on utilise donc
    une indexation entre 1 et n plutôt qu'entre 0 et n-1 pour
    y accéder.
  */

  nopo_lit_rec(fic_maillage, swap_end_lec, NOPO_TAILLE_E, (char**)(&nopo_e));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (ind = 0 ;
       ind < (int_32_t)(NOPO_TAILLE_E / sizeof(int_32_t)) ;
       ind++)
#endif

  /*
    D'après la documentation, la 5ème position de l'entête est réservée et
    vaut 0. En pratique, on constate qu'elle correspond à la dimension du
    tableau NOP3, qui existe lorsque l'on a des super-éléments ou des
    descriptions (contrairement toujours à la documentation Simail, qui
    ne semble pas à jour car elle indique que ce tableau est inutilisé).
  */

  lnop1 = nopo_e[3] ;
  lnop3 = nopo_e[5] ;
  lnop4 = nopo_e[6] ;
  m.lnop5 = nopo_e[7] ;

  if (   nopo_e[ 1] !=  6 || nopo_e[ 2] != 32 || nopo_e[ 4] != 27
      || nopo_e[ 8] !=  1 || nopo_e[ 9] !=  1 || nopo_e[10] !=  1
      || nopo_e[11] !=  1 || nopo_e[12] !=  2 || nopo_e[13] !=  1) {
    fprintf(stderr,
            "Erreur de format du fichier \"%s\" :\n"
            "Les valeurs réservées du descripteur ne correspondent\n"
            "pas aux valeurs attendues :\n"
            "ce fichier n'est probablement pas un fichier NOPO/Simail",
            nom_fic_maillage);
    exit(EXIT_FAILURE);
  }

  free(nopo_e);

  /* Lecture du tableau NOP0 */

  nopo_lit_rec(fic_maillage, swap_end_lec, NOPO_TAILLE_0, (char **)(&nop0));

  /*
    On rappelle que la première valeur d'un tableau NOP* donne la
    dimension du tableau, et le tableau est donc de taille n+1 ;
    On utilise donc une indexation entre 1 et n plutôt qu'entre 0 et n-1
    pour accéder aux autres valeurs.
  */

  /*
    NOP0 contient des chaînes de caractères codées sur des entiers,
    qu'il convient de remettre dans l'ordre si l'on a permuté les
    octets.
  */

  if (swap_end_lec == 1)
    nopo_permut_4((char *)(nop0 + 1), 29);

  /*
    Suppression blancs en fin de titre (80 caractères en général rarement
     utilisés) pour affichage sur ligne plus courte
  */
  for (ind = 80 ; ind > 0 && *((char *)(nop0) + ind) == ' ' ; ind--) {
    *((char *)(nop0) + ind) = '\0' ;
  }

  printf("  Titre     : %.80s\n", (char *)(nop0 + 1));
  printf("  Date      : %2.2s/%2.2s/%4.4s\n",
         (char *)(nop0 + 1 + 20),
         ((char *)(nop0 + 1 + 20)) + 2, ((char *)(nop0 + 1 + 20)) + 4);
  printf("  Créateur  : %24.24s\n", (char *)(nop0 + 1 + 22));

  if (strncmp("NOPO", (char *)(nop0 + 1 + 28), 4) != 0) {
    fprintf(stderr,
            "Erreur de format du fichier\"%s\" :\n"
            "Les caractères 'NOPO' n'apparaissent pas dans l'entête ;\n"
            "ce fichier n'est probablement pas un fichier NOPO/Simail",
            nom_fic_maillage);
    exit(EXIT_FAILURE);
  }


#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("nop0[   30] NIVEAU = %d\n",  nop0[30]) ;
  printf("nop0[   31] ETAT   = %d\n",  nop0[31]) ;
  printf("nop0[   32] NTACM  = %d\n",  nop0[32]) ;
#endif

  free(nop0) ;

  /* Lecture du tableau NOP1 */
  /* ----------------------- */

  /*
   * Ce tableau contient des tableaux auxiliaires, en général inutiles.
   */

  if (lnop1 != 0) {
    nopo_lit_rec(fic_maillage, swap_end_lec, (lnop1+1)*4, (char **)(&nop1));
    free(nop1);
  }
  nop1 = NULL;


  /* Lecture du tableau NOP2 */
  /* ----------------------- */

  nopo_lit_rec(fic_maillage, swap_end_lec, NOPO_TAILLE_2, (char **)(&nop2));

#if 1 && defined(DEBUG) && !defined(NDEBUG)
  for (ind = 0 ; ind < 28 ; ind++)
    printf("nop2[%2d] : % d\n", ind, nop2[ind]);
#endif

  /* Dimension du maillage */

  if (nop2[1] == 2)
    dim_e = 2;
  else
    dim_e = 3;

  printf("  Dimension : %d\n\n", (int) (nop2[1])) ;

  assert(nop2[1] == 2 || nop2[1] == 3) ;

  /* Autres dimensions et paramètres */

  m.ndsr   = (int_32_t) nop2[2] ;  /* numéro de référence maximal */
  m.ndsd   = (int_32_t) nop2[3] ;  /* numéro de sous-domaine maximal */
  m.ncopnp = (int_32_t) nop2[4] ;  /* 1 si sommets = noeuds, 0 sinon */
  m.ne     = (int_32_t) nop2[5] ;  /* nombre d'éléments */
  m.nepo   = (int_32_t) nop2[6] ;  /* nombre d'éléments point */
  m.nesg   = (int_32_t) nop2[7] ;  /* nombre de segments */
  m.ntri   = (int_32_t) nop2[8] ;  /* nombre de triangles */
  m.nqua   = (int_32_t) nop2[9] ;  /* nombre de quadrangles */
  m.ntet   = (int_32_t) nop2[10] ; /* nombre de tétraèdres */
  m.npen   = (int_32_t) nop2[11] ; /* nombre de pentaèdres */
  m.nhex   = (int_32_t) nop2[12] ; /* nombre d'hexaèdres */
  m.noe    = (int_32_t) nop2[15] ; /* nombre de noeuds */

  /*
   * Valeurs non utilisées ici :
   *
   * Remarque : on a d'après la documentation le nombre de mots dans NOP5
   *            en 27ème position de NOP2, et en 8ème position du tableau
   *            entête. Sur certains cas test (issus de Simail), la valeur
   *            fournie dans NOP2 n'est pas cohérente ; on l'ignore donc, et
   *            on conserve la valeur lue initialement. De la même manière,
   *            la présence de NOP3 dépend de LNOP3 (5ème valeur de l'entête)
   *            et non de NBEGM.
   *
   * nsup   = (int_32_t) nop2[13] ;  nombre de super éléments
   * nef    = (int_32_t) nop2[14] ;  nombre d'éléments de bord
   * n1     = (int_32_t) nop2[16] ;  nb. noeuds intérieurs segment ou arête
   * iset   = (int_32_t) nop2[17] ;  nb. noeuds intérieurs triangle ou face
   * iseq   = (int_32_t) nop2[18] ;  nb. noeuds intérieurs quadrangle ou face
   * isete  = (int_32_t) nop2[19] ;  nb. noeuds intérieurs tetraèdre
   * isepe  = (int_32_t) nop2[20] ;  nb. noeuds intérieurs pentaèdre
   * isehe  = (int_32_t) nop2[21] ;  nb. noeuds intérieurs hexaèdre
   * ntycoo = (int_32_t) nop2[23] ;  type de valeurs de coordonnées (2 ici)
   * lpgdn  = (int_32_t) nop2[24]    plus grande diff. num. noeuds même élé.
   * nbegm  = (int_32_t) nop2[25] ;  nombre de super-éléments dans NOP3
   */

  m.np     = (int_32_t) nop2[22] ; /* nombre de points */

  ntacoo = (int_32_t) nop2[27] ; /* type de système de coordonnées :
                                     1 : cartésien ; 2 : cyl ; 3 : sphérique */

  free(nop2) ;

  printf("  Données initiales : %10d points\n"
         "                      %10d noeuds\n"
         "                      %10d éléments points\n"
         "                      %10d segments\n"
         "                      %10d triangles\n"
         "                      %10d quadrangles\n",
         m.np, m.noe, m.nepo, m.nesg, m.ntri, m.nqua);

  if (dim_e == 2)
    printf("\n") ;
  else
    printf("                      %10d tétraèdres\n"
           "                      %10d pentaèdres\n"
           "                      %10d hexaèdres\n\n", 
           m.ntet, m.npen, m.nhex);

  if (ntacoo != 1) {
    fprintf(stderr,
            "Erreur à la lecture du fichier NOPO :\n\"%s\" ;\n"
            "Le système de coordonnés est cylindrique ou sphérique, ce\n"
            "qui n'est pas prévu dans la version actuelle",
            nom_fic_maillage);
    exit(EXIT_FAILURE);
  }

  /* "Blindages" (pas forcément nécessaires, mais on préfère être prudent) */

  assert(m.np <= m.noe);
  if (m.ncopnp == 1 && m.np == 0)
    m.np = m.noe;
  assert(m.np > 0);

  /* Lecture du tableau NOP3 */
  /* ----------------------- */

  /* Ce tableau contient des tableaux relatifs aux super-éléments */

  if (lnop3 != 0) {
    nopo_lit_rec(fic_maillage, swap_end_lec, (lnop3+1)*4, (char **)(&nop3));
    free(nop3);
  }
  nop3 = NULL;

  /* Lecture du tableau NOP4 */
  /* ----------------------- */

  nopo_lit_rec(fic_maillage, swap_end_lec, (lnop4+1)*4, (char **)(&nop4));

  /*
    Ce tableau contient des réels codés sur 4 octets ; de plus, il peut
    ne pas contenir de coordonnée "z" en 2D. Comme on ne s'intéresse
    qu'aux points et non aux noeuds intérieurs, on ne conserve que
    cette partie du tableau (les noeuds viennent après).
  */

  if (sizeof(real_32_t) != 4) {
    fprintf(stderr,
            "Ce sous-programme a été compilé avec un type \"real_32_t\"\n"
            "de taille différente de 4. (Erreur de portage bloquante pour\n"
            "la lecture et écriture des coordonnées simple précision NOPO).");
    exit(EXIT_FAILURE);
  }

  m.coord = malloc(m.np * 3 * sizeof(real_32_t));
  if (m.coord == NULL) {
    fprintf(stderr,
            "Impossible d'allouer %d octets\n",
            m.np * 3 * sizeof(real_32_t));
    exit(EXIT_FAILURE);
  }
  if (dim_e == 3) {
    for (ind = 0 ; ind < m.np * 3 ; ind++)
      m.coord[ind] = nop4[ind + 1] ;
  }
  else  if (dim_e == 2) {
    for (ind = 0 ; ind < m.np ; ind++) {
      m.coord[ind*3    ] = nop4[ind*3 + 1] ;
      m.coord[ind*3 + 1] = nop4[ind*3 + 2] ;
      m.coord[ind*3 + 2] = 0.0 ;
    }
  }
  free(nop4) ;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Coordonnées\n") ;
  for (ind = 0 ; ind < m.np ; ind++)
    printf("%d : % 10.5e % 10.5e % 10.5e\n",
           ind + 1, m.coord[ind*3    ],
           m.coord[ind*3 + 1], m.coord[ind*3 + 2]) ;
#endif

  /* Lecture du tableau NOP5 */
  /* ----------------------- */

  nopo_lit_rec(fic_maillage, swap_end_lec, (m.lnop5+1)*4, (char **)(&(m.nop5)));

  /* Fermeture du fichier de lecture du maillage */

  fclose(fic_maillage);

  /* renvoi du maillage */

  return m;

}


/*----------------------------------------------------------------------------
 * Traitement des modifications des sommets ; renvoie un tableau
 * de renumérotation des sommets (numérotation 1 à n) de dimension
 * maillage->
 *----------------------------------------------------------------------------*/

static int_32_t * nopo_traite_sommets
(
 nopo_maillage_t *m
)
{

  int_32_t  i, j, ret, taille;
  int_32_t  *renum_som;
  real_32_t coo_tmp[3];

  taille = m->np * 3 * sizeof(real_32_t);
  renum_som = malloc(taille);
  if (renum_som == NULL) {
    fprintf(stderr, "Impossible d'allouer %d octets\n", taille);
    exit(EXIT_FAILURE);
  }
  
  for (i = 0, j = 0 ; i < m->np ; i++) {

    coo_tmp[0] = m->coord[i*3];
    coo_tmp[1] = m->coord[i*3 + 1];
    coo_tmp[2] = m->coord[i*3 + 2];

    ret = nopo_modif_som(coo_tmp);

    if (ret != -1) {
      m->coord[j*3    ] = coo_tmp[0];
      m->coord[j*3 + 1] = coo_tmp[1];
      m->coord[j*3 + 2] = coo_tmp[2];
      j += 1;
      renum_som[i] = j;
    }
    else
      renum_som[i] = 0;

  }

  m->np = j;

  return renum_som;

}


/*----------------------------------------------------------------------------
 * Traitement des éléments (suppressions éventuelles, modification
 * des références des faces).
 *----------------------------------------------------------------------------*/

void nopo_traite_elements
(
 int_32_t *renum_som,
 nopo_maillage_t *m
)
{

  int_32_t    nbr_som ;
  int_32_t    nbr_sselt ;

  int_32_t    ient ;
  int_32_t    i, j;
  int_32_t    iel ;
  int_32_t    iloc ;
  int_32_t    ipos ;
  int_32_t    isom ;
  int_32_t    isselt ;
  int_32_t    issent ;

  int_32_t    tsselt ;

  int_32_t    ncge ;          /* Type d'élément */
  int_32_t    nmae ;
  int_32_t    ndsde ;
  int_32_t    nno ;
  int_32_t    npo ;
  int_32_t    ining ;
  int_32_t    iref ;
  int_32_t    tab_ref[NOPO_NBR_MAX_TAB_REF] ;

  int_32_t    sommets[8];
  int_32_t    references[6];
  int_32_t    nvol, lnop5, taille, diment;
  int_32_t    *nop5;

  real_32_t   cdg[3];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nvol = m->ntet + m->npen + m->nhex;

  /*
    Par type d'élément : surdimensionnement en prévoyant
    type + nombre de points (pour connectivité) + 1 (pour ining)
    + nombre de sous-éléments (pour références) + 1 (pour sous-domaine);
    On supprimera les sous-sous-éléments (i.e. les arêtes d'un maillage
    volumique) et les références correspondantes
  */

  lnop5 = 0;
  if (nvol > 0) {
    lnop5 += m->ntri * (1 + 1 + 3 + 0 + 1);
    lnop5 += m->nqua * (1 + 1 + 4 + 0 + 1);
    lnop5 += m->ntet * (1 + 1 + 4 + 4 + 1);
    lnop5 += m->npen * (1 + 1 + 6 + 5 + 1);
    lnop5 += m->nhex * (1 + 1 + 8 + 6 + 1);
  }
  else {
    lnop5 += m->nesg * (1 + 1 + 2 + 0 +1);
    lnop5 += m->ntri * (1 + 1 + 3 + 3 +1);
    lnop5 += m->nqua * (1 + 1 + 4 + 4 +1);
  }

  taille = (lnop5 + 1) * sizeof(int_32_t);
  nop5 = malloc(taille);
  if (nop5 == NULL) {
    fprintf(stderr, "Impossible d'allouer %d octets\n", taille);
    exit(EXIT_FAILURE);
  }

  /*========================================================*/
  /* Première boucle sur les éléments :                     */
  /* - comptages pour le dimensionnement des éléments       */
  /* - traitement des références des sommets                */
  /*========================================================*/

  i = 1; /* Indice ancien tableau */
  j = 1; /* Indice nouveau tableau */

  for (iel = 0 ; iel < m->ne ; iel++) {

    ncge   = m->nop5[i++] ;

    if (ncge == NOPO_NODE)
      diment = 0;
    else if (ncge == NOPO_SEGMENT)
      diment = 1;
    else if (ncge < NOPO_TETRAHEDRON)
      diment = 2;
    else if (ncge < NOPO_SUPER_ELEMENT)
      diment = 3;
    else
      diment = 0;

    /* Dimensions  */

    nmae   = m->nop5[i++] ; /* Nombre de mots pour réf. faces, arêtes, ... */
    ndsde  = m->nop5[i++] ; /* Numéro de sous-domaine */
    nno    = m->nop5[i++] ; /* Nombre de noeuds pour l'élément */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    printf("ele %d : t %d ; nmae %d ; ndsde %d ; nno %d\n",
           iel + 1, ncge, nmae, ndsde, nno) ;
#endif

    /* Connectivité de l'élément linéaire correspondant */

    nbr_som = nopo_elt_liste_c[ncge].nbr_som;

    for (isom = 0 ; isom < nno && isom < nbr_som ; isom++)
      sommets[isom] = m->nop5[i++];
    for (         ; isom < nno                   ; isom++)
      i++ ;

    /* Si sommets non confondus avec les noeuds,
       on prend les sommets (->éléments linéaires) */
    if (m->ncopnp == 0) {
      npo = m->nop5[i++] ;
      assert (npo == nbr_som) ;
      for (isom = 0 ; isom < npo ; isom++)
        sommets[isom] = m->nop5[i++] ;
    }

    /* Réinitialisation puis récupération éventuelle des références */

    for (iref = 0 ; iref < nopo_elt_liste_c[ncge].nbr_sselt ; iref++)
      tab_ref[iref] = m->nop5[i++];

    if (nmae != 0) {

      ining = m->nop5[i++] ;

      /* D'après la documention : boucle de 2 à NMAE (équiv. 0 à NMAE - 2) */

      for (iref = 0 ; iref < nmae - 1 ; iref++)
        tab_ref[iref] = m->nop5[i++];

    }

    /* renumérotation des sommets */

    for (isom = 0 ; isom < nbr_som ; isom++) {
      sommets[isom] = renum_som[isom];
      if (sommets[isom] < 1)
        diment = 0;
    }

    /*
      On ignore les éléments de dimension trop faibles (i.e. arêtes
      avec maillage volumique), les super-éléments, et les
      éléments dont un sommet a été supprimé
    */

    if (diment == 0 || (diment == 1 && nvol == 3))
      break;

    /*
      Boucle sur les sous-entités pour modification éventelle
      des références
    */

    if (nmae != 0) {

      /* En fonction de INING, on compte le nombre sous-éléments référencés */

      if (   (ining  < 3 && diment == 2)
          || (ining == 1 && diment == 3)) {

        nbr_sselt = nopo_elt_liste_c[ncge].nbr_sselt;

        for (isselt = 0 ; isselt < nbr_sselt ; isselt++) {

          /* CDG des sommets du sous-élément */

          tsselt = nopo_elt_liste_c[ncge].sous_elt[isselt].elt_typ;

          cdg[0] = 0.0;
          cdg[1] = 0.0;
          cdg[2] = 0.0;

          for (isom = 0 ;
               isom < nopo_elt_liste_c[tsselt].nbr_som;
               isom++) {
            iloc = nopo_elt_liste_c[ncge].sous_elt[isselt].som[isom] - 1 ;
            cdg[0] += m->coord[sommets[iloc]*3];
            cdg[0] += m->coord[sommets[iloc]*3 + 1];
            cdg[0] += m->coord[sommets[iloc]*3 + 2];
          }

          cdg[0] /= nopo_elt_liste_c[tsselt].nbr_som;
          cdg[1] /= nopo_elt_liste_c[tsselt].nbr_som;
          cdg[2] /= nopo_elt_liste_c[tsselt].nbr_som;

          references[isselt] = nopo_modif_ref(cdg,
                                              tab_ref[isselt]);

        }

      }

    }

  }

}


/*----------------------------------------------------------------------------
 *  Écriture d'un fichier NOPO (Format INRIA utilisé par Simail)
 *----------------------------------------------------------------------------*/

static void nopo_ecrit_maillage
(
 const nopo_maillage_t        m,                 /* --> Structure de maillage */
 const char            *const nom_fic_maillage   /* --> Nom du fichier écrit */
)
{
  FILE        *fic_maillage;         /* Descripteur du fichier */
  int_32_t    dim_e;                 /* Dimension spatiale */

  int_32_t     ind;
  int_32_t     ind_test;

  int_32_t     lnop1; /* nombre de mots dans NOP1 */
  int_32_t     lnop4; /* nombre de mots dans NOP4 */

  int_32_t  nopo_e[1 + 13];
  int_32_t  nop0[1 + 32];
  int_32_t  nop2[1 + 27];

  int_32_t  *nop1 = NULL;
  real_32_t *nop4 = NULL;

  int_32_t  ntacoo; /* type de système de coordonnées :
                       1 : cartésien ; 2 : cyl ; 3 : sphérique */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Affichage du titre */

  printf("\n\n"
         "Écriture du fichier de maillage au format NOPO (Simail)\n"
         "-------------------------------\n") ;

  printf("  Fichier de maillage : %s\n\n\n", nom_fic_maillage) ;


  /* Ouverture du fichier NOPO en écriture */

  fic_maillage = fopen(nom_fic_maillage, "wb");

  if (fic_maillage == NULL)
     nopo_err_fic(fic_maillage);

  /*  On écrit l'entête du fichier donnant les dimensions des tableaux */
  /*-------------------------------------------------------------------*/

  /*
    Attention, la première valeur d'un tableau NOP* donne la
    dimension du tableau, et le tableau est donc de taille n+1 ;
    les valeurs utiles étant aux positions 1 à n : on utilise donc
    une indexation entre 1 et n plutôt qu'entre 0 et n-1 pour
    y accéder.
  */

  lnop4 = m.np * 3;

  nopo_e[ 0] = 13;
  nopo_e[ 1] =  6;
  nopo_e[ 2] = 32;
  nopo_e[ 3] =  0;
  nopo_e[ 4] = 27;
  nopo_e[ 5] =  0;
  nopo_e[ 6] =  lnop4;
  nopo_e[ 7] =  m.lnop5;
  nopo_e[ 8] =  1;
  nopo_e[ 9] =  1;
  nopo_e[10] =  1;
  nopo_e[11] =  1;
  nopo_e[12] =  2;
  nopo_e[13] =  1;

  /* Écriture du tableau d'entête */

  nopo_ecr_rec(fic_maillage, NOPO_TAILLE_E, (char *)nopo_e);

  /* Écriture du tableau NOP0 */
  /* ------------------------ */

  /* Titre */

  strncpy((char *)(nop0 + 1), "Genere par mod_nopo", 80);
  for (ind = strlen((char *)(nop0 + 1)) ; ind < 80 ; ind++)
    ((char *)(nop0 + 1))[ind] = ' ';

  /* Date */

  {
    struct tm *dh;
    time_t    dh0;
    char   bufdh[5];

    if (time(&dh0) != -1) {
      dh = localtime(&dh0);
      sprintf(bufdh, "%02d", dh->tm_mday);
      strncpy(((char *)(nop0 + 1 +20)),     bufdh, 2);
      sprintf(bufdh, "%02d", dh->tm_mon + 1);
      strncpy(((char *)(nop0 + 1 +20)) + 2, bufdh, 2);
      sprintf(bufdh, "%04d", (dh->tm_year) % 100 + 2000);
      strncpy(((char *)(nop0 + 1 +20)) + 4, bufdh, 4);
    }
    else {
      strncpy(((char *)(nop0 + 1 +20)),     "00000000", 8);
    }
  }

  /* Nom de l'utilisateur */

  {
    char  str_user[L_cuserid];

    if ((char *)(cuserid(str_user)) == NULL)
      strcpy(str_user, "") ;

    strncpy((char *)(nop0 + 1 +22), str_user, 24);
    for (ind = strlen(str_user) ; ind < 24 ; ind++)
      ((char *)(nop0 + 1 +22))[ind] = ' ';
  }

  ((char *)(nop0 + 1 + 28))[0] = 'N' ;
  ((char *)(nop0 + 1 + 28))[1] = 'O' ;
  ((char *)(nop0 + 1 + 28))[2] = 'P' ;
  ((char *)(nop0 + 1 + 28))[3] = 'O' ;

  /* Autres enregistrements */

  nop0[30] = 4 ; /* NIVEAU (3 ou 4 sur exemples générés par Simail) */
  nop0[31] = 0 ; /* ETAT */
  nop0[32] = 0 ; /* NTACM */

  nopo_ecr_rec(fic_maillage, NOPO_TAILLE_0, (char *)nop0);

  /* Écriture du tableau NOP2 */
  /* ------------------------ */

  nop2[ 0] = 27;         /* dimension du tableau */
  nop2[ 1] =  3;         /* dimension du maillage */
  nop2[ 2] =  m.ndsr;    /* numéro de référence maximal */
  nop2[ 3] =  m.ndsd;    /* numéro de sous-domaine maximal */
  nop2[ 4] =  m.ncopnp;  /* 1 si sommets = noeuds, 0 sinon */
  nop2[ 5] =  m.ne;      /* nombre d'éléments */
  nop2[ 6] =  m.nepo;    /* nombre d'éléments point */
  nop2[ 7] =  m.nesg;    /* nombre de segments */
  nop2[ 8] =  m.ntri;    /* nombre de triangles */
  nop2[ 9] =  m.nqua;    /* nombre de quadrangles */
  nop2[10] =  m.ntet;    /* nombre de tétraèdres */
  nop2[11] =  m.npen;    /* nombre de pentaèdres */
  nop2[12] =  m.nhex;    /* nombre d'hexaèdres */
  nop2[13] =  0;         /* nombre de super éléments */
  nop2[14] =  0;         /* nombre d'éléments de bord */
  nop2[15] =  m.noe;     /* nombre de noeuds */
  nop2[16] =  0;         /* nb. noeuds intérieurs segment ou arête */
  nop2[17] =  0;         /* nb. noeuds intérieurs triangle ou face */
  nop2[18] =  0;         /* nb. noeuds intérieurs quadrangle ou face */
  nop2[19] =  0;         /* nb. noeuds intérieurs tétraèdre */
  nop2[20] =  0;         /* nb. noeuds intérieurs pentaèdre */
  nop2[21] =  0;         /* nb. noeuds intérieurs hexaèdre */
  nop2[22] =  m.np;      /* nombre de points */
  nop2[23] =  2;         /* type de valeurs de coordonnées */
  nop2[24] =  0;         /* plus grande diff. num. noeuds même élé. */
  nop2[25] =  0;         /* nombre de super-éléments dans NOP3 */
  nop2[26] =  m.lnop5;   /* nombre de mots dans NOP5 */
  nop2[27] =  1;         /* type de système de coordonnées :
                            1 : cartésien ; 2 : cyl ; 3 : sphérique */

  printf("  Données écrites : %10d points\n"
         "                    %10d segments\n"
         "                    %10d triangles\n"
         "                    %10d quadrangles\n"
         "                    %10d tétraèdres\n"
         "                    %10d pentaèdres\n"
         "                    %10d hexaèdres\n\n", 
         m.np, m.nesg, m.ntri, m.nqua, m.ntet, m.npen, m.nhex);

  nopo_ecr_rec(fic_maillage, NOPO_TAILLE_2, (char *)nop2);


  /* Écriture du tableau NOP4 */
  /* ------------------------ */

  /*
    Ce tableau contient des réels codés sur 4 octets ; de plus, il peut
    ne pas contenir de coordonnée "z" en 2D. Comme on ne s'intéresse
    qu'aux points et non aux noeuds intérieurs, on ne conserve que
    cette partie du tableau (les noeuds viennent après).
  */

  if (sizeof(real_32_t) != 4) {
    fprintf(stderr,
            "Ce sous-programme a été compilé avec un type \"real_32_t\"\n"
            "de taille différente de 4. (Erreur de portage bloquante pour\n"
            "la lecture et écriture des coordonnées simple précision NOPO).");
    exit(EXIT_FAILURE);
  }

  nop4 = malloc((lnop4 + 1) * sizeof(real_32_t));
  if (nop4 == NULL) {
    fprintf(stderr,
            "Impossible d'allouer %d octets\n",
            (lnop4 + 1) * sizeof(real_32_t));
    exit(EXIT_FAILURE);
  }

  ((int_32_t *) nop4)[0] = lnop4;
  for (ind = 0 ; ind < lnop4 ; ind++)
    nop4[ind + 1] = m.coord[ind];

  nopo_ecr_rec(fic_maillage, (lnop4+1)*4, (char *)nop4);

  free(nop4) ;

  /* Écriture du tableau NOP5 */
  /* ------------------------ */

  nopo_ecr_rec(fic_maillage, (m.lnop5+1)*4, (char *)(m.nop5));

  /* Fermeture du fichier de maillage */

  fclose(fic_maillage);

}


/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

int main
(
 int    argc,
 char * argv[]
)
{

  int_32_t *renum_som = NULL;

  nopo_maillage_t maillage;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if  (argc < 3) {
    printf("Utilisation :\n%s nom_fic nom_fic_conv\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (strcmp(argv[1], argv[2]) == 0) {
    printf("Les fichiers lus et écrits doivent être différents\n");
    exit(EXIT_FAILURE);
  }

  maillage = nopo_lit_maillage(argv[1]);

  renum_som = nopo_traite_sommets(&maillage);

  nopo_traite_elements(renum_som, &maillage);

  nopo_ecrit_maillage(maillage, argv[2]);

  exit(EXIT_SUCCESS);

}
