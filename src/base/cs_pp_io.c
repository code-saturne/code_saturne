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

/*============================================================================
 *  Low file I/O utility functions for Preprocessor output
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_file.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_pp_io.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 *  Structures locales
 *============================================================================*/

struct _cs_pp_io_t {

  char            *nom;          /* Nom du communicateur                      */

  bft_file_t      *fic;          /* Pointeur sur fichier associé              */

  cs_pp_io_mode_t  mode;         /* Mode de communication                     */
  cs_bool_t        swap_endian;  /* Permutation des octets ?                  */
  cs_int_t         echo;         /* Niveau d'impression des donnees           */

};

/*============================================================================
 *  Constantes et Macros
 *============================================================================*/


/*============================================================================
 *  Variables globales statiques
 *============================================================================*/

static char  cs_pp_io_nom_typ_elt_char[] = "c ";  /* Chaîne de caractères */
static char  cs_pp_io_nom_typ_elt_i4[] =   "i4";  /* Entier signé 32 bits */
static char  cs_pp_io_nom_typ_elt_i8[] =   "i8";  /* Entier signé 64 bits */
static char  cs_pp_io_nom_typ_elt_u4[] =   "u4";  /* Entier non signé 32 bits */
static char  cs_pp_io_nom_typ_elt_u8[] =   "u8";  /* Entier non signé 64 bits */
static char  cs_pp_io_nom_typ_elt_r4[] =   "r4";  /* Réel simple précision */
static char  cs_pp_io_nom_typ_elt_r8[] =   "r8";  /* Réel double précision */

/* Pointeurs globaux sur structure de lecture des données préprocesseur */
cs_pp_io_t  *cs_glob_pp_io = NULL;

/*============================================================================
 *  Prototypes de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui construit le descripteur du fichier d'interface et initialise
 *  ce fichier par l'envoi ou la lecture d'une éventuelle "chaîne magique"
 *  servant a vérifier le bon format des fichiers
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_fic_ouvre
(
       cs_pp_io_t  *const pp_io,
 const char        *const nom,
 const char        *const chaine_magique
);


/*----------------------------------------------------------------------------
 *  Fonction qui ferme le fichier d'interface
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_fic_ferme
(
 cs_pp_io_t  *pp_io
);


/*----------------------------------------------------------------------------
 *  Affichage de l'attente d'échange d'un message
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_echo_pre
(
 const cs_pp_io_t  *const pp_io
);


/*----------------------------------------------------------------------------
 *  Affichage de l'entete d'un message
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_echo_entete
(
 const char        *nom_rub,
       size_t       nbr_elt,
       cs_type_t    typ_elt,
       cs_type_t    typ_lu
);


/*----------------------------------------------------------------------------
 *  Affichage (partiel) du contenu d'un message
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_echo_donnees
(
 const cs_int_t          echo,
 const cs_int_t          nbr_elt,
 const cs_type_t         typ_elt,
 const void       *const elt_rub
);


/*----------------------------------------------------------------------------
 *  Conversion de données lues.
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_convert_read
(
 void            *buffer,
 void            *dest,
 size_t           n_elts,
 fvm_datatype_t   buffer_type,
 cs_type_t        dest_type
);


/*============================================================================
 *  Définitions de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui initialise une lecture ou écriture
 *----------------------------------------------------------------------------*/

cs_pp_io_t * cs_pp_io_initialize
(
 const char          *const nom_rep,        /* --> nom du répertoire associé  */
 const char          *const chaine_magique, /* --> Chaîne de vérif. de type   */
 const cs_pp_io_mode_t      mode,           /* --> Émission ou réception      */
 const cs_int_t             echo            /* --> Écho sur sortie principale
                                                    (< 0 si aucun, entête si 0,
                                                    n premiers et derniers
                                                    éléments si n)            */
)
{
  unsigned    int_endian;
  int         numero;

  cs_pp_io_t  *pp_io = NULL;

  BFT_MALLOC(pp_io, 1, cs_pp_io_t);

  /* Construction du nom */

  if (cs_glob_base_rang < 0)
    numero = 1;
  else
    numero = cs_glob_base_rang + 1;

  BFT_MALLOC(pp_io->nom, strlen(nom_rep) + strlen("/n00001") + 1, char);

  sprintf(pp_io->nom, "%s/n%05d", nom_rep, numero);


  /* Initialisation des autres champs */

  pp_io->mode = mode;
  pp_io->echo = echo;

  pp_io->fic  = NULL;


  /* Test si système "big-endian" ou "little-endian" */

  pp_io->swap_endian = CS_FALSE;

  int_endian = 0;
  *((char *)(&int_endian)) = '\1';

  if (int_endian == 1)
    pp_io->swap_endian = CS_TRUE;

#if defined(DEBUG) && !defined(NDEBUG)

  else {
    int_endian = 0;
    *((char *)(&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert(int_endian == 1);
  }

#endif

  /* Info sur la création de l'interface */

  bft_printf(_("\n  Lecture du pré traitement :  %s"), nom_rep);
  bft_printf_flush();


  /* Création du descripteur de fichier d'interface */

  _cs_pp_io_fic_ouvre(pp_io, pp_io->nom, chaine_magique);

  return pp_io;
}


/*----------------------------------------------------------------------------
 *  Fonction qui termine une lecture ou écriture
 *----------------------------------------------------------------------------*/

cs_pp_io_t * cs_pp_io_finalize
(
 cs_pp_io_t *pp_io
)
{

  /* Info sur la fermeture du fichier d'interface */

  bft_printf(_("\n  Fin de la lecture :  %s\n"), pp_io->nom);
  bft_printf_flush();

  _cs_pp_io_fic_ferme(pp_io);

  BFT_FREE(pp_io->nom);
  BFT_FREE(pp_io);

  return NULL;

}


/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur le nom d'un pré-traitement
 *----------------------------------------------------------------------------*/

const char * cs_pp_io_get_name
(
 const cs_pp_io_t  *const pp_io
)
{
  assert(pp_io != NULL);

  return(pp_io->nom);
}


/*----------------------------------------------------------------------------
 *  Lecture de l'entête d'un message.
 *----------------------------------------------------------------------------*/

void cs_pp_io_read_header
(
       cs_pp_io_msg_header_t  *const entete,
 const cs_pp_io_t             *const pp_io
)
{
  char    nom_typ_elt[] = "  ";

  assert(pp_io  != NULL);

  entete->nbr_elt = 0;

  if (pp_io->echo >= 0)
    _cs_pp_io_echo_pre(pp_io);


  /* Lecture sur fichier */
  /*---------------------*/

  /* nom de type de la rubrique */

  bft_file_read(entete->nom_rub, 1, CS_PP_IO_NAME_LEN, pp_io->fic);

  entete->nom_rub[CS_PP_IO_NAME_LEN] = '\0';

  /* nombre d'éléments */

  if (sizeof(unsigned long) == 8) {
    unsigned long rec;
    bft_file_read(&rec, 8, 1, pp_io->fic);
    entete->nbr_elt = rec;
  }
  else if (sizeof(unsigned long long) == 8) {
    unsigned long long rec;
    bft_file_read(&rec, 8, 1, pp_io->fic);
    entete->nbr_elt = rec;
  }

  assert(sizeof(unsigned long) == 8 || sizeof(unsigned long long) == 8);


  if (entete->nbr_elt != 0)
    bft_file_read(nom_typ_elt, 1, 2, pp_io->fic);


  entete->nom_rub[CS_PP_IO_NAME_LEN] = '\0';

  if (entete->nbr_elt != 0) {

    if (   strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_i4) == 0
        || strcmp(nom_typ_elt, "i ") == 0) {
      entete->typ_elt = CS_TYPE_cs_int_t;
      entete->typ_lu = FVM_INT32;
    }

    else if (strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_i8) == 0) {
      entete->typ_elt = CS_TYPE_cs_int_t;
      entete->typ_lu = FVM_INT64;
    }

    else if (strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_u4) == 0) {
      entete->typ_elt = CS_TYPE_cs_int_t;
      entete->typ_lu = FVM_UINT32;
    }

    else if (strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_u8) == 0) {
      entete->typ_elt = CS_TYPE_cs_int_t;
      entete->typ_lu = FVM_UINT64;
    }

    else if (strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_r4) == 0) {
      entete->typ_elt = CS_TYPE_cs_real_t;
      entete->typ_lu = FVM_FLOAT;
    }

    else if (strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_r8) == 0) {
      entete->typ_elt = CS_TYPE_cs_real_t;
      entete->typ_lu = FVM_DOUBLE;
    }

    else if (strcmp(nom_typ_elt, cs_pp_io_nom_typ_elt_char) == 0) {
      entete->typ_elt = CS_TYPE_char;
      entete->typ_lu = FVM_DATATYPE_NULL;
    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Erreur à la lecture du fichier de pré traitement : "
                  "\"%s\".\n"
                  "Le type de données \"%s\" n'est pas reconnu."),
                pp_io->nom, nom_typ_elt);

  }


  /* Affichage eventuel */

  if (pp_io->echo >= 0)
    _cs_pp_io_echo_entete(entete->nom_rub,
                          entete->nbr_elt,
                          entete->typ_elt,
                          entete->typ_lu);

}


/*----------------------------------------------------------------------------
 *  Lecture du corps d'un message.
 *
 *  Si la zone mémoire destinée à recevoir les données existe deja, on
 *  fournit un pointeur "elt" sur cette zone ; la fonction renvoie alors
 *  ce même pointeur. Sinon (si "elt" est à NULL), la mémoire est allouée
 *  ici, et la fonction renvoie un pointeur sur cette zone.
 *----------------------------------------------------------------------------*/

void * cs_pp_io_read_body
(
 const cs_pp_io_msg_header_t  *const entete,
       void                   *const elt,
 const cs_pp_io_t             *const pp_io
)
{

  size_t      type_size;
  cs_bool_t   convert_type = CS_FALSE;
  void      *_elt_rub = NULL;

  assert(pp_io  != NULL);

  switch(entete->typ_lu) {

  case FVM_INT32:
  case FVM_UINT32:
  case FVM_FLOAT:
    type_size = 4;
    break;

  case FVM_INT64:
  case FVM_UINT64:
  case FVM_DOUBLE:
    type_size = 8;
    break;

  default:
    type_size = 1;

  }


  _elt_rub = elt;

  if (_elt_rub == NULL && entete->nbr_elt != 0) {

    switch(entete->typ_elt) {

    case CS_TYPE_cs_int_t:
      {
        cs_int_t  *elt_rub_int;

        BFT_MALLOC(elt_rub_int, entete->nbr_elt, cs_int_t);
        _elt_rub = (void *) elt_rub_int;
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        cs_real_t  *elt_rub_rea;

        BFT_MALLOC(elt_rub_rea, entete->nbr_elt, cs_real_t);
        _elt_rub = (void *)elt_rub_rea;
      }
      break;

    case CS_TYPE_char:
      {
        char  *elt_rub_cha;

        BFT_MALLOC(elt_rub_cha, entete->nbr_elt + 1, char);
        _elt_rub = (void *)elt_rub_cha;
      }
      break;

    default:
      assert(0);
      break;

    }

  }

  /* valeurs des éléments */

  if (entete->nbr_elt != 0) {

    switch(entete->typ_elt) {

    case CS_TYPE_cs_int_t:
      {
        if (   (sizeof(cs_int_t) == 4 && entete->typ_lu == FVM_INT32)
            || (sizeof(cs_int_t) == 4 && entete->typ_lu == FVM_UINT32)
            || (sizeof(cs_int_t) == 8 && entete->typ_lu == FVM_INT64)
            || (sizeof(cs_int_t) == 8 && entete->typ_lu == FVM_UINT64))
          bft_file_read(_elt_rub, sizeof(cs_int_t), entete->nbr_elt, pp_io->fic);
        else
          convert_type = CS_TRUE;
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        if (   (sizeof(cs_real_t) == 4 && entete->typ_lu == FVM_FLOAT)
            || (sizeof(cs_real_t) == 8 && entete->typ_lu == FVM_DOUBLE))
          bft_file_read(_elt_rub, sizeof(cs_real_t), entete->nbr_elt, pp_io->fic);
        else
          convert_type = CS_TRUE;
      }
      break;

    case CS_TYPE_char:
      {
        size_t ii;
        bft_file_read(_elt_rub, 1, entete->nbr_elt, pp_io->fic);
        for (ii = 0 ;
             ii < entete->nbr_elt && ((char *)_elt_rub)[ii] != '\0' ;
             ii++);
        ((char *)_elt_rub)[ii] = '\0';
      }
      break;

    default:
      assert(0);
      break;

    }

  }

  /* Lecture avec conversion si nécessaire */

  if (convert_type == CS_TRUE) {

    void  *read_buf = NULL;
    size_t buf_size = CS_MIN(entete->nbr_elt, 131072);
    size_t n_remaining, read_size;

    BFT_MALLOC(read_buf, buf_size*type_size, char);

    for (n_remaining = entete->nbr_elt;
         n_remaining > 0;
         n_remaining -= buf_size) {

      read_size = CS_MIN(n_remaining, buf_size);

      bft_file_read(read_buf, type_size, read_size, pp_io->fic);

      _cs_pp_io_convert_read(read_buf,
                             _elt_rub,
                             read_size,
                             entete->typ_lu,
                             entete->typ_elt);

    }

    BFT_FREE(read_buf);

  }


  /* Affichage éventuel */

  if (entete->nbr_elt != 0 && pp_io->echo > 0)
    _cs_pp_io_echo_donnees(pp_io->echo,
                           entete->nbr_elt,
                           entete->typ_elt,
                           _elt_rub);


  /* Transmission des valeurs lues */

  return _elt_rub;

}


/*============================================================================
 *  Définitions de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui construit le descripteur du fichier d'interface et initialise
 *  ce fichier par l'envoi ou la lecture d'une eventuelle "chaîne magique"
 *  servant a vérifier le bon format des fichiers
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_fic_ouvre
(
       cs_pp_io_t  *const  pp_io,
 const char        *const  nom,
 const char        *const  chaine_magique
)
{

  bft_file_mode_t fic_mod_pp_io;


  /* Préparation de l'ouverture du fichier */

  switch(pp_io->mode) {

  case CS_PP_IO_MODE_READ:
    fic_mod_pp_io = BFT_FILE_MODE_READ;
    break;

  case CS_PP_IO_MODE_WRITE:
    fic_mod_pp_io = BFT_FILE_MODE_WRITE;
    break;

  default:
    assert(   pp_io->mode == CS_PP_IO_MODE_READ
           || pp_io->mode == CS_PP_IO_MODE_WRITE);

  }


  /* Création du descripteur du fichier d'interface */

  pp_io->fic = bft_file_open(nom, fic_mod_pp_io, BFT_FILE_TYPE_BINARY);
  bft_file_set_big_endian(pp_io->fic);


  /*-----------------------------------------------------*/
  /* Écriture ou lecture éventuelle d'une chaine magique */
  /*-----------------------------------------------------*/

  if (pp_io->mode == CS_PP_IO_MODE_READ) {

    char      *chaine_magique_lue;
    cs_int_t   lng_chaine_magique = strlen(chaine_magique);

    BFT_MALLOC(chaine_magique_lue, lng_chaine_magique + 1, char);

    bft_file_read(chaine_magique_lue, 1, strlen(chaine_magique), pp_io->fic);

    chaine_magique_lue[lng_chaine_magique] = '\0';

    /* Si la chaine magique ne correspond pas, on a une erreur */

    if (strcmp(chaine_magique_lue, chaine_magique) != 0) {

      bft_error(__FILE__, __LINE__, 0,
                _("Erreur à la lecture du fichier de pré traitement : "
                  "\"%s\".\n"
                  "Le format de l'interface n'est pas à la bonne version.\n"
                  "La chaîne magique repère la version du format "
                  "d'interface :\n"
                  "chaîne magique lue      : \"%s\"\n"
                  "chaîne magique actuelle : \"%s\"\n"),
                pp_io->nom, chaine_magique_lue, chaine_magique);

    }

    BFT_FREE(chaine_magique_lue);

  }
  else if (pp_io->mode == CS_PP_IO_MODE_WRITE) {

    bft_file_write(chaine_magique, 1, strlen(chaine_magique), pp_io->fic);

    bft_file_flush(pp_io->fic);

  }

}


/*----------------------------------------------------------------------------
 *  Fonction qui ferme le fichier d'interface
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_fic_ferme
(
 cs_pp_io_t  *pp_io
)
{

  pp_io->fic = bft_file_free(pp_io->fic);

}


/*----------------------------------------------------------------------------
 *  Affichage de l'attente d'échange d'un message
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_echo_pre
(
 const cs_pp_io_t  *const pp_io
)
{
  assert(pp_io != NULL);

  switch(pp_io->mode) {

  case CS_PP_IO_MODE_READ:
    bft_printf(_("\nMessage lu sur \"%s\" :\n"), pp_io->nom);
    break;

  case CS_PP_IO_MODE_WRITE:
    bft_printf(_("\nMessage écrit sur \"%s\" :\n"), pp_io->nom);
    break;

  default:
    assert(   pp_io->mode == CS_PP_IO_MODE_READ
           || pp_io->mode == CS_PP_IO_MODE_WRITE);
  }

  bft_printf_flush();

}


/*----------------------------------------------------------------------------
 *  Affichage de l'entete d'un message
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_echo_entete
(
 const char        *nom_rub,
       size_t       nbr_elt,
       cs_type_t    typ_elt,
       cs_type_t    typ_lu
)
{

  char nom_rub_ecr[CS_PP_IO_NAME_LEN + 1];

  /* instructions */

  strncpy(nom_rub_ecr, nom_rub,  CS_PP_IO_NAME_LEN);
  nom_rub_ecr[CS_PP_IO_NAME_LEN] = '\0';

  bft_printf(_("    nom de la rubrique    : \"%s\"\n"
               "    nombre d'éléments     : %d\n"),
             nom_rub_ecr, nbr_elt);

  if (nbr_elt > 0) {

    char *nom_typ;

    switch(typ_elt) {
    case CS_TYPE_char:
      nom_typ = cs_pp_io_nom_typ_elt_char;
      break;
    case CS_TYPE_cs_int_t:
      if (typ_lu == FVM_INT32)
        nom_typ = cs_pp_io_nom_typ_elt_i4;
      else if (typ_lu == FVM_INT64)
        nom_typ = cs_pp_io_nom_typ_elt_i8;
      else if (typ_lu == FVM_UINT32)
        nom_typ = cs_pp_io_nom_typ_elt_u4;
      else if (typ_lu == FVM_UINT64)
        nom_typ = cs_pp_io_nom_typ_elt_u8;
      break;
    case CS_TYPE_cs_real_t:
      if (typ_lu == FVM_FLOAT)
        nom_typ = cs_pp_io_nom_typ_elt_r4;
      else if (typ_lu == FVM_DOUBLE)
        nom_typ = cs_pp_io_nom_typ_elt_r8;
      break;
    default:
      assert(   typ_elt == CS_TYPE_char
             || typ_elt == CS_TYPE_cs_int_t
             || typ_elt == CS_TYPE_cs_real_t);
    }

    bft_printf(_("    nom du type d'élément : \"%s\"\n"), nom_typ);

  }

  bft_printf_flush();

}


/*----------------------------------------------------------------------------
 *  Affichage (partiel) du contenu d'un message
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_echo_donnees
(
 const cs_int_t          echo,
 const cs_int_t          nbr_elt,
 const cs_type_t         typ_elt,
 const void       *const elt_rub
)
{
  cs_int_t  echo_deb = 0;
  cs_int_t  echo_fin;
  cs_int_t  ind;

  /* Instructions */

  if (nbr_elt == 0) return;

  if (echo * 2 < nbr_elt) {
    echo_fin = echo;
    bft_printf(_("    %d premiers et derniers éléments :\n"), echo);
  }
  else {
    echo_fin = nbr_elt;
    bft_printf(_("    éléments :\n"));
  }

  do {

    switch (typ_elt) {

    case CS_TYPE_cs_int_t:
      {
        const cs_int_t *elt_rub_int = (const cs_int_t *) elt_rub;

        for (ind = echo_deb ; ind < echo_fin ; ind++)
          bft_printf("    %10d : %12d\n", ind + 1, *(elt_rub_int + ind));
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        const cs_real_t *elt_rub_real = (const cs_real_t *) elt_rub;

        for (ind = echo_deb ; ind < echo_fin ; ind++)
          bft_printf("    %10d : %12.5e\n", ind + 1, *(elt_rub_real + ind));
      }
      break;

    case CS_TYPE_char:
      {
        const char *elt_rub_char = (const char *) elt_rub;

        for (ind = echo_deb ; ind < echo_fin ; ind++) {
          if (*(elt_rub_char + ind) != '\0')
            bft_printf("    %10d : '%c'\n", ind + 1, *(elt_rub_char + ind));
          else
            bft_printf("    %10d : '\\0'\n", ind + 1);
        }
      }
      break;

    default:

      assert(   typ_elt == CS_TYPE_cs_int_t
             || typ_elt == CS_TYPE_cs_real_t
             || typ_elt == CS_TYPE_char);

    }

    if (echo_fin < nbr_elt) {
      bft_printf("    ..........   ............\n");
      echo_deb = nbr_elt - echo;
      echo_fin = nbr_elt;
    }
    else {
      assert(echo_fin == nbr_elt);
      echo_fin = nbr_elt + 1;
    }

  } while (echo_fin <= nbr_elt);

  bft_printf_flush();

}


/*----------------------------------------------------------------------------
 *  Conversion de données lues.
 *----------------------------------------------------------------------------*/

static void _cs_pp_io_convert_read
(
 void            *buffer,
 void            *dest,
 size_t           n_elts,
 fvm_datatype_t   buffer_type,
 cs_type_t        dest_type
)
{

  size_t  ii;

  switch(dest_type) {

    case CS_TYPE_cs_int_t:
      {
        cs_int_t *_dest = dest;

        switch(buffer_type) {

        case FVM_INT32:
          if (sizeof(int) == 4) {
            int * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          else if (sizeof(short) == 4) {
            short * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          assert(sizeof(int) == 4 || sizeof(short) == 4);
          break;

        case FVM_INT64:
          if (sizeof(long) == 8) {
            long * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          else if (sizeof(long long) == 8) {
            long long * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          assert(sizeof(long) == 8 || sizeof(long long) == 8);
          break;

        case FVM_UINT32:
          if (sizeof(unsigned) == 4) {
            unsigned * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          else if (sizeof(unsigned short) == 4) {
            unsigned short * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          assert(sizeof(unsigned) == 4 || sizeof(unsigned short) == 4);
          break;

        case FVM_UINT64:
          if (sizeof(unsigned long) == 8) {
            unsigned long * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          else if (sizeof(unsigned long long) == 8) {
            unsigned long long * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          assert(sizeof(unsigned long) == 8 || sizeof(unsigned long long) == 8);
          break;

        default:
          assert(0);

        }
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        cs_real_t *_dest = dest;

        switch(buffer_type) {

        case FVM_FLOAT:
          {
            float * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          break;

        case FVM_DOUBLE:
          {
            double * _buffer = buffer;
            for (ii = 0; ii < n_elts; ii++)
              _dest[ii] = _buffer[ii];
          }
          break;

        default:
          assert(0);

        }
      }
      break;

    default:
      assert(0);

  }

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
