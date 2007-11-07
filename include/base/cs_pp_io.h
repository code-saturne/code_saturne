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

#ifndef __CS_PP_IO_H__
#define __CS_PP_IO_H__

/*============================================================================
 *  Low file I/O utility functions for Preprocessor output
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_PP_IO_NAME_LEN   32    /* Section header nam length */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*  Input or output mode */

  typedef enum {

  CS_PP_IO_MODE_READ,
  CS_PP_IO_MODE_WRITE

} cs_pp_io_mode_t;

/* Structure associated with opaque pre-processing structure object */

typedef struct _cs_pp_io_t cs_pp_io_t;

/* Structure used to save section header data, so as to simplify
   passing this data to various functions */

typedef struct {

  char            nom_rub[CS_PP_IO_NAME_LEN + 1];     /* Nom  */
  fvm_gnum_t      nbr_elt;                            /* Nombre d'éléments */
  cs_type_t       typ_elt;                            /* Type si nbr_elt > 0 */
  fvm_datatype_t  typ_lu;                             /* Type dans le fichier */

} cs_pp_io_msg_header_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Global structure associated with pre-processor output (kernel input) file */

extern cs_pp_io_t  *cs_glob_pp_io;

/*============================================================================
 * Public function prototypes
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
);

/*----------------------------------------------------------------------------
 *  Fonction qui termine une lecture ou écriture
 *----------------------------------------------------------------------------*/

cs_pp_io_t * cs_pp_io_finalize
(
 cs_pp_io_t *pp_io
);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur le nom d'un pré traitemenent
 *----------------------------------------------------------------------------*/

const char * cs_pp_io_get_name
(
 const cs_pp_io_t  *const pp_io
);

/*----------------------------------------------------------------------------
 *  Lecture de l'entete d'un message.
 *----------------------------------------------------------------------------*/

void cs_pp_io_read_header
(
       cs_pp_io_msg_header_t  *const entete,
 const cs_pp_io_t             *const pp_io
);

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
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_PP_IO_H__ */
