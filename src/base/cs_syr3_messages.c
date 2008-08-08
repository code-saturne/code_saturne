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

/*============================================================================
 * Manage messages for Syrthes coupling: sending, receiving and interpolation
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_syr3_comm.h"
#include "cs_syr3_coupling.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_syr3_messages.h"

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

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if Syrthes coupling continues or if we must finalize communications.
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTSYR (IMSFIN)
 * *****************
 *
 * INTEGER          IMSFIN      : <-- : Indicates "end" message
 * INTEGER          NTMABS      : <-> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration numbern
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstsyr, TSTSYR)
(
 cs_int_t *const imsfin,
 cs_int_t *const ntmabs,
 cs_int_t *const ntcabs
)
{
  cs_int_t  ii, i_coupl;
  cs_syr3_comm_msg_entete_t  header;
  char  section_name[32];

  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_syr3_coupling_t *syr_coupling = NULL;
  cs_syr3_comm_t  *comm = NULL;

  *imsfin = 0;

  if (n_coupl > 0) {

    for (i_coupl = 0; i_coupl < n_coupl; i_coupl++) {

      *imsfin = 0;

      syr_coupling = cs_syr3_coupling_by_id(i_coupl);
      comm = cs_syr3_coupling_get_recv_comm(syr_coupling);

      cs_syr3_comm_recoit_entete(&header, comm);

      assert(header.num_rub == 0);

      for (ii = 0 ; ii < 32 ; ii++)
        section_name[ii] = header.nom_rub[ii];

      /* Treatment according to the received header */

      if (!strncmp(CS_SYR3_COMM_CMD_ARRET, section_name,
                   strlen(CS_SYR3_COMM_CMD_ARRET) )) {

        *ntmabs = *ntcabs;
        *imsfin = 1;

        cs_base_warn(__FILE__, __LINE__);
        bft_printf(_("========================================================\n"
                     "   ** ARRET PAR DEMANDE DE SYRTHES \n"
                     "      -------------------------------------- \n"
                     "      MESSAGE RECU    : \"%s\"\n"
                     "========================================================\n"),
                   section_name);

      }
      else if (!strncmp(CS_SYR3_COMM_FIN_FICHIER, section_name,
                        strlen(CS_SYR3_COMM_FIN_FICHIER) )) {

        *ntmabs = *ntcabs;
        *imsfin = 1;

        cs_base_warn(__FILE__, __LINE__);
        bft_printf
          (_("========================================================\n"
             "   ** ARRET PAR DEMANDE DE SYRTHES \n"
             "      -------------------------------------- \n"
             "      MESSAGE RECU    : \"%s\"\n"
             "========================================================\n"),
           section_name);

      }
      else if (strncmp(CS_SYR3_COMM_CMD_ITER_DEB, section_name,
                       strlen(CS_SYR3_COMM_CMD_ITER_DEB) )) {

        bft_error
          (__FILE__, __LINE__, 0,
           _("========================================================\n"
             "   ** MESSAGE INATTENDU DANS TSTSYI \n"
             "      -------------------------------------- \n"
             "      MESSAGE RECU    : \"%s\"\n"
             "      MESSAGE ATTENDU : cmd:iter:deb \n"
             "========================================================\n"),
           section_name);
      }

      assert(header.nbr_elt == 0);

    } /* End of loop on Syrthes couplings */

  } /* If there is at least one Syrthes coupling */

  else
    return;

}

/*----------------------------------------------------------------------------
 * Synchronize new time step message.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ITDSYR (NTCABS, NTMABS)
 * *****************
 *
 * INTEGER          NTMABS      : --> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration number
 *----------------------------------------------------------------------------*/

void CS_PROCF(itdsyr, ITDSYR)
(
 cs_int_t   *const ntcabs,
 cs_int_t   *const ntmabs
)
{
  cs_int_t  i_coupl;

  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_syr3_comm_t *comm = NULL;

  /* If there is at least one syrthes coupling */

  if (n_coupl > 0) {

    /*
       Code_Saturne tells Syrthes when we are ready to begin a new
       time step, also specifying if it is the last time step.
    */

    for (i_coupl = 0; i_coupl < n_coupl; i_coupl++) {

      cs_syr3_coupling_t *coupl = cs_syr3_coupling_by_id(i_coupl);

      comm = cs_syr3_coupling_get_send_comm(coupl);

      if (*ntcabs == *ntmabs)
        cs_syr3_comm_envoie_message(0,
                                    CS_SYR3_COMM_CMD_ITER_DEB_FIN,
                                    0,
                                    CS_TYPE_char,
                                    NULL,
                                    comm);

      else if (*ntcabs < *ntmabs)
        cs_syr3_comm_envoie_message(0,
                                    CS_SYR3_COMM_CMD_ITER_DEB,
                                    0,
                                    CS_TYPE_char,
                                    NULL,
                                    comm);

      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Le numéro d'itération courant \"%d\" est supérieur au "
                    "numéro d'itération maximal demandé \"%d\" \n"),
                  *ntcabs, *ntmabs);

    }

  } /* n_coupl > 0 */

}

/*----------------------------------------------------------------------------
 * Receive coupling variables from Syrthes
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYI (NUMSYR, NOMRUB, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of Syrthes coupling
 * INTEGER          NBRENT      : --> : Number of elements
 * DOUBLE PRECISION TWALL       : <-- : Wall temerature
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyi, VARSYI)
(
 cs_int_t   *const numsyr,
 cs_int_t   *const nbrent,
 cs_real_t  *const twall
)
{
  cs_int_t  ii;
  cs_syr3_comm_msg_entete_t  header;

  char  section_name[CS_SYR3_COMM_LNG_NOM_RUB + 1];

  cs_real_t  *syr_data = NULL;

  cs_int_t  n_vertices = 0;
  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_syr3_coupling_t  *syr_coupling = NULL;
  cs_syr3_comm_t  *comm = NULL;

  if (*numsyr < 1 || *numsyr > n_coupl)
    bft_error(__FILE__, __LINE__, 0,
              _("Numéro de couplage SYRTHES %d impossible ; "
                " on a %d couplages"),
              *numsyr, n_coupl);
  else {

    syr_coupling = cs_syr3_coupling_by_id(*numsyr - 1);
    comm = cs_syr3_coupling_get_recv_comm(syr_coupling);

    n_vertices = cs_syr3_coupling_get_n_vertices(syr_coupling);

  }

  if (n_vertices > 0) {

    /* Expected section name definition */

    sprintf(section_name, "coupl:b:tparoi:%04d", *numsyr);

    for (ii = strlen(section_name); ii < CS_SYR3_COMM_LNG_NOM_RUB; ii++)
      section_name[ii] = ' ';
    section_name[CS_SYR3_COMM_LNG_NOM_RUB] = '\0';

    /* Receive header and check consistency */

    cs_syr3_comm_recoit_entete(&header, comm);

    assert(header.num_rub == 0);

    if (   strncmp(header.nom_rub, section_name, CS_SYR3_COMM_LNG_NOM_RUB) != 0
        || (header.nbr_elt > 0 && header.typ_elt != CS_TYPE_cs_real_t)
        || header.nbr_elt != n_vertices)
      bft_error(__FILE__, __LINE__, 0,
                _("Message inatendu dans le couplage SYRTHES %d :\n"
                  " on attendait \"%s\" (%d éléments, type %d)\n"
                  " on a ici     \"%s\"(%d éléments, type %d)\n"),
                *numsyr, section_name, n_vertices, (int)CS_TYPE_cs_real_t,
                header.nom_rub, header.nbr_elt, (int)(header.typ_elt));

    /* Receive data */

    BFT_MALLOC(syr_data, header.nbr_elt, cs_real_t);

    cs_syr3_comm_recoit_corps(&header,
                              syr_data,
                              comm);

    /* Transfer received fields from vertices to faces */

    cs_syr3_coupling_post_var_update(syr_coupling, 0, syr_data);

    cs_syr3_coupling_vtx_to_elt(syr_coupling,
                                syr_data,
                                twall);
  }

  if (syr_data != NULL)
    BFT_FREE(syr_data);

}

/*----------------------------------------------------------------------------
 * Send coupling variables to Syrthes
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYO (NUMSYR, NOMRUB, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of Syrthes coupling
 * INTEGER          NBRENT      : --> : Number of elements
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION HWALL       : --> : Exchange coefficient
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyo, VARSYO)
(
 cs_int_t   *const numsyr,
 cs_int_t   *const nbrent,
 cs_real_t  *const tfluid,
 cs_real_t  *const hwall
)
{
  cs_int_t  ii, var_id;
  char  section_name[CS_SYR3_COMM_LNG_NOM_RUB + 1];

  cs_int_t  n_syr_values = 0;
  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_real_t  *src_data = NULL, *syr_data = NULL;
  cs_syr3_coupling_t  *syr_coupling = NULL;
  cs_syr3_comm_t  *comm = NULL;

  if (*numsyr < 1 || *numsyr > n_coupl)
    bft_error(__FILE__, __LINE__, 0,
              _("Numéro de couplage SYRTHES %d impossible ; "
                " on a %d couplages"),
              *numsyr, n_coupl);
  else {

    syr_coupling = cs_syr3_coupling_by_id(*numsyr - 1);
    comm = cs_syr3_coupling_get_send_comm(syr_coupling);

  }

  /* loop on fluid temperature then exchange coefficient */

  for (var_id = 0; var_id < 2; var_id++) {

    if (var_id == 0) {
      sprintf(section_name, "coupl:b:tfluid:%04d:", *numsyr);
      src_data = tfluid;
    }
    else { /* if (var_id == 1) */
      sprintf(section_name, "coupl:b:hparoi:%04d:", *numsyr);
      src_data = hwall;
    }

    for (ii = strlen(section_name); ii < CS_SYR3_COMM_LNG_NOM_RUB; ii++)
      section_name[ii] = ' ';
    section_name[CS_SYR3_COMM_LNG_NOM_RUB] = '\0';

    if (*nbrent > 0) { /* Send data to Syrthes */

      /* Number of vertices in coupled mesh */

      n_syr_values = cs_syr3_coupling_get_n_vertices(syr_coupling);

      /* Transfer Code_Saturne fields to vertices */

      BFT_MALLOC(syr_data, n_syr_values * 2, cs_real_t);

      cs_syr3_coupling_elt_to_vtx(syr_coupling,
                                  src_data,
                                  n_syr_values,
                                  syr_data);

      cs_syr3_coupling_post_var_update(syr_coupling, 1 + var_id, syr_data);

      /* Send data to Syrthes */

      cs_syr3_comm_envoie_message(0,
                                  section_name,
                                  n_syr_values,
                                  CS_TYPE_cs_real_t,
                                  syr_data,
                                  comm);

      BFT_FREE(syr_data);

    }

  } /* End of loop on var_id (fluid temperature then exchange coefficient) */

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
