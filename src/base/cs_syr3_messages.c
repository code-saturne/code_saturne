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

BEGIN_C_DECLS

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
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if Syrthes coupling continues or if we must finalize communications.
 *
 * parameters:
 *   is_end     --> "end" message indicator
 *   nt_cur_abs <-- current iteration number
 *   nt_max_abs <-> maximum iteration number
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_test_iter(int  *is_end,
                           int   nt_cur_abs,
                           int  *nt_max_abs)
{
  cs_int_t  ii, i_coupl;
  cs_syr3_comm_msg_entete_t  header;
  char  section_name[32 + 1];

  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_syr3_coupling_t *syr_coupling = NULL;
  cs_syr3_comm_t  *comm = NULL;

  *is_end = 0;

  if (n_coupl > 0) {

    for (i_coupl = 0; i_coupl < n_coupl; i_coupl++) {

      *is_end = 0;

      syr_coupling = cs_syr3_coupling_by_id(i_coupl);
      comm = cs_syr3_coupling_get_recv_comm(syr_coupling);

      cs_syr3_comm_recoit_entete(&header, comm);

      for (ii = 0 ; ii < 32 ; ii++)
        section_name[ii] = header.nom_rub[ii];
      section_name[32] = '\0';

      /* Treatment according to the received header */

      if (!strncmp("cmd:stop", section_name,
                   strlen("cmd:stop"))) {

        *nt_max_abs = nt_cur_abs;
        *is_end = 1;

        cs_base_warn(__FILE__, __LINE__);
        bft_printf
          (_("========================================================\n"
             "   ** Abort on SYRTHES request\n"
             "      ------------------------\n"
             "      received message: \"%s\"\n"
             "========================================================\n"),
                   section_name);

      }
      else if (!strncmp(CS_SYR3_COMM_FIN_FICHIER, section_name,
                        strlen(CS_SYR3_COMM_FIN_FICHIER))) {

        *nt_max_abs = nt_cur_abs;
        *is_end = 1;

        cs_base_warn(__FILE__, __LINE__);
        bft_printf
          (_("========================================================\n"
             "   ** Abort on SYRTHES request\n"
             "      ------------------------\n"
             "      received message: \"%s\"\n"
             "========================================================\n"),
           section_name);

      }
      else if (strncmp("cmd:iter:start", section_name,
                       strlen("cmd:iter:start"))) {

        bft_error
          (__FILE__, __LINE__, 0,
           _("========================================================\n"
             "   ** Unexpected message in cs_syr3_messages_test_iter\n"
             "      ------------------------------------------------\n"
             "      received message: \"%s\"\n"
             "      expected message: cmd:iter:start \n"
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
 * Synchronize new time step
 *
 * parameters:
 *   nt_cur_abs <-- current iteration number
 *   nt_max_abs --> maximum iteration number
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_new_time_step(int  nt_cur_abs,
                               int  nt_max_abs)
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

      if (nt_cur_abs == nt_max_abs)
        cs_syr3_comm_envoie_message("cmd:iter:start:last",
                                    0,
                                    CS_TYPE_char,
                                    NULL,
                                    comm);

      else if (nt_cur_abs < nt_max_abs)
        cs_syr3_comm_envoie_message("cmd:iter:start",
                                    0,
                                    CS_TYPE_char,
                                    NULL,
                                    comm);

      else
        bft_error(__FILE__, __LINE__, 0,
                  _("The current iteration number \"%d\" is greater than "
                    "the requested\n"
                    "maximum number iteration \"%d\" \n"),
                  nt_cur_abs, nt_max_abs);

    }

  } /* n_coupl > 0 */

}

/*----------------------------------------------------------------------------
 * Receive coupling variables from Syrthes
 *
 * parameters:
 *   syr_num <-- Syrthes 3 coupling number
 *   twall   --> wall temperature
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_recv_twall(cs_int_t   syr_num,
                            cs_real_t  twall[])
{
  cs_int_t  ii;
  cs_syr3_comm_msg_entete_t  header;

  char  section_name[CS_SYR3_COMM_H_LEN + 1];

  cs_real_t  *syr_data = NULL;

  cs_int_t  n_vertices = 0;
  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_syr3_coupling_t  *syr_coupling = NULL;
  cs_syr3_comm_t  *comm = NULL;

  if (syr_num < 1 || syr_num > n_coupl)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              syr_num, n_coupl);
  else {

    syr_coupling = cs_syr3_coupling_by_id(syr_num - 1);
    comm = cs_syr3_coupling_get_recv_comm(syr_coupling);

    n_vertices = cs_syr3_coupling_get_n_vertices(syr_coupling);

  }

  if (n_vertices > 0) {

    /* Expected section name definition */

    sprintf(section_name, "coupl:b:tparoi");

    for (ii = strlen(section_name); ii < CS_SYR3_COMM_H_LEN; ii++)
      section_name[ii] = ' ';
    section_name[CS_SYR3_COMM_H_LEN] = '\0';

    /* Receive header and check consistency */

    cs_syr3_comm_recoit_entete(&header, comm);

    if (   strncmp(header.nom_rub, section_name, CS_SYR3_COMM_H_LEN) != 0
        || (header.nbr_elt > 0 && header.typ_elt != CS_TYPE_cs_real_t)
        || header.nbr_elt != n_vertices)
      bft_error(__FILE__, __LINE__, 0,
                _("Unexpected message in the SYRTHES coupling %d:\n"
                  " expected \"%s\" (%d elements, type %d)\n"
                  " received \"%s\" (%d elements, type %d)\n"),
                syr_num, section_name, n_vertices, (int)CS_TYPE_cs_real_t,
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
 * parameters:
 *   syr_num <-- Syrthes 3 coupling number
 *   tfluid  <-- wall exchange coefficient
 *   hwall   <-- wall exchange coefficient
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_send_tf_hwall(cs_int_t   syr_num,
                               cs_real_t  tfluid[],
                               cs_real_t  hwall[])
{
  cs_int_t  ii, var_id;
  char  section_name[CS_SYR3_COMM_H_LEN + 1];

  cs_int_t  n_coupl = cs_syr3_coupling_n_couplings();
  cs_int_t  n_syr_values = 0;
  cs_real_t  *src_data = NULL, *syr_data = NULL;
  cs_syr3_coupling_t  *syr_coupling = NULL;
  cs_syr3_comm_t  *comm = NULL;

  if (syr_num < 1 || syr_num > n_coupl)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              syr_num, n_coupl);
  else {

    syr_coupling = cs_syr3_coupling_by_id(syr_num - 1);
    comm = cs_syr3_coupling_get_send_comm(syr_coupling);

    n_syr_values = cs_syr3_coupling_get_n_vertices(syr_coupling);

  }

  if (n_syr_values == 0)
    return;

  /* loop on fluid temperature then exchange coefficient */

  for (var_id = 0; var_id < 2; var_id++) {

    if (var_id == 0) {
      sprintf(section_name, "coupl:b:tfluid");
      src_data = tfluid;
    }
    else { /* if (var_id == 1) */
      sprintf(section_name, "coupl:b:hparoi");
      src_data = hwall;
    }

    for (ii = strlen(section_name); ii < CS_SYR3_COMM_H_LEN; ii++)
      section_name[ii] = ' ';
    section_name[CS_SYR3_COMM_H_LEN] = '\0';

    /* Transfer Code_Saturne fields to vertices */

    BFT_MALLOC(syr_data, n_syr_values * 2, cs_real_t);

    cs_syr3_coupling_elt_to_vtx(syr_coupling,
                                src_data,
                                n_syr_values,
                                syr_data);

    cs_syr3_coupling_post_var_update(syr_coupling, 1 + var_id, syr_data);

    /* Send data to Syrthes */

    cs_syr3_comm_envoie_message(section_name,
                                n_syr_values,
                                CS_TYPE_cs_real_t,
                                syr_data,
                                comm);

    BFT_FREE(syr_data);

  } /* End of loop on var_id (fluid temperature then exchange coefficient) */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
