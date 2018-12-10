/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"
#include "bft_mem.h"

#include "cs_log.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_thermal_model.h"
#include "cs_parameters.h"

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_rad_transfer_restart.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write Radiative restart file
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_write(void)
{
  cs_real_t tkelvi = 273.15;

  /* Open output */

  cs_log_printf(CS_LOG_DEFAULT,
                _("   ** Information on the radiative module\n"
                  "      -----------------------------------\n"
                  "    Writing a restart file\n"));

  cs_restart_t *rp = cs_restart_create ("radiative_transfer",
                                        NULL,
                                        CS_RESTART_MODE_WRITE);

  cs_log_printf(CS_LOG_DEFAULT,
                _("      Write start\n"));

  /* Headers */

  int ivers = 400000;
  cs_restart_write_section(rp,
                           "version_fichier_suite_rayonnement",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &ivers);

  cs_log_printf(CS_LOG_DEFAULT,
                _("      End of output for dimensions\n"));

  /* Temps (par securite) */

  cs_restart_write_section(rp,
                           "nbre_pas_de_temps",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &cs_glob_time_step->nt_cur);

  cs_restart_write_section(rp,
                           "instant_precedent",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_real_t,
                           &cs_glob_time_step->t_cur);

  /* Boundary values */

  cs_field_t *f_btemp = CS_F_(t_b);

  if (cs_glob_thermal_model->itpscl == 1)
    cs_restart_write_field_vals(rp, f_btemp->id, 0);

  else {
    cs_real_t *tb_save;
    BFT_MALLOC(tb_save, cs_glob_mesh->n_b_faces, cs_real_t);

    for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++)
      tb_save[ifac] = f_btemp->val[ifac] + tkelvi;

    cs_restart_write_section(rp,
                             "boundary_temperature::vals::0",
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             1,
                             CS_TYPE_cs_real_t,
                             tb_save);
    BFT_FREE(tb_save);
  }

  cs_restart_write_field_vals(rp, CS_F_(qinci)->id, 0);
  cs_restart_write_field_vals(rp, CS_F_(hconv)->id, 0);
  cs_restart_write_field_vals(rp, CS_F_(fconv)->id, 0);

  /* Cell values */

  cs_restart_write_field_vals(rp, CS_FI_(rad_est, 0)->id, 0);
  cs_restart_write_field_vals(rp, CS_FI_(rad_ist, 0)->id, 0);
  cs_restart_write_field_vals(rp, CS_F_(rad_lumin)->id, 0);

  cs_restart_write_fields(rp, CS_RESTART_RAD_TRANSFER);

  cs_log_printf(CS_LOG_DEFAULT,
                _("      End of output for data\n"));

  /* Close file */

  cs_restart_destroy(&rp);

  cs_log_printf(CS_LOG_DEFAULT,
                _("    End of output to restart file\n"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read Radiative restart file
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_read(void)
{
  if (cs_glob_rad_transfer_params->restart < 1)
    return;

  /* ====================================================================
   * 1. LECTURE DU FICHIER SUITE
   * ====================================================================   */

  /* Ouverture  */
  cs_log_printf(CS_LOG_DEFAULT,
                _("   ** INFORMATIONS SUR LE MODULE DE RAYONNEMENT\n"
                  "      ------------------------------------------  \n"
                  "    Lecture d''un fichier suite\n"));

  cs_restart_t *rp = cs_restart_create("radiative_transfer", NULL, CS_RESTART_MODE_READ);

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"));

  /* Type de fichier suite     */
  /*    Pourrait porter le numero de version si besoin. */
  /*    On ne se sert pas de IVERS pour le moment  */
  int ivers;
  {
    char rubriq[64];
    snprintf(rubriq, 63, "version_fichier_suite_rayonnement");
    rubriq[63] = '\0';
    int ierror = cs_restart_read_section(rp,
                                         rubriq,
                                         CS_MESH_LOCATION_NONE,
                                         1,
                                         CS_TYPE_cs_int_t,
                                         &ivers);
    if (ierror != 0)
      bft_error(__FILE__, __LINE__, 0,
                "@\n"
                "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                "@\n"
                "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
                "@    =========                                    RAYONNEMENT\n"
                "@      TYPE DE FICHIER INCORRECT\n"
                "@\n"
                "@    Le fichier %13s ne semble pas etre un fichier\n"
                "@      suite rayonnement.\n"
                "@\n"
                "@    Le calcul ne peut etre execute.\n"
                "@\n"
                "@    Verifier que le fichier suite utilise correspond bien\n"
                "@        a un fichier suite rayonnement.\n"
                "@\n"
                "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                "@\n",
                "radiative_transfer");
  }

  /*  Tests     */
  int iok  = 0;

  /* Dimensions des supports   */
  bool ncelok, nfaiok, nfabok, nsomok;
  cs_restart_check_base_location(rp, &ncelok, &nfaiok, &nfabok, &nsomok);

  if (! ncelok) {
    cs_log_printf(CS_LOG_DEFAULT,
                  "@\n"
                  "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                  "@\n"
                  "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
                  "@    =========   RAYONNEMENT\n"
                  "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
                  "@\n"
                  "@    Le nombre de cellules a ete modifie\n"
                  "@\n"
                  "@    Le calcul ne peut etre execute.\n"
                  "@\n"
                  "@\n"
                  "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                  "@\n");
    iok++;
  }

  if (! nfabok) {
    cs_log_printf(CS_LOG_DEFAULT,
                  "@\n"
                  "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                  "@\n"
                  "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
                  "@    =========   RAYONNEMENT\n"
                  "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
                  "@\n"
                  "@    Le nombre de faces de bord a ete modifie\n"
                  "@\n"
                  "@    Le calcul ne peut etre execute.\n"
                  "@\n"
                  "@\n"
                  "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                  "@\n");
    iok++;
  }
  if (iok != 0)
    cs_exit(1);

  /* ---> Pour test ulterieur si pb : arret   */
  int ierror;
  int nberro = 0;

  /* Boundary faces    */
  {
    cs_field_t *f_btemp = CS_F_(t_b);
    char rubriq[64];
    snprintf(rubriq, 63, "boundary_temperature::vals::0");
    rubriq[63] = '\0';
    char old_name[64];
    snprintf(old_name, 63, "wall_temperature");
    old_name[63] = '\0';

    ierror = cs_restart_read_section_compat(rp,
                                            rubriq,
                                            old_name,
                                            CS_MESH_LOCATION_BOUNDARY_FACES,
                                            1,
                                            CS_TYPE_cs_real_t,
                                            f_btemp->val);
    nberro += ierror;

    if (cs_glob_thermal_model->itpscl == 2) {
      for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++)
        f_btemp->val[ifac] -= 273.15;
    }
  }

  ierror = cs_restart_read_field_vals(rp, CS_F_(qinci)->id, 0);
  nberro += ierror;

  ierror = cs_restart_read_field_vals(rp, CS_F_(hconv)->id, 0);
  nberro += ierror;

  ierror = cs_restart_read_field_vals(rp, CS_F_(fconv)->id, 0);
  nberro += ierror;

  /* Cells */

  ierror = cs_restart_read_field_vals(rp, CS_FI_(rad_est, 0)->id, 0);
  nberro += ierror;
  ierror = cs_restart_read_field_vals(rp,  CS_FI_(rad_ist, 0)->id, 0);
  nberro += ierror;
  ierror = cs_restart_read_field_vals(rp, CS_F_(rad_lumin)->id, 0);
  nberro += ierror;

  cs_restart_read_fields(rp, CS_RESTART_RAD_TRANSFER);

  /* --> Si pb : arret    */

  if (nberro != 0)
    bft_error(__FILE__, __LINE__, 0,
              "Error(s) reading radiative restart.");

  cs_restart_destroy(&rp);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("    Finished reading radiative restart file.\n"));
  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       "-------------------------------------------------------------\n"));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
